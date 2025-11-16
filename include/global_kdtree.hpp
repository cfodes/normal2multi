#pragma once

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

#include "geometry.hpp"      // Point<double>
#include "partition_bvh.hpp" // AABB 定义在这里

/// --------------------------------------------------------------------------------------
/// 一个“动边界采样点”
///  - pos      : 全局坐标（来自 block_tree.grid(i)）
///  - block_id : 该点所属 mesh_block 的 block_id
///  - node_id  : 对应的全局节点 id（可选，仅用于调试或扩展）
/// --------------------------------------------------------------------------------------
struct MovingPoint
{
    Point<double> pos{};
    int           block_id = -1;
    int           node_id = -1;
    double        block_D = 0.0;
};

/// --------------------------------------------------------------------------------------
/// KD 树节点
///  - box       : 子树包围盒（用现成的 AABB）
///  - left/right: 子节点在 nodes_ 数组中的索引（-1 表示无）
///  - [begin,end): 该节点对应的 pts_ 下标范围（仅叶子节点有意义）
///  - axis/split: 划分维度 + 切分位置
///  - mono_block: 若该子树所有点的 block_id 相同，则为该 id；否则为 -1
///  - is_leaf   : 是否为叶子节点
/// --------------------------------------------------------------------------------------
struct KDNode
{
    AABB box{};
    int  left = -1;
    int  right = -1;

    int  begin = 0;
    int  end = 0;

    int    axis = 0;
    double split = 0.0;

    int  mono_block = -1;
    bool is_leaf = false;
    double min_block_D = std::numeric_limits<double>::infinity();
};

/// --------------------------------------------------------------------------------------
/// GlobalKDTree
/// 针对“当前层级的所有分区”的内部节点集合，构建一棵全局 KD 树：
///  - 采样点来自每个 mesh_block 的 block_tree（即你现在用来做距离搜索的那批点）
///  - 每个点带有所属 block_id
///  - 提供 query_point(...)：
///       给定查询点 q 和 du（到历史 unique 边界的距离）
///       返回：最近动边界距离 d1、静边界距离 d2，以及中间的 d_other/d_u
/// --------------------------------------------------------------------------------------
class GlobalKDTree
{
public:
    GlobalKDTree() : root_(-1) {}

    /// ------------------------------------------------------------------
    /// 从当前层级的 mesh_block 数组构建 KD 树
    ///
    /// 关键点：为了保证与原有 block_tree 搜索完全一致，采样点不从
    /// internal_nodes 等容器取，而是直接从
    ///   blk.block_tree.grid_size()
    ///   blk.block_tree.grid(i)
    ///   blk.block_tree.key(i)
    /// 读取。
    ///
    /// Blocks 可以是 std::vector<mesh_block> 或任何支持 range-for
    /// 且拥有成员：
    ///   int block_id;
    ///   <某个类型> block_tree;
    ///   block_tree.grid_size() -> int
    ///   block_tree.grid(i)     -> Point<double> const&
    ///   block_tree.key(i)      -> int
    /// ------------------------------------------------------------------
    template<class Blocks>
    void build_from_blocks(const Blocks& blks)
    {
        pts_.clear();

        // 预估总点数以减少 realloc
        std::size_t total = 0;
        for (const auto& b : blks) {
            const int n = b.block_tree.grid_size();
            if (n > 0) total += static_cast<std::size_t>(n);
        }
        pts_.reserve(total);

        for (const auto& b : blks) {
            const int n = b.block_tree.grid_size();
            if (n <= 0) continue;

            for (int i = 0; i < n; ++i) {
                MovingPoint mp;
                mp.pos = b.block_tree.grid(i); // 真实采样点坐标
                mp.block_id = b.block_id;
                mp.node_id = b.block_tree.key(i);   // 全局节点 id（可选）
                pts_.push_back(std::move(mp));
            }
        }

        nodes_.clear();
        if (pts_.empty()) {
            root_ = -1;
            return;
        }

        root_ = build_range(0, static_cast<int>(pts_.size()));
    }

    /// KD 树是否为空
    bool empty() const noexcept { return root_ < 0; }

    /// 采样点总数（仅用于统计/调试）
    std::size_t point_count() const noexcept { return pts_.size(); }

    /// ------------------------------------------------------------------
    /// 查询函数：给定查询点 q 和 du（到历史 unique 边界的距离）
    ///
    /// 输入：
    ///   q  : 查询点坐标
    ///   du : 到历史 unique 边界最近距离（线性距离，若不存在则传 +inf）
    ///
    /// 输出：
    ///   moving_blk_id : 最近动边界所在分区 id（block_id）
    ///   d1            : 最近动边界距离（线性）
    ///   d2            : 静边界距离 = min( d_other, d_u )（线性）
    ///   d_other       : 来自“其它分区”的最近距离（线性）
    ///   d_u           : du（直接返回，线性）
    ///
    /// 实现思路：
    ///   1) Pass1: 标准 KD-tree 最近邻，得到最近点 (d1, block_id = b1)
    ///   2) Pass2: 在同一棵树上再次遍历，但：
    ///        - 对 mono_block == b1 的子树整棵剪枝
    ///        - 仅用 block != b1 的点更新 d_other
    ///        - 同时使用 du^2 作为初始上界做剪枝
    ///   3) d2 = min( d_other, du )
    /// ------------------------------------------------------------------
    void query_point(
        const Point<double>& q,
        double du,
        int& moving_blk_id,
        double& d1,
        double& d2,
        double& d_other,
        double& d_u
    ) const
    {
        moving_blk_id = -1;
        d1 = d2 = d_other = d_u =
            std::numeric_limits<double>::infinity();

        if (root_ < 0 || pts_.empty()) {
            return; // 没有动边界点
        }

        // ============ Pass 1: 最近动边界 d1, b1 ======================
        double best1_d2 = std::numeric_limits<double>::infinity();
        int    b1 = -1;
        nearest_pass1(q, best1_d2, b1);

        if (b1 < 0 || !std::isfinite(best1_d2)) {
            // 理论上不应该发生，除非 KD 树构建失败
            return;
        }

        moving_blk_id = b1;
        d1 = std::sqrt(best1_d2);

        // du: 历史 unique_bndry 距离
        d_u = du;
        const double du2 = std::isfinite(du)
            ? du * du
            : std::numeric_limits<double>::infinity();

        // ============ Pass 2: 排除 b1，只看其它分区的最近距离 ============
        double best_other_d2 = std::numeric_limits<double>::infinity();
        double best_d2 = du2; // d2^2 的当前上界：初始取 du^2
        nearest_pass2_other(q, b1, best_other_d2, best_d2);

        d_other = std::sqrt(best_other_d2);
        const double d2_sq = std::min(du2, best_other_d2);
        d2 = std::sqrt(d2_sq);
    }

    bool has_neighbor_within_blockD(const Point<double>& q) const {
        if (root_ < 0) return false;
        struct Item { int node; double lb2; };
        Item stack[64]; int sp = 0;
        auto mind2 = [&](int idx){ return nodes_[idx].box.mindist2(q); };

        stack[sp++] = { root_, mind2(root_) };
        while (sp > 0) {
            Item it = stack[--sp];
            const KDNode& nd = nodes_[it.node];
            const double thr2 = nd.min_block_D * nd.min_block_D;
            if (it.lb2 > thr2) continue;

            if (nd.is_leaf) {
                for (int i = nd.begin; i < nd.end; ++i) {
                    const auto& mp = pts_[i];
                    const double r2 = mp.block_D * mp.block_D;
                    if (r2 <= 0.0) continue;
                    const double d2 = squared_distance(mp.pos, q);
                    if (d2 < r2) return true;
                }
            } else {
                const int L = nd.left, R = nd.right;
                if (L >= 0) {
                    const double lbL = mind2(L);
                    if (lbL <= nodes_[L].min_block_D * nodes_[L].min_block_D) {
                        stack[sp++] = { L, lbL };
                    }
                }
                if (R >= 0) {
                    const double lbR = mind2(R);
                    if (lbR <= nodes_[R].min_block_D * nodes_[R].min_block_D) {
                        stack[sp++] = { R, lbR };
                    }
                }
            }
        }
        return false;
    }

    bool has_neighbor_within_radius(const Point<double>& q, double radius) const {
        if (radius <= 0.0 || root_ < 0) {
            return false;
        }
        const double radius2 = radius * radius;
        struct Item { int node; double lb2; };
        Item stack[64];
        int sp = 0;
        auto mind2 = [&](int idx){ return nodes_[idx].box.mindist2(q); };

        stack[sp++] = { root_, mind2(root_) };
        while (sp > 0) {
            Item it = stack[--sp];
            if (it.lb2 > radius2) continue;
            const KDNode& nd = nodes_[it.node];
            if (nd.is_leaf) {
                for (int i = nd.begin; i < nd.end; ++i) {
                    const auto& mp = pts_[i];
                    const double d2 = squared_distance(mp.pos, q);
                    if (d2 < radius2) {
                        return true;
                    }
                }
            } else {
                const int L = nd.left, R = nd.right;
                if (L >= 0) {
                    const double lbL = mind2(L);
                    if (lbL <= radius2) {
                        stack[sp++] = { L, lbL };
                    }
                }
                if (R >= 0) {
                    const double lbR = mind2(R);
                    if (lbR <= radius2) {
                        stack[sp++] = { R, lbR };
                    }
                }
            }
        }
        return false;
    }

private:
    int                   root_;
    std::vector<KDNode>   nodes_;
    std::vector<MovingPoint> pts_;

    // ----------------------------------------------------------------------------------
    // 辅助函数：两点间平方距离
    // 使用 Point<double>::size() 和 operator[]，与 AABB::mindist2 一致。
    // ----------------------------------------------------------------------------------
    static double squared_distance(const Point<double>& a,
        const Point<double>& b)
    {
        const int dim = static_cast<int>(a.size());
        double s = 0.0;
        for (int d = 0; d < dim; ++d) {
            const double diff = a[d] - b[d];
            s += diff * diff;
        }
        return s;
    }

    // ----------------------------------------------------------------------------------
    // 递归构建 KD 树：在 [begin, end) 范围内构建子树，返回根节点索引
    // ----------------------------------------------------------------------------------
    int build_range(int begin, int end)
    {
        const int n = end - begin;
        const int dim = static_cast<int>(pts_[begin].pos.size());

        // 1) 统计该范围内所有点的 AABB
        KDNode node;
        node.box.bmin = pts_[begin].pos;
        node.box.bmax = pts_[begin].pos;

        for (int i = begin + 1; i < end; ++i) {
            for (int d = 0; d < dim; ++d) {
                const double v = pts_[i].pos[d];
                if (v < node.box.bmin[d]) node.box.bmin[d] = v;
                if (v > node.box.bmax[d]) node.box.bmax[d] = v;
            }
        }

        // 2) 叶子条件：点数较少（这里取 8，可根据需要调整）
        constexpr int kLeafSize = 8;
        if (n <= kLeafSize) {
            node.is_leaf = true;
            node.begin = begin;
            node.end = end;
            node.left = -1;
            node.right = -1;

            // 叶子 mono_block：检查该范围内 block_id 是否全相同
            int blk0 = pts_[begin].block_id;
            bool mono = true;
            for (int i = begin + 1; i < end; ++i) {
                if (pts_[i].block_id != blk0) {
                    mono = false;
                    break;
                }
            }
            node.mono_block = mono ? blk0 : -1;
            node.min_block_D = std::numeric_limits<double>::infinity();
            for (int i = begin; i < end; ++i) {
                node.min_block_D = std::min(node.min_block_D, pts_[i].block_D);
            }

            nodes_.push_back(node);
            return static_cast<int>(nodes_.size()) - 1;
        }

        // 3) 选择划分轴：AABB 最大展开方向
        int    axis = 0;
        double best_extent = node.box.bmax[0] - node.box.bmin[0];
        for (int d = 1; d < dim; ++d) {
            const double extent = node.box.bmax[d] - node.box.bmin[d];
            if (extent > best_extent) {
                best_extent = extent;
                axis = d;
            }
        }

        // 4) 中位数划分（以 axis 维度上的坐标为 key）
        const int mid = begin + n / 2;
        std::nth_element(
            pts_.begin() + begin,
            pts_.begin() + mid,
            pts_.begin() + end,
            [axis](const MovingPoint& a, const MovingPoint& b) {
                return a.pos[axis] < b.pos[axis];
            });

        // 5) 占位并递归构建左右子树
        const int my_index = static_cast<int>(nodes_.size());
        nodes_.push_back(KDNode{}); // 占位
        const int left = build_range(begin, mid);
        const int right = build_range(mid, end);

        // 6) 填充内部节点信息
        KDNode& me = nodes_[my_index];
        me = node;
        me.is_leaf = false;
        me.left = left;
        me.right = right;
        me.axis = axis;
        me.split = pts_[mid].pos[axis];
        me.begin = begin;
        me.end = end;

        const int monoL = nodes_[left].mono_block;
        const int monoR = nodes_[right].mono_block;
        me.mono_block = (monoL >= 0 && monoL == monoR) ? monoL : -1;

        return my_index;
    }

    // ----------------------------------------------------------------------------------
    // Pass1: 标准 KD-tree 最近邻搜索（只关心“最近点”和“其 block_id”）
    // 输出：
    //   best_d2      : 最近点的距离平方
    //   best_block_id: 最近点所属分区 id
    // ----------------------------------------------------------------------------------
    void nearest_pass1(
        const Point<double>& q,
        double& best_d2,
        int& best_block_id
    ) const
    {
        best_d2 = std::numeric_limits<double>::infinity();
        best_block_id = -1;

        if (root_ < 0) return;

        struct StackItem {
            int    node;
            double lb2;  // 该节点 AABB 到 q 的最小距离平方
        };

        // 简单的手写栈，避免递归
        StackItem stack[64];
        int       sp = 0;

        stack[sp++] = { root_, nodes_[root_].box.mindist2(q) };

        while (sp > 0) {
            StackItem it = stack[--sp];
            if (it.lb2 >= best_d2) continue; // AABB 下界已经不可能更优，剪枝

            const KDNode& nd = nodes_[it.node];

            if (nd.is_leaf) {
                for (int i = nd.begin; i < nd.end; ++i) {
                    const MovingPoint& mp = pts_[i];
                    const double d2 = squared_distance(mp.pos, q);
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best_block_id = mp.block_id;
                    }
                }
            }
            else {
                const int L = nd.left;
                const int R = nd.right;
                if (L < 0 && R < 0) continue;

                const double inf = std::numeric_limits<double>::infinity();
                double lbL = inf, lbR = inf;
                if (L >= 0) lbL = nodes_[L].box.mindist2(q);
                if (R >= 0) lbR = nodes_[R].box.mindist2(q);

                // 压栈次序：先压“远”的，再压“近”的，这样下次循环先处理更近的，剪枝效果更好
                if (L >= 0 && lbL < best_d2) {
                    if (R >= 0 && lbR < best_d2) {
                        if (lbL < lbR) {
                            stack[sp++] = { R, lbR };
                            stack[sp++] = { L, lbL };
                        }
                        else {
                            stack[sp++] = { L, lbL };
                            stack[sp++] = { R, lbR };
                        }
                    }
                    else {
                        stack[sp++] = { L, lbL };
                    }
                }
                else if (R >= 0 && lbR < best_d2) {
                    stack[sp++] = { R, lbR };
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------
    // Pass2: 在已知最近分区 b1 的条件下，只在其它分区中寻找最近距离
    //
    // 输入：
    //   b1              : 最近点所属分区 id
    //   best_d2 (inout) : d2^2 的当前上界。初始为 du^2，搜索中若找到更小的
    //                     其它分区距离则进一步收紧（用于剪枝）。
    //
    // 输出：
    //   best_other_d2   : 来自“其它分区”的最近距离平方（若不存在，则为 inf）
    // ----------------------------------------------------------------------------------
    void nearest_pass2_other(
        const Point<double>& q,
        int    b1,
        double& best_other_d2,
        double& best_d2
    ) const
    {
        best_other_d2 = std::numeric_limits<double>::infinity();
        if (root_ < 0) return;

        struct StackItem {
            int    node;
            double lb2;
        };

        StackItem stack[64];
        int       sp = 0;

        stack[sp++] = { root_, nodes_[root_].box.mindist2(q) };

        while (sp > 0) {
            StackItem it = stack[--sp];
            const KDNode& nd = nodes_[it.node];

            // 剪枝 1: 该子树 AABB 下界已经不可能改善当前 best_d2
            if (it.lb2 >= best_d2) continue;

            // 剪枝 2: 整棵子树的点都属于最近分区 b1
            if (nd.mono_block == b1) continue;

            if (nd.is_leaf) {
                for (int i = nd.begin; i < nd.end; ++i) {
                    const MovingPoint& mp = pts_[i];
                    if (mp.block_id == b1) continue; // 排除最近分区本身

                    const double d2 = squared_distance(mp.pos, q);
                    if (d2 < best_other_d2) {
                        best_other_d2 = d2;
                        if (d2 < best_d2) best_d2 = d2; // d2^2 上界收紧
                    }
                }
            }
            else {
                const int L = nd.left;
                const int R = nd.right;
                if (L < 0 && R < 0) continue;

                const double inf = std::numeric_limits<double>::infinity();
                double lbL = inf, lbR = inf;
                if (L >= 0) lbL = nodes_[L].box.mindist2(q);
                if (R >= 0) lbR = nodes_[R].box.mindist2(q);

                if (L >= 0 && lbL < best_d2) {
                    if (R >= 0 && lbR < best_d2) {
                        if (lbL < lbR) {
                            stack[sp++] = { R, lbR };
                            stack[sp++] = { L, lbL };
                        }
                        else {
                            stack[sp++] = { L, lbL };
                            stack[sp++] = { R, lbR };
                        }
                    }
                    else {
                        stack[sp++] = { L, lbL };
                    }
                }
                else if (R >= 0 && lbR < best_d2) {
                    stack[sp++] = { R, lbR };
                }
            }
        }
    }
};
