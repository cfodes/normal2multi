#pragma once
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include "geometry.hpp"  // 支持 Point

struct AABB
{
    // 轴对齐包围盒：用两个点描述 bmin / bmax
    Point<double> bmin{}, bmax{};

    // 合并两个盒（父盒 = 子盒并集）
    static AABB merge(const AABB& a, const AABB& b)
    {
        AABB r;
        for (int d = 0; d != Point<double>().size(); ++d) {
            r.bmin[d] = std::min(a.bmin[d], b.bmin[d]);
            r.bmax[d] = std::max(a.bmax[d], b.bmax[d]);
        }
        return r;
    }

    // 点到 AABB 的最小距离平方
    double mindist2(const Point<double>& p) const noexcept
    {
        double s = 0.0;
        for (int d = 0; d != Point<double>().size(); ++d) {
            double di = 0.0; // 点在盒内 → 该轴距离 0
            if (p[d] < bmin[d]) di = bmin[d] - p[d];   // 左侧
            else if (p[d] > bmax[d]) di = p[d] - bmax[d];   // 右侧
            s += di * di;
        }
        return s;
    }
};

struct PartitionLeaf
{
    int block_id = -1;     // 分区 id
    AABB box{};            // 分区包围盒
    Point<double> centroid{}; // 盒中心，用于“最大展宽维 + 中位切分”
};

class PartitionBVH
{
public:
    struct Node {
        AABB box{};
        int left = -1;      // 内部节点：左子
        int right = -1;      // 内部节点：右子
        int leaf_block_id = -1; // 叶子：分区 id（非叶为 -1）

        // 查询剪枝辅助
        int leaf_count = 0;   // 子树叶子数
        int mono_block_id = -1; // 若子树只含一个分区，则该分区 id；否则 -1
    };

    // 用叶子数组构建 BVH（扁平化节点数组，cache 友好）
    void build(std::vector<PartitionLeaf> leaves)
    {
        leaves_ = std::move(leaves);
        nodes_.clear();
        index_.resize(leaves_.size());
        for (int i = 0; i < static_cast<int>(leaves_.size()); ++i) index_[i] = i;
        root_ = leaves_.empty() ? -1 : build_range(0, static_cast<int>(index_.size()));
    }

    // 接口
    bool empty() const noexcept { return root_ < 0; }
    int  root()  const noexcept { return root_; }
    const Node& node(int idx) const { return nodes_[idx]; }
    std::size_t leaf_count() const noexcept { return leaves_.size(); };

private:
    int root_ = -1;
    std::vector<PartitionLeaf> leaves_;
    std::vector<int> index_;
    std::vector<Node> nodes_;

    // 在 [begin, end) 的叶子索引范围上建子树，返回该子树根在 nodes_ 的索引
    int build_range(int begin, int end)
    {
        const int n = end - begin;

        // 1) 合并包围盒，覆盖本范围全部叶子
        AABB box = leaves_[index_[begin]].box;
        for (int i = begin + 1; i < end; ++i) box = AABB::merge(box, leaves_[index_[i]].box);

        // 2) 单叶：构建叶子节点
        if (n == 1) {
            Node me;
            me.box = box;
            me.leaf_block_id = leaves_[index_[begin]].block_id;
            me.leaf_count = 1;
            me.mono_block_id = me.leaf_block_id;
            nodes_.push_back(me);
            return static_cast<int>(nodes_.size()) - 1;
        }

        // 3) 选择切分轴 = 最大展宽维（按 centroid）
        Point<double> cmin = leaves_[index_[begin]].centroid;
        Point<double> cmax = cmin;
        for (int i = begin + 1; i < end; ++i) {
            const auto& c = leaves_[index_[i]].centroid;
            for (int d = 0; d != Point<double>().size(); ++d) {
                cmin[d] = std::min(cmin[d], c[d]);
                cmax[d] = std::max(cmax[d], c[d]);
            }
        }
        int axis = 0;
        double best_extent = cmax[0] - cmin[0];
        for (int d = 1; d != Point<double>().size(); ++d) {
            double e = cmax[d] - cmin[d];
            if (e > best_extent) { best_extent = e; axis = d; }
        }

        // 4) 中位切分，平衡两侧规模 O(n)
        const int mid = begin + n / 2;
        std::nth_element(index_.begin() + begin, index_.begin() + mid, index_.begin() + end,
            [&](int ia, int ib) { return leaves_[ia].centroid[axis] < leaves_[ib].centroid[axis]; });

        // 5) 占位并递归建左右子
        const int my = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{}); // 占位
        const int L = build_range(begin, mid);
        const int R = build_range(mid, end);

        // 6) 填充内部节点
        Node me;
        me.left = L;
        me.right = R;
        me.box = AABB::merge(nodes_[L].box, nodes_[R].box);
        me.leaf_count = nodes_[L].leaf_count + nodes_[R].leaf_count;
        me.mono_block_id = (nodes_[L].mono_block_id == nodes_[R].mono_block_id)
            ? nodes_[L].mono_block_id : -1;

        nodes_[my] = me;
        return my;
    }
};
