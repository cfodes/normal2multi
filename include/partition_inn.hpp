#pragma once
#include <queue>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "partition_bvh.hpp"
#include "distance.hpp"   // GridBTree<int, Point<double>>

struct INNResult
{
    int    nearest_block_id = -1;   // 最邻近分区 id
    double d1 = std::numeric_limits<double>::infinity();   // 最近分区真实距离（动边界距离）
    double d2 = std::numeric_limits<double>::infinity();   // 静边界距离 = min{ d_other, d_u }
    double d_other = std::numeric_limits<double>::infinity();   // 除最邻近分区外的最小真实距离
    double d_u = std::numeric_limits<double>::infinity();   // 到 unique_bndry_nodes 的最小距离
};

class PartitionINN
    // INN 算法，参考 Hjaltason & Samet: Distance Browsing in Spatial Databases
    // 本项目特化：用分区级 BVH + Best-first（按 AABB 到 q 的 LB^2）遍历，
    // 在叶子上做一次真实最近邻以得到 d1 / d_other，再与 d_u 取 min 得到 d2。
{
public:
    template<class Blocks>
    static INNResult query(
        const PartitionBVH& bvh,
        const Blocks& blks,
        GridBTree<int, Point<double>> const* unique_tree, // 可为 nullptr
        const Point<double>& q
    )
    {
        INNResult out;
        if (bvh.empty()) {
            throw std::runtime_error("[PartitionINN::query] Error: BVH is empty, cannot perform query.");
        }

        // 小根堆：按下界距离的平方(LB^2)从小到大
        struct Item { int node; double key; };   // key = mindist2(q, node.box)
        struct Cmp { bool operator()(const Item& a, const Item& b) const { return a.key > b.key; } };
        std::priority_queue<Item, std::vector<Item>, Cmp> pq;

        // push root（避免重复 root() 调用）
        const int root = bvh.root();
        pq.push({ root, bvh.node(root).box.mindist2(q) });

        bool queried_du = false; // 到 unique_bndry 的距离先不查

        // help func: 跳过“只含最近分区”的条目，返回“其它分区”的 top LB^2
        auto top_rest_lb2 = [&](int best_blk) -> double {
            while (!pq.empty()) {
                const auto& nd = bvh.node(pq.top().node);
                if (best_blk >= 0 && nd.mono_block_id == best_blk) { pq.pop(); continue; }
                return pq.top().key; // 这是“其它分区”的最小 LB^2
            }
            return std::numeric_limits<double>::infinity();
        };

        // 主循环：INN / Best-first
        while (!pq.empty()) {
            // —— 早停检查（需要先知道 top_rest）
            const double Lrest = top_rest_lb2(out.nearest_block_id);

            // 若已有最近分区，且 du 尚未查询 && Lrest 已不小于 d1^2，则“懒”查询 du（只查一次）
            if (out.nearest_block_id >= 0 && !queried_du && unique_tree) { // 指针作布尔判断：非空即真
                if (Lrest >= out.d1 * out.d1) {
                    int    dummy_key = 0;
                    double du_tmp = std::numeric_limits<double>::infinity();
                    unique_tree->search(du_tmp, dummy_key, q);
                    out.d_u = du_tmp;
                    queried_du = true;
                }
            }

            // 再次计算当前“停止阈值”
            const double t2 = std::min(out.d_other, out.d_u); // 线性距离
            const double T = std::max(out.d1, t2);           // 线性距离（注意：若 d_u<t2 可能小于 d1）
            const double T2 = T * T;                          // 与 LB^2 同量纲

            if (out.nearest_block_id >= 0 && Lrest >= T2) {
                // 满足：L_rest >= max(d1^2, min(other^2, du^2))，可以停止
                out.d2 = t2;
                return out;
            }

            // 取一个真实需要展开的条目（其子树包含“非最近分区”的可能）
            if (pq.empty()) break; // 可能被 top_rest_lb2 清空
            auto it = pq.top(); pq.pop();
            const auto& nd = bvh.node(it.node);

            // 如果这个条目仍是“仅含最近分区”的节点（可能在 top_rest_lb2 之后产生），跳过
            if (out.nearest_block_id >= 0 && nd.mono_block_id == out.nearest_block_id) {
                continue;
            }

            if (nd.leaf_block_id >= 0) {
                // —— 命中叶子：对该分区做一次“真实最近邻”（线性距离）
                const int blk_id = nd.leaf_block_id;
                int    key_tmp = 0;
                double d_tmp = std::numeric_limits<double>::infinity();
                blks[blk_id].block_tree.search(d_tmp, key_tmp, q);

                // 维护滚动最小对：(d1, nearest_block_id) 与 d_other（来自“其它分区”）
                if (d_tmp < out.d1) {
                    out.d_other = std::min(out.d_other, out.d1); // 旧的 d1 成为 d_other 候选
                    out.d1 = d_tmp;
                    out.nearest_block_id = blk_id;
                }
                else if (blk_id != out.nearest_block_id && d_tmp < out.d_other) {
                    out.d_other = d_tmp;
                }

                // 有了 d_other 候选后，可立即查询一次 du（若还没查过）
                if (!queried_du && unique_tree && std::isfinite(out.d_other)) {
                    int    dummy_key = 0;
                    double du_tmp = std::numeric_limits<double>::infinity();
                    unique_tree->search(du_tmp, dummy_key, q);
                    out.d_u = du_tmp;
                    queried_du = true;
                }
            }
            else {
                // —— 内部节点：继续按 LB^2 展开
                const int L = nd.left, R = nd.right;
                if (L >= 0) pq.push({ L, bvh.node(L).box.mindist2(q) });
                if (R >= 0) pq.push({ R, bvh.node(R).box.mindist2(q) });
            }
        }

        // 堆耗尽：返回当前最优（通常不会走到这里）
        out.d2 = std::min(out.d_other, out.d_u);
        return out;
    }
};
