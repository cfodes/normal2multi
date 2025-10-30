#pragma once
#include <queue>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "partition_bvh.hpp"
#include "distance.hpp"   // GridBTree<int, Point<double>>

struct INNResult
{
    int    nearest_block_id = -1;   // ���ڽ����� id
    double d1 = std::numeric_limits<double>::infinity();   // ���������ʵ���루���߽���룩
    double d2 = std::numeric_limits<double>::infinity();   // ���߽���� = min{ d_other, d_u }
    double d_other = std::numeric_limits<double>::infinity();   // �����ڽ����������С��ʵ����
    double d_u = std::numeric_limits<double>::infinity();   // �� unique_bndry_nodes ����С����
};

class PartitionINN
    // INN �㷨���ο� Hjaltason & Samet: Distance Browsing in Spatial Databases
    // ����Ŀ�ػ����÷����� BVH + Best-first���� AABB �� q �� LB^2��������
    // ��Ҷ������һ����ʵ������Եõ� d1 / d_other������ d_u ȡ min �õ� d2��
{
public:
    template<class Blocks>
    static INNResult query(
        const PartitionBVH& bvh,
        const Blocks& blks,
        GridBTree<int, Point<double>> const* unique_tree, // ��Ϊ nullptr
        const Point<double>& q
    )
    {
        INNResult out;
        if (bvh.empty()) {
            throw std::runtime_error("[PartitionINN::query] Error: BVH is empty, cannot perform query.");
        }

        // С���ѣ����½�����ƽ��(LB^2)��С����
        struct Item { int node; double key; };   // key = mindist2(q, node.box)
        struct Cmp { bool operator()(const Item& a, const Item& b) const { return a.key > b.key; } };
        std::priority_queue<Item, std::vector<Item>, Cmp> pq;

        // push root�������ظ� root() ���ã�
        const int root = bvh.root();
        pq.push({ root, bvh.node(root).box.mindist2(q) });

        bool queried_du = false; // �� unique_bndry �ľ����Ȳ���

        // help func: ������ֻ���������������Ŀ�����ء������������� top LB^2
        auto top_rest_lb2 = [&](int best_blk) -> double {
            while (!pq.empty()) {
                const auto& nd = bvh.node(pq.top().node);
                if (best_blk >= 0 && nd.mono_block_id == best_blk) { pq.pop(); continue; }
                return pq.top().key; // ���ǡ���������������С LB^2
            }
            return std::numeric_limits<double>::infinity();
        };

        // ��ѭ����INN / Best-first
        while (!pq.empty()) {
            // ���� ��ͣ��飨��Ҫ��֪�� top_rest��
            const double Lrest = top_rest_lb2(out.nearest_block_id);

            // ����������������� du ��δ��ѯ && Lrest �Ѳ�С�� d1^2����������ѯ du��ֻ��һ�Σ�
            if (out.nearest_block_id >= 0 && !queried_du && unique_tree) { // ָ���������жϣ��ǿռ���
                if (Lrest >= out.d1 * out.d1) {
                    int    dummy_key = 0;
                    double du_tmp = std::numeric_limits<double>::infinity();
                    unique_tree->search(du_tmp, dummy_key, q);
                    out.d_u = du_tmp;
                    queried_du = true;
                }
            }

            // �ٴμ��㵱ǰ��ֹͣ��ֵ��
            const double t2 = std::min(out.d_other, out.d_u); // ���Ծ���
            const double T = std::max(out.d1, t2);           // ���Ծ��루ע�⣺�� d_u<t2 ����С�� d1��
            const double T2 = T * T;                          // �� LB^2 ͬ����

            if (out.nearest_block_id >= 0 && Lrest >= T2) {
                // ���㣺L_rest >= max(d1^2, min(other^2, du^2))������ֹͣ
                out.d2 = t2;
                return out;
            }

            // ȡһ����ʵ��Ҫչ������Ŀ��������������������������Ŀ��ܣ�
            if (pq.empty()) break; // ���ܱ� top_rest_lb2 ���
            auto it = pq.top(); pq.pop();
            const auto& nd = bvh.node(it.node);

            // ��������Ŀ���ǡ���������������Ľڵ㣨������ top_rest_lb2 ֮�������������
            if (out.nearest_block_id >= 0 && nd.mono_block_id == out.nearest_block_id) {
                continue;
            }

            if (nd.leaf_block_id >= 0) {
                // ���� ����Ҷ�ӣ��Ը÷�����һ�Ρ���ʵ����ڡ������Ծ��룩
                const int blk_id = nd.leaf_block_id;
                int    key_tmp = 0;
                double d_tmp = std::numeric_limits<double>::infinity();
                blks[blk_id].block_tree.search(d_tmp, key_tmp, q);

                // ά��������С�ԣ�(d1, nearest_block_id) �� d_other�����ԡ�������������
                if (d_tmp < out.d1) {
                    out.d_other = std::min(out.d_other, out.d1); // �ɵ� d1 ��Ϊ d_other ��ѡ
                    out.d1 = d_tmp;
                    out.nearest_block_id = blk_id;
                }
                else if (blk_id != out.nearest_block_id && d_tmp < out.d_other) {
                    out.d_other = d_tmp;
                }

                // ���� d_other ��ѡ�󣬿�������ѯһ�� du������û�����
                if (!queried_du && unique_tree && std::isfinite(out.d_other)) {
                    int    dummy_key = 0;
                    double du_tmp = std::numeric_limits<double>::infinity();
                    unique_tree->search(du_tmp, dummy_key, q);
                    out.d_u = du_tmp;
                    queried_du = true;
                }
            }
            else {
                // ���� �ڲ��ڵ㣺������ LB^2 չ��
                const int L = nd.left, R = nd.right;
                if (L >= 0) pq.push({ L, bvh.node(L).box.mindist2(q) });
                if (R >= 0) pq.push({ R, bvh.node(R).box.mindist2(q) });
            }
        }

        // �Ѻľ������ص�ǰ���ţ�ͨ�������ߵ����
        out.d2 = std::min(out.d_other, out.d_u);
        return out;
    }
};
