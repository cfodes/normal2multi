#pragma once
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include "geometry.hpp"  // ֧�� Point

struct AABB
{
    // ������Χ�У������������� bmin / bmax
    Point<double> bmin{}, bmax{};

    // �ϲ������У����� = �Ӻв�����
    static AABB merge(const AABB& a, const AABB& b)
    {
        AABB r;
        for (int d = 0; d != Point<double>().size(); ++d) {
            r.bmin[d] = std::min(a.bmin[d], b.bmin[d]);
            r.bmax[d] = std::max(a.bmax[d], b.bmax[d]);
        }
        return r;
    }

    // �㵽 AABB ����С����ƽ��
    double mindist2(const Point<double>& p) const noexcept
    {
        double s = 0.0;
        for (int d = 0; d != Point<double>().size(); ++d) {
            double di = 0.0; // ���ں��� �� ������� 0
            if (p[d] < bmin[d]) di = bmin[d] - p[d];   // ���
            else if (p[d] > bmax[d]) di = p[d] - bmax[d];   // �Ҳ�
            s += di * di;
        }
        return s;
    }
};

struct PartitionLeaf
{
    int block_id = -1;     // ���� id
    AABB box{};            // ������Χ��
    Point<double> centroid{}; // �����ģ����ڡ����չ��ά + ��λ�з֡�
};

class PartitionBVH
{
public:
    struct Node {
        AABB box{};
        int left = -1;      // �ڲ��ڵ㣺����
        int right = -1;      // �ڲ��ڵ㣺����
        int leaf_block_id = -1; // Ҷ�ӣ����� id����ҶΪ -1��

        // ��ѯ��֦����
        int leaf_count = 0;   // ����Ҷ����
        int mono_block_id = -1; // ������ֻ��һ����������÷��� id������ -1
    };

    // ��Ҷ�����鹹�� BVH����ƽ���ڵ����飬cache �Ѻã�
    void build(std::vector<PartitionLeaf> leaves)
    {
        leaves_ = std::move(leaves);
        nodes_.clear();
        index_.resize(leaves_.size());
        for (int i = 0; i < static_cast<int>(leaves_.size()); ++i) index_[i] = i;
        root_ = leaves_.empty() ? -1 : build_range(0, static_cast<int>(index_.size()));
    }

    // �ӿ�
    bool empty() const noexcept { return root_ < 0; }
    int  root()  const noexcept { return root_; }
    const Node& node(int idx) const { return nodes_[idx]; }
    std::size_t leaf_count() const noexcept { return leaves_.size(); };

private:
    int root_ = -1;
    std::vector<PartitionLeaf> leaves_;
    std::vector<int> index_;
    std::vector<Node> nodes_;

    // �� [begin, end) ��Ҷ��������Χ�Ͻ����������ظ��������� nodes_ ������
    int build_range(int begin, int end)
    {
        const int n = end - begin;

        // 1) �ϲ���Χ�У����Ǳ���Χȫ��Ҷ��
        AABB box = leaves_[index_[begin]].box;
        for (int i = begin + 1; i < end; ++i) box = AABB::merge(box, leaves_[index_[i]].box);

        // 2) ��Ҷ������Ҷ�ӽڵ�
        if (n == 1) {
            Node me;
            me.box = box;
            me.leaf_block_id = leaves_[index_[begin]].block_id;
            me.leaf_count = 1;
            me.mono_block_id = me.leaf_block_id;
            nodes_.push_back(me);
            return static_cast<int>(nodes_.size()) - 1;
        }

        // 3) ѡ���з��� = ���չ��ά���� centroid��
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

        // 4) ��λ�з֣�ƽ�������ģ O(n)
        const int mid = begin + n / 2;
        std::nth_element(index_.begin() + begin, index_.begin() + mid, index_.begin() + end,
            [&](int ia, int ib) { return leaves_[ia].centroid[axis] < leaves_[ib].centroid[axis]; });

        // 5) ռλ���ݹ齨������
        const int my = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{}); // ռλ
        const int L = build_range(begin, mid);
        const int R = build_range(mid, end);

        // 6) ����ڲ��ڵ�
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
