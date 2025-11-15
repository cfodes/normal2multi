#pragma once

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include "geometry.hpp"

// Simple AABB for 3D points
struct GK_AABB {
    Point<double> bmin{std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity()};
    Point<double> bmax{-std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()};

    inline void expand(const Point<double>& p) {
        for (int d = 0; d < 3; ++d) {
            if (p[d] < bmin[d]) bmin[d] = p[d];
            if (p[d] > bmax[d]) bmax[d] = p[d];
        }
    }

    inline double mindist2(const Point<double>& q) const {
        double acc = 0.0;
        for (int d = 0; d < 3; ++d) {
            const double v = q[d];
            const double lo = bmin[d];
            const double hi = bmax[d];
            double diff = 0.0;
            if (v < lo) diff = lo - v;
            else if (v > hi) diff = v - hi;
            acc += diff * diff;
        }
        return acc;
    }
};

struct MovingPoint {
    Point<double> pos;   // 坐标
    int           block_id = -1; // 所属 block
    int           node_id  = -1; // 原始 Node.id（可选）
};

struct KDNode {
    GK_AABB box;
    int  left  = -1;
    int  right = -1;

    int  begin = 0;    // 叶子区间 [begin, end)
    int  end   = 0;
    int  axis  = 0;
    double split = 0.0;

    int  mono_block = -1; // 子树是否全为同一 block；否则为 -1
    bool is_leaf = false;
};

class GlobalKDTree {
public:
    GlobalKDTree() = default;

    template<class Blocks>
    void build_from_blocks(const Blocks& blks) {
        pts_.clear();
        // 预估容量
        size_t estimate = 0;
        for (const auto& b : blks) estimate += b.internal_nodes.size();
        pts_.reserve(estimate);

        for (const auto& b : blks) {
            if (b.internal_nodes.empty()) continue;
            for (const auto& nd : b.internal_nodes) {
                MovingPoint mp;
                mp.pos = nd.point;
                mp.block_id = b.block_id;
                mp.node_id = nd.id;
                pts_.push_back(std::move(mp));
            }
        }

        nodes_.clear();
        if (pts_.empty()) { root_ = -1; return; }
        nodes_.reserve(pts_.size() * 2);
        root_ = build_range(0, static_cast<int>(pts_.size()));
    }

    bool empty() const noexcept { return root_ < 0; }

    void query_point(
        const Point<double>& q,
        double du2,
        double& d1, int& moving_block_id,
        double& d2, double& d_other, double& d_u
    ) const {
        double best1_d2; int b1;
        nearest_pass1(q, best1_d2, b1);

        if (b1 < 0 || !std::isfinite(best1_d2)) {
            d1 = d2 = d_other = d_u = std::numeric_limits<double>::infinity();
            moving_block_id = -1;
            return;
        }

        moving_block_id = b1;
        d1 = std::sqrt(best1_d2);
        d_u = std::sqrt(du2);

        double best_other_d2;
        double best_d2 = du2; // 初始由 du2 给出上界
        nearest_pass2_other(q, b1, best_other_d2, best_d2);

        d_other = std::sqrt(best_other_d2);
        const double d2_sq = std::min(du2, best_other_d2);
        d2 = std::sqrt(d2_sq);
    }

private:
    int root_ = -1;
    std::vector<KDNode>      nodes_;
    std::vector<MovingPoint> pts_;

    static inline double sqdist_pt_pt(const Point<double>& a, const Point<double>& b) {
        const double dx = a.x - b.x;
        const double dy = a.y - b.y;
        const double dz = a.z - b.z;
        return dx*dx + dy*dy + dz*dz;
    }

    int build_range(int begin, int end) {
        const int n = end - begin;

        GK_AABB box;
        for (int i = begin; i < end; ++i) box.expand(pts_[i].pos);

        constexpr int LEAF_SIZE = 8;
        if (n <= LEAF_SIZE) {
            KDNode me;
            me.box = box;
            me.is_leaf = true;
            me.begin = begin;
            me.end = end;

            int blk0 = pts_[begin].block_id;
            bool mono = true;
            for (int i = begin + 1; i < end; ++i) {
                if (pts_[i].block_id != blk0) { mono = false; break; }
            }
            me.mono_block = mono ? blk0 : -1;

            nodes_.push_back(me);
            return static_cast<int>(nodes_.size()) - 1;
        }

        // choose axis by max extent
        int axis = 0; double best_extent = box.bmax[0] - box.bmin[0];
        for (int d = 1; d < 3; ++d) {
            double e = box.bmax[d] - box.bmin[d];
            if (e > best_extent) { best_extent = e; axis = d; }
        }

        const int mid = begin + n / 2;
        std::nth_element(
            pts_.begin() + begin,
            pts_.begin() + mid,
            pts_.begin() + end,
            [axis](const MovingPoint& a, const MovingPoint& b){ return a.pos[axis] < b.pos[axis]; }
        );

        const int my = static_cast<int>(nodes_.size());
        nodes_.push_back(KDNode{}); // placeholder

        const int L = build_range(begin, mid);
        const int R = build_range(mid, end);

        KDNode me;
        me.is_leaf = false;
        me.left = L; me.right = R;
        me.axis = axis; me.split = pts_[mid].pos[axis];
        me.box = box;

        const int monoL = nodes_[L].mono_block;
        const int monoR = nodes_[R].mono_block;
        me.mono_block = (monoL >= 0 && monoL == monoR) ? monoL : -1;

        nodes_[my] = me;
        return my;
    }

    void nearest_pass1(const Point<double>& q, double& best_d2, int& best_block_id) const {
        best_d2 = std::numeric_limits<double>::infinity();
        best_block_id = -1;
        if (root_ < 0) return;

        struct Item { int node; double lb2; };
        Item stack[64]; int sp = 0;
        auto mind2 = [&](int idx){ return nodes_[idx].box.mindist2(q); };

        stack[sp++] = { root_, mind2(root_) };
        while (sp > 0) {
            Item it = stack[--sp];
            const KDNode& nd = nodes_[it.node];
            if (it.lb2 >= best_d2) continue;

            if (nd.is_leaf) {
                for (int i = nd.begin; i < nd.end; ++i) {
                    const auto& mp = pts_[i];
                    const double d2 = sqdist_pt_pt(mp.pos, q);
                    if (d2 < best_d2) { best_d2 = d2; best_block_id = mp.block_id; }
                }
            } else {
                const int L = nd.left, R = nd.right;
                if (L < 0 && R < 0) continue;
                const double lbL = (L >= 0) ? mind2(L) : std::numeric_limits<double>::infinity();
                const double lbR = (R >= 0) ? mind2(R) : std::numeric_limits<double>::infinity();

                if (L >= 0 && lbL < best_d2) {
                    if (R >= 0 && lbR < best_d2) {
                        if (lbL < lbR) { stack[sp++] = { R, lbR }; stack[sp++] = { L, lbL }; }
                        else            { stack[sp++] = { L, lbL }; stack[sp++] = { R, lbR }; }
                    } else {
                        stack[sp++] = { L, lbL };
                    }
                } else if (R >= 0 && lbR < best_d2) {
                    stack[sp++] = { R, lbR };
                }
            }
        }
    }

    void nearest_pass2_other(const Point<double>& q, int b1, double& best_other_d2, double& best_d2) const {
        best_other_d2 = std::numeric_limits<double>::infinity();
        if (root_ < 0) return;

        struct Item { int node; double lb2; };
        Item stack[64]; int sp = 0;
        auto mind2 = [&](int idx){ return nodes_[idx].box.mindist2(q); };

        stack[sp++] = { root_, mind2(root_) };
        while (sp > 0) {
            Item it = stack[--sp];
            const KDNode& nd = nodes_[it.node];
            if (it.lb2 >= best_d2) continue;         // 剪枝 1：下界
            if (nd.mono_block == b1) continue;       // 剪枝 2：整棵子树为最近分区

            if (nd.is_leaf) {
                for (int i = nd.begin; i < nd.end; ++i) {
                    const auto& mp = pts_[i];
                    if (mp.block_id == b1) continue;
                    const double d2 = sqdist_pt_pt(mp.pos, q);
                    if (d2 < best_other_d2) {
                        best_other_d2 = d2;
                        if (d2 < best_d2) best_d2 = d2;
                    }
                }
            } else {
                const int L = nd.left, R = nd.right;
                if (L < 0 && R < 0) continue;
                const double lbL = (L >= 0) ? mind2(L) : std::numeric_limits<double>::infinity();
                const double lbR = (R >= 0) ? mind2(R) : std::numeric_limits<double>::infinity();

                if (L >= 0 && lbL < best_d2) {
                    if (R >= 0 && lbR < best_d2) {
                        if (lbL < lbR) { stack[sp++] = { R, lbR }; stack[sp++] = { L, lbL }; }
                        else            { stack[sp++] = { L, lbL }; stack[sp++] = { R, lbR }; }
                    } else {
                        stack[sp++] = { L, lbL };
                    }
                } else if (R >= 0 && lbR < best_d2) {
                    stack[sp++] = { R, lbR };
                }
            }
        }
    }
};

