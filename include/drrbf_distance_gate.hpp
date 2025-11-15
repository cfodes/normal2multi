#pragma once

#include "global_kdtree.hpp"

namespace drrbf {

inline bool has_neighbor_within_Dm(const GlobalKDTree& tree,
                                   const Point<double>& q,
                                   double Dm)
{
    return tree.has_neighbor_within_radius(q, Dm);
}

} // namespace drrbf
