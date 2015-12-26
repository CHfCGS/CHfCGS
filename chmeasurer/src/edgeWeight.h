#pragma once

//#include "chgraph.h"
#include "nodes_and_edges.h"
#include "geoFunctions.h"

namespace chm {
    namespace ew {
        static double calcWeight(double dist, const CHNode& src, const CHNode& tgt) {
            double geo_dist = geo::geoDist(src, tgt);
            return geo_dist * geo_dist/dist;            
        }
    }
}