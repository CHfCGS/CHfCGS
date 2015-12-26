#pragma once

#include "chgraph.h"
#include "nodes_and_edges.h"
#include "geoFunctions.h"

namespace chc {
    namespace ew {
        static double calcWeight(double dist, const OSMNode& src, const OSMNode& tgt) {
            double geo_dist = geo::geoDist(src, tgt);
            return geo_dist * pow(geo_dist/dist, 3);            
        }
    }
}