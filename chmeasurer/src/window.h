#pragma once

#include <limits>

#include "nodes_and_edges.h"
#include "chgraph.h"

namespace chm {

struct Window {
    double MINLONGITUDE = std::numeric_limits<double>::max();
    double MAXLONGITUDE = 0;
    double MINLATITUDE = std::numeric_limits<double>::max();
    double MAXLATITUDE = 0;

    bool isIn(const CHGraph<CHNode, CHEdge> &graph, NodeID node_id) {
        CHNode node = graph.getNode(node_id);

        if (node.lat >= MINLATITUDE
                && node.lat <= MAXLATITUDE
                && node.lon >= MINLONGITUDE
                && node.lon <= MAXLONGITUDE) {
            return true;
        } else {
            return false;
        }
    }
    
    Window(const CHGraph<CHNode, CHEdge> &graph, const Chain &chain) {
        for (NodeID node_id : chain) {
            if (graph.getLon(node_id) > MAXLONGITUDE) {
                MAXLONGITUDE = graph.getLon(node_id);
            }
            if (graph.getLon(node_id) < MINLONGITUDE) {
                MINLONGITUDE = graph.getLon(node_id);
            }
            if (graph.getLat(node_id) > MAXLATITUDE) {
                MAXLATITUDE = graph.getLat(node_id);
            }
            if (graph.getLat(node_id) < MINLATITUDE) {
                MINLATITUDE = graph.getLat(node_id);
            }
        }        
    }    

};

}