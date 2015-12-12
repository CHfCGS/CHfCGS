#pragma once

namespace ls { 
    
namespace mc {

//functions used by different matching algorithms
    
void _setFollower (PrioNodeHandle node_h,
              PrioNodeHandle follower_h) {
    PrioNode2 &node = *node_h;
    PrioNode2 &follower = *follower_h;
    node.followerValid = true;
    node.follower_h = follower_h;            
    
    std::list<PrioNodeHandle> &guidelist = follower.guides;
    guidelist.push_back(node_h);
}

template <class GraphT>
double geoDist(const PrioNodeHandle pnh1, const PrioNodeHandle pnh2, const GraphT &graph) {
    const PrioNode2 &pn1 = *pnh1;
    const PrioNode2 &pn2 = *pnh2;
    return geo::geoDist(graph.getNode(pn1.node_id), graph.getNode(pn2.node_id));
}



template <class GraphT>
void getAndSetNearestFollower(PrioNodeHandle pnh,
            std::list<PrioNodeHandle> possibleFollowers,
            const GraphT &graph) {
    auto nearestFollowerIt = possibleFollowers.end();
    double minDist = std::numeric_limits<double>::max();
    for (auto it = possibleFollowers.begin(); it != possibleFollowers.end(); it++) {
        double dist = geoDist(pnh, *it, graph);
        if (dist < minDist) {
            minDist = dist;
            nearestFollowerIt = it;
        }            
    }
    if(nearestFollowerIt != possibleFollowers.end()) {
        _setFollower(pnh, *nearestFollowerIt);
    }

}
    
}

}