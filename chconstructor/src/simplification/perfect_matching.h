#pragma once

#include <math.h>
#include <vector>

#include "prio_nodes.h"
#include "../nodes_and_edges.h"
#include "match_functions.h"


namespace ls {          
    
namespace mc {
    
template <class GraphT>
class PerfectMatching {
    const GraphT &graph;
    
    struct HandledNode {
        const chc::OSMNode n;
        //const std::list<PrioNodeHandle>::const_iterator pn_h; 
        const PrioNodeHandle pn_h; 
        
        HandledNode(const chc::OSMNode node, PrioNodeHandle pn_h):
            n(node), pn_h(pn_h) {}
    };
    
    struct TableCoordinate {
        int i;
        int j;
        
        bool isLast() {
            if (i==-1 || j==-1) {
                assert (i == -1 && j== -1);
                return true;
            } else {
                return false;
            }            
        }  
        TableCoordinate(): i(-1), j(-1) {}
        TableCoordinate(int i, int j): i(i), j(j) {}
        
    };
    
    struct TableEntry {
        double c = -1.0; //maximal cost
        //TableEntry *prev = nullptr;          
        TableCoordinate tc_prev = {-1, -1};
        TableEntry(): c(-1.0), tc_prev(-1, -1) {}
    };
    
    typedef std::vector<HandledNode> NodeVec;
    NodeVec p, q;
    
    std::vector<std::vector<TableEntry> > ca;

    NodeVec ChainToVec(const std::list<PrioNodeHandle> &pnhs) const {
        NodeVec node_vec;                    
        for (std::list<PrioNodeHandle>::const_iterator it = pnhs.begin(); it != pnhs.end(); it++) {
            PrioNodeHandle pnh = *it;
            const PrioNode2 &pn  = *pnh;
            const chc::OSMNode &node = graph.getNode(pn.node_id);
            HandledNode hn(node, pnh);
            node_vec.push_back(hn);                    
        }        
        return node_vec;
    }
   
    double c(int i, int j) { 

        assert(i >= 0 && j >= 0);
        assert(i < p.size() && j < q.size());

        if (ca[i][j].c > -1) {
            return ca[i][j].c;
        } else if (i == 0 and j == 0) {
            ca[i][j].c = geo::geoDist(p[0].n, q[0].n);
        } else if (i > 0 and j == 0) {                                        
            ca[i][j].c = std::max(c(i - 1, 0), geo::geoDist(p[i].n, q[0].n)); 
            ca[i][j].tc_prev = {i-1, 0}; //&ca[i-1, 0];            
        } else if (i == 0 and j > 0) {            
            ca[i][j].c = std::max(c(0, j - 1), geo::geoDist(p[0].n, q[j].n));
            ca[i][j].tc_prev = {0, j-1}; //&ca[0, j-1];
        } else if (i > 0 and j > 0) {
            double dist1 = c(i - 1, j);
            double dist2 = c(i - 1, j - 1);
            double dist3 = c(i, j - 1);
            if (dist1 < dist2) {
                if (dist1 < dist3) {
                    ca[i][j].c = dist1;
                    ca[i][j].tc_prev = {i-1, j}; //&ca[i-1, j];
                } else {
                    ca[i][j].c = dist3;
                    ca[i][j].tc_prev = {i, j-1}; //&ca[i, j-1];
                }
            } else {
                if (dist2 < dist3) {
                    ca[i][j].c = dist2;
                    ca[i][j].tc_prev = {i-1, j-1}; //&ca[i-1, j-1];
                } else {
                    ca[i][j].c = dist3;
                    ca[i][j].tc_prev = {i, j-1}; //&ca[i, j-1];
                }
            }
            
            double dist4 = geo::geoDist(p[i].n, q[j].n);
            if (dist4 > ca[i][j].c) {
                ca[i][j].c = dist4;
            }
            /*
            ca[i][j].c = std::max(
                     std::min(c(i - 1, j), std::min(c(i - 1, j - 1), c(i, j - 1))),
                     geo::geoDist(p[i], q[j])
                     );*/
        } else {
            ca[i][j].c = std::numeric_limits<double>::max();
        }

        return ca[i][j].c;
    }

    void setFollowers(uint i, uint j) {
        std::vector<std::list<PrioNodeHandle> > possibleFollowersVector_p;
        std::vector<std::list<PrioNodeHandle> > possibleFollowersVector_q;
        possibleFollowersVector_p.resize(i+1);
        possibleFollowersVector_q.resize(j+1);
                
        TableCoordinate cur(i,j); //&ca[i, j];
        
        while (!cur.isLast()) {
            
            //PrioNodeHandle pnh = q[prev.j].pn_h;
            possibleFollowersVector_p[cur.i].push_back(q[cur.j].pn_h);
            possibleFollowersVector_q[cur.j].push_back(p[cur.i].pn_h);            
            
            cur = ca[cur.i][cur.j].tc_prev;            
        }
        for (uint i = 0; i < possibleFollowersVector_p.size(); i++) {            
            getAndSetNearestFollower(p[i].pn_h, possibleFollowersVector_p[i], graph);
        }
        for (uint i = 0; i < possibleFollowersVector_q.size(); i++) {            
            getAndSetNearestFollower(q[i].pn_h, possibleFollowersVector_q[i], graph);
        }
    }
    
    double getChainLength(const std::list<PrioNodeHandle> &list) const {
        PrioNode2 pn_start = *list.front();
        PrioNode2 pn_end = *list.back();
        return geo::geoDist(graph.getNode(pn_start.node_id), graph.getNode(pn_end.node_id));
    }
    
    double getCombinedChainLength(const std::list<PrioNodeHandle> &list1, const std::list<PrioNodeHandle> &list2) const {        
        return getChainLength(list1) + getChainLength(list2);
    }

public:

    PerfectMatching(const GraphT &graph) : graph(graph) {            
    }
    
    void match(const std::list<PrioNodeHandle> &list1, const std::list<PrioNodeHandle> &list2) {
        //chains are already going in the same direction
        p = ChainToVec(list1);
        q = ChainToVec(list2);
        std::vector<TableEntry> subvector(q.size(), TableEntry());                        
        ca.assign(p.size(), subvector);                        
        //fill dynamic programming table
        uint i = p.size()-1;
        uint j = q.size()-1;
        double discrete_frechet_dist = c(i, j);
        if (discrete_frechet_dist < getCombinedChainLength(list1, list2)/5) {
            setFollowers(i, j);
        }
        
                
    }        
};
    
}

}