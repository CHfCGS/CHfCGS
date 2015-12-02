#pragma once

#include <limits>
#include <algorithm>

#include "../self_intersection_checker.h"
#include "cdt_matching.h"
#include "../nodes_and_edges.h"
#include "prio_nodes.h"
#include "zip_order.h"

#undef matchChainPairNodesNDEBUG


//using namespace ls;

namespace ls {
        
    namespace mc {
        
        template <class GraphT>
        void match(const GraphT &graph, std::list<PrioNodeHandle> &list1, std::list<PrioNodeHandle> &list2) {
        debug_assert(list1.size()>=1 && list2.size()>=1); 
        
        
        Chain chain1;
        for (PrioNodeHandle pnh: list1) {
            //const PrioNodeHandle pnh = *it;
            const PrioNode2 &pn = *pnh;
            chain1.push_back(pn.node_id);
        }
        Chain chain2;
        for (PrioNodeHandle pnh: list2) {
            //const PrioNodeHandle pnh = *it;
            const PrioNode2 &pn = *pnh;
            chain2.push_back(pn.node_id);
        }
        //additionally to make sure cdt_matching constraints do not not intersect
        
        chain1.push_front(chain2.front());
        chain2.push_back(chain1.back());                
        
        
        SelfIntersectionChecker<GraphT> selfIC(graph);
        bool isSelfIntersecting = selfIC.isSelfIntersecting(chain1, chain2);
        //bool isSelfIntersecting = selfIC.isSelfIntersecting(chain);
        
        ZipOrder<GraphT> zipOrder(graph);
        bool useZipOrder = isSelfIntersecting;
        if (useZipOrder) {
            zipOrder.match(list1, list2);            
        } else {            
            //zipOrder.match(list1, list2);  
            CDTMatching2<GraphT> ctd_matcher(graph);
            ctd_matcher.match(list1, list2);
        }
        return;
    }    
        
        
//namespace MCPN {
//static class
        /*
template <class GraphT>
class matchChainPairNodes2 {
            
public:
    matchChainPairNodes2(const GraphT &base_graph, std::list<PrioNodeHandle> &list1, std::list<PrioNodeHandle> &list1)
        : base_graph(base_graph), list1(list1), list2(list2) {}
    
    void match() {
        debug_assert(list1.size()>=1 && list2.size()>=1); 
        //SelfIntersectionChecker<GraphT> selfIC(base_graph);
        Chain chainTo;
        for (auto it: list1) {
            chainTo.push_back(it->node_id);
        }
        Chain chainFrom;
        for (auto it: list2) {
            chainFrom.push_back(it->node_id);
        }
        //additionally to make sure cdt_matching constraints do not not intersect 
        chainTo.push_back(chainFrom.front());
        chainFrom.push_back(chainTo.front());
        bool isSelfIntersecting = true;//= selfIC.isSelfIntersecting(chainTo, chainFrom);
        
        bool useZipOrder = isSelfIntersecting;
        if (useZipOrder) {
            zipOrder();            
        } else {            
            zipOrder();  
            //CDTMatching<GraphT> ctd_matcher(base_graph);
            //ctd_matcher.match(toList, fromList);
        }
        return;
    }            
};*/

    }
}
