#pragma once

#include <limits>
#include <algorithm>

//#include "../self_intersection_checker.h"
//#include "cdt_matching.h"
#include "../nodes_and_edges.h"
#include "prio_nodes.h"
#include "zip_order.h"
#include "../s_options.h"

#undef matchChainPairNodesNDEBUG


//using namespace ls;

namespace ls {
        
    namespace mc {
        
        //chains are already going in the same direction
        template <class GraphT>
        void match(const GraphT &graph, std::list<PrioNodeHandle> &list1, std::list<PrioNodeHandle> &list2,
                PairMatchType pairMatch_type) {
            debug_assert(list1.size()>=1 && list2.size()>=1); 
        
            switch (pairMatch_type) {
                case PairMatchType::NONE:
                    break;
                case PairMatchType::ZO: {
                    ZipOrder<GraphT> zipOrder(graph);
                    zipOrder.match(list1, list2);  
                    break;
                }
                case PairMatchType::CD: {
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

                    //SelfIntersectionChecker<GraphT> selfIC(graph);
                    bool isSelfIntersecting = true;// selfIC.isSelfIntersecting(chain1, chain2);
                    ZipOrder<GraphT> zipOrder(graph);
                    if (isSelfIntersecting) {
                        zipOrder.match(list1, list2);            
                    } else {            
                        //CDTMatching2<GraphT> ctd_matcher(graph);
                        //ctd_matcher.match(list1, list2);
                    }                                

                    break;
                }
            }              
            return;
        }            

    }
}
