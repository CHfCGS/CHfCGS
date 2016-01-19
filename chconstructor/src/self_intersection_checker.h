#pragma once

// Computing intersection points among curves using the sweep line.
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <list>

#include "chains.h"
#include "chgraph.h"

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
//typedef CGAL::Simple_cartesian<double> Kernel;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel                       Kernel;
typedef CGAL::Cartesian<NT> Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Curve_2                               Segment_2;

template <class GraphT>
class SelfIntersectionChecker {
private:
    const GraphT &graph; 
public:
    SelfIntersectionChecker(const GraphT &graph) : graph(graph) {        
    }   
    
    bool isSelfIntersecting(Chain &chain) {        
        
        if (chain.size() > 0) {                        
            
            std::vector<Segment_2> chainSegments;
            for (auto it = chain.begin(); it != --chain.end(); it++) {
                auto next = it;
                next++;
                
                
                //Point_2 p_1(graph.getLon(*it), graph.getLat(*it));
                //Point_2 p_2(graph.getLon(*next), graph.getLat(*next));
                
                Segment_2 segment(Point_2(graph.getLon(*it), graph.getLat(*it)), Point_2(graph.getLon(*next), graph.getLat(*next)));
                chainSegments.push_back(segment);
                
                
            }
                     
            return CGAL::do_curves_intersect(chainSegments.begin(), chainSegments.end());
            
        } else {//chain of zero nodes can't intersect itself
            return false;
        }       
    }
    
    void getChainSegments(Chain &chain, std::vector<Segment_2> &chainSegments) {
        assert(chain.size() >= 1);
        //std::vector<Segment_2> chainSegments;
        for (auto it = chain.begin(); it != --chain.end(); it++) {
            auto next = it;
            next++;

            /*
            Point_2 p_1(graph.getLon(*it), graph.getLat(*it));
            Point_2 p_2(graph.getLon(*next), graph.getLat(*next));

            Segment_2 segment(p_1, p_2);
            chainSegments.push_back(segment);
             * */
            
            double lon1 = graph.getLon(*it);
            double lat1 = graph.getLat(*it);
            double lon2 = graph.getLon(*next);
            double lat2 = graph.getLat(*next);

            //assert(!(lon1==lon2 && lat1==lat2));
            if (!(lon1==lon2 && lat1==lat2)) {  //shit can happen
                Point_2 p_1(lon1, lat1);
                Point_2 p_2(lon2, lat2);
                Segment_2 segment(p_1, p_2);
                chainSegments.push_back(segment);
            }
            
        }
        return;
    }
    
    bool isSelfIntersecting(Chain &chain1, Chain &chain2) {        
        
        assert(chain1.size() > 0 && chain2.size() > 0);
        std::vector<Segment_2> chainSegments;
        getChainSegments(chain1, chainSegments);
        getChainSegments(chain2, chainSegments);
        
        // Compute all intersection points.
                
        
        return CGAL::do_curves_intersect(chainSegments.begin(), chainSegments.end());               
    }
    
};