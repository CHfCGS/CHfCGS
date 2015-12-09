#pragma once

#include "../geoFunctions.h"
#include "../chgraph.h"
#include "../nodes_and_edges.h"
#include <math.h>

class Regularity {
    CHGraph<CHNode, CHEdge> &graph;        
    
public:
    Regularity(CHGraph<CHNode, CHEdge>& graph) :
        graph(graph) {
    }
            
    double calcStandardDeviation(const Chain &chain, const double chain_length) {
        //assumption: redetectedChain is subsequence of expandedChain
                
        assert(chain.size() >= 2);        
        
        const uint edgeCount = chain.size()-1;
        /*
        double chain_length = 0;        
        for (auto it = chain.begin(); it!= --chain.end(); it++) {
            auto next = it;
            next++;
            chain_length += geo::geoDist(graph.getNode(*it), graph.getNode(*next));            
        }*/
        
        double avg_segment_length = chain_length/edgeCount;
        double quadratic_error_sum = 0;
        //standard deviation
        for (auto it = chain.begin(); it!= --chain.end(); it++) {
            auto next = it;
            next++;
            double edge_length = geo::geoDist(graph.getNode(*it), graph.getNode(*next));
            quadratic_error_sum += pow(avg_segment_length - edge_length, 2);
        }
        double variance = quadratic_error_sum/edgeCount;
        return sqrt(variance);        
    }
};