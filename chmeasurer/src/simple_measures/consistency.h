#pragma once

#include "../geoFunctions.h"
#include "../chgraph.h"
#include "../nodes_and_edges.h"
#include "../ch_measurer.h"
#include <math.h>

class ErrorConsistency {
    CHGraph<CHNode, CHEdge> &graph;        
    
public:
    ErrorConsistency(CHGraph<CHNode, CHEdge>& graph) :
        graph(graph) {
    }
            
    //mean absolute deviation, only useful for longer chains
    double calcMAD(const RedetectedChain &edge_chain, const double chain_length) {
        //calculate length of whole chain to know how to subdivide it in segments
        /*
        double edge_chain_length = 0;        
        for (EdgeID edge_id: edge_chain.edges) {
            auto edge  = graph.getEdge(edge_id);
            edge_chain_length += geo::geoDist(graph.getNode(edge.src), graph.getNode(edge.tgt));            
            
        }*/
        
        int nofSegments = abs(sqrt(edge_chain.edges.size()));    
        assert (nofSegments > 0);
        double pre_avg_segment_length = chain_length/nofSegments;
        
        //calculate summed error of segments
        std::list<double> summed_errors_per_segment;        
        
        double current_error_sum = 0;
        double cur_segment_length = 0;
        for (EdgeID edge_id: edge_chain.edges) {
            
            auto edge  = graph.getEdge(edge_id);
            current_error_sum += fabs(graph.calcEdgeError(edge_id));
            cur_segment_length += geo::geoDist(graph.getNode(edge.src), graph.getNode(edge.tgt));
            
            
            if (cur_segment_length > pre_avg_segment_length) {
                summed_errors_per_segment.push_back(current_error_sum);
                current_error_sum = 0;
                cur_segment_length -= pre_avg_segment_length;
            }            
        }
        
        //calc mean absolute vector error (mave)
        double sum = 0;
        for (double v_error: summed_errors_per_segment) {
            sum += v_error;
        }
        double mave = sum/nofSegments;
        
        //calc mean absolute deviation (mad)
        double quadratic_error_sum = 0;
        for (double v_error: summed_errors_per_segment) {
            quadratic_error_sum += pow(mave - v_error, 2);
        }
        double variance = quadratic_error_sum/nofSegments;
        return sqrt(variance);   
    }
};