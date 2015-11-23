/* 
 * File:   chainDetector.h
 * Author: tobias
 *
 * Created on 16. September 2015, 14:29
 */

#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "chains.h"
#include "graph.h"

template <class GraphT>
class DeadEndDetector {
    
public:
    //ChainDetector(const GraphT &base_graph);
    DeadEndDetector(const GraphT &base_graph): base_graph(base_graph), marked(base_graph.getNrOfNodes(), false) {                 
        //marked.resize(base_graph.getNrOfNodes());
        //marked.assign(base_graph.getNrOfNodes(), false);
    }
    
    //ChainDetector();
    //ChainDetector(const ChainDetector& orig);
    virtual ~DeadEndDetector(){        
    };    
    
    std::list<Chain> detectDeadEnds(const std::vector<chc::NodeID> &nodes) {
        //reset marked
        std::fill(marked.begin(), marked.end(), false);        

        std::list<Chain> deadEnds;
        //Chains_and_Remainder CaR;
        //std::vector<bool> marked(base_graph.getNrOfNodes(), false);

        //collect all chains
        for (chc::NodeID node_id : nodes) {
            debug_assert(0 <= node_id && node_id < marked.size());
            if (marked[node_id] == false) {
                //if (base_graph.getNrOfEdges(node_id, chc::EdgeType::OUT) <= 2
                    //    && base_graph.getNrOfEdges(node_id, chc::EdgeType::IN) <= 2) {
                    //std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(node_id);
                    //if (neighbours.size() <= 2) {                        
                    if (base_graph.degree_leq(node_id, 1)) {                        
                        deadEnds.push_back(collectDeadEnd(node_id));
                    }
                //}
            }
        }
        /*
        //collect remainder
        for (chc::NodeID node_id : nodes) {
            debug_assert(0 <= node_id && node_id < marked.size());
            if (marked[node_id] == false) {
                CaR.remainder.push_back(node_id);
            }
        }
         * */
        return deadEnds;
    }        

    Chain collectDeadEnd(const chc::NodeID node_id) {
        Chain deadEnd;

        debug_assert(!marked.at(node_id));
        deadEnd.emplace_back(node_id);
        marked[node_id] = true;
        
        //forward direction
        NodeID next = nextDeadEndElement(node_id);
        while (next != chc::c::NO_NID) {
            debug_assert(!marked.at(next));
            deadEnd.emplace_back(next);
            marked[next] = true;
            chc:: NodeID current = next;
            next = nextDeadEndElement(current);            
        }
        
        return deadEnd;
    }

    chc::NodeID nextDeadEndElement(const chc::NodeID current) {
        chc::NodeID next = chc::c::NO_NID;

        std::list<chc::NodeID> nodeNeighbours = base_graph.nodeNeighbours(current);
        debug_assert(nodeNeighbours.size() <= 2);
        for (auto it = nodeNeighbours.begin(); it != nodeNeighbours.end(); it++) {            
            if (marked.at(*it) == false) {
                if (base_graph.degree_leq(*it, 2)) {
                    next = *it;
                }
            }
        }
        return next;
    }
private:
    const GraphT &base_graph;
    std::vector<bool> marked;    
};