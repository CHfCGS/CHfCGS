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
#include "chgraph.h"

#include <algorithm>

namespace chm {

template <class GraphT>
class ChainDetector {
    
public:    
    ChainDetector(const GraphT &base_graph): base_graph(base_graph), marked(base_graph.getNrOfNodes(), false) {}
        
    virtual ~ChainDetector(){        
    };    
    
    Chains_and_Remainder detectChains(const std::vector<NodeID> &nodes) {
        //reset marked
        std::fill(marked.begin(), marked.end(), false);        
        
        Chains_and_Remainder CaR;    

        //collect all chains
        for (NodeID node_id : nodes) {        
            debug_assert(0 <= node_id && node_id < marked.size());
            if (marked[node_id] == false) {
                //if (base_graph.getNrOfEdges(node_id, chc::EdgeType::OUT) <= 2
                    //    && base_graph.getNrOfEdges(node_id, chc::EdgeType::IN) <= 2) {
                    //std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(node_id);
                    //if (neighbours.size() <= 2) {                        
                    if (base_graph.degree_leq_two(node_id)) {                        
                        collectChain(node_id, CaR);
                    }
                //}
            }
        }

        //collect remainder
        for (NodeID node_id : nodes) {
            debug_assert(0 <= node_id && node_id < marked.size());
            if (marked[node_id] == false) {
                CaR.remainder.push_back(node_id);
            }
        }
        return CaR;
    }        

    void collectChain(const NodeID node_id, Chains_and_Remainder &CaR) {
        Chain chain;

        debug_assert(!marked.at(node_id));
        chain.emplace_back(node_id);
        marked[node_id] = true;

        StreetType type = base_graph.getMaxStreetType(node_id);
        bool isOneway = base_graph.isOneway(node_id);
        
        //backward direction
        NodeID next = nextChainElement(node_id, type, EdgeType::IN, isOneway);
        while (next != c::NO_NID) {
            debug_assert(!marked.at(next));
            chain.emplace_front(next);
            marked[next] = true;
            NodeID current = next;
            next = nextChainElement(current, type, EdgeType::IN, isOneway);
        }
        
        //forward direction
        next = nextChainElement(node_id, type, EdgeType::OUT, isOneway);
        while (next != c::NO_NID) {
            debug_assert(!marked.at(next));
            chain.emplace_back(next);
            marked[next] = true;
            NodeID current = next;
            next = nextChainElement(current, type, EdgeType::OUT, isOneway);            
        }
        
        CaR.addChainOfType(chain, type, isOneway);
    }
    
    NodeID nextChainElement(const NodeID current, const StreetType streetType, const EdgeType edgeDirection, bool isOneway) {
        NodeID next = c::NO_NID;
        //TODO/Problem: chains cant share endpoints        
        //if (base_graph.getNrOfEdges(current, chc::EdgeType::OUT) <= 2 
          //      && base_graph.getNrOfEdges(current, chc::EdgeType::IN) <= 2) { //condtion optional but should increase performance
            //only neighbours which are reachable by a street of streettype
            //std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(current);
            //if (neighbours.size() <= 2) {
            if (base_graph.degree_leq_two(current)) {
                if (base_graph.isOneway(current) == isOneway) { //current node has to be oneway (twoway) if the chain is oneway (twoway)
                    std::list<NodeID> neighboursOnStreetType = base_graph.nodeNeighbours(current, streetType, edgeDirection);
                    debug_assert(neighboursOnStreetType.size() <=2);                                              
                    for (auto it = neighboursOnStreetType.begin(); it != neighboursOnStreetType.end(); it++) {
                        if (marked.at(*it) == false) {
                            next = *it;
                        }
                    }                
                }                
            }
        //}
        return next;
    }
    
    
    
    EdgeChain redetect(const Chain &old_chain, StreetType streettype) {
        
        Chain remaining_chain_nodes;
        
        for (NodeID node_id : old_chain) {
            if (base_graph.isValidNode(node_id)) {
                remaining_chain_nodes.push_back(node_id);
            }
        }
        
        
        bool valid = true;
        //first and last node have to remain
        if (old_chain.front() != remaining_chain_nodes.front()
                || old_chain.back() != remaining_chain_nodes.back()) {
            valid = false;
        }
        
        
        EdgeChain remaining_chain = EdgeChain(remaining_chain_nodes);
        
        if (valid) {
            //test connectivity
            for (auto it = remaining_chain_nodes.begin(); it != --remaining_chain_nodes.end(); it++) {
                //dont test for degree and oneway
                auto next = it;
                next++;
                
                //check if there is an edge connecting the two nodes and add the first
                bool found = false;
                //auto nodeEdges = base_graph.nodeEdges(*it, streettype);
                //Streettype streety
                NodeID node_id = 2;
                auto nodeEdges = base_graph.nodeEdges(node_id, streettype);
                for (auto edge: nodeEdges) {
                    if (edge.src == *next || edge.tgt == *next) {
                        remaining_chain.edges.push_back(edge.id);
                        found = true;
                        break;
                    }
                }
                valid = found;
                /*
                std::list<NodeID> nodeNeighbours =  base_graph.nodeNeighbours(*it, streettype);
                bool found = (std::find(nodeNeighbours.begin(), nodeNeighbours.end(), *next) != nodeNeighbours.end());
                if (!found) {
                    valid = false;
                    break;
                }
                 * */
            }        
        }
        
        if (valid) {
            return remaining_chain;
        } else {
            Chain emptyChain;
            return EdgeChain(emptyChain);
        }
    }
private:
    const GraphT &base_graph;
    std::vector<bool> marked;    
};

}