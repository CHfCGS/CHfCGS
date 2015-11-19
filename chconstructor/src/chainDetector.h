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
class ChainDetector {
    
public:
    //ChainDetector(const GraphT &base_graph);
    ChainDetector(const GraphT &base_graph): base_graph(base_graph), marked(base_graph.getNrOfNodes(), false) {                 
        //marked.resize(base_graph.getNrOfNodes());
        //marked.assign(base_graph.getNrOfNodes(), false);
    }
    
    //ChainDetector();
    //ChainDetector(const ChainDetector& orig);
    virtual ~ChainDetector(){        
    };    
    
    Chains_and_Remainder detectChains(const std::vector<chc::NodeID> &nodes) {
        //reset marked
        std::fill(marked.begin(), marked.end(), false);        

        //std::vector<chc::Chain> chains;
        Chains_and_Remainder CaR;
        //std::vector<bool> marked(base_graph.getNrOfNodes(), false);

        //collect all chains
        for (chc::NodeID node_id : nodes) {
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
        for (chc::NodeID node_id : nodes) {
            debug_assert(0 <= node_id && node_id < marked.size());
            if (marked[node_id] == false) {
                CaR.remainder.push_back(node_id);
            }
        }
        return CaR;
    }        

    void collectChain(const chc::NodeID node_id, Chains_and_Remainder &CaR) {
        Chain chain;

        debug_assert(!marked.at(node_id));
        chain.emplace_back(node_id);
        marked[node_id] = true;

        chc::StreetType type = base_graph.getMaxStreetType(node_id);
        bool isOneway = base_graph.isOneway(node_id);
        
        //backward direction
        chc::NodeID next = nextChainElement(node_id, type, chc::EdgeType::IN, isOneway);
        while (next != chc::c::NO_NID) {
            debug_assert(!marked.at(next));
            chain.emplace_front(next);
            marked[next] = true;
            chc:: NodeID current = next;
            next = nextChainElement(current, type, chc::EdgeType::IN, isOneway);
        }
        
        //forward direction
        next = nextChainElement(node_id, type, chc::EdgeType::OUT, isOneway);
        while (next != chc::c::NO_NID) {
            debug_assert(!marked.at(next));
            chain.emplace_back(next);
            marked[next] = true;
            chc:: NodeID current = next;
            next = nextChainElement(current, type, chc::EdgeType::OUT, isOneway);            
        }
        
        CaR.addChainOfType(chain, type, isOneway);
    }
    
    chc::NodeID nextChainElement(const chc::NodeID current, const chc::StreetType streetType, const chc::EdgeType edgeDirection, bool isOneway) {
        chc::NodeID next = chc::c::NO_NID;
        //TODO/Problem: chains cant share endpoints        
        //if (base_graph.getNrOfEdges(current, chc::EdgeType::OUT) <= 2 
          //      && base_graph.getNrOfEdges(current, chc::EdgeType::IN) <= 2) { //condtion optional but should increase performance
            //only neighbours which are reachable by a street of streettype
            //std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(current);
            //if (neighbours.size() <= 2) {
            if (base_graph.degree_leq_two(current)) {
                if (base_graph.isOneway(current) == isOneway) { //current node has to be oneway (twoway) if the chain is oneway (twoway)
                    std::list<chc::NodeID> neighboursOnStreetType = base_graph.nodeNeighbours(current, streetType, edgeDirection);
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
private:
    const GraphT &base_graph;
    std::vector<bool> marked;    
};