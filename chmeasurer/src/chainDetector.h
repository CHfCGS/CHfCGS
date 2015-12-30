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
    ChainDetector(const GraphT &base_graph): graph(base_graph), marked(base_graph.getNrOfNodes(), false) {}
        
    virtual ~ChainDetector(){        
    };  
        
    
    
    Chains_and_Remainder detectChains(const std::vector<NodeID> &nodes) {
        //reset marked
        std::fill(marked.begin(), marked.end(), false);        
        
        Chains_and_Remainder CaR(graph.getMaxStreetType());    

        //collect all chains
        for (NodeID node_id : nodes) {        
            debug_assert(0 <= node_id && node_id < (int) marked.size());
            if (marked[node_id] == false) {
                //if (base_graph.getNrOfEdges(node_id, chc::EdgeType::OUT) <= 2
                    //    && base_graph.getNrOfEdges(node_id, chc::EdgeType::IN) <= 2) {
                    //std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(node_id);
                    //if (neighbours.size() <= 2) {                        
                    if (graph.degree_leq_two(node_id)) {                        
                        collectChain(node_id, CaR);
                    }
                //}
            }
        }

        //collect remainder
        for (NodeID node_id : nodes) {
            debug_assert(0 <= node_id && node_id < (int) marked.size());
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

        StreetType type = graph.getMinStreetType(node_id);
        bool isOneway = graph.isOneway(node_id);
        
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
            if (graph.degree_leq_two(current)) {
                if (graph.isOneway(current) == isOneway) { //current node has to be oneway (twoway) if the chain is oneway (twoway)
                    std::list<NodeID> neighboursOnStreetType = graph.nodeNeighbours(current, streetType, edgeDirection);
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
    
    Chain getHullChain(const Chain &old_chain, const Chain &remaining_chain) {
        assert(old_chain.size() >= 2);
        assert(remaining_chain.size() >= 2);
        Chain full_chain;
        
        //get start
        Chain::const_iterator first = old_chain.end();
        for (auto it = old_chain.begin(); it != old_chain.end(); it++) {
            if (*it == remaining_chain.front()) {
                first = it;
                break;
            }
        }
        //get end
        Chain::const_iterator last = old_chain.end();
        for (Chain::const_reverse_iterator rit = old_chain.rbegin(); rit != old_chain.rend(); rit++) {
            if (*rit == remaining_chain.back()) {                
                last = --rit.base();
                assert(*last == remaining_chain.back());
                break;
            }
        }
        
        assert(first != old_chain.end());
        assert(last != old_chain.end());
        Chain::const_iterator end = last;
        end++;
        //copy part in between
        for (auto it = first; it != end; it++) {
            full_chain.push_back(*it);
        }
        return full_chain;
    }    
    
    RedetectedChain redetect(const Chain &old_chain, StreetType streettype) {     
        assert(old_chain.size() >=2);
        Chain remaining_chain_nodes;
               
        for (NodeID node_id : old_chain) {
            if (graph.isValidNode(node_id)) {
                remaining_chain_nodes.push_back(node_id);
            }
        }
                
        bool valid = true;
        /*
        //first and last node have to remain
        if (old_chain.front() != remaining_chain_nodes.front()
                || old_chain.back() != remaining_chain_nodes.back()) {
            valid = false;
        }*/
        if (remaining_chain_nodes.size() < 2) {
            valid = false;
        }
        
        
        RedetectedChain redetected_chain;
        redetected_chain.remaining_chain = remaining_chain_nodes;
        
        
        if (valid) {
            redetected_chain.hull_chain = getHullChain(old_chain, redetected_chain.remaining_chain);
            //test connectivity
            for (auto it = remaining_chain_nodes.begin(); it != --remaining_chain_nodes.end(); it++) {
                //dont test for degree and oneway
                auto next = it;
                next++;
                
                //check if there is an edge connecting the two nodes and add the first
                bool found = false;
                //auto nodeEdges = base_graph.nodeEdges(*it, streettype);
                
                
                auto nodeEdges = graph.nodeEdges(*it, streettype);
                for (auto edge: nodeEdges) {
                    if (edge.src == *next || edge.tgt == *next) {
                        redetected_chain.edges.push_back(edge.id);
                        found = true;
                        break;
                    }
                }
                if (found == false) {
                    valid = false;
                }                
            }                    
        }
        
        if (valid) {
            assert(redetected_chain.remaining_chain.size()>=2);
            assert(redetected_chain.edges.size()>=1);
            assert(redetected_chain.remaining_chain.size() == redetected_chain.edges.size()+1);
            return redetected_chain;
        } else {
            //Chain emptyChain;
           return RedetectedChain();
        }
    }
    RedetectedChain redetect(const Chain &old_chain) {
        assert(old_chain.size() >=2);
        Chain remaining_chain_nodes;
        
        for (NodeID node_id : old_chain) {
            if (graph.isValidNode(node_id)) {
                remaining_chain_nodes.push_back(node_id);
            }
        }
        
        
        bool valid = true;
        /*
        //first and last node have to remain
        if (old_chain.front() != remaining_chain_nodes.front()
                || old_chain.back() != remaining_chain_nodes.back()) {
            valid = false;
        }*/
        if (remaining_chain_nodes.size() < 2) {
            valid = false;
        }
        
        
        RedetectedChain redetected_chain;
        redetected_chain.remaining_chain = remaining_chain_nodes;
        
        if (valid) {
            redetected_chain.hull_chain = getHullChain(old_chain, redetected_chain.remaining_chain);
            //test connectivity
            for (auto it = remaining_chain_nodes.begin(); it != --remaining_chain_nodes.end(); it++) {
                
                auto next = it;
                next++;
                
                //check if there is an edge connecting the two nodes and add the first
                bool found = false;
                //auto nodeEdges = base_graph.nodeEdges(*it, streettype);
                                                               
                
                auto nodeEdges = graph.nodeEdges(*it);
                for (auto edge: nodeEdges) {
                    if (edge.src == *next || edge.tgt == *next) {
                        redetected_chain.edges.push_back(edge.id);
                        found = true;
                        break;
                    }
                }
                if (found == false) {
                    valid = false;
                }                
            }                    
        }
        
        if (valid) {
            assert(redetected_chain.remaining_chain.size()>=2);
            assert(redetected_chain.edges.size()>=1);
            assert(redetected_chain.remaining_chain.size() == redetected_chain.edges.size()+1);
            return redetected_chain;
        } else {
            //Chain emptyChain;
           return RedetectedChain();
        }
    }
    
    
    
    Chain expandChain (const Chain& chain) {
        Chain expanded_chain;
        assert(!chain.empty());
        assert(chain.size() >= 3);
        expanded_chain.push_back(chain.front());
        for (auto it = chain.begin(); it != --chain.end(); it++) {
            auto next_it = it;
            next_it++;
            NodeID next_node_id = *next_it;
            bool edge_found = false;            
            for (const CHEdge& edge: graph.valid_nodeEdges(*it, EdgeType::OUT)) {                
                //assumption: no loops
                if (edge.tgt == next_node_id) {// ||edge.src == next_node_id) {
                    Chain center_nodes = graph.getCenterNodes(edge.id);
                    for (NodeID node_id: center_nodes) {
                        expanded_chain.push_back(node_id);
                    }
                    //expanded_chain.splice(expanded_chain.end(), center_nodes, center_nodes.begin(), center_nodes.end());
                    edge_found = true;
                    break;
                }                
            }            
            assert(edge_found);
            expanded_chain.push_back(next_node_id);            
        } 
        /*
        expanded_chain.sort();
        uint sizeBefore = expanded_chain.size();
        expanded_chain.unique();
        assert(sizeBefore == expanded_chain.size());
        */
        return expanded_chain;
    }
    
private:
    const GraphT &graph;
    std::vector<bool> marked;    
};

}