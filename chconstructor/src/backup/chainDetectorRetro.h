/* 
 * File:   chainDetector.h
 * Author: tobias
 *
 * Created on 16. September 2015, 14:29
 */

//#ifndef CHAINDETECTOR_H
//#define	CHAINDETECTOR_H
#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "chains.h"

template <class GraphT>
class ChainDetector {
public:
    //ChainDetector(const GraphT &base_graph);
    ChainDetector(const GraphT &base_graph):base_graph(base_graph), marked(base_graph.getNrOfNodes(), false) {
        Print("Chaindetector constructed");
    }
    
    //ChainDetector();
    //ChainDetector(const ChainDetector& orig);
    virtual ~ChainDetector(){};
    
    Chains_and_Remainder detectChains(const std::vector<chc::NodeID> &nodes);
    //chc::Chain startChain(const chc::NodeID node, std::list<chc::NodeID> &neighbours);
    Chain startChain(const chc::NodeID node_id);
    chc::NodeID nextChainElement(const chc::NodeID previous, const chc::NodeID current, const chc::StreetType type);
private:
    const GraphT &base_graph;
    std::vector<bool> marked;    
};

template <class GraphT>
Chains_and_Remainder ChainDetector<GraphT>::detectChains(const std::vector<chc::NodeID> &nodes) {
    //reset marked
    std::fill(marked.begin(), marked.end(),false);
    
    //std::vector<chc::Chain> chains;
    Chains_and_Remainder CaR;
    //std::vector<bool> marked(base_graph.getNrOfNodes(), false);
    
    //collect all chains
    for (chc::NodeID node_id : nodes) {        
        if (marked[node_id]==false) {            
            if (base_graph.getNrOfEdges(node_id, chc::EdgeType::OUT)<=2 
                && base_graph.getNrOfEdges(node_id, chc::EdgeType::IN)<=2)
            {
                std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(node_id);    
                if (neighbours.size()<=2) {
                    Chain chain;
                    chain = startChain(node_id);
                    CaR.addChainOfType(chain, 1);
                }
            } 
        }                 
    }
    
    //collect remainder
    for (chc::NodeID node_id : nodes) {        
        if (marked[node_id]==false) {
            CaR.remainder.push_back(node_id);         
        }
    }
    return CaR;
}

template <class GraphT>
//chc::Chain ChainDetector<GraphT>::startChain(const chc::NodeID node_id, std::list<chc::NodeID> &neighbours) {
Chain ChainDetector<GraphT>::startChain(const chc::NodeID node_id) {
    Chain chain;        
    
    assert(!marked[node_id]);
    chain.emplace_back(node_id);    
    marked[node_id] = true;
    
    chc::StreetType type =  base_graph.getMaxStreetType(node_id);
    //std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(node_id);
    std::list<chc::NodeID> neighboursOnStreetType = base_graph.nodeNeighbours(node_id, type);
    
    debug_assert(neighboursOnStreetType.size()<=2);
    if (neighboursOnStreetType.size()==1) {
        //go in one direction
        chc::NodeID current = node_id;
        chc::NodeID next = neighboursOnStreetType.front();
        chc::NodeID oldnext;
        if (marked[next] == false) {
            while (next != chc::c::NO_NID) {
                assert(!marked[next]);
                chain.emplace_back(next);                

                oldnext = next;
                next = nextChainElement(current, next, type);
                marked[oldnext] = true;
                current = oldnext;
            }
        }

    } else if (neighboursOnStreetType.size()==2) {
        //go in both directions        
        
        //front direction
        chc::NodeID current = node_id;
        chc::NodeID next = neighboursOnStreetType.front();
        chc::NodeID oldnext;
        if (marked[next] == false) {
            while (next != chc::c::NO_NID) {
                assert(!marked[next]);
                chain.emplace_front(next);                                

                oldnext = next;
                next = nextChainElement(current, next, type);
                marked[oldnext] = true;
                current = oldnext;
            }
        }
        
        //back direction
        current = node_id;
        next = neighboursOnStreetType.back();        
        if (marked[next] == false) {
            while (next != chc::c::NO_NID) {
                assert(!marked[next]);
                chain.emplace_back(next);                                

                oldnext = next;
                next = nextChainElement(current, next, type);
                marked[oldnext] = true;
                current = oldnext;
            }       
        }

                 
    }
    return chain;           
}

template <class GraphT>
chc::NodeID ChainDetector<GraphT>::nextChainElement(const chc::NodeID previous, const chc::NodeID current, const chc::StreetType type) {
    chc::NodeID next = chc::c::NO_NID;
    //TODO/Problem: chains cant share endpoints
    //if (marked[current] == false) {
        if (base_graph.getNrOfEdges(current, chc::EdgeType::OUT) <= 2
                && base_graph.getNrOfEdges(current, chc::EdgeType::IN) <= 2) {
            //only neigbours which are reachable by a street of streettype
            std::list<chc::NodeID> neighbours = base_graph.nodeNeighbours(current);
            if (neighbours.size() == 2) {
                std::list<chc::NodeID> neighboursOnStreetType = base_graph.nodeNeighbours(current, type);
                if (neighboursOnStreetType.size() == 2) {
                    neighboursOnStreetType.remove(previous);
                    debug_assert(neighboursOnStreetType.size() == 1);
                    next = neighboursOnStreetType.front();
                    if (marked[next] == true) {
                        next = chc::c::NO_NID;
                    }
                }
            }
            
        }
    return next;
}


//#endif	/* CHAINDETECTOR_H */
