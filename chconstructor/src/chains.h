#pragma once

#include "nodes_and_edges.h"


using namespace chc;
/*
struct Chain {
    Chain() : node_ids() {}
    
    std::list<NodeID> node_ids;
    StreetType type = 0;
};
*/
typedef std::list<NodeID> Chain;

typedef std::list<Chain> ChainsOfType;
/*
struct PairChainNode;


typedef std::list<PairChainNode> PairChain; //only one of the chains


struct PairChainNode {
    NodeID nodeID;
    PairChain::iterator parter; //pointer to the element which is on the opposite way
    
};
 * */

struct ChainPair {
    Chain chainTo;
    Chain chainFrom;
    /*
    //function to reestablish a correct partner function
    remove(PairChain::iterator to_remove, bool direction) {
        //direction says whether a node is to removed of To(true) or From(false)
        if (direction) {
            
        }

    }
     * */
    ChainPair(): chainTo(), chainFrom() {}
};

/*
struct ChainsOfType {
    ChainsOfType() : chains() {}
    //StreetType type = 0;
    std::list<Chain> chains;
    
};
 * */

//typedef std::list<ChainsOfType> ChainsAccordingToType;

/*
struct ChainsAccordingToType {
    
};
 */

struct Chains_and_Remainder {

    Chains_and_Remainder(uint max_street_type) : remainder(), oneWayChainsAccordingToType() {
        oneWayChainsAccordingToType.resize(max_street_type + 1);
        twoWayChainsAccordingToType.resize(max_street_type + 1);//different street Types //TODO: calculate maxStreetType beforehand        
        //Print("Default-constructed capacity is " << chainsAccordingToType.capacity());
        //chainsAccordingToType.shrink_to_fit();
        //Print("after shrink " << chainsAccordingToType.capacity());
        
        
    }

    
    void addChainOfType(Chain chain, StreetType type, bool isOneWay) {
              
        if (isOneWay) {
            debug_assert(type < oneWayChainsAccordingToType.size());        
            oneWayChainsAccordingToType.at(type).push_back(chain); //use vector index as streettype
        }
        else {
            debug_assert(type < twoWayChainsAccordingToType.size());        
            twoWayChainsAccordingToType.at(type).push_back(chain); //use vector index as streettype
        }

        
    }
    
    uint getNrOfChains() {
        uint counter = 0;
        for (uint i = 0; i < oneWayChainsAccordingToType.size(); i++) {            
            counter += oneWayChainsAccordingToType.at(i).size();
        }
        for (uint i = 0; i < twoWayChainsAccordingToType.size(); i++) {            
            counter += twoWayChainsAccordingToType.at(i).size();
        }
        return counter;
    }
    
    uint getNrOfNodesInChains() {
        uint counter = 0;
        for (uint i = 0; i < oneWayChainsAccordingToType.size(); i++) {
            for (Chain &chain : oneWayChainsAccordingToType.at(i)) {
                counter += chain.size();
            }
        }
        for (uint i = 0; i < twoWayChainsAccordingToType.size(); i++) {
            for (Chain &chain : twoWayChainsAccordingToType.at(i)) {
                counter += chain.size();
            }
        }
        return counter;
    }

    std::vector<NodeID> remainder;
    std::vector<ChainsOfType> oneWayChainsAccordingToType;
    std::vector<ChainsOfType> twoWayChainsAccordingToType;
    std::vector<ChainPair> chainPairs;

};

#include "discreteFrechet.h"

namespace chains {
    template <class GraphT>
    static void cutDown(ChainPair& chain_pair, const GraphT& graph) {
        if (chain_pair.chainTo.size()<=3 || chain_pair.chainFrom.size()<=3) {
            return;
        } else {
            df::DiscreteFrechet<GraphT> df(graph);
            double df_dist = df.calc_dF(chain_pair.chainTo, chain_pair.chainFrom);
            double distFront =  geo::geoDist(graph.getNode(chain_pair.chainTo.front()), graph.getNode(chain_pair.chainFrom.back()));
            double distBack =  geo::geoDist(graph.getNode(chain_pair.chainTo.back()), graph.getNode(chain_pair.chainFrom.front()));
            if (df_dist <= distFront + std::numeric_limits<double>::epsilon()) {                
                double distToFrom = geo::geoDist(graph.getNode(*(chain_pair.chainTo.begin())),
                                             graph.getNode(*(++chain_pair.chainFrom.rbegin())));
                double distFromTo = geo::geoDist(graph.getNode(*(++chain_pair.chainTo.begin())),
                                             graph.getNode(*(chain_pair.chainFrom.rbegin())));
                if (distToFrom >= distFromTo) {
                    chain_pair.chainTo.pop_front();                    
                } else {
                    chain_pair.chainFrom.pop_back();   
                }
                cutDown(chain_pair, graph);
            } else if (df_dist <= distBack + std::numeric_limits<double>::epsilon()) {
                double distToFrom = geo::geoDist(graph.getNode(*(chain_pair.chainTo.rbegin())),
                                             graph.getNode(*(++chain_pair.chainFrom.begin())));
                double distFromTo = geo::geoDist(graph.getNode(*(++chain_pair.chainTo.rbegin())),
                                             graph.getNode(*(chain_pair.chainFrom.begin())));
                if (distToFrom >= distFromTo) {
                    chain_pair.chainTo.pop_back();                    
                } else {
                    chain_pair.chainFrom.pop_front();   
                }
                cutDown(chain_pair, graph);
            } else {
                return;
            }            
            
        }        
    }
}
