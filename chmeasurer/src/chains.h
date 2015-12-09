#pragma once

#include "nodes_and_edges.h"

using namespace chm;


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


bool chainsAreEqual(const Chain &chain1, const Chain &chain2) {
    if (chain1.size() == chain2.size()) {
        auto it2 = chain2.begin();
        for (auto it1 = chain1.begin(); it1 != chain1.end(); it1++) {
            if (*it1 == *it2) {
                it2++;
            } else {
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}


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

struct RedetectedChain {
    //EdgeChain(Chain &chain): chain(chain) {}
    Chain remaining_chain;
    //edges between the nodes of chain
    std::list<EdgeID> edges;
    //enclosing chain in full graph
    Chain hull_chain;
};


struct Chains_and_Remainder {

    Chains_and_Remainder() : remainder(), oneWayChainsAccordingToType() {
        oneWayChainsAccordingToType.resize(20);
        twoWayChainsAccordingToType.resize(20);//different street Types //TODO: calculate maxStreetType beforehand        
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