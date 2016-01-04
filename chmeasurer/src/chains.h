#pragma once

#include "nodes_and_edges.h"
#include "chgraph.h"

#include <limits>
//using namespace chm;


namespace chm {
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
typedef std::list<EdgeID> EdgeChain; 

namespace chains {
    
    
    static EdgeChain getChainEdges(const Chain& chain, const CHGraph<CHNode, CHEdge>& graph) {
        //assumption: connectivity is given
        EdgeChain edges;
        assert(!chain.empty());
        for (auto it = chain.begin(); it != --chain.end(); it++) {
            auto next_it = it;
            next_it++;
            NodeID next_node_id = *next_it;
            bool edge_found = false;
            for (const CHEdge& edge: graph.valid_nodeEdges(*it, EdgeType::OUT)) {                
                //assumption: no loops
                if (edge.tgt == next_node_id) {// || edge.tgt == next_node_id) {
                    edges.push_back(edge.id);
                    edge_found = true;
                    break;
                }                
            }
            assert(edge_found);
        }
        assert(edges.size() == chain.size()-1);
        return edges;
    }
}



struct EdgeChainPair {
    EdgeChain chainTo;
    EdgeChain chainFrom; //reversed Problem: get CenterNodes need to be reversed too
    
    EdgeChainPair() {}

    
    EdgeChainPair(const ChainPair& chain_pair, const CHGraph<CHNode, CHEdge>& graph) {
        this->chainTo = chains::getChainEdges(chain_pair.chainTo, graph);
        this->chainFrom = chains::getChainEdges(chain_pair.chainFrom, graph);                
    }    
};

struct RedetectedChain {
    //EdgeChain(Chain &chain): chain(chain) {}
    Chain remaining_chain;
    //edges between the nodes of chain
    std::list<EdgeID> edges;
    //enclosing chain in full graph
    Chain hull_chain;
};


struct Chains_and_Remainder {

    Chains_and_Remainder(uint max_street_type) : remainder(), oneWayChainsAccordingToType() {
        oneWayChainsAccordingToType.resize(max_street_type+1);
        twoWayChainsAccordingToType.resize(max_street_type+1);//different street Types //TODO: calculate maxStreetType beforehand        
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

namespace chains {
    
    bool areEqual(const Chain &chain1, const Chain &chain2) {
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
    
    //fill up an edge chain with the nodes in between (the center nodes of the connecting edges)
    static Chain expandEdgeChain(const RedetectedChain &redetected_edge_chain, const CHGraph<CHNode, CHEdge>& graph) {
        assert(redetected_edge_chain.remaining_chain.size()>=2 &&redetected_edge_chain.edges.size()>=1);
        assert(redetected_edge_chain.remaining_chain.size() == redetected_edge_chain.edges.size()+1);
        Chain expandedChain;

        auto nodeIt = redetected_edge_chain.remaining_chain.begin();
        expandedChain.push_back(*nodeIt);
        nodeIt++;
        for (auto edgeIt = redetected_edge_chain.edges.begin(); edgeIt!= redetected_edge_chain.edges.end(); edgeIt++) {
            std::list<NodeID> centerNodes = graph.getCenterNodes(*edgeIt);
                        
            for (NodeID node_id: centerNodes) {
                expandedChain.push_back(node_id);
            }                
            expandedChain.push_back(*nodeIt);
            nodeIt++;
        }
        assert(nodeIt == redetected_edge_chain.remaining_chain.end());

        return expandedChain;
    }

    static double calcChainGeoLength(const Chain &chain, const CHGraph<CHNode, CHEdge>& graph) {
        assert(chain.size() >= 2);                            
        double chain_length = 0;

        for (auto it = chain.begin(); it!= --chain.end(); it++) {
            auto next = it;
            next++;
            chain_length += geo::geoDist(graph.getNode(*it), graph.getNode(*next));            
        }
        return chain_length;            
    }
    
    static Chain toNodeChain(const EdgeChain& edge_chain, const CHGraph<CHNode, CHEdge>& graph) {
        assert(edge_chain.size() >= 1);
        Chain chain;
        NodeID src = graph.getEdge(edge_chain.front()).src;
        chain.push_back(src);
        
        for (EdgeID edge_id: edge_chain) {
            
            const CHEdge& edge = graph.getEdge(edge_id);                                    
            assert(edge.src == src);                                   
            chain.push_back(edge.tgt);
            src = edge.tgt;
        }
        return chain;
    }
    
    static Chain toExpandedNodeChain(const EdgeChain& edge_chain, const CHGraph<CHNode, CHEdge>& graph) {
        assert(edge_chain.size() >= 1);
        Chain chain;
        NodeID src = graph.getEdge(edge_chain.front()).src;
        chain.push_back(src);
        
        for (EdgeID edge_id: edge_chain) {
            
            const CHEdge& edge = graph.getEdge(edge_id);                                    
            assert(edge.src == src);
                        
            Chain center_nodes = graph.getCenterNodes(edge_id);            
            chain.splice(chain.end(), center_nodes);
            
            chain.push_back(edge.tgt);
            src = edge.tgt;
        }
        return chain;
    }
    
    static uint getExpandedNodeLength (const EdgeChain& chain, const CHGraph<CHNode, CHEdge>& graph) {
        uint length = 0; 
        for (const EdgeID edge_id: chain) {                        
            length += graph.getCenterNodes(edge_id).size() + 1;            
        }
        return length;
    }
    
    static std::list<EdgeChain> split (EdgeChain& chain, const CHGraph<CHNode, CHEdge>& graph) {
        std::list<EdgeChain> split_chains;
        const uint critical_size = 40;
        uint length = getExpandedNodeLength(chain, graph);        
        
        //only long chains are splitted
        if (length >= critical_size && chain.size() >=2) {
            uint collectedNofEdges = 0;
            uint toCollect = length/2;
                        
            //prevents degenerated splitting
            collectedNofEdges += graph.getCenterNodes(chain.front()).size() + 1;         
            auto split_it = ++chain.begin();
            
            for (auto it = ++chain.begin(); it != --chain.end(); it++) {
                collectedNofEdges += graph.getCenterNodes(*it).size() + 1;         
                if (collectedNofEdges >= toCollect) {
                    split_it = it;
                    break;
                }
            }
            EdgeChain front_part;
            EdgeChain back_part;
            uint chain_size_before = chain.size();
            front_part.splice(front_part.begin(), chain, chain.begin(), split_it); 
            assert(front_part.size()>=1);
            back_part.splice(back_part.begin(), chain); //rest
            assert(back_part.size()>=1);
            assert(front_part.size() + back_part.size() == chain_size_before);
                        
            std::list<EdgeChain> output_front = split(front_part, graph);
            std::list<EdgeChain> output_back = split(back_part, graph);
            
            split_chains.splice(split_chains.end(), output_front);
            split_chains.splice(split_chains.end(), output_back);                                                
            
        } else {
            split_chains.push_back(chain); 
        }   
        return split_chains;
    }   
          
    struct EdgeChainParts {
        EdgeChainPair front_part;
        EdgeChainPair back_part;
        
        EdgeChainParts(EdgeChain& edge_chain1, EdgeChain::iterator split_it1,
                        EdgeChain& edge_chain2, EdgeChain::iterator split_it2,
                        bool chain1_direction) {            
            EdgeChain* chainTo;
            EdgeChain* chainFrom;
            EdgeChain::iterator chainTo_it;
            EdgeChain::iterator chainFrom_it;
             
            if (chain1_direction) {
                chainTo = &edge_chain1;
                chainTo_it = split_it1;
                chainFrom = &edge_chain2;
                chainFrom_it = split_it2;
                
            } else {
                chainTo = &edge_chain2;
                chainTo_it = split_it2;
                chainFrom = &edge_chain1;
                chainFrom_it = split_it1;
            }
            front_part.chainTo.splice(front_part.chainTo.end(), *chainTo, chainTo->begin(), chainTo_it);
            back_part.chainTo.splice(back_part.chainTo.end(), *chainTo);
                        
            back_part.chainFrom.splice(back_part.chainFrom.end(), *chainFrom, chainFrom->begin(), chainFrom_it);
            front_part.chainFrom.splice(front_part.chainFrom.end(), *chainFrom);
            assert(edge_chain1.empty());
            assert(edge_chain2.empty());
        }
        
    };
    
    /*
    static bool uniqueElements (const Chain chain1, const Chain chain2) {
        std::list<NodeID> nodes;
        for (NodeID node_id: chain1) {
            nodes.push_back(node_id);
        }
        for (NodeID node_id: chain2) {
            nodes.push_back(node_id);
        }
        nodes.sort();
        uint size_before = nodes.size();
        nodes.unique();
        return (nodes.size()== size_before);
    }*/
    
    static bool sameLocation(NodeID node_id1, NodeID node_id2, const CHGraph<CHNode, CHEdge>& graph) {
        const CHNode& node1 = graph.getNode(node_id1);
        const CHNode& node2 = graph.getNode(node_id2);
        return node1.lat == node2.lat && node1.lon == node2.lon;
    }
    
    static bool uniqueLocations(const Chain& chain1, const Chain& chain2, const CHGraph<CHNode, CHEdge>& graph) {
        std::vector<NodeID> chain;
        for (NodeID node_id: chain1) {
            chain.push_back(node_id);
        }
        for (NodeID node_id: chain2) {
            chain.push_back(node_id);
        }
        //pairwise comparision
        bool unique = true;        
        for (uint i = 0; i < chain.size()-1; i++) {
            for (uint j = i+1; j < chain.size(); j++) {
                if (sameLocation(chain[i], chain[j], graph)) {
                    unique = false;
                }
            }                
        }
        return unique;
    }
    
    static bool uniqueElements(const Chain& chain1, const Chain& chain2) {
        std::vector<NodeID> chain;
        for (NodeID node_id: chain1) {
            chain.push_back(node_id);
        }
        for (NodeID node_id: chain2) {
            chain.push_back(node_id);
        }
        //pairwise comparision
        bool unique = true;        
        for (uint i = 0; i < chain.size()-1; i++) {
            for (uint j = i+1; j < chain.size(); j++) {
                if (chain[i]== chain[j]) {
                    unique = false;
                }
            }                
        }
        return unique;
    }
    
    static std::list<EdgeChainPair> split (EdgeChainPair& chain_pair, const CHGraph<CHNode, CHEdge>& graph) {        
        std::list<EdgeChainPair> edge_chain_pairs;
        const uint critical_size = 40;
        //get sizes
        EdgeChain* p_bigger_chain = nullptr;
        EdgeChain* p_smaller_chain = nullptr;
        bool bigger_chain_direction;
        if (chain_pair.chainTo.size() > chain_pair.chainFrom.size()) {
            p_bigger_chain = &chain_pair.chainTo;
            p_smaller_chain = &chain_pair.chainFrom;
            bigger_chain_direction = true;
        } else {
            p_bigger_chain = &chain_pair.chainFrom;
            p_smaller_chain = &chain_pair.chainTo;
            bigger_chain_direction = false;
        }
        EdgeChain& bigger_chain = *p_bigger_chain;
        EdgeChain& smaller_chain = *p_smaller_chain;
       
        //length in expanded edges
        uint bigger_chain_length = getExpandedNodeLength(bigger_chain, graph);
        uint smaller_chain_length = getExpandedNodeLength(smaller_chain, graph);
        
        //split in 2 edgechains
        if (bigger_chain.size()>=2// && smaller_chain.size()>=1
                && smaller_chain_length + bigger_chain_length >= critical_size) {
            //get split pos in bigger chain
            uint collectedNofEdges = 0;
            uint toCollect = bigger_chain_length/2;
                        
            //prevents degenerated splitting
            collectedNofEdges += graph.getCenterNodes(bigger_chain.front()).size() + 1;         
            auto bigger_split_it = ++bigger_chain.begin();
            
            for (auto it = ++bigger_chain.begin(); it != --bigger_chain.end(); it++) {
                collectedNofEdges += graph.getCenterNodes(*it).size() + 1;         
                if (collectedNofEdges >= toCollect) {
                    bigger_split_it = it;
                    break;
                }
            }
            
            //get nearest node in smaller chain
            double minDist = std::numeric_limits<double>::max();
            auto smaller_split_it = smaller_chain.end();
            for (auto it = smaller_chain.begin(); it != smaller_chain.end(); it++) {
                auto node1 = graph.getNode(graph.getEdge(*bigger_split_it).src);
                auto node2 = graph.getNode(graph.getEdge(*it).src);
                double dist = geo::geoDist(node1, node2);
                //double dist = geo::geoDist(graph.getNode(*bigger_split_it), graph.getNode(*it));
                if (dist < minDist) {
                    minDist = dist;
                    smaller_split_it = it;
                }                
            }
            
            EdgeChainParts edge_chain_parts(bigger_chain, bigger_split_it, smaller_chain, smaller_split_it, bigger_chain_direction);
            
            std::list<EdgeChainPair> output_front = split(edge_chain_parts.front_part, graph);
            std::list<EdgeChainPair> output_back = split(edge_chain_parts.back_part, graph);
            
            edge_chain_pairs.splice(edge_chain_pairs.end(), output_front);
            edge_chain_pairs.splice(edge_chain_pairs.end(), output_back);  
        } else {
            edge_chain_pairs.push_back(chain_pair);
        }
        return edge_chain_pairs;
    }
    
    
    
    /*
    static std::list<Chain> splitChain(Chain& chain, uint critical_size) {
        std::list<Chain> split_chains;
        //n-th chain
        int nof_segments = ceil( (double) chain.size()/(double) critical_size);
        int segment_length = ceil( (double) chain.size()/(double) nof_segments);
        for (int i = 0; i < nof_segments-1; i++) {
            Chain split_chain;
            for (int j = 0; j < segment_length; j++) {
                split_chain.push_back(chain.front()); //more efficient: splice
                chain.pop_front();
            }
            split_chains.push_back(split_chain);
        }
        //last part
        Chain last_split_chain;
        for (NodeID node_id: chain) {
            last_split_chain.push_back(node_id);
        }
        split_chains.push_back(last_split_chain);

        return split_chains;
    }*/
       
    /*
    static void splitChains(ChainsOfType &chains) {            
        for (auto it = chains.begin(); it != chains.end(); it++) {   
            Chain &chain = *it;
            //split chains larger than a critical size
            uint critical_size = 50;
            if (chain.size() > critical_size) {

                std::list<Chain> split_chains = splitChain(chain, critical_size);
                //delete old chain and hang in split chains instead
                it = chains.erase(it);
                chains.splice(it, split_chains, split_chains.begin(), split_chains.end());                                        
            }                                
        }
    }*/
}

}