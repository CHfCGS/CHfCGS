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
    
    //detect only normal dead ends
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
        
    struct DeadEndCandidate {
        Chain deadEnd;
        NodeID start;
        NodeID finish;
        bool isExtendedDeadEnd(const GraphT& base_graph) {
            //special case
            if (deadEnd.size()==1) {
                //front is only element
                if (base_graph.nodeNeighbours(deadEnd.front()).size() == 2) {
                    return false;
                }                
            }
            
            if (start == chc::c::NO_NID)
                return true;
            else if(finish == chc::c::NO_NID)
                return true;
            else if (start == finish) //circle
                return true;
            else
                return false;
        }
    };
    
    std::list<Chain> detectExtendedDeadEnds(const std::vector<chc::NodeID> &nodes) {
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
                    if (base_graph.degree_leq(node_id, 2)) {
                        DeadEndCandidate dee = collectDeadEndCandidate(node_id);
                        if (dee.isExtendedDeadEnd(base_graph)) {                            
                            deadEnds.push_back(dee.deadEnd);
                        }                        
                    }
                //}
            }
        }        
        return deadEnds;
    }
    
    DeadEndCandidate collectDeadEndCandidate(const chc::NodeID node_id) {
        DeadEndCandidate dee;
        //Chain deadEnd;

        debug_assert(!marked.at(node_id));
        dee.deadEnd.emplace_back(node_id);
        marked[node_id] = true;
        
        //backward direction
        NodeID next = nextDeadEndElement(node_id);
        while (next != chc::c::NO_NID) {
            debug_assert(!marked.at(next));
            dee.deadEnd.emplace_front(next);                        
            marked[next] = true;            
            chc:: NodeID current = next;
            next = nextDeadEndElement(current);            
        }
        //forward direction (the direction which isn't marked yet)
        next = nextDeadEndElement(node_id);
        while (next != chc::c::NO_NID) {
            debug_assert(!marked.at(next));
            dee.deadEnd.emplace_back(next);            
            marked[next] = true;            
            chc:: NodeID current = next;
            next = nextDeadEndElement(current);            
        }
        dee.start = getBorder(dee.deadEnd.front());
        dee.finish = getBorder(dee.deadEnd.back());
        
        return dee;
    }

    chc::NodeID nextDeadEndCandidateElement(const chc::NodeID current) {
        chc::NodeID next = chc::c::NO_NID;

        std::list<chc::NodeID> nodeNeighbours = base_graph.nodeNeighbours(current);
        if(nodeNeighbours.size() <= 2) {
            for (auto it = nodeNeighbours.begin(); it != nodeNeighbours.end(); it++) {            
                if (marked.at(*it) == false) {
                    if (base_graph.degree_leq(*it, 2)) {
                        next = *it;
                    }
                }
            }
        }
        return next;
    }
    
    chc::NodeID getBorder(const chc::NodeID current) {
        chc::NodeID border = chc::c::NO_NID;
        std::list<chc::NodeID> nodeNeighbours = base_graph.nodeNeighbours(current);
        debug_assert(nodeNeighbours.size() <= 2);
        for (auto it = nodeNeighbours.begin(); it != nodeNeighbours.end(); it++) {
            //nodes with degree bigger 2 are never marked
            if (marked.at(*it) == false) {
                //if (base_graph.degree_leq(*it, 2)) {
                    border = *it;
                //}
            }
        }
        return border;
    }
    
private:
    const GraphT &base_graph;
    std::vector<bool> marked;    
};

#pragma once

//#include "../s_options.h"

//struct SOptions;


//    namespace ls {
    //enum class LineSimplifierType { NONE = 0, DP };
//}
    
enum class DeadEndDetectType { NONE = 0, DE, EDE};

static constexpr DeadEndDetectType LastDeadEndDetectType = DeadEndDetectType::EDE;

DeadEndDetectType toDeadEndDetectType(std::string const& type)
{
	if (type == "NONE") {
            return DeadEndDetectType::NONE;
	}
        else if (type == "DE") {
            return DeadEndDetectType::DE;
	}
        else if (type == "EDE") {
            return DeadEndDetectType::EDE;
	}        
	return DeadEndDetectType::NONE;
}

std::string to_string(DeadEndDetectType type)
{
	switch (type) {
	case DeadEndDetectType::NONE:
		return "NONE";	
        case DeadEndDetectType::DE:
		return "DE";	
        case DeadEndDetectType::EDE:
		return "EDE";
        }
	std::cerr << "Unknown PairMatchType type: " << static_cast<int>(type) << "\n";
	return "NONE";
        
}