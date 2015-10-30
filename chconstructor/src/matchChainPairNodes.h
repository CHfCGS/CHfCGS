#pragma once

#include <limits>

#include "nodes_and_edges.h"

#undef matchChainPairNodesNDEBUG


using namespace DP;
        
//namespace MCPN {
//static class
template <class GraphT>
class matchChainPairNodes {
    
    struct zipNode {
        std::list<std::list<PPrioNode>::iterator>::iterator nodeIt; //iterator in an Intervall list
        //PPrioNodeRef &pnr;
        bool listDirection; //saves on which line the node was, true = To, false = From
        
        zipNode(std::list<std::list<PPrioNode>::iterator>::iterator nodeIt, bool listDirection): nodeIt(nodeIt), listDirection(listDirection) {};
    };
    
    double R = 6371.009; //Erdradius
    
    const GraphT &base_graph;   
    //PIntervall &intervallTo;
    //PIntervall &intervallFrom;
    std::list<std::list<PPrioNode>::iterator> &toList;
    std::list<std::list<PPrioNode>::iterator> &fromList;
    
    double geoDist(double lon1, double lat1, double lon2, double lat2) {
        double latmed, rdlon, rdlat;
        latmed = (lat1 - lat2) / 2;
        rdlon = (lon1 - lon2) * M_PI / 180;
        rdlat = (lat1 - lat2) * M_PI / 180;
        return R * sqrt((pow(rdlat, 2) + pow(cos(latmed) * rdlon, 2)));        
    }    
    
    double geoDist(NodeID nodeID1, NodeID nodeID2) {
        double lon1 = base_graph.getLon(nodeID1);
        double lat1 = base_graph.getLat(nodeID1);               
        double lon2 = base_graph.getLon(nodeID2);
        double lat2 = base_graph.getLat(nodeID2);
        return geoDist(lon1, lat1, lon2, lat2);
    }

    //on the active side it is pointed to the element which is already pushed on zipOrder and on the inactive side it is the first element which is still not pushed
    zipNode getNext(std::list<std::list<PPrioNode>::iterator>::iterator &toIt, std::list<std::list<PPrioNode>::iterator>::reverse_iterator &fromIt,
            bool &activeDirection) {
        
        //active direction To = true, From = false
        //std::list<PPrioNodeRef> &activeChain, &otherChain;
        //std::list<PPrioNodeRef>::iterator &activeChainIt, &otherChainIt;
        
        if (activeDirection) {
            std::list<std::list<PPrioNode>::iterator>::iterator nextIt = toIt;
            nextIt++;
            if (nextIt == toList.end()) { 
                activeDirection = false;
                toIt++;                  
                #ifdef matchChainPairNodesNDEBUG
                auto dbg = (*--(fromIt.base()))->node_id;  
                #endif
                return zipNode(--(fromIt.base()), false);
            } else {
                if (fromIt == fromList.rend()) {
                    toIt++;
                    #ifdef matchChainPairNodesNDEBUG
                    auto dbg = (*toIt)->node_id;
                    #endif
                    return zipNode(toIt, true);
                } else {
                    double currentDistance = geoDist((*toIt)->node_id, (*fromIt)->node_id);                                                
                    double afterleapedDistance = geoDist((*nextIt)->node_id, (*fromIt)->node_id);                            

                    //jump on the other chain if distance to its next node becomes worse
                    //i.e. stay on a side if the next node on the opposite chain comes closer by stepping forward
                    if (currentDistance > afterleapedDistance) {
                        toIt++;
                        #ifdef matchChainPairNodesNDEBUG
                        auto dbg = (*toIt)->node_id;
                        #endif
                        return zipNode(toIt, true);
                    } else {
                        activeDirection = false;
                        toIt++;
                        #ifdef matchChainPairNodesNDEBUG
                        auto dbg = (*--(fromIt.base()))->node_id;  
                        #endif
                        return zipNode(--(fromIt.base()), false);
                    }
                }                                
            }
        }
        else {                        
            std::list<std::list<PPrioNode>::iterator>::reverse_iterator nextIt = fromIt;
            nextIt++;
            if (nextIt == fromList.rend()) {
                activeDirection = true;
                fromIt++;
                #ifdef matchChainPairNodesNDEBUG
                auto dbg = (*toIt)->node_id;
                #endif
                return zipNode(toIt, true);
            } else {
                if (toIt == toList.end()) {
                    fromIt++;
                    #ifdef matchChainPairNodesNDEBUG
                    auto dbg = (*--(fromIt.base()))->node_id;  
                    #endif
                    return zipNode(--(fromIt.base()), false);
                } else {
                    double currentDistance = geoDist((*fromIt)->node_id, (*toIt)->node_id);                            
                    double afterleapedDistance = geoDist((*nextIt)->node_id, (*toIt)->node_id);                            

                    //jump on the other chain if distance to its next node becomes worse
                    //i.e. stay on a side if the next node on the opposite chain comes closer by stepping forward
                    if (currentDistance > afterleapedDistance) {
                        fromIt++;
                        #ifdef matchChainPairNodesNDEBUG
                        auto dbg = (*--(fromIt.base()))->node_id;  
                        #endif
                        return zipNode(--(fromIt.base()), false);
                    } else {
                        activeDirection = true;
                        fromIt++;
                        #ifdef matchChainPairNodesNDEBUG
                        auto dbg = (*toIt)->node_id;
                        #endif
                        return zipNode(toIt, true);
                    }
                }
            }
        }
    }

    int getPreviousNodeOnOtherSide(vector<zipNode> &zipOrder, int i) {
        debug_assert(0 <= i && i < zipOrder.size());
        bool activeSide = zipOrder[i].listDirection;
        
        while (zipOrder[i].listDirection == activeSide) {
            i--;
            if(i == -1) {
                break;
            }
        }
        return i;                
    }
    
    int getNextNodeOnOtherSide(vector<zipNode> &zipOrder, int i) {
        debug_assert(0 <= i && i < zipOrder.size());
        bool activeSide = zipOrder[i].listDirection;
        
        while (zipOrder[i].listDirection == activeSide) {
            i++;  
            if(i == zipOrder.size()) {
                break;
            }
        }
        return i;                
    }
    
    void setFollowers (vector<zipNode> &zipOrder) {
        for (uint i = 0; i < zipOrder.size(); i++) {      
            double distPrevious, distNext;
            //closest previous Node on the other chain according to zipOrder
            int previous = getPreviousNodeOnOtherSide(zipOrder, i);
            if (previous==-1) { //if there is no sooner node on the other side
                distPrevious = numeric_limits<double>::max();
            } else {
                distPrevious = geoDist((*zipOrder[i].nodeIt)->node_id, (*zipOrder[previous].nodeIt)->node_id);
            }
                        
            //closest next Node on the other chain according to zipOrder
            int next = getNextNodeOnOtherSide(zipOrder, i);
            if (next==zipOrder.size()) { //if there is no later node on the other side
                distNext = numeric_limits<double>::max();
            } else {
                distNext = geoDist((*zipOrder[i].nodeIt)->node_id, (*zipOrder[next].nodeIt)->node_id);
            }

            int follower;
            if (distPrevious<distNext) {
                follower = previous;
            } else {
                follower = next;
            }
            
            debug_assert(0<=follower && follower < zipOrder.size());
            
            (*zipOrder[i].nodeIt)->followerValid = true;
            (*zipOrder[i].nodeIt)->follower = *zipOrder[follower].nodeIt;            
            (*zipOrder[follower].nodeIt)->guides.push_back(*zipOrder[i].nodeIt);
        }

    }    
public:
    matchChainPairNodes(const GraphT &base_graph, PIntervall &intervallTo, PIntervall &intervallFrom)
        : base_graph(base_graph), toList(intervallTo.prioNodesIts), fromList(intervallFrom.prioNodesIts) {}
    
    void match() {
        //vector<list<PPrioNodeRef>::iterator> ordering;
        debug_assert(toList.size()>=1 && fromList.size()>=1);        
        
        vector<zipNode> zipOrder;
        size_t size = toList.size() + fromList.size();
        zipOrder.reserve(size);
        //construct an ordering to get follow-function
                
        bool activeDirection = true; //TODO: make it less arbitrary
        
        zipOrder.push_back(zipNode(toList.begin(), activeDirection));
                
        auto toIt = toList.begin();        
        auto fromIt = fromList.rbegin();
        
        for (int i = 0; i < size-1; i++) {
            zipOrder.push_back(getNext(toIt, fromIt, activeDirection));            
        }
        
        setFollowers(zipOrder);
        
        return;
    }            
};
