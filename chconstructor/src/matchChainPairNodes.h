#pragma once

#include <limits>
#include <algorithm>

#include "self_intersection_checker.h"
#include "cdt_matching.h"
#include "nodes_and_edges.h"

#undef matchChainPairNodesNDEBUG


using namespace DP;
        
//namespace MCPN {
//static class
template <class GraphT>
class matchChainPairNodes {
    
    struct ZipNode {
        std::list<std::list<PPrioNode>::iterator>::iterator nodeIt; //iterator in an Intervall list
        //PPrioNodeRef &pnr;
        bool listDirection = true; //saves on which line the node was, true = To, false = From
        ZipNode() {};
        ZipNode(std::list<std::list<PPrioNode>::iterator>::iterator nodeIt, bool listDirection): nodeIt(nodeIt), listDirection(listDirection) {};
    };
    
    struct ZipSegment {
        ZipNode prevOnOtherSide;
        ZipNode nextOnOtherSide;
        std::list<ZipNode> zipNodes;
    };
    
    double R = 6371.009; //Erdradius
    
    const GraphT &base_graph;   
    //PIntervall &intervallTo;
    //PIntervall &intervallFrom;
    std::list<std::list<PPrioNode>::iterator> &toList;
    std::list<std::list<PPrioNode>::iterator> &fromList;
    
    double pythagoras(double a, double b) {
	return sqrt(pow(a, 2) + pow(b, 2));
    }
          
    double geoDist(double lon1, double lat1, double lon2, double lat2) {
        double deltaLon = lon1 - lon2;
        double deltaLat = lat1 - lat2;
        return pythagoras(deltaLon, deltaLat);        
        /*
        double latmed, rdlon, rdlat;
        latmed = (lat1 - lat2) / 2;
        rdlon = (lon1 - lon2) * M_PI / 180;
        rdlat = (lat1 - lat2) * M_PI / 180;
        return R * sqrt((pow(rdlat, 2) + pow(cos(latmed) * rdlon, 2)));        
         * */
    }    
    
    double geoDist(NodeID nodeID1, NodeID nodeID2) {
        double lon1 = base_graph.getLon(nodeID1);
        double lat1 = base_graph.getLat(nodeID1);               
        double lon2 = base_graph.getLon(nodeID2);
        double lat2 = base_graph.getLat(nodeID2);
        return geoDist(lon1, lat1, lon2, lat2);
    }

    //on the active side it is pointed to the element which is already pushed on zipOrder and on the inactive side it is the first element which is still not pushed
    ZipNode getNext(std::list<std::list<PPrioNode>::iterator>::iterator &toIt, std::list<std::list<PPrioNode>::iterator>::reverse_iterator &fromIt,
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
                return ZipNode(--(fromIt.base()), false);
            } else {
                if (fromIt == fromList.rend()) {
                    toIt++;
                    #ifdef matchChainPairNodesNDEBUG
                    auto dbg = (*toIt)->node_id;
                    #endif
                    return ZipNode(toIt, true);
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
                        return ZipNode(toIt, true);
                    } else {
                        activeDirection = false;
                        toIt++;
                        #ifdef matchChainPairNodesNDEBUG
                        auto dbg = (*--(fromIt.base()))->node_id;  
                        #endif
                        return ZipNode(--(fromIt.base()), false);
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
                return ZipNode(toIt, true);
            } else {
                if (toIt == toList.end()) {
                    fromIt++;
                    #ifdef matchChainPairNodesNDEBUG
                    auto dbg = (*--(fromIt.base()))->node_id;  
                    #endif
                    return ZipNode(--(fromIt.base()), false);
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
                        return ZipNode(--(fromIt.base()), false);
                    } else {
                        activeDirection = true;
                        fromIt++;
                        #ifdef matchChainPairNodesNDEBUG
                        auto dbg = (*toIt)->node_id;
                        #endif
                        return ZipNode(toIt, true);
                    }
                }
            }
        }
    }
    
    int getPreviousNodeOnOtherSide(vector<ZipNode> &zipOrder, int i) {
        debug_assert(0 <= i && i < (int) zipOrder.size());
        bool activeSide = zipOrder[i].listDirection;
        
        while (zipOrder[i].listDirection == activeSide) {
            i--;
            if(i == -1) {
                break;
            }
        }
        return i;                
    }
    
    int getNextNodeOnOtherSide(vector<ZipNode> &zipOrder, int i) {
        debug_assert(0 <= i && i < (int) zipOrder.size());
        bool activeSide = zipOrder[i].listDirection;
        
        while (zipOrder[i].listDirection == activeSide) {
            i++;  
            if(i == (int) zipOrder.size()) {
                break;
            }
        }
        return i;                
    }
    
    
    void setFollower (std::list<std::list<PPrioNode>::iterator>::iterator nodeIt,
                      std::list<std::list<PPrioNode>::iterator>::iterator followerIt) {
        (*nodeIt)->followerValid = true;
        (*nodeIt)->follower = *followerIt;            
        (*followerIt)->guides.push_back(*nodeIt);
    }
    
    void setFollowers (list<ZipNode> &zipNodes, std::list<std::list<PPrioNode>::iterator>::iterator followerIt) {
        for (ZipNode &zn : zipNodes) {
            setFollower(zn.nodeIt, followerIt); 
        }
    }
    
    void processSegment (ZipSegment &zs) {         
        if (zs.zipNodes.empty()) {
            return; //Recursion end
        } else {
            double minDist = numeric_limits<double>::max();
            bool minSide; //true==previous false == next
            auto minIt = zs.zipNodes.begin();
            
            for (auto it = zs.zipNodes.begin(); it != zs.zipNodes.end(); it++) {
                #ifdef matchChainPairNodesNDEBUG
                NodeID node_id1 = (*it->nodeIt)->node_id;                
                NodeID node_id2 = (*zs.prevOnOtherSide.nodeIt)->node_id;
                #endif
                
                double distPrevious = geoDist((*it->nodeIt)->node_id, (*zs.prevOnOtherSide.nodeIt)->node_id);
                if (distPrevious < minDist) {
                    minDist = distPrevious;
                    minSide = true;
                    minIt = it;
                }
                double distNext = geoDist((*it->nodeIt)->node_id, (*zs.nextOnOtherSide.nodeIt)->node_id);
                if(distNext < minDist) {
                    minDist = distNext;
                    minSide = false;
                    minIt = it;
                }
            }
            //NodeID node_id  = (*minIt->nodeIt)->node_id;
            std::list<ZipNode> followerDeterminedZipNodes;
            if (minSide == true) { //cut away front part                
                followerDeterminedZipNodes.splice (followerDeterminedZipNodes.begin(),
                    zs.zipNodes, zs.zipNodes.begin(), ++minIt);
                setFollowers(followerDeterminedZipNodes, zs.prevOnOtherSide.nodeIt);                
                processSegment(zs);
                
            } else {  //cut away back part
                #ifdef matchChainPairNodesNDEBUG
                NodeID node_id_next = (*zs.nextOnOtherSide.nodeIt)->node_id;
                #endif
                
                followerDeterminedZipNodes.splice (followerDeterminedZipNodes.begin(),
                    zs.zipNodes, minIt, zs.zipNodes.end());
                setFollowers(followerDeterminedZipNodes, zs.nextOnOtherSide.nodeIt);                
                processSegment(zs);                
            }
            
        }
    }
    
    void setFollowers (vector<ZipNode> &zipOrder, bool careOrdering) {
        
        if (careOrdering && toList.size()>=1 && fromList.size()>=1) {
            
                        
            ZipNode last_zn = ZipNode();
            bool active_Direction;
            //first Segment
            uint i = 0;            
            ZipSegment firstSegment;
            active_Direction = zipOrder.at(i).listDirection;
            firstSegment.zipNodes.push_back(zipOrder.at(i));
            i++;
            while (zipOrder.at(i).listDirection == active_Direction && i < zipOrder.size()) {
                firstSegment.zipNodes.push_back(zipOrder.at(i));
                i++;
            }
            firstSegment.nextOnOtherSide = zipOrder.at(i);
            setFollowers(firstSegment.zipNodes, firstSegment.nextOnOtherSide.nodeIt);            
            uint lastOfFirst = i-1;            
            
            //last Segment
            uint j = zipOrder.size()-1;            
            ZipSegment lastSegment;
            active_Direction = zipOrder.at(j).listDirection;
            lastSegment.zipNodes.push_front(zipOrder.at(j));
            j--;
            while (zipOrder.at(j).listDirection == active_Direction && j >= 0) {
                lastSegment.zipNodes.push_front(zipOrder.at(j));
                j--;
            }
            lastSegment.prevOnOtherSide = zipOrder.at(j);
            setFollowers(lastSegment.zipNodes, lastSegment.prevOnOtherSide.nodeIt); 
            uint firstOfLast = j+1;
                 
            //middle Segments
            uint k = lastOfFirst;
            ZipSegment middleSegment;
            last_zn = zipOrder.at(k);
            middleSegment.prevOnOtherSide = last_zn;
            k++;            
            //middleSegment.zipNodes.push_back(zipOrder.at(k));
            active_Direction = zipOrder.at(k).listDirection;
            //k++;
            while ( k <= firstOfLast) {  
                ZipNode active_zn = zipOrder.at(k);
                if(active_zn.listDirection == active_Direction) {
                    //push back another node to current segment
                    middleSegment.zipNodes.push_back(active_zn);
                    last_zn = active_zn;
                } else {
                    //process closed segment
                    middleSegment.nextOnOtherSide = active_zn;                    
                    processSegment(middleSegment);
                    
                    //start new segment
                    active_Direction = active_zn.listDirection;
                    middleSegment = ZipSegment(); //old segment is cleared
                    middleSegment.prevOnOtherSide = last_zn;
                    middleSegment.zipNodes.push_back(active_zn);   
                    last_zn = active_zn;
                }                
                k++;
            }                
            
        }
        else {
            for (uint i = 0; i < zipOrder.size(); i++) {      
                double distPrevious, distNext;
                //closest previous Node on the other chain according to zipOrder
                int previous = getPreviousNodeOnOtherSide(zipOrder, i);
                if (previous == -1) { //if there is no sooner node on the other side
                    distPrevious = numeric_limits<double>::max();
                } else {
                    distPrevious = geoDist((*zipOrder[i].nodeIt)->node_id, (*zipOrder[previous].nodeIt)->node_id);
                }

                //closest next Node on the other chain according to zipOrder
                int next = getNextNodeOnOtherSide(zipOrder, i);
                if (next == (int) zipOrder.size()) { //if there is no later node on the other side
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

                debug_assert(0<=follower && follower < (int) zipOrder.size());
                
                std::list<std::list<PPrioNode>::iterator>::iterator nodeIt = zipOrder[i].nodeIt;
                std::list<std::list<PPrioNode>::iterator>::iterator followerIt = zipOrder[follower].nodeIt;
                setFollower(nodeIt, followerIt);
                /*
                (*zipOrder[i].nodeIt)->followerValid = true;
                (*zipOrder[i].nodeIt)->follower = *zipOrder[follower].nodeIt;            
                (*zipOrder[follower].nodeIt)->guides.push_back(*zipOrder[i].nodeIt);
                 * */
            }   
        }                 
    }

    void zipOrder() {
        vector<ZipNode> zipOrder;
        size_t size = toList.size() + fromList.size();
        zipOrder.reserve(size);
        //construct an ordering to get follow-function

        bool activeDirection = true; //TODO: make it less arbitrary

        zipOrder.push_back(ZipNode(toList.begin(), activeDirection));

        auto toIt = toList.begin();
        auto fromIt = fromList.rbegin();

        for (uint i = 0; i < size - 1; i++) {
            ZipNode zipNode = getNext(toIt, fromIt, activeDirection);
            zipOrder.push_back(zipNode);
        }

        setFollowers(zipOrder, true);
    }
    
public:
    matchChainPairNodes(const GraphT &base_graph, PIntervall &intervallTo, PIntervall &intervallFrom)
        : base_graph(base_graph), toList(intervallTo.prioNodesIts), fromList(intervallFrom.prioNodesIts) {}
    
    void match() {
        debug_assert(toList.size()>=1 && fromList.size()>=1); 
        SelfIntersectionChecker<GraphT> selfIC(base_graph);
        Chain chainTo;
        for (auto it: toList) {
            chainTo.push_back(it->node_id);
        }
        Chain chainFrom;
        for (auto it: fromList) {
            chainFrom.push_back(it->node_id);
        }
        //additionally to make sure cdt_matching constraints do not not intersect 
        chainTo.push_back(chainFrom.front());
        chainFrom.push_back(chainTo.front());
        bool isSelfIntersecting = selfIC.isSelfIntersecting(chainTo, chainFrom);
        
        bool useZipOrder = isSelfIntersecting;
        if (useZipOrder) {
            zipOrder();            
        } else {            
            CDTMatching<GraphT> ctd_matcher(base_graph);
            ctd_matcher.match(toList, fromList);
        }
        return;
    }            
};
