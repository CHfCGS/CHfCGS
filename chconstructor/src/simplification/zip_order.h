#pragma once

#include <limits>
#include <algorithm>

//#include "self_intersection_checker.h"
//#include "cdt_matching.h"
#include "../nodes_and_edges.h"
#include "../geoFunctions.h"
#include "prio_nodes.h"
#include "match_functions.h"

#undef matchChainPairNodesNDEBUG


//using namespace ls;

namespace ls { 
    
namespace mc {
//namespace MCPN {
//static class
template <class GraphT>
class ZipOrder {
    
    struct ZipNode {
        //std::list<PrioNodeHandle>::const_iterator handleIt; //iterator in an Intervall list
        PrioNodeHandle pnh;
        //PPrioNodeRef &pnr;
        bool listDirection = 0; //saves on which line the node was, true = To, false = From
        ZipNode() {}
        //ZipNode(std::list<PrioNodeHandle>::const_iterator handleIt, bool listDirection):
        //    handleIt(handleIt), listDirection(listDirection) {};
        ZipNode(PrioNodeHandle pnh, bool listDirection):
            pnh(pnh), listDirection(listDirection) {}
    };
    
    struct ZipSegment {
        ZipNode prevOnOtherSide;
        ZipNode nextOnOtherSide;
        std::list<ZipNode> zipNodes;
        
    };
    
    const GraphT &graph;           
       
    //double geoDist(NodeID nodeID1, NodeID nodeID2) {
    double geoDist(const ZipNode zn1, const ZipNode zn2) const {
        //const PrioNode2 &pn1 = *(*zn1.handleIt);
        //const PrioNode2 &pn2 = *(*zn2.handleIt);
        const PrioNode2 &pn1 = *zn1.pnh;
        const PrioNode2 &pn2 = *zn2.pnh;
        return geo::geoDist(graph.getNode(pn1.node_id), graph.getNode(pn2.node_id));
        /*
        double lon1 = base_graph.getLon(nodeID1);
        double lat1 = base_graph.getLat(nodeID1);               
        double lon2 = base_graph.getLon(nodeID2);
        double lat2 = base_graph.getLat(nodeID2);
        return geoDist(lon1, lat1, lon2, lat2);
         * */
    }
    /*
    double geoDist(const PrioNodeHandle pnh1, const PrioNodeHandle pnh2) {
        const PrioNode2 &pn1 = *pnh1;
        const PrioNode2 &pn2 = *pnh2;
        return geo::geoDist(graph.getNode(pn1.node_id), graph.getNode(pn2.node_id));
    }*/

    //on the active side it is pointed to the element which is already pushed on zipOrder and on the inactive side it is the first element which is still not pushed
    ZipNode getNext(const std::list<PrioNodeHandle> &list1, const std::list<PrioNodeHandle> &list2,
                    std::list<PrioNodeHandle>::const_iterator &it1, std::list<PrioNodeHandle>::const_iterator &it2,
                    bool &activeDirection) const {
        
        //active direction 1 = true, 2 = false
        const std::list<PrioNodeHandle> &activeList = activeDirection ?  list1 : list2;
        const std::list<PrioNodeHandle> &passiveList = activeDirection ?  list2 : list1;
        std::list<PrioNodeHandle>::const_iterator &activeIt = activeDirection ?  it1 : it2;
        std::list<PrioNodeHandle>::const_iterator &passiveIt = activeDirection ?  it2 : it1;
        
        //PrioNode2 &pn = *(*activeIt);
        std::list<PrioNodeHandle>::const_iterator nextIt = activeIt;
        nextIt++;
        if (nextIt == activeList.end()) {
            activeIt++;
            activeDirection = !activeDirection;

            assert(passiveIt != passiveList.end());
            return ZipNode(*passiveIt, activeDirection);
        } else {
            if (passiveIt == passiveList.end()) {
                activeIt++;

                return ZipNode(*activeIt, activeDirection);
            } else {                                   
                double currentDistance = ls::mc::geoDist(*activeIt, *passiveIt, graph);
                double afterleapedDistance = ls::mc::geoDist(*nextIt, *passiveIt, graph);

                //jump on the other chain if distance to its next node becomes worse
                //i.e. stay on a side if the next node on the opposite chain comes closer by stepping forward
                if (currentDistance > afterleapedDistance) {
                    activeIt++;

                    return ZipNode(*activeIt, activeDirection);
                } else {
                    activeIt++;
                    activeDirection = !activeDirection;

                    return ZipNode(*passiveIt, activeDirection);
                }
            }
        }            
    }
    
    static int getPreviousNodeOnOtherSide(const std::vector<ZipNode> &zipOrder, int i) {
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
    
    static int getNextNodeOnOtherSide(const std::vector<ZipNode> &zipOrder, int i) {
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
    
    /*
    static void _setFollowerZN (const std::list<PrioNodeHandle>::const_iterator nodeIt,
                      const std::list<PrioNodeHandle>::const_iterator followerIt) {
        
        PrioNode2 &node = *(*nodeIt);
        PrioNode2 &follower = *(*followerIt);
        node.followerValid = true;
        node.follower_h = *followerIt;            
        
        const PrioNodeHandle pnh = *nodeIt;
        std::list<PrioNodeHandle> &guidelist = follower.guides;
        guidelist.push_back(pnh);
        
        //_setFollower(*nodeIt, *followerIt);        
    }*/
    
    static void _setFollowers (const std::list<ZipNode> &zipNodes, PrioNodeHandle follower_pnh) {
        for (const ZipNode &zn : zipNodes) {
           _setFollower(zn.pnh, follower_pnh); 
        }
    }
    
    void processSegment (ZipSegment &zs) const {         
        if (zs.zipNodes.empty()) {
            return; //Recursion end
        } else {
            double minDist = std::numeric_limits<double>::max();
            bool minSide; //true==previous false == next
            auto minIt = zs.zipNodes.begin();
            
            for (auto it = zs.zipNodes.begin(); it != zs.zipNodes.end(); it++) {               
                
                double distPrevious = geoDist(*it, zs.prevOnOtherSide);
                if (distPrevious < minDist) {
                    minDist = distPrevious;
                    minSide = true;
                    minIt = it;
                }
                double distNext = geoDist(*it, zs.nextOnOtherSide);
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
                _setFollowers(followerDeterminedZipNodes, zs.prevOnOtherSide.pnh);                
                processSegment(zs);
                
            } else {  //cut away back part
                
                followerDeterminedZipNodes.splice (followerDeterminedZipNodes.begin(),
                    zs.zipNodes, minIt, zs.zipNodes.end());
                _setFollowers(followerDeterminedZipNodes, zs.nextOnOtherSide.pnh);                
                processSegment(zs);                
            }
            
        }
    }
    
    void setFollowers(const std::vector<ZipNode> &zipOrder, bool careOrdering, 
            const std::list<PrioNodeHandle> &list1, const std::list<PrioNodeHandle> &list2) const {
        
        if (careOrdering && list1.size()>=1 && list2.size()>=1) {
                  
            //more sophisticated with respect to ordering
            
            bool active_Direction;
            //first Segment
            uint i = 0;            
            ZipSegment firstSegment;
            active_Direction = zipOrder.at(i).listDirection;
            firstSegment.zipNodes.push_back(zipOrder.at(i));
            i++;
            while (i < zipOrder.size() && zipOrder.at(i).listDirection == active_Direction) {
                firstSegment.zipNodes.push_back(zipOrder.at(i));
                i++;
            }
            firstSegment.nextOnOtherSide = zipOrder.at(i);
            _setFollowers(firstSegment.zipNodes, firstSegment.nextOnOtherSide.pnh);            
            uint lastOfFirst = i-1;            
            
            //last Segment
            uint j = zipOrder.size()-1;            
            ZipSegment lastSegment;
            active_Direction = zipOrder.at(j).listDirection;
            lastSegment.zipNodes.push_front(zipOrder.at(j));
            j--;
            while (j >= 0 && zipOrder.at(j).listDirection == active_Direction) {
                lastSegment.zipNodes.push_front(zipOrder.at(j));
                j--;
            }
            lastSegment.prevOnOtherSide = zipOrder.at(j);
            _setFollowers(lastSegment.zipNodes, lastSegment.prevOnOtherSide.pnh); 
            uint firstOfLast = j+1;
                 
            //middle Segments
            uint k = lastOfFirst;
            ZipSegment middleSegment;
            ZipNode last_zn = zipOrder.at(k);
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
            //simple and robust, could lead to intersections with respect to ordering
            for (uint i = 0; i < zipOrder.size(); i++) {      
                double distPrevious, distNext;
                //closest previous Node on the other chain according to zipOrder
                int previous = getPreviousNodeOnOtherSide(zipOrder, i);
                if (previous == -1) { //if there is no sooner node on the other side
                    distPrevious = std::numeric_limits<double>::max();
                } else {
                    distPrevious = geoDist(zipOrder[i], zipOrder[previous]);
                }

                //closest next Node on the other chain according to zipOrder
                int next = getNextNodeOnOtherSide(zipOrder, i);
                if (next == (int) zipOrder.size()) { //if there is no later node on the other side
                    distNext = std::numeric_limits<double>::max();
                } else {
                    distNext = geoDist(zipOrder[i], zipOrder[next]);
                }

                int follower;
                if (distPrevious<distNext) {
                    follower = previous;
                } else {
                    follower = next;
                }

                debug_assert(0<=follower && follower < (int) zipOrder.size());
                
                //std::list<PrioNodeHandle>::const_iterator nodeIt = zipOrder[i].handleIt;
                //std::list<PrioNodeHandle>::const_iterator followerIt = zipOrder[follower].handleIt;                
                _setFollower(zipOrder[i].pnh, zipOrder[follower].pnh);                
            }   
        }                 
    }

    
    
public:
    ZipOrder(const GraphT &graph)
        : graph(graph) {}
    
    void match(const std::list<PrioNodeHandle> &list1, const std::list<PrioNodeHandle> &list2) const {
        std::vector<ZipNode> zipOrder;
        size_t size = list1.size() + list2.size();
        zipOrder.reserve(size);
        //construct an ordering to get follow-function

        bool activeDirection = true; //there cant happen that much if we start arbitrarily
        zipOrder.push_back(ZipNode(list1.front(), activeDirection));

        auto it1 = list1.begin();
        auto it2 = list2.begin();

        for (uint i = 0; i < size - 1; i++) {
            ZipNode zipNode = getNext(list1, list2, it1, it2, activeDirection);
            zipOrder.push_back(zipNode);
        }

        setFollowers(zipOrder, true, list1, list2);
    }   
};

}

}
