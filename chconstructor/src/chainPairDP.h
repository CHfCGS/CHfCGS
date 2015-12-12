#pragma once

#include "matchChainPairNodes.h"

#include "nodes_and_edges.h"
#include "simplification/prio_nodes.h"
#include "grid.h"
#include "chains.h"

#include <limits>
#include <math.h>

namespace ls {


template <class GraphT>
class chainPairDP{
    private:       
        std::list<PPrioNode> nodePrios;
        std::list<PIntervall> intervalls;        
        const GraphT &base_graph;
        const Grid<GraphT> &grid;
        std::list<std::vector<NodeInBox>> boxesWithNodes;
        
        //naiv TODO: heap implementation;
        std::list<PPrioNode>::iterator getMax() {
            
            assert(!nodePrios.empty());
            nodePrios.sort();
            std::list<PPrioNode>::iterator maxIt = --(nodePrios.end());            

            //PPrioNode r(*maxIt);            
            //nodePrios.erase(maxIt);
            return maxIt;
        }
    
    std::list<simplePrioNode> simplify() {
        if (nodePrios.empty()) {
            return std::list<simplePrioNode>(); //empty list
        }
        else {
            std::list<PPrioNode>::iterator maxIt = getMax();
            simplePrioNode spnMax(*maxIt);
            //PPrioNode max = getMax(nodePrios);            
            split(maxIt);
                                              
            //update guidance of max, i.e. it can't follow any other node anymore
            for (std::list<PPrioNode>::iterator guide : maxIt->guides) {
                guide->followerValid = false;
            }
            
            //check if there is a follower
            if (maxIt->followerValid) {
                std::list<PPrioNode>::iterator followerIt = maxIt->follower;
                simplePrioNode spnFollower(*followerIt);
                //nodePrios.erase(max.follower);
                split(followerIt);
                  
                 //update guidance of follower
                for (std::list<PPrioNode>::iterator guide : followerIt->guides) {
                    guide->followerValid = false;
                }
                
                if (followerIt->followerValid) {
                    //the follower of follower has one guide less now
                    size_t nofFollowerFollowerguides = followerIt->follower->guides.size();
                    followerIt->follower->guides.remove(followerIt);
                    debug_assert(nofFollowerFollowerguides-1 == followerIt->follower->guides.size());
                }

                nodePrios.erase(maxIt);
                nodePrios.erase(followerIt);
                std::list<simplePrioNode> r = simplify();
                //PrioNodes have to be pushed in that order such that most important is last
                r.push_back(simplePrioNode(spnFollower));                
                r.push_back(simplePrioNode(spnMax));
                
                return r;
            } else {
                nodePrios.erase(maxIt);
                std::list<simplePrioNode> r = simplify();
                //r.push_back(simplePrioNode(max));
                //std::list<simplePrioNode> r = split();//(max.intervallIt); //, max.posInIntervall);                    
                r.push_back(spnMax);
                return r;
            }
        }
    }
                
    void split(std::list<PPrioNode>::iterator splitNodeIt) {                
            
            std::list<PIntervall>::iterator splitIt = splitNodeIt->intervallIt;            
            std::list<std::list<PPrioNode>::iterator>::iterator splitPos = splitNodeIt->posInIntervallIt;                                                                            
            
            PIntervall leftIntervall(splitIt->start, (*splitPos)->node_id);
            PIntervall rightIntervall((*splitPos)->node_id, splitIt->end);                        
            
            std::list<PIntervall>::iterator leftIt = intervalls.insert(splitIt, leftIntervall);           
            std::list<PIntervall>::iterator rightIt = intervalls.insert(splitIt, rightIntervall);                                                                                                
            
            //left Intervall
            leftIt->prioNodesIts.splice(leftIt->prioNodesIts.begin(), splitIt->prioNodesIts, splitIt->prioNodesIts.begin(), splitPos);
            //leftIt->prioNodes = std::list<PrioNodeRef>(splitIt->prioNodes.begin(), splitPos);
            //cout << "left: \n";            
            for (std::list<std::list<PPrioNode>::iterator>::iterator it = leftIt->prioNodesIts.begin(); it!= leftIt->prioNodesIts.end(); it++) {
                //cout << nodes[it->get().node_id].coord.lon << " " << nodes[it->get().node_id].coord.lat << "\n";
                (*it)->intervallIt = leftIt;
                (*it)->posInIntervallIt = it;                
            }

            //right Intervall  
            rightIt->prioNodesIts.splice(rightIt->prioNodesIts.begin(), splitIt->prioNodesIts, ++splitPos, splitIt->prioNodesIts.end());
            //rightIt->prioNodes = std::list<PrioNodeRef>(++splitPos, splitIt->prioNodes.end());
            //cout << "right: \n";            
            for (std::list<std::list<PPrioNode>::iterator>::iterator it = rightIt->prioNodesIts.begin(); it!= rightIt->prioNodesIts.end(); it++) {                
                //cout << nodes[it->get().node_id].coord.lon << " " << nodes[it->get().node_id].coord.lat << "\n";
                (*it)->intervallIt = rightIt;
                (*it)->posInIntervallIt = it;                 
            }                                    
                        
            intervalls.erase(splitIt);
            
                        
            update (*leftIt);
            update (*rightIt);                                               
            
            return;               
    }        
    
    inline double OrientationTest(NodeID aSource, NodeID bTarget, NodeID cOutlier) {
        double ax = base_graph.getLon(aSource);
        double ay = base_graph.getLat(aSource);
        double bx = base_graph.getLon(bTarget);
        double by = base_graph.getLat(bTarget);
        double cx = base_graph.getLon(cOutlier);
        double cy = base_graph.getLat(cOutlier);
        
        double acx = ax - cx; // pa[0] - pc[0];
        double bcx = bx - cx; //pb[0] - pc[0];
        double acy = ay - cy; //pa[1] - pc[1];
        double bcy = by - cy; //pb[1] - pc[1];
        return acx * bcy - acy * bcx;
    }
    
    inline bool intersect(NodeID a, NodeID b, NodeID c, NodeID d) {
        double sideOfc = OrientationTest(a, b, c);
        double sideOfd = OrientationTest(a, b, d);
        return (sideOfc * sideOfd) < 0;
    }
    
    /*
     * calculates how many edges cross the lines from start to split and from split to end
     */
    uint calcNoFIntersections(NodeID start, NodeID split, NodeID end){
        uint intersectionCount = 0;
        
        //neighbors in the sense of being geographically nearby, not necessarily connected by edges
        std::vector<chc::NodeID> geoNeighbours(grid.nodeGeoNeighbours(split));        
        for (chc::NodeID geoNeighbour: geoNeighbours) {
            //neigbours by edges
            auto neighbours = base_graph.nodeNeighbours(geoNeighbour);
            for (chc::NodeID neighbour: neighbours) {
                if (intersect(start, split, geoNeighbour, neighbour)) {
                    intersectionCount++;
                }
                if (intersect(split, end, geoNeighbour, neighbour)) {
                    intersectionCount++;
                }
            }                       
        }
        
        return intersectionCount;
    }
    
    /*
    calculates how many nodes are correctly oriented after the intervall
    is spanned by the lines from start to split and from split to end
    */
    uint calcOrientationMisses(const PPrioNode &split, const PIntervall &intervall){
        uint misscount = 0;
        uint assertcounter = 0;
        auto splitPosIntervallincr= split.posInIntervallIt;
        splitPosIntervallincr++;
        for (auto it = intervall.prioNodesIts.begin(); it != splitPosIntervallincr; it++) {                                    
            for(NodeInBox nodeInBox : (*it)->leftBox) {                
                bool sign = (OrientationTest(intervall.start, split.node_id, nodeInBox.nodeID) >= 0) ? true : false;  
                if (sign != nodeInBox.sign) {
                    misscount++;
                }
            }                  
            assertcounter++;
        }
        for (auto it = split.posInIntervallIt; it != intervall.prioNodesIts.end(); it++) {
            for(NodeInBox nodeInBox : (*it)->rightBox) {
                bool sign = (OrientationTest(split.node_id, intervall.end, nodeInBox.nodeID) >= 0) ? true : false;  
                if (sign != nodeInBox.sign) {
                    misscount++;
                }
            }                                
            assertcounter++;
        }        
        debug_assert(assertcounter == intervall.prioNodesIts.size()+1);
        return misscount;
    }
    
    /*
     * updates the priority data for all prionodes in an Intervall
     */    
    void update(PIntervall &intervall){
        
        struct vec{
            double x;
            double y;            
        };
        
        //as if the earth would be flat, i.e lat and lon orthogonal, fails at the poles
        double x1 = base_graph.getLon(intervall.start); //prioNodes.front().get().node.coord.lon;
        double y1 = base_graph.getLat(intervall.start); //prioNodes.front().get().node.coord.lat;
        double x2 = base_graph.getLon(intervall.end); //prioNodes.back().get().node.coord.lon;
        double y2 = base_graph.getLat(intervall.end); //prioNodes.back().get().node.coord.lat;
        
        vec parallelVec = {x2-x1, y2-y1};
        //vec orthogonalVec = {parallelVec.y, -parallelVec.x};
        
        
        /*
        //line of the intervall endpoints
        double gradient_Intervall = (y2-y1)/(x2-x1);
        double c_Intervall = y1-gradient_Intervall*x1;
        */
        for (std::list<PPrioNode>::iterator it: intervall.prioNodesIts) {
            PPrioNode &prioNode = *it;
            //PrioNode &prioNode = prioNodeData.ref;
            //outlier
            double x_o = base_graph.getLon(prioNode.node_id);
            double y_o = base_graph.getLat(prioNode.node_id);
            
            double dist_x;
            double dist_y;
            if (parallelVec.x == 0 && parallelVec.y ==0) { //exception case: divide by zero
                dist_x = x_o - x1;
                dist_y = y_o - y1;
            }else{ //regular case
                double divisor = pow(parallelVec.x, 2)+pow(parallelVec.y, 2);
                debug_assert(divisor != 0);
                //solve linear system for intersection point
                //i.e. (x1, y1) + c * parallelVec = (x_o, y_o) + d * orthogonalVec        
                double c= -((y1-y_o)*parallelVec.y + (x1-x_o)*parallelVec.x)/divisor;
                vec intersection = {x1+c*parallelVec.x, y1+c*parallelVec.y};
                dist_x = x_o - intersection.x;
                dist_y = y_o - intersection.y;
            }
                                    
            //calculate distance between outlier and intersection point                    
            double dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
            
            prioNode.perpendicularLength = dist;
            //prioNode.nOfIntersections = calcNoFIntersections(intervall.start, prioNode.node_id, intervall.end);
            prioNode.cross_diff = calcOrientationMisses(prioNode, intervall);
        }
        return;
    }
        
    //calculates for a line all nodes in the lines' bounding box and also save their orientation with respect to the line
    std::vector<NodeInBox> calcNodesInBoundingBox(NodeID node_id_src, NodeID node_id_tgt) {
        BoundingBox bb = BoundingBox(
                base_graph.getLat(node_id_src),
                base_graph.getLon(node_id_src),
                base_graph.getLat(node_id_tgt),
                base_graph.getLon(node_id_tgt));
        
        std::vector<NodeID> nodesInGrid = grid.nodesInNeigbourhood(bb.midLat, bb.midLon);
        std::vector<NodeInBox> nodesInBox;
        for (NodeID node_id: nodesInGrid) {                    
            //consider only the nodes in the bounding box
            if (bb.contains(base_graph.getLat(node_id), base_graph.getLon(node_id))) {
                //int i1   =       3 > 4 ? 0 : 1;
                //save on which side of the line the node lies
                bool sign = (OrientationTest(node_id_src, node_id_tgt, node_id) >= 0) ? true : false;                
                nodesInBox.push_back(NodeInBox(node_id, sign));
            }            
        }
        return nodesInBox;
    }    
    
    void initializeChain(const Chain &chain) {
            debug_assert(chain.size() >= 3);
            PIntervall initialIntervall(chain.front(), chain.back());
            intervalls.emplace_back(initialIntervall);
            const auto intervallIt = --(intervalls.end());
                        
            //first box of a chain
            std::vector<NodeInBox> nodesInFirstBox(calcNodesInBoundingBox(*(chain.begin()), *(++(chain.begin()))));
            boxesWithNodes.push_back(nodesInFirstBox);            
            
            //initialize lists
            debug_assert(boxesWithNodes.size()>=1);
            auto leftBox = --(boxesWithNodes.end()); //Iterator which runs in tandem with chainNodeIterator                                    
            
            for (auto it = ++(chain.begin()); it != --(chain.end()); it++) {
                
                //compute orientations in bounding box
                auto nextIt = it;
                nextIt++;
                std::vector<NodeInBox> nodesInBox(calcNodesInBoundingBox(*it, *nextIt));
                
                boxesWithNodes.push_back(nodesInBox);
                
                //getting iterators on the left and right boxes of a node
                auto rightBox = leftBox;
                rightBox++;
                
                nodePrios.emplace_back(PPrioNode(*it, intervallIt, *leftBox, *rightBox));                                
                leftBox++;

                //reference each other
                intervallIt->prioNodesIts.emplace_back(--(nodePrios.end()));
                std::list<std::list<PPrioNode>::iterator>::iterator pos = --(intervallIt->prioNodesIts.end());
                nodePrios.back().posInIntervallIt = pos;


            }
            update(intervalls.back());   
    }
    
    public:
        chainPairDP(const GraphT &base_graph, const Grid<GraphT> &grid):
        nodePrios(), intervalls(), base_graph(base_graph), grid(grid), boxesWithNodes() {}
        
        std::list<simplePrioNode> process(const ChainPair &chainpair) {
            debug_assert(chainpair.chainTo.size() + chainpair.chainFrom.size() > 6);
            debug_assert(chainpair.chainTo.size() >=3 && chainpair.chainFrom.size() >= 3);
            nodePrios.clear();            
            intervalls.clear();                                    
            boxesWithNodes.clear();
            
            initializeChain(chainpair.chainTo);
            initializeChain(chainpair.chainFrom);
            
            matchChainPairNodes<GraphT> matcher(base_graph, intervalls.front(), intervalls.back());
            matcher.match();
                        
            return simplify();
        }                     
};

}


