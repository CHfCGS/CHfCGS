/* 
 * File:   DouglasPeucker.h
 * Author: tobias
 *
 * Created on 4. September 2015, 14:14
 */

#include "nodes_and_edges.h"
#include "grid.h"

#include <limits>
#include <math.h>

#ifndef DOUGLASPEUCKER_H
#define	DOUGLASPEUCKER_H

namespace DP {


template <class GraphT>
class DouglasPeucker{
    private:       
        std::list<PrioNode> nodePrios;
        std::list<Intervall> intervalls;
        //vector<Node> &nodes;
        const GraphT &base_graph;
        const Grid<GraphT> &grid;
        std::vector<std::vector<NodeInBox>> boxesWithNodes;
        
        //naiv TODO: heap implementation; or make_heap
         std::list<PrioNode>::iterator getMax() {// (std::list<PrioNode> &nodePrios) {
            
            assert(!nodePrios.empty());
            nodePrios.sort();
            std::list<PrioNode>::iterator maxIt = --(nodePrios.end());            
            /*
            for (auto it = nodePrios.begin(); it != nodePrios.end(); ++it) {
                Print("\n"<<it->node_id);
                Print("\n"<<it->perpendicularLength);
                Print("\n"<<it->intervallIt);
                Print("\n"<<it->intervallIt)

            }*/
                        
            //nodePrios.erase(maxIt);            
            //PrioNode r(*maxIt);
                                    
            return maxIt;
        }
    
    
                
    std::list<simplePrioNode> split() {        
        if (nodePrios.empty()) {
            return std::list<simplePrioNode>(); //empty list
        }
        else {
            
             std::list<PrioNode>::iterator maxIt = getMax();             
            
            
            //cout << "extracting Node:";
            //cout << nodes[max.node_id].coord.lon << " " <<  nodes[max.node_id].coord.lat << "\n";
            
            std::list<Intervall>::iterator splitIt = maxIt->intervallIt;            
            std::list<PrioNodeRef>::iterator splitPos = maxIt->posInIntervallIt;                                                                            
            
            Intervall leftIntervall(splitIt->start, splitPos->get().node_id);
            Intervall rightIntervall(splitPos->get().node_id, splitIt->end);                        
            
            std::list<Intervall>::iterator leftIt = intervalls.insert(splitIt, leftIntervall);           
            std::list<Intervall>::iterator rightIt = intervalls.insert(splitIt, rightIntervall);                                                                                                
            
            //left Intervall
            leftIt->prioNodesRefs.splice(leftIt->prioNodesRefs.begin(), splitIt->prioNodesRefs, splitIt->prioNodesRefs.begin(), splitPos);
            //leftIt->prioNodes = std::list<PrioNodeRef>(splitIt->prioNodes.begin(), splitPos);
            //cout << "left: \n";            
            for (std::list<PrioNodeRef>::iterator it = leftIt->prioNodesRefs.begin(); it!= leftIt->prioNodesRefs.end(); it++) {
                //cout << nodes[it->get().node_id].coord.lon << " " << nodes[it->get().node_id].coord.lat << "\n";
                it->get().intervallIt = leftIt;
                it->get().posInIntervallIt = it;                
            }

            //right Intervall  
            rightIt->prioNodesRefs.splice(rightIt->prioNodesRefs.begin(), splitIt->prioNodesRefs, ++splitPos, splitIt->prioNodesRefs.end());
            //rightIt->prioNodes = std::list<PrioNodeRef>(++splitPos, splitIt->prioNodes.end());
            //cout << "right: \n";            
            for (std::list<PrioNodeRef>::iterator it = rightIt->prioNodesRefs.begin(); it!= rightIt->prioNodesRefs.end(); it++) {                
                //cout << nodes[it->get().node_id].coord.lon << " " << nodes[it->get().node_id].coord.lat << "\n";
                it->get().intervallIt = rightIt;
                it->get().posInIntervallIt = it;                 
            }                                    
            
            intervalls.erase(splitIt);
                        
            update (*leftIt);
            update (*rightIt);
            
            simplePrioNode spn(*maxIt);
            nodePrios.erase(maxIt);
            //cout << "\n";
            
            std::list<simplePrioNode> r = split();//(max.intervallIt); //, max.posInIntervall);                    
            r.push_back(spn);            
            return r;
        }           
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
    uint calcOrientationMisses(const PrioNode &split, const Intervall &intervall){
        uint misscount = 0;
        uint assertcounter = 0;
        auto splitPosIntervallincr= split.posInIntervallIt;
        splitPosIntervallincr++;
        for (auto it = intervall.prioNodesRefs.begin(); it != splitPosIntervallincr; it++) {
            for(NodeInBox nodeInBox : it->get().leftBox) {
                bool sign = (OrientationTest(intervall.start, split.node_id, nodeInBox.nodeID) >= 0) ? true : false;  
                if (sign != nodeInBox.sign) {
                    misscount++;
                }
            }                  
            assertcounter++;
        }
        for (auto it = split.posInIntervallIt; it != intervall.prioNodesRefs.end(); it++) {
            for(NodeInBox nodeInBox : it->get().rightBox) {
                bool sign = (OrientationTest(split.node_id, intervall.end, nodeInBox.nodeID) >= 0) ? true : false;  
                if (sign != nodeInBox.sign) {
                    misscount++;
                }
            }                                
            assertcounter++;
        }        
        debug_assert(assertcounter == intervall.prioNodesRefs.size()+1);
        return misscount;
    }
    
    /*
     * updates the priority data for all prionodes in an Intervall
     */    
    void update(Intervall &intervall){
        
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
        for (PrioNode &prioNode: intervall.prioNodesRefs) {
          
            //PrioNode &prioNode = prioNodeData.ref;
            //outlier
            double x_o = base_graph.getLon(prioNode.node_id);
            double y_o = base_graph.getLat(prioNode.node_id);
            
            double dist_x;
            double dist_y;
            if (parallelVec.x == 0 && parallelVec.y == 0) { //exception case: divide by zero
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
            prioNode.nOfIntersections = calcOrientationMisses(prioNode, intervall);
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
    
    public:        
        DouglasPeucker(const GraphT &base_graph, const Grid<GraphT> &grid):
        nodePrios(), intervalls(), base_graph(base_graph), grid(grid), boxesWithNodes() {}
        
        std::list<simplePrioNode> process(const Chain &chain) {                

            debug_assert(chain.size() >= 3);

            nodePrios.clear();
            intervalls.clear();
            boxesWithNodes.clear();
                    
            Intervall initialIntervall(chain.front(), chain.back());
            intervalls.emplace_back(initialIntervall);
            const auto intervallIt = intervalls.begin();

            //precompute orientations in bounding boxes
            for (auto it = chain.begin(); it != --(chain.end()); it++) {
                auto nextIt = it;
                nextIt++;                
                std::vector<NodeInBox> nodesInBox(calcNodesInBoundingBox(*it, *nextIt));                
                boxesWithNodes.push_back(nodesInBox);                
            }
            
            //initialize lists
            debug_assert(chain.size() == boxesWithNodes.size() +1);
            debug_assert(boxesWithNodes.size() >= 1);
            auto leftBox = boxesWithNodes.begin(); //Iterator which runs in tandem with chainNodeIterator
            for (auto it = ++(chain.begin()); it != --(chain.end()); it++) {
                //getting iterators on the left and right boxes of a node
                auto rightBox = leftBox;
                rightBox++;
                
                nodePrios.emplace_back(PrioNode(*it, intervallIt, *leftBox, *rightBox));
                leftBox++;

                //reference each other
                intervallIt->prioNodesRefs.emplace_back(nodePrios.back());
                std::list<PrioNodeRef>::iterator pos = --(intervallIt->prioNodesRefs.end());
                nodePrios.back().posInIntervallIt = pos;


            }
            update(intervalls.front());     
            
            return split();
        }                     
};

}

#endif	/* DOUGLASPEUCKER_H */

