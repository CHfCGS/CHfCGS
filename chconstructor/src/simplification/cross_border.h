#pragma once

#include "../chgraph.h"
#include "../grid.h"
#include "prio_nodes.h"
#include "../geoFunctions.h"
#include "../bounding_box.h"

namespace ls {
    namespace cb {

        /*
        inline double OrientationTest(NodeID aSource, NodeID bTarget, NodeID cOutlier) {
            double ax = graph.getLon(aSource);
            double ay = graph.getLat(aSource);
            double bx = graph.getLon(bTarget);
            double by = graph.getLat(bTarget);
            double cx = graph.getLon(cOutlier);
            double cy = graph.getLat(cOutlier);

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

        
        //calculates how many nodes are correctly oriented after the intervall
        //is spanned by the lines from start to split and from split to end
        
        
        //calculates for a line all nodes in the lines' bounding box and also save their orientation with respect to the line
        std::vector<NodeInBox> calcNodesInBoundingBox(NodeID node_id_src, NodeID node_id_tgt) {
            BoundingBox bb = BoundingBox(
                    graph.getLat(node_id_src),
                    graph.getLon(node_id_src),
                    graph.getLat(node_id_tgt),
                    graph.getLon(node_id_tgt));

            std::vector<NodeID> nodesInGrid = grid.nodesInNeigbourhood(bb.midLat, bb.midLon);
            std::vector<NodeInBox> nodesInBox;
            for (NodeID node_id: nodesInGrid) {                    
                //consider only the nodes in the bounding box
                if (bb.contains(graph.getLat(node_id), graph.getLon(node_id))) {
                    //int i1   =       3 > 4 ? 0 : 1;
                    //save on which side of the line the node lies
                    bool sign = (OrientationTest(node_id_src, node_id_tgt, node_id) >= 0) ? true : false;                
                    nodesInBox.push_back(NodeInBox(node_id, sign));
                }            
            }
            return nodesInBox;
        } 
        */
        /*
        inline bool intersect(NodeID a, NodeID b, NodeID c, NodeID d) {
            double sideOfc = geo::calcArea(a, b, c);
            double sideOfd = geo::calcArea(a, b, d);
            return (sideOfc * sideOfd) < 0;
        }*/
        
        //calculates for a line all nodes in the lines' bounding box and also save their orientation with respect to the line
        NodeBox calcNodesInBoundingBox(NodeID node_id_src, NodeID node_id_tgt,
                                        const CHGraph<OSMNode, OSMEdge> &graph,
                                        const Grid<CHGraph<OSMNode, OSMEdge> > &grid) {
            const OSMNode src = graph.getNode(node_id_src);
            const OSMNode tgt = graph.getNode(node_id_tgt);            
            BoundingBox bb = BoundingBox(
                    src.lat,
                    src.lon,
                    tgt.lat,
                    tgt.lon);

            std::list<NodeID> nodesInGrid = grid.nodesInBoundingBox(bb);
            std::vector<NodeInBox> nodesInBox;
            for (NodeID node_id: nodesInGrid) {                  
                //if(graph.isActive(node_id)) {
                    bool sign = (geo::calcSignedArea(src, tgt, graph.getNode(node_id)) >= 0) ? true : false;                
                    nodesInBox.push_back(NodeInBox(node_id, sign));
                //}                
            }
            /*
            std::vector<NodeID> nodesInGrid = grid.nodesInNeigbourhood(bb.midLat, bb.midLon);
            std::vector<NodeInBox> nodesInBox;
            for (NodeID node_id: nodesInGrid) {                    
                //consider only the nodes in the bounding box
                if (bb.contains(graph.getLat(node_id), graph.getLon(node_id))) {                    
                    //save on which side of the line the node lies
                    bool sign = (geo::calcArea(src, tgt, graph.getNode(node_id)) >= 0) ? true : false;                
                    nodesInBox.push_back(NodeInBox(node_id, sign));
                }            
            }
             * */
            return nodesInBox;
        }                       
    }
}