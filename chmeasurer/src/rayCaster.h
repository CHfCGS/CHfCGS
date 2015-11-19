#pragma once

#include "nodes_and_edges.h"
#include "grid.h"
#include "geoFunctions.h"
#include "chgraph.h"

namespace chm {

class RayCaster {
    
    CHGraph<CHNode, CHEdge> &graph;
    Grid<CHGraph<CHNode, CHEdge> > &grid;
    const CHNode outerRayPoint = {1000.0, 1000.0, 0, false};   
    /*
    double calcArea(CHNode aSource, CHNode bTarget, CHNode cOutlier) {
        double ax = aSource.lon;
        double ay = aSource.lat;
        double bx = bTarget.lon;
        double by = bTarget.lat;
        double cx = cOutlier.lon;
        double cy = cOutlier.lat;

        double acx = ax - cx; // pa[0] - pc[0];
        double bcx = bx - cx; //pb[0] - pc[0];
        double acy = ay - cy; //pa[1] - pc[1];
        double bcy = by - cy; //pb[1] - pc[1];
        double determinant = acx * bcy - acy * bcx;
        return 0.5 * determinant;
    }

    bool differentSign(double d1, double d2) {
        return (d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0);
    }

    bool testIntersection(Line line1, Line line2) {
        //4 Orientation tests
        double line2StartArea = calcArea(line1.start, line1.end, line2.start);
        double line2EndArea = calcArea(line1.start, line1.end, line2.end);
        bool line1BetweenLine2 = differentSign(line2StartArea, line2EndArea);

        double line1StartArea = calcArea(line2.start, line2.end, line1.start);
        double line1EndArea = calcArea(line2.start, line2.end, line1.end);
        bool line2BetweenLine1 = differentSign(line1StartArea, line1EndArea);

        return line1BetweenLine2 && line2BetweenLine1;
    }
    */
    void subtractLists(std::list<NodeID> &minuend, std::list<NodeID> subtrahend) {
        minuend.sort();
        minuend.unique();
        subtrahend.sort();
        subtrahend.unique();
        
        auto mIt = minuend.begin();
        auto sIt = subtrahend.begin();        
        while (mIt != minuend.end() && sIt != subtrahend.end()) {
            if (*mIt == *sIt) {
                mIt = minuend.erase(mIt);                
                sIt++;
            } else if (*mIt > *sIt) {
                sIt++;
            } else { //(*mIt < *sIt)
                mIt++;
            }
        }        
    }
    
    std::list<NodeID> getOutliers(Chain &chain) {
        
        std::list<NodeID> outliers;
        
        for (NodeID node_id : chain) {
            outliers.splice(outliers.end(), grid.nodeGeoNeighbours(node_id));
        }
        /*
        assert(!chain.empty());
        outliers.splice(outliers.end(), grid.nodeGeoNeighbours(chain.front()));
        outliers.splice(outliers.end(), grid.nodeGeoNeighbours(chain.back()));
         * */
        
        outliers.sort();
        outliers.unique();
        subtractLists(outliers, chain);
        
        return outliers;
    }
    
    bool rayIntersects(NodeID src_id, NodeID tgt_id, NodeID outlier_id) {
        const CHNode src = graph.getNode(src_id);
        const CHNode tgt = graph.getNode(tgt_id);
        CHLine borderline = CHLine(src, tgt);
        const CHNode outlier = graph.getNode(outlier_id);
        CHLine ray = CHLine(outlier, outerRayPoint);
        return geo::testIntersection(borderline, ray);        
    }
    
public:
    RayCaster(CHGraph<CHNode, CHEdge> &graph, Grid<CHGraph<CHNode, CHEdge> > &grid) : graph(graph), grid(grid) {}
    
    uint getNofCrossings(Chain &chain) {
        assert(chain.size() >=3);
        uint nofCrossings = 0;
        std::list<NodeID> outliers = getOutliers(chain);         
        for (NodeID outlier: outliers) {
            uint nofIntersections = 0;            
            //shortcut
            if (rayIntersects(chain.front(), chain.back(), outlier)) {
                nofIntersections++;
            }            
            //former line
            for (auto it = chain.begin(); it != --chain.end(); it++) {                
                auto nextIt = it;
                nextIt++;
                if (rayIntersects(*it, *nextIt, outlier)) {
                    nofIntersections++;
                }
            }
            if (nofIntersections % 2 == 1) {
                nofCrossings++;
            }
        }
        return nofCrossings;        
    }
};

}