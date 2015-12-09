#pragma once

#include "nodes_and_edges.h"
#include "grid.h"
#include "geoFunctions.h"
#include "chgraph.h"
#include "CGAL/range_tree.h"
#include "limits"

namespace chm {
    
    namespace unit_tests {
        void testIsBetween();
    }
    

class RayCaster {
    
    CHGraph<CHNode, CHEdge> &graph;
    Grid<CHGraph<CHNode, CHEdge> > &grid;
    RangeTree rangeTree;
    
    const CHNode outerRayPoint;   
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
    
    struct twoDvec {
        double x;
        double y;
        twoDvec getOrthogonal() {
            twoDvec orthogonal;
            orthogonal.x = -y;
            orthogonal.y = x;
            return orthogonal;
        }
    };
    
    twoDvec calculateDirectionVec(CHLine line) {
        twoDvec tdv;
        tdv.x = line.tgt.lon - line.src.lon;
        tdv.y = line.tgt.lat - line.src.lat;
        return tdv;
    }    
    
        
      
    
    
    
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
    
    struct Window {
        double MINLONGITUDE = numeric_limits<double>::max();
        double MAXLONGITUDE = 0;
        double MINLATITUDE = numeric_limits<double>::max();
        double MAXLATITUDE = 0;
        
        
    };
    
    bool isInWindow(Window window, NodeID node_id) {
        CHNode node = graph.getNode(node_id);
        
        if(node.lat > window.MINLATITUDE
                && node.lat < window.MAXLATITUDE
                && node.lat > window.MINLONGITUDE
                && node.lat < window.MAXLONGITUDE) {
            return true;
        } else {
            return false;
        }
        
    }
    
    Window calculateWindow(const Chain &chain) {
        Window window;
        for (NodeID node_id: chain) {
            if (graph.getLon(node_id) > window.MAXLONGITUDE) {
                window.MAXLONGITUDE = graph.getLon(node_id);
            }
            if (graph.getLon(node_id) < window.MINLONGITUDE) {
                window.MINLONGITUDE = graph.getLon(node_id);
            }
            if (graph.getLat(node_id) > window.MAXLATITUDE) {
                window.MAXLATITUDE = graph.getLat(node_id);
            }
            if (graph.getLat(node_id) < window.MINLATITUDE) {
                window.MINLATITUDE = graph.getLat(node_id);
            }
        }
        return window;
    }
    
    std::list<NodeID> getOutliers(Chain &chain) {
        
        Window window = calculateWindow(chain);
        std::list<NodeID> outliers = rangeTree.rectangleQuery(window.MINLATITUDE, window.MINLONGITUDE, window.MAXLATITUDE, window.MAXLONGITUDE);
        //assert(outliers.empty());
        /*
        std::list<NodeID> outliers;
        for (NodeID node_id : chain) {
            outliers.splice(outliers.end(), grid.nodeGeoNeighbours(node_id));
        }
        for (auto it = outliers.begin(); it != outliers.end();) {
            if(isInWindow(window, *it)){
                it = outliers.erase(it);
            }else {
                it++;
            }               
        }
        */
        /*
        assert(!chain.empty());
        outliers.splice(outliers.end(), grid.nodeGeoNeighbours(chain.front()));
        outliers.splice(outliers.end(), grid.nodeGeoNeighbours(chain.back()));
         * */
        
        //outliers.sort();
        //outliers.unique();
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
    uint nofBehinds;
    uint nofIns;
    
    RayCaster(CHGraph<CHNode, CHEdge> &graph, Grid<CHGraph<CHNode, CHEdge> > &grid):
        graph(graph), grid(grid), rangeTree(RangeTree(graph)), outerRayPoint(1000, 1000) {
        nofBehinds = 0;
        nofIns = 0;
        
    }
    
    uint getNofCrossings(Chain &chain) {
        assert(chain.size() >=3);
        uint nofCrossings = 0;
        std::list<NodeID> outliers = getOutliers(chain); 
        CHLine shortcutLine(graph.getNode(chain.front()), graph.getNode(chain.back()));        
        
        /*
        for (auto it = ++chain.begin(); it != --chain.end(); it++) {
                if (!isBetween(shortcutLine, graph.getNode(*it))) {
                    nofBehinds++;                    
                }
                else {
                    nofIns++;
                }
                //Print("nofBehinds" << nofBehinds);
                //Print("nofIns" << nofIns);
            }

        */
        for (NodeID outlier: outliers) {
            
            //isBetween(shortcutLine, graph.getNode(outlier));
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
    
    bool isBetween(const CHLine line, CHNode outlier) {
        //calculate direction vector
        twoDvec directionVector = calculateDirectionVec(line);
                
        //get orthogonal vector
        twoDvec orthogonal = directionVector.getOrthogonal();
        //create imaginary point for src
        CHNode iSrc;
        iSrc.lon = line.src.lon + orthogonal.x;
        iSrc.lat = line.src.lat + orthogonal.y;
        //create line through src and its imaginary point
        CHLine line1(line.src, iSrc);
        
        //create imaginary point for tgt
        CHNode iTgt;
        iTgt.lon = line.tgt.lon + orthogonal.x;
        iTgt.lat = line.tgt.lat + orthogonal.y;
        //create line through tgt and its imaginary point
        CHLine line2(line.tgt, iTgt);
        //test if outlier lies between these lines
        double line1area = geo::calcArea(line1.src, line1.tgt, outlier);
        double line2area = geo::calcArea(line2.src, line2.tgt, outlier);
        //cant happen: (line1area > 0 && line2area < 0)                     
        return (line1area < 0 && line2area > 0);
    }
};

}