#pragma once

#include <cmath>
#include "nodes_and_edges.h"

using namespace chc;

//utility functions for geometric calculations (in a plane)
namespace geo {
    
    //static double R = 6371.009; //Erdradius 
    
    double calcArea(OSMNode aSource, OSMNode bTarget, OSMNode cOutlier) {
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
    
    double calcArea(NodeID src_id, NodeID tgt_id, NodeID outlier_id, const CHGraph<OSMNode, OSMEdge> &graph) {
        return calcArea(graph.getNode(src_id), graph.getNode(tgt_id), graph.getNode(outlier_id));
    }
    
    double pythagoras(double a, double b) {
	return sqrt(pow(a, 2) + pow(b, 2));
    }
    
    double geoDist(double lon1, double lat1, double lon2, double lat2) {        
        double deltaLon = lon1 - lon2;
        double deltaLat = lat1 - lat2;
        return pythagoras(deltaLon, deltaLat);        
    }

    double geoDist(OSMNode node1, OSMNode node2) {
        double lon1 = node1.lon;
        double lat1 = node1.lat;
        double lon2 = node2.lon;
        double lat2 = node2.lat;
        return geoDist(lon1, lat1, lon2, lat2);
    }

    double calcPerpendicularLength(OSMNode source, OSMNode target, OSMNode outlier) {
        //calculate perpendicular length as height of a triangle
        double area = calcArea(source, target, outlier);
        double baselength = geoDist(source, target);
        double positiveRectArea = std::abs(2.0 * area);
        return positiveRectArea / baselength;        
    }
    
    //used with calcArea for orientation tests
    bool differentSign(double d1, double d2) {
        return (d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0);
    }

    bool testIntersection(OSMLine line1, OSMLine line2) {
        //4 Orientation tests
        double line2StartArea = calcArea(line1.src, line1.tgt, line2.src);
        double line2EndArea = calcArea(line1.src, line1.tgt, line2.tgt);
        bool line1BetweenLine2 = differentSign(line2StartArea, line2EndArea);

        double line1StartArea = calcArea(line2.src, line2.tgt, line1.src);
        double line1EndArea = calcArea(line2.src, line2.tgt, line1.tgt);
        bool line2BetweenLine1 = differentSign(line1StartArea, line1EndArea);

        return line1BetweenLine2 && line2BetweenLine1;
    }  
    /*
    //originally from curve evolution
    struct Kink {
        const OSMNode src;
        const OSMNode peak;
        const OSMNode tgt;
        Kink(OSMNode src, OSMNode peak, OSMNode tgt): src(src), peak(peak), tgt(tgt) {}
    };*/
    
    struct twoDvector {
        double x;
        double y;
        double length;
        twoDvector(OSMNode src, OSMNode tgt) {
            x = tgt.lon - src.lon;
            y = tgt.lat - src.lat;
            length = geo::geoDist(src, tgt);
        }        
    };
    
    double dotProduct (twoDvector s1, twoDvector s2) {        
        double xProduct = s1.x * s2.x;
        double yProduct = s1.y * s2.y;
        return xProduct + yProduct;
    }
    
    double calcTurnAngle(twoDvector s1, twoDvector s2){
        //double length1 = geo::geoDist(graph.getNode(kink.src), graph.getNode(kink.peak));
        //double length2 = geo::geoDist(graph.getNode(kink.peak), graph.getNode(kink.tgt));        
        double divisor = s1.length * s2.length;   
        assert(divisor > 0);
        double toAcos = fabs(dotProduct(s1, s2))/divisor;
        
        return acos(toAcos);        
    }       
    
}
