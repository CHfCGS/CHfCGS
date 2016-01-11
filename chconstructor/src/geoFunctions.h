#pragma once

#include <cmath>
#include <math.h>
#include "nodes_and_edges.h"

using namespace chc;

//utility functions for geometric calculations (in a plane)
namespace geo {
    
    //static double R = 6371.009; //Erdradius 
    
    double calcSignedArea(OSMNode aSource, OSMNode bTarget, OSMNode cOutlier) {
        
        /*
        double avgLat = (aSource.lat + bTarget.lat + cOutlier.lat) / 3;
        double scale_factor_x = std::cos(avgLat * M_PI / 180.0);
        Print(scale_factor_x);
        double x_1 = aSource.lon * scale_factor_x;
        double y_1 = aSource.lat;
        double x_2 = bTarget.lon * scale_factor_x;
        double y_2 = bTarget.lat;
        double x_3 = cOutlier.lon * scale_factor_x;
        double y_3 = cOutlier.lat;
        
        double sum = (x_3 * y_2 - x_2 * y_3)
                           - (x_3 * y_1 - x_1 * y_3)
                           + (x_2 * y_1 - x_1 * y_2);
        return 0.5 * sum;
        */
        /*
        return 0.5 * determinant = (x_2 * y_3 - x_3 * y_2)
                           - (x_1 * y_3 - x_3 * y_1)
                           - (x_1 * y_2 - x_2 * y_1);
         * */
       
        
        
        double ax = aSource.lon;
        double ay = aSource.lat;
        double bx = bTarget.lon;
        double by = bTarget.lat;
        double cx = cOutlier.lon;
        double cy = cOutlier.lat;
        
        /*
        double avgLat = (aSource.lat + bTarget.lat + cOutlier.lat) / 3;
        double scale_factor_x = std::cos(avgLat * M_PI / 180.0);
        Print(scale_factor_x);
         * */
        
        double avgLat_ac = (ay + cy)/2;
        double avgLat_bc = (by + cy)/2;                
        
        double acx = (ax - cx)  * std::cos(avgLat_ac * M_PI / 180.0); // pa[0] - pc[0];
        double bcx = (bx - cx)  * std::cos(avgLat_bc * M_PI / 180.0); //pb[0] - pc[0];
        double acy = ay - cy; //pa[1] - pc[1];
        double bcy = by - cy; //pb[1] - pc[1];
        

        
        
        //double acx = ax - cx; // pa[0] - pc[0];
        //double bcx = bx - cx; //pb[0] - pc[0];
        //double acy = ay - cy; //pa[1] - pc[1];
        //double bcy = by - cy; //pb[1] - pc[1];
        
        
        double determinant = acx * bcy - acy * bcx;
        return 0.5 * determinant; 
        
       
    }
    
    double calcSignedArea(NodeID src_id, NodeID tgt_id, NodeID outlier_id, const CHGraph<OSMNode, OSMEdge> &graph) {
        return calcSignedArea(graph.getNode(src_id), graph.getNode(tgt_id), graph.getNode(outlier_id));
    }
    
    double pythagoras(double a, double b) {
	return sqrt(pow(a, 2) + pow(b, 2));
    }
    
    double geoDist(double lon1, double lat1, double lon2, double lat2) {                
        double deltaLon = lon1 - lon2;
        double deltaLat = lat1 - lat2;
        
        
        double avgLat = (lat1 + lat2)/2;
       
        double x = (deltaLon) * std::cos(avgLat * M_PI / 180.0);
        double y = (deltaLat);
        return sqrt(x*x + y*y); //R only scaling
        
        //return pythagoras(deltaLon, deltaLat);        
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
        double signed_area = calcSignedArea(source, target, outlier);
        double baselength = geoDist(source, target);
        double positiveRectArea = std::fabs(2.0 * signed_area);
        //return positiveRectArea / baselength;    
        if (baselength == 0) {            
            return std::numeric_limits<double>::max();
        } else {
            return positiveRectArea / baselength;        
        } 
    }
    
    //height divided by baselength
    double getTriangleProportion(OSMNode source, OSMNode target, OSMNode outlier) {
        double perpendicularLength = calcPerpendicularLength(source, target, outlier);
        double base_length = geoDist(source, target);
        if (base_length == 0) {
            return std::numeric_limits<double>::max();
        } else {
            double proportion = perpendicularLength/base_length;
            return proportion;
        }   
        //return calcPerpendicularLength(source, target, outlier)/geoDist(source, target);
    }
    
    //used with calcArea for orientation tests
    bool differentSign(double d1, double d2) {
        return (d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0);
    }

    bool testIntersection(OSMLine line1, OSMLine line2) {
        //4 Orientation tests
        double line2StartArea = calcSignedArea(line1.src, line1.tgt, line2.src);
        double line2EndArea = calcSignedArea(line1.src, line1.tgt, line2.tgt);
        bool line1BetweenLine2 = differentSign(line2StartArea, line2EndArea);

        double line1StartArea = calcSignedArea(line2.src, line2.tgt, line1.src);
        double line1EndArea = calcSignedArea(line2.src, line2.tgt, line1.tgt);
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
        //assert(divisor > 0);
        if (divisor != 0.0) {
            double toAcos = dotProduct(s1, s2)/divisor;
            
            if (toAcos < -1.0) toAcos = -1.0;
            else if (toAcos > 1.0) toAcos = 1.0;

            double turnAngle = acos(toAcos) * 180.0 / M_PI;
            return turnAngle;                        
        } else {
            return 0; //because dotProduct would also be 0
        }
        
    }       
    
}
