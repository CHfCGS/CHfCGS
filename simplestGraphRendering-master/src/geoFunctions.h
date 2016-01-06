#pragma once

#include <cmath>
#include "nodesAndEdges.h"


//utility functions for geometric calculations (in a plane)
namespace geo {
    
    //static const double PI = 3.14159265;
    //static double R = 6371.009; //Erdradius 
    
    double calcSignedArea(CHNode aSource, CHNode bTarget, CHNode cOutlier) {
        double ax = aSource.lon;
        double ay = aSource.lat;
        double bx = bTarget.lon;
        double by = bTarget.lat;
        double cx = cOutlier.lon;
        double cy = cOutlier.lat;

	double avgLat_ac = (ay + cy)/2;
        double avgLat_bc = (by + cy)/2;
        
        double acx = (ax - cx) * std::cos(avgLat_ac * M_PI / 180.0); // pa[0] - pc[0];
        double bcx = (bx - cx) * std::cos(avgLat_bc * M_PI / 180.0); //pb[0] - pc[0];
        double acy = ay - cy; //pa[1] - pc[1];
        double bcy = by - cy; //pb[1] - pc[1];
	/*
        double acx = ax - cx; // pa[0] - pc[0];
        double bcx = bx - cx; //pb[0] - pc[0];
        double acy = ay - cy; //pa[1] - pc[1];
        double bcy = by - cy; //pb[1] - pc[1];
	*/
        double determinant = acx * bcy - acy * bcx;
        return 0.5 * determinant;
    }
    
    
    /*
    double pythagoras(double a, double b) {
	return sqrt(pow(a, 2) + pow(b, 2));
    }*/
    
    double geoDist(double lon1, double lat1, double lon2, double lat2) {        
        double deltaLon = lon1 - lon2;
        double deltaLat = lat1 - lat2;
        
        double avgLat = (lat1 + lat2)/2;
        
        double x = (deltaLon) * std::cos(avgLat * M_PI / 180.0);
        double y = (deltaLat);
        return sqrt(x*x + y*y); //R only scaling
        
        //return pythagoras(deltaLon, deltaLat);        
    }

    double geoDist(CHNode node1, CHNode node2) {
        double lon1 = node1.lon;
        double lat1 = node1.lat;
        double lon2 = node2.lon;
        double lat2 = node2.lat;
        return geoDist(lon1, lat1, lon2, lat2);
    }  

    double calcPerpendicularLength(CHNode source, CHNode target, CHNode outlier) {
        //calculate perpendicular length as height of a triangle
        double signed_area = calcSignedArea(source, target, outlier);
        double baselength = geoDist(source, target);
        double positiveRectArea = std::fabs(2.0 * signed_area);
        if (baselength == 0) {
            return std::numeric_limits<double>::max();
        } else {
            return positiveRectArea / baselength;        
        }       
    }

    double getTriangleProportion(CHNode source, CHNode target, CHNode outlier) {
	double perpendicularLength = calcPerpendicularLength(source, target, outlier);
        double base_length = geoDist(source, target);
        if (base_length == 0) {
            return std::numeric_limits<double>::max();
        } else {
            double proportion = perpendicularLength/base_length;
            return proportion;
        }   
        return calcPerpendicularLength(source, target, outlier)/geoDist(source, target);
    }
}

