#pragma once

struct BoundingBox {
    double upperLat;
    double LeftLon;   
    double lowerLat;
    double RightLon;
    double midLat;
    double midLon;
    
    BoundingBox(double lat1, double lon1, double lat2, double lon2) {
        if (lat1 > lat2) {
            upperLat = lat1;
            lowerLat = lat2;
        } else {
            upperLat = lat2;
            lowerLat = lat1;
        }
        if (lon1 > lon2) {
            RightLon = lon1;
            LeftLon = lon2;
        } else {
            RightLon = lon2;
            LeftLon = lon1;
        }
        midLat = (upperLat + lowerLat) /2;
        midLon = (LeftLon + RightLon) /2;
    }
    
    bool contains(double lat, double lon) {
        return ((lowerLat < lat) && (lat < upperLat)
                && (LeftLon < lon) && (lon < RightLon));
    }
    
};