#pragma once

#include <cmath>
#include "nodes_and_edges.h"

using namespace chm;

//utility functions for geometric calculations
namespace geo
{

    //static const double PI = 3.14159265;
    static const double R = 6371; //Erdradius in kilometern


    //flag to calculate with equirectangular approximation
    //http://www.movable-type.co.uk/scripts/latlong.html
    //simpler is euclidean which is also used by unit tests
    static const bool equirectangular = false;

    /*
    double toRadians (double degree) {
        return degree * M_PI / 180.0;
    } */

    double square(double x)
    {
        return x * x;
    }

    double calcSignedArea(CHNode aSource, CHNode bTarget, CHNode cOutlier)
    {
        double ax = aSource.lon;
        double ay = aSource.lat;
        double bx = bTarget.lon;
        double by = bTarget.lat;
        double cx = cOutlier.lon;
        double cy = cOutlier.lat;

        double acx = (ax - cx);
        double bcx = (bx - cx);

        if (equirectangular)
        {
            double avgLat_ac = (ay + cy) / 2;
            double avgLat_bc = (by + cy) / 2;

            acx *= std::cos(avgLat_ac * M_PI / 180.0);
            bcx *= std::cos(avgLat_bc * M_PI / 180.0);
        }

        double acy = ay - cy;
        double bcy = by - cy;
        double determinant = acx * bcy - acy * bcx;

        double signedArea = 0.5 * determinant;
        if (equirectangular)
        {
            signedArea *= square(R) * square((M_PI / 180.0));
        }
        return signedArea;
    }

    double pythagoras(double a, double b)
    {
        return sqrt(square(a) + square(b));
    }

    double geoDist(double lon1, double lat1, double lon2, double lat2)
    {
        double deltaLon = lon1 - lon2;
        double deltaLat = lat1 - lat2;

        if (equirectangular)
        {
            double avgLat = (lat1 + lat2) / 2;

            double x = (deltaLon) * std::cos(avgLat * M_PI / 180.0);
            double y = (deltaLat);
            return pythagoras(x, y) * R * (M_PI / 180.0); //R only scaling
        } else
        {
            return pythagoras(deltaLon, deltaLat);
        }
    }

    double geoDist(CHNode node1, CHNode node2)
    {
        double lon1 = node1.lon;
        double lat1 = node1.lat;
        double lon2 = node2.lon;
        double lat2 = node2.lat;
        return geoDist(lon1, lat1, lon2, lat2);
    }

    double geoDist(const Line line)
    {
        return geoDist(line.start.ch_node, line.end.ch_node);
    }

    double calcPerpendicularLength(const CHNode source, const CHNode target, const CHNode outlier)
    {
        //calculate perpendicular length as height of a triangle
        double signed_area = calcSignedArea(source, target, outlier);
        double baselength = geoDist(source, target);
        double positiveRectArea = std::fabs(2.0 * signed_area);
        if (baselength == 0)
        {
            return std::numeric_limits<double>::max();
        } else
        {
            return positiveRectArea / baselength;
        }
    }

    double calcPerpendicularLength(const Line line, const CHNode outlier)
    {
        return calcPerpendicularLength(line.start.ch_node, line.end.ch_node, outlier);
    }

    //used with calcArea for orientation tests

    bool differentSign(double d1, double d2)
    {
        return (d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0);
    }

    //height divided by baselength

    double getTriangleProportion(CHNode source, CHNode target, CHNode outlier)
    {
        double perpendicularLength = calcPerpendicularLength(source, target, outlier);
        double base_length = geoDist(source, target);
        if (base_length == 0)
        {
            return std::numeric_limits<double>::max();
        } else
        {
            double proportion = perpendicularLength / base_length;
            return proportion;
        }
    }

    bool testIntersection(CHLine line1, CHLine line2)
    {
        //4 Orientation tests
        double line2StartArea = calcSignedArea(line1.src, line1.tgt, line2.src);
        double line2EndArea = calcSignedArea(line1.src, line1.tgt, line2.tgt);
        bool line1BetweenLine2 = differentSign(line2StartArea, line2EndArea);

        double line1StartArea = calcSignedArea(line2.src, line2.tgt, line1.src);
        double line1EndArea = calcSignedArea(line2.src, line2.tgt, line1.tgt);
        bool line2BetweenLine1 = differentSign(line1StartArea, line1EndArea);

        return line1BetweenLine2 && line2BetweenLine1;
    }

    struct twoDvector
    {
        double x;
        double y;
        double length;

        twoDvector(CHNode src, CHNode tgt)
        {
            x = tgt.lon - src.lon;
            y = tgt.lat - src.lat;
            length = geo::geoDist(src, tgt);
        }
    };

    double dotProduct(twoDvector s1, twoDvector s2)
    {
        double xProduct = s1.x * s2.x;
        double yProduct = s1.y * s2.y;
        return xProduct + yProduct;
    }

    double calcTurnAngle(twoDvector s1, twoDvector s2)
    {
        //double length1 = geo::geoDist(graph.getNode(kink.src), graph.getNode(kink.peak));
        //double length2 = geo::geoDist(graph.getNode(kink.peak), graph.getNode(kink.tgt));
        double divisor = s1.length * s2.length;
        assert(divisor > 0);
        double toAcos = dotProduct(s1, s2) / divisor;

        //safe for acos
        if (toAcos < -1.0) toAcos = -1.0;
        else if (toAcos > 1.0) toAcos = 1.0;

        double turnAngle = acos(toAcos) * 180.0 / M_PI;
        return turnAngle;
    }

    struct twoDvec
    {
        double x;
        double y;

        twoDvec getOrthogonal()
        {
            twoDvec orthogonal;
            orthogonal.x = -y;
            orthogonal.y = x;
            return orthogonal;
        }
    };

    twoDvec calculateDirectionVec(CHLine line)
    {
        twoDvec tdv;
        tdv.x = line.tgt.lon - line.src.lon;
        tdv.y = line.tgt.lat - line.src.lat;
        return tdv;
    }

    //Deprecated: Will not work in non euclidean space probably

    bool isBetweenByVector(const CHLine line, CHNode outlier)
    {
        assert(!equirectangular);

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
        double line1area = geo::calcSignedArea(line1.src, line1.tgt, outlier);
        double line2area = geo::calcSignedArea(line2.src, line2.tgt, outlier);
        //cant happen: (line1area > 0 && line2area < 0)
        return (line1area < 0 && line2area > 0);
    }

    //Deprecated: Wont work for FrechetDistance test: TODO find out why

    bool isBetweenByAngles(const CHLine line, CHNode outlier)
    {
        double perpendicularLength = calcPerpendicularLength(line.src, line.tgt, outlier);

        //The perpendicular on the line is dividing the triangle in two right-angled triangles
        //We now consider the angles alpha and beta at the corners of our given line

        //The opposite side is the same for both angles
        double distOppositeSide = perpendicularLength;

        double alpha = asin(distOppositeSide / geoDist(line.src, outlier));
        double beta = asin(distOppositeSide / geoDist(line.tgt, outlier));

        return alpha < 90 && beta < 90;
    }

    bool isBetweenByThales(const CHLine line, CHNode outlier)
    {
        double a = geoDist(line.tgt, outlier);
        double b = geoDist(line.src, outlier);
        double c = geoDist(line.src, line.tgt);

        //alpha belongs to src
        bool alphaIsAcute = square(a) < square(b) + square(c);

        //beta belongs to tgt
        bool betaIsAcute = square(b) < square(a) + square(c);

        return alphaIsAcute && betaIsAcute;
    }
}
