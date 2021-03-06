#pragma once

#include "../defs.h"
#include "../nodes_and_edges.h"
#include "../chains.h"
#include "../geoFunctions.h"

#include <memory>

namespace ls {

    /*
     * Interface of a ErrorMeasure used in line simplfication algorithms.
     */
    class ErrorMeasure {
    public:
                   
        virtual double calcError(OSMNode source, OSMNode target, OSMNode outlier) = 0;
        
        //virtual destructor                
        virtual ~ErrorMeasure() {            
        }
    };

}

// class defitions
namespace ls {

//template <class GraphT>
class EpsilonError : public ErrorMeasure{
    double calcError(OSMNode source, OSMNode target, OSMNode outlier) {
        return geo::calcPerpendicularLength(source, target, outlier);    
    }
    
    ~EpsilonError() {}    
};

//template <class GraphT>
class AreaError : public ErrorMeasure{
    double calcError(OSMNode source, OSMNode target, OSMNode outlier) {
        return fabs(geo::calcSignedArea(source, target, outlier)); 
    }    
    ~AreaError() {}
};

class KinkError : public ErrorMeasure{
    double calcError(OSMNode source, OSMNode target, OSMNode outlier) {
        geo::twoDvector s1(source, outlier);
        geo::twoDvector s2(outlier, target);
        double divisor = s1.length + s2.length;
        if (divisor == 0) {
            return 0;
        } else {
            //double turn_angle = geo::calcTurnAngle(s1, s2);
            double turn_angle = geo::calcTurnAngle2(source, target, outlier);
            return (turn_angle * s1.length * s2.length) / divisor;        
        }
        
    } 
    ~KinkError() {}
};

}

//utitility functions
namespace ls {

enum class ErrorMeasureType { NONE = 0, AE, KE, VE };
static constexpr ErrorMeasureType LastErrorMeasureType = ErrorMeasureType::VE;

ErrorMeasureType toErrorMeasureType(std::string const& type)
{
	if (type == "NONE") {
		return ErrorMeasureType::NONE;
	}else if (type == "AE") {
		return ErrorMeasureType::AE;
	}
        else if (type == "KE") {
		return ErrorMeasureType::KE;
	}	
        else if (type == "VE") {
		return ErrorMeasureType::VE;
	}	

	return ErrorMeasureType::NONE;
}

std::string to_string(ErrorMeasureType type)
{
	switch (type) {
	case ErrorMeasureType::NONE:
		return "NONE";
        case ErrorMeasureType::AE:
		return "AE";
        case ErrorMeasureType::KE:
		return "KE";
        case ErrorMeasureType::VE:
		return "VE";	
        }
	std::cerr << "Unknown ErrorMeasureType type: " << static_cast<int>(type) << "\n";
	return "NONE";
        
}

//template <class GraphT>
    std::unique_ptr<ErrorMeasure> createErrorMeasure(ErrorMeasureType errorMeasure_type) {
        switch (errorMeasure_type) {
            case ErrorMeasureType::NONE:
                return nullptr;
            case ErrorMeasureType::AE:
                return std::unique_ptr<ErrorMeasure>(new AreaError());
            case ErrorMeasureType::KE:
                return std::unique_ptr<ErrorMeasure>(new KinkError());
            case ErrorMeasureType::VE:
                return std::unique_ptr<ErrorMeasure>(new EpsilonError());
        }
        return nullptr;
    }
}