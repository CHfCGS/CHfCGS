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
        return geo::calcArea(source, target, outlier); 
    }    
    ~AreaError() {}
};

class KinkError : public ErrorMeasure{
    double calcError(OSMNode source, OSMNode target, OSMNode outlier) {
        geo::twoDvector s1(source, outlier);
        geo::twoDvector s2(outlier, target);
        return (calcTurnAngle(s1, s2) * s1.length * s2.length) / (s1.length + s2.length);        
    } 
    ~KinkError() {}
};

}

//utitility functions
namespace ls {

enum class ErrorMeasureType { NONE = 0, AE, KE, EF };
static constexpr ErrorMeasureType LastErrorMeasureType = ErrorMeasureType::EF;

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
        else if (type == "DP") {
		return ErrorMeasureType::EF;
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
		return "AE";
        case ErrorMeasureType::EF:
		return "DP";	
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
            case ErrorMeasureType::EF:
                return std::unique_ptr<ErrorMeasure>(new EpsilonError());
        }
        return nullptr;
    }
}