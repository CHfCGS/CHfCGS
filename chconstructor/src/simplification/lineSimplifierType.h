#pragma once

#include "nth_point_simplifier.h"
#include "dp_simplifier.h"
#include "lineSimplifier.h"
#include "../s_options.h"
#include "discreteCurveEvolution.h"


#include <memory>

struct SOptions;

namespace ls {

//    namespace ls {
    //enum class LineSimplifierType { NONE = 0, DP };
//}


static constexpr LineSimplifierType LastLineSimpliferType = LineSimplifierType::BU;

LineSimplifierType toLineSimplifierType(std::string const& type)
{
	if (type == "NONE") {
            return LineSimplifierType::NONE;
	}
        else if (type == "NTH") {
            return LineSimplifierType::NTH;
	}
        else if (type == "DP") {
            return LineSimplifierType::DP;
	}
        else if (type == "BU") {
            return LineSimplifierType::BU;
	}	

	return LineSimplifierType::NONE;
}

std::string to_string(LineSimplifierType type)
{
	switch (type) {
	case LineSimplifierType::NONE:
		return "NONE";	
        case LineSimplifierType::NTH:
		return "NTH";	
        case LineSimplifierType::DP:
		return "DP";	
        case LineSimplifierType::BU:
		return "BU";	
        }
	std::cerr << "Unknown LineSimplifier type: " << static_cast<int>(type) << "\n";
	return "NONE";
        
}



template <class GraphT>
    std::unique_ptr<LineSimplifier> createLineSimplifier(SOptions s_options, LineSimplifierType lineSimplifier_type, GraphT const& graph, Grid<GraphT> const& grid) {        
        switch (lineSimplifier_type) {
            case LineSimplifierType::NONE:
                return nullptr;
            case LineSimplifierType::NTH:                
                return std::unique_ptr<LineSimplifier>(new NthPointSimplifier<GraphT>());
            case LineSimplifierType::DP:
                return std::unique_ptr<LineSimplifier>(new DPSimplifier<GraphT>(s_options, graph, grid));
            case LineSimplifierType::BU:                
                return std::unique_ptr<LineSimplifier>(new DiscreteCurveEvolution<GraphT>(s_options, graph, grid));
        }
        return nullptr;
    }
}