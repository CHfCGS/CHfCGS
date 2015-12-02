#pragma once

#include "dp_simplifier.h"
#include "lineSimplifier.h"
#include <memory>

namespace ls {

enum class LineSimplifierType { NONE = 0, DP };
static constexpr LineSimplifierType LastLineSimpliferType = LineSimplifierType::DP;

LineSimplifierType toLineSimpliferType(std::string const& type)
{
	if (type == "NONE") {
		return LineSimplifierType::NONE;
	}
        else if (type == "DP") {
		return LineSimplifierType::DP;
	}	

	return LineSimplifierType::NONE;
}

std::string to_string(LineSimplifierType type)
{
	switch (type) {
	case LineSimplifierType::NONE:
		return "NONE";	
        case LineSimplifierType::DP:
		return "DP";	
        }
	std::cerr << "Unknown LineSimplifier type: " << static_cast<int>(type) << "\n";
	return "NONE";
        
}

template <class GraphT>
    std::unique_ptr<LineSimplifier> createLineSimplifer(LineSimplifierType lineSimpliferType_type, GraphT const& graph, Grid<GraphT> const& grid) {
        switch (lineSimpliferType_type) {
            case LineSimplifierType::NONE:
                return nullptr;
            case LineSimplifierType::DP:
                return std::unique_ptr<LineSimplifier>(new DPSimplifier<GraphT>(graph, grid));
        }
        return nullptr;
    }
}