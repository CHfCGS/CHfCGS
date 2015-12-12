#pragma once

//#include "../s_options.h"

//struct SOptions;

namespace ls {

//    namespace ls {
    //enum class LineSimplifierType { NONE = 0, DP };
//}
    
enum class PairMatchType { NONE = 0, ZO, CD, P};

static constexpr PairMatchType LastPairMatchType = PairMatchType::P;

PairMatchType toPairMatchType(std::string const& type)
{
	if (type == "NONE") {
            return PairMatchType::NONE;
	}
        else if (type == "ZO") {
            return PairMatchType::ZO;
	}
        else if (type == "CD") {
            return PairMatchType::CD;
	}
        else if (type == "P") {
            return PairMatchType::P;
	}

	return PairMatchType::NONE;
}

std::string to_string(PairMatchType type)
{
	switch (type) {
	case PairMatchType::NONE:
		return "NONE";	
        case PairMatchType::ZO:
		return "ZO";	
        case PairMatchType::CD:
		return "CD";
        case PairMatchType::P:
		return "P";
        }
	std::cerr << "Unknown PairMatchType type: " << static_cast<int>(type) << "\n";
	return "NONE";
        
}

}