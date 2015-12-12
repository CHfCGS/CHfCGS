#pragma once

#include "simplification/lineSimplifierType.h"
#include "simplification/error_measure.h"
#include "simplification/pairMatchType.h"
#include "dead_end_detector.h"

namespace ls {    
    enum class LineSimplifierType { NONE = 0, BU, DP}; 
//    enum class PairMatchType { NONE = 0, ZO, CD};
}

struct SOptions {                
    
    ls::LineSimplifierType lineSimplifier_type = ls::LineSimplifierType::DP;        
    ls::ErrorMeasureType errorMeasure_type = ls::ErrorMeasureType::VE;
    ls::PairMatchType pairMatch_type = ls::PairMatchType::NONE;  
    DeadEndDetectType deadEndDetect_type = DeadEndDetectType::NONE;
    bool checkBorderCrossing = false;   
    bool taggingViusallyUnpleasantSh = false;
};