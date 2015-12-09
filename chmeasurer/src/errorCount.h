#pragma once
#include "defs.h"

struct ErrorCounts {
    uint chainNotRedetectedCounter = 0;
    uint chainRedetectedCounter = 0;
    uint chainPairRedetectedCounter = 0;
            
    uint epsilonerrorEqualszero = 0;
    uint epsilonerrorGreaterzero = 0;    
    uint selfIntersecting = 0;
    uint addedDiffs = 0;
    
    uint chain_counter = 0;
    double chain_length_sum = 0;
    double angular_change_sum = 0;
    uint nofNodesWithAngelChange = 0;
        
    double weighted_regularity_sum = 0;
    double weighted_variance_sum = 0;
    
    
    void print () {
        Print("");
        Print("errorCountPrint:");
        Print("chainNotRedetectedCounter:" << chainNotRedetectedCounter);
        Print("chainRedetectedCounter:" << chainRedetectedCounter);
        Print("chainPairRedetectedCounter:" << chainPairRedetectedCounter);
                
        Print("epsilonerrorEqualszero:" << epsilonerrorEqualszero);
        Print("epsilonerrorGreaterzero:" << epsilonerrorGreaterzero);
        Print("selfIntersecting:" << selfIntersecting);
        Print("addedDiffs:" << addedDiffs);
        
        Print("chain_counter:" << chain_counter);
        if (nofNodesWithAngelChange==0) {
            Print("avg_angular_change:" << 0);        
        } else {
            Print("avg_angular_change:" << angular_change_sum/nofNodesWithAngelChange);        
        }
        if (chain_length_sum == 0) {
            Print("regularity_sum:" << "NAN");
            Print("error_consistency_sum:" << "NAN");
        } else {
            Print("avg_regularity:" << weighted_regularity_sum/chain_length_sum);
            Print("avg_error_consistency:" << weighted_variance_sum/chain_length_sum);
        }
        Print("");


        
    }
};