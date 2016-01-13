#pragma once
#include "defs.h"

struct ErrorCounts {
    uint chainNotRedetectedCounter = 0;
    uint chainRedetectedCounter = 0;
    uint chainPairRedetectedCounter = 0;
            
    
    
    uint nofCrossings = 0;
    
    double weighted_epsilon = 0;
    
    uint epsilonerrorEqualszero = 0;
    uint epsilonerrorGreaterzero = 0;    
    uint nofAnalyzedPolygons = 0;
    uint selfIntersecting = 0;
    uint sameLocation = 0;
    double addedDiffs = 0; 
    double weighted_eta = 0;
    double weighted_eta2 = 0;
    
    uint chain_counter = 0;
    uint node_length_sum = 0;
    double length_sum = 0;
    double length_sum2 = 0;
    double angular_change_sum = 0;
    uint nofNodesWithAngelChange = 0;
        
    double weighted_regularity_sum = 0;
    double weighted_variance_sum = 0;
    
    
    void print () {
        Print("");
        Print("errorCountPrint:");
        /*
        Print("chainNotRedetectedCounter:" << chainNotRedetectedCounter);
        Print("chainRedetectedCounter:" << chainRedetectedCounter);
        Print("chainPairRedetectedCounter:" << chainPairRedetectedCounter);
        
        Print("epsilonerrorEqualszero:" << epsilonerrorEqualszero);
        Print("epsilonerrorGreaterzero:" << epsilonerrorGreaterzero);
        */
                
        Print("length_sum" << length_sum);
        
        Print("sameLocation:" << sameLocation);
        Print("selfIntersecting:" << selfIntersecting);
        Print("addedDiffs:" << addedDiffs);
        
        /*
        Print("chain_counter:" << chain_counter);
        if (nofNodesWithAngelChange==0) {
            Print("avg_angular_change:" << 0);        
        } else {
            Print("avg_angular_change:" << angular_change_sum/nofNodesWithAngelChange);        
        }
        if (length_sum == 0) {
            Print("regularity_sum:" << "NAN");
            Print("error_consistency_sum:" << "NAN");
        } else {
            Print("avg_regularity:" << weighted_regularity_sum/length_sum);
            Print("avg_error_consistency:" << weighted_variance_sum/length_sum);
        }
         * */
        Print("");


        
    }
};