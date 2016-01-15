#pragma once
#include "defs.h"

#include <iomanip> 
#include <sstream>

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


struct OutData {
        /*    
        & t
        & $\frac{|E_{n_z, \epsilon_z}|}{|E'_i|} $
        & $\sharp gw$
        & $\Delta n_{\Pslash}$
        & $\% \times_{\Pslash}$   
        & $\Delta n_{\Pslash_t}$
        & $\% \times_{\Pslash_t}$
        & $\delta_{dF}({\Pslash_t})$
        & $\delta_{dF}(P_{t})$ \\ \hline
        % |E'_i| |
        %$dist(\Plslash)$
        %| $dist(\Plslash_t)$ 
        */ 
        double dijkstra_time;
        
        double fraction;
        
        uint crossings;
        
        double deltapslash;
        double intersection_plash;
        
        double deltapslasht;
        double intersection_plasht;
        double etapslasht;
        double etapt;
        
        double E_i_size;
        double dist_Plslash;
        double dist_Plslash_t;
        double E_plus_size;
        
        /*
        double pi = 3.14159265359;
        stringstream stream;
        stream << fixed << setprecision(2) << pi;
        string s = stream.str();
        */
        
                
        
        
        
        void print() {
            const uint nofVariables = 13;
            std::vector<std::stringstream> streams(nofVariables);
            //std::vector<std::string> strings(nofVariables);
            
            streams[0] << std::fixed << std::setprecision(0) << dijkstra_time;
            
            streams[1] << std::fixed << std::setprecision(2) << fraction;
            
            streams[2] << crossings;
            
            streams[3] << std::fixed << std::setprecision(3) << deltapslash;
            streams[4] << std::fixed << std::setprecision(3) << intersection_plash * 100;
            
            streams[5] << std::fixed << std::setprecision(3) << deltapslasht;
            streams[6] << std::fixed << std::setprecision(3) << intersection_plasht * 100;
            streams[7] << std::fixed << std::setprecision(3) << etapslasht;            
            streams[8] << std::fixed << std::setprecision(3) << etapt;
            
            streams[9] << E_i_size;            
            streams[10] << dist_Plslash;
            streams[11] << dist_Plslash_t;  
            streams[12] << E_plus_size; 
            
            
            Print("--------------------");
            Print("Outdata:");
            for (uint i = 0; i < 9; i++) {                
                std::cout << " & " << streams[i].str() ;
            }
            std::cout << "\\\\ %";
            for (uint i = 9; i < nofVariables; i++) {                
                std::cout << " | " << streams[i].str() ;
            }                                    
            std::cout << std::endl;
                                    
                    
                    
            //streams.resize(nofVariables);
            //streams.resize(nofVariables);
            
            
            
            
            
            //for (uint i = 0, i)
            
            //std::cout << x << std::endl;
            
            
            
            Print("--------------------");
        }
        
    };