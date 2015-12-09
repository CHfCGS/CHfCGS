#pragma once

#include "../geoFunctions.h"
#include "../chgraph.h"
#include "../nodes_and_edges.h"
#include <math.h>

class AngularChange {
    CHGraph<CHNode, CHEdge> &graph;
    
    double calcTurnAngle(Chain::const_iterator cur) {
        Chain::const_iterator prev = cur;
        prev--;
        Chain::const_iterator next = cur;
        next++;
        geo::twoDvector s1(graph.getNode(*prev), graph.getNode(*cur));
        geo::twoDvector s2(graph.getNode(*cur), graph.getNode(*next));
        
        return geo::calcTurnAngle(s1, s2);        
    }
    
public:
    AngularChange(CHGraph<CHNode, CHEdge>& graph) :
        graph(graph) {
    }
            
    double calcACsum(const Chain &redetectedChain, const Chain &expandedChain) {
        //assumption: redetectedChain is subsequence of expandedChain
                
        assert(redetectedChain.size() >= 2);
        assert(expandedChain.size() >= 2);
        assert(expandedChain.front() == redetectedChain.front());
        assert(expandedChain.back() == redetectedChain.back());
        
        uint count = 0;
        double sum = 0;
        
        auto rIt = ++redetectedChain.begin();
        auto eIt = ++expandedChain.begin();        
        while (rIt != --redetectedChain.end() && eIt != --expandedChain.end()) {
            if (*rIt == *eIt) {
                double turnAngleR = calcTurnAngle(rIt);
                double turnAngleE = calcTurnAngle(eIt);
                double delta = turnAngleR - turnAngleE;
                sum += fabs(delta);
                count++;
                rIt++;
                eIt++;
            } else {
                eIt++;
            }
        }
        
        return sum;

        
    }
        
    

};