#pragma once

#include <math.h>

#include "../nodes_and_edges.h"
#include "../geoFunctions.h"
#include "prio_nodes.h"

namespace ls {


//Bottom Up generalisation

template<GraphT>
class DiscreteCurveEvolution {
    
    PrioNodeHeap prioNodeHeap;
    std::list<Intervall2> intervalls;
    std::list<std::vector<NodeInBox> > boxesWithNodes;
    
    GraphT &graph;    
    

    struct Kink {
        const NodeID src;
        const NodeID peak;
        const NodeID tgt;
        Kink(NodeID src, NodeID peak, NodeID tgt): src(src), peak(peak), tgt(tgt) {}
    };
    
    struct twoDvector {
        double x;
        double y;
        double length;
        twoDvector(NodeID src, NodeID tgt) {
            x = graph.getNode(tgt).lon - graph.getNode(src).lon;
            y = graph.getNode(tgt).lat - graph.getNode(src).lat;
            length = geo::geoDist(graph.getNode(src), graph.getNode(tgt));
        }        
    };
    
    double dotProduct (twoDvector s1, twoDvector s2) {        
        double xProduct = s1.x * s2.x;
        double yProduct = s1. * s2.y;
        return xProduct + yProduct;
    }
    
    double calcTurnAngle(twoDvector s1, twoDvector s2){
        //double length1 = geo::geoDist(graph.getNode(kink.src), graph.getNode(kink.peak));
        //double length2 = geo::geoDist(graph.getNode(kink.peak), graph.getNode(kink.tgt));        
        double toAcos = abs(dotProduct(s1, s2))/(s1.length * s2.length);
        return acos(toAcos);        
    }
    
    double calcRelevanceMeasure(const Kink kink) {
        twoDvector s1(kink.src, kink.peak);
        twoDvector s2(kink.peak, kink.tgt);
        return (calcTurnAngle(s1, s2) * s1.length * s2.length) / (s1.length + s2.length);
    }
            
    void setIntervall(Intervall2 &intervall){
        if(intervall.prioNodeHandles.empty()) {
            return;
        } else if(intervall.prioNodeHandles.size() == 1) {
            PrioNodeHandle prioNodeH = *intervall.prioNodeHandles.begin();
            PrioNode2 &prioNode = *prioNodeH;
            Kink kink(intervall.start, prioNode.node_id, intervall.finish);
            prioNode.perpendicularLength = calcRelevanceMeasure(kink);
            prioNodeHeap.update(prioNodeH);
        } else {
            PrioNodeHandle prev_h;
            PrioNodeHandle cur_h;
            PrioNodeHandle next_h;
            //front
            cur_h = *intervall.prioNodeHandles.begin();
            PrioNode2 &cur = *cur_h;
            next_h = *++(intervall.prioNodeHandles.begin());
            PrioNode2 &next = *next_h;
            Kink kink(intervall.start, cur.node_id, next.node_id);
            cur.perpendicularLength = calcRelevanceMeasure(kink);
            prioNodeHeap.update(cur);
            //middle
            for(std::list<PrioNodeHandle>::iterator it = ++intervall.prioNodeHandles.begin();
                it != --intervall.prioNodeHandles.end();
                it++) {
                auto lastIt = it;
                lastIt--;
                auto nextIt = it;
                nextIt++;
                
                prev_h = *lastIt;
                PrioNode2 &prev = *prev_h;
                cur_h = *it;
                PrioNode2 &cur = *cur_h;
                next_h = *nextIt;
                PrioNode2 &next = *next_h;
                
                Kink kink(prev.node_id, cur.node_id, next.node_id);
                cur.perpendicularLength = calcRelevanceMeasure(kink);
                prioNodeHeap.update(cur);                
            }
            //back            
            prev_h = *--(intervall.prioNodeHandles.end());
            PrioNode2 &prev = *prev_h;
            cur_h = *intervall.prioNodeHandles.end();
            PrioNode2 &cur = *cur_h;
            Kink kink(prev.node_id, cur.node_id, intervall.finish);
            cur.perpendicularLength = calcRelevanceMeasure(kink);
            prioNodeHeap.update(cur);
        }        
        
        for (PrioNodeHandle prioNodeH: intervall.prioNodeHandles) {
            PrioNode2 &prioNode = *prioNodeH;  
            
            prioNode.perpendicularLength = geo::calcRelevanceMeasure(graph.getNode(intervall.start),
                                                                        graph.getNode(intervall.finish),
                                                                        graph.getNode(prioNode.node_id));
            
            prioNode.nOfIntersections = calcOrientationMisses(prioNode, intervall);
            prioNodeHeap.update(prioNodeH);
        }
        return;
    }
    
    void initializeChain(const Chain &chain) {                   
        
            debug_assert(chain.size() >= 3);            
            
            Intervall2 initialIntervall(chain.front(), chain.back());
            intervalls.emplace_back(initialIntervall);
            const auto intervallIt = --(intervalls.end());
            
            //empty boxes for now
            std::vector<NodeInBox> nodesInFirstBox();
            boxesWithNodes.push_back(nodesInFirstBox);  
            
            //initialize lists
            
            std::list<std::vector<NodeInBox> >::iterator leftBox = boxesWithNodes.begin();
            for (std::list<NodeID>::const_iterator it = ++(chain.begin()); it != --(chain.end()); it++) {                
                
                
                //Print("boxesWithNodes" << boxesWithNodes.size());
                
                
                
                
                PrioNodeHeap::handle_type lastInserted
                        = prioNodeHeap.push(PrioNode2(*it, intervallIt, *leftBox, *leftBox);                

                //reference each other                
                intervallIt->prioNodeHandles.emplace_back(lastInserted);
                std::list<PrioNodeHandle>::iterator pos = --(intervallIt->prioNodeHandles.end());                
                PrioNode2 &prioNodeCopy = (*lastInserted);                
                prioNodeCopy.posInIntervallIt = pos;                                
                
                PrioNode2 prioNodeAssert(*lastInserted);
                assert(prioNodeAssert.posInIntervallIt == pos);                              
            }
            setIntervall(intervalls.back());             
    }
    
    
public:
    DiscreteCurveEvolution (GraphT &graph): graph(graph) {
        
    }
    
    ~DiscreteCurveEvolution() {        
    }
    
    std::list<ls::simplePrioNode> match(Chain chain1, Chain chain2) {
        assert(chain1.size() >= 3);
        prioNodeHeap.clear();
        
        //initialize
        std::list<>
        Intervall2 intervall
        
        
        //
        while (!prioNodeHeap.empty()) {
            
        }
        
        
    }
};

}