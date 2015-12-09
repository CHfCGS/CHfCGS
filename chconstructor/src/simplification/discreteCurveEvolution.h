#pragma once

#include "../nodes_and_edges.h"
#include "prio_nodes.h"
#include "cross_border.h"
#include "../grid.h"
//#include "lineSimplifierType.h"
#include "lineSimplifier.h"
#include "../geoFunctions.h"
#include "matchChainPairNodes.h"
#include "error_measure.h"
#include "../s_options.h"

#include <limits>
#include <math.h>

namespace ls {


//Bottom Up generalization

template<class GraphT>
class DiscreteCurveEvolution: public LineSimplifier {
private:
    
    PrioNodeHeap prioNodeHeap;
    std::list<Intervall2> intervalls;
    
    const GraphT &graph;
    const Grid<GraphT> &grid;
    
    
    std::list<NodeBox> node_boxes;
    
    std::unique_ptr<ErrorMeasure> errorMeasure;
    const SOptions s_options;
    
     
    
    std::list<ls::simplePrioNode> simplify () {        
        std::list<simplePrioNode> priolist;       
        
        while (!prioNodeHeap.empty()) {
            const PrioNode2 max =  prioNodeHeap.top();            
            //std::list<PPrioNode>::iterator maxIt = getMax();
            simplePrioNode spnMax(max);
            //PPrioNode max = getMax(nodePrios);   

            //update guidance of max, i.e. it can't follow any other node anymore
            for (PrioNodeHandle guide_h : max.guides) {
                PrioNode2 &guide = *guide_h;
                guide.followerValid = false;
            }
            
            //check if there is a follower
            if (max.followerValid) {
            //if(false) {
                //std::list<PrioNode2>::iterator followerIt = max.follower;
                const PrioNodeHandle follower_h = max.follower_h;
                const PrioNode2 follower = *follower_h;
                assert(follower.node_id != max.node_id);
                simplePrioNode spnFollower(follower);
                //nodePrios.erase(max.follower);
                                  
                //update guidance of follower
                for (PrioNodeHandle guide_h : follower.guides) {
                    PrioNode2 &guide = *guide_h;
                    guide.followerValid = false;
                }
                
                if (follower.followerValid) {
                    //the follower of follower has one guide less now
                    PrioNode2 &follower_follower = *follower.follower_h;
                    size_t nofFollowerFollowerguides = follower_follower.guides.size();
                    follower_follower.guides.remove(follower_h);
                    debug_assert(nofFollowerFollowerguides-1 == follower_follower.guides.size());
                }
                
                const PrioNode2 maxAssert =  prioNodeHeap.top();     
                assert(max.node_id == maxAssert.node_id);
                prioNodeHeap.pop();   //remove max, weights didnt change yet                                        
                prioNodeHeap.erase(follower_h);                
                
                assert(max.intervallIt != follower.intervallIt);
                
                smoothOut(max);
                smoothOut(follower);
                                
                //PrioNodes have to be pushed in that order such that most important is last
                priolist.push_back(simplePrioNode(spnMax));
                priolist.push_back(simplePrioNode(spnFollower));                                
            } else {
                prioNodeHeap.pop();  
                smoothOut(max);
                
                //r.push_back(simplePrioNode(max));
                //std::list<simplePrioNode> r = split();//(max.intervallIt); //, max.posInIntervall);                    
                priolist.push_back(spnMax);
            }         
        }
        
        return priolist;
    }       
    
    void smoothOut(PrioNode2 pn) {
        assert(!pn.intervallIt->prioNodeHandles.empty());
        if (pn.intervallIt->prioNodeHandles.size() == 1) {
            pn.intervallIt->prioNodeHandles.erase(pn.posInIntervallIt);
        }
        else if (pn.posInIntervallIt == pn.intervallIt->prioNodeHandles.begin()) {
            auto nextIt = pn.posInIntervallIt;
            nextIt++;
            PrioNodeHandle next_h = *nextIt;
            PrioNode2 &next = *next_h;
            next.left_node_boxes_its.splice(next.left_node_boxes_its.begin(),
                                            pn.left_node_boxes_its,
                                            pn.left_node_boxes_its.begin(),
                                            pn.left_node_boxes_its.end());
            pn.intervallIt->prioNodeHandles.erase(pn.posInIntervallIt);
            update(*pn.intervallIt, nextIt);
        } else if (pn.posInIntervallIt == --pn.intervallIt->prioNodeHandles.end()) {
            auto prevIt = pn.posInIntervallIt;
            prevIt--;
            PrioNodeHandle prev_h = *prevIt;
            PrioNode2 &prev = *prev_h;
            prev.right_node_boxes_its.splice(prev.right_node_boxes_its.end(),
                                            pn.right_node_boxes_its,
                                            pn.right_node_boxes_its.begin(),
                                            pn.right_node_boxes_its.end());
            pn.intervallIt->prioNodeHandles.erase(pn.posInIntervallIt);
            update(*pn.intervallIt, prevIt);            
        } else {
            auto prevIt = pn.posInIntervallIt;
            prevIt--;
            PrioNodeHandle prev_h = *prevIt;
            PrioNode2 &prev = *prev_h;
            prev.right_node_boxes_its.splice(prev.right_node_boxes_its.end(),
                                            pn.right_node_boxes_its,
                                            pn.right_node_boxes_its.begin(),
                                            pn.right_node_boxes_its.end());     
            auto nextIt = pn.posInIntervallIt;
            nextIt++;
            PrioNodeHandle next_h = *nextIt;
            PrioNode2 &next = *next_h;
            next.left_node_boxes_its.splice(next.left_node_boxes_its.begin(),
                                            pn.left_node_boxes_its,
                                            pn.left_node_boxes_its.begin(),
                                            pn.left_node_boxes_its.end());
            pn.intervallIt->prioNodeHandles.erase(pn.posInIntervallIt);
            update(*pn.intervallIt, prevIt);
            update(*pn.intervallIt, nextIt);
        } 
    }
    
    void update(Intervall2 &intervall, const std::list<PrioNodeHandle>::iterator posInIntervallIt) {
        PrioNodeHandle prev_h;
        PrioNodeHandle cur_h;
        PrioNodeHandle next_h;
        assert(!intervall.prioNodeHandles.empty());
        if (intervall.prioNodeHandles.size() == 1) {
            cur_h = *intervall.prioNodeHandles.begin();
            PrioNode2 &cur = *cur_h;            
            cur.perpendicularLength = errorMeasure->calcError(graph.getNode(intervall.start),
                                                              graph.getNode(intervall.finish),
                                                              graph.getNode(cur.node_id));            
        } else if (posInIntervallIt == intervall.prioNodeHandles.begin()) {
            //front
            cur_h = *intervall.prioNodeHandles.begin();
            PrioNode2 &cur = *cur_h;
            next_h = *++(intervall.prioNodeHandles.begin());
            PrioNode2 &next = *next_h;            
            cur.perpendicularLength = errorMeasure->calcError(graph.getNode(intervall.start),
                                                              graph.getNode(next.node_id),
                                                              graph.getNode(cur.node_id));
            
        } else if (posInIntervallIt == --intervall.prioNodeHandles.end()) {
            //back            
            prev_h = *--(--intervall.prioNodeHandles.end());
            PrioNode2 &prev = *prev_h;
            cur_h = *--intervall.prioNodeHandles.end();
            PrioNode2 &cur = *cur_h;            
            cur.perpendicularLength = errorMeasure->calcError(graph.getNode(prev.node_id),
                    graph.getNode(intervall.finish),
                    graph.getNode(cur.node_id));            
        } else {
            //middle
            auto lastIt = posInIntervallIt;
            lastIt--;
            auto nextIt = posInIntervallIt;
            nextIt++;

            prev_h = *lastIt;
            PrioNode2 &prev = *prev_h;
            cur_h = *posInIntervallIt;
            PrioNode2 &cur = *cur_h;
            next_h = *nextIt;
            PrioNode2 &next = *next_h;
            
            cur.perpendicularLength = errorMeasure->calcError(graph.getNode(prev.node_id),
                        graph.getNode(next.node_id),
                        graph.getNode(cur.node_id)); 
            
        }        
        prioNodeHeap.update(cur_h);  
    }
    
    void setIntervall(Intervall2 &intervall){
        for(std::list<PrioNodeHandle>::iterator it = intervall.prioNodeHandles.begin();
                it != intervall.prioNodeHandles.end();
                it++) {
            update(intervall, it);
        }
    }
    
    /*
    void setIntervall(Intervall2 &intervall){
        if(intervall.prioNodeHandles.empty()) {
            return;
        } else if(intervall.prioNodeHandles.size() == 1) {
            PrioNodeHandle prioNodeH = *intervall.prioNodeHandles.begin();
            PrioNode2 &prioNode = *prioNodeH;
            //Kink kink(intervall.start, prioNode.node_id, intervall.finish);
            prioNode.perpendicularLength = errorMeasure->calcError(intervall.start, intervall.finish, prioNode.node_id)
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
            //Kink kink(intervall.start, cur.node_id, next.node_id);
            cur.perpendicularLength = errorMeasure->calcError(intervall.start, next.node_id, cur.node_id);
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
                
                //Kink kink(prev.node_id, cur.node_id, next.node_id);
                cur.perpendicularLength = errorMeasure->calcError(prev.node_id, next.node_id, cur.node_id);
                prioNodeHeap.update(cur);                
            }
            //back            
            prev_h = *--(intervall.prioNodeHandles.end());
            PrioNode2 &prev = *prev_h;
            cur_h = *intervall.prioNodeHandles.end();
            PrioNode2 &cur = *cur_h;
            //Kink kink(prev.node_id, cur.node_id, intervall.finish);
            cur.perpendicularLength = errorMeasure->calcError(prev.node_id, intervall.finish, cur.node_id);
            prioNodeHeap.update(cur);
        }                        
    }*/
    
    void initializeChain(const Chain &chain) {                   
        
        debug_assert(chain.size() >= 3);            

        Intervall2 initialIntervall(chain.front(), chain.back());
        intervalls.emplace_back(initialIntervall);
        const auto intervallIt = --(intervalls.end());

        //first box of a chain
        std::vector<NodeInBox> nodesInFirstBox(cb::calcNodesInBoundingBox(*(chain.begin()), *(++(chain.begin())), graph, grid));
        node_boxes.push_back(nodesInFirstBox);  

        //initialize lists                        
        std::list<std::vector<NodeInBox> >::iterator leftBox = --node_boxes.end(); //Iterator which runs in tandem with chainNodeIterator
        for (std::list<NodeID>::const_iterator it = ++(chain.begin()); it != --(chain.end()); it++) {                
            auto nextIt = it;
            nextIt++;
            
            //getting iterators on the left and right boxes of a node
            auto rightBox = leftBox;
            if (s_options.checkBorderCrossing) {
                
                std::vector<NodeInBox> nodesInBox(cb::calcNodesInBoundingBox(*it, *nextIt, graph, grid));                
                //uint sizeBefore = boxesWithNodes.size();
                node_boxes.push_back(nodesInBox);
                //Print("boxesWithNodes" << boxesWithNodes.size());
                rightBox++;
                assert(rightBox == --(node_boxes.end()));  
            }

            PrioNode2 prio_node(*it, intervallIt, *leftBox, *rightBox);
    
            PrioNodeHeap::handle_type lastInserted
                    = prioNodeHeap.push(prio_node);
            
            if (s_options.checkBorderCrossing) {
                PrioNode2 &prioNodeRef = (*lastInserted);                                     
                prioNodeRef.left_node_boxes_its.push_back(leftBox);
                prioNodeRef.right_node_boxes_its.push_back(rightBox);
                leftBox++;
            }
            
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
    DiscreteCurveEvolution (SOptions s_options, const GraphT &graph, const Grid<GraphT> &grid):
                prioNodeHeap(), intervalls(), graph(graph), grid(grid), node_boxes(), s_options(s_options) {
        errorMeasure = createErrorMeasure(s_options.errorMeasure_type);
        assert(errorMeasure != nullptr);
    }    
    
    ~DiscreteCurveEvolution() {        
    }
    
    std::list<ls::simplePrioNode> process(const Chain &chain1, const Chain &chain2) {
        assert(chain1.size() >= 3);
        prioNodeHeap.clear();
        intervalls.clear();
        node_boxes.clear(); 
                
        
        if (chain2.empty()) {
            initializeChain(chain1);
            intervalls.back().prioNodeHandles.reverse();
            return simplify();
        } else {
            debug_assert(chain1.size() + chain2.size() > 6);
            debug_assert(chain1.size() >= 3 && chain2.size() >= 3);

            initializeChain(chain1);
            initializeChain(chain2);

            intervalls.back().prioNodeHandles.reverse();

            mc::match<GraphT> (graph, intervalls.front().prioNodeHandles, intervalls.back().prioNodeHandles, s_options.pairMatch_type);
            //matchChainPairNodes2<GraphT> matcher(base_graph, intervalls.front(), intervalls.back());
            //matcher.match();                
            return simplify();
        }
    }
};

}
/*
std::list<simplePrioNode> process(const Chain &chain1, const Chain &chain2) {                
            debug_assert(chain1.size() >= 3);
            
            prioNodeHeap.clear();
            intervalls.clear();
            node_boxes.clear();  
            
            if(chain2.empty()) {
                initializeChain(chain1);
                intervalls.back().prioNodeHandles.reverse();
                return simplify();
            } else {
                debug_assert(chain1.size() + chain2.size() > 6);
                debug_assert(chain1.size() >=3 && chain2.size() >= 3);              
                                                
                initializeChain(chain1);
                initializeChain(chain2);

                intervalls.back().prioNodeHandles.reverse();
                
                mc::match<GraphT> (graph, intervalls.front().prioNodeHandles, intervalls.back().prioNodeHandles, s_options.pairMatch_type);
                //matchChainPairNodes2<GraphT> matcher(base_graph, intervalls.front(), intervalls.back());
                //matcher.match();                
                return simplify();
            }
            
        }    */