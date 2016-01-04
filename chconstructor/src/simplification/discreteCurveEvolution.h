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
#include "../bounding_box.h"

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
            update(next);
            //update(*pn.intervallIt, nextIt);
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
            update(prev);
            //update(*pn.intervallIt, prevIt);            
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
            update(prev);
            update(next);
            //update(*pn.intervallIt, prevIt);
            //update(*pn.intervallIt, nextIt);            
        } 
    }
    
    
    uint _calcOrientationMisses(const std::list<std::list<NodeBox>::iterator> &node_boxes_its, NodeID src, NodeID tgt) {
        uint crossings_counter = 0;
        for (std::list<NodeBox>::iterator it: node_boxes_its) {
            NodeBox &node_box = *it;                   
            for (NodeInBox nodeInBox : node_box) {
                bool sign = (geo::calcArea(src, tgt, nodeInBox.nodeID, graph) >= 0) ? true : false;
                if (sign != nodeInBox.sign) {
                    crossings_counter++;
                }
            }
        }
        return crossings_counter;
    }
    
    NofCrossings calcOrientationMisses(const PrioNode2 &pn, NodeID src, NodeID tgt) {        
        NofCrossings nof_crossings;                        
        /*
        for (std::list<NodeBox>::iterator it: pn.left_node_boxes_its) {
            NodeBox &node_box = *it;                   
            for (NodeInBox nodeInBox : node_box) {
                bool sign = (geo::calcArea(src, tgt, nodeInBox.nodeID, graph) >= 0) ? true : false;
                if (sign != nodeInBox.sign) {
                    nof_crossings.left++;
                }
            }
        }
                    
        for (std::list<NodeBox>::iterator it: pn.right_node_boxes_its) {
            NodeBox &node_box = *it;
            for (NodeInBox nodeInBox : node_box) {
                bool sign = (geo::calcArea(src, tgt, nodeInBox.nodeID, graph) >= 0) ? true : false;
                if (sign != nodeInBox.sign) {
                    nof_crossings.right++;
                }
            }            
        }*/
        nof_crossings.left = _calcOrientationMisses(pn.left_node_boxes_its, src, tgt);
        nof_crossings.right = _calcOrientationMisses(pn.right_node_boxes_its, src, tgt);
        return nof_crossings;
    }
    
    void _updateErrors(PrioNode2 &cur, NodeID prev, NodeID next) {
        //negative weight since heap is a max heap it should always eliminate the least critical value
        cur.error = -errorMeasure->calcError(graph.getNode(prev),
                        graph.getNode(next),
                        graph.getNode(cur.node_id));            
        cur.cross_diff = calcOrientationMisses(cur, prev, next).getSum();
    }
    
    //void update(Intervall2 &intervall, const std::list<PrioNodeHandle>::iterator posInIntervallIt) {
    void update(PrioNode2 &pn) {
        Intervall2 &intervall = *pn.intervallIt;
        const std::list<PrioNodeHandle>::iterator &posInIntervallIt = pn.posInIntervallIt;
        PrioNodeHandle prev_h;
        PrioNodeHandle cur_h;
        PrioNodeHandle next_h;
        assert(!intervall.prioNodeHandles.empty());
        if (intervall.prioNodeHandles.size() == 1) {
            cur_h = *intervall.prioNodeHandles.begin();
            PrioNode2 &cur = *cur_h;
            _updateErrors(cur, intervall.start, intervall.finish);   
            /*
            cur.error = -errorMeasure->calcError(graph.getNode(intervall.start),
                        graph.getNode(intervall.finish),
                        graph.getNode(cur.node_id));            
            cur.cross_diff = calcOrientationMisses(cur, intervall.start, intervall.finish).getSum();*/
        } else if (posInIntervallIt == intervall.prioNodeHandles.begin()) {
            //front
            cur_h = *intervall.prioNodeHandles.begin();
            PrioNode2 &cur = *cur_h;
            next_h = *++(intervall.prioNodeHandles.begin());
            PrioNode2 &next = *next_h;   
            _updateErrors(cur, intervall.start, next.node_id);
            /*
            cur.error = -errorMeasure->calcError(graph.getNode(intervall.start),
                        graph.getNode(next.node_id),
                        graph.getNode(cur.node_id));
            cur.cross_diff = calcOrientationMisses(cur, intervall.start, next.node_id).getSum();*/
            
        } else if (posInIntervallIt == --intervall.prioNodeHandles.end()) {
            //back            
            prev_h = *--(--intervall.prioNodeHandles.end());
            PrioNode2 &prev = *prev_h;
            cur_h = *--intervall.prioNodeHandles.end();
            PrioNode2 &cur = *cur_h; 
            _updateErrors(cur, prev.node_id, intervall.finish);
            /*
            cur.error = -errorMeasure->calcError(graph.getNode(prev.node_id),
                        graph.getNode(intervall.finish),
                        graph.getNode(cur.node_id));
            cur.cross_diff = calcOrientationMisses(cur, prev.node_id, intervall.finish).getSum();  */          
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
        
            _updateErrors(cur, prev.node_id, next.node_id);
            /*
            cur.error = -errorMeasure->calcError(graph.getNode(prev.node_id),
                        graph.getNode(next.node_id),
                        graph.getNode(cur.node_id)); 
            cur.cross_diff = calcOrientationMisses(cur, prev.node_id, next.node_id).getSum(); */
        }                
        
        prioNodeHeap.update(cur_h);  
    }
    
    void setIntervall(Intervall2 &intervall){
        for(std::list<PrioNodeHandle>::iterator it = intervall.prioNodeHandles.begin();
                it != intervall.prioNodeHandles.end();
                it++) {
            PrioNodeHandle pn_h = *it;
            PrioNode2 &pn = *pn_h;
            update(pn);
            //update(intervall, it);
        }
    }        
    
    void initializeChain(const Chain &chain) {                   
        
        debug_assert(chain.size() >= 3);            

        Intervall2 initialIntervall(chain.front(), chain.back());
        intervalls.emplace_back(initialIntervall);
        const auto intervallIt = --(intervalls.end());

        //first box of a chain
        std::vector<NodeInBox> nodesInFirstBox;
        if (s_options.checkBorderCrossing) {
            nodesInFirstBox = cb::calcNodesInBoundingBox(*(chain.begin()), *(++(chain.begin())), graph, grid);
        }            
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
            return simplify();
        } else {
            debug_assert(chain1.size() + chain2.size() > 6);
            debug_assert(chain1.size() >= 3 && chain2.size() >= 3);

            initializeChain(chain1);
            initializeChain(chain2);

            //one chain need to reversed for matching
            intervalls.back().prioNodeHandles.reverse();
            mc::match<GraphT> (graph, intervalls.front().prioNodeHandles, intervalls.back().prioNodeHandles, s_options.pairMatch_type);
            //matchChainPairNodes2<GraphT> matcher(base_graph, intervalls.front(), intervalls.back());
            //matcher.match();                
            return simplify();
        }
    }
};

}