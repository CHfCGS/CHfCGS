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

template <class GraphT>
class DPSimplifier : public LineSimplifier{
    private:   
        PrioNodeHeap prioNodeHeap;
        //std::list<PrioNode2> PrioNodes;
        std::list<Intervall2> intervalls;
        //vector<Node> &nodes;
        const GraphT &graph;
        const Grid<GraphT> &grid;
        std::list<NodeBox> node_boxes;
        std::unique_ptr<ErrorMeasure> errorMeasure;
        const SOptions s_options;        
     
    std::list<simplePrioNode> simplify() {
        if (prioNodeHeap.empty()) {
            return std::list<simplePrioNode>(); //empty list
        }
        else {            
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
                split(max);
                split(follower);
                
                std::list<simplePrioNode> r = simplify();
                //PrioNodes have to be pushed in that order such that most important is last
                r.push_back(simplePrioNode(spnFollower));                
                r.push_back(simplePrioNode(spnMax));
                
                return r;
            } else {
                prioNodeHeap.pop();  
                split(max);
                std::list<simplePrioNode> r = simplify();
                //r.push_back(simplePrioNode(max));
                //std::list<simplePrioNode> r = split();//(max.intervallIt); //, max.posInIntervall);                    
                r.push_back(spnMax);
                return r;
        
            }
            
            //nodePrios.erase(maxIt);
        }     
    }
        
    void split(PrioNode2 max) {   
        /*
        if (prioNodeHeap.empty()) {
            return std::list<simplePrioNode>(); //empty list
        }
        else {
                        
            const PrioNode2 max =  prioNodeHeap.top();
            */
            //
            
            //std::list<PrioNode>::iterator maxIt = PrioNodes.top();// getMax();                         
            
            std::list<Intervall2>::iterator splitIt = max.intervallIt;            
            std::list<PrioNodeHandle>::iterator splitPos = max.posInIntervallIt;
            
            /*
            const PrioNodeHandle splitPh = *splitPos;
            const PrioNode2 &splitP = *splitPh; //equal to max
            const NodeID split_node_id = splitP.node_id; 
            */
            
            //Intervall2 leftIntervall(splitIt->start, split_node_id);
            Intervall2 leftIntervall(splitIt->start, max.node_id);
            //Intervall2 rightIntervall(split_node_id, splitIt->end);                                    
            Intervall2 rightIntervall(max.node_id, splitIt->finish);                                    
            std::list<Intervall2>::iterator leftIt = intervalls.insert(splitIt, leftIntervall);           
            std::list<Intervall2>::iterator rightIt = intervalls.insert(splitIt, rightIntervall);                                                                                                
            
            //left Intervall
            leftIt->prioNodeHandles.splice(leftIt->prioNodeHandles.begin(), splitIt->prioNodeHandles, splitIt->prioNodeHandles.begin(), splitPos);            
            for (std::list<PrioNodeHandle>::iterator it = leftIt->prioNodeHandles.begin(); it!= leftIt->prioNodeHandles.end(); it++) {                
                PrioNode2 &p = (*(*it));
                p.intervallIt = leftIt;
                p.posInIntervallIt = it;  
                //prioNodes.update(*it, p);
            }

            //right Intervall  
            rightIt->prioNodeHandles.splice(rightIt->prioNodeHandles.begin(), splitIt->prioNodeHandles, ++splitPos, splitIt->prioNodeHandles.end());            
            for (std::list<PrioNodeHandle>::iterator it = rightIt->prioNodeHandles.begin(); it!= rightIt->prioNodeHandles.end(); it++) {                                
                PrioNode2 &p = (*(*it));
                p.intervallIt = rightIt;
                p.posInIntervallIt = it;                 
                //prioNodes.update(*it, p);
            }                                    
            
            intervalls.erase(splitIt);
              
            
            
            update (*leftIt, max.nof_crossings.left);
            update (*rightIt, max.nof_crossings.right);
            
            /*
            simplePrioNode spn(max);
            
            
            
            std::list<simplePrioNode> r = split();
            r.push_back(spn);            
            return r;
             * */
    }
    
    /*
    //calculates for a line all nodes in the lines' bounding box and also save their orientation with respect to the line
    std::vector<NodeInBox> calcNodesInBoundingBox(NodeID node_id_src, NodeID node_id_tgt) {
        BoundingBox bb = BoundingBox(
                graph.getLat(node_id_src),
                graph.getLon(node_id_src),
                graph.getLat(node_id_tgt),
                graph.getLon(node_id_tgt));

        std::vector<NodeID> nodesInGrid = grid.nodesInNeigbourhood(bb.midLat, bb.midLon);
        std::vector<NodeInBox> nodesInBox;
        for (NodeID node_id: nodesInGrid) {                    
            //consider only the nodes in the bounding box
            if (bb.contains(graph.getLat(node_id), graph.getLon(node_id))) {
                //int i1   =       3 > 4 ? 0 : 1;
                //save on which side of the line the node lies
                bool sign = (geo::calcArea(node_id_src, node_id_tgt, node_id, graph) >= 0) ? true : false;                
                nodesInBox.push_back(NodeInBox(node_id, sign));
            }            
        }
        return nodesInBox;
    }      */     
           
    
    NofCrossings calcOrientationMisses(const PrioNode2 &split, const Intervall2 &intervall) {
        NofCrossings nof_crossings;        
        uint assertcounter = 0;
        auto splitPosIntervallincr = split.posInIntervallIt;
        splitPosIntervallincr++;
        for (auto it = intervall.prioNodeHandles.begin(); it != splitPosIntervallincr; it++) {
            const PrioNode2 &p = (*(*it));
            for (NodeInBox nodeInBox : p.leftBox) {
                bool sign = (geo::calcArea(intervall.start, split.node_id, nodeInBox.nodeID, graph) >= 0) ? true : false;
                if (sign != nodeInBox.sign) {
                    nof_crossings.left++;
                }
            }
            assertcounter++;
        }
        for (auto it = split.posInIntervallIt; it != intervall.prioNodeHandles.end(); it++) {
            const PrioNode2 &p = (*(*it));
            for (NodeInBox nodeInBox : p.rightBox) {
                bool sign = (geo::calcArea(split.node_id, intervall.finish, nodeInBox.nodeID, graph) >= 0) ? true : false;
                if (sign != nodeInBox.sign) {
                    nof_crossings.right++;
                }
            }
            assertcounter++;
        }
        debug_assert(assertcounter == intervall.prioNodeHandles.size() + 1);
        return nof_crossings;
    }
    
    /*
     * updates the priority data for all prionodes in an Intervall
     */    
    void update(Intervall2 &intervall, uint pre_nof_crossings){
        
        for (PrioNodeHandle prioNodeH: intervall.prioNodeHandles) {
            PrioNode2 &prioNode = *prioNodeH;            
            prioNode.error = errorMeasure->calcError(graph.getNode(intervall.start),
                                                                        graph.getNode(intervall.finish),
                                                                        graph.getNode(prioNode.node_id));
            
            prioNode.nof_crossings = calcOrientationMisses(prioNode, intervall);
                                    
            prioNode.cross_diff = (int) prioNode.nof_crossings.getSum() - (int) pre_nof_crossings;
            prioNodeHeap.update(prioNodeH);
        }
        return;
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

                PrioNodeHeap::handle_type lastInserted
                        = prioNodeHeap.push(PrioNode2(*it, intervallIt, *leftBox, *rightBox));
                if (s_options.checkBorderCrossing) {
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
            update(intervalls.back(), 0);             
    }
    
        
    public:        
        DPSimplifier(const SOptions s_options, const GraphT &graph, const Grid<GraphT> &grid):
                prioNodeHeap(), intervalls(), graph(graph), grid(grid), node_boxes(), s_options(s_options)  {
            errorMeasure = createErrorMeasure(s_options.errorMeasure_type);            
            assert(errorMeasure != nullptr);
        }
        
        ~DPSimplifier() {            
        }
        
        std::list<simplePrioNode> process(const Chain &chain1, const Chain &chain2) {                
            debug_assert(chain1.size() >= 3);
            
            prioNodeHeap.clear();
            intervalls.clear();
            node_boxes.clear();  
            
            if(chain2.empty()) {
                initializeChain(chain1);                
                return simplify();
            } else {
                debug_assert(chain1.size() + chain2.size() > 6);
                debug_assert(chain1.size() >=3 && chain2.size() >= 3);              
                                                
                initializeChain(chain1);
                initializeChain(chain2);

                //one chain needs to reversed for matching
                intervalls.back().prioNodeHandles.reverse();                
                mc::match<GraphT> (graph, intervalls.front().prioNodeHandles, intervalls.back().prioNodeHandles, s_options.pairMatch_type);
                //matchChainPairNodes2<GraphT> matcher(base_graph, intervalls.front(), intervalls.back());
                //matcher.match();                
                return simplify();
            }
            
        }                            
};

}


