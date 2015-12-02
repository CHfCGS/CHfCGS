/* 
 * File:   dpPrioritizer.h
 * Author: tobias
 *
 * Created on 6. Oktober 2015, 13:10
 */

#pragma once

#include "prioritizers.h"
#include "simplification/lineSimplifierType.h"


#include "simplification/lineSimplifier.h"
#include "simplification/dp_simplifier.h"

#include "chains.h"
#include "grid.h"
#include "chainDetector.h"
#include "DouglasPeucker.h"
#include "indexed_container.h"
#include "dpPrioritizer.h"
#include "dead_end_detector.h"

#include "4Dgrid.h"
#include "chainPairDP.h"

#include "simplification/prio_nodes.h"
#include "nodes_and_edges.h"

/*
 * Prioritizer that prioritizes chain nodes by DouglasPeucker.
 */
namespace chc {
    /*
    template <class GraphT, class CHConstructorT>
    class EdgeDiffPrioritizer;
    */
    template <class GraphT, class CHConstructorT>
    class DPPrioritizer : public Prioritizer  {// : public EdgeDiffPrioritizer<GraphT, CHConstructorT> {

        enum class State {
            start, removingDeadEnds, removingChains//, removingRemaining, done
        };

    private:
        
        struct CompInOutProduct {
            GraphT const& g;

            CompInOutProduct(GraphT const& g)
                    : g(g) {}

            bool operator()(NodeID node1, NodeID node2) const
            {
                    uint edge_product1(g.getNrOfEdges(node1, EdgeType::IN)
                                    * g.getNrOfEdges(node1, EdgeType::OUT));
                    uint edge_product2(g.getNrOfEdges(node2, EdgeType::IN)
                                    * g.getNrOfEdges(node2, EdgeType::OUT));

                    return edge_product1 < edge_product2;
            }
        };

	GraphT const& _base_graph;
        CHConstructorT const& _chc;
        std::vector<NodeID> _prio_vec;
                        
        State state;
        Grid<GraphT> grid;
        FourDGrid<GraphT> fourDGrid;
        
        DeadEndDetector<GraphT> deadEndDetector;
        ChainDetector<GraphT> chaindetector;        
        ls::DouglasPeucker<GraphT> dp;
        ls::chainPairDP<GraphT> cpdp;
                
        std::unique_ptr<ls::LineSimplifier> lineSimplifier; 
        //node prio list
        //std::vector<NodeID> _chooseIndependentSet();		
        //std::vector<chc::Chain> chains;

        std::list<Chain> deadEnds;
        //_prio_vec = chains + remainder        
        Chains_and_Remainder CaR;
        std::vector<std::list<ls::simplePrioNode>> priolists;

        double epsilon;
        int roundcounter;
        
        void init(std::vector<NodeID>& node_ids) {
            _prio_vec = std::move(node_ids);
        }
        
        void _remove(std::vector<NodeID> const& nodes) {
            std::vector<bool> to_remove(_base_graph.getNrOfNodes(), false);
            for (auto node : nodes) {
                debug_assert(0 <= node && node < to_remove.size());
                to_remove[node] = true;
            }

            size_t remaining_nodes(_prio_vec.size());
            size_t i(0);
            while (i < remaining_nodes) {
                debug_assert(0 <= i && i < _prio_vec.size());
                NodeID node(_prio_vec[i]);
                if (to_remove[node]) {
                    remaining_nodes--;
                    debug_assert(0 <= remaining_nodes && remaining_nodes < _prio_vec.size());
                    _prio_vec[i] = _prio_vec[remaining_nodes];
                    _prio_vec[remaining_nodes] = node;
                } else {
                    i++;
                }
            }

            _prio_vec.resize(remaining_nodes);
        }
        
        void _removeFromRemainder(std::vector<NodeID> const& nodes) {
            std::vector<bool> to_remove(this->_base_graph.getNrOfNodes(), false);
            for (auto node : nodes) {     
                debug_assert(0 <= node && node < to_remove.size());
                to_remove[node] = true;                
            }

            size_t remaining_nodes(CaR.remainder.size());
            size_t i(0);
            while (i < remaining_nodes) {                
                debug_assert(0 <= i && i < CaR.remainder.size());
                NodeID node(CaR.remainder[i]);                
                if (to_remove[node]) {
                    remaining_nodes--;
                    debug_assert(0 <= remaining_nodes && remaining_nodes < CaR.remainder.size());
                    CaR.remainder[i] = CaR.remainder[remaining_nodes];
                    CaR.remainder[remaining_nodes] = node;
                } else {
                    i++;
                }
            }

            CaR.remainder.resize(remaining_nodes);
        }

        //TODO slowdown
        std::vector<NodeID> _chooseIndependentSetFromRemainder() {
            //CompInOutProduct ciop(this->_base_graph);            
            if (CaR.remainder.empty()) {
                std::vector<NodeID> empytList;
                return empytList;
            } else {
                std::sort(CaR.remainder.begin(), CaR.remainder.end(), CompInOutProduct(_base_graph));
                //auto independent_set(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc.calcIndependentSet(CaR.remainder));
                //auto edge_diffs(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc.calcEdgeDiffs(independent_set));

                auto independent_set(_chc.calcIndependentSet(CaR.remainder));
                Print("Calc EdgeDiffs");
                auto edge_diffs(_chc.calcEdgeDiffs(independent_set));
                double edge_diff_mean(0);
                for (size_t i(0); i < edge_diffs.size(); i++) {                    
                    edge_diff_mean += edge_diffs[i];
                }
                debug_assert(independent_set.size() != 0);
                edge_diff_mean /= independent_set.size();

                std::vector<NodeID> low_edge_diff_nodes;
                for (size_t i(0); i < independent_set.size(); i++) {
                    debug_assert(0 <= i && i < edge_diffs.size());
                    if (edge_diffs[i] <= edge_diff_mean) {
                        NodeID node(independent_set[i]);
                        low_edge_diff_nodes.push_back(node);
                    }
                }

                return low_edge_diff_nodes;
            }
        }        
        
        void FillChainsInPriolists(std::vector<ChainsOfType> &chainsaccordingToType) {
            for (ChainsOfType &chainsOfType: chainsaccordingToType) {                                                                
                
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                for (Chain &chain: chainsOfType) {                                                                
                    //big chains are generalized
                    if (chain.size() >= 3) {                                      
                        
                        priolists.push_back(lineSimplifier->process(chain, Chain()));
                        //priolists.push_back(dp.process(chain));
                        
                    //small chains are assigned to the remainder
                    } else {
                        for (NodeID node_id: chain) {
                            CaR.remainder.push_back(node_id);
                        }                    
                    }   
                }
            }
        }
        
        void FillPriolists() {
            priolists.clear();        
            FillChainsInPriolists(CaR.oneWayChainsAccordingToType);
            FillChainsInPriolists(CaR.twoWayChainsAccordingToType);                        
                        
            for (ChainPair &chainPair: CaR.chainPairs) {    
                    if (chainPair.chainTo.size() + chainPair.chainFrom.size() >= 7
                            && chainPair.chainTo.size() >=3 && chainPair.chainFrom.size() >= 3) { 
                        priolists.push_back(lineSimplifier->process(chainPair.chainTo, chainPair.chainFrom)); 
                        //priolists.push_back(cpdp.process(chainPair));                        
                        //Print("Length of Priolist: " << pl.size());                                        
                    //small chains are assigned to the remainder
                    } else {
                        for (NodeID node_id: chainPair.chainTo) {
                            CaR.remainder.push_back(node_id);
                        }                    
                        for (NodeID node_id: chainPair.chainFrom) {
                            CaR.remainder.push_back(node_id);
                        }                    
                    }    
            }
        }
        
        std::vector<NodeID> removeDeadEnds(std::list<Chain> &deadEnds) {
            std::vector<NodeID> next_nodes;
            for (auto it = deadEnds.begin(); it != deadEnds.end();) {
                if (it->empty()) {
                    it = deadEnds.erase(it);
                }else {                    
                    next_nodes.push_back(it->front());
                    it->pop_front();
                    it++;
                }
            }
            return next_nodes;
        }
        
        /*
        std::vector<NodeID> removeRemainder() {
            std::vector<NodeID> next_nodes;
            Print("Getting Independent set from Remainder");
            //independent set from remainder
            
           
            next_nodes = _chooseIndependentSetFromRemainder();            
            //remove from remainder
            _removeFromRemainder(next_nodes);
            return next_nodes;
        }
         * */
        
        
    public:
        /*
        DPPrioritizer(GraphT const& base_graph, CHConstructorT const& chc)
        : EdgeDiffPrioritizer<GraphT, CHConstructorT>(base_graph, chc), state(State::start),
        grid(1000, base_graph), fourDGrid(1, base_graph), chaindetector(base_graph), dp(this->_base_graph, grid),
        cpdp(this->_base_graph, grid), CaR(), priolists(), epsilon(0.0001), roundcounter(1) {
        }*/
        
        DPPrioritizer(GraphT const& base_graph, CHConstructorT const& chc)
        : _base_graph(base_graph), _chc(chc), state(State::removingDeadEnds),
        grid(1000, base_graph), fourDGrid(1, base_graph), deadEndDetector(base_graph), chaindetector(base_graph), dp(this->_base_graph, grid),
        cpdp(this->_base_graph, grid), CaR(), priolists(), epsilon(10000), roundcounter(1) {            
            lineSimplifier = createLineSimplifer(ls::LineSimplifierType::DP, base_graph, grid);
        }        
        //epsilon(0.0001)
        ~DPPrioritizer() {
            //Print("DPPrioritizer is destructed");
        }

        std::vector<NodeID> extractNextNodes() {
            bool empty = _prio_vec.empty();
            debug_assert(!empty);
            std::vector<NodeID> next_nodes;
            
                        
            switch (state) {
                case State::start: {
                    deadEnds = deadEndDetector.detectDeadEnds(_prio_vec);
                    state = State::removingDeadEnds;
                    break;
                }
                case State::removingDeadEnds: {
                    next_nodes = removeDeadEnds(deadEnds);
                    if (!next_nodes.empty()) {
                        break;
                    } else {
                        state = State::removingChains;
                        break;
                    }
                }
                case State::removingChains: {
                    if ((roundcounter - 1) % 5 == 0) {
                        if (epsilon < 100000) { //prevents epsilon overflow
                            epsilon *= 2.0;
                        }

                        Print("Detecting chains");
                        
                        CaR = chaindetector.detectChains(_prio_vec);
                        Print("Number of chains: " << CaR.getNrOfChains());
                        debug_assert(CaR.getNrOfNodesInChains() + CaR.remainder.size() == this->_prio_vec.size());

                        Print("IdentifyingChainPairs ");
                        fourDGrid.identifyPairs(CaR);
                        Print("Number of chain pairs: " << CaR.chainPairs.size());

                        FillPriolists();
                    }

                    Print("Getting Independent set from Remainder");
                    //independent set from remainder

                    if (roundcounter > 20) {
                        next_nodes = _chooseIndependentSetFromRemainder();
                    }
                    //remove from remainder
                    _removeFromRemainder(next_nodes);
                    //std::vector<NodeID> next_nodes;

                    Print("Getting Independent set from Chains");
                    //calc independent set of nodes from priolists
                    std::vector<bool> marked(this->_base_graph.getNrOfNodes(), false);

                    //if (marked.at(marked.size()+1)) {Print("Getting Independent set from Chains");}

                    //extraction from chains
                    for (std::list<ls::simplePrioNode> &priolist : priolists) {
                        for (auto it = priolist.begin(); it != priolist.end();) {
                            if (epsilon > it->perpendicularLength) {
                                NodeID node_id = it->node_id;
                                if (!marked[node_id]) {

                                    //EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc._markNeighbours(node_id, marked); 
                                    _chc._markNeighbours(node_id, marked);
                                    next_nodes.push_back(node_id);
                                    it = priolist.erase(it);
                                } else {
                                    break;
                                }
                            } else {
                                break;
                            }


                        }
                    }
                    roundcounter++;
                    break;
                }
                default: {
                    Print("DefaultCase");
                    break;
                }
            }
            
        _remove(next_nodes); //remove from priovector               

        /*
        for (NodeID node_id: next_nodes) {
            Print(node_id);
        }
         * */
        return next_nodes;
    }
        
        bool hasNodesLeft() {
            return !_prio_vec.empty();
        }

    };
}