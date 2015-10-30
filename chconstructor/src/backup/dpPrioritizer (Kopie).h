/* 
 * File:   dpPrioritizer.h
 * Author: tobias
 *
 * Created on 6. Oktober 2015, 13:10
 */

#pragma once

#include "prioritizers.h"
/*
#include "nodes_and_edges.h"
#include "grid.h"
#include "chainDetector.h"
#include "DouglasPeucker.h"
 * */
/*
 * Prioritizer that prioritizes chain nodes by DouglasPeucker.
 */
namespace chc {
    template <class GraphT, class CHConstructorT>
    class EdgeDiffPrioritizer;

    template <class GraphT, class CHConstructorT>
    class DPPrioritizer : public EdgeDiffPrioritizer<GraphT, CHConstructorT> {

        enum class State {
            start, removingDeadEnds, removingChains, removingRemaining, done
        };

    private:
        State state;
        Grid<GraphT> grid;
        ChainDetector<GraphT> chaindetector;
        DP::DouglasPeucker<GraphT> dp;
        //node prio list
        //std::vector<NodeID> _chooseIndependentSet();		
        //std::vector<chc::Chain> chains;

        //_prio_vec = chains + remainder
        chc::Chains_and_Remainder CaR;
        std::vector<DP::PrioList> priolists;

        double epsilon;
        uint roundcounter;

        void _removeFromRemainder(std::vector<NodeID> const& nodes) {
            std::vector<bool> to_remove(this->_base_graph.getNrOfNodes(), false);
            for (auto node : nodes) {
                to_remove[node] = true;
            }

            size_t remaining_nodes(CaR.remainder.size());
            size_t i(0);
            while (i < remaining_nodes) {
                NodeID node(CaR.remainder[i]);
                if (to_remove[node]) {
                    remaining_nodes--;
                    CaR.remainder[i] = CaR.remainder[remaining_nodes];
                    CaR.remainder[remaining_nodes] = node;
                } else {
                    i++;
                }
            }

            CaR.remainder.resize(remaining_nodes);
        }

        std::vector<NodeID> _chooseIndependentSetFromRemainder() {
            typename EdgeDiffPrioritizer<GraphT, CHConstructorT>::CompInOutProduct ciop(this->_base_graph);
            std::sort(CaR.remainder.begin(), CaR.remainder.end(), ciop);
            auto independent_set(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc.calcIndependentSet(CaR.remainder));
            auto edge_diffs(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc.calcEdgeDiffs(independent_set));

            double edge_diff_mean(0);
            for (size_t i(0); i < edge_diffs.size(); i++) {
                edge_diff_mean += edge_diffs[i];
            }
            edge_diff_mean /= independent_set.size();

            std::vector<NodeID> low_edge_diff_nodes;
            for (size_t i(0); i < independent_set.size(); i++) {
                if (edge_diffs[i] <= edge_diff_mean) {
                    NodeID node(independent_set[i]);
                    low_edge_diff_nodes.push_back(node);
                }
            }

            return low_edge_diff_nodes;
        }
        
        bool isDeadEnd(Chain &chain) {
            NodeID front = chain.node_ids.front();
            if (_base_graph.nodeNeighbours(front).size()==1) {
                return true;
            }
            NodeID back = chain.node_ids.back();
            if (_base_graph.nodeNeighbours(back).size()==1) {
                return true;
            }
            return false;            
        }
        
    public:
        DPPrioritizer(GraphT const& base_graph, CHConstructorT const& chc)
        : EdgeDiffPrioritizer<GraphT, CHConstructorT>(base_graph, chc), state(State::start),
        grid(1000, base_graph), chaindetector(base_graph), dp(this->_base_graph, grid),
        CaR(), priolists(), epsilon(0.0001), roundcounter(1) {
        }

        std::vector<NodeID> extractNextNodes() {
            bool empty = EdgeDiffPrioritizer<GraphT, CHConstructorT>::_prio_vec.empty();
            assert(!empty);
            
            std::vector<NodeID> next_nodes;
            
            switch (state) {
                case State::start:
                    Print("Detecting chains");
                    CaR = chaindetector.detectChains(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_prio_vec);
                    Print("Number of chains: " << CaR.chains.size());

                    //test
                    uint counter = 0;
                    for (Chain &chain : CaR.chains) {
                        counter += chain.node_ids.size();
                    }
                    counter += CaR.remainder.size();
                    //size_t priovecsize = this->_prio_vec.size();
                    debug_assert(counter == this->_prio_vec.size());

                    priolists.clear();
                    for (Chain &chain : CaR.chains) {
                        //big chains are generalized
                        if (chain.node_ids.size() >= 3) {                            
                            priolists.push_back(dp.process(chain));                        
                        } else { //small chains are assigned to the remainder
                            for (NodeID node_id : chain.node_ids) {
                                CaR.remainder.push_back(node_id);
                            }
                        }
                    }

                    state = State::removingChains;
                    break;
                    
                case State::removingDeadEnds:                    
                    break;
                    
                case State::removingChains:
                    //calc independent set of nodes in priolists
                    std::vector<bool> marked(this->_base_graph.getNrOfNodes(), false);
                    //extraction from chains
                    for (DP::PrioList &priolist : priolists) {
                        for (auto it = priolist.begin(); it != priolist.end();) {
                            if (epsilon > it->perpendicularLength) {
                                NodeID node_id = it->node_id;
                                if (!marked[node_id]) {

                                    //mark neighbors
                                    EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc._markNeighbours(node_id, marked);

                                    it = priolist.erase(it);
                                    next_nodes.push_back(node_id);
                                } else {
                                    break;
                                }
                            } else {
                                break;
                            }
                        }
                    }
                    break;
                case State::removingRemaining:

                    break;
                case State::done:
                    
                    break;
                default:
                    break;
            }

            
            

            if ((roundcounter - 1) % 5 == 0) {
                epsilon *= 1.5;
                


                
            }
            roundcounter++;



            //independent set from remainder
            std::vector<NodeID> next_nodes(_chooseIndependentSetFromRemainder());
            //remove from remainder
            _removeFromRemainder(next_nodes);
            //std::vector<NodeID> next_nodes;


            //calc independent set of nodes in priolists
            std::vector<bool> marked(this->_base_graph.getNrOfNodes(), false);
            //extraction from chains
            for (DP::PrioList &priolist : priolists) {                
                for (auto it = priolist.begin(); it != priolist.end();) {
                    if (epsilon > it->perpendicularLength) {
                        NodeID node_id = it->node_id;
                        if (!marked[node_id]) {

                            //mark neighbors
                            EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chc._markNeighbours(node_id, marked);

                            it = priolist.erase(it);
                            next_nodes.push_back(node_id);
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
            }


            //auto next_nodes(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_chooseIndependentSet());

            EdgeDiffPrioritizer<GraphT, CHConstructorT>::_remove(next_nodes); //remove from priovector

            return next_nodes;
        }

    };
}