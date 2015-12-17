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
#include "simplification/discreteCurveEvolution.h"

#include "chains.h"
#include "grid.h"
#include "chainDetector.h"

#include "indexed_container.h"
#include "dead_end_detector.h"

#include "4Dgrid.h"

#include "simplification/prio_nodes.h"
#include "nodes_and_edges.h"
#include "s_options.h"


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
        
        struct CompWeight {
            const std::vector<double> &comp_weights;

            CompWeight(const std::vector<double> &comp_weights)
                    : comp_weights(comp_weights) {}

            bool operator()(NodeID node1, NodeID node2) const
            {                    
                    return comp_weights[node1] < comp_weights[node2];
            }                        
        }; 
        
        struct CompMST {
            GraphT const& g;

            CompMST(GraphT const& g)
                    : g(g) {}

            bool operator()(NodeID node1, NodeID node2) const
            {               
                return g.getMaxStreetType(node1) > g.getMaxStreetType(node1);                    
            }                        
        }; 
        
        
        
        /*
        struct CompWithTable {
            const std::vector<int> &table;
            const std::vector<int> &map;

            CompWithTable(const std::vector<int> &table, const std::vector<int> &map)
                    : table(table), map(map) {}

            bool operator()(NodeID node1, NodeID node2) const
            {                    
                    return table[map[node1]] < table[map[node2]];
            }                        
        };*/

	GraphT const& _base_graph;
        CHConstructorT const& _chc;
        std::vector<NodeID> _prio_vec;
        std::vector<double> weights;
                        
        State state;
        Grid<GraphT> grid;
        FourDGrid<GraphT> fourDGrid;
        
        DeadEndDetector<GraphT> deadEndDetector;
        ChainDetector<GraphT> chaindetector;                
                
        std::unique_ptr<ls::LineSimplifier> lineSimplifier;         

        std::list<Chain> deadEnds;  
        
        //_prio_vec = chains + remainder        
        Chains_and_Remainder CaR;
        
        typedef std::list<ls::simplePrioNode> PrioList;
        struct ChainPriolist{
            PrioList priolist;
            std::list<NodeID> outer_nodes;
            ChainPriolist(PrioList &&priolist, const Chain &chain): priolist(priolist) {
                outer_nodes.push_back(chain.front());
                outer_nodes.push_back(chain.back());
            }
            ChainPriolist(PrioList &&priolist, const ChainPair &chain_pair): priolist(priolist) {
                outer_nodes.push_back(chain_pair.chainTo.front());
                outer_nodes.push_back(chain_pair.chainTo.back());
                outer_nodes.push_back(chain_pair.chainFrom.front());
                outer_nodes.push_back(chain_pair.chainFrom.back());
            }
            ChainPriolist(std::list<NodeID> outer_nodes): outer_nodes(outer_nodes){}
        };
        std::list<ChainPriolist> chain_prio_lists;
        std::list<ChainPriolist> deadEndPrioLists;
        //std::list<std::list<ls::simplePrioNode>> priolists;

        double epsilon;
        int roundcounter = 1;
        uint deadEndPhaseCounter = 0;
        const SOptions s_options;
        
        void init(std::vector<NodeID>& node_ids) {
            _prio_vec = std::move(node_ids);
        }        
        
        void _removeFrom(std::vector<NodeID> const& nodes, std::vector<NodeID> &from_nodes) {
            std::vector<bool> to_remove(this->_base_graph.getNrOfNodes(), false);
            for (auto node : nodes) {     
                debug_assert(0 <= node && node < to_remove.size());
                to_remove[node] = true;                
            }

            size_t remaining_nodes(from_nodes.size());
            size_t i(0);
            while (i < remaining_nodes) {                
                debug_assert(0 <= i && i < from_nodes.size());
                NodeID node(from_nodes[i]);                
                if (to_remove[node]) {
                    remaining_nodes--;
                    debug_assert(0 <= remaining_nodes && remaining_nodes < from_nodes.size());
                    from_nodes[i] = from_nodes[remaining_nodes];
                    from_nodes[remaining_nodes] = node;
                } else {
                    i++;
                }
            }

            from_nodes.resize(remaining_nodes);
        }

        std::vector<NodeID> _chooseIndependentSetFromRemainderWeights(std::vector<NodeID> &remainder) {
            std::sort(remainder.begin(), remainder.end(), CompInOutProduct(_base_graph));        
            auto independent_set(_chc.calcIndependentSet(remainder));
            //auto edge_diffs(_chc.calcEdgeDiffs(independent_set));
            //auto edge_diffs(_chc.calcWeightedEdgeDiffs(independent_set));
            auto edge_diffs(_chc.calcGeoImportance(independent_set));

            double edge_diff_mean(0);
            for (size_t i(0); i<edge_diffs.size(); i++) {
                    edge_diff_mean += edge_diffs[i];
            }
            edge_diff_mean /= independent_set.size();

            std::vector<NodeID> low_edge_diff_nodes;
            for (size_t i(0); i<independent_set.size(); i++) {
                    //if (edge_diffs[i] <= edge_diff_mean) {
                if (edge_diffs[i] <= edge_diff_mean) {
                            NodeID node(independent_set[i]);
                            low_edge_diff_nodes.push_back(node);
                    }
            }
            //low_edge_diff_nodes.resize((low_edge_diff_nodes.size()/4) +1);
            return low_edge_diff_nodes;
        }  
        
        std::vector<NodeID> _chooseIndependentSetFromRemainderEDTable(std::vector<NodeID> &remainder) {
            std::sort(remainder.begin(), remainder.end(), CompInOutProduct(_base_graph));        
            auto independent_set(_chc.calcIndependentSet(remainder));
            //auto edge_diffs(_chc.calcEdgeDiffs(independent_set));
            //auto edge_diffs(_chc.calcWeightedEdgeDiffs(independent_set));
            auto edge_diffs(_chc.calcEdgeDiffs(independent_set));

            double edge_diff_mean(0);
            for (size_t i(0); i<edge_diffs.size(); i++) {
                    edge_diff_mean += edge_diffs[i];
            }
            edge_diff_mean /= independent_set.size();

            std::vector<NodeID> low_edge_diff_nodes;
            for (size_t i(0); i<independent_set.size(); i++) {
                    //if (edge_diffs[i] <= edge_diff_mean) {
                if (edge_diffs[i] <= edge_diff_mean) {
                    NodeID node(independent_set[i]);
                    low_edge_diff_nodes.push_back(node);
                }
            }
            
            //std::sort(low_edge_diff_nodes.begin(), low_edge_diff_nodes.end(), CompWithTable(edge_diffs, map));
            
            low_edge_diff_nodes.resize((low_edge_diff_nodes.size()/4) +1);
            return low_edge_diff_nodes;
        }   
        
        std::vector<NodeID> _chooseIndependentSetFromRemainderED(std::vector<NodeID> &remainder) {
            std::sort(remainder.begin(), remainder.end(), CompInOutProduct(_base_graph));        
            auto independent_set(_chc.calcIndependentSet(remainder));
            //auto edge_diffs(_chc.calcEdgeDiffs(independent_set));
            //auto edge_diffs(_chc.calcWeightedEdgeDiffs(independent_set));
            auto edge_diffs(_chc.calcEdgeDiffs(independent_set));

            double edge_diff_mean(0);
            for (size_t i(0); i<edge_diffs.size(); i++) {
                    edge_diff_mean += edge_diffs[i];
            }
            edge_diff_mean /= independent_set.size();

            std::vector<NodeID> low_edge_diff_nodes;
            for (size_t i(0); i<independent_set.size(); i++) {
                    //if (edge_diffs[i] <= edge_diff_mean) {
                if (edge_diffs[i] <= edge_diff_mean) {
                            NodeID node(independent_set[i]);
                            low_edge_diff_nodes.push_back(node);
                    }
            }
            //low_edge_diff_nodes.resize((low_edge_diff_nodes.size()/4) +1);
            return low_edge_diff_nodes;
        } 
        
        std::vector<NodeID> _chooseIndependentSetFromRemainder2(std::vector<NodeID> &remainder) {
            //first step: sort after InOutProduct and get independent Set
            std::sort(remainder.begin(), remainder.end(), CompInOutProduct(_base_graph));        
            auto independent_set(_chc.calcIndependentSet(remainder));
            
            std::vector<NodeWeight> node_weights(_chc.calcGeoImportance2(independent_set));            

            //second step: get geo-unimportant nodes from here
            double geo_measure_mean = 0;
            for (size_t i(0); i<node_weights.size(); i++) {
                geo_measure_mean += node_weights[i].geo_measure;                    
            }
            geo_measure_mean /= independent_set.size();
            
            double edge_diff_mean = 0;
            std::vector<uint> low_geo_measure_nodes; //saves index of independent set
            for (size_t i(0); i<independent_set.size(); i++) {                    
                if (node_weights[i].geo_measure <= geo_measure_mean) {                                        
                    edge_diff_mean += node_weights[i].edge_diff;                   
                    low_geo_measure_nodes.push_back(i);
                }
            }
            edge_diff_mean /= low_geo_measure_nodes.size();
            
            //third step: get nodes with low edge difference
            std::vector<NodeID> low_edge_diff_nodes;
            for (size_t i(0); i<low_geo_measure_nodes.size(); i++) {                    
                uint index = low_geo_measure_nodes[i];
                if (node_weights[index].edge_diff <= edge_diff_mean) {
                    NodeID node_id(independent_set[index]);
                    low_edge_diff_nodes.push_back(node_id);
                }
            }
            
            //low_edge_diff_nodes.resize((low_edge_diff_nodes.size()/4) +1);
            return low_edge_diff_nodes;
        } 
        
        
        //TODO slowdown
        std::vector<NodeID> _chooseIndependentSetFromRemainder(std::vector<NodeID> &remainder) {
            //CompInOutProduct ciop(this->_base_graph);            
            if (remainder.empty()) {
                std::vector<NodeID> emptyList;
                return emptyList;
            } else {
                Print("Calc Weights");
                _chc.calcNodeWeights(remainder, weights);
                std::sort(remainder.begin(), remainder.end(), CompWeight(weights));                
                //auto weightNodes(_chc.calcWeightNodes(CaR.remainder));                
                auto independent_set(_chc.calcIndependentSet(remainder));
                
                //auto edge_diffs(_chc.calcEdgeDiffs(independent_set));
                //auto edge_diffs(_chc.calcWeightedEdgeDiffs(independent_set));
                double weight_mean(0);
                for (NodeID node_id: independent_set) {                    
                    weight_mean += weights[node_id];
                }
                debug_assert(independent_set.size() != 0);
                weight_mean /= independent_set.size();

                std::vector<NodeID> low_weight_nodes;
                for (NodeID node_id: independent_set) {                    
                    if (weights[node_id] <= weight_mean) {                        
                        low_weight_nodes.push_back(node_id);
                    }
                }

                return low_weight_nodes;
            }
        }   
        
        std::vector<NodeID> _chooseIndependentSetFromRemainderMST(std::vector<NodeID> &remainder) {
            std::sort(remainder.begin(), remainder.end(), CompInOutProduct(_base_graph));        
            auto independent_set(_chc.calcIndependentSet(remainder));
            //auto edge_diffs(_chc.calcEdgeDiffs(independent_set));
            //auto edge_diffs(_chc.calcWeightedEdgeDiffs(independent_set));
            auto max_street_types(_chc.calcMaxStreetTypes(independent_set));

            int highest_mst(-3);
            for (size_t i(0); i<max_street_types.size(); i++) {
                if (highest_mst <= max_street_types[i]) {
                    highest_mst = max_street_types[i];
                }
                //highest_mst = std::max(max_street_types[i], highest_mst);
            }
            Print("independent_set.size(): " << independent_set.size());
            Print("highest_mst: " << highest_mst);
            //highest_mst /= independent_set.size();

            std::vector<NodeID> high_mst_nodes;
            for (size_t i(0); i<independent_set.size(); i++) {
                    //if (edge_diffs[i] <= edge_diff_mean) {
                if (max_street_types[i] >= highest_mst) {
                    NodeID node(independent_set[i]);
                    high_mst_nodes.push_back(node);
                }
            }
            //low_edge_diff_nodes.resize((low_edge_diff_nodes.size()/4) +1);
            return high_mst_nodes;
        } 
        
        
        
        void FillChainsInPriolists(std::vector<ChainsOfType> &chainsaccordingToType) {
            for (ChainsOfType &chainsOfType: chainsaccordingToType) {                                                                
                
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                for (Chain &chain: chainsOfType) {                                                                
                    //big chains are generalized
                    if (chain.size() >= 3) {                          
                        chain_prio_lists.push_back(ChainPriolist(lineSimplifier->process(chain, Chain()), chain));
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
            chain_prio_lists.clear();        
            FillChainsInPriolists(CaR.oneWayChainsAccordingToType);
            FillChainsInPriolists(CaR.twoWayChainsAccordingToType);                        
                        
            for (ChainPair &chainPair: CaR.chainPairs) {    
                    if (chainPair.chainTo.size() + chainPair.chainFrom.size() >= 7
                            && chainPair.chainTo.size() >=3 && chainPair.chainFrom.size() >= 3) {                        
                        chain_prio_lists.push_back(ChainPriolist(lineSimplifier->process(chainPair.chainTo, chainPair.chainFrom), chainPair)); 
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
        
        void FillDeadEndPrioLists (std::list<Chain> deadEnds) {
            deadEndPrioLists.clear();
            for (Chain &chain: deadEnds) {
                if (chain.size() >= 3) {
                    deadEndPrioLists.push_back(ChainPriolist(lineSimplifier->process(chain, Chain()), chain));
                } else {                    
                    //pseudo remainder
                    deadEndPrioLists.push_back(ChainPriolist(chain));                    
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
        std::vector<NodeID> removeDeadEnds {
            
        }*/
        
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
        
        //, dp(this->_base_graph, grid), cpdp(this->_base_graph,
        
        DPPrioritizer(SOptions s_options, GraphT const& base_graph, CHConstructorT const& chc)
                : _base_graph(base_graph), _chc(chc), weights(_base_graph.getNrOfNodes()), state(State::start),
                grid(1000, base_graph), fourDGrid(1, base_graph), deadEndDetector(base_graph),
                chaindetector(base_graph), epsilon(10000),        
                s_options(s_options) {                        
                    
//            node_weights.resize();
            lineSimplifier = createLineSimplifier(s_options, s_options.lineSimplifier_type, base_graph, grid);
            assert(lineSimplifier != nullptr);
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
                    //deadEnds = deadEndDetector.detectDeadEnds(_prio_vec);
                    state = State::removingDeadEnds;                    
                }
                case State::removingDeadEnds: {
                    if(deadEndPrioLists.empty()) {                        
                        Print("Detecting Dead Ends");
                        switch (s_options.deadEndDetect_type) {
                            case DeadEndDetectType::NONE: {
                                break;
                            }
                            case DeadEndDetectType::DE: {
                                FillDeadEndPrioLists(deadEndDetector.detectDeadEnds(_prio_vec));
                                break;
                            }
                            case DeadEndDetectType::EDE: {
                                FillDeadEndPrioLists(deadEndDetector.detectExtendedDeadEnds(_prio_vec));
                                break;
                            }
                        }                        
                        deadEndPhaseCounter++;
                    }
                                        
                    Print("Removing Dead Ends");
                    if (!deadEndPrioLists.empty()) {
                        std::vector<bool> marked(this->_base_graph.getNrOfNodes(), false);    
                        for (auto p_it = deadEndPrioLists.begin(); p_it != deadEndPrioLists.end();) {
                            PrioList &priolist = p_it->priolist;
                            if (!priolist.empty()) {
                                for (auto it = priolist.begin(); it != priolist.end();) {                                    
                                    NodeID node_id = it->node_id;
                                    if (!marked[node_id]) {                                   
                                        _chc._markNeighbours(node_id, marked);
                                        next_nodes.push_back(node_id);
                                        it = priolist.erase(it);                                        
                                    } else {
                                        break;
                                    }                                    
                                }
                                p_it++;
                            } else {
                                //remove remainder of chain
                                if (!p_it->outer_nodes.empty()) {
                                    NodeID node_id = p_it->outer_nodes.front();
                                    if (!marked[node_id]) {  
                                        _chc._markNeighbours(node_id, marked);
                                        next_nodes.push_back(node_id);
                                        p_it->outer_nodes.pop_front();
                                    }
                                    p_it++;                                    
                                } else {
                                    p_it = deadEndPrioLists.erase(p_it);
                                }
                            }                        
                        }
                        //dont contract peanuts
                        double threshold = log(_base_graph.getNrOfNodes()) / log (1.01);
                        if(next_nodes.size() < threshold) {
                            deadEndPrioLists.clear();
                            deadEndPhaseCounter++;
                            if (deadEndPhaseCounter > 5) {
                                state = State::removingChains;  
                            }

                        }
                        break;
                    } else {
                        state = State::removingChains;                        
                    }
                    
                    //if(deadEndPrioLists.empty && deadEndPhaseCounter < 3)
                        /*
                    next_nodes = removeDeadEnds(deadEnds);
                    if (!next_nodes.empty()) {
                        break;
                    } else {
                        state = State::removingChains;
                        break;
                    }*/    
                }
                case State::removingChains: {                    
                    
                    if ((roundcounter - 1) % 5 == 0) {
                        if (epsilon < 100000) { //prevents epsilon overflow
                            epsilon *= 2.0;
                        }

                        Print("Detecting chains");
                        
                        //CaR = chaindetector.detectChains(_prio_vec);                        
                        CaR.remainder = _prio_vec;
                        
                        
                        Print("Number of chains: " << CaR.getNrOfChains());
                        debug_assert(CaR.getNrOfNodesInChains() + CaR.remainder.size() == this->_prio_vec.size());
                        if (s_options.pairMatch_type != ls::PairMatchType::NONE) {
                            Print("IdentifyingChainPairs ");
                            fourDGrid.identifyPairs(CaR);
                        }
                        
                        Print("Number of chain pairs: " << CaR.chainPairs.size());

                        Print("FillPriolists");
                        FillPriolists();
                    }

                    Print("Getting Independent set from Remainder");
                    //independent set from remainder
                    //if (roundcounter > 20) {
                    //next_nodes = _chooseIndependentSetFromRemainderMST(CaR.remainder);
                    
                    //next_nodes = _chooseIndependentSetFromRemainderEDTable(CaR.remainder);
                    next_nodes = _chooseIndependentSetFromRemainder2(CaR.remainder);
                    //next_nodes = _chooseIndependentSetFromRemainderED(CaR.remainder);
                    //}
                    //remove from remainder
                    //_removeFromRemainder(next_nodes);
                    
                    _removeFrom(next_nodes, CaR.remainder);
                    
                    
                    Print("Getting Independent set from Chains");
                    //calc independent set of nodes from priolists
                    std::vector<bool> marked(this->_base_graph.getNrOfNodes(), false);                    
                    //extraction from chains
                    //for (std::list<ls::simplePrioNode> &priolist : priolists) {
                    for (auto p_it = chain_prio_lists.begin(); p_it != chain_prio_lists.end();) {
                        PrioList &priolist = p_it->priolist;
                        if (!priolist.empty()) {
                            for (auto it = priolist.begin(); it != priolist.end();) {
                                if (epsilon > it->perpendicularLength) {
                                    NodeID node_id = it->node_id;
                                    if (!marked[node_id]) {                                   
                                        _chc._markNeighbours(node_id, marked);
                                        next_nodes.push_back(node_id);
                                        it = priolist.erase(it);
                                        break; //DEBUG
                                    } else {
                                        break;
                                    }
                                } else {
                                    break;
                                }
                            }
                            p_it++;
                        } else {
                            //outer nodes are free to be contracted now
                            for (NodeID outer_node: p_it->outer_nodes) {
                                CaR.remainder.push_back(outer_node);
                            }
                            p_it = chain_prio_lists.erase(p_it);
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
            
        _removeFrom(next_nodes, _prio_vec);        
        //_remove(next_nodes); //remove from priovector               

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