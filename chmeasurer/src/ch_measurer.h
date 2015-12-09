#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "chgraph.h"
#include "geoFunctions.h"
#include "rayCaster.h"
#include "ILP/lineSimplificationILP.h"
#include "ILP/parallelLineSimplficationILP.h"
#include "chainDetector.h"
#include "chains.h"
#include "4Dgrid.h"
#include "CGAL/self_intersection_checker.h"
#include "dijkstra.h"
#include "errorCount.h"
#include "discreteFrechet.h"
#include "dijkstra.h"
#include "track_time.h"
#include "CGAL/cdthp_cross.h"
#include "measure_options.h"
#include "simple_measures/angular_change.h"
#include "simple_measures/regularity.h"
#include "simple_measures/consistency.h"


#include <chrono>
#include <queue>
#include <mutex>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <limits>

namespace chm {
    
    namespace unit_tests
    {
	void testLineSimplfication();
        void testParallelLineSimplfication();    
        void testCDTHCross();
    }

    class CHMeasurer {
    private:
        CHGraph<CHNode, CHEdge> &graph;
        Grid<CHGraph<CHNode, CHEdge> > grid;        
        RangeTree rangeTree;
        ChainDetector<CHGraph<CHNode, CHEdge> > chaindetector;
        FourDGrid<CHGraph<CHNode, CHEdge> > fourDGrid;
        SelfIntersectionChecker selfIntersectionChecker;
        
        
        void _make_chain_measure(Chain &chain, int streettype, ErrorCounts &error_counts) {
            //only big chains are measured
            if (chain.size() >= 3) {    
                AngularChange ac(graph);
                Regularity rg(graph);
                ErrorConsistency ec(graph);

                RedetectedChain redetected_edge_chain = chaindetector.redetect(chain, streettype);
                if (!redetected_edge_chain.remaining_chain.empty()) {
                    Chain expandedChain = graph.expandEdgeChain(redetected_edge_chain);
                    
                    if (chainsAreEqual(expandedChain, redetected_edge_chain.hull_chain)) {
                        double chain_length = graph.calcChainLength(redetected_edge_chain.remaining_chain);
                        
                        //weight with length of redetected(visible) sum
                        error_counts.weighted_regularity_sum +=
                                chain_length * rg.calcStandardDeviation(redetected_edge_chain.remaining_chain, chain_length);
                        error_counts.weighted_variance_sum += chain_length * ec.calcMAD(redetected_edge_chain, chain_length);
                        error_counts.chain_length_sum += chain_length;
                        
                        error_counts.nofNodesWithAngelChange += redetected_edge_chain.remaining_chain.size() - 2;
                        error_counts.angular_change_sum += ac.calcACsum(redetected_edge_chain.remaining_chain, expandedChain);
                        error_counts.chainRedetectedCounter++;
                    } else {
                        error_counts.chainNotRedetectedCounter++;
                    }
                    error_counts.chainNotRedetectedCounter++;
                }
            }
        }
        
        void make_chain_measure(Chains_and_Remainder &CaR) {                        
            ErrorCounts error_counts;                                
            for (uint i = 0; i < CaR.oneWayChainsAccordingToType.size(); i++) {
                ChainsOfType &chainsOfType = CaR.oneWayChainsAccordingToType.at(i);                          
                for (Chain &chain: chainsOfType) {
                    _make_chain_measure(chain, i, error_counts);                    
                }
            }
            for (uint i = 0; i < CaR.twoWayChainsAccordingToType.size(); i++) {
                ChainsOfType &chainsOfType = CaR.twoWayChainsAccordingToType.at(i);                          
                for (Chain &chain: chainsOfType) {  
                    _make_chain_measure(chain, i, error_counts);
                }
            }
            error_counts.print();
        }
        
        void _make_ilp_measure(Chain &chain, int streettype, ErrorCounts &errorcounts) {
            lineSimplificationILP ilp(graph);
            //only big chains are measured
            if (chain.size() >= 3) {        
                //get epsilon error by expanding the used edges
                double epsilon_error = 0;
                RedetectedChain redetected_edge_chain = chaindetector.redetect(chain, streettype);
                if (!redetected_edge_chain.remaining_chain.empty()) {
                    errorcounts.chainRedetectedCounter++;
                    for (EdgeID edge_id: redetected_edge_chain.edges) {
                        double error = graph.calcEdgeError(edge_id);
                        if (error > epsilon_error) {
                            epsilon_error = error;
                        }
                    }

                    //expandedChain should be equal to chain
                    Chain expandedChain = graph.expandEdgeChain(redetected_edge_chain);

                    if(chainsAreEqual(expandedChain, redetected_edge_chain.hull_chain)) {
                        //neither the original chain, nor the simplfication may intersect to use the ilp
                        if(!selfIntersectionChecker.isSelfIntersecting(expandedChain)
                                && !selfIntersectionChecker.isSelfIntersecting(expandedChain) ) {
                            if (epsilon_error > 0) {
                                double ilpNeededNumberOfEdges
                                    = ilp.solve(expandedChain, epsilon_error + std::numeric_limits<double>::epsilon());
                                //assert(redetected_edge_chain.edges.size() >= ilpNeededNumberOfEdges);
                                /*
                                if (!(redetected_edge_chain.edges.size() >= ilpNeededNumberOfEdges)) {

                                    for (NodeID node_id: expandedChain) {
                                        Print(node_id);
                                    }
                                }*/

                                assert(redetected_edge_chain.edges.size() >= ilpNeededNumberOfEdges);
                                double diff = redetected_edge_chain.edges.size() - ilpNeededNumberOfEdges;
                                //Debug("ilpNeededNumberOfEdges:" << ilpNeededNumberOfEdges);
                                //Debug("diff:" << diff);
                                errorcounts.addedDiffs += diff;                                    
                                errorcounts.epsilonerrorGreaterzero++;
                            }else {
                                errorcounts.epsilonerrorEqualszero++;
                            }
                        }else {
                            errorcounts.selfIntersecting++;
                        }       
                    }else {
                        errorcounts.chainNotRedetectedCounter++;
                    }        

                } else {
                    errorcounts.chainNotRedetectedCounter++;
                }          
            }
        }
        
        void make_ilp_measure(Chains_and_Remainder &CaR) {
            Print("make_ilp_measure");  
            
            ErrorCounts errorcounts;                                    
            
            
            for (uint i = 0; i < CaR.oneWayChainsAccordingToType.size(); i++) {
                ChainsOfType &chainsOfType = CaR.oneWayChainsAccordingToType.at(i);            
                for (Chain &chain: chainsOfType) {                                                                
                    _make_ilp_measure(chain, i, errorcounts);                    
                }
            }            
            for (uint i = 0; i < CaR.twoWayChainsAccordingToType.size(); i++) {
                ChainsOfType &chainsOfType = CaR.twoWayChainsAccordingToType.at(i);            
                for (Chain &chain: chainsOfType) {                                                                
                    _make_ilp_measure(chain, i, errorcounts);                    
                }
            }

            errorcounts.print();            
                        
        }
        
        
        
        void make_p_ilp_measure(Chains_and_Remainder &CaR) {
            Print("make_p_ilp_measure");  
            //ParallelLineSimplificationILP p_ilp(graph);
            ParallelLineSimplificationILP p_ilp(graph);
            DiscreteFrechet dF(graph);
            
            ErrorCounts errorcounts;
            
                        
            
            //for (ChainsOfType &chainsOfType: CaR.oneWayChainsAccordingToType) {                                                                
            for (ChainPair chainpair: CaR.chainPairs) {
                
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                //for (Chain &chain: chainsOfType) {                                                                
                
                    //only big chains are measured
                    if (chainpair.chainTo.size() >= 3 && chainpair.chainFrom.size() >= 3) {     
                       
                        //get epsilon error by expanding the used edges
                        double epsilon_error = 0;
                        RedetectedChain redetected_chainTo = chaindetector.redetect(chainpair.chainTo);
                        RedetectedChain redetected_chainFrom = chaindetector.redetect(chainpair.chainFrom);
                        if (!redetected_chainTo.remaining_chain.empty() && !redetected_chainFrom.remaining_chain.empty()) {                            
                            for (EdgeID edge_id: redetected_chainTo.edges) {
                                double error = graph.calcEdgeError(edge_id);
                                if (error > epsilon_error) {
                                    epsilon_error = error;
                                }
                            }
                            for (EdgeID edge_id: redetected_chainFrom.edges) {
                                double error = graph.calcEdgeError(edge_id);
                                if (error > epsilon_error) {
                                    epsilon_error = error;
                                }
                            }
                            
                            
                            //expandedChain should be equal to chain
                            Chain expandedChainTo = graph.expandEdgeChain(redetected_chainTo);                            
                            Chain expandedChainFrom = graph.expandEdgeChain(redetected_chainFrom);/*
                            Chain testchain;
                            testchain.push_back(1);
                            testchain.push_back(2);
                            testchain.push_back(3);
                             * */   
                            
                            
                            if(chainsAreEqual(expandedChainTo, redetected_chainTo.hull_chain)
                                    && chainsAreEqual(expandedChainFrom, redetected_chainFrom.hull_chain)) {
                                errorcounts.chainPairRedetectedCounter++;
                                //neither the original chains, nor the simplfications may intersect to use the ilp
                                if(!selfIntersectionChecker.isSelfIntersecting(expandedChainTo, expandedChainFrom)
                                        && !selfIntersectionChecker.isSelfIntersecting(redetected_chainTo.remaining_chain,
                                                                                        redetected_chainFrom.remaining_chain)) {
                                    double eta1 = dF.calc_dF(chainpair.chainTo, chainpair.chainFrom);
                                    double eta2 = dF.calc_dF(redetected_chainTo.remaining_chain, redetected_chainFrom.remaining_chain);
                                    double eta = std::max(eta1, eta2);
                                    double p_ilpNeededNumberOfEdges
                                        = p_ilp.solve(expandedChainTo, expandedChainFrom,
                                            epsilon_error + std::numeric_limits<double>::epsilon(),
                                            eta + std::numeric_limits<double>::epsilon());
                                    size_t redetectedSize = redetected_chainTo.edges.size() + redetected_chainFrom.edges.size();

                                    /*
                                    if (!(redetectedSize >= p_ilpNeededNumberOfEdges)) {
                                        Print("difference for p_ilp" << redetectedSize -p_ilpNeededNumberOfEdges);
                                        for (NodeID node_id: expandedChainTo) {
                                            Print(node_id);
                                        }
                                        Print(" ");
                                        for (NodeID node_id: expandedChainFrom) {
                                            Print(node_id);
                                        }
                                        Print("----");
                                        for (NodeID node_id: redetected_edge_chainTo.chain) {
                                            Print(node_id);
                                        }
                                        Print(" ");
                                        for (NodeID node_id: redetected_edge_chainFrom.chain) {
                                            Print(node_id);
                                        }

                                        int i = 3;
                                    }*/
                                    assert(redetectedSize >= p_ilpNeededNumberOfEdges);
                                    double diff = redetectedSize - p_ilpNeededNumberOfEdges;
                                    //Debug("ilpNeededNumberOfEdges:" << ilpNeededNumberOfEdges);
                                    //Debug("diff:" << diff);
                                    
                                    errorcounts.addedDiffs += diff;                                    

                                } else {
                                    errorcounts.selfIntersecting++;
                                }
                            } else {
                                errorcounts.chainNotRedetectedCounter++;
                            }
                        } else {
                            errorcounts.chainNotRedetectedCounter++;
                        }
                    }
                //}
            }

            errorcounts.print();
            
            //SelfIntersectionChecker selfIntersectionChecker;
            
            /*
            lineSimplificationILP ilp(graph);
            if (chain.size()>2 && epsilon_error > 0) {
                        ilp.solve(chain, epsilon_error);
            }
             * */
        }
        
        void makeDijkstraMeasure() {            
            Print("make_dijkstra_measure");  
            //random Dijkstras
            std::string queryFile = std::to_string(graph.getNrOfNodes()) + ".txt";
        
            std::string line;
            std::ifstream file;
            std::vector<NodeID> srcList;
            std::vector<NodeID> tgtList;
            file.open(queryFile.c_str(), std::ios::in);
            uint queryCount = 0;
            if (file.is_open()) { //check if there exist a query file
                file.seekg(0, std::ios::beg);

                getline(file, line, '\n');
                queryCount = (uint) atoi(line.c_str());

                for (uint i = 0; i < queryCount; i++) {
                    getline(file, line, '\n');
                    std::stringstream ss(line);
                    std::string buffer;

                    ss >> buffer;
                    srcList.push_back(atoi(buffer.c_str()));

                    ss >> buffer;
                    tgtList.push_back(atoi(buffer.c_str()));
                }            
            } else { //create one 
                std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
                std::uniform_int_distribution<uint> dist(0,graph.getNrOfNodes()-1);
                auto rand_node = std::bind (dist, gen);
                ofstream outFileStream;
                outFileStream.open(queryFile);
                queryCount = (graph.getNrOfNodes()-1)/10000;
                outFileStream << queryCount << std::endl;
                for (uint i = 0; i < queryCount; i++) {
                    NodeID src = rand_node();
                    NodeID tgt = rand_node();
                    srcList.push_back(src);
                    tgtList.push_back(tgt);
                    outFileStream << src << " " << tgt << std::endl;
                }
                outFileStream.close();
            }
            
            CHDijkstra<CHNode, CHEdge> chdij(graph);
            TrackTime tt =VerboseTrackTime();
            tt.track("Starting Dijkstras");
            for (uint i = 0; i < queryCount; i++) {
                std::vector<EdgeID> path;
                chdij.calcShopa(srcList.at(i), tgtList.at(i), path);
            }
            tt.track("Finished Dijkstras");
            
            
        }
        
        void makeEdgeMeasure() {
            Print("makeEdgeMeasure");
            Print("Accumulating Error and Crossings Counting");
            CDTHPCross cdthpC(graph, grid, rangeTree);
            //Grid grid = Grid(1, graph);
            double getNofCrossings = 0;
            double accumulatedError = 0;
            //RayCaster raycaster(graph, grid);       
            for (uint edge_id = 0; edge_id < graph.getNrOfEdges(); edge_id++) {
                
                if (graph.isValidEdge(edge_id) && graph.isShortcut(edge_id)) {
                                        
                    Chain chain = graph.getCenterNodes(edge_id);
                    chain.push_front(graph.getEdge(edge_id).src);
                    chain.push_back(graph.getEdge(edge_id).tgt);                    
                                        
                    //getNofCrossings += raycaster.getNofCrossings(chain);
                    getNofCrossings += cdthpC.getNofCrossings(chain);
                    double epsilon_error = graph.calcEdgeError(edge_id);
                    accumulatedError += epsilon_error;                        
                }
                
            }
            Print("number of Crossings: " << getNofCrossings);
            Print("accumulated error: " << accumulatedError);
        }
        
        
    public:

        CHMeasurer(CHGraph<CHNode, CHEdge> &graph):
            graph(graph), grid(1000, graph), rangeTree(graph), chaindetector(graph), fourDGrid(10, graph), selfIntersectionChecker(graph) {
            
        }
            
        void makeMeasurement(MeasureOptions m_options) {
            Print("make measure");    
            //measurements without path-finding info
            graph.zoom(100, false, 0);
            for (uint i = 0; i < graph.getNrOfNodes(); i++) {
                assert(graph.isValidNode(i));
            }

            Chains_and_Remainder CaR;
            
            
            Print("Detecting chains");            
            //CaR = chaindetector.detectChains(EdgeDiffPrioritizer<GraphT, CHConstructorT>::_prio_vec);  
            std::vector<NodeID> activeNodeIDs;
            for (uint i = 0; i < graph.getNrOfNodes(); i++) {                
                if(graph.isValidNode(i)) {
                    activeNodeIDs.push_back(i);                
                }
            }           

            CaR = chaindetector.detectChains(activeNodeIDs);
            Print("Number of chains: " << CaR.getNrOfChains());
            
            Print("IdentifyingChainPairs ");                                      
            fourDGrid.identifyPairs(CaR);
            Print("Number of chain pairs: " << CaR.chainPairs.size());
            
            
            Print("Zooming ");
            graph.zoom(50, false, 0);               
            
            //makeEdgeMeasure();
            
            make_chain_measure(CaR);
            
            if (m_options.dijkstra) {
                makeDijkstraMeasure();
            }
            if (m_options.ilp) {
                make_ilp_measure(CaR);
            }
            if (m_options.p_ilp) {
                make_p_ilp_measure(CaR);
            }
           
        }

    };
}
