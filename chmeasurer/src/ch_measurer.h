#pragma once


#include "defs.h"
#include "nodes_and_edges.h"
#include "chgraph.h"
#include "geoFunctions.h"
//#include "rayCaster.h"

#include "ILP/lineSimplificationILP.h"
#include "ILP/parallelLineSimplficationILP.h"
#include "ILP/calc_frechet.h"

#include "chainDetector.h"
#include "chains.h"
#include "4Dgrid.h"

#include "dijkstra.h"
#include "errorCount.h"
#include "discreteFrechet.h"
#include "dijkstra.h"
#include "track_time.h"

#include "CGAL/cdthp_cross.h"
#include "CGAL/self_intersection_checker.h"

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
    
    

    class CHMeasurer {
    private:
        CHGraph<CHNode, CHEdge> &graph;
        Grid<CHGraph<CHNode, CHEdge> > grid;                
        ChainDetector<CHGraph<CHNode, CHEdge> > chaindetector;
        FourDGrid<CHGraph<CHNode, CHEdge> > fourDGrid;
        ErrorCounts sadf;
        OutData out_data;       
       
        //SelfIntersectionChecker selfIntersectionChecker;
        
        /*
        void _make_chain_measureOLD(Chain &chain, int streettype, ErrorCounts &error_counts) {
            //only big chains are measured
            if (chain.size() >= 3) {    
                AngularChange ac(graph);
                Regularity rg(graph);
                ErrorConsistency ec(graph);

                RedetectedChain redetected_edge_chain = chaindetector.redetect(chain, streettype);
                if (!redetected_edge_chain.remaining_chain.empty()) {
                    Chain expandedChain = chains::expandEdgeChain(redetected_edge_chain, graph);
                    
                    if (chains::areEqual(expandedChain, redetected_edge_chain.hull_chain)) {
                        double chain_length = chains::calcChainGeoLength(redetected_edge_chain.remaining_chain, graph);
                        
                        //weight with length of redetected(visible) sum
                        error_counts.weighted_regularity_sum +=
                                chain_length * rg.calcStandardDeviation(redetected_edge_chain.remaining_chain, chain_length);
                        error_counts.weighted_variance_sum += chain_length * ec.calcMAD(redetected_edge_chain, chain_length);
                        error_counts.length_sum += chain_length;
                        
                        error_counts.nofNodesWithAngelChange += redetected_edge_chain.remaining_chain.size() - 2;
                        error_counts.angular_change_sum += ac.calcACsum(redetected_edge_chain.remaining_chain, expandedChain);
                        error_counts.chainRedetectedCounter++;
                    } else {
                        error_counts.chainNotRedetectedCounter++;
                    }
                    error_counts.chainNotRedetectedCounter++;
                }
            }
        }*/
        
        void _make_chain_measure(const Chain &chain, int streettype, ErrorCounts &error_counts) {
            //only big chains are measured
            if (chain.size() >= 3) {                    
                
                EdgeChain edge_chain = chains::getChainEdges(chain, graph);
                CDTHPCross cdthpC(graph, grid);                                              
                
                for (EdgeID edge_id: edge_chain) {
                    
                    if (graph.isValidEdge(edge_id) && graph.isShortcut(edge_id)) {

                        Chain center_nodes = graph.getCenterNodes(edge_id);
                        center_nodes.push_front(graph.getEdge(edge_id).src);
                        center_nodes.push_back(graph.getEdge(edge_id).tgt);                    
               
                        error_counts.nofCrossings += cdthpC.getNofCrossings(center_nodes);                                                
                    }
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
            for (const ChainPair& chainpair: CaR.chainPairs) {
                //TODO Streettype needs to be removed 
                _make_chain_measure(chainpair.chainTo, 1, error_counts);                   
                _make_chain_measure(chainpair.chainFrom, 1, error_counts);                   
            }
            out_data.crossings = error_counts.nofCrossings;
            //error_counts.print();
                        
            //Print("length: " << error_counts.length_sum);
            Print("----------------------------");
            Print("number of Crossings: " << error_counts.nofCrossings);
            Print("----------------------------");
            //Print("number of Crossings/length: " << nofCrossings/error_counts.length_sum);
            Print("number of Crossings * NrOfValidEdges: " << error_counts.nofCrossings * graph.getNrOfValidEdges());
            //Print("accumulated error/length: " << accumulatedError/error_counts.length_sum);
        }
        
        
        void _make_ilp_measure(const Chain &chain, int streettype, ErrorCounts &errorcounts) {                       
            //only big chains are measured
            if (chain.size() >= 3) {
                //EdgeChain edge_chain = chaindetector.getChainEdges(chain);
                EdgeChain edge_chain = chains::getChainEdges(chain, graph);

                for (EdgeChain split_edge_chain : chains::split(edge_chain, graph)) {
                    if (split_edge_chain.size() >= 2) {        
                        //get epsilon error by expanding the used edges
                        double epsilon_error = 0;                                                    
                        for (EdgeID edge_id: split_edge_chain) {
                            double error = graph.calcEdgeError(edge_id);
                            if (error > epsilon_error) {
                                epsilon_error = error;
                            }
                        }                    
                        Chain split_chain = chains::toNodeChain(split_edge_chain, graph);
                        Chain expanded_split_chain = chains::toExpandedNodeChain(split_edge_chain, graph);  
                        //Print("split_chain.size()" << split_chain.size());
                        //Print("expanded_split_chain.size()" << expanded_split_chain.size());

                        if(chains::uniqueLocations(expanded_split_chain, Chain(), graph)) { //can happen by having same centernode or nodes in the same place
                            assert(chains::uniqueLocations(split_chain, Chain(), graph)); //because they are subsets of expanded
                            SelfIntersectionChecker selfIntersectionChecker(graph);
                            if(!selfIntersectionChecker.isSelfIntersecting(split_chain)
                                    && !selfIntersectionChecker.isSelfIntersecting(expanded_split_chain)) {                           

                                lineSimplificationILP ilp(graph);
                                double ilpNeededNumberOfEdges
                                    = ilp.solve(expanded_split_chain, epsilon_error + std::numeric_limits<double>::epsilon() * geo::R);                    

                                int nof_edges = split_chain.size()-1;
                                assert(nof_edges >= ilpNeededNumberOfEdges);
                                double diff = nof_edges - ilpNeededNumberOfEdges;
                                //Debug("ilpNeededNumberOfEdges:" << ilpNeededNumberOfEdges);
                                //Debug("diff:" << diff);
                                double chain_length = chains::calcChainGeoLength(expanded_split_chain, graph);
                                errorcounts.weighted_epsilon += epsilon_error * chain_length;
                                errorcounts.length_sum += chain_length;
                                errorcounts.addedDiffs += diff;                                                   

                            }else {
                                errorcounts.selfIntersecting++;
                            }
                            errorcounts.nofAnalyzedPolygons++;
                        } else {
                            errorcounts.sameLocation++;
                        }
                    }                    
                }
            }
        }
        
        
        
        void make_ilp_measure(Chains_and_Remainder &CaR) {
            Print("make_ilp_measure");              
            ErrorCounts errorcounts;                                                
            
            for (uint i = 0; i < CaR.oneWayChainsAccordingToType.size(); i++) {
                const ChainsOfType &chainsOfType = CaR.oneWayChainsAccordingToType.at(i);                            
                for (const Chain &chain: chainsOfType) {                     
                    _make_ilp_measure(chain, i, errorcounts);                    
                }
            }            
            for (uint i = 0; i < CaR.twoWayChainsAccordingToType.size(); i++) {
                const ChainsOfType &chainsOfType = CaR.twoWayChainsAccordingToType.at(i);            
                for (const Chain &chain: chainsOfType) {                                                                                    
                    _make_ilp_measure(chain, i, errorcounts);                   
                }              
            }
            for (const ChainPair& chainpair: CaR.chainPairs) {
                //TODO Streettype needs to be removed 
                _make_ilp_measure(chainpair.chainTo, 1, errorcounts);                   
                _make_ilp_measure(chainpair.chainFrom, 1, errorcounts);                   
            }

            
            errorcounts.print(); 
            Print("length" << errorcounts.length_sum);
            //assert(errorcounts.length_sum != 0);
            //Print("avg_epsilon: " << errorcounts.weighted_epsilon/errorcounts.length_sum);
            Print("-------------------");
            Print("diff/length ratio: " << errorcounts.addedDiffs/errorcounts.length_sum);
            Print("self intersecting ratio: " << (double) errorcounts.selfIntersecting/ (double) errorcounts.nofAnalyzedPolygons );
            Print("-------------------");
            Print(" ");
            
            out_data.dist_Plslash = errorcounts.length_sum;
            out_data.deltapslash = errorcounts.addedDiffs/errorcounts.length_sum;
            out_data.intersection_plash = (double) errorcounts.selfIntersecting/ (double) errorcounts.nofAnalyzedPolygons;                                    
        }                
        
        
        void make_p_ilp_measure(Chains_and_Remainder &CaR) {
            Print("make_p_ilp_measure");              
            //CalcFrechetILP cf_ilp(graph);            
            
            ErrorCounts errorcounts;
            
            //for (ChainsOfType &chainsOfType: CaR.oneWayChainsAccordingToType) {                                                                
            for (ChainPair chainpair: CaR.chainPairs) {
                
                
                
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                //for (Chain &chain: chainsOfType) {                                                                
                    
                    //only big chains are measured
                    if (chainpair.chainTo.size() >= 3 && chainpair.chainFrom.size() >= 3) {
                        
                        DiscreteFrechet dF(graph);
                        double eta_whole = dF.calc_dF(chainpair.chainTo, chainpair.chainFrom);
                        EdgeChain w_edge_chain_to = chains::getChainEdges(chainpair.chainTo, graph);
                        EdgeChain w_edge_chain_from = chains::getChainEdges(chainpair.chainFrom, graph);
                        Chain w_expandedChainTo = chains::toExpandedNodeChain(w_edge_chain_to, graph);                            
                        Chain w_expandedChainFrom = chains::toExpandedNodeChain(w_edge_chain_from, graph);
                        double combined_whole_chain_length = chains::calcChainGeoLength(w_expandedChainTo, graph)
                                                            + chains::calcChainGeoLength(w_expandedChainFrom, graph);                                    
                        errorcounts.weighted_eta2 += eta_whole * combined_whole_chain_length;
                        errorcounts.length_sum2 += combined_whole_chain_length;
                        
                        EdgeChainPair edge_chain_pair(chainpair, graph);                                                
                        
                        //process the splitted parts of the chains
                        for (EdgeChainPair split_edge_chain_pair : chains::split(edge_chain_pair, graph)) {
                            if (split_edge_chain_pair.chainTo.size() >= 2 && split_edge_chain_pair.chainFrom.size() >= 2) { 
                                double epsilon_error = 0;
                                for (EdgeID edge_id: split_edge_chain_pair.chainTo) {
                                    double error = graph.calcEdgeError(edge_id);
                                    if (error > epsilon_error) {
                                        epsilon_error = error;
                                    }
                                }
                                for (EdgeID edge_id: split_edge_chain_pair.chainFrom) {
                                    double error = graph.calcEdgeError(edge_id);
                                    if (error > epsilon_error) {
                                        epsilon_error = error;
                                    }
                                }
                                
                                
                                Chain chainTo = chains::toNodeChain(split_edge_chain_pair.chainTo, graph);
                                Chain chainFrom = chains::toNodeChain(split_edge_chain_pair.chainFrom, graph);
                                Chain expandedChainTo = chains::toExpandedNodeChain(split_edge_chain_pair.chainTo, graph);                            
                                Chain expandedChainFrom = chains::toExpandedNodeChain(split_edge_chain_pair.chainFrom, graph);
                                
                                
                                if(chains::uniqueLocations(expandedChainTo, expandedChainFrom, graph)) { //can happen by having same centernode or nodes in the same place
                                    assert(chains::uniqueLocations(chainTo, chainFrom, graph)); //because they are subsets of expanded
                                    SelfIntersectionChecker selfIntersectionChecker(graph);
                                    if(!selfIntersectionChecker.isSelfIntersecting(chainTo, chainFrom)                                      
                                        && !selfIntersectionChecker.isSelfIntersecting(expandedChainTo, expandedChainFrom)) {  

                                        //Print("split_chain.size()" << chainTo.size() + chainFrom.size());
                                        //Print("expanded_split_chain.size()" << expandedChainTo.size() + expandedChainFrom.size());                                                                

                                        DiscreteFrechet dF(graph);
                                        double eta = dF.calc_dF(chainTo, chainFrom);
                                        //double eta2 = dF.calc_dF(redetected_chainTo.remaining_chain, redetected_chainFrom.remaining_chain);
                                        //double eta = std::max(eta1, eta2);
                                        //Print("eta " << eta);                                                                                                


                                        ParallelLineSimplificationILP p_ilp(graph);
                                        const double epsilon_relaxation = 1.2;
                                        double p_ilpNeededNumberOfEdges
                                                = p_ilp.solve(expandedChainTo, expandedChainFrom,
                                                    epsilon_relaxation * epsilon_error + std::numeric_limits<double>::epsilon() * geo::R,
                                                    eta + std::numeric_limits<double>::epsilon() * geo::R);

                                        size_t preSize = split_edge_chain_pair.chainTo.size() + split_edge_chain_pair.chainFrom.size();

                                        assert(preSize >= p_ilpNeededNumberOfEdges);
                                        double diff = preSize - p_ilpNeededNumberOfEdges;
                                        //Debug("ilpNeededNumberOfEdges:" << ilpNeededNumberOfEdges);
                                        //Debug("diff:" << diff);
                                        double combined_chain_length
                                            = chains::calcChainGeoLength(expandedChainTo, graph)
                                            + chains::calcChainGeoLength(expandedChainFrom, graph);                                    
                                        //errorcounts.weighted_eta += eta * combined_chain_length;
                                        errorcounts.length_sum += combined_chain_length;      
                                        uint node_chain_size = expandedChainTo.size() + expandedChainFrom.size();
                                        errorcounts.node_length_sum += node_chain_size; 
                                        errorcounts.length_sum += combined_chain_length;
                                        errorcounts.addedDiffs += diff * node_chain_size;
                                        errorcounts.weighted_eta += eta * combined_chain_length;
                                    } else {
                                        errorcounts.selfIntersecting++;
                                    }
                                    errorcounts.nofAnalyzedPolygons++;
                                } else {
                                    errorcounts.sameLocation++;
                                }                                        
                                
                            }                                                
                        }
                    }
            }
            errorcounts.print();   
            Print("length: " << errorcounts.length_sum);
            Print("node_length_sum: " << errorcounts.node_length_sum);
            //assert(errorcounts.length_sum != 0);
            Print("-------------------");
            Print("avg_diff_node_size: " << errorcounts.addedDiffs/errorcounts.node_length_sum);
            //Print("avg_diff_length" << errorcounts.addedDiffs/errorcounts.length_sum);            
            Print("avg_eta: " << errorcounts.weighted_eta/errorcounts.length_sum);
            Print("self intersecting ratio: " << (double) errorcounts.selfIntersecting/(double) errorcounts.nofAnalyzedPolygons );                        
            Print("-------------------");
            Print("avg_eta_whole: " << errorcounts.weighted_eta2/  errorcounts.length_sum2);
            
            out_data.dist_Plslash_t = errorcounts.length_sum;
            out_data.deltapslasht = errorcounts.addedDiffs/errorcounts.length_sum;
            out_data.intersection_plasht = (double) errorcounts.selfIntersecting/ (double) errorcounts.nofAnalyzedPolygons;  
            out_data.etapslasht =  errorcounts.weighted_eta/errorcounts.length_sum;
            out_data.etapt = errorcounts.weighted_eta2/  errorcounts.length_sum2;
           
            
        }
        
        void makeDijkstraMeasure() {            
            Print("make_dijkstra_measure");  
            //random Dijkstras
            std::string queryFile = "query_nof_nodes_" + std::to_string(graph.getNrOfNodes()) + ".txt";
        
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
                queryCount = 50000;
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
            //TrackTime tt =VerboseTrackTime();
            using namespace std::chrono;
            steady_clock::time_point t1 =  steady_clock::now();
            Print("nof Dijkstras:" << queryCount);
            Print("Starting Dijkstras");            
            for (uint i = 0; i < queryCount; i++) {
                std::vector<EdgeID> path;
                chdij.calcShopa(srcList.at(i), tgtList.at(i), path);
            }
            //tt.track("Finished Dijkstras");
            
            duration<double> time_span =  duration_cast<duration<double>>(steady_clock::now() - t1);
            double time = time_span.count();
            Print("Dijkstras took " << time << " seconds.\n");
            out_data.dijkstra_time = time;
            Unused(time_span);                        
        }
        
        
        void makeEdgeMeasure() {
            using namespace std::chrono;
            steady_clock::time_point t1 =  steady_clock::now();
            Print("makeEdgeMeasure");
            Print("Accumulating Error and Crossings Counting");
            ErrorCounts error_counts;
            CDTHPCross cdthpC(graph, grid);
            //Grid grid = Grid(1, graph);
            double nofCrossings = 0;
            double accumulatedError = 0;
            //RayCaster raycaster(graph, grid);       
            for (uint edge_id = 0; edge_id < graph.getNrOfEdges(); edge_id++) {
                
                if (graph.isValidEdge(edge_id) && graph.isShortcut(edge_id)) {
                                        
                    Chain chain = graph.getCenterNodes(edge_id);
                    chain.push_front(graph.getEdge(edge_id).src);
                    chain.push_back(graph.getEdge(edge_id).tgt);                    
                                        
                    //getNofCrossings += raycaster.getNofCrossings(chain);
                    nofCrossings += cdthpC.getNofCrossings(chain);
                    double epsilon_error = graph.calcEdgeError(edge_id);
                    accumulatedError += epsilon_error;
                    const CHEdge& edge = graph.getEdge(edge_id);
                    error_counts.length_sum += geo::geoDist(
                                graph.getNode(edge.src),
                                graph.getNode(edge.tgt));
                }
                
            }        
            error_counts.print();        
            Print("length: " << error_counts.length_sum);
            Print("number of Crossings: " << nofCrossings);
            Print("number of Crossings/length: " << nofCrossings/error_counts.length_sum);
            Print("number of Crossings * NrOfValidEdges: " << nofCrossings * graph.getNrOfValidEdges());
            Print("accumulated error/length: " << accumulatedError/error_counts.length_sum);
            
            duration<double> time_span =  duration_cast<duration<double>>(steady_clock::now() - t1);
            Print("makeEdgeMeasure took " << time_span.count() << " seconds.\n");
            Unused(time_span);
        }
        
        
    public:

        CHMeasurer(CHGraph<CHNode, CHEdge> &graph):
            graph(graph), grid(1000, graph), chaindetector(graph), fourDGrid(5, graph) {
            
        }
            
        void makeMeasurement(MeasureOptions m_options) {
            Print("make measure");    
            
            if (m_options.dijkstra) {
                makeDijkstraMeasure();
            }
            /*
            //measurements without path-finding info
            graph.zoom(100, false, false, 0, 0);
            for (uint i = 0; i < graph.getNrOfNodes(); i++) {
                assert(graph.isValidNode(i));
            }*/

            Chains_and_Remainder CaR(graph.getMaxStreetType());
            
            Print("Zooming ");                                               
            //graph.zoom(50, false, 800); //time dist   
            //graph.zoom(0.02, true, true, 800);               
            //graph.zoom(0.02, true, true, 0.005);
            //graph.zoom(0.02, true, false, 0.02, 0); //hat gut funktioniert für ilp
            //graph.zoom(0.02, true, true, 0.0075, 0.02); //hat gut funktioniert für ilp
            //graph.zoom(0.02, true, true, 0.0075, 0.05);
            //graph.zoom(0.02, true, true, 0.0075, 0.0002);
            //graph.zoom(0.02, true, true, 0.03, 0.0002);
            //graph.zoom(0.002, true, true, 2000, 0.00002);
            //graph.zoom(0.02, true, true, 2000, 0.00002); //gut für ilp auch gut 0.002 repspektive 0.00002
            //graph.zoom(0.02, true, false, 0.002, 0.00002); //dist
            //graph.zoom(22.644, true, true, 2000, 1.0);
            //graph.zoom(0.02, true, false, 0.002, 0.0);
            //graph.zoom(0.02, true, true, 0.0075, 0.1); 
            //graph.zoom(45, true, true, 10000.0);  
            
            //km
            //graph.zoom(0.02, false, true, 0, 0.12);
            //graph.zoom(0.02, false, true, 0, 0.00);
            //graph.zoom(0.02, false, true, 0, 0.05);//gut
            //graph.zoom(0.02, false, true, 0, 0.005);//gut
            //graph.zoom(0.2, false, true, 0, 0.002);
            //graph.zoom(0.2, false, true, 0, 0.0);
            graph.zoom(0.2, false, true, 0, 0.01, out_data); //gut für zip
            //graph.zoom(0.2, false, true, 0, 0.1); //sehr gut, aber kleine Daten
            //graph.zoom(0.02, false, true, 0, 0.06);
            
            std::vector<NodeID> visibleNodeIDs;
            for (uint i = 0; i < graph.getNrOfNodes(); i++) {                
                if(graph.isVisibleNode(i)) {
                    visibleNodeIDs.push_back(i);                
                }
            }
            Print("visibleNodeIDs.size()" << visibleNodeIDs.size());
            
            Print("Detecting chains");            
            CaR = chaindetector.detectChains(visibleNodeIDs);            
            
            //CaR = chaindetector.detectChains(activeNodeIDs);
            Print("Number of chains: " << CaR.getNrOfChains());
            
            Print("IdentifyingChainPairs ");                                      
            fourDGrid.identifyPairs(CaR);
            Print("Number of chain pairs: " << CaR.chainPairs.size());
                        
            /*
            if(m_options.edge) {
                //graph.zoom(0.02, true, false, 0.005, 0);
                makeEdgeMeasure();    
            }*/
            if(m_options.chain) {
                make_chain_measure(CaR);
            }                        
            
            if (m_options.ilp) {                
                make_ilp_measure(CaR);
            }
            
            if (m_options.p_ilp) {
                make_p_ilp_measure(CaR);
            }
            
            /*
            //affords only zooming
            if (m_options.expand_diff) {
                graph.zoom(0.02, true, true, 0.2, 0.1);
            }
            */

            out_data.print();
            //beep for test end            
            //cout << '\a' << flush;
            
           
        }

    };
}
