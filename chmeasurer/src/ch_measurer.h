#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "chgraph.h"
#include "geoFunctions.h"
#include "rayCaster.h"
#include "lineSimplificationILP.h"
#include "parallelLineSimplficationILP.h"
#include "chainDetector.h"
#include "chains.h"
#include "4Dgrid.h"
#include "self_intersection_checker.h"
#include "dijkstra.h"
#include "errorCount.h"
#include "discreteFrechet.h"
#include "dijkstra.h"
#include "track_time.h"


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
    }

    class CHMeasurer {
    private:
        CHGraph<CHNode, CHEdge> &graph;
        Grid<CHGraph<CHNode, CHEdge> > grid;
        ChainDetector<CHGraph<CHNode, CHEdge> > chaindetector;
        FourDGrid<CHGraph<CHNode, CHEdge> > fourDGrid;
        SelfIntersectionChecker selfIntersectionChecker;

        std::list<NodeID> getCenterNodes(const EdgeID edge_id) {
            std::list<NodeID> centerNodes;
            const CHEdge &edge = graph.getEdge(edge_id);
            if (graph.isShortcut(edge_id)) {
                centerNodes.splice(centerNodes.end(), getCenterNodes(edge.child_edge1));

                assert(graph.getEdge(edge.child_edge1).tgt == graph.getEdge(edge.child_edge2).src);
                const NodeID centerNode_id = graph.getEdge(edge.child_edge1).tgt;
                centerNodes.push_back(centerNode_id);

                centerNodes.splice(centerNodes.end(), getCenterNodes(edge.child_edge2));
            }
            return centerNodes;
        }
        
        double calcEdgeError(EdgeID edge_id) {
            std::list<NodeID> centerNodes = getCenterNodes(edge_id);
            double maxError = 0;
            const CHEdge &edge = graph.getEdge(edge_id);
            for (NodeID center_node_id : centerNodes) {
                double error = geo::calcPerpendicularLength(graph.getNode(edge.src), graph.getNode(edge.tgt), graph.getNode(center_node_id));
                if (error > maxError) maxError = error;
            }
            return maxError;
        }
        
        //fill up an edge chain with the nodes in between (the center nodes of the connecting edges)
        Chain expandEdgeChain(const EdgeChain &redetected_edge_chain) {
            assert(redetected_edge_chain.chain.size()>=2 &&redetected_edge_chain.edges.size()>=1);
            assert(redetected_edge_chain.chain.size() == redetected_edge_chain.edges.size()+1);
            Chain expandedChain;
            
            auto nodeIt = redetected_edge_chain.chain.begin();
            expandedChain.push_back(*nodeIt);
            nodeIt++;
            for (auto edgeIt = redetected_edge_chain.edges.begin(); edgeIt!= redetected_edge_chain.edges.end(); edgeIt++) {
                std::list<NodeID> centerNodes = getCenterNodes(*edgeIt);
                for (NodeID node_id: centerNodes) {
                    expandedChain.push_back(node_id);
                }                
                expandedChain.push_back(*nodeIt);
                nodeIt++;
            }
            assert(nodeIt == redetected_edge_chain.chain.end());

            return expandedChain;
        }
        
        void make_ilp_measure(Chains_and_Remainder &CaR) {
            Print("make_ilp_measure");  
            lineSimplificationILP ilp(graph);
            ErrorCounts errorcounts;
            /*
            uint chainNotRedetectedCounter = 0;
            uint chainRedetectedCounter = 0;
            uint epsilonerrorEqualszero = 0;
            uint epsilonerrorGreaterzero = 0;
            uint selfIntersecting = 0;
            double addedDiffs = 0;
             * */
                        
            
            //for (ChainsOfType &chainsOfType: CaR.oneWayChainsAccordingToType) {                                                                
            for (uint i = 0; i < CaR.oneWayChainsAccordingToType.size(); i++) {
                ChainsOfType &chainsOfType = CaR.oneWayChainsAccordingToType.at(i);
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                for (Chain &chain: chainsOfType) {                                                                
                    //only big chains are measured
                    if (chain.size() >= 3) {        
                        //get epsilon error by expanding the used edges
                        double epsilon_error = 0;
                        EdgeChain redetected_edge_chain = chaindetector.redetect(chain, i);
                        if (!redetected_edge_chain.chain.empty()) {
                            errorcounts.chainRedetectedCounter++;
                            for (EdgeID edge_id: redetected_edge_chain.edges) {
                                double error = calcEdgeError(edge_id);
                                if (error > epsilon_error) {
                                    epsilon_error = error;
                                }
                            }
                            
                            //expandedChain should be equal to chain
                            Chain expandedChain = expandEdgeChain(redetected_edge_chain);
                            
                            if(chainsAreEqual(expandedChain, chain)) {
                                if(!selfIntersectionChecker.isSelfIntersecting(expandedChain)) {
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
            }

            std::cout << "chainNotRedetectedCounter: " << errorcounts.chainNotRedetectedCounter << std::endl;
            std::cout << "chainRedetectedCounter: " << errorcounts.chainRedetectedCounter << std::endl;
            
            std::cout << "epsilonerrorEqualszero: " << errorcounts.epsilonerrorEqualszero << std::endl;            
            std::cout << "epsilonerrorGreaterzero: " << errorcounts.epsilonerrorGreaterzero << std::endl;
            
            std::cout << "selfIntersecting: " << errorcounts.selfIntersecting << std::endl;
            
            std::cout << "addedDiffs: " << errorcounts.addedDiffs << std::endl;
            
            //SelfIntersectionChecker selfIntersectionChecker;
            
            /*
            lineSimplificationILP ilp(graph);
            if (chain.size()>2 && epsilon_error > 0) {
                        ilp.solve(chain, epsilon_error);
            }
             * */
        }
        
        
        
        void make_p_ilp_measure(Chains_and_Remainder &CaR) {
            Print("make_p_ilp_measure");  
            //ParallelLineSimplificationILP p_ilp(graph);
            ParallelLineSimplificationILP p_ilp(graph);
            DiscreteFrechet dF(graph);
            uint chainNotRedetectedCounter = 0;
            uint chainPairRedetectedCounter = 0;
            uint epsilonerrorEqualszero = 0;
            uint epsilonerrorGreaterzero = 0;
            uint selfIntersecting = 0;
            double addedDiffs = 0;
                        
            
            //for (ChainsOfType &chainsOfType: CaR.oneWayChainsAccordingToType) {                                                                
            for (ChainPair chainpair: CaR.chainPairs) {
                
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                //for (Chain &chain: chainsOfType) {                                                                
                
                    //only big chains are measured
                    if (chainpair.chainTo.size() >= 3 && chainpair.chainFrom.size() >= 3) {     
                       
                        //get epsilon error by expanding the used edges
                        double epsilon_error = 0;
                        EdgeChain redetected_edge_chainTo = chaindetector.redetect(chainpair.chainTo);
                        EdgeChain redetected_edge_chainFrom = chaindetector.redetect(chainpair.chainFrom);
                        if (!redetected_edge_chainTo.chain.empty() && !redetected_edge_chainFrom.chain.empty()) {
                            chainPairRedetectedCounter++;
                            for (EdgeID edge_id: redetected_edge_chainTo.edges) {
                                double error = calcEdgeError(edge_id);
                                if (error > epsilon_error) {
                                    epsilon_error = error;
                                }
                            }
                            for (EdgeID edge_id: redetected_edge_chainFrom.edges) {
                                double error = calcEdgeError(edge_id);
                                if (error > epsilon_error) {
                                    epsilon_error = error;
                                }
                            }
                            
                            
                            //expandedChain should be equal to chain
                            Chain expandedChainTo = expandEdgeChain(redetected_edge_chainTo);                            
                            Chain expandedChainFrom = expandEdgeChain(redetected_edge_chainFrom);/*
                            Chain testchain;
                            testchain.push_back(1);
                            testchain.push_back(2);
                            testchain.push_back(3);
                             * */               
                            if(!selfIntersectionChecker.isSelfIntersecting(expandedChainTo, expandedChainFrom)) {
                                    double eta1 = dF.calc_dF(chainpair.chainTo, chainpair.chainFrom);
                                    double eta2 = dF.calc_dF(redetected_edge_chainTo.chain, redetected_edge_chainFrom.chain);
                                    double eta = std::max(eta1, eta2);
                                    double p_ilpNeededNumberOfEdges
                                        = p_ilp.solve(expandedChainTo, expandedChainFrom,
                                            epsilon_error + std::numeric_limits<double>::epsilon(),
                                            eta + std::numeric_limits<double>::epsilon());
                                    size_t redetectedSize = redetected_edge_chainTo.edges.size() + redetected_edge_chainFrom.edges.size();
                                    
                                    /*
                                    if (!(redetectedSize >= p_ilpNeededNumberOfEdges)) {
                                        
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
                                    addedDiffs += diff;                                    
                                
                            }else {
                                selfIntersecting++;
                            }                           
                        } else {
                            chainNotRedetectedCounter++;
                        }

                        
                    }
                //}
            }

            std::cout << "chainNotRedetectedCounter: " << chainNotRedetectedCounter << std::endl;
            std::cout << "chainRedetectedCounter: " << chainPairRedetectedCounter << std::endl;
            
            std::cout << "epsilonerrorEqualszero: " << epsilonerrorEqualszero << std::endl;            
            std::cout << "epsilonerrorGreaterzero: " << epsilonerrorGreaterzero << std::endl;
            
            std::cout << "selfIntersecting: " << selfIntersecting << std::endl;
            
            std::cout << "addedDiffs: " << addedDiffs << std::endl;
            
            //SelfIntersectionChecker selfIntersectionChecker;
            
            /*
            lineSimplificationILP ilp(graph);
            if (chain.size()>2 && epsilon_error > 0) {
                        ilp.solve(chain, epsilon_error);
            }
             * */
        }
        
        void makeDijkstraMeasure() {            
            
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
            /*
        
            std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<uint> dist(0,g.getNrOfNodes()-1);
	auto rand_node = std::bind (dist, gen);
	std::vector<EdgeID> path;
	for (uint i(0); i<nr_of_dij; i++) {
		NodeID src = rand_node();
		NodeID tgt = rand_node();
        }
        */
        
    public:

        CHMeasurer(CHGraph<CHNode, CHEdge> &graph):
            graph(graph), grid(1000, graph), chaindetector(graph), fourDGrid(10, graph), selfIntersectionChecker(graph) {
            
        }
            
        void makeMeasurement() {
                       
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
            Print("Accumulating Error and Raycasting");
            //Grid grid = Grid(1, graph);
            double getNofCrossings = 0;
            double accumulatedError = 0;
            RayCaster raycaster(graph, grid);
            for (uint edge_id = 0; edge_id < graph.getNrOfEdges(); edge_id++) {
                //const CHEdge &edge = graph.getEdge(edge_id);
                
                if (graph.isValidEdge(edge_id) && graph.isShortcut(edge_id)) {
                    
                    
                    Chain chain = getCenterNodes(edge_id);
                    chain.push_front(graph.getEdge(edge_id).src);
                    chain.push_back(graph.getEdge(edge_id).tgt);                    
                                        
                    getNofCrossings += raycaster.getNofCrossings(chain);
                    double epsilon_error = calcEdgeError(edge_id);
                    accumulatedError += epsilon_error;                        
  
                }
                
            }
            std::cout << "number of Crossings: " << getNofCrossings << std::endl;
            std::cout << "accumulated error: " << accumulatedError << std::endl;
            
            //make_ilp_measure(CaR);  
            //make_p_ilp_measure(CaR);
            
            //makeDijkstraMeasure();
        }

    };
}
