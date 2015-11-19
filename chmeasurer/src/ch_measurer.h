#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "chgraph.h"
#include "geoFunctions.h"
#include "rayCaster.h"
#include "lineSimplificationILP.h"
#include "chainDetector.h"
#include "chains.h"
#include "4Dgrid.h"


#include <chrono>
#include <queue>
#include <mutex>
#include <vector>
#include <omp.h>
#include <algorithm>

namespace chm {

    class CHMeasurer {
    private:
        CHGraph<CHNode, CHEdge> &graph;
        Grid<CHGraph<CHNode, CHEdge>> grid;
        ChainDetector<CHGraph<CHNode, CHEdge> > chaindetector;
        FourDGrid<CHGraph<CHNode, CHEdge> > fourDGrid;

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
        
    public:

        CHMeasurer(CHGraph<CHNode, CHEdge> &graph) : graph(graph), grid(1000, graph), chaindetector(graph), fourDGrid(10, graph) {

        }

        void makeMeasurement() {            
            //measurements without path-finding info
            graph.zoom(0, false, 0);
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
            
            
            //Grid grid = Grid(1, graph);
            double getNofCrossings = 0;
            double accumulatedError = 0;
            RayCaster raycaster(graph, grid);
            
            
            
            for (uint edge_id = 0; edge_id < graph.getNrOfEdges(); edge_id++) {
                //const CHEdge &edge = graph.getEdge(edge_id);
                
                if (graph.isValidEdge(edge_id) && graph.isShortcut(edge_id)) {
                    
                    /*
                    Chain chain = getCenterNodes(edge_id);
                    chain.push_front(edge.src);
                    chain.push_back(edge.tgt);
                     * */
                                        
                    //getNofCrossings += raycaster.getNofCrossings(chain);
                    double epsilon_error = calcEdgeError(edge_id);
                    accumulatedError += epsilon_error;                        
  
                }
                
            }
            std::cout << "number of Intersections: " << getNofCrossings << std::endl;
            std::cout << "accumulated error: " << accumulatedError << std::endl;
            
            
            //for (ChainsOfType &chainsOfType: CaR.oneWayChainsAccordingToType) {                                                                
            for (uint i = 0; i < CaR.oneWayChainsAccordingToType.size(); i++) {
                ChainsOfType &chainsOfType = CaR.oneWayChainsAccordingToType.at(i);
                //for (auto it = chainsOfType.begin(); it!=chainsOfType.end(); ++it) {                       
                for (Chain &chain: chainsOfType) {                                                                
                    //only big chains are measured
                    if (chain.size() >= 3) {                                                                                      
                        EdgeChain edge_chain = chaindetector.redetect(chain, i);
                    }
                }
            }

            
            /*
            lineSimplificationILP ilp(graph);
            if (chain.size()>2 && epsilon_error > 0) {
                        ilp.solve(chain, epsilon_error);
            }
             * */
            
        }

    };
}
