/* always Print() stuff */
#ifdef NVERBOSE
# undef NVERBOSE
#endif

#include "unit_tests.h"

#include "nodes_and_edges.h"
#include "chains.h"
#include "ch_parser.h"
#include "chgraph.h"
#include "ch_measurer.h"
#include "lineSimplificationILP.h"
#include "parallelLineSimplficationILP.h"
#include "dijkstra.h"
//#include "prioritizers.h"

#include <map>
#include <iostream>
#include <random>
#include <chrono>

namespace chm
{

void unit_tests::testAll()
{
    unit_tests::testLineSimplfication();
    unit_tests::testParallelLineSimplfication();
    unit_tests::testDijkstra();
    unit_tests::testCHDijkstra();
    unit_tests::testIsBetween();
}

void unit_tests::testLineSimplfication()
{
	Print("\n=======================");
	Print("TEST: LineSimplfication test.");
	Print("=======================\n");

        GraphInData<CHNode, CHEdge> graphInData;
    
        CH_Parser ch_Parser = CH_Parser();    

        string filepath = "../test_data/ProDPExample_ch_out.graph";
        Test(ch_Parser.parseTxtGraphFile(filepath, graphInData.nodes, graphInData.edges));

        CHGraph<CHNode, CHEdge> graph;
        graph.init(std::move(graphInData));
        
        
        //ProDPxample
        //epsilon = 0.3 => 9 nodes 
        //epsilon = 0.5 => 6 nodes
        
        
        //whole graph is only one chain
        Chain chain;    
        for (uint i = 0; i < 16; i++) {
            chain.push_back(i);
        }
        
        lineSimplificationILP ilp(graph);
        Test(ilp.solve(chain, 0.3) == 9);
        Test(ilp.solve(chain, 0.5) == 6);
                
	Print("\n============================");
	Print("TEST: LineSimplfication test successful.");
	Print("============================\n");
}

void unit_tests::testParallelLineSimplfication()
{
	Print("\n=======================");
	Print("TEST: ParallelLineSimplfication test.");
	Print("=======================\n");

        GraphInData<CHNode, CHEdge> graphInData;    
        CH_Parser ch_Parser = CH_Parser();    

        string filepath = "../test_data/ProP_ILPChainPair_ch_out.graph";
        Test(ch_Parser.parseTxtGraphFile(filepath, graphInData.nodes, graphInData.edges));

        CHGraph<CHNode, CHEdge> graph;
        graph.init(std::move(graphInData));
        
        
        //ProP_ILPChainPair:
        //epsilon = 5 , eta = 2000 => bad result
        //epsilon = ? , eta = ? => good result
                        
        
        //whole graph consists of 2 chains
        Chain chain1;
        for (uint i = 0; i < 4; i++) {
            chain1.push_back(i);
        }
        Chain chain2;
        for (uint i = 4; i < 8; i++) {
            chain2.push_front(i);
        }
        
        //p_ILP_data p_ilp_data = set_p_ILP_data(chain1, chain2, epsilon, eta, nodes);
        
        //parallel
        
        ParallelLineSimplificationILP p_ilp(graph);
        
        Test(p_ilp.solve(chain1, chain2, 4.5 , 2000) == 2);        
        Test(p_ilp.solve(chain1, chain2, 2.9 , 2000) == 4);
        
        Test(p_ilp.solve(chain1, chain2, 2.9 , 2.0001)==4);
        
                
        
	Print("\n============================");
	Print("TEST: ParallelLineSimplfication test successful.");
	Print("============================\n");
}

void unit_tests::testCHDijkstra()
{
	Print("\n============================");
	Print("TEST: Start CHDijkstra test.");
	Print("============================\n");

	//typedef CHEdge Shortcut;
	//typedef CHGraph<CHNode, CHEdge> CHGraphOSM;

        GraphInData<CHNode, CHEdge> graphInData;    
        CH_Parser ch_Parser = CH_Parser();    

        string filepath = "../test_data/15kSZHK_ch_out.graph";
        Test(ch_Parser.parseTxtGraphFile(filepath, graphInData.nodes, graphInData.edges));

        CHGraph<CHNode, CHEdge> g;
        g.init(std::move(graphInData));
        g.zoom(100, false, 0);
	/*
	Graph<CHNode, CHEdge> g;
	g.init(FormatSTD::Reader::readGraph<CHNode, CHEdge>("../test_data/15kSZHK.txt"));

	
	CHGraphOSM chg;
	chg.init(FormatSTD::Reader::readGraph<CHNode, Shortcut>("../test_data/15kSZHK.txt"));

	
	CHConstructor<CHNode, CHEdge> chc(chg, 2);
	std::vector<NodeID> all_nodes(g.getNrOfNodes());
	for (NodeID i(0); i<all_nodes.size(); i++) {
		all_nodes[i] = i;
	}
	chc.quickContract(all_nodes, 4, 5);
	chc.contract(all_nodes);
	chc.rebuildCompleteGraph();
         * */

	/* Random Dijkstras */
	Print("\nStarting random Dijkstras.");
	uint nr_of_dij(10);
	Dijkstra<CHNode, CHEdge> dij(g);
	CHDijkstra<CHNode, CHEdge> chdij(g);

	std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<uint> dist(0,g.getNrOfNodes()-1);
	auto rand_node = std::bind (dist, gen);
	std::vector<EdgeID> path;
	for (uint i(0); i<nr_of_dij; i++) {
		NodeID src = rand_node();
		NodeID tgt = rand_node();
		Debug("From " << src << " to " << tgt << ".");
		Test(dij.calcShopa(src,tgt,path) == chdij.calcShopa(src,tgt,path));
	}

	// Export (destroys graph data)
	//writeCHGraphFile<FormatSTD::Writer>("../out/ch_15kSZHK.txt", chg.exportData());

	Print("\n=================================");
	Print("TEST: CHDijkstra test successful.");
	Print("=================================\n");
}


void unit_tests::testDijkstra()
{
	Print("\n============================");
	Print("TEST: Start Dijkstra test.");
	Print("============================\n");

        GraphInData<CHNode, CHEdge> graphInData;    
        CH_Parser ch_Parser = CH_Parser();    

        string filepath = "../test_data/DijkstraTest_ch_out.graph";
        Test(ch_Parser.parseTxtGraphFile(filepath, graphInData.nodes, graphInData.edges));

        CHGraph<CHNode, CHEdge> g;
        g.init(std::move(graphInData));
        g.zoom(100, false, 0);
        
	

	Dijkstra<CHNode, CHEdge> dij(g);
	std::vector<EdgeID> path;
	NodeID tgt(g.getNrOfNodes() - 1);
	uint dist = dij.calcShopa(0, tgt, path);

	Print("Dist of Dijkstra from 0 to " << tgt << ": " << dist);
	Test(dist == 9);

	Print("Shortest path from 0 to " << tgt << ":");
	for (auto edge_id: path) {
		CHEdge const& edge(g.getEdge(edge_id));
		Print("EdgeID: " << edge_id << ", src: " << edge.src << ", tgt: " << edge.tgt);
	}

	Print("Test if shortest paths are the same from both sides for the 'test' graph.");
	for (NodeID src(0); src< (int) g.getNrOfNodes(); src++) {
		for (NodeID tgt(src); tgt< (int) g.getNrOfNodes(); tgt++) {
			Test(dij.calcShopa(src, tgt, path) == dij.calcShopa(tgt, src, path));
		}
	}

	Print("\n=================================");
	Print("TEST: Dijkstra test successful.");
	Print("=================================\n");
}

void unit_tests::testIsBetween()
{
	Print("\n============================");
	Print("TEST: Start IsBetween test.");
	Print("============================\n");
        
        
        CHGraph<CHNode, CHEdge> graph;
        Grid<CHGraph<CHNode, CHEdge> > grid(1, graph);
        RayCaster raycaster(graph, grid);
        CHNode chlinenodeSrc;
        chlinenodeSrc.lon = 3;
        chlinenodeSrc.lat = 2;
        CHNode chlinenodeTgt;
        chlinenodeTgt.lon = 6;
        chlinenodeTgt.lat = 7;
        CHLine chline(chlinenodeSrc, chlinenodeTgt);
                
        CHNode in1;
        in1.lon = 1;
        in1.lat = 9;
        Test(raycaster.isBetween(chline, in1));
        
        
        CHNode in2;
        in2.lon = 8;
        in2.lat = 0;
        Test(raycaster.isBetween(chline, in2));
        
        CHNode out1;
        out1.lon = 4;
        out1.lat = 1;
        Test(!raycaster.isBetween(chline, out1));
        
        CHNode out2;
        out2.lon = 10;
        out2.lat = 5;
        Test(!raycaster.isBetween(chline, out2));                                

	Print("\n=================================");
	Print("TEST: IsBetween test successful.");
	Print("=================================\n");
}

}
