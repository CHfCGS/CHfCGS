/* always Print() stuff */
#ifdef NVERBOSE
# undef NVERBOSE
#endif

#include "unit_tests.h"

#include "chgraph.h"
#include "defs.h"

#include "nodes_and_edges.h"
#include "chains.h"
#include "ch_parser.h"
#include "grid.h"


#include "ch_measurer.h"
#include "ILP/lineSimplificationILP.h"
#include "ILP/parallelLineSimplficationILP.h"
#include "dijkstra.h"
#include "discreteFrechet.h"

//#include "CGAL/range_tree.h"
//#include "CGAL/cdthp_cross.h"

#include <map>
#include <iostream>
#include <random>
#include <chrono>

namespace chm
{

    namespace unit_tests
    {
	void testLineSimplfication();
        void testParallelLineSimplfication();    
        void testCDTHCross();
    }
    
void unit_tests::testAll()
{
    
    unit_tests::testLineSimplfication();
    unit_tests::testParallelLineSimplfication();
    
    unit_tests::testDijkstra();
    unit_tests::testCHDijkstra();
    unit_tests::testIsBetween();
     
    unit_tests::testdiscreteFrechet();
    unit_tests::testCDTHCross();
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
	

	// Random Dijkstras
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
        CHNode chlinenodeSrc(3, 2);        
        CHNode chlinenodeTgt(6, 7);        
        CHLine chline(chlinenodeSrc, chlinenodeTgt);
                
        CHNode in1(1, 9);        
        Test(raycaster.isBetween(chline, in1));        
        
        CHNode in2(8, 0);        
        Test(raycaster.isBetween(chline, in2));
        
        CHNode out1(4, 1);        
        Test(!raycaster.isBetween(chline, out1));
        
        CHNode out2(10, 5);        
        Test(!raycaster.isBetween(chline, out2));                                

	Print("\n=================================");
	Print("TEST: IsBetween test successful.");
	Print("=================================\n");
}

void unit_tests::testdiscreteFrechet()
{
	Print("\n============================");
	Print("TEST: Start discreteFrechet test.");
	Print("============================\n");
        
        
        CHGraph<CHNode, CHEdge> graph;
        GraphInData<CHNode, CHEdge> graphInData;
        
        
        graphInData.nodes.push_back(CHNode(0, 0));
        graphInData.nodes.push_back(CHNode(1, 0));
        graphInData.nodes.push_back(CHNode(2, 0));
        graphInData.nodes.push_back(CHNode(3, 0));
        graphInData.nodes.push_back(CHNode(4, 0));
        graphInData.nodes.push_back(CHNode(5, 0));
        graphInData.nodes.push_back(CHNode(8, 0));
        
        Chain chainTo;
        for (int i = 0; i < 7; i++) {
            chainTo.push_back(i);
        }

        graphInData.nodes.push_back(CHNode(0, 2));
        graphInData.nodes.push_back(CHNode(2, 1));
        graphInData.nodes.push_back(CHNode(5, 3));
        graphInData.nodes.push_back(CHNode(8, 2)); 
        Chain chainFrom;
        for (int i = 7; i < 11; i++) {
            chainFrom.push_front(i);
        }
        graph.init(std::move(graphInData));
        

                
                
        DiscreteFrechet df(graph);                                       
        Test(df.calc_dF(chainTo, chainFrom)==3);      
                
        

	Print("\n=================================");
	Print("TEST: discreteFrechet test successful.");
	Print("=================================\n");
}

void unit_tests::testCDTHCross()
{
	Print("\n============================");
	Print("TEST: Start CDTHCross test.");
	Print("============================\n");
        
        
        CHGraph<CHNode, CHEdge> graph;
        GraphInData<CHNode, CHEdge> graphInData;
        
        
        graphInData.nodes.push_back(CHNode(0, 0));
        graphInData.nodes.push_back(CHNode(1, 1));
        graphInData.nodes.push_back(CHNode(2, 0));
        graphInData.nodes.push_back(CHNode(0, -2));
        graphInData.nodes.push_back(CHNode(-2, 0));
        graphInData.nodes.push_back(CHNode(1, 3));
        graphInData.nodes.push_back(CHNode(6, -2));
        graphInData.nodes.push_back(CHNode(8, 0));
        
        Chain chain;
        for (int i = 0; i < 8; i++) {
            chain.push_back(i);
        }

        //outliers
        graphInData.nodes.push_back(CHNode(1, 0.5));
        graphInData.nodes.push_back(CHNode(0, -1));
        graphInData.nodes.push_back(CHNode(-3, -1));
        graphInData.nodes.push_back(CHNode(1.5, 1.5)); 
        graphInData.nodes.push_back(CHNode(3, -1)); 
        graphInData.nodes.push_back(CHNode(4, 0.5)); 
        graphInData.nodes.push_back(CHNode(7, -0.5)); 
        
        graph.init(std::move(graphInData));
        Grid<CHGraph<CHNode, CHEdge> > grid(1, graph);
        
        
        //RangeTree rangeTree(graph);
        CDTHPCross cdthpC(graph, grid);
        for (uint i= 0; i< 100; i++) {
            Print("cdthpC.getNofCrossings(chain)" << cdthpC.getNofCrossings(chain));
        }
                
        Test(cdthpC.getNofCrossings(chain) == 4);
                   
	Print("\n=================================");
	Print("TEST: CDTHCross test successful.");
	Print("=================================\n");
}

}
