/* always Print() stuff */
#ifdef NVERBOSE
# undef NVERBOSE
#endif

#include "unit_tests.h"

#include "nodes_and_edges.h"
#include "graph.h"
#include "file_formats.h"
#include "chgraph.h"
#include "ch_constructor.h"
#include "dijkstra.h"
#include "prioritizers.h"
#include "simplification/cdt_matching.h"
#include "self_intersection_checker.h"


#include <map>
#include <iostream>
#include <random>
#include <chrono>

namespace chc
{

void unit_tests::testAll()
{
	unit_tests::testNodesAndEdges();
	//unit_tests::testGraph();
	unit_tests::testCHConstructor();
	unit_tests::testCHDijkstra();
	//unit_tests::testDijkstra();
	//unit_tests::testPrioritizers();
        unit_tests::testCDTMatching();
}

void unit_tests::testNodesAndEdges()
{
	Print("\n=================================");
	Print("TEST: Start Nodes And Edges test.");
	Print("=================================\n");

	typedef MetricEdge<Edge> MetricEdge;
	typedef CHEdge<MetricEdge> CHEdge;

	CHNode<Node> node0(Node(0), 0);
	CHNode<Node> node1(Node(1), 1);
	CHNode<Node> node2(Node(2), 2);

	MetricEdge edge0(Edge(0, node0.id, node1.id, 42), 23);
	MetricEdge edge1(Edge(1, node1.id, node2.id, 24), 32);
	CHEdge ch_edge(make_shortcut(edge0, edge1));

	UnitTest(otherNode(edge0, EdgeType::IN) == 0);
	UnitTest(otherNode(ch_edge, EdgeType::OUT) == 2);

	Print("\n======================================");
	Print("TEST: Nodes and edges test successful.");
	Print("======================================\n");
}
/*
void unit_tests::testGraph()
{
	Print("\n=======================");
	Print("TEST: Start Graph test.");
	Print("=======================\n");

	Graph<OSMNode, Edge> g;
	g.init(FormatSTD::Reader::readGraph<OSMNode, Edge>("../test_data/15kSZHK.txt"));

	// Test the normal iterator.
	for (NodeID node_id(0); node_id<g.getNrOfNodes(); node_id++) {
		UnitTest(g.getNode(node_id).id == node_id);

		// Find for every out_edge the corresponding in edge.
		for (auto const& out_edge: g.nodeEdges(node_id, EdgeType::OUT)) {
			bool found(false);

			for (auto const& in_edge: g.nodeEdges(out_edge.tgt, EdgeType::IN)) {
				if (equalEndpoints(in_edge, out_edge)) {
					found = true;
					break;
				}
			}

			UnitTest(found);
		}
	}

	Print("\n============================");
	Print("TEST: Graph test successful.");
	Print("============================\n");
}*/

void unit_tests::testCHConstructor()
{
	Print("\n===============================");
	Print("TEST: Start CHConstructor test.");
	Print("===============================\n");

	typedef CHEdge<OSMEdge> Shortcut;
	typedef CHGraph<OSMNode, OSMEdge> CHGraphOSM;

	CHGraphOSM g;
	g.init(FormatSTD::Reader::readGraph<OSMNode, Shortcut>("../test_data/test"));

	CHConstructor<OSMNode, OSMEdge> chc(g, 2);

	/*
	 * Test the independent set construction.
	 */

	std::vector<NodeID> all_nodes(g.getNrOfNodes());
	for (NodeID i(0); i<all_nodes.size(); i++) {
		all_nodes[i] = i;
	}

	std::vector<NodeID> remaining_nodes(all_nodes);
	std::vector<bool> is_in_ind_set(g.getNrOfNodes(), false);

	auto independent_set = chc.calcIndependentSet(remaining_nodes);
	Print("Size of the independent set of all nodes: " << independent_set.size());

	for (NodeID node: independent_set) {
		is_in_ind_set[node] = true;
	}

	for (NodeID node: independent_set) {
		for (auto const& edge: g.nodeEdges(node, EdgeType::OUT)) {
			UnitTest(!is_in_ind_set[edge.tgt]);
		}
	}

	/*
	 * Test the contraction.
	 */
	chc.contract(all_nodes);

	// Export
	writeCHGraphFile<FormatSTD::Writer>("../out/ch_test", g.exportData());

	Print("\n====================================");
	Print("TEST: CHConstructor test successful.");
	Print("====================================\n");
}

void unit_tests::testCHDijkstra()
{
	Print("\n============================");
	Print("TEST: Start CHDijkstra test.");
	Print("============================\n");

	typedef CHEdge<OSMEdge> Shortcut;
	typedef CHGraph<OSMNode, OSMEdge> CHGraphOSM;

	/* Init normal graph */
	Graph<OSMNode, OSMEdge> g;
	g.init(FormatSTD::Reader::readGraph<OSMNode, OSMEdge>("../test_data/15kSZHK.txt"));

	/* Init CH graph */
	CHGraphOSM chg;
	chg.init(FormatSTD::Reader::readGraph<OSMNode, Shortcut>("../test_data/15kSZHK.txt"));

	/* Build CH */
	CHConstructor<OSMNode, OSMEdge> chc(chg, 2);
	std::vector<NodeID> all_nodes(g.getNrOfNodes());
	for (NodeID i(0); i<all_nodes.size(); i++) {
		all_nodes[i] = i;
	}
	chc.quickContract(all_nodes, 4, 5);
	chc.contract(all_nodes);
	chc.rebuildCompleteGraph();

	/* Random Dijkstras */
	Print("\nStarting random Dijkstras.");
	uint nr_of_dij(10);
	Dijkstra<OSMNode, OSMEdge> dij(g);
	CHDijkstra<OSMNode, OSMEdge> chdij(chg);

	std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<uint> dist(0,g.getNrOfNodes()-1);
	auto rand_node = std::bind (dist, gen);
	std::vector<EdgeID> path;
	for (uint i(0); i<nr_of_dij; i++) {
		NodeID src = rand_node();
		NodeID tgt = rand_node();
		Debug("From " << src << " to " << tgt << ".");
		UnitTest(dij.calcShopa(src,tgt,path) == chdij.calcShopa(src,tgt,path));
	}

	// Export (destroys graph data)
	writeCHGraphFile<FormatSTD::Writer>("../out/ch_15kSZHK.txt", chg.exportData());

	Print("\n=================================");
	Print("TEST: CHDijkstra test successful.");
	Print("=================================\n");
}
/*
void unit_tests::testDijkstra()
{
	Print("\n============================");
	Print("TEST: Start Dijkstra test.");
	Print("============================\n");

	Graph<OSMNode, Edge> g;
	g.init(FormatSTD::Reader::readGraph<OSMNode, Edge>("../test_data/test"));

	Dijkstra<OSMNode, Edge> dij(g);
	std::vector<EdgeID> path;
	NodeID tgt(g.getNrOfNodes() - 1);
	uint dist = dij.calcShopa(0, tgt, path);

	Print("Dist of Dijkstra from 0 to " << tgt << ": " << dist);
	UnitTest(dist == 18);

	Print("Shortest path from 0 to " << tgt << ":");
	for (auto edge_id: path) {
		Edge const& edge(g.getEdge(edge_id));
		Print("EdgeID: " << edge_id << ", src: " << edge.src << ", tgt: " << edge.tgt);
	}

	Print("Test if shortest paths are the same from both sides for the 'test' graph.");
	for (NodeID src(0); src<g.getNrOfNodes(); src++) {
		for (NodeID tgt(src); tgt<g.getNrOfNodes(); tgt++) {
			UnitTest(dij.calcShopa(src, tgt, path) == dij.calcShopa(tgt, src, path));
		}
	}

	Print("\n=================================");
	Print("TEST: Dijkstra test successful.");
	Print("=================================\n");
}*/

void unit_tests::testPrioritizers()
{
	Print("\n=============================");
	Print("TEST: Start Prioritizer test.");
	Print("=============================\n");

	typedef CHEdge<OSMEdge> Shortcut;
	typedef CHGraph<OSMNode, OSMEdge> CHGraphOSM;

	/* Init normal graph */
	Graph<OSMNode, OSMEdge> g;
	g.init(FormatSTD::Reader::readGraph<OSMNode, OSMEdge>("../test_data/test"));

	size_t const last = from_enum(LastPrioritizerType);
	for (size_t t = 0; t <= last; ++t) {
		PrioritizerType type(static_cast<PrioritizerType>(t));
		Print("\n------------------------------------");
		Print("Testing prioritizer type: " << to_string(type));
		Print("------------------------------------\n");

		/* Init CH graph */
		CHGraphOSM chg;
		chg.init(FormatSTD::Reader::readGraph<OSMNode, Shortcut>("../test_data/test"));

		/* Build CH */
		CHConstructor<OSMNode, OSMEdge> chc(chg, 2);
		std::vector<NodeID> all_nodes(g.getNrOfNodes());
		for (NodeID i(0); i<all_nodes.size(); i++) {
			all_nodes[i] = i;
		}
		std::random_shuffle(all_nodes.begin(), all_nodes.end()); /* random initial node order */
		auto prioritizer(createPrioritizer(type, chg, chc));
		if (prioritizer == nullptr && type == PrioritizerType::NONE) { continue; }
		chc.contract(all_nodes, *prioritizer, false);
		chc.rebuildCompleteGraph();

		/* Random Dijkstras */
		Print("\nStarting random Dijkstras.");
		uint nr_of_dij(1000);
		Dijkstra<OSMNode, OSMEdge> dij(g);
		CHDijkstra<OSMNode, OSMEdge> chdij(chg);

		std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_int_distribution<uint> dist(0,g.getNrOfNodes()-1);
		auto rand_node = std::bind (dist, gen);
		std::vector<EdgeID> path;
		for (uint i(0); i<nr_of_dij; i++) {
			NodeID src = rand_node();
			NodeID tgt = rand_node();
			Debug("From " << src << " to " << tgt << ".");
			UnitTest(dij.calcShopa(src,tgt,path) == chdij.calcShopa(src,tgt,path));
		}
	}

	Print("\n==================================");
	Print("TEST: Prioritizer test successful.");
	Print("==================================\n");
}

void unit_tests::testCDTMatching()
{
	Print("\n============================");
	Print("TEST: Start CDTMatching test.");
	Print("============================\n");

	typedef CHEdge<OSMEdge> Shortcut;
	typedef CHGraph<OSMNode, OSMEdge> CHGraphOSM;
	
	/* Init CH graph */
	CHGraphOSM chg;
	chg.init(FormatSTD::Reader::readGraph<OSMNode, Shortcut>("../test_data/cdt_matchingTest2.txt"));

        using namespace ls;
        
        std::vector<NodeInBox> nodesInBox;
        std::list<Intervall2>::iterator intervallIt;
        PrioNodeHeap heap;
             
        Chain chain1;
        std::list<PrioNodeHandle> list1;        
        for (int i = 0; i < 9; i++) {            
            chain1.push_back(i);
            PrioNode2 p(i, intervallIt, nodesInBox, nodesInBox);
            PrioNodeHandle h = heap.push(p);
            list1.push_back(h);
        }        
                
        Chain chain2;
        std::list<PrioNodeHandle> list2;        
        for (int i = 9; i < 14; i++) {
            chain2.push_back(i);
            PrioNode2 p(i, intervallIt, nodesInBox, nodesInBox);
            PrioNodeHandle h = heap.push(p);
            list2.push_back(h);
        }   
        
        //constrains space between chains
        chain1.push_front(chain2.front());
        chain2.push_back(chain1.back());   
        
        Print("self Intersection_test");
        SelfIntersectionChecker<CHGraphOSM> sic(chg);
        UnitTest(!sic.isSelfIntersecting(chain1, chain2));
        
        Print("CDTMatching");
        ls::mc::CDTMatching2<CHGraphOSM> cdtm(chg);        
        cdtm.match(list1, list2);                  
        Print("CDTmatched");
        
        for (auto it = list1.begin(); it!= list1.end(); it++) {
            const PrioNode2 pn = *(*it);
                        
            switch (pn.node_id) {
                case 0: {                    
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 9);
                    break;
                }
                case 1: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 9);
                    break;
                }
                case 3: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 10);
                    break;
                }                
                case 5: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 11);
                    break;
                }
                case 6: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 11);
                    break;
                }
                case 7: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 12);
                    break;
                }
                case 8: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 13);
                    break;
                }
                default: {             
                    break;
                }                    
            }            
        }
        
        for (auto it = list2.begin(); it!= list2.end(); it++) {
            const PrioNode2 pn = *(*it);   
            
            switch (pn.node_id) {                                
                case 12: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 7);
                    break;
                }
                case 13: {
                    const PrioNode2 follower = *pn.follower_h;
                    UnitTest(follower.node_id == 8);
                    break;
                }      
                default: {             
                    break;
                }  
            }            
        }
        
	Print("\n=================================");
	Print("TEST: CDTMatching test successful.");
	Print("=================================\n");
}

}
