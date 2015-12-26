#pragma once

#include <assert.h>

#include "widthAndColors.h"
#include "nodesAndEdges.h"

#define DEBUG_BUILD

#ifdef DEBUG_BUILD
#   define DEBUG(x) do { std::cout << x << std::endl; } while (0)
#else
#   define DEBUG(x)
#endif

using namespace std;

//class Zoomer {	
namespace Zoomer {

	struct ValidLevelInfo {
	    Level allValidUntilLevel;
	    size_t validNofNodesOnThatLevel;
	};			
	
	
	static void removeShortcut(const vector<CHNode> &nodes, vector<CHEdge> &edges, const EdgeID edgeID){
		CHEdge &edge = edges[edgeID];
		if (edge.remaining && edge.is_shortcut()) {        
			if (edges[edge.child_edge1].is_valid(nodes)) {
				edge.remaining = false;
				removeShortcut(nodes, edges, edge.child_edge1);
				assert(edges[edge.child_edge2].is_valid(nodes));
				removeShortcut(nodes, edges, edge.child_edge2);
			}
			/*
			if (edges[edge.child_edge2].is_valid(nodes)) {
				edge.remaining = false;
				removeShortcut(nodes, edges, edge.child_edge2);
			}*/
		}   
	}
	
	static void expandEdge(const vector<CHNode> &nodes, vector<CHEdge> &edges, const EdgeID edgeID, double expandSize) {
		CHEdge &edge = edges[edgeID];
		if (edge.is_shortcut() && edge.remaining) {
			//if (edge.getDist(nodes) > expandSize) {
			if (edge.dist > expandSize) {
			//double BendingRatio = calcBendingRatio(nodes [edge.src], nodes[edge.getCenterPoint(edges)], nodes[edge.tgt]);
			//if (BendingRatio > expandSize) {
				edge.remaining = false;
				edges[edge.child_edge1].remaining = true;
				edges[edge.child_edge2].remaining = true;
				expandEdge(nodes, edges, edge.child_edge1, expandSize);
				expandEdge(nodes, edges, edge.child_edge2, expandSize);
			}
		}
	}
	
	static Level calculateZoomlvl(vector<CHNode> &nodes, double percent_of_showed_nodes) {		
		Level maxLevel = 0;
		for (CHNode node: nodes) {
			if (node.lvl> maxLevel) {
				maxLevel = node.lvl;
			}
		}
			
		vector<uint> nofNodesPerLevel(maxLevel+1, 0); //one field for level 0 is also needed
		for (CHNode node: nodes) {
			nofNodesPerLevel[node.lvl]++;
		}
				
		assert(percent_of_showed_nodes >= 0 && percent_of_showed_nodes <= 100);
		double neededNofNodes = (nodes.size() * (percent_of_showed_nodes/100.0));
		uint nofCollectedNodes = 0;
		int i = maxLevel;
		while (nofCollectedNodes < neededNofNodes) {
			nofCollectedNodes += nofNodesPerLevel[i];
			i--;
		}    
		return i+1;
	}	 	
	
	static ValidLevelInfo calcValidLevel(vector<CHNode> &nodes, size_t numberOfValidNodes) {
		ValidLevelInfo vli;
		Level maxLevel = 0;
		for (CHNode node : nodes) {
		    if (node.lvl > maxLevel) {
		        maxLevel = node.lvl;
		    }
		}

		std::vector<uint> nofNodesPerLevel(maxLevel + 1, 0); //one field for level 0 is also needed
		for (CHNode node : nodes) {
		    nofNodesPerLevel[node.lvl]++;
		}

		assert(numberOfValidNodes >= 0 && numberOfValidNodes <= nodes.size());        
		uint nofCollectedNodes = 0;
		uint oldNofCollectedNodes = 0;
		Level l = maxLevel;
		while (nofCollectedNodes < numberOfValidNodes) {
		    oldNofCollectedNodes = nofCollectedNodes;
		    nofCollectedNodes += nofNodesPerLevel[l];
		    l--;
		}
		vli.allValidUntilLevel = l + 1;
		vli.validNofNodesOnThatLevel = numberOfValidNodes - oldNofCollectedNodes;
		return vli;
	}

	static void markValidNodes(vector<CHNode> &nodes, double percent_of_valid_nodes) { 
	    assert(percent_of_valid_nodes >= 0 && percent_of_valid_nodes <= 100);
	    size_t numberOfValidNodes = (uint) trunc((percent_of_valid_nodes/100.0) * nodes.size());
	    ValidLevelInfo vli = calcValidLevel(nodes, numberOfValidNodes);
//	    _validNodes.resize(_nodes.size());
	    size_t nofValidatedNodesOnCriticalLevel = 0;
	    for (NodeID nodeID = 0; nodeID < (int) nodes.size(); nodeID++) {
		CHNode &node = nodes[nodeID];
		if (node.lvl > vli.allValidUntilLevel) {
		    nodes[nodeID].remaining = true;
		}else if(node.lvl == vli.allValidUntilLevel) {
		    if(nofValidatedNodesOnCriticalLevel < vli.validNofNodesOnThatLevel) {
			nofValidatedNodesOnCriticalLevel++;
			nodes[nodeID].remaining = true;
		    }else {
			nodes[nodeID].remaining = false;
		    }        
		}else {
		    nodes[nodeID].remaining = false;
		}
	    }
	    assert(nofValidatedNodesOnCriticalLevel == vli.validNofNodesOnThatLevel);    
	}

	static void zoom(std::vector<CHNode> &ch_nodes, std::vector<CHEdge> &ch_edges, std::vector<Node> &nodes, std::vector<Edge> &edges,
		double percent_of_showed_nodes = 2, //normal: 2
		double expandSize = 0.002, //normal: 0.002
		bool expand = false,
                bool spareShortcuts = false
		) {	

		DEBUG("markValidNodes");
		markValidNodes(ch_nodes, percent_of_showed_nodes);					
			
		/*
		DEBUG("Calculating zoomlvl");
		Level zoomlvl = calculateZoomlvl(ch_nodes, percent_of_showed_nodes);    
		
		DEBUG("Processing zoomlvl");
		//process zoomlvl    
		//build graph of level zoomlvl
		//mark valid nodes
		for (uint nodeID = 0; nodeID < ch_nodes.size(); nodeID++){
			ch_nodes[nodeID].remaining = (ch_nodes[nodeID].lvl >= zoomlvl);
		}        
		*/

		DEBUG("Processing edges");
		//mark valid edges
		for (uint edgeID = 0; edgeID < ch_edges.size(); edgeID++){
			ch_edges[edgeID].remaining = ch_edges[edgeID].is_valid(ch_nodes);
		}    
		
		//remove high-level shortcuts
		//not all valid shortcuts will remain
		for (uint edgeID = 0; edgeID < ch_edges.size(); edgeID++){        
			removeShortcut(ch_nodes, ch_edges, edgeID);             
		}
                if (spareShortcuts) {
                    //do not draw visually unpleasing shortcuts
                    for (uint edgeID = 0; edgeID < ch_edges.size(); edgeID++){
                        if (ch_edges[edgeID].speed == -1000) {
                            ch_edges[edgeID].remaining = false;       
                        }				
                    }
                }
                
		
		if (expand) {
			DEBUG("Processing expandSize");
			//process expandSize
			for (uint edgeID = 0; edgeID < ch_edges.size(); edgeID++){
				//if (ch_edges[edgeID].remaining) {
				expandEdge(ch_nodes, ch_edges, edgeID, expandSize);
				//}        
			}  
		}
		
		/*			
		EdgeID nofRemainingEdges = 0;
		for (auto it = ch_edges.begin(); it != ch_edges.end(); it++){
			if (it->remaining) {
				nofRemainingEdges++;
			}       
		}
		*/
		
		//TODO omit this part, nodes don't need to be refreshed, maybe no 2 node structures at all
		nodes.clear();
		for (auto it = ch_nodes.begin(); it != ch_nodes.end(); it++){
			Node node = Node(it->lat, it->lon);
			nodes.push_back(node);
		}
		
		edges.clear();
		for (auto it = ch_edges.begin(); it != ch_edges.end(); it++){
			if (it->remaining) {					
				Edge edge = Edge(it->src, it->tgt, 1 /*width(it->type)*/, color(it->type));
				edges.push_back(edge);
			}       
		}												
	}		
}
