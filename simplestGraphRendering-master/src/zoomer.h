#pragma once

#include <assert.h>
#include <list>

#include "widthAndColors.h"
#include "nodesAndEdges.h"
#include "geoFunctions.h"

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

	/*
	static double CHGraph<NodeT, EdgeT>::edgeError(const vector<CHNode> &nodes, const vector<CHEdge> &edges const EdgeID edge_id) const {
		const EdgeT& edge = getEdge(edge_id);
		const NodeT& src_node = getNode(edge.src);
		const NodeT& tgt_node = getNode(edge.tgt);
		const EdgeT& child_edge1 = getEdge(edge.child_edge1);
		const EdgeT& child_edge2 = getEdge(edge.child_edge2);
		assert(child_edge1.tgt == child_edge2.src);
		const NodeT center_node = getNode(child_edge1.tgt);
		return geo::calcPerpendicularLength(src_node, tgt_node, center_node);
	}
	*/
 
	static std::list<NodeID> getCenterNodes(const vector<CHNode> &nodes, const vector<CHEdge> &edges, const EdgeID edge_id) {
            std::list<NodeID> centerNodes;
            const CHEdge &edge = edges[edge_id];
            if (edge.is_shortcut()) {
                centerNodes.splice(centerNodes.end(), getCenterNodes(nodes, edges, edge.child_edge1));

                assert(edges[edge.child_edge1].tgt == edges[edge.child_edge2].src);
                const NodeID centerNode_id = edges[edge.child_edge1].tgt;
                centerNodes.push_back(centerNode_id);

                centerNodes.splice(centerNodes.end(), getCenterNodes(nodes, edges, edge.child_edge2));
            }
            return centerNodes;
        }

	static double otherMeasure(const vector<CHNode> &nodes, const vector<CHEdge> &edges, const EdgeID edge_id) {
		const CHEdge& edge = edges[edge_id];

		std::list<NodeID> center_nodes = getCenterNodes(nodes, edges, edge_id);
		const CHNode& src_node = nodes[edge.src];
		const CHNode& tgt_node = nodes[edge.tgt];		 

		double maxProportion = 0;
		for (NodeID node_id: center_nodes) {
		    const CHNode& center_node = nodes[node_id];
		    double proportion = geo::calcPerpendicularLength(src_node, tgt_node, center_node);
		    //double proportion = geo::getTriangleProportion(src_node, tgt_node, center_node);
		    if (proportion > maxProportion) {
		        maxProportion = proportion;
		    }
		}
		return maxProportion;		
    	}
	
	static void expandEdge(const vector<CHNode> &nodes, vector<CHEdge> &edges, const EdgeID edgeID, bool expandDist, double expandSizeDist, bool expandOther, double expandSizeOther) {
		CHEdge &edge = edges[edgeID];
		if (edge.is_shortcut() && edge.remaining) {
			//const CHNode& src_node = nodes[edge.src];
			//const CHNode& tgt_node = nodes[edge.tgt];
			//const CHEdge& child_edge1 = edges[edge.child_edge1];
//			const CHNode center_node = nodes[child_edge1.tgt];
			//DEBUG("geo::getTriangleProportion)" <<  geo::getTriangleProportion(src_node, tgt_node, center_node));
			if ((expandDist && geo::geoDist(nodes[edge.src], nodes[edge.tgt]) > expandSizeDist)
				|| (expandOther && otherMeasure(nodes, edges, edgeID) > expandSizeOther)) {
				//if (edge.getDist(nodes) > expandSize) {
				//if (edge.dist > expandSize) {			
				
				//if (geo::calcPerpendicularLength(src_node, tgt_node, center_node) > expandSize) {
				
			
				//double BendingRatio = calcBendingRatio(nodes [edge.src], nodes[edge.getCenterPoint(edges)], nodes[edge.tgt]);
				//if (BendingRatio > expandSize) {
				edge.remaining = false;
				edges[edge.child_edge1].remaining = true;
				edges[edge.child_edge2].remaining = true;
				expandEdge(nodes, edges, edge.child_edge1, expandDist, expandSizeDist, expandOther, expandSizeOther);
				expandEdge(nodes, edges, edge.child_edge2, expandDist, expandSizeDist, expandOther, expandSizeOther);
				
			}
		}
	}
	/*
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
	}*/	 	
	
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
	    DEBUG("numberOfValidNodes" << numberOfValidNodes);
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
		double expandSizeDist = 0.002, //normal: 0.002
		bool expandDist = false,
		double expandSizeOther = 0.002, //normal: 0.002
		bool expandOther = false,
                bool spareShortcuts = false
		) {	

		DEBUG("markValidNodes");
		markValidNodes(ch_nodes, percent_of_showed_nodes);

		/*
		size_t nofValidNodes2 = 0;
		for(NodeID nodeID = 0; nodeID< (int) nodes.size(); nodeID++) {
			if(ch_nodes[nodeID].remaining) {
				nofValidNodes2++;
			}
		}
		DEBUG("numberOfValidNodes2" << nofValidNodes2);
		*/
			
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
                        if (ch_edges[edgeID].vis_unpleasing) {
                            ch_edges[edgeID].remaining = false;       
                        }				
                    }
                }
                
		
		if (expandDist || expandOther) {
			DEBUG("Processing expandSize");
			//process expandSize
			for (uint edgeID = 0; edgeID < ch_edges.size(); edgeID++){
				//if (ch_edges[edgeID].remaining) {
				expandEdge(ch_nodes, ch_edges, edgeID, expandDist, expandSizeDist, expandOther, expandSizeOther);
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
				Edge edge = Edge(it->src, it->tgt, width(it->type), color(it->type));
				edges.push_back(edge);
			}       
		}												
	}		
}
