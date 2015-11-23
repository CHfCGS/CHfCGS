#pragma once

#include "defs.h"
#include "indexed_container.h"
#include "geoFunctions.h"
#include "nodes_and_edges.h"

#include <vector>
#include <list>
#include <algorithm>

namespace chm
{

template <typename NodeT, typename EdgeT>
class CHGraph
{
	protected:
		Metadata _meta_data;

		std::vector<NodeT> _nodes;

		std::vector<uint> _out_offsets;
		std::vector<uint> _in_offsets;
		std::vector<EdgeT> _out_edges;
		std::vector<EdgeT> _in_edges;
                
                //marks all valid elements right now
                std::vector<bool> _validNodes;
                std::vector<bool> _validEdges;

		/* Maps edge id to index in the _out_edge vector. */
		//std::vector<uint> _id_to_index;

		EdgeID edge_count = 0;

		void sortInEdges();
		void sortOutEdges();
		void initOffsets();
		void initIdToIndex();

		void update();
                
                void setCenterNodesLists();
                
                //zoom helper functions
                bool isRemaining(EdgeID edge_id) const;
                void invalidateShortcut(const EdgeID edge_id);
                void expandEdge(const EdgeID edgeID, double expandSize);  
                void markValidNodes(double percent_of_valid_nodes);
                ValidLevelInfo calcValidLevel(size_t numberOfValidNodes);
                //void removeShortcut(const EdgeID edgeID);
                //void expandEdge(const EdgeID edgeID, double expandSize);
                

	public:
		typedef NodeT node_type;
		typedef EdgeT edge_type;

		typedef EdgeSortSrcTgt<EdgeT> OutEdgeSort;
		typedef EdgeSortTgtSrc<EdgeT> InEdgeSort;

		/* Init the graph from file 'filename' and sort
		 * the edges according to OutEdgeSort and InEdgeSort. */
		void init(GraphInData<NodeT,EdgeT>&& data);
                
                void zoom(double percent_of_valid_nodes, bool expand, double expandsize);

		void printInfo() const;
		template<typename Range>
		void printInfo(Range&& nodes) const;
                              
                bool isUp(EdgeT const& edge, EdgeType direction) const;
                
		uint getNrOfNodes() const { return _nodes.size(); }
		uint getNrOfEdges() const { return _out_edges.size(); }
		//EdgeT const& getEdge(EdgeID edge_id) const { return _out_edges[_id_to_index[edge_id]]; }
                EdgeT const& getEdge(EdgeID edge_id) const { return _out_edges[edge_id]; }
		NodeT const& getNode(NodeID node_id) const { return _nodes[node_id]; }
		Metadata const& getMetadata() const { return _meta_data; }
                                                                
		uint getNrOfEdges(NodeID node_id) const;
		uint getNrOfEdges(NodeID node_id, EdgeType type) const;
                double getEdgeDist(EdgeID edge_id) const;
                bool isValidNode(NodeID node_id) const;
                bool isValidEdge(EdgeID edge_id) const;
                bool isShortcut(EdgeID edge_id) const;                
                bool isHigh(EdgeID edge_id) const;
                
                
                double getLat(NodeID node_id) const;
                double getLon(NodeID node_id) const;
                std::list<NodeID> nodeNeighbours(NodeID node_id) const;                
                std::list<NodeID> nodeNeighbours(NodeID node_id, StreetType type) const;
                std::list<NodeID> nodeNeighbours(NodeID node_id, StreetType streetType, EdgeType edgeDirection) const;
                StreetType getMaxStreetType(NodeID node_id) const;
                bool degree_leq_two(NodeID node_id) const;
                bool isOneway(NodeID node_id) const;

		typedef range<typename std::vector<EdgeT>::const_iterator> node_edges_range;
		node_edges_range _nodeEdges(NodeID node_id, EdgeType type) const;
                std::vector<EdgeT> nodeEdges(NodeID node_id, EdgeType type) const;
                std::vector<EdgeT> nodeEdges(NodeID node_id, StreetType streetType) const;
                
		
};

/*
 * Graph member functions.
 */

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::init(GraphInData<NodeT, EdgeT>&& data)
{
	_meta_data.swap(data.meta_data);
	_nodes.swap(data.nodes);
	_out_edges.swap(data.edges);
	_in_edges = _out_edges;
	edge_count = _out_edges.size();

	update();
        setCenterNodesLists();
        zoom(0, false, 0);
	Print("Graph info:");
	Print("===========");
	printInfo();
}

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::printInfo() const
{
	printInfo(counting_iteration(range<NodeID>(0, _nodes.size())));
}

template <typename NodeT, typename EdgeT>
template <typename Range>
void CHGraph<NodeT, EdgeT>::printInfo(Range&& nodes) const
{
#ifdef NVERBOSE
	(void) nodes;
#else
	uint active_nodes(0);

	double avg_out_deg(0);
	double avg_in_deg(0);
	double avg_deg(0);

	std::vector<uint> out_deg;
	std::vector<uint> in_deg;
	std::vector<uint> deg;

	for (auto node: nodes) {
		uint out(getNrOfEdges(node, EdgeType::OUT));
		uint in(getNrOfEdges(node, EdgeType::IN));

		if (out != 0 || in != 0) {
			++active_nodes;

			out_deg.push_back(out);
			in_deg.push_back(in);
			deg.push_back(out + in);

			avg_out_deg += out;
			avg_in_deg += in;
			avg_deg += out+in;
		}
	}

	Print("#nodes: " << nodes.size() << ", #active nodes: " << active_nodes << ", #edges: " << _out_edges.size() << ", maximal edge id: " << edge_count - 1);

	if (active_nodes != 0) {
		auto mm_out_deg = std::minmax_element(out_deg.begin(), out_deg.end());
		auto mm_in_deg = std::minmax_element(in_deg.begin(), in_deg.end());
		auto mm_deg = std::minmax_element(deg.begin(), deg.end());

		avg_out_deg /= active_nodes;
		avg_in_deg /= active_nodes;
		avg_deg /= active_nodes;

		Print("min/max/avg degree:"
			<< " out "   << *mm_out_deg.first << " / " << *mm_out_deg.second << " / " << avg_out_deg
			<< ", in "   << *mm_in_deg.first  << " / " << *mm_in_deg.second  << " / " << avg_in_deg
			<< ", both " << *mm_deg.first     << " / " << *mm_deg.second     << " / " << avg_deg);
	}
	else {
		Debug("(no degree info is provided as there are no active nodes)");
	}
#endif
}

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::sortInEdges()
{
	Debug("Sort the incomming edges.");

	std::sort(_in_edges.begin(), _in_edges.end(), InEdgeSort());
	debug_assert(std::is_sorted(_in_edges.begin(), _in_edges.end(), InEdgeSort()));
}

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::sortOutEdges()
{
	Debug("Sort the outgoing edges.");

	std::sort(_out_edges.begin(), _out_edges.end(), OutEdgeSort());
	debug_assert(std::is_sorted(_out_edges.begin(), _out_edges.end(), OutEdgeSort()));
}

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::initOffsets()
{
	Debug("Init the offsets.");
	debug_assert(std::is_sorted(_out_edges.begin(), _out_edges.end(), OutEdgeSort()));
	debug_assert(std::is_sorted(_in_edges.begin(), _in_edges.end(), InEdgeSort()));

	uint nr_of_nodes(_nodes.size());

	_out_offsets.assign(nr_of_nodes + 1, 0);
	_in_offsets.assign(nr_of_nodes + 1, 0);

	/* assume "valid" edges are in _out_edges and _in_edges */
	for (auto const& edge: _out_edges) {
		_out_offsets[edge.src]++;
		_in_offsets[edge.tgt]++;
	}

	uint out_sum(0);
	uint in_sum(0);
	for (NodeID i(0); i<(int) nr_of_nodes; i++) {
		auto old_out_sum = out_sum, old_in_sum = in_sum;
		out_sum += _out_offsets[i];
		in_sum += _in_offsets[i];
		_out_offsets[i] = old_out_sum;
		_in_offsets[i] = old_in_sum;
	}
	assert(out_sum == _out_edges.size());
	assert(in_sum == _in_edges.size());
	_out_offsets[nr_of_nodes] = out_sum;
	_in_offsets[nr_of_nodes] = in_sum;
}
/*
template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::initIdToIndex()
{
	Debug("Renew the index mapper.");

	_id_to_index.resize(edge_count);
	for (uint i(0), size(_out_edges.size()); i<size; i++) {
            assert(_out_edges[i].id == i);
            _id_to_index[_out_edges[i].id] = i;
            assert(_id_to_index[_out_edges[i].id] == _out_edges[i].id);
	}
}
*/
template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::update()
{
	//sortOutEdges(); //out edges should already be sorted
	sortInEdges();
	initOffsets();
	//initIdToIndex();        
}

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::setCenterNodesLists() {
    for (EdgeID edge_id = 0; edge_id < edge_count; edge_id++) {
        if(isShortcut(edge_id)) {
            EdgeT child_edge1 = getEdge(getEdge(edge_id).child_edge1);
            EdgeT child_edge2 = getEdge(getEdge(edge_id).child_edge2);
            assert(child_edge1.tgt == child_edge2.src);
            _nodes[child_edge1.tgt].shortcuts.push_back(edge_id);
        }      
    }
}


template <typename NodeT, typename EdgeT>
ValidLevelInfo CHGraph<NodeT, EdgeT>::calcValidLevel(size_t numberOfValidNodes) {
        ValidLevelInfo vli;
        Level maxLevel = 0;
        for (NodeT node : _nodes) {
            if (node.lvl > maxLevel) {
                maxLevel = node.lvl;
            }
        }

        std::vector<uint> nofNodesPerLevel(maxLevel + 1, 0); //one field for level 0 is also needed
        for (NodeT node : _nodes) {
            nofNodesPerLevel[node.lvl]++;
        }

        assert(numberOfValidNodes >= 0 && numberOfValidNodes <= _nodes.size());        
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

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::markValidNodes(double percent_of_valid_nodes) { 
    assert(percent_of_valid_nodes >= 0 && percent_of_valid_nodes <= 100);
    size_t numberOfValidNodes = (uint) trunc((percent_of_valid_nodes/100.0) * _nodes.size());
    ValidLevelInfo vli = calcValidLevel(numberOfValidNodes);
    _validNodes.resize(_nodes.size());
    size_t nofValidatedNodesOnCriticalLevel = 0;
    for (NodeID nodeID = 0; nodeID < (int) _validNodes.size(); nodeID++) {
        NodeT &node = _nodes[nodeID];
        if (node.lvl > vli.allValidUntilLevel) {
            _validNodes[nodeID] = true;
        }else if(node.lvl == vli.allValidUntilLevel) {
            if(nofValidatedNodesOnCriticalLevel < vli.validNofNodesOnThatLevel) {
                nofValidatedNodesOnCriticalLevel++;
                _validNodes[nodeID] = true;
            }else {
                _validNodes[nodeID] = false;
            }        
        }else {
            _validNodes[nodeID] = false;
        }
    }
    assert(nofValidatedNodesOnCriticalLevel == vli.validNofNodesOnThatLevel);    
}

template <typename NodeT, typename EdgeT>
void CHGraph<NodeT, EdgeT>::zoom(double percent_of_valid_nodes, bool expand, double expandSize)
{
    markValidNodes(percent_of_valid_nodes);

    //mark valid edges
    _validEdges.resize(edge_count);
    for (EdgeID edgeID = 0; edgeID < (int) _validEdges.size(); edgeID++) {
        _validEdges[edgeID] = isHigh(edgeID);
    }
    
    //remove too high-level shortcuts    
    for (EdgeID edgeID = 0; edgeID < (int) edge_count; edgeID++) {
        invalidateShortcut(edgeID);
    }
    
    
    if (expand) {
        Debug("Processing expandSize");
        //process expandSize
        for (uint edgeID = 0; edgeID < _validEdges.size(); edgeID++) {
            if (_validEdges[edgeID]) {
                expandEdge(edgeID, expandSize);
            }
        }
    }

    #ifndef NDEBUG
    if (percent_of_valid_nodes==100) {
        for (int i = 0; i < _validNodes.size(); i++) {
            assert(_validNodes[i] == true);
        }
        for (int i = 0; i < edge_count; i++) {            
            assert(isValidEdge(i) != isShortcut(i));            
        }
    }        
    #endif
}

    template <typename NodeT, typename EdgeT>
    void CHGraph<NodeT, EdgeT>::invalidateShortcut(const EdgeID edge_id) {
        //CHEdge &edge = _out_edges[edgeID];
        if (isValidEdge(edge_id) && isShortcut(edge_id)) {
            
            if (isHigh(_out_edges[edge_id].child_edge1)) {
                _validEdges[edge_id] = false;
                invalidateShortcut(_out_edges[edge_id].child_edge1);
                debug_assert(isHigh(_out_edges[edge_id].child_edge2));
                invalidateShortcut(_out_edges[edge_id].child_edge2);
            }
            /*
            if (isHigh(_out_edges[edge_id].childedge2)) {
                _validEdges[edge_id] = false;
                invalidateShortcut(_out_edges[edge_id].childedge2);
            }
             * */
        }
    }
    
    template <typename NodeT, typename EdgeT>
    void CHGraph<NodeT, EdgeT>::expandEdge(const EdgeID edge_id, double expandSize) {
        //
        if (isShortcut(edge_id)) {
            if (getEdgeDist(edge_id) > expandSize) {
                //double BendingRatio = calcBendingRatio(nodes [edge.src], nodes[edge.getCenterPoint(edges)], nodes[edge.tgt]);
                //if (BendingRatio > expandSize) {
                
                _validEdges[edge_id] = false;
                
                EdgeID child_edge1_id = _out_edges[edge_id].child_edge1;
                EdgeID child_edge2_id = _out_edges[edge_id].child_edge2;
                _validEdges[child_edge1_id] = true;
                _validEdges[child_edge2_id] = true;
                expandEdge(child_edge1_id, expandSize);
                expandEdge(child_edge2_id, expandSize);
            }
        }
    }
    

template <typename NodeT, typename EdgeT>
uint CHGraph<NodeT, EdgeT>::getNrOfEdges(NodeID node_id) const
{
	return getNrOfEdges(node_id, EdgeType::OUT) + getNrOfEdges(node_id, EdgeType::IN);
}

template <typename NodeT, typename EdgeT>
uint CHGraph<NodeT, EdgeT>::getNrOfEdges(NodeID node_id, EdgeType type) const
{
	if (type == EdgeType::IN) {
		return _in_offsets[node_id+1] - _in_offsets[node_id];
	}
	else {
		return _out_offsets[node_id+1] - _out_offsets[node_id];
	}
}

template <typename NodeT, typename EdgeT>
double CHGraph<NodeT, EdgeT>:: getEdgeDist(EdgeID edge_id) const
{
    return geo::geoDist(_nodes[_out_edges[edge_id].src], _nodes[_out_edges[edge_id].tgt]);
}

template <typename NodeT, typename EdgeT>
std::list<NodeID> CHGraph<NodeT, EdgeT>::nodeNeighbours(NodeID node_id) const
{
        std::list<NodeID> neighbors;
        for (auto const& edge : nodeEdges(node_id, EdgeType::IN)) {
            neighbors.emplace_back(edge.src);
        }
        for (auto const& edge : nodeEdges(node_id, EdgeType::OUT)) {
            neighbors.emplace_back(edge.tgt);
        }
        neighbors.sort();
        neighbors.unique();
        return neighbors;
}

template <typename NodeT, typename EdgeT>
std::list<NodeID> CHGraph<NodeT, EdgeT>::nodeNeighbours(NodeID node_id, StreetType type) const
{
        std::list<NodeID> neighbors;
        for (auto const& edge : nodeEdges(node_id, EdgeType::IN)) {
            if (edge.type == type) {
                neighbors.emplace_back(edge.src);
            }            
        }
        for (auto const& edge : nodeEdges(node_id, EdgeType::OUT)) {
            if (edge.type == type) {
                neighbors.emplace_back(edge.tgt);
            }
        }
        neighbors.sort();
        neighbors.unique();
        return neighbors;
}

template <typename NodeT, typename EdgeT>
std::list<NodeID> CHGraph<NodeT, EdgeT>::nodeNeighbours(NodeID node_id, StreetType streetType, EdgeType edgeDirection) const
{
        std::list<NodeID> neighbors;
        for (auto const& edge : nodeEdges(node_id, edgeDirection)) {
            if (edge.type == streetType) {
                if (edgeDirection ==  EdgeType::IN) {
                    neighbors.emplace_back(edge.src);
                }
                else {
                    neighbors.emplace_back(edge.tgt);
                }      
            }            
        }
        
        //neighbors.sort();
        //neighbors.unique();
        return neighbors;
}

template <typename NodeT, typename EdgeT>
StreetType CHGraph<NodeT, EdgeT>::getMaxStreetType(NodeID node_id) const
{
        StreetType MaxStreetType = 0;
        for (auto const& edge : nodeEdges(node_id, EdgeType::IN)) {
            if (edge.type > MaxStreetType) {
                MaxStreetType = edge.type;
            }
        }
        for (auto const& edge : nodeEdges(node_id, EdgeType::OUT)) {
            if (edge.type > MaxStreetType) {
                MaxStreetType = edge.type;
            }
        }
        return MaxStreetType;
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::degree_leq_two(NodeID node_id) const
{
    if (getNrOfEdges(node_id, EdgeType::OUT)<=2 
            && getNrOfEdges(node_id, EdgeType::IN)<=2) { //condition optional but should increase performance       
        if (nodeNeighbours(node_id).size()<=2) {
            return true;
        }
    }
    return false;
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isOneway(NodeID node_id) const
{
    std::list<NodeID> neighbors;
    for (auto const& edge : nodeEdges(node_id, EdgeType::IN)) {        
        neighbors.emplace_back(edge.src);
    }
    for (auto const& edge : nodeEdges(node_id, EdgeType::OUT)) {        
        neighbors.emplace_back(edge.tgt);        
    }
    neighbors.sort();
    size_t sizeWithDuplicates = neighbors.size();
    neighbors.unique();
    size_t sizeWithOutDuplicates = neighbors.size();
    size_t nofDuplicates = sizeWithDuplicates - sizeWithOutDuplicates;
    return (nofDuplicates == 0);
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isValidNode(NodeID node_id) const
{
    debug_assert(0 <= node_id && node_id < (int) _validNodes.size());
    return _validNodes[node_id];
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isValidEdge(EdgeID edge_id) const
{
    debug_assert(0 <= edge_id && edge_id < (int) _validEdges.size());
    return _validEdges[edge_id];
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isUp(EdgeT const& edge, EdgeType direction) const
{
	uint src_lvl = _nodes[edge.src].lvl;
	uint tgt_lvl = _nodes[edge.tgt].lvl;

	if (src_lvl > tgt_lvl) {
		return direction == EdgeType::IN ? true : false;
	}
	else if (src_lvl < tgt_lvl) {
		return direction == EdgeType::OUT ? true : false;
	}

	/* should never reach this: */
	assert(src_lvl != tgt_lvl);
	return false;
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isShortcut(EdgeID edge_id) const
{
    return (_out_edges[edge_id].child_edge1 != -1 && _out_edges[edge_id].child_edge2 != -1);
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isHigh(EdgeID edge_id) const {
    return isValidNode(_out_edges[edge_id].src) && isValidNode(_out_edges[edge_id].tgt);
}

template <typename NodeT, typename EdgeT>
bool CHGraph<NodeT, EdgeT>::isRemaining(EdgeID edge_id) const
{
    debug_assert(0 <= edge_id && edge_id < _out_edges.size());    
    return isValidNode(_out_edges[edge_id].src) && isValidNode(_out_edges[edge_id].tgt);
}

template <typename NodeT, typename EdgeT>
double CHGraph<NodeT, EdgeT>::getLat(NodeID node_id) const
{
    debug_assert(0 <= node_id && node_id < (int) _nodes.size());
    return _nodes[node_id].lat;
}

template <typename NodeT, typename EdgeT>
double CHGraph<NodeT, EdgeT>::getLon(NodeID node_id) const
{
    debug_assert(0 <= node_id && node_id < (int) _nodes.size());
    return _nodes[node_id].lon;
}

template <typename NodeT, typename EdgeT>
auto CHGraph<NodeT, EdgeT>::_nodeEdges(NodeID node_id, EdgeType type) const -> node_edges_range {
	if (EdgeType::OUT == type) {
		return node_edges_range(_out_edges.begin() + _out_offsets[node_id], _out_edges.begin() + _out_offsets[node_id+1]);
	} else {
		return node_edges_range(_in_edges.begin() + _in_offsets[node_id], _in_edges.begin() + _in_offsets[node_id+1]);
	}
}

template <typename NodeT, typename EdgeT>
std::vector<EdgeT> CHGraph<NodeT, EdgeT>::nodeEdges(NodeID node_id, EdgeType type) const  {
        std::vector<EdgeT> nodeEdges;
        for (auto const edge : _nodeEdges(node_id, type)) {
            if (isValidEdge(edge.id)) {
                nodeEdges.push_back(edge);
            }
        }
        return nodeEdges;
}

template <typename NodeT, typename EdgeT>
std::vector<EdgeT> CHGraph<NodeT, EdgeT>::nodeEdges(NodeID node_id, StreetType streetType) const  {       
        std::vector<EdgeT> edges;
        for (auto const edge : nodeEdges(node_id, EdgeType::IN)) {
            if (edge.type == streetType) {
                edges.push_back(edge);                
            }            
        }        
        for (auto const edge : nodeEdges(node_id, EdgeType::OUT)) {
            if (edge.type == streetType) {
                edges.push_back(edge);                
            }            
        }                
        return edges;        
}

}
