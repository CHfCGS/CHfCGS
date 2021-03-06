#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "indexed_container.h"
//#include "geoFunctions.h"

#include <limits>
#include <vector>
#include <list>
#include <algorithm>

namespace chc
{

namespace unit_tests
{
	void testGraph();
}

template <typename NodeT, typename EdgeT>
class Graph
{
	protected:
		Metadata _meta_data;

		std::vector<NodeT> _nodes;

		std::vector<uint> _out_offsets;
		std::vector<uint> _in_offsets;
		std::vector<EdgeT> _out_edges;
		std::vector<EdgeT> _in_edges;
                
		/* Maps edge id to index in the _out_edge vector. */
		std::vector<uint> _id_to_index;

		EdgeID edge_count = 0;
                uint max_street_type = 0;

                void fixEdgeType();
		void sortInEdges();
		void sortOutEdges();
		void initOffsets();
		void initIdToIndex();

		void update();

	public:
		typedef NodeT node_type;
		typedef EdgeT edge_type;

		typedef EdgeSortSrcTgt<EdgeT> OutEdgeSort;
		typedef EdgeSortTgtSrc<EdgeT> InEdgeSort;

		/* Init the graph from file 'filename' and sort
		 * the edges according to OutEdgeSort and InEdgeSort. */
		void init(GraphInData<NodeT,EdgeT>&& data);

		void printInfo() const;
		template<typename Range>
		void printInfo(Range&& nodes) const;
                
                void setVisualFlag(EdgeID edge_id, bool tag);
		uint getNrOfNodes() const { return _nodes.size(); }
		uint getNrOfEdges() const { return _out_edges.size(); }
		EdgeT const& getEdge(EdgeID edge_id) const { return _out_edges[_id_to_index[edge_id]]; }
		NodeT const& getNode(NodeID node_id) const { return _nodes[node_id]; }
                std::vector<EdgeT> const& getAllEdges() {return _out_edges;};
		Metadata const& getMetadata() const { return _meta_data; }
                
                uint getMaxStreetType() const {return max_street_type;} 
		uint getNrOfEdges(NodeID node_id) const;
		uint getNrOfEdges(NodeID node_id, EdgeType type) const;                                
                double getLat(NodeID node_id) const;
                double getLon(NodeID node_id) const;
                std::list<NodeID> nodeNeighbours(NodeID node_id) const;                
                std::list<NodeID> nodeNeighbours(NodeID node_id, StreetType type) const;
                std::list<NodeID> nodeNeighbours(NodeID node_id, StreetType streetType, EdgeType edgeDirection) const;
                StreetType getMinStreetType(NodeID node_id) const;   
                std::vector<StreetType> getStreetTypeVector() const;
                bool degree_leq(NodeID node_id , uint degree) const;
                bool isActive(NodeID node_id) const;
                bool isOneway(NodeID node_id) const;

		typedef range<typename std::vector<EdgeT>::const_iterator> node_edges_range;
		node_edges_range nodeEdges(NodeID node_id, EdgeType type) const;

		friend void unit_tests::testGraph();
};

/*
 * Graph member functions.
 */

template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::init(GraphInData<NodeT, EdgeT>&& data)
{
	_meta_data.swap(data.meta_data);
	_nodes.swap(data.nodes);
	_out_edges.swap(data.edges);
	_in_edges = _out_edges;
	edge_count = _out_edges.size();

	update();

	Print("Graph info:");
	Print("===========");
	printInfo();
        
        fixEdgeType();
}

template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::printInfo() const
{
	printInfo(counting_iteration(range<NodeID>(0, _nodes.size())));
}

template <typename NodeT, typename EdgeT>
template <typename Range>
void Graph<NodeT, EdgeT>::printInfo(Range&& nodes) const
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
void Graph<NodeT, EdgeT>::setVisualFlag(EdgeID edge_id, bool tag) { 
    _out_edges[_id_to_index[edge_id]].vis_unpleasing = tag;
    /*
        if (tag) {
            _out_edges[_id_to_index[edge_id]].speed = 1000; 
        }else {
            _out_edges[_id_to_index[edge_id]].speed = -1000; 
        }    */

}
    
//change edge types so that smaller numbers indicate higher importance
template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::fixEdgeType() { 
    for (EdgeT& edge: _out_edges) {                
        if (edge.type < 1) {
            edge.type = 15;
        }
        //red lanes
        else if (edge.type == 9 || edge.type == 10) {
            edge.type = 4; 
        }
        max_street_type = std::max(max_street_type, edge.type);
    }    
    Print("maxStreetType: " << max_street_type);
}


template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::sortInEdges()
{
	Debug("Sort the incomming edges.");

	std::sort(_in_edges.begin(), _in_edges.end(), InEdgeSort());
	debug_assert(std::is_sorted(_in_edges.begin(), _in_edges.end(), InEdgeSort()));
}

template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::sortOutEdges()
{
	Debug("Sort the outgoing edges.");

	std::sort(_out_edges.begin(), _out_edges.end(), OutEdgeSort());
	debug_assert(std::is_sorted(_out_edges.begin(), _out_edges.end(), OutEdgeSort()));
}

template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::initOffsets()
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
	for (NodeID i(0); i<nr_of_nodes; i++) {
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

template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::initIdToIndex()
{
	Debug("Renew the index mapper.");
	_id_to_index.resize(edge_count);               
	for (uint i(0), size(_out_edges.size()); i<size; i++) {             
		_id_to_index[_out_edges[i].id] = i;
	}
}

template <typename NodeT, typename EdgeT>
void Graph<NodeT, EdgeT>::update()
{
	sortOutEdges();
	sortInEdges();
	initOffsets();
	initIdToIndex();
}

template <typename NodeT, typename EdgeT>
uint Graph<NodeT, EdgeT>::getNrOfEdges(NodeID node_id) const
{
	return getNrOfEdges(node_id, EdgeType::OUT) + getNrOfEdges(node_id, EdgeType::IN);
}

template <typename NodeT, typename EdgeT>
uint Graph<NodeT, EdgeT>::getNrOfEdges(NodeID node_id, EdgeType type) const
{
	if (type == EdgeType::IN) {
		return _in_offsets[node_id+1] - _in_offsets[node_id];
	}
	else {
		return _out_offsets[node_id+1] - _out_offsets[node_id];
	}
}

template <typename NodeT, typename EdgeT>
std::list<NodeID> Graph<NodeT, EdgeT>::nodeNeighbours(NodeID node_id) const
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
std::list<NodeID> Graph<NodeT, EdgeT>::nodeNeighbours(NodeID node_id, StreetType type) const
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
std::list<NodeID> Graph<NodeT, EdgeT>::nodeNeighbours(NodeID node_id, StreetType streetType, EdgeType edgeDirection) const
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
StreetType Graph<NodeT, EdgeT>::getMinStreetType(NodeID node_id) const
{
        StreetType MinStreetType = getMaxStreetType();
        for (auto const& edge : nodeEdges(node_id, EdgeType::IN)) {
            if (edge.type < MinStreetType) {
                MinStreetType = edge.type;
            }
        }
        for (auto const& edge : nodeEdges(node_id, EdgeType::OUT)) {
            if (edge.type < MinStreetType) {
                MinStreetType = edge.type;
            }
        }
        return MinStreetType;
}

template <typename NodeT, typename EdgeT>
std::vector<StreetType> Graph<NodeT, EdgeT>::getStreetTypeVector() const
{
    std::vector<StreetType> streetTypes;
    streetTypes.resize(getNrOfNodes());
    for(NodeID node_id = 0; node_id < getNrOfNodes(); node_id++) {
        streetTypes[node_id] = getMinStreetType(node_id);
    }
    return streetTypes;
}

template <typename NodeT, typename EdgeT>
bool Graph<NodeT, EdgeT>::isActive(NodeID node_id) const
{
    return getNrOfEdges(node_id)!=0;
}

template <typename NodeT, typename EdgeT>
bool Graph<NodeT, EdgeT>::degree_leq(NodeID node_id , uint degree) const
{
    if (getNrOfEdges(node_id, EdgeType::OUT)<=degree 
            && getNrOfEdges(node_id, EdgeType::IN)<=degree) { //condition optional but should increase performance       
        if (nodeNeighbours(node_id).size()<=degree) {
            return true;
        }
    }
    return false;
}

template <typename NodeT, typename EdgeT>
bool Graph<NodeT, EdgeT>::isOneway(NodeID node_id) const
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
double Graph<NodeT, EdgeT>::getLat(NodeID node_id) const
{
    debug_assert(0 <= node_id && node_id < _nodes.size());
    return _nodes[node_id].lat;
}

template <typename NodeT, typename EdgeT>
double Graph<NodeT, EdgeT>::getLon(NodeID node_id) const
{
    debug_assert(0 <= node_id && node_id < _nodes.size());
    return _nodes[node_id].lon;
}

template <typename NodeT, typename EdgeT>
auto Graph<NodeT, EdgeT>::nodeEdges(NodeID node_id, EdgeType type) const -> node_edges_range {
	if (EdgeType::OUT == type) {
		return node_edges_range(_out_edges.begin() + _out_offsets[node_id], _out_edges.begin() + _out_offsets[node_id+1]);
	} else {
		return node_edges_range(_in_edges.begin() + _in_offsets[node_id], _in_edges.begin() + _in_offsets[node_id+1]);
	}
}

}
