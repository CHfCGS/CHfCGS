#pragma once

#include "defs.h"
#include "enum_helpers.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <map>
#include <string>
#include <list>
#include <functional>

namespace chc
{

namespace unit_tests
{
	void testNodesAndEdges();
}

typedef uint NodeID;
typedef uint EdgeID;
typedef uint StreetType;

namespace c
{
	uint const NO_NID(std::numeric_limits<NodeID>::max());
	uint const NO_EID(std::numeric_limits<EdgeID>::max());
	uint const NO_DIST(std::numeric_limits<uint>::max());
	uint const NO_LVL(std::numeric_limits<uint>::max());
}

typedef std::map<std::string, std::string> Metadata;

enum class EdgeType : uint8_t {OUT = 0, IN = 1};
inline EdgeType operator!(EdgeType type) { return to_enum<EdgeType>(1 - from_enum(type)); }

//ways
typedef uint WayID;

struct WayInfo {    
    bool undirected;
    uint waylength;
};

template <typename NodeT>
struct DirectedWay {
    DirectedWay() : nodes() {}
    //uint length;
    std::list<std::reference_wrapper<NodeT>> nodes;
    //std::list<NodeT> nodes;
    //uint32_t startEdge;
    //uint32_t endEdge;
    //uint32_t startNode;
    //uint32_t endNode;
    //int32_t oppositeWay;
    //double distToOpp;
};
/*
struct Chain {
    Chain() : node_ids() {}
    
    std::list<NodeID> node_ids;
    StreetType type = 0;
};

struct Chains_and_Remainder {
    Chains_and_Remainder() : chains(), remainder() {}
        
    std::list<Chain> chains;
    std::vector<NodeID> remainder;
};
*/

/*
 * Data required to construct or write out graphs.
 */
template <typename NodeT, typename EdgeT>
struct GraphInData {
	/* graph will "steal" data */
	std::vector<NodeT> nodes;
	std::vector<EdgeT> edges;
        std::vector<WayInfo> wayInfos;
        std::vector<DirectedWay<NodeT> > directedWays;
	Metadata meta_data;
};

template <typename NodeT, typename EdgeT>
struct GraphOutData {
	std::vector<NodeT> const& nodes;
	std::vector<EdgeT> const& edges;
	Metadata meta_data;
};

template <typename NodeT, typename EdgeT>
struct GraphCHOutData {
	std::vector<NodeT> const& nodes;
	std::vector<uint> const& node_levels;
	std::vector<EdgeT> const& edges;
	Metadata meta_data;
};


/*
 * Nodes
 */

template<typename NodeT>
struct CHNode : NodeT
{
	typedef NodeT base_node_type;

	uint lvl = c::NO_LVL;

	explicit CHNode() { }
	explicit CHNode(NodeT const& node) : NodeT(node) { }
	explicit CHNode(NodeT const& node, uint lvl) : NodeT(node), lvl(lvl){}
};

template<typename NodeT>
inline CHNode<NodeT> makeCHNode(NodeT const& node, uint lvl) {
	return CHNode<NodeT>(node, lvl);
}

template <typename NodeT>
struct _MakeCHNode
{
	typedef CHNode<NodeT> type;
};

template <typename NodeT>
struct _MakeCHNode<CHNode<NodeT>>
{
	typedef CHNode<NodeT> type;
};

template <typename NodeT>
using MakeCHNode = typename _MakeCHNode<typename std::remove_reference<NodeT>::type>::type;
/* remove possible CHNode<> with: typename MakeCHNode<T>::base_node_type */

struct Node
{
	NodeID id = c::NO_NID;

	explicit Node() { }
	explicit Node(NodeID id) : id(id) { }

	bool operator<(Node const& node) const { return id < node.id; }
};

struct GeoNode
{
	NodeID id = c::NO_NID;
	double lat = 0;
	double lon = 0;
	int elev = 0;

	explicit GeoNode() { }
	explicit GeoNode(NodeID id, double lat, double lon, int elev)
		: id(id), lat(lat), lon(lon), elev(elev) { }

	explicit operator Node() const
	{
		return Node(id);
	}

	bool operator<(GeoNode const& other) const { return id < other.id; }
};

struct StefanNode
{
	NodeID id = c::NO_NID;
	uint64_t osm_id = std::numeric_limits<uint>::max();
	double lat = 0;
	double lon = 0;

	StefanNode() { }
	StefanNode(NodeID id, uint64_t osm_id, double lat, double lon)
		: id(id), osm_id(osm_id), lat(lat), lon(lon) {}

	explicit operator Node() const
	{
		return Node(id);
	}

	operator GeoNode() const
	{
		return GeoNode(id, lat, lon, 0);
	}

	bool operator<(StefanNode const& other) const { return id < other.id; }
};

struct WayOSMNode
{       typedef std::list<std::reference_wrapper<DirectedWay<WayOSMNode>>> waylist;
        //std::list<std::reference_wrapper<DirectedWay<WayOSMNode>>> nodes;
        
        
        NodeID id = c::NO_NID;
	uint64_t osm_id = std::numeric_limits<uint>::max();
	double lat = 0;
	double lon = 0;
	int elev = 0;
        
        //osmways
        waylist ways = waylist();
        /*
        //2-degree chains
        bool degree_leq_two = false;
        NodeID neighbor1 = c::NO_NID;
        NodeID neighbor2 = c::NO_NID;
        bool collected = false;
         * */
        
        WayOSMNode() { }
        WayOSMNode(NodeID id, uint64_t osm_id, double lat, double lon, int elev)
		: id(id), osm_id(osm_id), lat(lat), lon(lon), elev(elev) {}
	WayOSMNode(NodeID id, uint64_t osm_id, double lat, double lon, int elev, waylist ways)
		: id(id), osm_id(osm_id), lat(lat), lon(lon), elev(elev), ways(ways) {}
        

	explicit operator Node() const
	{
		return Node(id);
	}

	operator GeoNode() const
	{
		return GeoNode(id, lat, lon, elev);
	}

	operator StefanNode() const
	{
		return StefanNode(id, osm_id, lat, lon);
	}
        /*
        operator OSMNode() const
	{
		return OSMNode(id, osm_id, lat, lon, elev);
	}
         * */
        
        //std::list<DirectedWay<WayOSMNode> > ways;
};

struct OSMNode
{
	NodeID id = c::NO_NID;
	uint64_t osm_id = std::numeric_limits<uint>::max();
	double lat = 0;
	double lon = 0;
	int elev = 0;
        
        OSMNode() { }
	OSMNode(NodeID id, uint64_t osm_id, double lat, double lon, int elev)
		: id(id), osm_id(osm_id), lat(lat), lon(lon), elev(elev) {}

	explicit operator Node() const
	{
		return Node(id);
	}

	operator GeoNode() const
	{
		return GeoNode(id, lat, lon, elev);
	}
        
	operator StefanNode() const
	{
		return StefanNode(id, osm_id, lat, lon);
	}    
        operator WayOSMNode() const
	{
		return WayOSMNode(id, osm_id, lat, lon, elev);
	}   

	bool operator<(OSMNode const& other) const { return id < other.id; }
};

/*
 * Edges
 */

template<typename EdgeT>
inline NodeID otherNode(EdgeT const& edge, EdgeType edge_type) {
	switch (edge_type) {
	case EdgeType::OUT:
		return edge.tgt;
	case EdgeType::IN:
		break;
	}
	return edge.src;
}

template<typename EdgeT1, typename EdgeT2>
inline bool equalEndpoints(EdgeT1 const& edge1, EdgeT2 const& edge2) {
	return edge1.src == edge2.src && edge1.tgt == edge2.tgt;
}

template<typename EdgeT>
struct CHEdge : EdgeT
{
	typedef EdgeT base_edge_type;
        
	EdgeID child_edge1 = c::NO_EID;
	EdgeID child_edge2 = c::NO_EID;
	NodeID center_node = c::NO_NID;

	CHEdge() { }
	CHEdge(EdgeT const& edge): EdgeT(edge) { }
	CHEdge(EdgeT const& edge, EdgeID child_edge1, EdgeID child_edge2, NodeID center_node)
		: EdgeT(edge), child_edge1(child_edge1),
		child_edge2(child_edge2), center_node(center_node){}        
};

template<typename EdgeT>
struct _MakeCHEdge
{
	typedef CHEdge<EdgeT> type;
};

template<typename EdgeT>
struct _MakeCHEdge<CHEdge<EdgeT>>
{
	typedef CHEdge<EdgeT> type;
};

template<typename EdgeT>
using MakeCHEdge = typename _MakeCHEdge<typename std::remove_reference<EdgeT>::type>::type;
/* remove possible CHEdge<> with: typename MakeCHEdge<T>::base_edge_type */

template<typename EdgeT>
MakeCHEdge<EdgeT> make_shortcut(EdgeT const& edge1, EdgeT const& edge2) {
	typedef typename MakeCHEdge<EdgeT>::base_edge_type base;
	assert(edge1.tgt == edge2.src);                
	return MakeCHEdge<EdgeT>(concat(static_cast<base const&>(edge1), static_cast<base const&>(edge2)), edge1.id, edge2.id, edge1.tgt);
}

struct Edge
{
	EdgeID id = c::NO_EID;
	NodeID src = c::NO_NID;
	NodeID tgt = c::NO_NID;
	uint dist = c::NO_DIST;

	Edge() { }
	Edge(EdgeID id, NodeID src, NodeID tgt, uint dist)
		: id(id), src(src), tgt(tgt), dist(dist) { }

	uint distance() const { return dist; }
};
Edge concat(Edge const& edge1, Edge const& edge2);

template <typename EdgeT>
struct MetricEdge : EdgeT
{
	uint metric = 0;

	MetricEdge() {}
	MetricEdge(EdgeT const& edge) : EdgeT(edge) {}
	MetricEdge(EdgeT const& edge, uint metric) : EdgeT(edge), metric(metric){}

	friend MetricEdge concat(MetricEdge const& edge1, MetricEdge const& edge2) {
		return MetricEdge(concat(static_cast<EdgeT const&>(edge1), static_cast<EdgeT const&>(edge2)), edge1.metric + edge2.metric);
	}
};

struct StefanEdge
{
	EdgeID id = c::NO_EID;
	NodeID src = c::NO_NID;
	NodeID tgt = c::NO_NID;
	uint dist = c::NO_DIST;

	StefanEdge() { }
	StefanEdge(EdgeID id, NodeID src, NodeID tgt, uint dist)
	: id(id), src(src), tgt(tgt), dist(dist) { }

	uint distance() const { return dist; }

	operator Edge() const
	{
		return Edge(id, src, tgt, dist);
	}
};
StefanEdge concat(StefanEdge const& edge1, StefanEdge const& edge2);

struct OSMEdge
{
	EdgeID id = c::NO_EID;
	NodeID src = c::NO_NID;
	NodeID tgt = c::NO_NID;
	uint dist = c::NO_DIST;
	StreetType type = 0;
	int speed = -1;

	OSMEdge() { }
	OSMEdge(EdgeID id, NodeID src, NodeID tgt, uint dist, StreetType type, int speed)
	: id(id), src(src), tgt(tgt), dist(dist), type(type), speed(speed) { }

	uint distance() const { return dist; }

	operator Edge() const
	{
		return Edge(id, src, tgt, dist);
	}

	operator StefanEdge() const
	{
		return StefanEdge(id, src, tgt, dist);
	}
};
OSMEdge concat(OSMEdge const& edge1, OSMEdge const& edge2);

struct OSMDistEdge : OSMEdge { };

struct EuclOSMEdge : OSMEdge
{
	uint eucl_dist = c::NO_DIST;

	EuclOSMEdge() { }
	EuclOSMEdge(EdgeID id, NodeID src, NodeID tgt, uint dist, uint type, int speed, uint eucl_dist)
	: OSMEdge(id, src, tgt, dist, type, speed), eucl_dist(eucl_dist) { }
};
EuclOSMEdge concat(EuclOSMEdge const& edge1, EuclOSMEdge const& edge2);

/*
 * EdgeSort
 */

template <typename EdgeT>
struct EdgeSortSrcTgt
{
	bool operator()(EdgeT const& edge1, EdgeT const& edge2) const
	{
		return edge1.src < edge2.src ||
		       (edge1.src == edge2.src && edge1.tgt < edge2.tgt);
	}
};

template <typename EdgeT>
struct EdgeSortTgtSrc
{
	bool operator()(EdgeT const& edge1, EdgeT const& edge2) const
	{
		return edge1.tgt < edge2.tgt ||
		       (edge1.tgt == edge2.tgt && edge1.src < edge2.src);
	}
};

template <typename EdgeT>
struct EdgeSortSrcTgtDist
{
	bool operator()(EdgeT const& edge1, EdgeT const& edge2) const
	{
		return edge1.src < edge2.src ||
		       (edge1.src == edge2.src && edge1.tgt < edge2.tgt) ||
		       (edge1.src == edge2.src && edge1.tgt == edge2.tgt && edge1.dist < edge2.dist);
	}
};

}

//similar to an Edge, but can be handled easier for geometric calculations
struct OSMLine {
    OSMLine(chc::OSMNode src, chc::OSMNode tgt):
        src(src), tgt(tgt) {}
    chc::OSMNode src;
    chc::OSMNode tgt;
};

/*
namespace ls {

using namespace chc;
//DouglasPeucker related types


struct simplePrioNode {
    const NodeID node_id;
    double perpendicularLength;
    uint nOfIntersections = 0;
    
    simplePrioNode(NodeID node_id, double perpendicularLength, uint nOfIntersections):
        node_id(node_id), perpendicularLength(perpendicularLength), nOfIntersections(nOfIntersections) {}
};

struct PrioNode;
//struct Intervall;

//typedef std::reference_wrapper<Node> NodeRef;
//typedef std::list<NodeRef> NodeRefList;

//Liste von Zeigern auf Einträge in anderer Liste
//typedef list<list<PrioNode>::iterator> QueueNodeIteratorList;
typedef std::reference_wrapper<PrioNode> PrioNodeRef;

struct BoundingBox {
    double upperLat;
    double LeftLon;   
    double lowerLat;
    double RightLon;
    double midLat;
    double midLon;
    
    BoundingBox(double lat1, double lon1, double lat2, double lon2) {
        if (lat1 > lat2) {
            upperLat = lat1;
            lowerLat = lat2;
        } else {
            upperLat = lat2;
            lowerLat = lat1;
        }
        if (lon1 > lon2) {
            RightLon = lon1;
            LeftLon = lon2;
        } else {
            RightLon = lon2;
            LeftLon = lon1;
        }
        midLat = (upperLat + lowerLat) /2;
        midLon = (LeftLon + RightLon) /2;
    }
    
    bool contains(double lat, double lon) {
        return ((lowerLat < lat) && (lat < upperLat)
                && (LeftLon < lon) && (lon < RightLon));
    }
    
};

struct NodeInBox {
    NodeInBox():nodeID(0), sign(false) {}
    NodeInBox(NodeID nodeID, bool sign): nodeID(nodeID), sign(sign) {}
    NodeID nodeID;
    bool sign; // true: >= 0 false: < 0
};


struct Intervall {    

    Intervall(const NodeID start, const NodeID end): start(start), end(end), prioNodesRefs() {        
    }    
    const NodeID start;
    const NodeID end;
    std::list<PrioNodeRef> prioNodesRefs; //all nodes between start and end (excluding))    
};


struct PrioNode {
    const NodeID node_id;
    std::list<Intervall>::iterator intervallIt;    
    std::list<PrioNodeRef>::iterator posInIntervallIt;
        
    std::vector<NodeInBox> &leftBox;
    std::vector<NodeInBox> &rightBox;
    
    double perpendicularLength;
    uint nOfIntersections = 0; //used as orientation misses
            
    PrioNode(const NodeID node, std::list<Intervall>::iterator intervallIt, std::vector<NodeInBox> &leftBox, std::vector<NodeInBox> &rightBox)
    : node_id(node), intervallIt(intervallIt), leftBox(leftBox), rightBox(rightBox), perpendicularLength(0) {        
    }

    explicit operator simplePrioNode() const {
        return simplePrioNode(node_id, perpendicularLength, nOfIntersections);
    }
    
    //< means later processed in DP and sooner contracted
    bool operator <(const PrioNode &rhs) const {
        
        if (this->nOfIntersections > rhs.nOfIntersections) {
            return true;
        } else if (this->nOfIntersections == rhs.nOfIntersections) {            
            return this->perpendicularLength < rhs.perpendicularLength;     
        }
        else {
            return false;
        }
        
        //return this->perpendicularLength < rhs.perpendicularLength; 
    }
};


struct Intervall2 {    

    Intervall2(const NodeID start, const NodeID end): start(start), end(end), prioNodesRefs() {        
    }    
    const NodeID start;
    const NodeID end;
    std::list<PrioNodeRef> prioNodesRefs; //all nodes between start and end (excluding))    
};


struct PrioNode2 {
    const NodeID node_id;
    std::list<Intervall2>::iterator intervallIt;    
    std::list<PrioNodeRef>::iterator posInIntervallIt;
        
    std::vector<NodeInBox> &leftBox;
    std::vector<NodeInBox> &rightBox;
    
    double perpendicularLength;
    uint nOfIntersections = 0; //used as orientation misses
            
    PrioNode2(const NodeID node, std::list<Intervall2>::iterator intervallIt, std::vector<NodeInBox> &leftBox, std::vector<NodeInBox> &rightBox)
    : node_id(node), intervallIt(intervallIt), leftBox(leftBox), rightBox(rightBox), perpendicularLength(0) {        
    }

    explicit operator simplePrioNode() const {
        return simplePrioNode(node_id, perpendicularLength, nOfIntersections);
    }
    
    //< means later processed in DP and sooner contracted
    bool operator <(const PrioNode2 &rhs) const {
        
        if (this->nOfIntersections > rhs.nOfIntersections) {
            return true;
        } else if (this->nOfIntersections == rhs.nOfIntersections) {            
            return this->perpendicularLength < rhs.perpendicularLength;     
        }
        else {
            return false;
        }
        
        //return this->perpendicularLength < rhs.perpendicularLength; 
    }
};

//data types used for ChainPair version of DouglasPeucker ("P" = "Pair")

struct PPrioNode;

struct PIntervall {    
    PIntervall(const NodeID start, const NodeID end): start(start), end(end), prioNodesIts() {        
    }
    
    const NodeID start;
    const NodeID end;
    //std::list<PPrioNodeRef> prioNodesIts; //all nodes between start and end (excluding))    
    std::list<std::list<PPrioNode>::iterator> prioNodesIts; //all nodes between start and end (excluding)), iterators in chain-global Priolist    
};

struct PPrioNode {
    std::list<PPrioNode>::iterator follower; //if a PPrioNode is chosen to subdivide one way the follower is chosen to subdivide the other in the DP algorithm
    bool followerValid = false;
    std::list<std::list<PPrioNode>::iterator> guides; // this object would follow all PPrioNodes in the list
    
    const NodeID node_id;
    std::list<PIntervall>::iterator intervallIt;    
    std::list<std::list<PPrioNode>::iterator>::iterator posInIntervallIt;
        
    std::vector<NodeInBox> &leftBox;
    std::vector<NodeInBox> &rightBox;
    
    double perpendicularLength;
    uint nOfIntersections = 0; //used as orientation misses
            
    PPrioNode(const NodeID node, std::list<PIntervall>::iterator intervallIt, std::vector<NodeInBox> &leftBox, std::vector<NodeInBox> &rightBox)
    : node_id(node), intervallIt(intervallIt), leftBox(leftBox), rightBox(rightBox), perpendicularLength(0) {        
    }    
    
    explicit operator simplePrioNode() const {
        return simplePrioNode(node_id, perpendicularLength, nOfIntersections);
    }
    
    //< means later processed in DP and sooner contracted
    bool operator <(const PPrioNode &rhs) const {
        //TODO also take follower values into account                
        if (this->nOfIntersections > rhs.nOfIntersections) {
            return true;
        } else if (this->nOfIntersections == rhs.nOfIntersections) {            
            return this->perpendicularLength < rhs.perpendicularLength;     
        }
        else {
            return false;
        }        
    }    
};

}
*/