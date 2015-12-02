#pragma once
#include "defs.h"
#include "enum_helpers.h"
#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <limits>

//#include "geoFunctions.h"

namespace chm
{
    
enum class EdgeType : uint8_t {OUT = 0, IN = 1};
inline EdgeType operator!(EdgeType type) { return to_enum<EdgeType>(1 - from_enum(type)); }

typedef std::map<std::string, std::string> Metadata;

typedef int NodeID;
typedef int EdgeID;
typedef uint StreetType;

typedef uint Level;

struct Node {

    Node() : lat(0), lon(0) {
    }

    Node(double la, double lo) :
    lat(la), lon(lo) {
    }

    double lat;
    double lon;
};

struct Edge {

    Edge() : source(0), target(0), width(0.0), color(0) {
    }

    Edge(uint s, uint t, uint w, uint c) :
    source(s), target(t), width(w), color(c) {
    }

    uint source;
    uint target;
    uint width;
    int color;
};

struct CHNode {
    CHNode(): lat(0), lon(0), lvl(0), remaining(false) {}
    
    //care: lat and lon parameter is switched 
    CHNode(double lon, double lat): lat(lat), lon(lon), lvl(0), remaining(false), shortcuts(std::list<EdgeID>()) {} 
    //NodeID ID;
    //OSMID osmID;
    double lat;
    double lon;
    //int elev;
    Level lvl;
    bool remaining;
    std::list<EdgeID> shortcuts; //all shortcuts edges which use this node as centerNode
};
/*
double pythagoras(double a, double b) {
    return sqrt(pow(a, 2) + pow(b, 2));
}

double geoDist(const CHNode &node1, const CHNode &node2) {
    //2D plane
    double dlat = node1.lat - node2.lat;
    double dlon = node1.lon - node2.lon;
    return pythagoras(dlat, dlon);
}
*/
struct CHEdge {  
    EdgeID id;
    NodeID src;
    NodeID tgt;
    uint dist;
    uint type;
    //int speed;
    EdgeID child_edge1;
    EdgeID child_edge2;
    bool remaining;
    
    uint distance() const { return dist; }
/*
    bool is_shortcut() {
        return (child_edge1 != -1 && child_edge2 != -1);
    }
 * */
    /*
    bool is_valid(const std::vector<CHNode> &nodes) {
        return (nodes[src].remaining && nodes[tgt].remaining);
    }
    */
    /*
    double getDist(CHEdge ch_edge) {
        return geo::geoDist(nodes[src], nodes[tgt]);
    }    
     * */
    
};

inline NodeID otherNode(CHEdge const& edge, EdgeType edge_type) {
	switch (edge_type) {
	case EdgeType::OUT:
		return edge.tgt;
	case EdgeType::IN:
		break;
	}
	return edge.src;
}

namespace c
{
	uint const NO_NID(std::numeric_limits<NodeID>::max());
	uint const NO_EID(std::numeric_limits<EdgeID>::max());
	uint const NO_DIST(std::numeric_limits<uint>::max());
	uint const NO_LVL(std::numeric_limits<uint>::max());
}

struct ValidLevelInfo {
    Level allValidUntilLevel;
    size_t validNofNodesOnThatLevel;
};

template <typename NodeT, typename EdgeT>
struct GraphInData {
    std::vector<NodeT> nodes;
    std::vector<EdgeT> edges;
    Metadata meta_data;
};


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


typedef std::list<NodeID> Chain;


//similar to CHEdge, but can be handled easier for geometric calculations
struct CHLine {
    CHLine(CHNode src, CHNode tgt):
        src(src), tgt(tgt) {}
    CHNode src;
    CHNode tgt;
};

//ILP data structures

typedef uint LineID;

//wrapper for normal ch_node
struct ILP_Node {
    ILP_Node(const CHNode &ch_node, const NodeID ch_node_id, const uint posInChain, bool side):
        ch_node(ch_node), ch_node_id(ch_node_id), posInChain(posInChain), side(side) {}
    const CHNode &ch_node;
    const NodeID ch_node_id; //NodeID of ch_node in graph
    const uint posInChain;
    bool side; //true == to, false == from     
};

typedef std::list<ILP_Node> ILP_Chain;

//potential edges
struct Line {
    Line(const ILP_Node start, const ILP_Node end, LineID id):
        start(start), end(end), id(id) {}    
    const ILP_Node start;
    const ILP_Node end;
    LineID id;    
};

struct Intersection {
    Intersection(const Line line1, const Line line2):
        line1(line1), line2(line2) {}    
    const Line line1;
    const Line line2;
};


}
