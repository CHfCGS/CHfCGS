#pragma once

typedef int NodeID;
typedef int EdgeID;

typedef uint Level;

struct Node
{
	Node() : lat(0), lon(0) {}
	Node(double la, double lo) :
	lat(la), lon(lo) {}

	double lat;
	double lon;
};

struct Edge
{
	Edge() : source(0), target(0), width(0.0), color(0) {}
	Edge(uint s, uint t, uint w, uint c) :
	source(s), target(t), width(w), color(c) {}

	uint source;
	uint target;
	uint width;
	int color;
};

struct CHNode{
	//NodeID ID;
    //OSMID osmID;
    double lat;
    double lon;
    //int elev;
    Level lvl;
    bool remaining;
};

double pythagoras(double a, double b) {
	return sqrt(pow(a, 2) + pow(b, 2));
}

double geoDist(const CHNode &node1, const CHNode &node2) {
	//2D plane
	double dlat = node1.lat- node2.lat;
	double dlon = node1.lon -node2.lon;
	return pythagoras(dlat, dlon);    
}

struct CHEdge {
	NodeID src;
    NodeID tgt;
    uint dist;
    uint type;
    int speed;
    EdgeID child_edge1;
    EdgeID child_edge2;
    bool remaining;
    
    bool is_shortcut(){return (child_edge1 !=-1 && child_edge2 != -1);}	                   
    bool is_valid(const std::vector<CHNode> &nodes) {
        return (nodes[src].remaining && nodes[tgt].remaining);
    }    
    double getDist(const std::vector<CHNode> &nodes) {return geoDist(nodes[src], nodes[tgt]);}
};
