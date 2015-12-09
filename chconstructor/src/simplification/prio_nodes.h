#pragma once

#include "../nodes_and_edges.h"
#include "boost/heap/binomial_heap.hpp"

namespace ls {

using namespace chc;

//data types for DouglasPeucker ----------------------------------------------------------------------------------------------
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

/*
struct PrioNodeData {
    PrioNodeData(PrioNodeRef ref):ref(ref), nodesInBox(){}
    PrioNodeRef ref;    
    std::vector<NodeInBox> nodesInBox;
};*/

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

//data used for more general simplification--------------------------------------------------------------------------------------
struct PrioNode2;
struct PrioNodeBase;



//typedef typename boost::heap::binomial_heap<PrioNodeBase> PrioNodeHeap;
typedef typename boost::heap::binomial_heap<PrioNode2> PrioNodeHeap;
typedef typename PrioNodeHeap::handle_type PrioNodeHandle;    


struct Intervall2 {    

    Intervall2(const NodeID start, const NodeID end): start(start), finish(end), prioNodeHandles() {        
    }    
    const NodeID start;
    const NodeID finish;
    std::list<PrioNodeHandle> prioNodeHandles; //all nodes between start and end (excluding))    
};

//base class which handles follower information and can be put in a heap
struct PrioNodeBase {    
    const NodeID node_id;
    
    double perpendicularLength = 0;
    uint nOfIntersections = 0; //used as orientation misses
    
    PrioNodeHandle follower_h; //if a PPrioNode is chosen to subdivide one way the follower is chosen to subdivide the other in the DP algorithm
    bool followerValid = false;
    std::list<PrioNodeHandle> guides; // this object would follow all PPrioNodes in the list    
    
    PrioNodeBase(const NodeID node_id): node_id(node_id) {}
    
    explicit operator simplePrioNode() const {
        return simplePrioNode(node_id, perpendicularLength, nOfIntersections);
    } 
    
    /*
    PrioNode2 operator =(const PrioNode2 &p_IN) {
        PrioNode2 p_out(p_IN.node_id, p_IN.intervallIt, p_IN.leftBox, p_IN.rightBox);
        p_out.posInIntervallIt = p_IN.posInIntervallIt;
        p_out.perpendicularLength = p_IN.perpendicularLength;
        p_out.nOfIntersections = p_IN.nOfIntersections;
        p_out.follower_h = p_IN.follower_h;
        p_out.followerValid = p_IN.followerValid;
        p_out.guides = p_IN.guides;
        
        return p_out;        
    }*/
    
};

struct PrioNodeTD: public PrioNodeBase {
    std::list<Intervall2>::iterator intervallIt;    
    std::list<PrioNodeHandle>::iterator posInIntervallIt;
    
    std::vector<NodeInBox> &leftBox;
    std::vector<NodeInBox> &rightBox;
    //PrioNodeTD(): PrioNode2()
     PrioNodeTD(const NodeID node, std::list<Intervall2>::iterator intervallIt,
                std::vector<NodeInBox> &leftBox, std::vector<NodeInBox> &rightBox
                //, std::list<PrioNode2>::iterator follower, bool followerValid,
                //std::list<std::list<PrioNode2>::iterator> guides
                    ): PrioNodeBase(node_id), intervallIt(intervallIt), leftBox(leftBox), rightBox(rightBox)
                    //, follower(follower), followerValid(followerValid), guides(guides)
    {}
      
};


typedef std::vector<NodeInBox> NodeBox;

struct PrioNode2 {        
    const NodeID node_id;
    std::list<Intervall2>::iterator intervallIt;    
    std::list<PrioNodeHandle>::iterator posInIntervallIt;
        
    NodeBox &leftBox;
    NodeBox &rightBox;
    
    //iterators in a list of all node boxes
    std::list<std::list<NodeBox>::iterator> left_node_boxes_its;
    std::list<std::list<NodeBox>::iterator> right_node_boxes_its;
    
    double perpendicularLength;
    uint nOfIntersections = 0; //used as orientation misses
            
    PrioNodeHandle follower_h; //if a PPrioNode is chosen to subdivide one way the follower is chosen to subdivide the other in the DP algorithm
    bool followerValid = false;
    std::list<PrioNodeHandle> guides; // this object would follow all PPrioNodes in the list
    
    PrioNode2(const NodeID node, std::list<Intervall2>::iterator intervallIt,
                std::vector<NodeInBox> &leftBox, std::vector<NodeInBox> &rightBox
                //, std::list<PrioNode2>::iterator follower, bool followerValid,
                //std::list<std::list<PrioNode2>::iterator> guides
                    ): node_id(node), intervallIt(intervallIt), leftBox(leftBox), rightBox(rightBox), perpendicularLength(0)
                    //, follower(follower), followerValid(followerValid), guides(guides)
    {}

     explicit operator simplePrioNode() const {
        return simplePrioNode(node_id, perpendicularLength, nOfIntersections);
    } 
    
    PrioNode2 operator =(const PrioNode2 &p_IN) {
        PrioNode2 p_out(p_IN.node_id, p_IN.intervallIt, p_IN.leftBox, p_IN.rightBox);
        p_out.posInIntervallIt = p_IN.posInIntervallIt;
        p_out.perpendicularLength = p_IN.perpendicularLength;
        p_out.nOfIntersections = p_IN.nOfIntersections;
        p_out.follower_h = p_IN.follower_h;
        p_out.followerValid = p_IN.followerValid;
        p_out.guides = p_IN.guides;
        
        return p_out;        
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

struct MatchNode {
    std::list<MatchNode>::iterator follower; //if a PPrioNode is chosen to subdivide one way the follower is chosen to subdivide the other in the DP algorithm
    bool followerValid = false;
    std::list<std::list<MatchNode>::iterator> guides; // this object would follow all PPrioNodes in the list
    
    NodeID node_id = c::NO_NID;
    
    MatchNode(NodeID node_id): node_id(node_id) {}
};

struct MatchPair {
    std::list<MatchNode> matchNodes1;
    std::list<MatchNode> matchNodes2;
    
    //MatchPair(Chain chain1, Chain chain2) {     
    MatchPair(std::list<NodeID> &chain1, std::list<NodeID> &chain2) {                
        for (chc::NodeID node_id: chain1) {            
            matchNodes1.push_back(MatchNode(node_id));
        }
        for (chc::NodeID node_id: chain2) {            
            matchNodes2.push_front(MatchNode(node_id));
        }
    }
};

//data types used for ChainPair version of DouglasPeucker ("P" = "Pair")--------------------------------------------------------
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