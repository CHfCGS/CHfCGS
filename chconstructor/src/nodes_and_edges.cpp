#include "nodes_and_edges.h"

namespace chc
{
	Edge concat(Edge const& edge1, Edge const& edge2)
	{
		assert(edge1.tgt == edge2.src);
		return Edge(constant::NO_EID, edge1.src, edge2.tgt, edge1.dist + edge2.dist);
	}

	OSMEdge concat(OSMEdge const& edge1, OSMEdge const& edge2)
	{
		assert(edge1.tgt == edge2.src);
                uint type=0;
                /*
                if (edge1.type == edge2.type) {
                    type = edge1.type;
                }*/                
                type = std::max(edge1.type, edge2.type); 
		return OSMEdge(constant::NO_EID, edge1.src, edge2.tgt, edge1.dist + edge2.dist, type, -1);
	}

	EuclOSMEdge concat(EuclOSMEdge const& edge1, EuclOSMEdge const& edge2)
	{
		assert(edge1.tgt == edge2.src);
                uint type=0;
                if (edge1.type == edge2.type) {
                    type = edge1.type;
                }
		return EuclOSMEdge(constant::NO_EID, edge1.src, edge2.tgt, edge1.dist + edge2.dist, type, -1,
				edge1.eucl_dist + edge2.eucl_dist);
	}

	StefanEdge concat(StefanEdge const& edge1, StefanEdge const& edge2)
	{
		assert(edge1.tgt == edge2.src);
		return StefanEdge(constant::NO_EID, edge1.src, edge2.tgt, edge1.dist + edge2.dist);
	}
}
