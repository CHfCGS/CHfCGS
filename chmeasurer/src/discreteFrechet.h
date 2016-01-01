#pragma once

#include "chgraph.h"
#include "nodes_and_edges.h"
#include "chains.h"
#include "limits"

#include <math.h>
#include <vector>



namespace chm {
    
    namespace unit_tests
    {
	void testdiscreteFrechet();        
    }
    

    class DiscreteFrechet {
        const CHGraph<CHNode, CHEdge> &graph;
        typedef std::vector<CHNode> ChainNodeVec;
        ChainNodeVec p, q;
        std::vector<std::vector<double> > ca;
        
        ChainNodeVec ChainToChainNodeVecChain(const Chain &chain, bool sameDirection) const {
            ChainNodeVec chain_node_vec;            
            if (sameDirection) {
                for (Chain::const_iterator it = chain.begin(); it != chain.end(); it++) {                      
                    chain_node_vec.push_back(graph.getNode(*it));                    
                }
            } else {
                for (Chain::const_reverse_iterator it = chain.rbegin(); it != chain.rend(); it++) {                    
                    chain_node_vec.push_back(graph.getNode(*it));                       
                }
            }
            return chain_node_vec;
        }
        
        double c(uint i, uint j) { 
            
            assert(i < p.size() && j < q.size());
            
            if (ca[i][j] > -1) {
                return ca[i][j];
            } else if (i == 0 and j == 0) {
                ca[i][j] = geo::geoDist(p[0], q[0]);
            } else if (i > 0 and j == 0) {                            
                ca[i][j] = std::max(c(i - 1, 0), geo::geoDist(p[i], q[0])); 
            } else if (i == 0 and j > 0) {
                ca[i][j] = std::max(c(0, j - 1), geo::geoDist(p[0], q[j]));
            } else if (i > 0 and j > 0) {
                ca[i][j] = std::max(
                         std::min(c(i - 1, j), std::min(c(i - 1, j - 1), c(i, j - 1))),
                         geo::geoDist(p[i], q[j])
                         );
            } else {
                ca[i][j] = std::numeric_limits<double>::max();
            }
            
            return ca[i][j];
        }
        
        

    public:
        
        
        DiscreteFrechet(const CHGraph<CHNode, CHEdge> &graph) : graph(graph) {            
        }
        
        double calc_dF (const Chain &chain1, const Chain &chain2) {
            p = ChainToChainNodeVecChain(chain1, true);
            q = ChainToChainNodeVecChain(chain2, false);
            std::vector<double> subvector(q.size(), -1.0);                        
            ca.assign(p.size(), subvector);                        
            return c(p.size()-1, q.size()-1);
        }                
    };
    

}