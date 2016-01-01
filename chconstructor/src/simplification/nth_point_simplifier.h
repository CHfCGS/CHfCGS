#pragma once

#include "../nodes_and_edges.h"
#include "prio_nodes.h"
#include "lineSimplifier.h"

namespace ls {

template <class GraphT>
class NthPointSimplifier : public LineSimplifier{
    

public:
    NthPointSimplifier() {}
    ~NthPointSimplifier() {}
    std::list<simplePrioNode> process(const Chain &chain1, const Chain &chain2) { 
        assert(chain2.empty());
        assert(chain1.size() >= 3);
        Chain chain_copy = chain1;
        chain_copy.pop_front();
        chain_copy.pop_back();
        std::list<simplePrioNode> priolist;
        while (!chain_copy.empty()) {
            uint mod_counter = 0;
        
            for (auto it = chain_copy.begin(); it != chain_copy.end();) {
                if (mod_counter%2 == 0) {
                    priolist.push_back(simplePrioNode(*it, 0, 0));
                    it = chain_copy.erase(it);                    
                } else {
                    it++;
                }
                mod_counter++;                               
            }
        }
        assert(chain1.size() == priolist.size()+2);
        return priolist;
    }
    
};

}