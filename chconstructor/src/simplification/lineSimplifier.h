#pragma once

#include "../defs.h"
#include "../nodes_and_edges.h"
#include "../chains.h"

#include <vector>

namespace chc {

    /*
     * Interface of a LineSimplifier used in own Prioritizer.
     */
    class LineSimplifier {
    public:
        /*
         * processes a chain and 
         */
        virtual std::list<ls::simplePrioNode> process(const Chain &chain1, const Chain &chain2) = 0;
        

        //virtual destructor                
        virtual ~LineSimplifier() {            
        }
    };

}