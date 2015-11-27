#include "defs.h"
//#include "track_time.h"
#include <getopt.h>
#include "nodes_and_edges.h"
#include "ch_parser.h"
#include "chgraph.h"
#include "ch_measurer.h"



using namespace chm;
//using namespace std::chrono;





int main(int argc, char** argv) {

    
    
    /////////////////////////////////////
    // overly simple command line parsing
    /////////////////////////////////////

    std::string filepath;

    int i = 0;
    if (argc == 1) {
        std::cout << "Supply a graph with -gf <graph.gl>" << std::endl;
        return 0;
    }
    while (i < argc) {
        if (argv[i] == (std::string) "-gf") {
            i++;
            if (i < argc) {
                filepath = argv[i];
            } else {
                std::cout << "Missing parameter for -gf" << std::endl;
                return 0;
            }
        } else {
            i++;
        }
    }

    GraphInData<CHNode, CHEdge> graphInData;
    
    CH_Parser ch_Parser = CH_Parser();    
    
    
    if (!ch_Parser.parseTxtGraphFile(filepath, graphInData.nodes, graphInData.edges)) {
        std::cout << "Couln't open file " << filepath << std::endl;
        return 1;
    }
        
    CHGraph<CHNode, CHEdge> g;
    g.init(std::move(graphInData));
    
    CHMeasurer ch_measurer = CHMeasurer(g);
    ch_measurer.makeMeasurement();
    
    
    
    
return 0;
}
