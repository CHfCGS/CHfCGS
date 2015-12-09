#include "defs.h"
//#include "track_time.h"
#include <getopt.h>
#include "nodes_and_edges.h"
#include "ch_parser.h"
#include "chgraph.h"
#include "ch_measurer.h"
#include "measure_options.h"



using namespace chm;
//using namespace std::chrono;


void printHelp() {
    std::cout
            << "Usage: ./ch_measurer [ARGUMENTS]\n"
            << "Mandatory arguments are:\n"
            << "  -i, --infile <path>        Read graph from <path>\n"
            << "Optional arguments are:\n"
            << "  -d, --dijkstra <type>"
            << "  -l, --ilp <type>"
            << "  -p, --P_ilp <type>";
}


int main(int argc, char** argv) {

    std::string filepath;
    MeasureOptions m_options;
    
    const struct option longopts[] = {        
        {"infile", required_argument, 0, 'i'}, 
        {"dijkstra", no_argument, 0, 'd'},   
        {"ilp", no_argument, 0, 'l'},   
        {"p_ilp", no_argument, 0, 'p'},   
        {0, 0, 0, 0},
    };

    int index(0);
    int iarg(0);
    opterr = 1;

    while ((iarg = getopt_long(argc, argv, "i:dlp", longopts, &index)) != -1) {
        switch (iarg) {            
            case 'i':
                filepath = optarg;
                break;                                    
            case 'd':
                m_options.dijkstra = true;
                break;
            case 'l':
                m_options.ilp = true;
                break;
            case 'p':
                m_options.p_ilp = true;
                break;            
            default:
                printHelp();
                return 1;
                break;
        }
    }
    /*
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
     * */
    GraphInData<CHNode, CHEdge> graphInData;
    
    CH_Parser ch_Parser = CH_Parser();    
    
    
    if (!ch_Parser.parseTxtGraphFile(filepath, graphInData.nodes, graphInData.edges)) {
        std::cout << "Couln't open file " << filepath << std::endl;
        return 1;
    }
        
    CHGraph<CHNode, CHEdge> g;
    g.init(std::move(graphInData));
    
    CHMeasurer ch_measurer = CHMeasurer(g);
    ch_measurer.makeMeasurement(m_options);
    
    
    
    
return 0;
}
