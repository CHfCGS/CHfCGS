#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "nodes_and_edges.h"
#include "ILP_data.h"

class lineSimplificationILP {
    
    CHGraph<CHNode, CHEdge> &graph;
    glp_prob *lp;
    
    const static size_t size = 40000;
    
    
    //numbers of needed sizes beforehand 
    int preNofRows;
    int preNofCols;
    int needed_size;
    
    //whole ILP is invalid if size is not big enough
    bool enoughSpace;
    
    int ia[1+size], ja[1+size];
    //std::vector<double> ar;
    double ar[1+size];        
    
    size_t nofNonZeros = 0;
    
    
    double objective_value = 0;
    
    void fillArrays() {
        
        for (int i = 0; i < 1 + needed_size; i++) {
            ia[i] = 0;
            ja[i] = 0;
            ar[i] = 0;
        }
        return;
    }
    
    bool isInRange(size_t index) {
        return (1 <= index && index <= size); //0 is not valid for ILP
    }
    
    void setInOutDegreeCoefficients(std::vector<Line> &lines, NodeID node_id) {   
                
        for (Line &line : lines) { 
            nofNonZeros++;
            size_t arrayIndex = nofNonZeros;
            assert(isInRange(arrayIndex));
            ia[arrayIndex] = glp_get_num_rows(lp);
            ja[arrayIndex] = line.id + 1;
            assert(line.start.ch_node_id != line.end.ch_node_id);                                 
            
            if (line.start.ch_node_id == node_id) {
                ar[arrayIndex] = -1.0;                
            } else if(line.end.ch_node_id == node_id) {
                ar[arrayIndex] = 1.0;                
            } else {
                //ar[arrayIndex] = 0;
            }
        }
    }
    
    void addDegreeRows(ILP_Chain ilp_chain, std::vector<Line> &lines) {        
        //nodes 0 and n-1 should have degree 1
                
        
        assert(ilp_chain.size() >1);
        //nodes 0 and n-1 should have degree 1        
        glp_add_rows(lp, 1);            
        glp_set_row_name(lp, glp_get_num_rows(lp), "node: 0");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, -1.0, 0);
        setInOutDegreeCoefficients(lines, ilp_chain.front().ch_node_id);         
        
        //all other nodes should have indegree = outdegree
        //for (uint node_id = 1; node_id < graph.nodes.size()-1; node_id++) {                    
        for (auto it = ++ilp_chain.begin(); it != --ilp_chain.end(); it++) {              
            glp_add_rows(lp, 1);
            
            NodeID node_id = it->ch_node_id;
            std::stringstream ss("");            
            ss << "node: " << node_id;            
            std::string s = ss.str();
            char const *row_name = s.c_str();
                        
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name);            
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);                                                 
            setInOutDegreeCoefficients(lines, node_id);            
        }
        
        
        glp_add_rows(lp, 1);
        glp_set_row_name(lp, glp_get_num_rows(lp), "node: n-1");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);        
        setInOutDegreeCoefficients(lines, ilp_chain.back().ch_node_id); 
        
        
    }
    
    void addIntersectionRows(std::vector<Line> &lines, std::vector<Intersection> &intersections) {
        for (Intersection &intersection : intersections ) {            
            //glp_set_row_name(lp, nofRows + i, "("<< intersection.line.start << ", " << intersection.line.end << ")");
            glp_add_rows(lp, 1);                              
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_DB, 0.0, 1.0);
            /*
            //first zeroing out the whole row
            for (Line &line : lines) {        
                size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
                if(!isInRange(arrayIndex)) {
                    int i = glp_get_num_cols(lp);
                    auto j = glp_get_num_rows(lp)-1;
                    assert(isInRange(arrayIndex));
                }
                assert(isInRange(arrayIndex));
                
                ia[arrayIndex] = glp_get_num_rows(lp);
                ja[arrayIndex] = line.id + 1;
                assert(line.start.ch_node_id != line.end.ch_node_id);
                ar[arrayIndex] = 0;
            }
            */
            
            //then the two entries corresponding to the intersecting lines are set
            Line &line1 = intersection.line1;
            nofNonZeros++;
            size_t Line1Index = nofNonZeros;
            assert(isInRange(Line1Index));
            ia[Line1Index] = glp_get_num_rows(lp);
            ja[Line1Index] = line1.id + 1;
            ar[Line1Index] = 1.0;                       
            
            Line &line2 = intersection.line2;
            nofNonZeros++;
            size_t Line2Index = nofNonZeros;
            assert(isInRange(Line2Index));
            ia[Line2Index] = glp_get_num_rows(lp);
            ja[Line2Index] = line2.id + 1;
            ar[Line2Index] = 1.0;              
            
            std::stringstream ss("");            
            ss << "I: (" << line1.start.ch_node_id << "," << line1.end.ch_node_id << "),(" << line2.start.ch_node_id << "," << line2.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *row_name = s.c_str();
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name); 
            
            assert(line1.start.ch_node_id != line2.end.ch_node_id);
        }
    }
    
    //number of Columns is equal to the number of edges
    void setColumns(std::vector<Line> &lines) {
        debug_assert(lines.size() > 0);        

        glp_add_cols(lp, lines.size());        
        for (uint i = 0; i < lines.size(); i++) {
            Line &line = lines.at(i);
            std::stringstream ss("");            
            ss << "Line: (" << line.start.ch_node_id << "," << line.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *col_name = s.c_str();
            glp_set_col_name(lp, i+1, col_name);
            //glp_set_col_name(lp, i, (line.id << " "<< line.start << ", " << line.end << ")");
            glp_set_col_kind(lp, i+1, GLP_BV);
            //glp_set_obj_coef(lp, i+1, line.maxError);
            glp_set_obj_coef(lp, i+1, 1); //objective is number of used edges
        }
    }
    /*
    void addNofLinesRow(std::vector<Node> &nodes, std::vector<Line> &lines, LineID constraint) {
        glp_add_rows(lp, 1);
        //glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_DB, 0.0, constraint);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, constraint, 0.0);
        for (Line &line: lines) {
            size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
            assert(isInRange(arrayIndex));
            ia[arrayIndex] = glp_get_num_rows(lp);
            ja[arrayIndex] = line.id + 1;
            assert(line.start.ch_node_id != line.end.ch_node_id);
            ar[arrayIndex] = 1.0;
            nofNonZeros++;
        }

    }
     * */
    void prepareArray(ILP_data &ilp_data) {
        assert(ilp_data.ilp_chain1.size() > 1);
        
        preNofCols = ilp_data.potEdges1.size();        
        preNofRows = (ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size()) //in out degree equation for each node                
                + ilp_data.edgeIntersections.size();               
        //needed_size = preNofCols *preNofRows;
        nofNonZeros = 0;
        enoughSpace = true;// needed_size < (int) size;        
        /*
        if (enoughSpace) {
            fillArrays();            
        }
         * */
    }
    Chain::iterator calculateHalfListIt(Chain &chain) {
        int i = 0;
        int size = chain.size();
        Chain::iterator halfListIt = chain.begin();        
        while (size > 2 * i) {
            halfListIt++;
            i++;
        }
        return halfListIt;
    }
    
public:

    lineSimplificationILP(CHGraph<CHNode, CHEdge> &graph): graph(graph) {    
        
    }
    
    ~lineSimplificationILP() {
        
    }
    
    void solve(Chain chain, double epsilon) {
        
        
        ILP_data ilp_data = ILP_data(graph, chain, Chain(), epsilon, 0, false); 
        prepareArray(ilp_data); //==nofRows calculated beforehand       
        
        
        lp = glp_create_prob();
        glp_term_out(GLP_OFF); //make glp not verbose

        glp_set_prob_name(lp, "lineGeneralization");
        glp_set_obj_dir(lp, GLP_MIN);

        setColumns(ilp_data.potEdges1);

        addDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1);
        //addNofLinesRow(nodes, lines, 9); //useful when optimizing for epsilon
        addIntersectionRows(ilp_data.potEdges1, ilp_data.edgeIntersections);

        assert(preNofRows == glp_get_num_rows(lp));

        //auto numCols = glp_get_num_cols(lp)
        //nofNonZeros scheint doch wenig mit ne zu tun zu haben, glp rechnet sich das selber aus, sieht aber doch so aus
        //glp_load_matrix(lp, glp_get_num_cols(lp) * glp_get_num_rows(lp), ia, ja, ar);
        glp_load_matrix(lp, nofNonZeros, ia, ja, ar);
        //glp_write_lp(lp, NULL, "lp.txt");


        //solve ILP
        glp_simplex(lp, NULL);

        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        glp_intopt(lp, &parm);



        //glp_print_mip(lp, "solution.txt");
        //glp_intopt(lp, NULL);
        //
        // recover and display results 
        objective_value = glp_get_obj_val(lp);
        //std::cout << "Needed number of Edges: " << objective_value;
        //x1 = glp_get_col_prim(lp, 1);
        //x2 = glp_get_col_prim(lp, 2);

        /*
        printf("objective_value = %g\n", objective_value);
        for (int i = 1; i < glp_get_num_cols(lp)+1; i++) {
            printf("line %d: %g \n", i, glp_get_col_prim(lp, i));
        }
         * */
        glp_delete_prob(lp);
        glp_free_env();

        /*
        //cut Chain in half            
        auto halfListIt = calculateHalfListIt(chain);
        Chain firstHalfChain;            
        firstHalfChain.splice(firstHalfChain.begin(), chain, halfListIt);
        solve(firstHalfChain, epsilon);
        solve(chain, epsilon);
 
        Print("too big");
         */       
        
    }
    
};