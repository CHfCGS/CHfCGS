#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "nodes_and_edges.h"
#include "ILP_data.h"

class ParallelLineSimplificationILP {
    
    
  
    
    Graph &graph;
    glp_prob *lp;
    
    const static size_t size = 400000;
    
    int ia[1+size], ja[1+size];
    double ar[1+size];        
    
    
    double objective_value = 0;
    
    void fillArrays(uint nofCols, uint nofRows) {
        
        for (uint i = 1; i < 1 + nofRows; i++) {
            for (uint j = 1; j < 1 + nofCols; j++) {
                uint index = (i-1) * nofCols + j;
                assert(isInRange(index));
                ia[index] = i;
                ja[index] = j;
                ar[index] = 0;
            }
        }                
    }
    
    bool isInRange(size_t index) {
        return (1 <= index && index <= size); //0 is not valid for ILP
    }
    
    void setInOutDegreeCoefficients(std::vector<Line> &lines, NodeID node_id) {   
                
        for (Line &line : lines) {        
            size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
            assert(isInRange(arrayIndex));
            ia[arrayIndex] = glp_get_num_rows(lp);
            ja[arrayIndex] = line.id + 1;
            assert(line.start.ch_node_id != line.end.ch_node_id);                                    
            
            if (line.start.ch_node_id == node_id) {
                ar[arrayIndex] = -1.0;                
            } else if(line.end.ch_node_id == node_id) {
                ar[arrayIndex] = 1.0;                
            } else {
                ar[arrayIndex] = 0;
            }
        }
    }
    
    void addInOutDegreeRows(ILP_Chain ilp_chain, std::vector<Line> &lines) {        
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
    
    
    
    void setOutDegreeCoefficients(std::vector<Line> &potEdges, std::vector<Line> &followerLines, NodeID node_id) {                   
        for (Line &line : potEdges) {        
            size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
            assert(isInRange(arrayIndex));            
            assert(line.start.ch_node_id != line.end.ch_node_id);                                                
            if (line.start.ch_node_id == node_id) {
                ar[arrayIndex] = -1.0;                
            }
        }        
        for (Line &line : followerLines) {        
            size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
            assert(isInRange(arrayIndex));            
            assert(line.start.ch_node_id != line.end.ch_node_id);                                                
            if (line.start.ch_node_id == node_id) {
                ar[arrayIndex] = 1.0;                
            }
        }
    }
    
    //outdegrees has too be equal except for ends
    void addOutDegreeRows(ILP_Chain ilp_chain, std::vector<Line> &potEdges, std::vector<Line> &followerLines) {                
        //all other nodes should have edge_outdegree = follower_outdegree
        //for (uint node_id = 0; node_id < graph.nodes.size()-1; node_id++) {   
        for (auto it = ilp_chain.begin(); it != --ilp_chain.end(); it++) {       
            glp_add_rows(lp, 1);    
            NodeID node_id = it->ch_node_id;
            std::string s = "node: " + node_id;
            char const *row_name = s.c_str();
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name);            
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);                                                 
            setOutDegreeCoefficients(potEdges, followerLines, node_id);            
        }
        
        //except the last one
        glp_add_rows(lp, 1);
        glp_set_row_name(lp, glp_get_num_rows(lp), "node: n-1");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);        
        setOutDegreeCoefficients(potEdges, followerLines, ilp_chain.back().ch_node_id); 
    }
    
    void setInDegreeCoefficients(std::vector<Line> &potEdges, std::vector<Line> &followerLines, NodeID node_id) {                   
        for (Line &line : potEdges) {        
            size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
            assert(isInRange(arrayIndex));            
            assert(line.start.ch_node_id != line.end.ch_node_id);                                                
            if (line.end.ch_node_id == node_id) {
                ar[arrayIndex] = 1.0 * followerLines.size();                
            }
        }        
        for (Line &line : followerLines) {        
            size_t arrayIndex = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line.id + 1 ;
            assert(isInRange(arrayIndex));            
            assert(line.start.ch_node_id != line.end.ch_node_id);                                                
            if (line.end.ch_node_id == node_id) {
                ar[arrayIndex] = -1.0;                
            }
        }
    }
    
    //outdegrees has too be equal except for ends
    void addInDegreeRows(ILP_Chain ilp_chain, std::vector<Line> &potEdges, std::vector<Line> &followerLines) {   
        //first node is always in the simplfication so follower_indegree is unconstrained
        /*
        glp_add_rows(lp, 1);            
        glp_set_row_name(lp, glp_get_num_rows(lp), "node: 0");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0, 0);
        setInOutDegreeCoefficients(potEdges, followerLines, 0);         
         * */
        
        //all other nodes should have const * edge_indegree > follower_indegree
        //=> that means: if a node has at least edge_indegree, follower_indgree is practically unconstrained
        //for (uint node_id = 1; node_id < graph.nodes.size(); node_id++) {            
        for (auto it = ++ilp_chain.begin(); it != ilp_chain.end(); it++) {   
            glp_add_rows(lp, 1); 
            NodeID node_id = it->ch_node_id;
            std::string s = "node: " + node_id;
            char const *row_name = s.c_str();
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name);            
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, 0.0, 0.0);                                                 
            setInDegreeCoefficients(potEdges, followerLines, node_id);            
        }
        
        //except the last one
         
    }
    
    void addIntersectionRows(std::vector<Intersection> &intersections) {        
        
        for (Intersection &intersection : intersections ) {            
            //glp_set_row_name(lp, nofRows + i, "("<< intersection.line.start << ", " << intersection.line.end << ")");
            glp_add_rows(lp, 1);            
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_DB, 0.0, 1.0);
                                           
            Line &line1 = intersection.line1;
            size_t Line1Index = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line1.id + 1;
            assert(isInRange(Line1Index));
            //ia[Line1Index] = glp_get_num_rows(lp);
            //ja[Line1Index] = line1.id;
            ar[Line1Index] = 1.0;                   
            
            Line &line2 = intersection.line2;
            size_t Line2Index = glp_get_num_cols(lp)*(glp_get_num_rows(lp)-1) + line2.id + 1;
            assert(isInRange(Line2Index));
            //ia[Line2Index] = glp_get_num_rows(lp);
            //ja[Line2Index] = line2.id;
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
    void addColumns(std::vector<Line> &lines, double obj_coef) {
        debug_assert(lines.size() > 0);        
        glp_add_cols(lp, lines.size());        
        for (uint i = 0; i < lines.size(); i++) {            
            Line &line = lines.at(i);
            uint col_id = line.id + 1;
            std::string colType = obj_coef == 1 ? "Edge" : "followerLine";            
            std::stringstream ss("");           
            ss << colType << ": (" << line.start.ch_node_id << "," << line.end.ch_node_id << ")";            
            std::string s = ss.str();    
            char const *col_name = s.c_str();            
            glp_set_col_name(lp, col_id, col_name);
            //glp_set_col_name(lp, i, (line.id << " "<< line.start << ", " << line.end << ")");
            glp_set_col_kind(lp, col_id, GLP_BV);
            //glp_set_obj_coef(lp, i+1, line.maxError);
            glp_set_obj_coef(lp, col_id, obj_coef); //objective is number of used edges
        }
    }
    
    void setAllColumns(ILP_data &ilp_data) {
        addColumns(ilp_data.potEdges1, 1.0);
        addColumns(ilp_data.potEdges2, 1.0);
        addColumns(ilp_data.followerLines1, 0);
        addColumns(ilp_data.followerLines2, 0);        
    }
    
    uint calcNumberOfEquations(ILP_data &ilp_data) {
        assert(ilp_data.ilp_chain1.size() > 1 && ilp_data.ilp_chain2.size() > 1);
        uint nofEquations = 3 * (ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size()) //3 degree equations for each node
                - 2 //except for the last nodes which need no indegree-equation
                + ilp_data.edgeIntersections.size()
                + ilp_data.followerLinesUnorderings.size();
        return nofEquations;
    }
   
public:
    
    ParallelLineSimplificationILP(Graph &graph, Chain chain1, Chain chain2, double epsilon, double eta): graph(graph) {
        ILP_data ilp_data = ILP_data(graph, chain1, chain2, epsilon, eta, true); 
                
        lp = glp_create_prob();
        
        glp_set_prob_name(lp, "lineGeneralization");
        
        glp_set_obj_dir(lp, GLP_MIN); //minimize number of used edges
        
        setAllColumns(ilp_data);
        int nofEquations = calcNumberOfEquations(ilp_data); //==nofRows calculated beforehand
        fillArrays(glp_get_num_cols(lp), nofEquations);        
        
        addInOutDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1);
        addInOutDegreeRows(ilp_data.ilp_chain2, ilp_data.potEdges2);
        
        addOutDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1, ilp_data.followerLines1);
        addOutDegreeRows(ilp_data.ilp_chain2, ilp_data.potEdges2, ilp_data.followerLines2);
        
        addInDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1, ilp_data.followerLines2);
        addInDegreeRows(ilp_data.ilp_chain2, ilp_data.potEdges2, ilp_data.followerLines1);
        
        
        addIntersectionRows(ilp_data.edgeIntersections);
        addIntersectionRows(ilp_data.followerLinesUnorderings);                
        
        //auto numCols = glp_get_num_cols(lp)
        assert(nofEquations == glp_get_num_rows(lp));
        
        //nofNonZeros scheint doch wenig mit ne zu tun zu haben, glp rechnet sich das selber aus
        glp_load_matrix(lp, glp_get_num_cols(lp)*glp_get_num_rows(lp), ia, ja, ar);         
        glp_write_lp(lp, NULL, "lp.txt");
    }
    
    ~ParallelLineSimplificationILP() {
        glp_delete_prob(lp);
        glp_free_env();
    }
    
    void solve() {
        /* solve problem */
        
        glp_simplex(lp, NULL);
        
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        glp_intopt(lp, &parm);
        
        
        
        glp_print_mip(lp, "solution.txt");
        //glp_intopt(lp, NULL);
        //
        /* recover and display results */
        objective_value = glp_get_obj_val(lp);
        //x1 = glp_get_col_prim(lp, 1);
        //x2 = glp_get_col_prim(lp, 2);
        /*
        printf("objective_value = %g\n", objective_value);
        for (int i = 1; i < glp_get_num_cols(lp)+1; i++) {
            printf("line %d: %g \n", i, glp_get_col_prim(lp, i));
        }
         * */

    }
};