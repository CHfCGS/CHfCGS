#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "nodes_and_edges.h"
#include "ILP_data.h"

class ParallelLineSimplificationILP {
    
    CHGraph<CHNode, CHEdge> &graph;
    glp_prob *lp;
    
    const static size_t size = 40000;
    size_t nofNonZeros = 0;
    
    //numbers of needed sizes beforehand 
    int preNofRows;
    int preNofCols;
    int preNofNonZeros = 0;
    
    //whole ILP is invalid if size is not big enough
    bool enoughSpace;
    
    int ia[1+size], ja[1+size];
    double ar[1+size];        
    
    
    double objective_value = 0;    
    
    bool isInRange(size_t index) {
        return (1 <= index && index <= size); //0 is not valid for ILP
    }
    
    void setRowName(NodeID node_id) {        
        std::stringstream ss("");
        ss << "node: " << node_id;
        std::string s = ss.str();
        char const *row_name = s.c_str();     
        glp_set_row_name(lp, glp_get_num_rows(lp), row_name);
    }
    
    void setDegreeCoefficient(uint i, uint j, double r) {        
        nofNonZeros++;
        size_t index = nofNonZeros;
        assert(isInRange(index));
        ia[index] = i;
        ja[index] = j;
        ar[index] = r;
    }
        
    /*
    void setDegreeCoefficients(std::vector<Line> &lines, size_t srcOffset1, size_t tgtOffset1, size_t srcOffset2, size_t tgtOffset2) {
        for (Line &line : ilp_data.potEdges1) { 
            assert(line.start.ch_node_id != line.end.ch_node_id);                                 
            setDegreeCoefficient(srcOffset1 + line.start.posInChain + 1, line.id + 1 , -1.0);           
            setDegreeCoefficient(tgtOffset1 + line.end.posInChain + 1, line.id + 1 , 1.0);
            setDegreeCoefficient(srcOffset1 + line.start.posInChain + 1, line.id + 1 , -1.0);           
            setDegreeCoefficient(tgtOffset1 + line.end.posInChain + 1, line.id + 1 , 1.0);
        }
    }
     * */
    
    void setAllDegreeCoefficients(ILP_data ilp_data) {
        size_t nofNodes = ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size();
        size_t chainPosOffset = ilp_data.ilp_chain1.size();
        
        for (Line &line : ilp_data.potEdges1) { 
            assert(line.start.ch_node_id != line.end.ch_node_id);                                 
            setDegreeCoefficient(line.start.posInChain + 1, line.id + 1 , -1.0);           
            setDegreeCoefficient(line.end.posInChain + 1, line.id + 1 , 1.0);                                 
            setDegreeCoefficient(nofNodes + line.start.posInChain + 1, line.id + 1, -1.0); //out
            setDegreeCoefficient(2*nofNodes + line.end.posInChain + 1, line.id + 1, 1.0); //in
        }
                
        for (Line &line : ilp_data.potEdges2) { 
            assert(line.start.ch_node_id != line.end.ch_node_id);                                 
            setDegreeCoefficient(chainPosOffset + line.start.posInChain + 1, line.id + 1 , -1.0);           
            setDegreeCoefficient(chainPosOffset + line.end.posInChain + 1, line.id + 1 , 1.0);                                 
            setDegreeCoefficient(nofNodes + chainPosOffset + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(2*nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, 1.0);
        }
                
        for (Line &line : ilp_data.followerLines1) { 
            assert(line.start.ch_node_id != line.end.ch_node_id);                                 
            //out
            setDegreeCoefficient(nofNodes + line.start.posInChain + 1, line.id + 1 , 1.0);           
            //setDegreeCoefficient(nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1 , 1.0);                                 
            //in
            setDegreeCoefficient(2*nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, -1.0);
            //setDegreeCoefficient(2*nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, -1.0);
        }                
        
        for (Line &line : ilp_data.followerLines2) { 
            assert(line.start.ch_node_id != line.end.ch_node_id);                                 
            //out
            setDegreeCoefficient(nofNodes + chainPosOffset + line.start.posInChain + 1, line.id + 1 , 1.0);           
            //setDegreeCoefficient(nofNodes + line.end.posInChain + 1, line.id + 1 , 1.0);                                 
            //in
            setDegreeCoefficient(2*nofNodes + line.end.posInChain + 1, line.id + 1, -1.0);
            //setDegreeCoefficient(2*nofNodes + line.end.posInChain + 1, line.id + 1, -1.0);
        }   
        
    }    
    void addInOutDegreeRows(ILP_Chain ilp_chain) {        
        assert(ilp_chain.size() >1);
        //nodes 0 and n-1 should have degree 1        
        glp_add_rows(lp, 1);            
        setRowName(ilp_chain.front().ch_node_id);
        //glp_set_row_name(lp, glp_get_num_rows(lp), "node: 0");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, -1.0, 0);
        //setInOutDegreeCoefficients(lines, ilp_chain.front().ch_node_id);         
        
        //all other nodes should have indegree = outdegree
        //for (uint node_id = 1; node_id < graph.nodes.size()-1; node_id++) {                    
        for (auto it = ++ilp_chain.begin(); it != --ilp_chain.end(); it++) {              
            glp_add_rows(lp, 1);            
            setRowName(it->ch_node_id);
            //glp_set_row_name(lp, glp_get_num_rows(lp), row_name);            
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);                                                 
            //setInOutDegreeCoefficients(lines, node_id);            
        }
        
        
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id);        
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);        
        //setInOutDegreeCoefficients(lines, ilp_chain.back().ch_node_id);         
    }
    
       
    //outdegrees has too be equal except for ends
    void addOutDegreeRows(ILP_Chain ilp_chain) { 
        assert(ilp_chain.size() >1);
        //all other nodes should have edge_outdegree = follower_outdegree
        //for (uint node_id = 0; node_id < graph.nodes.size()-1; node_id++) {   
        for (auto it = ilp_chain.begin(); it != --ilp_chain.end(); it++) {       
            glp_add_rows(lp, 1);    
            setRowName(it->ch_node_id);          
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);                                                 
            //setOutDegreeCoefficients(potEdges, followerLines, node_id);            
        }
        
        //except the last one
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id); 
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);        
        //setOutDegreeCoefficients(potEdges, followerLines, ilp_chain.back().ch_node_id); 
    }        
        
    void addInDegreeRows(ILP_Chain ilp_chain) {   
        assert(ilp_chain.size() >1);
        //first node is always in the simplfication so follower_indegree is unconstrained        
        glp_add_rows(lp, 1);            
        setRowName(ilp_chain.front().ch_node_id); 
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FR, 0, 0);
        //glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, -1000, 0);
        //setInOutDegreeCoefficients(potEdges, followerLines, 0);         
                
        //all other nodes should have const * edge_indegree > follower_indegree
        //=> that means: if a node has at least edge_indegree, follower_indgree is practically unconstrained
        //for (uint node_id = 1; node_id < graph.nodes.size(); node_id++) {            
        for (auto it = ++ilp_chain.begin(); it != ilp_chain.end(); it++) {   
            glp_add_rows(lp, 1); 
            setRowName(it->ch_node_id);      
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, 0.0, 0.0);                                                 
            //setInDegreeCoefficients(potEdges, followerLines, node_id);            
        }                         
    }

    void addAllDegreeRows(ILP_data &ilp_data) {
        addInOutDegreeRows(ilp_data.ilp_chain1);
        addInOutDegreeRows(ilp_data.ilp_chain2);
        addOutDegreeRows(ilp_data.ilp_chain1);
        addOutDegreeRows(ilp_data.ilp_chain2);
        addInDegreeRows(ilp_data.ilp_chain1);
        addInDegreeRows(ilp_data.ilp_chain2);
    }    
    
    void setIntersectionEntry(LineID line_id) {
            nofNonZeros++;
            size_t index = nofNonZeros;            
            ia[index] = glp_get_num_rows(lp);
            ja[index] = line_id + 1;
            ar[index] = 1.0;
    }   
    
    void addIntersectionRows(const std::vector<Intersection> &intersections) {
        for (const Intersection &intersection : intersections ) {                        
            glp_add_rows(lp, 1);                              
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_DB, 0.0, 1.0);            
            
            //then the two entries corresponding to the intersecting lines are set
            
            const Line &line1 = intersection.line1;            
            setIntersectionEntry(line1.id);           
            
            const Line &line2 = intersection.line2;            
            setIntersectionEntry(line2.id);              
            
            std::stringstream ss("");            
            ss << "I: (" << line1.start.ch_node_id << "," << line1.end.ch_node_id << "),(" << line2.start.ch_node_id << "," << line2.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *row_name = s.c_str();
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name); 
            
            assert(line1.start.ch_node_id != line2.end.ch_node_id);
        }
    }
    
    void addAllIntersectionRows(ILP_data &ilp_data) {        

        addIntersectionRows(ilp_data.edgeIntersections);
        addIntersectionRows(ilp_data.followerLinesUnorderings);
    }  
              
        
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
            glp_set_col_kind(lp, col_id, GLP_BV);            
            glp_set_obj_coef(lp, col_id, obj_coef); //objective is number of used edges
        }
    }
    
    void setAllColumns(ILP_data &ilp_data) {
        addColumns(ilp_data.potEdges1, 1.0);
        addColumns(ilp_data.potEdges2, 1.0);
        addColumns(ilp_data.followerLines1, 0);
        addColumns(ilp_data.followerLines2, 0);        
    }
        
    
    void prepareArray(const ILP_data &ilp_data) {
        assert(ilp_data.ilp_chain1.size() > 1);
        assert(ilp_data.ilp_chain2.size() > 1);
        
        preNofCols = ilp_data.potEdges1.size() + ilp_data.potEdges2.size()
                + ilp_data.followerLines1.size() + ilp_data.followerLines2.size();        
        preNofRows = 3 * (ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size()) //3 degree equations for each node                
                + ilp_data.edgeIntersections.size()
                + ilp_data.followerLinesUnorderings.size();
        preNofNonZeros = 4*ilp_data.potEdges1.size() + 4*ilp_data.potEdges2.size()
                + 2*ilp_data.followerLines1.size() + 2*ilp_data.followerLines2.size()
                + 2*ilp_data.edgeIntersections.size()
                + 2*ilp_data.followerLinesUnorderings.size(); 
        nofNonZeros = 0;
        enoughSpace =  preNofNonZeros < (int) size;                     
    }
    
   
public:
    
    ParallelLineSimplificationILP(CHGraph<CHNode, CHEdge> &graph): graph(graph) {
        
    }
    
    ~ParallelLineSimplificationILP() {
        
    }
    
    double solve(Chain chain1, Chain chain2, double epsilon, double eta) {
        
        ILP_data ilp_data = ILP_data(graph, chain1, chain2, epsilon, eta, true); 
        prepareArray(ilp_data);
        lp = glp_create_prob();
        glp_term_out(GLP_OFF); //make glp not verbose
        
        glp_set_prob_name(lp, "parallellineGeneralization");
        
        glp_set_obj_dir(lp, GLP_MIN); //minimize number of used edges
        
        setAllColumns(ilp_data);        
        
        addAllDegreeRows(ilp_data);
        setAllDegreeCoefficients(ilp_data);
        addAllIntersectionRows(ilp_data);
                
        //addIntersectionRows(ilp_data)
        
        /*
        addInOutDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1);
        addInOutDegreeRows(ilp_data.ilp_chain2, ilp_data.potEdges2);
        
        addOutDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1, ilp_data.followerLines1);
        addOutDegreeRows(ilp_data.ilp_chain2, ilp_data.potEdges2, ilp_data.followerLines2);
        
        addInDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1, ilp_data.followerLines2);
        addInDegreeRows(ilp_data.ilp_chain2, ilp_data.potEdges2, ilp_data.followerLines1);
        
        
        addIntersectionRows(ilp_data.edgeIntersections);
        addIntersectionRows(ilp_data.followerLinesUnorderings);                
        */
        
        assert(preNofRows == glp_get_num_rows(lp));
        assert(preNofNonZeros == (int) nofNonZeros);
        //auto numCols = glp_get_num_cols(lp)
        //assert(nofEquations == glp_get_num_rows(lp));
        
        
        glp_load_matrix(lp, nofNonZeros, ia, ja, ar);         
        glp_write_lp(lp, NULL, "lp.txt");
        
        
        
        
        /* solve problem */
        
        glp_simplex(lp, NULL);
        
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        glp_intopt(lp, &parm);
        
        
        
        glp_print_mip(lp, "assdfAS.txt");
        
        /* recover and display results */
        objective_value = glp_mip_obj_val(lp);        
        /*
        printf("objective_value = %g\n", objective_value);
        for (int i = 1; i < glp_get_num_cols(lp)+1; i++) {
            printf("line %d: %g \n", i, glp_get_col_prim(lp, i));
        }
         * */
        
        glp_delete_prob(lp);
        glp_free_env();
        
        return objective_value;

    }
};