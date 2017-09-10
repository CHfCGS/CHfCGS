#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "../nodes_and_edges.h"
#include "linearProgram.h"
#include "ILP_data.h"

class lineSimplificationILP : LinearProgram
{

    void setInOutDegreeCoefficients(const std::vector<Line> &lines, NodeID node_id)
    {
        for (const Line &line : lines)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            //src
            setDegreeCoefficient(line.start.posInChain + 1, line.id + 1, -1.0);
            //tgt
            setDegreeCoefficient(line.end.posInChain + 1, line.id + 1, 1.0);
        }
    }

    void addDegreeRows(const ILP_Chain &ilp_chain, const std::vector<Line> &lines)
    {
        assert(ilp_chain.size() > 1);

        //nodes 0 and n-1 should have degree 1
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.front().ch_node_id, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, -1.0, 0);
        //setInOutDegreeCoefficients(lines, ilp_chain.front().ch_node_id);

        //all other nodes should have indegree = outdegree
        //for (uint node_id = 1; node_id < graph.nodes.size()-1; node_id++) {
        for (auto it = ++ilp_chain.begin(); it != --ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            /*
            NodeID node_id2 = it->ch_node_id;
            std::stringstream ss2("");
            ss2 << "node: " << node_id2;
            std::string s2 = ss2.str();
            char const *row_name2 = s2.c_str();
             * */
            setRowName(it->ch_node_id, 1);
            //glp_set_row_name(lp, glp_get_num_rows(lp), it->ch_node_id);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);
            //setInOutDegreeCoefficients(lines, it->ch_node_id);
        }
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);

        assert(glp_get_num_rows(lp) == (int) ilp_chain.size());
        setInOutDegreeCoefficients(lines, ilp_chain.back().ch_node_id);
    }


    //number of Columns is equal to the number of edges

    void setColumns(const ILP_data &ilp_data)
    {
        const std::vector<Line> &lines = ilp_data.potEdges1;
        debug_assert(lines.size() > 0);

        glp_add_cols(lp, lines.size());
        for (uint i = 0; i < lines.size(); i++)
        {
            const Line &line = lines.at(i);
            std::stringstream ss("");
            ss << "Line: (" << line.start.ch_node_id << "," << line.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *col_name = s.c_str();
            glp_set_col_name(lp, i + 1, col_name);
            glp_set_col_kind(lp, i + 1, GLP_BV);
            glp_set_obj_coef(lp, i + 1, 1); //objective is number of used edges
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
    void prepareArray(const ILP_data &ilp_data)
    {
        assert(ilp_data.ilp_chain1.size() > 1);

        preNofCols = ilp_data.potEdges1.size();
        preNofRows = (ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size()) //in out degree equation for each node
                + ilp_data.edgeIntersections.size();
        preNofNonZeros = 2 * ilp_data.potEdges1.size() + 2 * ilp_data.edgeIntersections.size();
        nofNonZeros = 0;
        enoughSpace = preNofNonZeros < (int) size;

        //Print("preNofNonZeros" << preNofNonZeros);
        assert(enoughSpace);
    }
    /*
    static ILP_Chain::iterator calculateHalfListIt(const ILP_Chain &chain) {
        int i = 0;
        int size = chain.size();
        ILP_Chain::iterator halfListIt = chain.begin();
        while (size > 2 * i) {
            halfListIt++;
            i++;
        }
        return halfListIt;
    }*/

public:

    lineSimplificationILP(const CHGraph<CHNode, CHEdge> &graph) : LinearProgram(graph) { }

    ~lineSimplificationILP() { }

    double simplify(const Chain &chain, double epsilon)
    {

        Chain emptyChain;
        ILP_data ilp_data = ILP_data(graph, chain, emptyChain, epsilon, 0, false);
        prepareArray(ilp_data);


        lp = glp_create_prob();
        glp_term_out(GLP_OFF); //make glp not verbose

        glp_set_prob_name(lp, "lineGeneralization");
        glp_set_obj_dir(lp, GLP_MIN);

        setColumns(ilp_data);

        addDegreeRows(ilp_data.ilp_chain1, ilp_data.potEdges1);
        //addNofLinesRow(nodes, lines, 9); //useful when optimizing for epsilon
        addIntersectionRows(ilp_data.edgeIntersections);

        assert(preNofRows == glp_get_num_rows(lp));
        assert(preNofNonZeros == (int) nofNonZeros);

        return solve();
    }

};