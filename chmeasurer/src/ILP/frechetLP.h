#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "../nodes_and_edges.h"
#include "linearProgram.h"
#include "frechet_test_data.h"

class FrechetLP : public LinearProgram
{
protected:
    virtual void setAllDegreeCoefficients(const FrechetTest_data& ilp_data) = 0;

    void addOutDegreeRows(ILP_Chain ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        for (auto it = ilp_chain.begin(); it != ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 2);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0.0);
        }
    }

    virtual void addAllDegreeRows(FrechetTest_data &ilp_data) = 0;

    void addLinkUnorderingRows(const std::vector<CrossLinkUnordering> &unorderings)
    {
        for (const CrossLinkUnordering &unordering : unorderings)
        {
            glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_DB, 0.0, 1.0);

            //then the two entries corresponding to the intersecting lines are set

            const CrossLink &link1 = unordering.link1;
            setIntersectionEntry(link1.id);

            const CrossLink &link2 = unordering.link2;
            setIntersectionEntry(link2.id);

            std::stringstream ss("");
            ss << "CI: (" << link1.src.ch_node_id << "->" << link1.line.start.ch_node_id << "," << link1.line.end.ch_node_id << "),"
                    "(" << link2.src.ch_node_id << "->" << link2.line.start.ch_node_id << "," << link2.line.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *row_name = s.c_str();
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name);

            //assert(line1.start.ch_node_id != line2.end.ch_node_id);
        }
    }

    virtual void addAllIntersectionRows(const FrechetTest_data &ilp_data) = 0;

    void addLinkColumns(const std::vector<CrossLink> &links)
    {
        debug_assert(links.size() > 0);
        glp_add_cols(lp, links.size());
        for (uint i = 0; i < links.size(); i++)
        {
            const CrossLink &link = links.at(i);
            uint col_id = link.id + 1;
            std::string colType = "cross link";
            std::stringstream ss("");
            ss << colType << ": (" << link.src.ch_node_id << "->" << link.line.start.ch_node_id << "," << link.line.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *col_name = s.c_str();
            glp_set_col_name(lp, col_id, col_name);
            glp_set_col_kind(lp, col_id, GLP_BV);
            glp_set_obj_coef(lp, col_id, 0);
        }
    }

    void addColumns(const std::vector<Line> &lines, double obj_coef)
    {
        debug_assert(lines.size() > 0);
        glp_add_cols(lp, lines.size());
        for (uint i = 0; i < lines.size(); i++)
        {
            const Line &line = lines.at(i);
            uint col_id = line.id + 1;
            std::string colType = "line"; // obj_coef == 1 ? "Edge" : "followerLine";
            std::stringstream ss("");
            ss << colType << ": (" << line.start.ch_node_id << "," << line.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *col_name = s.c_str();
            glp_set_col_name(lp, col_id, col_name);
            glp_set_col_kind(lp, col_id, GLP_BV);
            glp_set_obj_coef(lp, col_id, obj_coef); //objective is number of used edges
        }
    }

    virtual void setColumns(const FrechetTest_data &ilp_data) = 0;

    FrechetLP(const CHGraph<CHNode, CHEdge> &graph) : LinearProgram(graph) { }

    virtual ~FrechetLP() { }
};