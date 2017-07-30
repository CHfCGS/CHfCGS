#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "../nodes_and_edges.h"
#include "calculation.h"
#include "ILP_data.h"

class ParallelLineSimplificationILP : Calculation
{

    void setDegreeCoefficient(uint i, uint j, double r)
    {
        nofNonZeros++;
        size_t index = nofNonZeros;
        assert(isInRange(index));
        ia[index] = i;
        ja[index] = j;
        ar[index] = r;
    }

    void setAllDegreeCoefficients(const ILP_data &ilp_data)
    {
        size_t nofNodes = ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size();
        size_t chainPosOffset = ilp_data.ilp_chain1.size();
        const double INDEGREEMULTIPLIER = 1000.0;

        for (const Line &line : ilp_data.potEdges1)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            setDegreeCoefficient(line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(line.end.posInChain + 1, line.id + 1, 1.0);
            setDegreeCoefficient(nofNodes + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(2 * nofNodes + line.end.posInChain + 1, line.id + 1, INDEGREEMULTIPLIER);
        }

        for (const Line &line : ilp_data.potEdges2)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            setDegreeCoefficient(chainPosOffset + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(chainPosOffset + line.end.posInChain + 1, line.id + 1, 1.0);
            setDegreeCoefficient(nofNodes + chainPosOffset + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(2 * nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, INDEGREEMULTIPLIER);
        }

        for (const Line &line : ilp_data.followerLines1)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id); //can happen by having the same center nodes
            //out
            setDegreeCoefficient(nofNodes + line.start.posInChain + 1, line.id + 1, 1.0);
            //setDegreeCoefficient(nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1 , 1.0);
            //in
            setDegreeCoefficient(2 * nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, -1.0);
            //setDegreeCoefficient(2*nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, -1.0);
        }

        for (const Line &line : ilp_data.followerLines2)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id); //can happen by having the same center nodes
            //out
            setDegreeCoefficient(nofNodes + chainPosOffset + line.start.posInChain + 1, line.id + 1, 1.0);
            //setDegreeCoefficient(nofNodes + line.end.posInChain + 1, line.id + 1 , 1.0);
            //in
            setDegreeCoefficient(2 * nofNodes + line.end.posInChain + 1, line.id + 1, -1.0);
            //setDegreeCoefficient(2*nofNodes + line.end.posInChain + 1, line.id + 1, -1.0);
        }

    }

    void addInOutDegreeRows(const ILP_Chain& ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        //nodes 0 and n-1 should have degree 1
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.front().ch_node_id, 1);
        //glp_set_row_name(lp, glp_get_num_rows(lp), "node: 0");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, -1.0, 0);
        //setInOutDegreeCoefficients(lines, ilp_chain.front().ch_node_id);

        //all other nodes should have indegree = outdegree
        //for (uint node_id = 1; node_id < graph.nodes.size()-1; node_id++) {
        for (auto it = ++ilp_chain.begin(); it != --ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 1);
            //glp_set_row_name(lp, glp_get_num_rows(lp), row_name);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);
            //setInOutDegreeCoefficients(lines, node_id);
        }


        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);
        //setInOutDegreeCoefficients(lines, ilp_chain.back().ch_node_id);
    }


    //outdegrees has too be equal except for ends

    void addOutDegreeRows(const ILP_Chain& ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        //all other nodes should have edge_outdegree = follower_outdegree
        //for (uint node_id = 0; node_id < graph.nodes.size()-1; node_id++) {
        for (auto it = ilp_chain.begin(); it != --ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 2);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);
            //setOutDegreeCoefficients(potEdges, followerLines, node_id);
        }

        //except the last one
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id, 2);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);
        //setOutDegreeCoefficients(potEdges, followerLines, ilp_chain.back().ch_node_id);
    }

    void addInDegreeRows(const ILP_Chain& ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        //first node is always in the simplfication so follower_indegree is unconstrained
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.front().ch_node_id, 3);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FR, 0, 0);
        //glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, -1000, 0);
        //setInOutDegreeCoefficients(potEdges, followerLines, 0);

        //all other nodes should have const * edge_indegree > follower_indegree
        //=> that means: if a node has at least edge_indegree, follower_indgree is practically unconstrained
        //for (uint node_id = 1; node_id < graph.nodes.size(); node_id++) {
        for (auto it = ++ilp_chain.begin(); it != ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 3);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, 0.0, 0.0);
            //setInDegreeCoefficients(potEdges, followerLines, node_id);
        }
        //last node could be handled like first one
    }

    void addAllDegreeRows(const ILP_data &ilp_data)
    {
        addInOutDegreeRows(ilp_data.ilp_chain1);
        addInOutDegreeRows(ilp_data.ilp_chain2);
        addOutDegreeRows(ilp_data.ilp_chain1);
        addOutDegreeRows(ilp_data.ilp_chain2);
        addInDegreeRows(ilp_data.ilp_chain1);
        addInDegreeRows(ilp_data.ilp_chain2);
    }

    void addIntersectionRows(const std::vector<Intersection> &intersections)
    {
        for (const Intersection &intersection : intersections)
        {
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

    void addAllIntersectionRows(const ILP_data &ilp_data)
    {

        addIntersectionRows(ilp_data.edgeIntersections);
        addIntersectionRows(ilp_data.followerLinesUnorderings);
    }

    void addColumns(const std::vector<Line> &lines, double obj_coef)
    {
        debug_assert(lines.size() > 0);
        glp_add_cols(lp, lines.size());
        for (uint i = 0; i < lines.size(); i++)
        {
            const Line &line = lines.at(i);
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

    void setAllColumns(const ILP_data &ilp_data)
    {
        addColumns(ilp_data.potEdges1, 1.0);
        addColumns(ilp_data.potEdges2, 1.0);
        addColumns(ilp_data.followerLines1, 0);
        addColumns(ilp_data.followerLines2, 0);
    }

    void prepareArray(const ILP_data &ilp_data)
    {
        assert(ilp_data.ilp_chain1.size() > 1);
        assert(ilp_data.ilp_chain2.size() > 1);

        preNofCols = ilp_data.potEdges1.size() + ilp_data.potEdges2.size()
                + ilp_data.followerLines1.size() + ilp_data.followerLines2.size();
        preNofRows = 3 * (ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size()) //3 degree equations for each node
                + ilp_data.edgeIntersections.size()
                + ilp_data.followerLinesUnorderings.size();
        preNofNonZeros = 4 * ilp_data.potEdges1.size() + 4 * ilp_data.potEdges2.size()
                + 2 * ilp_data.followerLines1.size() + 2 * ilp_data.followerLines2.size()
                + 2 * ilp_data.edgeIntersections.size()
                + 2 * ilp_data.followerLinesUnorderings.size();
        nofNonZeros = 0;
        //Print("preNofNonZeros" << preNofNonZeros);
        enoughSpace = preNofNonZeros < (int) size;
        assert(enoughSpace);
    }


public:

    ParallelLineSimplificationILP(const CHGraph<CHNode, CHEdge> &graph) : Calculation(graph) { }

    ~ParallelLineSimplificationILP() { }

    double solve(const Chain& chain1, const Chain& chain2, double epsilon, double eta)
    {

        ILP_data ilp_data = ILP_data(graph, chain1, chain2, epsilon, eta, true);
        prepareArray(ilp_data);
        lp = glp_create_prob();
        glp_term_out(GLP_OFF); //make glp not verbose

        glp_set_prob_name(lp, "parallel line generalization");

        glp_set_obj_dir(lp, GLP_MIN); //minimize number of used edges

        setAllColumns(ilp_data);

        addAllDegreeRows(ilp_data);
        setAllDegreeCoefficients(ilp_data);
        addAllIntersectionRows(ilp_data);

        assert(preNofRows == glp_get_num_rows(lp));
        assert(preNofNonZeros == (int) nofNonZeros);

        glp_load_matrix(lp, nofNonZeros, &ia[0], &ja[0], &ar[0]);
        //glp_load_matrix(lp, nofNonZeros, ia, ja, ar);
        //glp_write_lp(lp, NULL, "lp.txt");

        //Print("factorize" << glp_factorize(lp));
        /* solve problem */
        //Print("solving");
        //glp_simplex(lp, NULL);

        //WorkingSign();
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        glp_intopt(lp, &parm);


        switch (glp_mip_status(lp))
        {
            case GLP_UNDEF:
            {
                objective_value = std::numeric_limits<double>::max();
                break;
            }
            case GLP_OPT:
            {
                objective_value = glp_mip_obj_val(lp);
                break;
            }
            case GLP_FEAS:
            {
                objective_value = std::numeric_limits<double>::max();
                break;
            }
            case GLP_NOFEAS:
            {
                objective_value = std::numeric_limits<double>::max();
                break;
            }
            default:
                break;
        }

        //glp_print_mip(lp, "p_ilp_solution.txt");

        /* recover and display results */

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