/*
 * File:  calculation.h
 * Author: tobias
 *
 * Created on July 30, 2017, 1:23 PM
 */

#ifndef CALCULATION_H
#define	CALCULATION_H

class Calculation
{
protected:
    const CHGraph<CHNode, CHEdge> &graph;

    glp_prob *lp;

    const static size_t size = 4000000;
    size_t nofNonZeros = 0;

    //numbers of needed sizes beforehand
    int preNofRows;
    int preNofCols;
    int preNofNonZeros = 0;

    //whole ILP is invalid if size is not big enough
    bool enoughSpace;

    std::vector<int> ia;
    std::vector<int> ja;
    std::vector<double> ar;

    double objective_value = std::numeric_limits<double>::max();

    bool isInRange(size_t index)
    {
        return (1 <= index && index <= size); //0 is not valid for ILP
    }

    void setDegreeCoefficient(uint i, uint j, double r)
    {
        nofNonZeros++;
        size_t index = nofNonZeros;
        assert(isInRange(index));
        ia[index] = i;
        ja[index] = j;
        ar[index] = r;
    }

    void setRowName(NodeID node_id, uint nr)
    {
        std::stringstream ss("");
        ss << "node: " << nr << ". " << node_id;
        std::string s = ss.str();
        char const *row_name = s.c_str();
        glp_set_row_name(lp, glp_get_num_rows(lp), row_name);
    }

    void setIntersectionEntry(LineID line_id)
    {
        nofNonZeros++;
        size_t index = nofNonZeros;
        assert(isInRange(index));
        ia[index] = glp_get_num_rows(lp);
        ja[index] = line_id + 1;
        ar[index] = 1.0;
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

    double solve()
    {
        glp_load_matrix(lp, nofNonZeros, &ia[0], &ja[0], &ar[0]);
        //glp_load_matrix(lp, nofNonZeros, ia, ja, ar);
        glp_write_lp(lp, NULL, "../logs/ILP/lp.txt");

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

        glp_print_mip(lp, "../logs/ILP/solution.txt");

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

    Calculation(const CHGraph<CHNode, CHEdge> &graph) : graph(graph)
    {
        ia.resize(size + 1);
        ja.resize(size + 1);
        ar.resize(size + 1);
    }

    virtual ~Calculation() { }
public:

};


#endif	/* CALCULATION_H */

