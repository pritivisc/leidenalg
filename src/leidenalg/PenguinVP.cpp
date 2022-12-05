#include "PenguinVP.h"
#include <fstream>

#ifdef DEBUG
#include <iostream>
#include <math.h>
#include <fstream>
using std::ofstream;
using std::cerr;
using std::endl;
using std::cout;
using std::basic_ofstream<char>;
#endif

PenguinVP::PenguinVP(Graph* graph,
      vector<size_t> const& membership) :
        MutableVertexPartition(graph,
        membership)
{ }

PenguinVP::PenguinVP(Graph* graph) :
        MutableVertexPartition(graph)
{ }

PenguinVP::~PenguinVP()
{ }

PenguinVP* PenguinVP::create(Graph* graph)
{
  return new PenguinVP(graph);
}

PenguinVP* PenguinVP::create(Graph* graph, vector<size_t> const& membership)
{
  return new PenguinVP(graph, membership);
}

double PenguinVP::l2_norm(double u[], int size) {
    double accum = 0.;
    for (int i = 0; i < size; ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
*****************************************************************************/

double PenguinVP::diff_move(size_t v, size_t new_comm)
{
    double K_EXP = 0.5;
    double K_COEFFICIENT = 100;
    double B_EXP = 70;

    int k = this->n_communities();
    int num_nodes = this->graph->vcount();

    int old_comm = this->_membership[v];
    if (old_comm == new_comm) {
        return 0.0;
    }

    if (this->cnodes(old_comm) == 1) {
        int nk = k - 1;
        int counts[nk];
        int index = 0;
        for (int i = 0; i < k; i++) {
            if (i != old_comm) {
                if (i == new_comm) {
                    counts[index] = this->cnodes(i) + 1;
                } else {
                    counts[index] = this->cnodes(i);
                }
                index++;
            }
        }

        double counts_processed[nk];

        for (int i = 0; i < nk; i++) {
            counts_processed[i] = (((double)counts[i]/(double)num_nodes) - ((double)1/(double)nk));
        }

        double b = l2_norm(counts_processed, nk);

        double C_w = this->total_weight_in_all_comms();
        C_w += this->weight_to_comm(v, new_comm);

        double new_cost = (C_w + K_COEFFICIENT * exp(K_EXP * nk) + exp(B_EXP * b));

        return new_cost - this->quality();

    } else {
        int counts[k];
        for (int i = 0; i < k; i++) {
            if (i == new_comm) {
                counts[i] = this->cnodes(i) + 1;
            } else if (i == old_comm) {
                counts[i] = this->cnodes(i) - 1;
            } else {
                counts[i] = this->cnodes(i);
            }
        }

        double counts_processed[k];

        for (int i = 0; i < k; i++) {
            counts_processed[i] = (((double)counts[i]/(double)num_nodes) - ((double)1/(double)k));
        }

        double b = l2_norm(counts_processed, k);

        double C_w = this->total_weight_in_all_comms();

        double w_to_new = this->weight_to_comm(v, new_comm); //sum of weights of edges from v to the new community
        double w_to_old = this->weight_to_comm(v, old_comm); //sum of weights from v to the other edges in to comm

        C_w = C_w - w_to_old + w_to_new;

        double new_cost = (C_w + K_COEFFICIENT * exp(K_EXP * k) + exp(B_EXP * b));

        return new_cost - this->quality();

    }
}

// double PenguinVP::diff_move(size_t v, size_t new_comm)
// {
    
//     int p_vector[k];
//     for (int i = 0; i < k; i++) 
//     {
//         p_vector[i] = this->cnodes(i);
//     }

//     double b_vector[k];

//     for (int i = 0; i < k; i++)
//     {
//         b_vector[i] = (((double)p_vector[i]/(double)num_nodes) - ((double)1/(double)k)); //make sure this returns a fpdivision not integer division
//     }

//     double l2n_b = l2_norm(b_vector, k);


//     double C_subp = exp(70 * l2n_b);

//     double C_subp_new = 0;

//     if (this->cnodes(old_comm) == 1) {
//         int new_p_vector[k-1];
//         int index = 0;
//         for (int i = 0; i < k; i++) {
//             if (i != old_comm) {
//                 if (i == new_comm) {
//                     new_p_vector[index] = this->cnodes(i) + 1;
//                 } else {
//                     new_p_vector[index] = this->cnodes(i);
//                 }
//                 index++;
//             }
//         }

//         double new_b_vector[k-1];

//         for (int i = 0; i < k-1; i++) {
//             new_b_vector[i] = (((double)new_p_vector[i]/(double)num_nodes) - ((double)1/(double)(k-1)));
//         }

//         double new_l2n_b = l2_norm(new_b_vector, k-1);
        
//         C_subp_new = exp(70 * new_l2n_b);

//     } else {
//         double b_i2 = (b_vector[old_comm] * b_vector[old_comm]);
//         double b_j2 = (b_vector[new_comm] * b_vector[new_comm]);
//         double b_l22 = (l2n_b*l2n_b);
//         double b_ix = (b_vector[old_comm] - (double)1/(double)num_nodes) * (b_vector[old_comm] - (double)1/(double)num_nodes);
//         double b_jx = (b_vector[new_comm] + (double)1/(double)num_nodes) * (b_vector[new_comm] + (double)1/(double)num_nodes);

//         C_subp_new = exp(70 * sqrt(b_l22 - b_i2 - b_j2 + b_ix + b_jx));
//     }

    

//     double diffC_subp = C_subp_new - C_subp;//difference in the C_subp component

//     double old_2ndterm = 100 * exp(0.5 * k);

//     double new_2ndterm = 0;
//     if (this->cnodes(old_comm) == 1)
//     {
//         new_2ndterm = 100 * exp(0.5 * k-1);
//     }
//     else
//     {
//         new_2ndterm = old_2ndterm;
//     }

//     double diff_secondterm = new_2ndterm - old_2ndterm;





//     // double C_w = this->total_weight_in_all_comms();
//     // double C_w_n = C_w - this.weight_to_comm(old)
//     // C_w_n = C_w + this.weight_to_comm(new)

//     // C_w_n - old = C_w - old + new - C_w

//     double w_to_new = this->weight_to_comm(v, new_comm); //sum of weights of edges from v to the new community
//     double w_to_old = this->weight_to_comm(v, old_comm); //sum of weights from v to the other edges in to comm


//     double diff_thirdterm = w_to_new - w_to_old;

//     return (diff_secondterm + diff_thirdterm + diffC_subp);

// }

/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/

double PenguinVP::t1() {

    std::ofstream myfile;
    myfile.open ("example1.txt");

    double sum = 0;
    int k = this->n_communities();
    for (int i = 0; i < k; i++) {
        myfile << this->total_weight_in_comm(i) << "\n";
        sum += this->total_weight_in_comm(i);
    }

    myfile.close();
    
    // return this->total_weight_in_all_comms();
    return sum;
}

double PenguinVP::t2() {
    int k = this->n_communities();
    double K_COEFFICIENT = 100;
    double K_EXP = 0.5;
    double t2 = K_COEFFICIENT * exp(K_EXP * k);
    return t2;
}

double PenguinVP::t3() {
    double B_EXP = 70;
    int num_nodes = this->graph->vcount();
    int k = this->n_communities();
    int counts[k];
    for (int i = 0; i < k; i++) 
    {
        counts[i] = this->cnodes(i);
    }
    double counts_processed[k];
    for (int i = 0; i < k; i++)
    {
        counts_processed[i] = (((double)counts[i]/(double)num_nodes) - ((double)1/(double)k)); //make sure this returns a fpdivision not integer division
    }
    double b = l2_norm(counts_processed, k);

    return exp(B_EXP * b);

}

double PenguinVP::quality()
{
    // double K_EXP = 0.5;
    // double K_COEFFICIENT = 100;
    // double B_EXP = 70;

    // int k = this->n_communities();
    // int num_nodes = this->graph->vcount();

    // int counts[k];
    // for (int i = 0; i < k; i++) 
    // {
    //     counts[i] = this->cnodes(i);
    // }

    // double counts_processed[k];
    // for (int i = 0; i < k; i++)
    // {
    //     counts_processed[i] = (((double)counts[i]/(double)num_nodes) - ((double)1/(double)k)); //make sure this returns a fpdivision not integer division
    // }

    // double b = l2_norm(counts_processed, k);

    // double C_w = this->total_weight_in_all_comms();

    // double t1 = C_w;
    // double t2 = K_COEFFICIENT * exp(K_EXP * k);
    // double t3 = exp(B_EXP * b); 
    std::ofstream myfile;
    myfile.open ("example.txt");
    myfile << t1() << "\n";
    myfile << t2() << "\n";
    myfile << t3() << "\n";
    myfile.close();



    return t1() + t2() + t3();
}
