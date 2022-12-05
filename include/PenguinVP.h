#ifndef PENGUINVP_H
#define PENGUINVP_H

#include <MutableVertexPartition.h>

class PenguinVP : public MutableVertexPartition
{
  public:
    PenguinVP(Graph* graph,
        vector<size_t> const& membership);
    PenguinVP(Graph* graph);
    virtual ~PenguinVP();
    virtual PenguinVP* create(Graph* graph);
    virtual PenguinVP* create(Graph* graph, vector<size_t> const& membership);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();
    virtual double l2_norm(double u[], int size);
    virtual double t1();
    virtual double t2();
    virtual double t3();

  protected:
  private:
};

#endif // MODULARITYVERTEXPARTITION_H