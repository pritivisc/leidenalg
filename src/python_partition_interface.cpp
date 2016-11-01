#include "python_partition_interface.h"

Graph* create_graph_from_py(PyObject* py_obj_graph)
{
  return create_graph_from_py(py_obj_graph, NULL, NULL);
}

Graph* create_graph_from_py(PyObject* py_obj_graph, PyObject* py_weights)
{
  return create_graph_from_py(py_obj_graph, py_weights, NULL);
}

Graph* create_graph_from_py(PyObject* py_obj_graph, PyObject* py_weights, PyObject* py_node_sizes)
{
  #ifdef DEBUG
    cerr << "create_graph_from_py" << endl;
  #endif

  #if PY_MAJOR_VERSION >= 3
    igraph_t* py_graph = (igraph_t*) PyCapsule_GetPointer(py_obj_graph, NULL);
  #else
    igraph_t* py_graph = (igraph_t*) PyCObject_AsVoidPtr(py_obj_graph);
  #endif
  #ifdef DEBUG
    cerr << "Got igraph_t " << py_graph << endl;
  #endif

  // If necessary create a weighted graph
  Graph* graph = NULL;
  #ifdef DEBUG
    cerr << "Creating graph."<< endl;
  #endif

  size_t n = igraph_vcount(py_graph);
  size_t m = igraph_ecount(py_graph);

  vector<size_t> node_sizes;
  vector<double> weights;
  if (py_node_sizes != NULL && py_node_sizes != Py_None)
  {
    #ifdef DEBUG
      cerr << "Reading node_sizes." << endl;
    #endif

    size_t nb_node_size = PyList_Size(py_node_sizes);
    node_sizes.resize(n);
    for (size_t v = 0; v < nb_node_size; v++)
      node_sizes[v] = PyLong_AsLong(PyList_GetItem(py_node_sizes, v));
  }

  if (py_weights != NULL && py_weights != Py_None)
  {
    #ifdef DEBUG
      cerr << "Reading weights." << endl;
    #endif
    size_t nb_weights = PyList_Size(py_weights);
    weights.resize(m);
    for (size_t e = 0; e < nb_weights; e++)
      weights[e] = PyFloat_AsDouble(PyList_GetItem(py_weights, e));
  }

  // TODO: Pass correct_for_self_loops as parameter
  int correct_self_loops = false;
  if (node_sizes.size() == n)
  {
    if (weights.size() == m)
      graph = new Graph(py_graph, weights, node_sizes, correct_self_loops);
    else
      graph = new Graph(py_graph, node_sizes, correct_self_loops);
  }
  else
  {
    if (weights.size() == m)
      graph = new Graph(py_graph, weights, correct_self_loops);
    else
      graph = new Graph(py_graph, correct_self_loops);
  }

  #ifdef DEBUG
    cerr << "Created graph " << graph << endl;
    cerr << "Number of nodes " << graph->vcount() << endl;
    cerr << "Number of edges " << graph->ecount() << endl;
    cerr << "Total weight " << graph->total_weight() << endl;
  #endif

  return graph;
}

PyObject* capsule_MutableVertexPartition(MutableVertexPartition* partition)
{
  PyObject* py_partition = PyCapsule_New(partition, "louvain.VertexPartition.MutableVertexPartition", del_MutableVertexPartition);
  return py_partition;
}

MutableVertexPartition* decapsule_MutableVertexPartition(PyObject* py_partition)
{
  MutableVertexPartition* partition = (MutableVertexPartition*) PyCapsule_GetPointer(py_partition, "louvain.VertexPartition.MutableVertexPartition");
  return partition;
}

void del_MutableVertexPartition(PyObject* py_partition)
{
  MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);
  delete partition;
}

#ifdef __cplusplus
extern "C"
{
#endif
  PyObject* _new_ModularityVertexPartition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_obj_graph = NULL;
    PyObject* py_initial_membership = NULL;
    PyObject* py_weights = NULL;

    static char* kwlist[] = {"graph", "initial_membership", "weights", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OO", kwlist,
                                     &py_obj_graph, &py_initial_membership, &py_weights))
        return NULL;

    try
    {

      Graph* graph = create_graph_from_py(py_obj_graph, py_weights);

      ModularityVertexPartition* partition = NULL;

      // If necessary create an initial partition
      if (py_initial_membership != NULL && py_initial_membership != Py_None)
      {

        vector<size_t> initial_membership;

        #ifdef DEBUG
          cerr << "Reading initial membership." << endl;
        #endif
        size_t n = PyList_Size(py_initial_membership);
        initial_membership.resize(n);
        for (size_t v = 0; v < n; v++)
          initial_membership[v] = PyLong_AsLong(PyList_GetItem(py_initial_membership, v));

        partition = new ModularityVertexPartition(graph, initial_membership);
      }
      else
        partition = new ModularityVertexPartition(graph);

      // Do *NOT* forget to remove the graph upon deletion
      partition->destructor_delete_graph = true;

      PyObject* py_partition = capsule_MutableVertexPartition(partition);
      #ifdef DEBUG
        cerr << "Created capsule partition at address " << py_partition << endl;
      #endif

      return py_partition;
    }
    catch (std::exception const & e )
    {
      PyErr_SetString(PyExc_BaseException, "Could not constuct partition.");
      return NULL;
    }
  }

  PyObject* _new_SignificanceVertexPartition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_obj_graph = NULL;
    PyObject* py_initial_membership = NULL;

    static char* kwlist[] = {"graph", "initial_membership", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|O", kwlist,
                                     &py_obj_graph, &py_initial_membership))
        return NULL;

    try
    {

      Graph* graph = create_graph_from_py(py_obj_graph);

      SignificanceVertexPartition* partition = NULL;

      // If necessary create an initial partition
      if (py_initial_membership != NULL && py_initial_membership != Py_None)
      {

        vector<size_t> initial_membership;

        #ifdef DEBUG
          cerr << "Reading initial membership." << endl;
        #endif
        size_t n = PyList_Size(py_initial_membership);
        initial_membership.resize(n);
        for (size_t v = 0; v < n; v++)
          initial_membership[v] = PyLong_AsLong(PyList_GetItem(py_initial_membership, v));

        partition = new SignificanceVertexPartition(graph, initial_membership);
      }
      else
        partition = new SignificanceVertexPartition(graph);

      // Do *NOT* forget to remove the graph upon deletion
      partition->destructor_delete_graph = true;

      PyObject* py_partition = capsule_MutableVertexPartition(partition);
      #ifdef DEBUG
        cerr << "Created capsule partition at address " << py_partition << endl;
      #endif

      return py_partition;
    }
    catch (std::exception const & e )
    {
      PyErr_SetString(PyExc_BaseException, "Could not constuct partition.");
      return NULL;
    }
  }

  PyObject* _new_SurpriseVertexPartition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_obj_graph = NULL;
    PyObject* py_initial_membership = NULL;
    PyObject* py_weights = NULL;

    static char* kwlist[] = {"graph", "initial_membership", "weights", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OO", kwlist,
                                     &py_obj_graph, &py_initial_membership, &py_weights))
        return NULL;

    try
    {

      Graph* graph = create_graph_from_py(py_obj_graph, py_weights);

      SurpriseVertexPartition* partition = NULL;

      // If necessary create an initial partition
      if (py_initial_membership != NULL && py_initial_membership != Py_None)
      {

        vector<size_t> initial_membership;

        #ifdef DEBUG
          cerr << "Reading initial membership." << endl;
        #endif
        size_t n = PyList_Size(py_initial_membership);
        initial_membership.resize(n);
        for (size_t v = 0; v < n; v++)
          initial_membership[v] = PyLong_AsLong(PyList_GetItem(py_initial_membership, v));

        partition = new SurpriseVertexPartition(graph, initial_membership);
      }
      else
        partition = new SurpriseVertexPartition(graph);

      // Do *NOT* forget to remove the graph upon deletion
      partition->destructor_delete_graph = true;

      PyObject* py_partition = capsule_MutableVertexPartition(partition);
      #ifdef DEBUG
        cerr << "Created capsule partition at address " << py_partition << endl;
      #endif

      return py_partition;
    }
    catch (std::exception const & e )
    {
      PyErr_SetString(PyExc_BaseException, "Could not constuct partition.");
      return NULL;
    }
  }

  PyObject* _new_CPMVertexPartition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_obj_graph = NULL;
    PyObject* py_initial_membership = NULL;
    PyObject* py_weights = NULL;
    PyObject* py_node_sizes = NULL;
    double resolution_parameter = 1.0;

    static char* kwlist[] = {"graph", "initial_membership", "weights", "node_sizes", "resolution_parameter", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOOd", kwlist,
                                     &py_obj_graph, &py_initial_membership, &py_weights, &py_node_sizes, &resolution_parameter))
        return NULL;

    try
    {

      Graph* graph = create_graph_from_py(py_obj_graph, py_weights, py_node_sizes);

      CPMVertexPartition* partition = NULL;

      // If necessary create an initial partition
      if (py_initial_membership != NULL && py_initial_membership != Py_None)
      {

        vector<size_t> initial_membership;

        #ifdef DEBUG
          cerr << "Reading initial membership." << endl;
        #endif
        size_t n = PyList_Size(py_initial_membership);
        initial_membership.resize(n);
        for (size_t v = 0; v < n; v++)
          initial_membership[v] = PyLong_AsLong(PyList_GetItem(py_initial_membership, v));

        partition = new CPMVertexPartition(graph, initial_membership, resolution_parameter);
      }
      else
        partition = new CPMVertexPartition(graph, resolution_parameter);

      // Do *NOT* forget to remove the graph upon deletion
      partition->destructor_delete_graph = true;

      PyObject* py_partition = capsule_MutableVertexPartition(partition);
      #ifdef DEBUG
        cerr << "Created capsule partition at address " << py_partition << endl;
      #endif

      return py_partition;
    }
    catch (std::exception const & e )
    {
      PyErr_SetString(PyExc_BaseException, "Could not constuct partition.");
      return NULL;
    }
  }

  PyObject* _new_RBERVertexPartition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_obj_graph = NULL;
    PyObject* py_initial_membership = NULL;
    PyObject* py_weights = NULL;
    PyObject* py_node_sizes = NULL;
    double resolution_parameter = 1.0;

    static char* kwlist[] = {"graph", "initial_membership", "weights", "node_sizes", "resolution_parameter", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOOd", kwlist,
                                     &py_obj_graph, &py_initial_membership, &py_weights, &py_node_sizes, &resolution_parameter))
        return NULL;

    try
    {

      Graph* graph = create_graph_from_py(py_obj_graph, py_weights, py_node_sizes);

      RBERVertexPartition* partition = NULL;

      // If necessary create an initial partition
      if (py_initial_membership != NULL && py_initial_membership != Py_None)
      {

        vector<size_t> initial_membership;

        #ifdef DEBUG
          cerr << "Reading initial membership." << endl;
        #endif
        size_t n = PyList_Size(py_initial_membership);
        initial_membership.resize(n);
        for (size_t v = 0; v < n; v++)
          initial_membership[v] = PyLong_AsLong(PyList_GetItem(py_initial_membership, v));

        partition = new RBERVertexPartition(graph, initial_membership, resolution_parameter);
      }
      else
        partition = new RBERVertexPartition(graph, resolution_parameter);

      // Do *NOT* forget to remove the graph upon deletion
      partition->destructor_delete_graph = true;

      PyObject* py_partition = capsule_MutableVertexPartition(partition);
      #ifdef DEBUG
        cerr << "Created capsule partition at address " << py_partition << endl;
      #endif

      return py_partition;
    }
    catch (std::exception const & e )
    {
      PyErr_SetString(PyExc_BaseException, "Could not constuct partition.");
      return NULL;
    }
  }

  PyObject* _new_RBConfigurationVertexPartition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_obj_graph = NULL;
    PyObject* py_initial_membership = NULL;
    PyObject* py_weights = NULL;
    double resolution_parameter = 1.0;

    static char* kwlist[] = {"graph", "initial_membership", "weights", "resolution_parameter", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOOd", kwlist,
                                     &py_obj_graph, &py_initial_membership, &py_weights, &resolution_parameter))
        return NULL;

    try
    {

      Graph* graph = create_graph_from_py(py_obj_graph, py_weights);

      RBConfigurationVertexPartition* partition = NULL;

      // If necessary create an initial partition
      if (py_initial_membership != NULL && py_initial_membership != Py_None)
      {

        vector<size_t> initial_membership;

        #ifdef DEBUG
          cerr << "Reading initial membership." << endl;
        #endif
        size_t n = PyList_Size(py_initial_membership);
        initial_membership.resize(n);
        for (size_t v = 0; v < n; v++)
          initial_membership[v] = PyLong_AsLong(PyList_GetItem(py_initial_membership, v));

        partition = new RBConfigurationVertexPartition(graph, initial_membership, resolution_parameter);
      }
      else
        partition = new RBConfigurationVertexPartition(graph, resolution_parameter);

      // Do *NOT* forget to remove the graph upon deletion
      partition->destructor_delete_graph = true;

      PyObject* py_partition = capsule_MutableVertexPartition(partition);
      #ifdef DEBUG
        cerr << "Created capsule partition at address " << py_partition << endl;
      #endif

      return py_partition;
    }
    catch (std::exception const & e )
    {
      PyErr_SetString(PyExc_BaseException, "Could not constuct partition.");
      return NULL;
    }
  }

  PyObject* _MutableVertexPartition_get_py_igraph(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;

    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "get_py_igraph();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    Graph* graph = partition->get_graph();

    size_t n = graph->vcount();
    size_t m = graph->ecount();

    PyObject* edges = PyList_New(m);
    for (size_t e = 0; e < m; e++)
    {
      vector<size_t> edge = graph->edge(e);
      PyList_SetItem(edges, e, Py_BuildValue("(KK)", edge[0], edge[1]));
    }

    PyObject* weights = PyList_New(m);
    for (size_t e = 0; e < m; e++)
    {
      PyObject* item = PyFloat_FromDouble(graph->edge_weight(e));
      PyList_SetItem(weights, e, item);
    }

    PyObject* node_sizes = PyList_New(n);
    for (size_t v = 0; v < n; v++)
    {
      #if PY_MAJOR_VERSION >= 3
        PyObject* item = PyLong_FromSize_t(graph->node_size(v));
      #else
        PyObject* item = PyInt_FromSize_t(graph->node_size(v));
      #endif
      PyList_SetItem(node_sizes, v, item);
    }

    return Py_BuildValue("lOOO", n, edges, weights, node_sizes);
  }

  PyObject* _MutableVertexPartition_from_coarse_partition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    PyObject* py_membership = NULL;
    PyObject* py_coarse_node = NULL;

    static char* kwlist[] = {"partition", "membership", "coarse_node", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    // TODO : Instead of simply returning NULL, we should also set an error.
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|O", kwlist,
                                     &py_partition, &py_membership, &py_coarse_node))
        return NULL;

    #ifdef DEBUG
      cerr << "from_coarse_partition();" << endl;
    #endif

    size_t n = PyList_Size(py_membership);
    vector<size_t> membership;
    membership.resize(n);
    for (size_t v = 0; v < n; v++)
      membership[v] = PyLong_AsLong(PyList_GetItem(py_membership, v));

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    if (py_coarse_node != NULL && py_coarse_node != Py_None)
    {
      cerr << "Get coarse node list" << endl;
      size_t n = PyList_Size(py_coarse_node);
      vector<size_t> coarse_node;
      membership.resize(n);
      for (size_t v = 0; v < n; v++)
        coarse_node[v] = PyLong_AsLong(PyList_GetItem(py_coarse_node, v));

    cerr << "Got coarse node list" << endl;
      partition->from_coarse_partition(membership, coarse_node);
    }
    else
      partition->from_coarse_partition(membership);

    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* _MutableVertexPartition_renumber_communities(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;

    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "renumber_communities();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    partition->renumber_communities();

    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* _MutableVertexPartition_diff_move(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t v;
    size_t new_comm;

    static char* kwlist[] = {"partition", "v", "new_comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OKK", kwlist,
                                     &py_partition, &v, &new_comm))
        return NULL;

    #ifdef DEBUG
      cerr << "diff_move(" << v << ", " << new_comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double diff = partition->diff_move(v, new_comm);
    return PyFloat_FromDouble(diff);
  }

  PyObject* _MutableVertexPartition_move_node(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t v;
    size_t new_comm;

    static char* kwlist[] = {"partition", "v", "new_comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OKK", kwlist,
                                     &py_partition, &v, &new_comm))
        return NULL;

    #ifdef DEBUG
      cerr << "move_node(" << v << ", " << new_comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    partition->move_node(v, new_comm);

    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* _MutableVertexPartition_quality(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;

    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "quality();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double q = partition->quality();
    return PyFloat_FromDouble(q);
  }

  PyObject* _MutableVertexPartition_aggregate_partition(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;

    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "aggregate_partition();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    // First collapse graph (i.e. community graph)
    Graph* collapsed_graph = partition->get_graph()->collapse_graph(partition);

    // Create collapsed partition (i.e. default partition of each node in its own community).
    MutableVertexPartition* collapsed_partition = partition->create(collapsed_graph);
    collapsed_partition->destructor_delete_graph = true;

    PyObject* py_collapsed_partition = capsule_MutableVertexPartition(collapsed_partition);
    return py_collapsed_partition;
  }

  PyObject* _MutableVertexPartition_total_weight_in_comm(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t comm;

    static char* kwlist[] = {"partition", "comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OK", kwlist,
                                     &py_partition, &comm))
        return NULL;

    #ifdef DEBUG
      cerr << "total_weight_in_comm(" << comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    if (comm >= partition->nb_communities())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of communities.");
      return NULL;
    }

    double w = partition->total_weight_in_comm(comm);
    return PyFloat_FromDouble(w);
  }

  PyObject* _MutableVertexPartition_total_weight_from_comm(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t comm;

    static char* kwlist[] = {"partition", "comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OK", kwlist,
                                     &py_partition, &comm))
        return NULL;

    #ifdef DEBUG
      cerr << "total_weight_from_comm(" << comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    if (comm >= partition->nb_communities())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of communities.");
      return NULL;
    }

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double w = partition->total_weight_from_comm(comm);
    return PyFloat_FromDouble(w);
  }

  PyObject* _MutableVertexPartition_total_weight_to_comm(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t comm;

    static char* kwlist[] = {"partition", "comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OK", kwlist,
                                     &py_partition, &comm))
        return NULL;

    #ifdef DEBUG
      cerr << "total_weight_to_comm(" << comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    if (comm >= partition->nb_communities())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of communities.");
      return NULL;
    }

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double w = partition->total_weight_to_comm(comm);
    return PyFloat_FromDouble(w);
  }

  PyObject* _MutableVertexPartition_total_weight_in_all_comms(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;

    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "total_weight_in_all_comms();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double w = partition->total_weight_in_all_comms();
    return PyFloat_FromDouble(w);
  }

  PyObject* _MutableVertexPartition_total_possible_edges_in_all_comms(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;

    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "total_possible_edges_in_all_comms();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    size_t e = partition->total_possible_edges_in_all_comms();
    return PyLong_FromSize_t(e);
  }

  PyObject* _MutableVertexPartition_weight_to_comm(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t v;
    size_t comm;

    static char* kwlist[] = {"partition", "v", "comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OKK", kwlist,
                                     &py_partition, &v, &comm))
        return NULL;

    #ifdef DEBUG
      cerr << "weight_to_comm(" << v << ", " << comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    if (comm >= partition->nb_communities())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of communities.");
      return NULL;
    }

    if (v >= partition->get_graph()->vcount())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of nodes.");
      return NULL;
    }

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double diff = partition->weight_to_comm(v, comm);
    return PyFloat_FromDouble(diff);
  }

  PyObject* _MutableVertexPartition_weight_from_comm(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    size_t v;
    size_t comm;

    static char* kwlist[] = {"partition", "v", "comm", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OKK", kwlist,
                                     &py_partition, &v, &comm))
        return NULL;

    #ifdef DEBUG
      cerr << "weight_to_comm(" << v << ", " << comm << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    if (comm >= partition->nb_communities())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of communities.");
      return NULL;
    }

    if (v >= partition->get_graph()->vcount())
    {
      PyErr_SetString(PyExc_IndexError, "Try to index beyond the number of nodes.");
      return NULL;
    }

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double diff = partition->weight_to_comm(v, comm);
    return PyFloat_FromDouble(diff);
  }

  PyObject* _MutableVertexPartition_get_membership(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "get_membership();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    size_t n = partition->get_graph()->vcount();
    PyObject* py_membership = PyList_New(n);
    for (size_t v = 0; v < n; v++)
    {
      #if PY_MAJOR_VERSION >= 3
        PyObject* item = PyLong_FromSize_t(partition->membership(v));
      #else
        PyObject* item = PyInt_FromSize_t(partition->membership(v));
      #endif
      PyList_SetItem(py_membership, v, item);
    }
    return py_membership;
  }

  PyObject* _MutableVertexPartition_set_membership(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    PyObject* py_membership = NULL;

    static char* kwlist[] = {"partition", "membership", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO", kwlist,
                                     &py_partition, &py_membership))
        return NULL;

    #ifdef DEBUG
      cerr << "set_membership();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    MutableVertexPartition* partition = decapsule_MutableVertexPartition(py_partition);

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    size_t n = PyList_Size(py_membership);
    vector<size_t> membership;
    membership.resize(n);
    for (size_t v = 0; v < n; v++)
      membership[v] = PyLong_AsLong(PyList_GetItem(py_membership, v));

    partition->set_membership(membership);

    #ifdef DEBUG
      cerr << "Exiting set_membership();" << endl;
    #endif

    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* _ResolutionParameterVertexPartition_get_resolution(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    static char* kwlist[] = {"partition", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O", kwlist,
                                     &py_partition))
        return NULL;

    #ifdef DEBUG
      cerr << "get_resolution();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule ResolutionParameterVertexPartition at address " << py_ResolutionParameterVertexPartition << endl;
    #endif
    ResolutionParameterVertexPartition* partition = (ResolutionParameterVertexPartition*)decapsule_MutableVertexPartition(py_partition);
    #ifdef DEBUG
      cerr << "Using ResolutionParameterVertexPartition at address " << ResolutionParameterVertexPartition << endl;
    #endif

    return PyFloat_FromDouble(partition->resolution_parameter);
  }

  PyObject* _ResolutionParameterVertexPartition_set_resolution(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    double resolution_parameter = 1.0;
    static char* kwlist[] = {"partition", "resolution_parameter", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "Od", kwlist,
                                     &py_partition, &resolution_parameter))
        return NULL;

    #ifdef DEBUG
      cerr << "set_resolution(" << resolution << ");" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule ResolutionParameterVertexPartition at address " << py_ResolutionParameterVertexPartition << endl;
    #endif
    ResolutionParameterVertexPartition* partition = (ResolutionParameterVertexPartition*)decapsule_MutableVertexPartition(py_partition);
    #ifdef DEBUG
      cerr << "Using ResolutionParameterVertexPartition at address " << ResolutionParameterVertexPartition << endl;
    #endif

    partition->resolution_parameter = resolution_parameter;

    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* _ResolutionParameterVertexPartition_quality(PyObject *self, PyObject *args, PyObject *keywds)
  {
    PyObject* py_partition = NULL;
    PyObject* py_res = NULL;
    double resolution_parameter = 0.0;

    static char* kwlist[] = {"partition", "resolution_parameter", NULL};

    #ifdef DEBUG
      cerr << "Parsing arguments..." << endl;
    #endif

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|O", kwlist,
                                     &py_partition, &py_res))
        return NULL;

    #ifdef DEBUG
      cerr << "quality();" << endl;
    #endif

    #ifdef DEBUG
      cerr << "Capsule partition at address " << py_partition << endl;
    #endif

    ResolutionParameterVertexPartition* partition = (ResolutionParameterVertexPartition*)decapsule_MutableVertexPartition(py_partition);

    if (py_res != NULL && py_res != Py_None)
      resolution_parameter = PyFloat_AsDouble(py_res);
    else
      resolution_parameter = partition->resolution_parameter;

    #ifdef DEBUG
      cerr << "Using partition at address " << partition << endl;
    #endif

    double q = partition->quality(resolution_parameter);
    return PyFloat_FromDouble(q);
  }

#ifdef __cplusplus
}
#endif
