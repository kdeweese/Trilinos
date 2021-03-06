/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Functions contained in this file:
 *	main()
 *	print_input()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <cstddef> // for size_t
#include <cstdio>  // for printf, nullptr, fprintf, etc
#include <cstdlib> // for free, exit, malloc
#include <cstring> // for strcmp
#include <iostream>
#include <stdexcept>

#include "add_to_log.h" // for add_to_log
#include "elb.h"        // for LB_Description<INT>, get_time, etc
#include "elb_allo.h"   // for array_alloc
#include "elb_elem.h"   // for E_Type, ::NULL_EL
#include "elb_err.h"    // for error_report, Gen_Error, etc
#include "elb_exo.h"    // for init_weight_struct, etc
#include "elb_format.h"
#include "elb_graph.h"   // for generate_graph
#include "elb_inp.h"     // for check_inp_specs, etc
#include "elb_loadbal.h" // for generate_loadbal, etc
#include "elb_output.h"  // for write_nemesis, write_vis

#ifdef USE_ZOLTAN
#include <mpi.h> // for MPI_Finalize, etc
#endif

#ifdef SGI10K
#include <sys/resource.h>
#endif

namespace {
  template <typename INT>
  void print_input(Machine_Description * /*machine*/, LB_Description<INT> * /*lb*/,
                   Problem_Description * /*prob*/, Solver_Description * /*solver*/,
                   Weight_Description<INT> * /*weight*/);
} // namespace

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Routine which reads in a EXODUS II mesh database and load balances the
 * mesh by either a nodal or element assignment.
 *----------------------------------------------------------------------------
 * Functions called:
 *	cmd_line_arg_parse()
 *	read_cmd_file()
 *	check_inp_specs()
 *	print_input()
 *	read_exo_weights()
 *	read_mesh_params()
 *      read_mesh()
 *	generate_gaph()
 *	generate_loadbal()
 *	write_vis()
 *	generate_maps()
 *	write_nemesis()
 *----------------------------------------------------------------------------
 * Variable Index:
 *	machine       - structure containing information about the machine
 *			for which the load balance is to be generated.
 *	lb	      - structure containing information about what type
 *			of load balance is to be performed.
 *	problem	      - structure containing information about the problem
 *			type. Currently whether the problem is a nodal or
 *			elemental based decomposition.
 *	solver	      - structure containing various parameters for the
 *			eigensolver used by Chaco.
 *	weight	      - structure used to store parameters about how to
 *			weight the resulting graph before it is fed into
 *			Chaco.
 *	graph	      - structure containing the graph of the problem.
 *	mesh	      - structure containing a description of the FEM mesh.
 *	exoII_inp_file - name of the ExodusII file containing the problem
 *			 geometry.
 *	ascii_inp_file - name of the input command file.
 *	nemI_out_file  - name of the output NemesisI file.
 *
 *****************************************************************************/
template <typename INT> int internal_main(int argc, char *argv[], INT /* dummy */);

int main(int argc, char *argv[])
{
  std::cerr << "Beginning nem_slice execution." << '\n';

  double start_time = get_time();

  /* Initialize to just reporting the error */
  error_lev = 1;

  /* check if the user just wants to know the version (or forcing 64-bit mode)*/
  bool int64com = false;
  bool int32com = false;
  int  int64db  = 0;

  for (int cnt = 0; cnt < argc; cnt++) {
    if (strcmp(argv[cnt], "-V") == 0) {
      printf("%s version %s\n", UTIL_NAME, ELB_VERSION);
      exit(0);
    }
    if (strcmp(argv[cnt], "-64") == 0) {
      int64com = true;
    }

    if (strcmp(argv[cnt], "-32") == 0) {
      int32com = true;
    }
  }

  if (argc > 2) {
    /* Get the input mesh file so we can determine the integer size... */
    /* Should be the last item on the command line */
    char *mesh_file_name = argv[argc - 1];

    /* Open file and get the integer size... */
    int   cpu_ws = 0;
    int   io_ws  = 0;
    float vers   = 00.0;
    std::cerr << "Input Mesh File = '" << mesh_file_name << "'" << '\n';
    int exoid = ex_open(mesh_file_name, EX_READ, &cpu_ws, &io_ws, &vers);
    if (exoid < 0) {
      std::string error("fatal: unable to open input ExodusII file ");
      error += mesh_file_name;
      Gen_Error(0, error);
      return 0;
    }

    int64db = ex_int64_status(exoid) & EX_ALL_INT64_DB;
    ex_close(exoid);

    ex_opts(EX_VERBOSE);
  }

  int status;
  if (int32com && (int64db != 0)) {
    std::cerr << "Forcing 32-bit integer mode for decomposition even though database is 64-bit.\n";
    status = internal_main(argc, argv, int(0));
  }
  else if ((int64db != 0) || int64com) {
    std::cerr << "Using 64-bit integer mode for decomposition...\n"
              << "NOTE: Only 'linear' and 'scattered' methods are supported for 64-bit models\n";

    status = internal_main(argc, argv, int64_t(0));
  }
  else {
    std::cerr << "Using 32-bit integer mode for decomposition...\n";
    status = internal_main(argc, argv, int(0));
  }

  /* Report any non-fatal errors that may have occurred */
  error_report();

  /* Get ending time */
  double end_time = get_time();
  std::cerr << "The entire load balance took " << end_time - start_time << " seconds.\n";
  add_to_log(argv[0], end_time - start_time);
  return status;
}

template <typename INT> int internal_main(int argc, char *argv[], INT /* dummy */)
{
  /* Local variables */
  std::string exoII_inp_file;
  std::string ascii_inp_file;
  std::string nemI_out_file;

  Machine_Description     machine;
  LB_Description<INT>     lb;
  Problem_Description     problem;
  Solver_Description      solver;
  Weight_Description<INT> weight;
  Mesh_Description<INT>   mesh;
  Sphere_Info             sphere;
  Graph_Description<INT>  graph;

/*-----------------------------Execution Begins------------------------------*/
#ifdef USE_ZOLTAN
  MPI_Init(&argc, &argv);
#endif

  mesh.title[0] = '\0';

  if (sizeof(INT) == 8) {
    problem.int64api = EX_ALL_INT64_API;
  }

  /* Parse the command line */
  if (!cmd_line_arg_parse(argc, argv, exoII_inp_file, ascii_inp_file, nemI_out_file, &machine, &lb,
                          &problem, &solver, &weight)) {
    fprintf(stderr, "error parsing command line\n");
    error_report();
    exit(1);
  }

  /*
   *If the information is to be taken from an ASCII input file then
   * read that file.
   */
  if (!ascii_inp_file.empty()) {
    if (!read_cmd_file(ascii_inp_file, exoII_inp_file, nemI_out_file, &machine, &lb, &problem,
                       &solver, &weight)) {
      fprintf(stderr, "error parsing command file\n");
      error_report();
      exit(1);
    }
  }

  /* make sure that this type is set */
  if (weight.type < 0) {
    weight.type = NO_WEIGHT;
  }

  /*
   * Perform at least some rudimentary error checks on the user
   * specified input.
   */
  if (!check_inp_specs(exoII_inp_file, nemI_out_file, &machine, &lb, &problem, &solver, &weight)) {
    fprintf(stderr, "Error in user specified parameters\n");
    error_report();
    exit(1);
  }

  /* Output the parameters for the run to the screen */
  print_input(&machine, &lb, &problem, &solver, &weight);

  /* Read in mesh parameters */
  double time1 = get_time();
  if (!read_mesh_params(exoII_inp_file, &problem, &mesh, &sphere)) {
    fprintf(stderr, "Error reading mesh parameters\n");
    error_report();
    exit(1);
  }
  double time2 = get_time();
  printf("Time to read mesh parameters: %fs\n", time2 - time1);

  /* Check for error conditions */
  if ((mesh.num_nodes) / (machine.num_procs) < 1) {
    Gen_Error(0, "fatal: problem divided among too many processors");
    error_report();
    exit(1);
  }
  else if ((problem.type == ELEMENTAL) && ((mesh.num_elems) / (machine.num_procs) < 1)) {
    Gen_Error(0, "fatal: problem divided among too many processors");
    error_report();
    exit(1);
  }

  /* If vertex weighting is turned on, prepare weight struct */
  if ((weight.type & READ_EXO) || (weight.type & EL_BLK)) {
    time1 = get_time();
    if (!init_weight_struct(&problem, &mesh, &weight)) {
      fprintf(stderr, "Error during initialization of weight struct\n");
      error_report();
      exit(1);
    }
    time2 = get_time();
    printf("Time to initialize weight struct: %fs\n", time2 - time1);
  }

  /* If desired, read in the weighting factors from the ExodusII file */
  if (weight.type & READ_EXO) {
    time1 = get_time();
    if (!read_exo_weights(&problem, &weight)) {
      fprintf(stderr, "Error during read of ExodusII weights\n");
      error_report();
      exit(1);
    }
    time2 = get_time();
    printf("Time to read ExodusII weights: %fs\n", time2 - time1);
  }

  /* Initialize various parameters */
  if (lb.type == INERTIAL || lb.type == ZPINCH || lb.type == BRICK || lb.type == ZOLTAN_RCB ||
      lb.type == ZOLTAN_RIB || lb.type == ZOLTAN_HSFC || problem.vis_out == 1 ||
      problem.vis_out == 2) {
    problem.read_coords = ELB_TRUE;
  }
  else {
    problem.read_coords = ELB_FALSE;
  }

  if (lb.type != SPECTRAL) {
    problem.coarse_flag = ELB_FALSE;
  }
  else {
    problem.coarse_flag = ELB_TRUE;
  }

  if (lb.refine == KL_REFINE) {
    problem.alloc_graph = ELB_TRUE;
  }
  else if (lb.type == SPECTRAL) {
    problem.alloc_graph = ELB_TRUE;
  }
  else {
    problem.alloc_graph = ELB_FALSE;
  }

  /* if fix_columns is on, we need the face adjacency graph. So if
   * nothing else has asked for the full adjacency graph, ask for the
   * face adjacency graph. If something else did ask for the adjacency
   * graph, we don't know if its full or face adjacency only, so leave
   * the option as is */

  if (problem.fix_columns) {
    if (problem.alloc_graph == ELB_FALSE) {
      problem.alloc_graph = ELB_TRUE;
      problem.face_adj    = ELB_TRUE;
    }
  }

  /* Allocate necessary memory */
  if (problem.type == NODAL) {
    problem.num_vertices = mesh.num_nodes;
  }
  else if (problem.type == ELEMENTAL) {
    problem.num_vertices = (mesh.num_elems - sphere.num);
  }

  if (problem.read_coords == ELB_TRUE) {
    size_t mem_req = (size_t)(mesh.num_dims) * (mesh.num_nodes) * sizeof(float);
    mesh.coords    = reinterpret_cast<float *>(malloc(mem_req));
    if (!(mesh.coords)) {
      Gen_Error(0, "fatal: insufficient memory for coordinates");
      error_report();
      exit(1);
    }
  }
  else {
    mesh.coords = nullptr;
  }

  mesh.elem_type = static_cast<E_Type *>(array_alloc(1, mesh.num_elems, sizeof(E_Type)));
  mesh.connect = static_cast<INT **>(array_alloc(2, mesh.num_elems, mesh.max_np_elem, sizeof(INT)));
  if (!(mesh.elem_type) || !(mesh.connect)) {
    Gen_Error(0, "fatal: insufficient memory");
    error_report();
    exit(1);
  }

  /* Read the mesh */
  time1 = get_time();
  if (!read_mesh(exoII_inp_file, &problem, &mesh, &weight)) {
    Gen_Error(0, "fatal: could not read mesh");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to read mesh: %fs\n", time2 - time1);

  /* free unneeded memory */
  vec_free(weight.ow);
  vec_free(weight.elemblk);
  vec_free(weight.elemblk_wgt);

  /* Generate the graph for the mesh */
  time1 = get_time();
  if (!generate_graph(&problem, &mesh, &graph, &weight, &sphere)) {
    Gen_Error(0, "fatal: could not generate graph");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to generate graph: %fs\n", time2 - time1);

  /* Generate load balance */
  try {
    time1 = get_time();
    if (!generate_loadbal(&machine, &problem, &mesh, &lb, &solver, &graph, &weight, &sphere, argc,
                          argv)) {
      Gen_Error(0, "fatal: could not generate load balance");
      error_report();
      exit(1);
    }
    time2 = get_time();
    printf("Time to generate load balance: %fs\n", time2 - time1);
  }
  catch (const std::exception &e) {
    std::cerr << "NEM_SLICE: Exception in generate_loadbal: " << e.what();
  }

  /* free up memory */
  if (sphere.adjust) {
    free(sphere.adjust);
  }

#ifdef PRINT_VERT
  for (size_t cnt = 0; cnt < problem.num_vertices; cnt++)
    printf("element = %lu, proc = %i\n", cnt, lb.vertex2proc[cnt]);
#endif

  /*
   * NOTE: in Chaco, if FREE_GRAPH is set to 1, the following arrays
   * are freed: graph.start
   *		graph.adj
   *		weight.vertices
   *		weight.edges
   *
   * Need to take into account special case where the mesh contains
   * only spheres. In this case Chaco is not called, and the above
   * arrays are not freed
   */

  if (sphere.num >= mesh.num_elems) {
    vec_free(graph.start);
    vec_free(graph.adj);
    vec_free(weight.vertices);
    vec_free(weight.edges);
  }

  /* Generate the load balance maps */
  time1 = get_time();
  if (!generate_maps(&machine, &problem, &mesh, &lb, &graph)) {
    Gen_Error(0, "fatal: could not generate load-balance maps");
    error_report();
    exit(1);
  }
  time2 = get_time();
  printf("Time to generate load-balance maps: %fs\n", time2 - time1);

  /* Output the visualization file */
  if (problem.vis_out == 1 || problem.vis_out == 2) {
    time1 = get_time();
    if (!write_vis(nemI_out_file, exoII_inp_file, &machine, &problem, &mesh, &lb)) {
      Gen_Error(0, "warning: unable to output visualization file");
    }
    time2 = get_time();
    printf("Time to output visualization file: %fs\n", time2 - time1);
  }

  /* Free up memory no longer needed */
  if (problem.read_coords == ELB_TRUE) {
    free(mesh.coords);
  }

  free(mesh.elem_type);
  free(mesh.connect);

  if (!graph.sur_elem.empty()) {
    for (size_t cnt = 0; cnt < mesh.num_nodes; cnt++) {
      vec_free(graph.sur_elem[cnt]);
    }
    vec_free(graph.sur_elem);
  }

  /* Output a Nemesis load balance file */
  time1 = get_time();
  if (!write_nemesis(nemI_out_file, &machine, &problem, &mesh, &lb, &sphere)) {
    Gen_Error(0, "fatal: could not output Nemesis file");
    error_report();
    exit(1);
  }
  time2 = get_time() - time1;
  printf("Time to write Nemesis file: %fs\n", time2);

  /* Free up unused memory for leak checking */
  for (int cnt = 0; cnt < machine.num_procs; cnt++) {
    vec_free(lb.int_nodes[cnt]);
    vec_free(lb.bor_nodes[cnt]);
    if (problem.type == NODAL) {
      vec_free(lb.ext_nodes[cnt]);
      vec_free(lb.ext_procs[cnt]);
    }

    vec_free(lb.int_elems[cnt]);
    if (problem.type == ELEMENTAL) {
      vec_free(lb.bor_elems[cnt]);
    }
  }

  vec_free(lb.int_nodes);

  if (problem.type == NODAL) {
    vec_free(lb.ext_nodes);
    vec_free(lb.ext_procs);
  }

  vec_free(lb.int_elems);

  if (problem.type == ELEMENTAL) {
    vec_free(lb.bor_elems);
    for (int cnt = 0; cnt < machine.num_procs; cnt++) {
      for (size_t cnt1 = 0; cnt1 < lb.bor_nodes[cnt].size(); cnt1++) {
        vec_free(lb.born_procs[cnt][cnt1]);
      }

      if (!lb.born_procs[cnt].empty()) {
        vec_free(lb.born_procs[cnt]);
      }
      vec_free(lb.ext_procs[cnt]);
      vec_free(lb.e_cmap_elems[cnt]);
      vec_free(lb.e_cmap_sides[cnt]);
      vec_free(lb.e_cmap_procs[cnt]);
      vec_free(lb.e_cmap_neigh[cnt]);
    }
    vec_free(lb.e_cmap_elems);
    vec_free(lb.born_procs);
    vec_free(lb.ext_procs);
  }

  vec_free(lb.bor_nodes);
  free(lb.vertex2proc);

#ifdef USE_ZOLTAN
  MPI_Finalize();
#endif

  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This function prints out parameters as read from the command line and/or
 * the input ASCII load-balance file.
 *---------------------------------------------------------------------------*/
namespace {
  template <typename INT>
  void print_input(Machine_Description *machine, LB_Description<INT> *lb, Problem_Description *prob,
                   Solver_Description *solver, Weight_Description<INT> *weight)
  {
    printf("%s version %s\n", UTIL_NAME, ELB_VERSION);

    printf("Performing ");
    switch (prob->type) {
    case NODAL: printf("a nodal "); break;

    case ELEMENTAL: printf("an elemental "); break;
    }

    printf("load balance with the following parameters...\n");

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                         Machine_Description PARAMETERS                                */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    printf("Machine Description\n");
    if (machine->num_boxes > 1) {
      printf("\tCluster of %d boxes\n", machine->num_boxes);
      printf("\tarchitecture of each box: ");
    }
    else {
      printf("\tarchitecture: ");
    }
    switch (machine->type) {
    case HCUBE: printf("hypercube\n"); break;

    case MESH: printf("mesh\n"); break;
    }
    if (machine->num_boxes > 1) {
      printf("\tdimension(s) of each box: ");
    }
    else {
      printf("\tdimension(s): ");
    }
    switch (machine->type) {
    case HCUBE: printf("%d\n", machine->dim[0]); break;

    case MESH:
      for (int cnt = 0; cnt < (machine->num_dims) - 1; cnt++) {
        printf("%dx", machine->dim[cnt]);
      }

      printf("%d\n", machine->dim[(machine->num_dims) - 1]);
      break;
    }
    printf("\ttotal number of processors: %d\n", machine->num_procs);

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                       LOAD BALANCE PARAMETERS                             */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    printf("Load Balance Parameters\n");
    switch (lb->type) {
    case MULTIKL:
      printf("\ttype: multilevel\n");
      printf("\tnumber of sections: %d\n", lb->num_sects);
      break;

    case SPECTRAL:
      printf("\ttype: spectral\n");
      printf("\tnumber of sections: %d\n", lb->num_sects);
      break;

    case INERTIAL: printf("\ttype: inertial\n"); break;

    case ZPINCH: printf("\ttype: zpinch\n"); break;

    case BRICK: printf("\ttype: brick\n"); break;

    case ZOLTAN_RCB: printf("\ttype: rcb\n"); break;

    case ZOLTAN_RIB: printf("\ttype: rib\n"); break;

    case ZOLTAN_HSFC: printf("\ttype: hsfc\n"); break;

    case LINEAR: printf("\ttype: linear\n"); break;

    case RANDOM: printf("\ttype: random\n"); break;

    case SCATTERED: printf("\ttype: scattered\n"); break;

    case INFILE:
      printf("\ttype: input from file\n");
      printf("\tfile name: %s\n", lb->file.c_str());
      break;
    }
    printf("\trefinement: ");
    switch (lb->refine) {
    case KL_REFINE: printf("Kernighan-Lin\n"); break;

    case NO_REFINE: printf("none\n"); break;
    }
    if (lb->cnctd_dom) {
      printf("\tConnected Domain enforced\n");
    }
    if (lb->outfile) {
      printf("\toutput assignment vector\n");
      printf("\tfile name: %s\n", lb->file.c_str());
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                         EIGENSOLVER PARAMETERS                            */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    if (lb->type == MULTIKL || lb->type == SPECTRAL) {
      printf("Eigensolver Parameters\n");
      printf("\teignsolver tolerance: %f\n", solver->tolerance);
      if (solver->rqi_flag == USE_RQI) {
        printf("\tusing RQI/Symmlq eigensolver\n");
        printf("\tnumber of vertices to coarsen down to: %d\n", solver->vmax);
      }
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                          WEIGHTING PARAMETERS                             */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    printf("Weighting Parameters\n");
    if (weight->type == NO_WEIGHT) {
      printf("\tno weighting\n");
    }
    else if (weight->type & READ_EXO) {
      printf("\tweights from: ExodusII file\n");
      printf("\ttime index: %d\n", weight->exo_tindx);
      printf("\tvariable index: %d\n", weight->exo_vindx);
      printf("\tvariable name: %s\n", weight->exo_varname.c_str());
    }
    else if (weight->type & EL_BLK) {
      printf("\tElement Block weights specified\n");
      for (size_t cnt = 0; cnt < weight->elemblk.size(); cnt++) {
        printf("\telement block: " ST_ZU ", weight: " ST_ZU "\n", (size_t)weight->elemblk[cnt],
               (size_t)weight->elemblk_wgt[cnt]);
      }
    }
    else if (weight->type & EDGE_WGT) {
      printf("\tedge weights turned on\n");
    }

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                          WEIGHTING PARAMETERS                             */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    printf("Miscellaneous Options\n");
    if (prob->face_adj == 1) {
      printf("\tusing face definition of adjacency\n");
    }
    if (prob->global_mech == 1) {
      printf("\tlooking for global mechanisms\n");
    }
    if (prob->local_mech == 1) {
      printf("\tlooking for local mechanisms\n");
    }
    if (prob->find_cnt_domains == 1) {
      printf("\tidentifying the number of disconnected element blocks on a subdomain\n");
    }
    if (prob->mech_add_procs == 1) {
      printf("\tincreasing the number of processors if mechanisms are found\n");
    }
    if (prob->dsd_add_procs == 1) {
      printf("\tincreasing the number of processors if disconnected sudomains are found\n");
    }
    if (prob->no_sph == 1) {
      printf("\tSPHERES are being treated as concentrated mass - connectivity exists\n");
    }
    if (prob->skip_checks == 1) {
      printf("\tWARNING: side id error checks turned off\n");
    }
    if (prob->groups != nullptr) {
      printf("\telement block groups defined\n");
      printf("\tgroup string: \"%s\"\n", prob->groups);
    }
    return;
  }
} // namespace
