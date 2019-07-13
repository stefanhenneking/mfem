// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443211. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability see http://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// #define DEFINITE

#ifndef MFEM_USE_PETSC
#error This example requires that MFEM is built with MFEM_USE_PETSC=YES
#endif

int main(int argc, char *argv[])
{
   StopWatch chrono;

   // 1. Initialise MPI
   int num_procs, myid;
   MPI_Init(&argc, &argv); // Initialise MPI
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs); //total number of processors available
   MPI_Comm_rank(MPI_COMM_WORLD, &myid); // Determine process identifier

   const char *mesh_file = "ellipse_origin_h0pt0625_o3.mesh";
   int order = 1;
   int ser_ref_levels = 0; // number of serial refinements
   double freq = 1.0e6; //frequecy in herz
   const char *petscrc_file = "petscrc_mult_options";


   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh","Mesh file to use.");
   args.AddOption(&ser_ref_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree)."); 
   args.AddOption(&freq, "-f", "--frequency",
                  "Frequency in Hertz (of course...)");        


   args.Parse();
   // check if the inputs are correct
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      MPI_Finalize();
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }


// 2b. We initialize PETSc
   MFEMInitializePetsc(NULL, NULL, petscrc_file, NULL);


   MFEMFinalizePetsc();
   MPI_Finalize();

}