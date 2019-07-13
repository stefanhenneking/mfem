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
#include <cmath>
#include <complex>

using namespace std;
using namespace mfem;

// #define DEFINITE

#ifndef MFEM_USE_PETSC
#error This example requires that MFEM is built with MFEM_USE_PETSC=YES
#endif

// Define exact solution
void E_exact_Re(const Vector & x, Vector & E);
void f_exact_Re(const Vector & x, Vector & f);
void get_maxwell_solution_Re(const Vector & x, double E[], double curl2E[]);

void E_exact_Im(const Vector & x, Vector & E);
void f_exact_Im(const Vector & x, Vector & f);
void get_maxwell_solution_Im(const Vector & x, double E[], double curl2E[]);


// Mesh Size
Vector mesh_dim_(0); // x, y, z dimensions of mesh
int dim;
double omega;
double complex_shift;
int isol = 1;


int main(int argc, char *argv[])
{
   StopWatch chrono;

   // 1. Initialise MPI
   int num_procs, myid;
   MPI_Init(&argc, &argv); // Initialise MPI
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs); //total number of processors available
   MPI_Comm_rank(MPI_COMM_WORLD, &myid); // Determine process identifier

   int order = 1;
   double freq = 1.0e9; //frequecy in hertz
   // number of wavelengths
   double k = 0.5;
   //
   const char *petscrc_file = "petscrc_mult_options";
   int num_elements_x = 8;
   int num_elements_y = 8;
   int num_elements_z = 8;
   // visualization flag
   bool visualization = 1;
   // number of initial ref
   int initref = 1;
   complex_shift = 0.0;


   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree)."); 
   args.AddOption(&freq, "-f", "--frequency",
                  "Frequency in Hertz (of course...)");        
   args.AddOption(&num_elements_x, "-nex", "--num-elements",
                  "The number of mesh elements in x");
   args.AddOption(&num_elements_y, "-ney", "--num-elements",
                  "The number of mesh elements in y");
   args.AddOption(&num_elements_z, "-nez", "--num-elements",
                  "The number of mesh elements in z");        
   args.AddOption(&k, "-k", "--wavelengths",
                  "Number of wavelengths"); 
   args.AddOption(&complex_shift, "-cs", "--complex_shift",
                  "Complex shift");                
   args.AddOption(&isol, "-isol", "--exact",
                  "Exact solution flag - 0:polynomial, 1: plane wave");
   args.AddOption(&initref, "-initref", "--initref",
                  "Number of initial refinements.");
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

   mesh_dim_.SetSize(3);
   mesh_dim_[0] = 1.0;
   mesh_dim_[1] = 1.0;
   mesh_dim_[2] = 1.0;

   // Angular frequency
   omega = 2.0*k*M_PI;

// 2b. We initialize PETSc
   MFEMInitializePetsc(NULL, NULL, petscrc_file, NULL);

   // Create serial mesh 
   Mesh * mesh = new Mesh(num_elements_x, num_elements_y, num_elements_z, Element::HEXAHEDRON, 1,
                          mesh_dim_(0), mesh_dim_(1), mesh_dim_(2));
   dim = mesh->Dimension();

   // 3. Executing uniform h-refinement
   for (int i = 0; i < initref; i++ )
   {
      mesh->UniformRefinement();
   }

   // create parallel mesh and delete the serial one
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;


   // Create H(curl) (Nedelec) Finite element space
   FiniteElementCollection *fec   = new ND_FECollection(order, dim);
   ParFiniteElementSpace *ND_fespace = new ParFiniteElementSpace(pmesh, fec);

   // Offsets for block vectors (real and imag part)
   Array<int> block_offsets(3);
   block_offsets[0] = 0;
   block_offsets[1] = ND_fespace->GetVSize();
   block_offsets[2] = ND_fespace->GetVSize();
   block_offsets.PartialSum();

   Array<int> block_trueOffsets(3);
   block_trueOffsets[0] = 0;
   block_trueOffsets[1] = ND_fespace->TrueVSize();
   block_trueOffsets[2] = ND_fespace->TrueVSize();
   block_trueOffsets.PartialSum();

   // Solution vector and right-hand side
   BlockVector x(block_offsets), rhs(block_offsets);
   BlockVector trueX(block_trueOffsets), trueRhs(block_trueOffsets);

   x = 0.0;  rhs = 0.0;
   trueX = 0.0;  trueRhs = 0.0;



   Array<int> ess_tdof_list;
   Array<int> ess_bdr(pmesh->bdr_attributes.Max());
   ess_bdr = 1;
   ND_fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);



   // Solution grid functions
   VectorFunctionCoefficient E_Re(dim, E_exact_Re);
   VectorFunctionCoefficient E_Im(dim, E_exact_Im);


   ParComplexGridFunction E_gf(ND_fespace);

   E_gf.real().MakeRef(ND_fespace, x.GetBlock(0));
   E_gf.imag().MakeRef(ND_fespace, x.GetBlock(1));
   E_gf.ProjectCoefficient(E_Re,E_Im);

   // 7. Linear form b(.) (Right hand side)
   ParLinearForm b_Re(ND_fespace);
   VectorFunctionCoefficient f_Re(dim, f_exact_Re);
   b_Re.AddDomainIntegrator(new VectorFEDomainLFIntegrator(f_Re));
   b_Re.Assemble();


   // 7. Linear form b(.) (Right hand side)
   ParLinearForm b_Im(ND_fespace);
   VectorFunctionCoefficient f_Im(dim, f_exact_Im);
   b_Im.AddDomainIntegrator(new VectorFEDomainLFIntegrator(f_Im));
   b_Im.Assemble();


   // 7. Bilinear form a(.,.) on the finite element space
   ConstantCoefficient muinv(1.0);
   ConstantCoefficient sigma(-pow(omega, 2));
   ParBilinearForm a_Re(ND_fespace);
   a_Re.AddDomainIntegrator(new CurlCurlIntegrator(muinv)); 
   a_Re.AddDomainIntegrator(new VectorFEMassIntegrator(sigma));
   a_Re.Assemble();
   HypreParMatrix A_Re;

   a_Re.FormLinearSystem(ess_tdof_list, x.GetBlock(0), b_Re, A_Re, trueX.GetBlock(0), trueRhs.GetBlock(0));
   a_Re.FormLinearSystem(ess_tdof_list, x.GetBlock(1), b_Im, A_Re, trueX.GetBlock(1), trueRhs.GetBlock(1));



   // zero imaginary part just for testing
   ConstantCoefficient zero(0.0);
   ParBilinearForm a_Im(ND_fespace);
   ConstantCoefficient alpha(complex_shift);
   a_Im.AddDomainIntegrator(new VectorFEMassIntegrator(alpha));
   a_Im.Assemble();
   a_Im.EliminateEssentialBC(ess_tdof_list,mfem::Matrix::DIAG_ZERO);
   a_Im.Finalize();
   HypreParMatrix *A_Im = a_Im.ParallelAssemble();

   ComplexHypreParMatrix chpm(&A_Re, A_Im, false, false);
   HypreParMatrix *A = chpm.GetSystemMatrix();

   E_gf.real()=0.0;
   E_gf.imag()=0.0;

   PetscLinearSolver * invA = new PetscLinearSolver(MPI_COMM_WORLD, "direct");
   invA->SetOperator(PetscParMatrix(A, Operator::PETSC_MATAIJ));
   trueX = 0.0;
   invA->Mult(trueRhs,trueX);

   a_Re.RecoverFEMSolution(trueX.GetBlock(0), rhs.GetBlock(0), E_gf.real());
   a_Re.RecoverFEMSolution(trueX.GetBlock(1), rhs.GetBlock(1), E_gf.imag());

   // Compute error
   int order_quad = max(2, 2 * order + 1);
   const IntegrationRule *irs[Geometry::NumGeom];
   for (int i = 0; i < Geometry::NumGeom; ++i)
   {
      irs[i] = &(IntRules.Get(i, order_quad));
   }

   double L2Error_Re = E_gf.real().ComputeL2Error(E_Re, irs);
   double norm_E_Re = ComputeGlobalLpNorm(2, E_Re, *pmesh, irs);

   double L2Error_Im = E_gf.imag().ComputeL2Error(E_Im, irs);
   double norm_E_Im = ComputeGlobalLpNorm(2, E_Im, *pmesh, irs);


   if (myid == 0)
   {
      cout << " Real Part: || E_h - E || / ||E|| = " << L2Error_Re / norm_E_Re << '\n' << endl;
      cout << " Imag Part: || E_h - E || / ||E|| = " << L2Error_Im / norm_E_Im << '\n' << endl;

      cout << " Real Part: || E_h - E || = " << L2Error_Re << '\n' << endl;
      cout << " Imag Part: || E_h - E || = " << L2Error_Im << '\n' << endl;
   }


   // visualization   
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << E_gf.real() << "window_title 'Imaginary part'" << flush;

      socketstream sol_sock_Im(vishost, visport);
      sol_sock_Im << "parallel " << num_procs << " " << myid << "\n";
      sol_sock_Im.precision(8);
      sol_sock_Im << "solution\n" << *pmesh << E_gf.imag() << "window_title 'Imaginary part'" << flush;
   }

   delete invA;
   delete fec;
   delete ND_fespace;
   delete pmesh;

   MFEMFinalizePetsc();
   MPI_Finalize();

}
//define exact solution
void E_exact_Re(const Vector &x, Vector &E)
{
   double curl2E[3];
   get_maxwell_solution_Re(x, E, curl2E);
}

//calculate RHS from exact solution
void f_exact_Re(const Vector &x, Vector &f)
{
   double E_Re[3], curl2E_Re[3];
   double E_Im[3], curl2E_Im[3];

   get_maxwell_solution_Re(x, E_Re, curl2E_Re);
   get_maxwell_solution_Im(x, E_Im, curl2E_Im);

   // curl ( curl E) - omega^2 E = f
   double coeff;
   coeff = -omega * omega;
   f(0) = curl2E_Re[0] + coeff * E_Re[0];
   f(1) = curl2E_Re[1] + coeff * E_Re[1];
   f(2) = curl2E_Re[2] + coeff * E_Re[2];

   // Acount for the complex shift
   f(0) += -complex_shift*E_Im[0];
   f(1) += -complex_shift*E_Im[1];
   f(2) += -complex_shift*E_Im[2];
}

void get_maxwell_solution_Re(const Vector & x, double E[], double curl2E[])
{

   if (isol == 0) // polynomial
   {
      E[0] = x[1] * x[2]      * (1.0 - x[1]) * (1.0 - x[2]);
      E[1] = x[0] * x[1] * x[2] * (1.0 - x[0]) * (1.0 - x[2]);
      E[2] = x[0] * x[1]      * (1.0 - x[0]) * (1.0 - x[1]);
      curl2E[0] = 2.0 * x[1] * (1.0 - x[1]) - (2.0 * x[0] - 3.0) * x[2] * (1 - x[2]);
      curl2E[1] = 2.0 * x[1] * (x[0] * (1.0 - x[0]) + (1.0 - x[2]) * x[2]);
      curl2E[2] = 2.0 * x[1] * (1.0 - x[1]) + x[0] * (3.0 - 2.0 * x[2]) * (1.0 - x[0]);
   }
   else
   {
      double alpha = omega / sqrt(3);
      E[0] = cos(alpha*(x(0) + x(1) + x(2)));
      E[1] = 0.0;
      E[2] = 0.0;

      curl2E[0] = 2.0 * alpha * alpha * E[0];
      curl2E[1] = -alpha * alpha * E[0];
      curl2E[2] = -alpha * alpha * E[0];
   }
}



//define exact solution
void E_exact_Im(const Vector &x, Vector &E)
{
   double curl2E[3];
   get_maxwell_solution_Im(x, E, curl2E);
}

//calculate RHS from exact solution
void f_exact_Im(const Vector &x, Vector &f)
{
   double E_Re[3], curl2E_Re[3];
   double E_Im[3], curl2E_Im[3];

   get_maxwell_solution_Im(x, E_Im, curl2E_Im);
   get_maxwell_solution_Re(x, E_Re, curl2E_Re);

   // curl ( curl E) - omega^2 E = f
   double coeff;
   coeff = -omega * omega;
   f(0) = curl2E_Im[0] + coeff * E_Im[0];
   f(1) = curl2E_Im[1] + coeff * E_Im[1];
   f(2) = curl2E_Im[2] + coeff * E_Im[2];

   // Acount for the complex shift
   f(0) += complex_shift*E_Re[0];
   f(1) += complex_shift*E_Re[1];
   f(2) += complex_shift*E_Re[2];
}

void get_maxwell_solution_Im(const Vector & x, double E[], double curl2E[])
{
   if (isol == 0) // polynomial
   {
      E[0] = x[1] * x[2]      * (1.0 - x[1]) * (1.0 - x[2]);
      E[1] = x[0] * x[1] * x[2] * (1.0 - x[0]) * (1.0 - x[2]);
      E[2] = x[0] * x[1]      * (1.0 - x[0]) * (1.0 - x[1]);
      curl2E[0] = 2.0 * x[1] * (1.0 - x[1]) - (2.0 * x[0] - 3.0) * x[2] * (1 - x[2]);
      curl2E[1] = 2.0 * x[1] * (x[0] * (1.0 - x[0]) + (1.0 - x[2]) * x[2]);
      curl2E[2] = 2.0 * x[1] * (1.0 - x[1]) + x[0] * (3.0 - 2.0 * x[2]) * (1.0 - x[0]);
   }
   else
   {
      double alpha = omega / sqrt(3);   
      E[0] = sin(alpha * (x(0) + x(1) + x(2)));
      E[1] = 0.0;
      E[2] = 0.0;
      curl2E[0] = 2.0 * alpha * alpha * E[0];
      curl2E[1] = -alpha * alpha * E[0];
      curl2E[2] = -alpha * alpha * E[0];
   }   
}