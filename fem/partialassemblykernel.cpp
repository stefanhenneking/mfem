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

#include "fem.hpp"

namespace mfem
{

void MultBtDB1(FiniteElementSpace* fes, DenseMatrix const & shape1d,
   DummyTensor & D, const Vector &V, Vector &U)
{
   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads = quads1d;

   Vector Q(quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const Vector Vmat(V.GetData() + offset, dofs1d);
      Vector Umat(U.GetData() + offset, dofs1d);

      // Q_k1 = dshape_j1_k1 * V_i1
      shape1d.MultTranspose(Vmat, Q);

      double *data_q = Q.GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k) { data_q[k] *= data_d[k]; }

      // Q_k1 = dshape_j1_k1 * Q_k1
      shape1d.AddMult(Q, Umat);
   }
}

void MultBtDB2(FiniteElementSpace* fes, DenseMatrix const & shape1d,
   DummyTensor & D, const Vector &V, Vector &U)
{
   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d*quads1d;

   DenseMatrix QQ(quads1d, quads1d);
   DenseMatrix DQ(dofs1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseMatrix Vmat(V.GetData() + offset, dofs1d, dofs1d);
      DenseMatrix Umat(U.GetData() + offset, dofs1d, dofs1d);

      // DQ_j2_k1   = E_j1_j2  * dshape_j1_k1 -- contract in x direction
      // QQ_0_k1_k2 = DQ_j2_k1 * shape_j2_k2  -- contract in y direction
      MultAtB(Vmat, shape1d, DQ);
      MultAtB(DQ, shape1d, QQ);

      // QQ_c_k1_k2 = Dmat_c_d_k1_k2 * QQ_d_k1_k2
      // NOTE: (k1, k2) = k -- 1d index over tensor product of quad points
      double *data_qq = QQ.GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k) { data_qq[k] *= data_d[k]; }

      // DQ_i2_k1   = shape_i2_k2  * QQ_0_k1_k2
      // U_i1_i2   += dshape_i1_k1 * DQ_i2_k1
      MultABt(shape1d, QQ, DQ);
      AddMultABt(shape1d, DQ, Umat);

      // increment offset
      offset += dofs;
   }
}

void MultBtDB3(FiniteElementSpace* fes, DenseMatrix const & shape1d,
   DummyTensor & D, const Vector &V, Vector &U)
{
   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d * quads1d;

   Vector Q(quads1d);
   DenseMatrix QQ(quads1d, quads1d);
   DenseTensor QQQ(quads1d, quads1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseTensor Vmat(V.GetData() + offset, dofs1d, dofs1d, dofs1d);
      DenseTensor Umat(U.GetData() + offset, dofs1d, dofs1d, dofs1d);

      // QQQ_k1_k2_k3 = shape_j1_k1 * shape_j2_k2  * shape_j3_k3  * Vmat_j1_j2_j3
      QQQ = 0.;
      for (int j3 = 0; j3 < dofs1d; ++j3)
      {
         QQ = 0.;
         for (int j2 = 0; j2 < dofs1d; ++j2)
         {
            Q = 0.;
            for (int j1 = 0; j1 < dofs1d; ++j1)
            {
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  Q(k1) += Vmat(j1, j2, j3) * shape1d(j1, k1);
               }
            }
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQ(k1, k2) += Q(k1) * shape1d(j2, k2);
               }
         }
         for (int k3 = 0; k3 < quads1d; ++k3)
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQQ(k1, k2, k3) += QQ(k1, k2) * shape1d(j3, k3);
               }
      }

      // QQQ_k1_k2_k3 = Dmat_k1_k2_k3 * QQQ_k1_k2_k3
      // NOTE: (k1, k2, k3) = q -- 1d quad point index
      double *data_qqq = QQQ.GetData(0);
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k) { data_qqq[k] *= data_d[k]; }

      // Apply transpose of the first operator that takes V -> QQQ -- QQQ -> U
      for (int k3 = 0; k3 < quads1d; ++k3)
      {
         QQ = 0.;
         for (int k2 = 0; k2 < quads1d; ++k2)
         {
            Q = 0.;
            for (int k1 = 0; k1 < quads1d; ++k1)
            {
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Q(i1) += QQQ(k1, k2, k3) * shape1d(i1, k1);
               }
            }
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  QQ(i1, i2) += Q(i1) * shape1d(i2, k2);
               }
         }
         for (int i3 = 0; i3 < dofs1d; ++i3)
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Umat(i1, i2, i3) += shape1d(i3, k3) * QQ(i1, i2);
               }
      }

      // increment offset
      offset += dofs;
   }  
}

void MultGtDG1(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads = quads1d;

   Vector Q(quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const Vector Vmat(V.GetData() + offset, dofs1d);
      Vector Umat(U.GetData() + offset, dofs1d);

      // Q_k1 = dshape_j1_k1 * V_i1
      dshape1d.MultTranspose(Vmat, Q);

      double *data_q = Q.GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         data_q[k] *= data_d[k];
      }

      // Q_k1 = dshape_j1_k1 * Q_k1
      dshape1d.AddMult(Q, Umat);
   }   
}

void MultGtDG2(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dim = 2;
   const int terms = dim*dim;

   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d;

   DenseTensor QQ(quads1d, quads1d, dim);
   DenseMatrix DQ(dofs1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseMatrix Vmat(V.GetData() + offset, dofs1d, dofs1d);
      DenseMatrix Umat(U.GetData() + offset, dofs1d, dofs1d);

      // DQ_j2_k1   = E_j1_j2  * dshape_j1_k1 -- contract in x direction
      // QQ_0_k1_k2 = DQ_j2_k1 * shape_j2_k2  -- contract in y direction
      MultAtB(Vmat, dshape1d, DQ);
      MultAtB(DQ, shape1d, QQ(0));

      // DQ_j2_k1   = E_j1_j2  * shape_j1_k1  -- contract in x direction
      // QQ_1_k1_k2 = DQ_j2_k1 * dshape_j2_k2 -- contract in y direction
      MultAtB(Vmat, shape1d, DQ);
      MultAtB(DQ, dshape1d, QQ(1));

      // QQ_c_k1_k2 = Dmat_c_d_k1_k2 * QQ_d_k1_k2
      // NOTE: (k1, k2) = k -- 1d index over tensor product of quad points
      double *data_qq = QQ(0).GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         const double D00 = data_d[terms*k + 0];
         const double D10 = data_d[terms*k + 1];
         const double D01 = data_d[terms*k + 2];
         const double D11 = data_d[terms*k + 3];

         const double q0 = data_qq[0*quads + k];
         const double q1 = data_qq[1*quads + k];

         data_qq[0*quads + k] = D00 * q0 + D01 * q1;
         data_qq[1*quads + k] = D10 * q0 + D11 * q1;
      }

      // DQ_i2_k1   = shape_i2_k2  * QQ_0_k1_k2
      // U_i1_i2   += dshape_i1_k1 * DQ_i2_k1
      MultABt(shape1d, QQ(0), DQ);
      AddMultABt(dshape1d, DQ, Umat);

      // DQ_i2_k1   = dshape_i2_k2 * QQ_1_k1_k2
      // U_i1_i2   += shape_i1_k1  * DQ_i2_k1
      MultABt(dshape1d, QQ(1), DQ);
      AddMultABt(shape1d, DQ, Umat);

      // increment offset
      offset += dofs;
   }   
}

void MultGtDG3(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dim = 3;
   const int terms = dim*dim;

   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d * quads1d;

   DenseMatrix Q(quads1d, dim);
   DenseTensor QQ(quads1d, quads1d, dim);

   Array<double> QQQmem(quads1d * quads1d * quads1d * dim);
   double *data_qqq = QQQmem.GetData();
   DenseTensor QQQ0(data_qqq + 0*quads, quads1d, quads1d, quads1d);
   DenseTensor QQQ1(data_qqq + 1*quads, quads1d, quads1d, quads1d);
   DenseTensor QQQ2(data_qqq + 2*quads, quads1d, quads1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseTensor Vmat(V.GetData() + offset, dofs1d, dofs1d, dofs1d);
      DenseTensor Umat(U.GetData() + offset, dofs1d, dofs1d, dofs1d);

      // QQQ_0_k1_k2_k3 = dshape_j1_k1 * shape_j2_k2  * shape_j3_k3  * Vmat_j1_j2_j3
      // QQQ_1_k1_k2_k3 = shape_j1_k1  * dshape_j2_k2 * shape_j3_k3  * Vmat_j1_j2_j3
      // QQQ_2_k1_k2_k3 = shape_j1_k1  * shape_j2_k2  * dshape_j3_k3 * Vmat_j1_j2_j3
      QQQ0 = 0.; QQQ1 = 0.; QQQ2 = 0.;
      for (int j3 = 0; j3 < dofs1d; ++j3)
      {
         QQ = 0.;
         for (int j2 = 0; j2 < dofs1d; ++j2)
         {
            Q = 0.;
            for (int j1 = 0; j1 < dofs1d; ++j1)
            {
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  Q(k1, 0) += Vmat(j1, j2, j3) * dshape1d(j1, k1);
                  Q(k1, 1) += Vmat(j1, j2, j3) * shape1d(j1, k1);
               }
            }
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQ(k1, k2, 0) += Q(k1, 0) * shape1d(j2, k2);
                  QQ(k1, k2, 1) += Q(k1, 1) * dshape1d(j2, k2);
                  QQ(k1, k2, 2) += Q(k1, 1) * shape1d(j2, k2);
               }
         }
         for (int k3 = 0; k3 < quads1d; ++k3)
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQQ0(k1, k2, k3) += QQ(k1, k2, 0) * shape1d(j3, k3);
                  QQQ1(k1, k2, k3) += QQ(k1, k2, 1) * shape1d(j3, k3);
                  QQQ2(k1, k2, k3) += QQ(k1, k2, 2) * dshape1d(j3, k3);
               }
      }

      // QQQ_c_k1_k2_k3 = Dmat_c_d_k1_k2_k3 * QQQ_d_k1_k2_k3
      // NOTE: (k1, k2, k3) = q -- 1d quad point index
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         const double D00 = data_d[terms*k + 0];
         const double D10 = data_d[terms*k + 1];
         const double D20 = data_d[terms*k + 2];
         const double D01 = data_d[terms*k + 3];
         const double D11 = data_d[terms*k + 4];
         const double D21 = data_d[terms*k + 5];
         const double D02 = data_d[terms*k + 6];
         const double D12 = data_d[terms*k + 7];
         const double D22 = data_d[terms*k + 8];

         const double q0 = data_qqq[0*quads + k];
         const double q1 = data_qqq[1*quads + k];
         const double q2 = data_qqq[2*quads + k];

         data_qqq[0*quads + k] = D00 * q0 + D01 * q1 + D02 * q2;
         data_qqq[1*quads + k] = D10 * q0 + D11 * q1 + D12 * q2;
         data_qqq[2*quads + k] = D20 * q0 + D21 * q1 + D22 * q2;
      }

      // Apply transpose of the first operator that takes V -> QQQd -- QQQd -> U
      for (int k3 = 0; k3 < quads1d; ++k3)
      {
         QQ = 0.;
         for (int k2 = 0; k2 < quads1d; ++k2)
         {
            Q = 0.;
            for (int k1 = 0; k1 < quads1d; ++k1)
            {
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Q(i1, 0) += QQQ0(k1, k2, k3) * dshape1d(i1, k1);
                  Q(i1, 1) += QQQ1(k1, k2, k3) * shape1d(i1, k1);
                  Q(i1, 2) += QQQ2(k1, k2, k3) * shape1d(i1, k1);
               }
            }
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  QQ(i1, i2, 0) += Q(i1, 0) * shape1d(i2, k2);
                  QQ(i1, i2, 1) += Q(i1, 1) * dshape1d(i2, k2);
                  QQ(i1, i2, 2) += Q(i1, 2) * shape1d(i2, k2);
               }
         }
         for (int i3 = 0; i3 < dofs1d; ++i3)
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Umat(i1, i2, i3) +=
                     QQ(i1, i2, 0) * shape1d(i3, k3) +
                     QQ(i1, i2, 1) * shape1d(i3, k3) +
                     QQ(i1, i2, 2) * dshape1d(i3, k3);
               }
      }

      // increment offset
      offset += dofs;
   }
}

void MultBtDG1(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads = quads1d;

   Vector Q(quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const Vector Vmat(V.GetData() + offset, dofs1d);
      Vector Umat(U.GetData() + offset, dofs1d);

      // Q_k1 = dshape_j1_k1 * V_i1
      dshape1d.MultTranspose(Vmat, Q);

      double *data_q = Q.GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         data_q[k] *= data_d[k];
      }

      // Q_k1 = dshape_j1_k1 * Q_k1
      shape1d.AddMult(Q, Umat);
   }   
}

void MultBtDG2(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dim = 2;
   const int terms = dim*dim;

   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d;

   DenseTensor QQ(quads1d, quads1d, dim);
   DenseMatrix DQ(dofs1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseMatrix Vmat(V.GetData() + offset, dofs1d, dofs1d);
      DenseMatrix Umat(U.GetData() + offset, dofs1d, dofs1d);

      // DQ_j2_k1   = E_j1_j2  * dshape_j1_k1 -- contract in x direction
      // QQ_0_k1_k2 = DQ_j2_k1 * shape_j2_k2  -- contract in y direction
      MultAtB(Vmat, dshape1d, DQ);
      MultAtB(DQ, shape1d, QQ(0));

      // DQ_j2_k1   = E_j1_j2  * shape_j1_k1  -- contract in x direction
      // QQ_1_k1_k2 = DQ_j2_k1 * dshape_j2_k2 -- contract in y direction
      MultAtB(Vmat, shape1d, DQ);
      MultAtB(DQ, dshape1d, QQ(1));

      // QQ_c_k1_k2 = Dmat_c_d_k1_k2 * QQ_d_k1_k2
      // NOTE: (k1, k2) = k -- 1d index over tensor product of quad points
      double *data_qq = QQ(0).GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         const double D00 = data_d[terms*k + 0];
         const double D10 = data_d[terms*k + 1];
         const double D01 = data_d[terms*k + 2];
         const double D11 = data_d[terms*k + 3];

         const double q0 = data_qq[0*quads + k];
         const double q1 = data_qq[1*quads + k];

         data_qq[0*quads + k] = D00 * q0 + D01 * q1;
         data_qq[1*quads + k] = D10 * q0 + D11 * q1;
      }

      // DQ_i2_k1   = shape_i2_k2  * QQ_0_k1_k2
      // U_i1_i2   += dshape_i1_k1 * DQ_i2_k1
      MultABt(shape1d, QQ(0), DQ);
      AddMultABt(shape1d, DQ, Umat);

      // DQ_i2_k1   = dshape_i2_k2 * QQ_1_k1_k2
      // U_i1_i2   += shape_i1_k1  * DQ_i2_k1
      MultABt(shape1d, QQ(1), DQ);
      AddMultABt(shape1d, DQ, Umat);

      // increment offset
      offset += dofs;
   }   
}

void MultBtDG3(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dim = 3;
   const int terms = dim*dim;

   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d * quads1d;

   DenseMatrix Q(quads1d, dim);
   DenseTensor QQ(quads1d, quads1d, dim);

   Array<double> QQQmem(quads1d * quads1d * quads1d * dim);
   double *data_qqq = QQQmem.GetData();
   DenseTensor QQQ0(data_qqq + 0*quads, quads1d, quads1d, quads1d);
   DenseTensor QQQ1(data_qqq + 1*quads, quads1d, quads1d, quads1d);
   DenseTensor QQQ2(data_qqq + 2*quads, quads1d, quads1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseTensor Vmat(V.GetData() + offset, dofs1d, dofs1d, dofs1d);
      DenseTensor Umat(U.GetData() + offset, dofs1d, dofs1d, dofs1d);

      // QQQ_0_k1_k2_k3 = dshape_j1_k1 * shape_j2_k2  * shape_j3_k3  * Vmat_j1_j2_j3
      // QQQ_1_k1_k2_k3 = shape_j1_k1  * dshape_j2_k2 * shape_j3_k3  * Vmat_j1_j2_j3
      // QQQ_2_k1_k2_k3 = shape_j1_k1  * shape_j2_k2  * dshape_j3_k3 * Vmat_j1_j2_j3
      QQQ0 = 0.; QQQ1 = 0.; QQQ2 = 0.;
      for (int j3 = 0; j3 < dofs1d; ++j3)
      {
         QQ = 0.;
         for (int j2 = 0; j2 < dofs1d; ++j2)
         {
            Q = 0.;
            for (int j1 = 0; j1 < dofs1d; ++j1)
            {
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  Q(k1, 0) += Vmat(j1, j2, j3) * dshape1d(j1, k1);
                  Q(k1, 1) += Vmat(j1, j2, j3) * shape1d(j1, k1);
               }
            }
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQ(k1, k2, 0) += Q(k1, 0) * shape1d(j2, k2);
                  QQ(k1, k2, 1) += Q(k1, 1) * dshape1d(j2, k2);
                  QQ(k1, k2, 2) += Q(k1, 1) * shape1d(j2, k2);
               }
         }
         for (int k3 = 0; k3 < quads1d; ++k3)
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQQ0(k1, k2, k3) += QQ(k1, k2, 0) * shape1d(j3, k3);
                  QQQ1(k1, k2, k3) += QQ(k1, k2, 1) * shape1d(j3, k3);
                  QQQ2(k1, k2, k3) += QQ(k1, k2, 2) * dshape1d(j3, k3);
               }
      }

      // QQQ_c_k1_k2_k3 = Dmat_c_d_k1_k2_k3 * QQQ_d_k1_k2_k3
      // NOTE: (k1, k2, k3) = q -- 1d quad point index
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         const double D00 = data_d[terms*k + 0];
         const double D10 = data_d[terms*k + 1];
         const double D20 = data_d[terms*k + 2];
         const double D01 = data_d[terms*k + 3];
         const double D11 = data_d[terms*k + 4];
         const double D21 = data_d[terms*k + 5];
         const double D02 = data_d[terms*k + 6];
         const double D12 = data_d[terms*k + 7];
         const double D22 = data_d[terms*k + 8];

         const double q0 = data_qqq[0*quads + k];
         const double q1 = data_qqq[1*quads + k];
         const double q2 = data_qqq[2*quads + k];

         data_qqq[0*quads + k] = D00 * q0 + D01 * q1 + D02 * q2;
         data_qqq[1*quads + k] = D10 * q0 + D11 * q1 + D12 * q2;
         data_qqq[2*quads + k] = D20 * q0 + D21 * q1 + D22 * q2;
      }

      // Apply transpose of the first operator that takes V -> QQQd -- QQQd -> U
      for (int k3 = 0; k3 < quads1d; ++k3)
      {
         QQ = 0.;
         for (int k2 = 0; k2 < quads1d; ++k2)
         {
            Q = 0.;
            for (int k1 = 0; k1 < quads1d; ++k1)
            {
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Q(i1, 0) += QQQ0(k1, k2, k3) * shape1d(i1, k1);
                  Q(i1, 1) += QQQ1(k1, k2, k3) * shape1d(i1, k1);
                  Q(i1, 2) += QQQ2(k1, k2, k3) * shape1d(i1, k1);
               }
            }
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  QQ(i1, i2, 0) += Q(i1, 0) * shape1d(i2, k2);
                  QQ(i1, i2, 1) += Q(i1, 1) * shape1d(i2, k2);
                  QQ(i1, i2, 2) += Q(i1, 2) * shape1d(i2, k2);
               }
         }
         for (int i3 = 0; i3 < dofs1d; ++i3)
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Umat(i1, i2, i3) +=
                     QQ(i1, i2, 0) * shape1d(i3, k3) +
                     QQ(i1, i2, 1) * shape1d(i3, k3) +
                     QQ(i1, i2, 2) * shape1d(i3, k3);
               }
      }

      // increment offset
      offset += dofs;
   }
}

void MultGtDB1(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads = quads1d;

   Vector Q(quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const Vector Vmat(V.GetData() + offset, dofs1d);
      Vector Umat(U.GetData() + offset, dofs1d);

      // Q_k1 = shape_j1_k1 * V_i1
      shape1d.MultTranspose(Vmat, Q);

      double *data_q = Q.GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         data_q[k] *= data_d[k];
      }

      // Q_k1 = dshape_j1_k1 * Q_k1
      dshape1d.AddMult(Q, Umat);
   }    
}

void MultGtDB2(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dim = 2;
   const int terms = dim*dim;

   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d;

   DenseTensor QQ(quads1d, quads1d, dim);
   DenseMatrix DQ(dofs1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseMatrix Vmat(V.GetData() + offset, dofs1d, dofs1d);
      DenseMatrix Umat(U.GetData() + offset, dofs1d, dofs1d);

      //TODO One QQ should be enough
      // DQ_j2_k1   = E_j1_j2  * shape_j1_k1 -- contract in x direction
      // QQ_k1_k2 = DQ_j2_k1 * shape_j2_k2  -- contract in y direction
      MultAtB(Vmat, shape1d, DQ);
      MultAtB(DQ, shape1d, QQ(0));

      // DQ_j2_k1   = E_j1_j2  * shape_j1_k1  -- contract in x direction
      // QQ_1_k1_k2 = DQ_j2_k1 * dshape_j2_k2 -- contract in y direction
      //Can be optimized since this is the same data
      MultAtB(Vmat, shape1d, DQ);
      MultAtB(DQ, shape1d, QQ(1));

      // QQ_c_k1_k2 = Dmat_c_d_k1_k2 * QQ_d_k1_k2
      // NOTE: (k1, k2) = k -- 1d index over tensor product of quad points
      double *data_qq = QQ(0).GetData();
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         const double D00 = data_d[terms*k + 0];
         const double D10 = data_d[terms*k + 1];
         const double D01=  data_d[terms*k + 2];
         const double D11 = data_d[terms*k + 3];

         const double q0  =data_qq[0*quads + k];
         const double q1 = data_qq[1*quads + k];

         data_qq[0*quads + k] = D00 * q0 + D01 * q1;
         data_qq[1*quads + k] = D10 * q0 + D11 * q1;
      }

      // DQ_i2_k1   = shape_i2_k2  * QQ_0_k1_k2
      // U_i1_i2   += dshape_i1_k1 * DQ_i2_k1
      MultABt(shape1d, QQ(0), DQ);
      AddMultABt(dshape1d, DQ, Umat);

      // DQ_i2_k1   = dshape_i2_k2 * QQ_1_k1_k2
      // U_i1_i2   += shape_i1_k1  * DQ_i2_k1
      MultABt(dshape1d, QQ(1), DQ);
      AddMultABt(shape1d, DQ, Umat);

      // increment offset
      offset += dofs;
   }    
}

void MultGtDB3(FiniteElementSpace* fes, DenseMatrix const& shape1d,
   DenseMatrix const& dshape1d, DummyTensor & D, const Vector &V, Vector &U)
{
   const int dim = 3;
   const int terms = dim*(dim+1)/2;

   const FiniteElement *fe = fes->GetFE(0);
   const int dofs   = fe->GetDof();

   const int dofs1d = shape1d.Height();
   const int quads1d = shape1d.Width();
   const int quads  = quads1d * quads1d * quads1d;

   DenseMatrix Q(quads1d, dim);
   DenseTensor QQ(quads1d, quads1d, dim);

   Array<double> QQQmem(quads1d * quads1d * quads1d * dim);
   double *data_qqq = QQQmem.GetData();
   DenseTensor QQQ0(data_qqq + 0*quads, quads1d, quads1d, quads1d);
   DenseTensor QQQ1(data_qqq + 1*quads, quads1d, quads1d, quads1d);
   DenseTensor QQQ2(data_qqq + 2*quads, quads1d, quads1d, quads1d);

   int offset = 0;
   for (int e = 0; e < fes->GetNE(); e++)
   {
      const DenseTensor Vmat(V.GetData() + offset, dofs1d, dofs1d, dofs1d);
      DenseTensor Umat(U.GetData() + offset, dofs1d, dofs1d, dofs1d);

      // TODO One QQQ should be enough
      // QQQ_0_k1_k2_k3 = shape_j1_k1 * shape_j2_k2 * shape_j3_k3 * Vmat_j1_j2_j3
      // QQQ_1_k1_k2_k3 = shape_j1_k1 * shape_j2_k2 * shape_j3_k3 * Vmat_j1_j2_j3
      // QQQ_2_k1_k2_k3 = shape_j1_k1 * shape_j2_k2 * shape_j3_k3 * Vmat_j1_j2_j3
      QQQ0 = 0.; QQQ1 = 0.; QQQ2 = 0.;
      for (int j3 = 0; j3 < dofs1d; ++j3)
      {
         QQ = 0.;
         for (int j2 = 0; j2 < dofs1d; ++j2)
         {
            Q = 0.;
            for (int j1 = 0; j1 < dofs1d; ++j1)
            {
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  Q(k1, 0) += Vmat(j1, j2, j3) * shape1d(j1, k1);
                  Q(k1, 1) += Vmat(j1, j2, j3) * shape1d(j1, k1);
               }
            }
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQ(k1, k2, 0) += Q(k1, 0) * shape1d(j2, k2);
                  QQ(k1, k2, 1) += Q(k1, 1) * shape1d(j2, k2);
                  QQ(k1, k2, 2) += Q(k1, 1) * shape1d(j2, k2);
               }
         }
         for (int k3 = 0; k3 < quads1d; ++k3)
            for (int k2 = 0; k2 < quads1d; ++k2)
               for (int k1 = 0; k1 < quads1d; ++k1)
               {
                  QQQ0(k1, k2, k3) += QQ(k1, k2, 0) * shape1d(j3, k3);
                  QQQ1(k1, k2, k3) += QQ(k1, k2, 1) * shape1d(j3, k3);
                  QQQ2(k1, k2, k3) += QQ(k1, k2, 2) * shape1d(j3, k3);
               }
      }

      //TODO insert the three QQQ only here
      // QQQ_c_k1_k2_k3 = Dmat_c_d_k1_k2_k3 * QQQ_d_k1_k2_k3
      // NOTE: (k1, k2, k3) = q -- 1d quad point index
      const double *data_d = D.GetElmtData(e);
      for (int k = 0; k < quads; ++k)
      {
         const double D00 = data_d[terms*k + 0];
         const double D10 = data_d[terms*k + 1];
         const double D20 = data_d[terms*k + 2];
         const double D01 = data_d[terms*k + 3];
         const double D11 = data_d[terms*k + 4];
         const double D21 = data_d[terms*k + 5];
         const double D02 = data_d[terms*k + 6];
         const double D12 = data_d[terms*k + 7];
         const double D22 = data_d[terms*k + 8];

         const double q0 = data_qqq[0*quads + k];
         const double q1 = data_qqq[1*quads + k];
         const double q2 = data_qqq[2*quads + k];

         data_qqq[0*quads + k] = D00 * q0 + D01 * q1 + D02 * q2;
         data_qqq[1*quads + k] = D10 * q0 + D11 * q1 + D12 * q2;
         data_qqq[2*quads + k] = D20 * q0 + D21 * q1 + D22 * q2;
      }

      // Apply transpose of the first operator that takes V -> QQQd -- QQQd -> U
      for (int k3 = 0; k3 < quads1d; ++k3)
      {
         QQ = 0.;
         for (int k2 = 0; k2 < quads1d; ++k2)
         {
            Q = 0.;
            for (int k1 = 0; k1 < quads1d; ++k1)
            {
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Q(i1, 0) += QQQ0(k1, k2, k3) * dshape1d(i1, k1);
                  Q(i1, 1) += QQQ1(k1, k2, k3) * shape1d(i1, k1);
                  Q(i1, 2) += QQQ2(k1, k2, k3) * shape1d(i1, k1);
               }
            }
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  QQ(i1, i2, 0) += Q(i1, 0) * shape1d(i2, k2);
                  QQ(i1, i2, 1) += Q(i1, 1) * dshape1d(i2, k2);
                  QQ(i1, i2, 2) += Q(i1, 2) * shape1d(i2, k2);
               }
         }
         for (int i3 = 0; i3 < dofs1d; ++i3)
            for (int i2 = 0; i2 < dofs1d; ++i2)
               for (int i1 = 0; i1 < dofs1d; ++i1)
               {
                  Umat(i1, i2, i3) +=
                     QQ(i1, i2, 0) * shape1d(i3, k3) +
                     QQ(i1, i2, 1) * shape1d(i3, k3) +
                     QQ(i1, i2, 2) * dshape1d(i3, k3);
               }
      }

      // increment offset
      offset += dofs;
   }   
}

}