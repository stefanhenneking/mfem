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
//
// Created on: Jan 22, 2019
// Author: William F Godoy godoywf@ornl.gov
// adios2: Adaptable Input/Output System https://github.com/ornladios/ADIOS2

#ifndef MFEM_ADIOS2STREAM
#define MFEM_ADIOS2STREAM

#include <adios2.h>

#include <string>
#include <memory> // std::shared_ptr
#include <array>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace mfem
{

//forward declaring friends
class Vector;
class FiniteElementSpace;
class GridFunction;
class Mesh;

#ifdef MFEM_USE_MPI
class ParGridFunction;
class ParMesh;
#endif


class adios2stream
{
   friend class Vector;
   friend class FiniteElementSpace;
   friend class GridFunction;
   friend class Mesh;

#ifdef MFEM_USE_MPI
   friend class ParGridFunction;
   friend class ParMesh;
#endif


public:

   /**
    * Open modes for adios2stream (from fstream)
    * out: write
    * in:  read
    * app: append
    */
   enum class openmode {out,in,app};


#ifdef MFEM_USE_MPI

   /**
     * adios2stream MPI constructor, allows for passing parameters in source
     * code (compile-time) only.
     * @param name stream name
     * @param mode adios2stream::openmode::in (Read), adios2stream::openmode::out (Write)
     * @param comm MPI communicator establishing domain for fstream
     * @param engineType available adios2 engine, default is BP3 file
     *        see https://adios2.readthedocs.io/en/latest/engines/engines.html
     * @throws std::invalid_argument (user input error) or std::runtime_error
     *         (system error)
    */
   adios2stream(const std::string &name, const openmode mode,
                MPI_Comm comm, const std::string engineType = "BPFile");

   /**
    * adios2stream MPI constructor, allows for passing parameters in
    * a configuration file (run-time).
    * @param name stream name
    * @param mode adios2stream::openmode::in (Read), adios2stream::openmode::out (Write)
    * @param comm MPI communicator establishing domain for fstream
    * @param configFile adios2 xml runtime configuration (must have .xml extension)
    *        see https://adios2.readthedocs.io/en/latest/components/components.html#runtime-configuration-files
    * @param ioInConfigFile identifies an io field in the configFile to relate for this stream parameters
    * @throws std::invalid_argument (user input error) or std::runtime_error
    *         (system error)
    */
   adios2stream(const std::string &name, const openmode mode,
                MPI_Comm comm, const std::string &configFile,
                const std::string ioInConfigFile);

#else
   /**
    * adios2stream Non-MPI serial constructor, allows for passing parameters in source
    * code (compile-time) only.
    * @param name stream name
    * @param mode adios2stream::openmode::in (Read), adios2stream::openmode::out (Write)
    * @param engineType available adios2 engine, default is BP3 file
    * @throws std::invalid_argument (user input error) or std::runtime_error
    *         (system error)
    */
   adios2stream(const std::string &name, const openmode mode,
                const std::string engineType = "BPFile");

   /**
    * adios2stream Non-MPI serial constructor, allows for passing parameters in
    * a configuration file (run-time).
    * @param name stream name
    * @param mode adios2stream::openmode::in (Read), adios2stream::openmode::out (Write)
    * @param configFile adios2 xml runtime configuration (must have .xml extension)
    *        see https://adios2.readthedocs.io/en/latest/components/components.html#runtime-configuration-files
    * @param ioInConfigFile identifies an io field in the configFile to relate for this stream parameters
    * @throws std::invalid_argument (user input error) or std::runtime_error
    *         (system error)
    */
   adios2stream(const std::string &name, const openmode mode,
                const std::string &configFile,
                const std::string ioInConfigFile);

#endif

   /** using RAII components, nothing to be deallocated **/
   virtual ~adios2stream() = default;

   /**
    * Set parameters for a particular adios2stream Engine
    * See https://adios2.readthedocs.io/en/latest/engines/engines.html#bp3-default
    * @param parameters map of key/value string elements
    */
   void SetParameters(const adios2::Params& parameters = adios2::Params() );

   /**
    * Single parameter version of SetParameters passing a key/value pair
    * See https://adios2.readthedocs.io/en/latest/engines/engines.html#bp3-default
    * @param key input parameter key
    * @param value input parameter value
    */
   void SetParameter(const std::string key, const std::string value);

private:

   /** placeholder for engine name */
   const std::string name;

   /** placeholder for engine openmode */
   const adios2::Mode adios2openmode;

   /** main adios2 object that owns all the io and engine components */
   std::shared_ptr<adios2::ADIOS> adios;

   /** io object to set parameters, variables and engines */
   adios2::IO io;

   /** heavy object doing system-level I/O operations */
   adios2::Engine engine;

   static const std::array<uint32_t,7> glvis2vtk;
   static const std::vector<std::string> glvis_types;

   /**
    * convert openmode input to adios2::Mode format for adios2openmode placeholder
    * @param mode input
    * @return adios2::Mode format
    */
   adios2::Mode ToADIOS2Mode(const openmode mode);
};


} //end namespace mfem

#endif /* MFEM_ADIOS2STREAM */
