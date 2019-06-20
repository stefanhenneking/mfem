/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2stream.cpp : implementation of adios2stream functions
 *
 *  Created on: Feb 4, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "adios2stream.hpp"

namespace mfem
{

//PUBLIC
#ifdef MFEM_USE_MPI
adios2stream::adios2stream(const std::string &name, const openmode mode,
                           MPI_Comm comm, const std::string engineType):
   name(name), adios2openmode(ToADIOS2Mode(mode)),
   adios(std::make_shared<adios2::ADIOS>(comm)),
   io( adios->DeclareIO(name) )
{
   io.SetEngine(engineType);
}

adios2stream::adios2stream(const std::string &name, const openmode mode,
                           MPI_Comm comm, const std::string &configFile,
                           const std::string ioInConfigFile) :
   name(name), adios2openmode(ToADIOS2Mode(mode)),
   adios(std::make_shared<adios2::ADIOS>(configFile, comm)),
   io( adios->DeclareIO(ioInConfigFile) )
{
   int rank;
   MPI_Comm_rank(comm, &rank);

   if (rank == 0 && !io.InConfigFile())
   {
      std::cout << "WARNING: adios2stream io: " << ioInConfigFile
                << " not found in config file: " << configFile
                << " assuming defaults, "
                << " in call to adios2stream " << name << " constructor\n";
   }
}
#else
adios2stream::adios2stream(const std::string &name, const openmode mode,
                           const std::string engineType) :
   name(name), adios2openmode(ToADIOS2Mode(mode)),
   adios(std::make_shared<adios2::ADIOS>()),
   io( adios->DeclareIO(name) )
{
   io.SetEngine(engineType);
}


adios2stream::adios2stream(const std::string &name, const openmode mode,
                           const std::string &configFile,
                           const std::string ioInConfigFile):
   name(name), adios2openmode(ToADIOS2Mode(mode)),
   adios(std::make_shared<adios2::ADIOS>(configFile)),
   io( adios->DeclareIO(ioInConfigFile) )
{
   if (!io.InConfigFile())
   {
      std::cout << "WARNING: adios2stream io: " << ioInConfigFile
                << " not found in config file: " << configFile
                << " assuming defaults, "
                << " in call to adios2stream " << name << " constructor\n";
   }
}
#endif

adios2stream::~adios2stream()
{
}

void adios2stream::SetParameters(const adios2::Params& parameters)
{
   io.SetParameters(parameters);
}

void adios2stream::SetParameter(const std::string key,
                                const std::string value)
{
   io.SetParameter(key, value);
}

// PRIVATE
adios2::Mode adios2stream::ToADIOS2Mode(const adios2stream::openmode mode)
{
   adios2::Mode adios2Mode = adios2::Mode::Undefined;
   switch (mode)
   {
      case openmode::out : adios2Mode = adios2::Mode::Write; break;
      case openmode::in : adios2Mode = adios2::Mode::Read; break;
      default:
         throw std::invalid_argument("ERROR: invalid adios2stream, "
                                     "only openmode::out and "
                                     "openmode::in are valid, "
                                     "in call to adios2stream constructor\n");
   }
   return adios2Mode;
}


} //end namespace mfem
