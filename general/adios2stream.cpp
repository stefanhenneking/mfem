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

#include <adios2.h>

namespace mfem
{

namespace{

/**
    * convert openmode input to adios2::Mode format for adios2openmode placeholder
    * @param mode input
    * @return adios2::Mode format
    */
adios2::Mode ToADIOS2Mode(const adios2stream::openmode mode)
{
   adios2::Mode adios2Mode = adios2::Mode::Undefined;
   switch (mode)
   {
      case adios2stream::openmode::out : adios2Mode = adios2::Mode::Write; break;
      case adios2stream::openmode::in : adios2Mode = adios2::Mode::Read; break;
      default:
         throw std::invalid_argument("ERROR: invalid adios2stream, "
                                     "only openmode::out and "
                                     "openmode::in are valid, "
                                     "in call to adios2stream constructor\n");
   }
   return adios2Mode;
}

} //end empty namespace


//PUBLIC
#ifdef MFEM_USE_MPI
adios2stream::adios2stream(const std::string &name, const openmode mode,
                           MPI_Comm comm, const std::string engineType):
   name(name), adios2_openmode(mode),
   adios(std::make_shared<adios2::ADIOS>(comm))
{
   *io = adios->DeclareIO(name);
   io->SetEngine(engineType);
}
#else
adios2stream::adios2stream(const std::string &name, const openmode mode,
                           const std::string engineType) :
   name(name), adios2_openmode(mode),
   adios(std::make_shared<adios2::ADIOS>())
{
  *io = adios->DeclareIO(name);
   io->SetEngine(engineType);
}
#endif


void adios2stream::SetParameters(const std::map<std::string,std::string>& parameters)
{
   io->SetParameters(parameters);
}

void adios2stream::SetParameter(const std::string key,
                                const std::string value)
{
   io->SetParameter(key, value);
}

} //end namespace mfem
