/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    EXACT                                                                  */
/*                                                                           */
/* Copyright (C) 2024             Zuse Institute Berlin                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with MIP-DD; see the file LICENSE. If not visit scipopt.org.       */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _EXACT_MISC_VERSION_LOGGER_HPP_
#define _EXACT_MISC_VERSION_LOGGER_HPP_

#include "exact/io/MpsParser.hpp"
#include "exact/io/MpsWriter.hpp"
#include "exact/io/SolParser.hpp"
#include "exact/io/SolWriter.hpp"
#include "exact/misc/NumericalStatistics.hpp"
#include "exact/misc/OptionsParser.hpp"
#include "exact/misc/tbb.hpp"
#include "scip/scipgithash.h"
#include "soplex.h"
#include "soplex/spxgithash.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <string>
#include <utility>

namespace exact
{

void
join( const Vec<std::string>& v, char c, std::string& s )
{

   s.clear();

   for( auto p = v.begin(); p != v.end(); ++p )
   {
      s += *p;
      if( p != v.end() - 1 )
         s += c;
   }
}

void
print_header()
{
   std::string mode = "optimized";
#ifndef NDEBUG
   mode = "debug";
#endif

#ifdef EXACT_GITHASH_AVAILABLE
   fmt::print( "exact version {}.{}.{} [mode: {}][GitHash: {}]\n",
               EXACT_VERSION_MAJOR, EXACT_VERSION_MINOR, EXACT_VERSION_PATCH,
               mode, EXACT_GITHASH );
#else
   fmt::print( "exact version {}.{}.{} [mode: {}][GitHash: ]\n",
               EXACT_VERSION_MAJOR, EXACT_VERSION_MINOR, EXACT_VERSION_PATCH,
               mode );
#endif
   fmt::print( "Copyright (C) 2024 Zuse Institute Berlin (ZIB)\n" );
   fmt::print( "\n" );

   fmt::print( "External libraries: \n" );

#ifdef BOOST_FOUND
   fmt::print( "  Boost    {}.{}.{} \t (https://www.boost.org/)\n",
               BOOST_VERSION_NUMBER_MINOR( BOOST_VERSION ),
               BOOST_VERSION_NUMBER_PATCH( BOOST_VERSION ) / 100,
               BOOST_VERSION_NUMBER_MAJOR( BOOST_VERSION ) );
#endif

#ifdef EXACT_HAVE_PAPILO
   fmt::print( "  PaPILO     {}.{}.{} \t Parallel Presolving Integer und Linear Optimization "
               "developed at Zuse "
               "Institute Berlin (scip.zib.de) [GitHash: {}]\n",
               PAPILO_VERSION_MAJOR, PAPILO_VERSION_MINOR, PAPILO_VERSION_PATCH,
               SCIPgetGitHash() );
#endif

   fmt::print( "  TBB            \t Thread building block https://github.com/oneapi-src/oneTBB developed by Intel\n");


   fmt::print( "  SCIP     {}.{}.{} \t Mixed Integer Programming Solver "
               "developed at Zuse "
               "Institute Berlin (scip.zib.de) [GitHash: {}]\n",
               SCIP_VERSION_MAJOR, SCIP_VERSION_MINOR, SCIP_VERSION_PATCH,
               SCIPgetGitHash() );
   fmt::print(
       "  SoPlex   {}.{}.{} \t Linear Programming Solver developed at Zuse "
       "Institute Berlin (soplex.zib.de) [GitHash: {}]\n",
       SOPLEX_VERSION / 100, ( SOPLEX_VERSION % 100 ) / 10, SOPLEX_VERSION % 10,
       soplex::getGitHash() );

   fmt::print( "\n" );

}



} // namespace exact

#endif
