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

#ifndef _EXACT_MISC_WRAPPERS_HPP_
#define _EXACT_MISC_WRAPPERS_HPP_


#include "exact/io/MpsParser.hpp"
#include "exact/io/MpsWriter.hpp"
#include "exact/io/SolParser.hpp"
#include "exact/io/SolWriter.hpp"
#include "exact/misc/NumericalStatistics.hpp"
#include "exact/misc/OptionsParser.hpp"
#include "exact/misc/tbb.hpp"
#include "Timer.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <string>
#include <utility>

namespace exact
{

enum class ResultStatus
{
   kOk = 0,
   kUnbndOrInfeas,
   kError
};

template <typename REAL>
ResultStatus
presolve_and_solve(const OptionsInfo& opts)
{
   try
   {
      double readtime = 0;
      Problem<REAL> problem;
      boost::optional<Problem<REAL>> prob;

      {
         Timer t( readtime );
         prob = MpsParser<REAL>::loadProblem( opts.instance_file );
      }

      // Check whether reading was successful or not
      if( !prob )
      {
         fmt::print( "error loading problem {}\n", opts.instance_file );
         return ResultStatus::kError;
      }
      problem = *prob;

      fmt::print( "reading took {:.3} seconds\n", readtime );

      NumericalStatistics<REAL> nstats( problem );
      nstats.printStatistics();

      if( !opts.param_settings_file.empty() || !opts.unparsed_options.empty() ||
          opts.print_params )
      {
         ParameterSet paramSet { };

         if( !opts.param_settings_file.empty() && !opts.print_params )
         {
            std::ifstream input( opts.param_settings_file );
            if( input )
            {
               String theoptionstr;
               String thevaluestr;
               for( String line; getline( input, line ); )
               {
                  std::size_t pos = line.find_first_of( '#' );
                  if( pos != String::npos )
                     line = line.substr( 0, pos );

                  pos = line.find_first_of( '=' );

                  if( pos == String::npos )
                     continue;

                  theoptionstr = line.substr( 0, pos - 1 );
                  thevaluestr = line.substr( pos + 1 );

                  boost::algorithm::trim( theoptionstr );
                  boost::algorithm::trim( thevaluestr );

                  try
                  {
                     paramSet.parseParameter( theoptionstr.c_str(),
                                              thevaluestr.c_str() );
                     fmt::print( "set {} = {}\n", theoptionstr, thevaluestr );
                  }
                  catch( const std::exception& e )
                  {
                     fmt::print( "parameter '{}' could not be set: {}\n", line,
                                 e.what() );
                  }
               }
            }
            else
            {
               fmt::print( "could not read parameter file '{}'\n",
                           opts.param_settings_file );
            }
         }

         if( !opts.unparsed_options.empty() )
         {
            String theoptionstr;
            String thevaluestr;

            for( const auto& option : opts.unparsed_options )
            {
               std::size_t pos = option.find_first_of( '=' );
               if( pos != String::npos && pos > 2 )
               {
                  theoptionstr = option.substr( 2, pos - 2 );
                  thevaluestr = option.substr( pos + 1 );
                  try
                  {
                     paramSet.parseParameter( theoptionstr.c_str(),
                                              thevaluestr.c_str() );
                     fmt::print( "set {} = {}\n", theoptionstr, thevaluestr );
                  }
                  catch( const std::exception& e )
                  {
                     fmt::print( "parameter '{}' could not be set: {}\n",
                                 option, e.what() );
                  }
               }
               else
               {
                  fmt::print(
                      "parameter '{}' could not be set: value expected\n",
                      option );
               }
            }
         }

         if( opts.print_params )
         {
            if( !opts.param_settings_file.empty() )
            {
               std::ofstream outfile( opts.param_settings_file );

               if( outfile )
               {
                  std::ostream_iterator<char> out_it( outfile );
                  paramSet.printParams( out_it );
               }
               else
               {
                  fmt::print( "could not write to parameter file '{}'\n",
                              opts.param_settings_file );
               }
            }
            else
            {
               String paramDesc;
               paramSet.printParams( std::back_inserter( paramDesc ) );
               puts( paramDesc.c_str() );
            }
         }
      }

   }
   catch( std::bad_alloc& ex )
   {
      fmt::print( "Memory out exception occured! Please assign more memory\n" );
      return ResultStatus::kError;
   }

   return ResultStatus::kOk;
}


} // namespace exact

#endif
