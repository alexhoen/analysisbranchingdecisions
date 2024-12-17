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


#include "exact/misc/MultiPrecision.hpp"
#include "exact/misc/OptionsParser.hpp"
#include "exact/misc/VersionLogger.hpp"
#include "exact/misc/Timer.hpp"
#include "exact/interfaces/ScipInterface.hpp"
#include "exact/interfaces/ViprInterface.hpp"
#include "exact/interfaces/EventCatcher.hpp"
#include "exact/core/PropagationView.hpp"
#include "exact/ExactSolver.hpp"

#include <boost/program_options.hpp>
#include <fstream>

void start(double readtime, exact::Problem<exact::Rational> &problem);

int
main( int argc, char* argv[] )
{
   using namespace exact;

   print_header();



   // get the options passed by the user
   OptionsInfo optionsInfo;
   try
   {
      optionsInfo = parseOptions( argc, argv );
   }
   catch( const boost::program_options::error& ex )
   {
      std::cerr << "Error while parsing the options.\n" << '\n';
      std::cerr << ex.what() << '\n';
      return 1;
   }

   if( !optionsInfo.is_complete )
      return 0;


   ExactOptions options{};
   ParameterSet paramSet{};
   options.addParameters(paramSet);
   if( !optionsInfo.param_settings_file.empty() )
   {
      std::ifstream input( optionsInfo.param_settings_file );
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
                     optionsInfo.param_settings_file );
      }
   }

   double readtime = 0;
   Problem<exact::Rational> problem;
   boost::optional<Problem<exact::Rational>> prob;

   {
      Timer t( readtime );
      prob = MpsParser<exact::Rational>::loadProblem( optionsInfo.instance_file, options.randomseed );
   }

   // Check whether reading was successful or not
   if( !prob )
   {
      fmt::print( "error loading problem {}\n", optionsInfo.instance_file );
      return 0;
   }
   problem = *prob;
   fmt::print("reading took {:.3} seconds\n", readtime );



   ExactSolver exact{problem};
   exact.start(options);


   return 0;
}

