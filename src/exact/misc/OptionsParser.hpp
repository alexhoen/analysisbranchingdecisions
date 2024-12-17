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

#ifndef _EXACT_MISC_OPTIONS_PARSER_HPP_
#define _EXACT_MISC_OPTIONS_PARSER_HPP_

#include "exact/misc/fmt.hpp"
#include <boost/program_options.hpp>
#include <fstream>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

namespace exact
{

using namespace boost::program_options;



struct OptionsInfo
{
   std::string instance_file;
   std::string param_settings_file;
   std::vector<std::string> unparsed_options;
   double tlim = std::numeric_limits<double>::max();
   int nthreads;
   bool print_stats;
   bool print_params;
   bool is_complete;

   bool
   checkFiles()
   {
      if( existsFile(instance_file))
      {
         fmt::print( "file {} is not valid\n", instance_file );
         return false;
      }


      if( !print_params && existsFile( param_settings_file ))
      {
         fmt::print( "file {} is not valid\n", param_settings_file );
         return false;
      }

      return true;
   }

   bool
   existsFile( std::string& filename ) const
   {
      return !filename.empty() && !std::ifstream( filename );
   }

   void
   parse( const std::vector<std::string>& opts = std::vector<std::string>() )
   {
      is_complete = false;


      options_description desc( fmt::format( ""));

      desc.add_options()( "file,f", value( &instance_file ), "instance file" );

      desc.add_options( )("parameter-settings,p",
                          value(&param_settings_file),
                          "filename for presolve parameter settings");

      desc.add_options( )("threads,t",
                          value(&nthreads)->default_value(0));



      if( opts.empty() )
      {
         fmt::print( "\n{}\n", desc );
         return;
      }

      variables_map vm;

      parsed_options parsed = command_line_parser( opts )
                                  .options( desc )
                                  .allow_unregistered()
                                  .run();
      store( parsed, vm );
      notify( vm );

      if( !checkFiles() )
         return;

      unparsed_options =
          collect_unrecognized( parsed.options, exclude_positional );

      is_complete = true;
   }
};

OptionsInfo
parseOptions( int argc, char* argv[] )
{
   OptionsInfo optionsInfo;
   using namespace boost::program_options;
//   using boost::none;
   using boost::optional;
   std::string usage =
       fmt::format( "usage:\n {} [ARGUMENTS]\n", argv[0] );

   // global description.
   // will capture the command and arguments as unrecognised
   options_description global{};
   global.add_options()( "help,h", "produce help message" );
   global.add_options()( "args", value<std::vector<std::string>>(),
                         "arguments for the command" );

   positional_options_description pos;
   pos.add( "args", -1 );

   parsed_options parsed = command_line_parser( argc, argv )
                               .options( global )
                               .positional( pos )
                               .allow_unregistered()
                               .run();

   variables_map vm;
   store( parsed, vm );

   if( vm.count( "help" ) || vm.empty() )
   {
      fmt::print( "{}\n{}", usage, global );
      optionsInfo.parse( );

      return optionsInfo;
   }

   // we branch on each command
   // and parse the arguments passed with the command


   std::vector<std::string> opts =
         collect_unrecognized(parsed.options, include_positional);

//   opts.erase(opts.begin( ));

   optionsInfo.parse(opts);

   return optionsInfo;
}

} // namespace exact

#endif
