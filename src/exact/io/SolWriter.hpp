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

#ifndef _EXACT_IO_SOL_WRITER_HPP_
#define _EXACT_IO_SOL_WRITER_HPP_

#include "exact/Config.hpp"
#include "exact/misc/Vec.hpp"
#include "exact/misc/fmt.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

namespace exact
{

/// Writer to write problem structures into an mps file
template <typename REAL>
struct SolWriter
{
   static void
   writePrimalSol( const std::string& filename, const Vec<REAL>& sol,
                   const Vec<REAL>& objective, const REAL& solobj,
                   const Vec<std::string>& colnames )
   {
      std::ofstream file( filename, std::ofstream::out );
      boost::iostreams::filtering_ostream out;

#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );
#endif

      out.push( file );

      fmt::print( out, "{: <50} {: <18.15}\n", "=obj=", double( solobj ) );

      for( int i = 0; i != (int) sol.size(); ++i )
      {
         if( sol[i] != 0.0 )
         {
            fmt::print( out, "{: <50} {: <18.15}   obj({:.15})\n", colnames[i],
                        double( sol[i] ), double( objective[i] ) );
         }
      }
   }

   static void
   writeDualSol( const std::string& filename, const Vec<REAL>& sol,
                 const Vec<REAL>& rhs, const Vec<REAL>& lhs,
                 const REAL& obj_value, const Vec<std::string>& row_names )
   {
      std::ofstream file( filename, std::ofstream::out );
      boost::iostreams::filtering_ostream out;

#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );
#endif

      out.push( file );

      fmt::print( out, "{: <50} {: <18.15}\n", "=obj=", double( obj_value ) );

      for( int i = 0; i < (int) sol.size(); ++i )
      {
         if( sol[i] != 0.0 )
         {
            REAL objective = lhs[i];
            if( sol[i] < 0 )
               objective = rhs[i];
            fmt::print( out, "{: <50} {: <18.15}   obj({:.15})\n", row_names[i],
                        double( sol[i] ), double( objective ) );
         }
      }
   }

   static void
   writeReducedCostsSol( const std::string& filename, const Vec<REAL>& sol,
                         const Vec<REAL>& ub, const Vec<REAL>& lb,
                         const REAL& solobj, const Vec<std::string>& col_names )
   {
      std::ofstream file( filename, std::ofstream::out );
      boost::iostreams::filtering_ostream out;

#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );
#endif

      out.push( file );

      fmt::print( out, "{: <50} {: <18.15}\n", "=obj=", double( solobj ) );

      for( int i = 0; i < (int) sol.size(); ++i )
      {
         if( sol[i] != 0.0 )
         {
            REAL objective = lb[i];
            if( sol[i] < 0 )
               objective = ub[i];
            fmt::print( out, "{: <50} {: <18.15}   obj({:.15})\n", col_names[i],
                        double( sol[i] ), double( objective ) );
         }
      }
   }

   static void
   writeBasis( const std::string& filename, const Vec<VarBasisStatus>& colBasis,
               const Vec<VarBasisStatus>& rowBasis, const Vec<std::string>& col_names, const Vec<std::string>& row_names )
   {
      std::ofstream file( filename, std::ofstream::out );
      boost::iostreams::filtering_ostream out;

#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
         out.push( boost::iostreams::gzip_compressor() );
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
         out.push( boost::iostreams::bzip2_compressor() );
#endif

      int rowSize = (int) rowBasis.size();
      assert(colBasis.size() == col_names.size());
      assert( rowSize == row_names.size());


      out.push( file );
      int row = 0;
      fmt::print( out, "NAME  papilo.bas\n");

      for( int col = 0; col < (int) colBasis.size(); ++col )
      {
         if( colBasis[col] == VarBasisStatus::BASIC )
         {
            /* Find non basic row */
            for(; row < rowSize; row++)
            {
               if(rowBasis[row] != VarBasisStatus::BASIC)
                  break;
            }
            assert( rowBasis[row] == VarBasisStatus::ON_UPPER ||
                    rowBasis[row] == VarBasisStatus::ON_LOWER ||
                    rowBasis[row] == VarBasisStatus::ZERO ||
                    rowBasis[row] == VarBasisStatus::FIXED
                    );
            if( colBasis[col] == VarBasisStatus::ON_UPPER )
               fmt::print( out, "  XU {: <50} {: <50}\n", col_names[col],
                           row_names[row]);
            else
               fmt::print( out, "  XL {: <50} {: <50}\n", col_names[col],
                           row_names[row]);
            row++;
         }
         else if( colBasis[col] == VarBasisStatus::ON_UPPER )
            fmt::print( out, "  UL {: <50}\n", col_names[col]);
         else if( colBasis[col] == VarBasisStatus::ON_LOWER ||
                  colBasis[col] == VarBasisStatus::ZERO )
         {
            /* Default is all non-basic variables on lower bound (if finite) or
             * at zero (if free). nothing to do in this case.
             */
         }
         else
            assert( false );
      }
      fmt::print( out, "ENDDATA\n");

//      assert(check_if_remaining_rows_are_basic( rowBasis, rowSize, row ));
   }

   static bool
   check_if_remaining_rows_are_basic( const Vec<VarBasisStatus>& rowBasis,
                                      int rowSize, int row )
   {
      // Check that we covered all nonbasic rows - the remaining should be basic.
      for(; row < rowSize; row++)
      {
         if(rowBasis[row] != VarBasisStatus::BASIC)
            break;
      }

      return row == rowSize ;
   }
};

} // namespace exact

#endif
