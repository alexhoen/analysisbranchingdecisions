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

#ifndef _EXACT_IO_SOL_PARSER_HPP_
#define _EXACT_IO_SOL_PARSER_HPP_

#include "exact/misc/Hash.hpp"
#include "exact/misc/String.hpp"
#include "exact/misc/Vec.hpp"
#include "exact/data/Solution.hpp"
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>

namespace exact
{

template <typename REAL>
struct SolParser
{

   static bool
   read( const std::string& filename, const Vec<String>& colnames, Vec<Rational>& solution_vector )
   {
      std::ifstream file( filename, std::ifstream::in );
      boost::iostreams::filtering_istream in;

      if( !file )
         return false;

      in.push( file );

      HashMap<String, int> nameToCol;

      for( size_t i = 0; i != colnames.size(); ++i )
      {
         nameToCol.emplace( colnames[i], i );
      }

      solution_vector.resize( colnames.size(), Rational{ 0 } );
      String strline;

      std::pair<bool, Rational> result;
      do
      {
         auto tokens = split( strline.c_str() );
         assert( !tokens.empty() );

         auto it = nameToCol.find( tokens[0] );
         if( it != nameToCol.end() )
         {
            assert( tokens.size() > 1 );
            tokens[1].erase( std::remove(tokens[1].begin(), tokens[1].end(), '\r'), tokens[1].end() );
            result = parse_number( tokens[1] );
            if( result.first )
            {
               fmt::print("Could not parse solution {}\n", tokens[1]);
               return false;
            }
            solution_vector[it->second] = result.second;
         }
         else if(strline.empty()){}
         else
         {
            fmt::print( stderr,
                        "WARNING: skipping unknown column {} in solution\n",
                        tokens[0] );
         }
      } while( getline( in, strline ) );

      return true;
   }

   Vec<String> static split( const char* str )
   {
      Vec<String> tokens;
      char c1 = ' ';
      char c2 = '\t';

      do
      {
         const char* begin = str;

         while( *str != c1 && *str != c2 && *str )
            str++;

         tokens.emplace_back( begin, str );

         while( ( *str == c1 || *str == c2 ) && *str )
            str++;

      } while( 0 != *str );

      return tokens;
   }

   std::pair<bool, Rational> static
   parse_number( const String& s )
   {
      Rational number;
      try
      {
         std::stringstream string_stream;
         string_stream.str( s );
         string_stream >> number;
         if( !string_stream.fail() && string_stream.eof() )
            return { false, number };
      }
      catch( ... ) { }
      using Integral = boost::multiprecision::mpz_int;

      bool failure = false;
      Integral numerator = 0;
      Integral denominator = 1;
      unsigned exponent = 0;
      unsigned phase = 0;
      bool num_negated = false;
      bool exp_negated = false;
      for( char c : s )
      {
         int digit = '0' <= c && c <= '9' ? c - '0' : -1;
         switch( phase )
         {
            // number sign
            case 0:
               ++phase;
               if( c == '+' )
                  break;
               else if( c == '-' )
               {
                  num_negated = true;
                  break;
               }
               // before delimiter
            case 1:
               if( digit >= 0 )
               {
                  numerator *= 10;
                  numerator += digit;
                  break;
               }
               else
               {
                  ++phase;
                  if( num_traits<Rational>::is_rational )
                  {
                     if( c == '.' )
                        break;
                  }
                  else
                  {
                     if( c == '/' )
                     {
                        denominator = 0;
                        break;
                     }
                  }
               }
               // after delimiter
            case 2:
               if( digit >= 0 )
               {
                  if( num_traits<Rational>::is_rational )
                  {
                     numerator *= 10;
                     numerator += digit;
                     denominator *= 10;
                  }
                  else
                  {
                     denominator *= 10;
                     denominator += digit;
                  }
               }
               else if( ( c == 'e' || c == 'E' ) && num_traits<Rational>::is_rational )
                  ++phase;
               else
                  failure = true;
               break;
               // exponent sign
            case 3:
               ++phase;
               if( c == '+' )
                  break;
               else if( c == '-' )
               {
                  exp_negated = true;
                  break;
               }
               // exponent value
            case 4:
               if( digit >= 0 )
               {
                  exponent *= 10;
                  exponent += digit;
                  break;
               }
               else
                  ++phase;
            default:
               failure = true;
               break;
         }
         if( failure )
            break;
      }
      if( denominator == 0 )
      {
         denominator = 1;
         failure = true;
      }
      if( num_negated )
         numerator *= -1;
      if( exp_negated )
         denominator *= boost::multiprecision::pow( Integral( 10 ), exponent );
      else
         numerator *= boost::multiprecision::pow( Integral( 10 ), exponent );
      number = Rational( Rational( numerator, denominator ) );

      return { failure, number };
   }
};

} // namespace exact

#endif
