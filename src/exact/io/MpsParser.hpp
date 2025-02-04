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

#ifndef _EXACT_IO_MPS_PARSER_HPP_
#define _EXACT_IO_MPS_PARSER_HPP_

#include "boost/random.hpp"
#include "exact/Config.hpp"
#include "exact/data/ConstraintMatrix.hpp"
#include "exact/data/Objective.hpp"
#include "exact/data/Problem.hpp"
#include "exact/data/VariableDomains.hpp"
#include "exact/misc/Flags.hpp"
#include "exact/misc/Hash.hpp"
#include "exact/misc/Num.hpp"
#include "exact/io/ParseKey.hpp"
#include "exact/io/BoundType.hpp"
#include "exact/external/pdqsort/pdqsort.h"
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/optional.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>

#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif
#include <random>

namespace exact
{

/// Parser for mps files in fixed and free format
   template <typename REAL>
   class MpsParser
   {
   public:
      static boost::optional<Problem<REAL>>
      loadProblem( const std::string& filename, int seed )
      {
         MpsParser<REAL> parser;

         Problem<REAL> problem;

         if( !parser.parseFile( filename ) )
            return boost::none;

         assert( parser.nnz >= 0 );

         Vec<REAL> obj_vec( size_t( parser.nCols ), REAL{ 0.0 } );
         Vec<int> random_col_perm{};
         Vec<int> random_row_perm{};
         random_col_perm.resize( parser.nCols );
         random_row_perm.resize( parser.nRows );
         for( int i = 0; i < parser.nCols; ++i )
            random_col_perm[i] = i;
         for( int i = 0; i < parser.nRows; ++i )
            random_row_perm[i] = i;

         std::ranlux24 randgen( seed );

         if(seed != 0 )
         {
            parser.shuffle(randgen, random_col_perm);
            parser.shuffle(randgen, random_row_perm);
         }
         for( auto i : parser.coeffobj )
            obj_vec[random_col_perm[i.first]] = i.second;

         problem.setObjective( std::move( obj_vec ), parser.objoffset );



         Vec<REAL> lb(size_t( parser.nCols ), REAL{ 0.0 } );
         Vec<REAL> ub(size_t( parser.nCols ), REAL{ 0.0 } );
         Vec<REAL> lhs(size_t( parser.nRows ), REAL{ 0.0 } );
         Vec<REAL> rhs(size_t( parser.nRows ), REAL{ 0.0 } );
         Vec<ColFlags> flags(size_t( parser.nCols ), ColFlag::kNone );
         Vec<RowFlags> rflags(size_t( parser.nRows ), RowFlag::kLhsInf );
         Vec<String> colnames(size_t( parser.nCols ), "" );

         for( int i = 0; i < parser.nCols; i++ )
            lb[ random_col_perm[ i ]] = parser.lb4cols[ i ];
         for( int i = 0; i < parser.nCols; i++ )
            ub[ random_col_perm[ i ]] = parser.ub4cols[ i ];
         for( int i = 0; i < parser.nCols; i++ )
            flags[ random_col_perm[ i ]] = parser.col_flags[ i ];
         for( int i = 0; i < parser.nCols; i++ )
            colnames[ random_col_perm[ i ]] = parser.colnames[ i ];
         for( int i = 0; i < parser.nRows; i++ )
            rhs[ random_row_perm[ i ]] = parser.rowrhs[ i ];
         for( int i = 0; i < parser.nRows; i++ )
            lhs[ random_row_perm[ i ]] = parser.rowlhs[ i ];
         for( int i = 0; i < parser.nRows; i++ )
            rflags[ random_row_perm[ i ]] = parser.row_flags[ i ];

         for( auto& [x,y,z]: parser.entries )
         {
            x = random_col_perm[ x ];
            y = random_row_perm[ y ];
         }
         pdqsort( parser.entries.begin(), parser.entries.end(),
                  []( Triplet<REAL> a, Triplet<REAL> b ) {
                     if(std::get<0>( b ) != std::get<0>( a ))
                        return std::get<0>( b ) > std::get<0>( a );
                     return std::get<1>( b ) > std::get<1>( a );
                  } );

         problem.setConstraintMatrix(
               SparseStorage<REAL>{ std::move( parser.entries ), parser.nCols,
                                    parser.nRows, true },
               std::move( lhs ), std::move( rhs ),
               std::move( rflags ), true );
         problem.setVariableDomains( std::move( lb ),
                                     std::move( ub ),
                                     std::move( flags ) );
         problem.setVariableNames( std::move( colnames ) );
         problem.setName( std::move( filename ) );
         problem.setConstraintNames( std::move( parser.rownames ) );
         problem.initflags( );
         return problem;
      }

   private:
      MpsParser() {}

      /// load LP from MPS file as transposed triplet matrix
      bool
      parseFile( const std::string& filename );

      void shuffle( std::ranlux24& random_generator, Vec<int>& array );

      bool
      parse( boost::iostreams::filtering_istream& file );

      void
      printErrorMessage( ParseKey keyword )
      {
         switch( keyword )
         {
            case ParseKey::kRows:
               std::cerr << "read error in section ROWS " << std::endl;
               break;
            case ParseKey::kCols:
               std::cerr << "read error in section COLUMNS " << std::endl;
               break;
            case ParseKey::kRhs:
               std::cerr << "read error in section RHS " << std::endl;
               break;
            case ParseKey::kBounds:
               std::cerr << "read error in section BOUNDS " << std::endl;
               break;
            case ParseKey::kRanges:
               std::cerr << "read error in section RANGES " << std::endl;
               break;
            default:
               std::cerr << "undefined read error " << std::endl;
               break;
         }
      };

      /*
       * data for mps problem
       */

      Vec<Triplet<REAL>> entries;
      Vec<std::pair<int, REAL>> coeffobj;
      Vec<REAL> rowlhs;
      Vec<REAL> rowrhs;
      Vec<std::string> rownames;
      Vec<std::string> colnames;

      HashMap<std::string, int> rowname2idx;
      HashMap<std::string, int> colname2idx;
      Vec<REAL> lb4cols;
      Vec<REAL> ub4cols;
      Vec<BoundType> row_type;
      Vec<RowFlags> row_flags;
      Vec<ColFlags> col_flags;
      REAL objoffset = 0;

      int nCols = 0;
      int nRows = 0;
      int nnz = -1;

      /// checks first word of strline and wraps it by it_begin and it_end
      ParseKey
      checkFirstWord( std::string& strline, std::string::iterator& it,
                      boost::string_ref& word_ref ) const;

      ParseKey
      parseDefault( boost::iostreams::filtering_istream& file ) const;

      ParseKey
      parseRows( boost::iostreams::filtering_istream& file,
                 Vec<BoundType>& rowtype );

      ParseKey
      parseCols( boost::iostreams::filtering_istream& file,
                 const Vec<BoundType>& rowtype );

      ParseKey
      parseRhs( boost::iostreams::filtering_istream& file );

      ParseKey
      parseRanges( boost::iostreams::filtering_istream& file );

      ParseKey
      parseBounds( boost::iostreams::filtering_istream& file );

      std::pair<bool, REAL>
      read_number( const std::string& s );

      REAL
      pow( int base, int exponent );
   };

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::checkFirstWord( std::string& strline,
                                    std::string::iterator& it,
                                    boost::string_ref& word_ref ) const
   {
      using namespace boost::spirit;

      it = strline.begin() + strline.find_first_not_of( " " );
      std::string::iterator it_start = it;

      // TODO: Daniel
      qi::parse( it, strline.end(), qi::lexeme[+qi::graph] );

      const std::size_t length = std::distance( it_start, it );

      boost::string_ref word( &( *it_start ), length );

      word_ref = word;

      if( word.front() == 'R' ) // todo
      {
         if( word == "ROWS" )
            return ParseKey::kRows;
         else if( word == "RHS" )
            return ParseKey::kRhs;
         else if( word == "RANGES" )
            return ParseKey::kRanges;
         else
            return ParseKey::kNone;
      }
      else if( word == "COLUMNS" )
         return ParseKey::kCols;
      else if( word == "BOUNDS" )
         return ParseKey::kBounds;
      else if( word == "ENDATA" )
         return ParseKey::kEnd;
      else
         return ParseKey::kNone;
   }

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::parseDefault( boost::iostreams::filtering_istream& file ) const
   {
      std::string strline;
      getline( file, strline );

      std::string::iterator it;
      boost::string_ref word_ref;
      return checkFirstWord( strline, it, word_ref );
   }

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::parseRows( boost::iostreams::filtering_istream& file,
                               Vec<BoundType>& rowtype )
   {
      using namespace boost::spirit;

      std::string strline;
      size_t nrows = 0;
      bool hasobj = false;

      while( getline( file, strline ) )
      {
         bool isobj = false;
         std::string::iterator it;
         boost::string_ref word_ref;
         ParseKey key = checkFirstWord( strline, it, word_ref );

         // start of new section?
         if( key != ParseKey::kNone )
         {
            nRows = int( nrows );
            if( !hasobj )
            {
               std::cout << "WARNING: no objective row found" << std::endl;
               rowname2idx.emplace( "artificial_empty_objective", -1 );
            }

            return key;
         }

         if( word_ref.front() == 'G' )
         {
            rowlhs.push_back( REAL{ 0.0 } );
            rowrhs.push_back( REAL{ 0.0 } );
            row_flags.emplace_back( RowFlag::kRhsInf );
            rowtype.push_back( BoundType::kGE );
         }
         else if( word_ref.front() == 'E' )
         {
            rowlhs.push_back( REAL{ 0.0 } );
            rowrhs.push_back( REAL{ 0.0 } );
            row_flags.emplace_back( RowFlag::kEquation );
            rowtype.push_back( BoundType::kEq );
         }
         else if( word_ref.front() == 'L' )
         {
            rowlhs.push_back( REAL{ 0.0 } );
            rowrhs.push_back( REAL{ 0.0 } );
            row_flags.emplace_back( RowFlag::kLhsInf );
            rowtype.push_back( BoundType::kLE );
         }
            // todo properly treat multiple free rows
         else if( word_ref.front() == 'N' )
         {
            if( hasobj )
            {
               rowlhs.push_back( REAL{ 0.0 } );
               rowrhs.push_back( REAL{ 0.0 } );
               RowFlags rowf;
               rowf.set( RowFlag::kLhsInf, RowFlag::kRhsInf );
               row_flags.emplace_back( rowf );
               rowtype.push_back( BoundType::kLE );
            }
            else
            {
               isobj = true;
               hasobj = true;
            }
         }
         else if( word_ref.empty() ) // empty line
            continue;
         else
            return ParseKey::kFail;

         std::string rowname = ""; // todo use ref

         // get row name
         qi::phrase_parse( it, strline.end(), qi::lexeme[+qi::graph], ascii::space,
                           rowname ); // todo use ref

         // todo whitespace in name possible?
         auto ret = rowname2idx.emplace( rowname, isobj ? ( -1 ) : ( nrows++ ) );

         if( !isobj )
            rownames.push_back( rowname );

         if( !ret.second )
         {
            std::cerr << "duplicate row " << rowname << std::endl;
            return ParseKey::kFail;
         }
      }

      return ParseKey::kFail;
   }

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::parseCols( boost::iostreams::filtering_istream& file,
                               const Vec<BoundType>& rowtype )
   {
      using namespace boost::spirit;

      std::string colname = "";
      std::string strline;
      int rowidx;
      int ncols = 0;
      int colstart = 0;
      bool integral_cols = false;

      auto parsename = [&rowidx, this]( std::string name ) {
         auto mit = rowname2idx.find( name );

         assert( mit != rowname2idx.end() );
         rowidx = mit->second;

         if( rowidx >= 0 )
            this->nnz++;
         else
            assert( -1 == rowidx );
      };

      auto addtuple = [&rowidx, &ncols, this]( std::string sval) {
         auto result = read_number(sval);
         if(result.first)
         {
            fmt::print("could not parse {}\n", sval);
            exit(0);
         }
         REAL coeff = result.second;
         if( rowidx >= 0 )
            entries.push_back(
                  std::make_tuple( ncols - 1, rowidx, coeff ) );
         else
            coeffobj.push_back( std::make_pair( ncols - 1,  coeff ) );
      };

      while( getline( file, strline ) )
      {
         std::string::iterator it;
         boost::string_ref word_ref;
         ParseKey key = checkFirstWord( strline, it, word_ref );

         // start of new section?
         if( key != ParseKey::kNone )
         {
            if( ncols > 1 )
               pdqsort( entries.begin() + colstart, entries.end(),
                        []( Triplet<REAL> a, Triplet<REAL> b ) {
                           return std::get<1>( b ) > std::get<1>( a );
                        } );

            return key;
         }

         // check for integrality marker
         std::string marker = ""; // todo use ref
         std::string::iterator it2 = it;

         qi::phrase_parse( it2, strline.end(), qi::lexeme[+qi::graph],
                           ascii::space, marker );

         if( marker == "'MARKER'" )
         {
            marker = "";
            qi::phrase_parse( it2, strline.end(), qi::lexeme[+qi::graph],
                              ascii::space, marker );

            if( ( integral_cols && marker != "'INTEND'" ) ||
                ( !integral_cols && marker != "'INTORG'" ) )
            {
               std::cerr << "integrality marker error " << std::endl;
               return ParseKey::kFail;
            }
            integral_cols = !integral_cols;

            continue;
         }

         // new column?
         if( !( word_ref == colname ) )
         {
            if( word_ref.empty() ) // empty line
               continue;

            colname = word_ref.to_string();
            auto ret = colname2idx.emplace( colname, ncols++ );
            colnames.push_back( colname );

            if( !ret.second )
            {
               std::cerr << "duplicate column " << std::endl;
               return ParseKey::kFail;
            }

            assert( lb4cols.size() == col_flags.size() );

            col_flags.emplace_back( integral_cols ? ColFlag::kIntegral
                                                  : ColFlag::kNone );

            // initialize with default bounds
            if( integral_cols )
            {
               lb4cols.push_back( REAL{ 0.0 } );
               ub4cols.push_back( REAL{ 1.0 } );
            }
            else
            {
               lb4cols.push_back( REAL{ 0.0 } );
               ub4cols.push_back( REAL{ 0.0 } );
               col_flags.back().set( ColFlag::kUbInf );
            }

            assert( col_flags.size() == lb4cols.size() );

            if( ncols > 1 )
               pdqsort( entries.begin() + colstart, entries.end(),
                        []( Triplet<REAL> a, Triplet<REAL> b ) {
                           return std::get<1>( b ) > std::get<1>( a );
                        } );

            colstart = entries.size();
         }

         assert( ncols > 0 );

         std::istringstream is( strline );
         std::vector<std::string> tokens;
         std::string tmp;
         while( is >> tmp )
            tokens.push_back( tmp );
         if( tokens.size() != 3 && tokens.size() != 5 )
            return ParseKey::kFail;
         parsename( tokens[1] );
         addtuple( tokens[2] );
         if( tokens.size() == 5 )
         {
            parsename( tokens[3] );
            addtuple( tokens[4] );
         }
      }

      return ParseKey::kFail;
   }

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::parseRanges( boost::iostreams::filtering_istream& file )
   {
      using namespace boost::spirit;
      std::string strline;
      assert( rowrhs.size() == rowlhs.size() );

      while( getline( file, strline ) )
      {
         std::string::iterator it;
         boost::string_ref word_ref;
         ParseKey key = checkFirstWord( strline, it, word_ref );

         // start of new section?
         if( key != ParseKey::kNone && key != ParseKey::kRanges )
            return key;

         if( word_ref.empty() )
            continue;

         int rowidx;

         auto parsename = [&rowidx, this]( std::string name ) {
            auto mit = rowname2idx.find( name );

            assert( mit != rowname2idx.end() );
            rowidx = mit->second;

            assert( rowidx >= 0 && rowidx < nRows );
         };

         auto addrange = [&rowidx, this]( std::string sval ) {
            auto result = read_number(sval);
            if(result.first)
            {
               fmt::print("could not parse {}\n", sval);
               exit(0);
            }
            REAL val = result.second;
            assert( size_t( rowidx ) < rowrhs.size() );

            if( row_type[rowidx] == BoundType::kGE )
            {
               row_flags[rowidx].unset( RowFlag::kRhsInf );
               rowrhs[rowidx] = rowlhs[rowidx] + REAL(abs( val ));
            }
            else if( row_type[rowidx] == BoundType::kLE )
            {
               row_flags[rowidx].unset( RowFlag::kLhsInf );
               rowlhs[rowidx] = rowrhs[rowidx] - REAL(abs( val ));
            }
            else
            {
               assert( row_type[rowidx] == BoundType::kEq );
               assert( rowrhs[rowidx] == rowlhs[rowidx] );
               assert( row_flags[rowidx].test(RowFlag::kEquation) );

               if( val > REAL{ 0.0 } )
               {
                  row_flags[rowidx].unset( RowFlag::kEquation );
                  rowrhs[rowidx] = rowrhs[rowidx] + REAL( val );
               }
               else if( val < REAL{ 0.0 } )
               {
                  rowlhs[rowidx] = rowlhs[rowidx] + REAL( val );
                  row_flags[rowidx].unset( RowFlag::kEquation );
               }
            }
         };

         std::istringstream is( strline );
         std::vector<std::string> tokens;
         std::string tmp;
         while( is >> tmp )
            tokens.push_back( tmp );
         if( tokens.size() != 3 && tokens.size() != 5 )
            return ParseKey::kFail;
         parsename( tokens[1] );
         addrange( tokens[2] );
         if( tokens.size() == 5 )
         {
            parsename( tokens[3] );
            addrange( tokens[4] );
         }
      }

      return ParseKey::kFail;
   }

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::parseRhs( boost::iostreams::filtering_istream& file )
   {
      using namespace boost::spirit;
      std::string strline;

      while( getline( file, strline ) )
      {
         std::string::iterator it;
         boost::string_ref word_ref;
         ParseKey key = checkFirstWord( strline, it, word_ref );

         // start of new section?
         if( key != ParseKey::kNone && key != ParseKey::kRhs )
            return key;

         if( word_ref.empty() )
            continue;

         int rowidx;

         auto parsename = [&rowidx, this]( std::string name ) {
            auto mit = rowname2idx.find( name );

            assert( mit != rowname2idx.end() );
            rowidx = mit->second;

            assert( rowidx >= -1 );
            assert( rowidx < nRows );
         };

         auto addrhs = [&rowidx, this]( std::string sval ) {
            auto result = read_number(sval);
            if(result.first)
            {
               fmt::print("could not parse {}\n", sval);
               exit(0);
            }
            REAL val = result.second;
            if( rowidx == -1 )
            {
               objoffset = -REAL{ val };
               return;
            }
            if( row_type[rowidx] == BoundType::kEq ||
                row_type[rowidx] == BoundType::kLE )
            {
               assert( size_t( rowidx ) < rowrhs.size() );
               rowrhs[rowidx] = REAL{ val };
               row_flags[rowidx].unset( RowFlag::kRhsInf );
            }

            if( row_type[rowidx] == BoundType::kEq ||
                row_type[rowidx] == BoundType::kGE )
            {
               assert( size_t( rowidx ) < rowlhs.size() );
               rowlhs[rowidx] = REAL{ val };
               row_flags[rowidx].unset( RowFlag::kLhsInf );
            }
         };

         std::istringstream is( strline );
         std::vector<std::string> tokens;
         std::string tmp;
         while( is >> tmp )
            tokens.push_back( tmp );
         if( tokens.size() != 3 && tokens.size() != 5 )
            return ParseKey::kFail;
         parsename( tokens[1] );
         addrhs( tokens[2] );
         if( tokens.size() == 5 )
         {
            parsename( tokens[3] );
            addrhs( tokens[4] );
         }
      }

      return ParseKey::kFail;
   }

   template <typename REAL>
   ParseKey
   MpsParser<REAL>::parseBounds( boost::iostreams::filtering_istream& file )
   {
      using namespace boost::spirit;
      std::string strline;

      Vec<bool> ub_is_default( lb4cols.size(), true );
      Vec<bool> lb_is_default( lb4cols.size(), true );

      while( getline( file, strline ) )
      {
         std::string::iterator it;
         boost::string_ref word_ref;
         ParseKey key = checkFirstWord( strline, it, word_ref );

         // start of new section?
         if( key != ParseKey::kNone )
            return key;

         if( word_ref.empty() )
            continue;

         bool islb = false;
         bool isub = false;
         bool isintegral = false;
         bool isdefaultbound = false;

         if( word_ref == "UP" ) // lower bound
            isub = true;
         else if( word_ref == "LO" ) // upper bound
            islb = true;
         else if( word_ref == "FX" ) // fixed
         {
            islb = true;
            isub = true;
         }
         else if( word_ref == "MI" ) // infinite lower bound
         {
            islb = true;
            isdefaultbound = true;
         }
         else if( word_ref == "PL" ) // infinite upper bound (redundant)
         {
            isub = true;
            isdefaultbound = true;
         }
         else if( word_ref == "BV" ) // binary
         {
            isintegral = true;
            isdefaultbound = true;
            islb = true;
            isub = true;
         }
         else if( word_ref == "LI" ) // integer lower bound
         {
            islb = true;
            isintegral = true;
         }
         else if( word_ref == "UI" ) // integer upper bound
         {
            isub = true;
            isintegral = true;
         }
         else if( word_ref == "FR" ) // free variable
         {
            islb = true;
            isub = true;
            isdefaultbound = true;
         }
         else
         {
            if( word_ref == "INDICATORS" )
               std::cerr << "PaPILO does not support INDICATORS in the MPS file!!"<< std::endl;
            else
               std::cerr << "unknown bound type " << word_ref << std::endl;
            return ParseKey::kFail;
         }

         // parse over next word
         qi::phrase_parse( it, strline.end(), qi::lexeme[+qi::graph], ascii::space );

         int colidx;

         auto parsename = [&colidx, this]( std::string name ) {
            auto mit = colname2idx.find( name );
            assert( mit != colname2idx.end() );
            colidx = mit->second;
            assert( colidx >= 0 );
         };

         if( isdefaultbound )
         {
            if( !qi::phrase_parse(
                  it, strline.end(),
                  ( qi::lexeme[qi::as_string[+qi::graph][( parsename )]] ),
                  ascii::space ) )
               return ParseKey::kFail;

            if( isintegral ) // binary
            {
               if( islb )
                  lb4cols[colidx] = REAL{ 0.0 };
               if( isub )
               {
                  col_flags[colidx].unset( ColFlag::kUbInf );
                  ub4cols[colidx] = REAL{ 1.0 };
               }
               col_flags[colidx].set( ColFlag::kIntegral );
            }
            else
            {
               if( islb )
                  col_flags[colidx].set( ColFlag::kLbInf );
               if( isub )
                  col_flags[colidx].set( ColFlag::kUbInf );
            }
            continue;
         }

         auto adddomains = [&ub_is_default, &lb_is_default, &colidx, &islb, &isub, &isintegral, this]
               ( std::string sval )
         {
            auto result = read_number(sval);
            if(result.first)
            {
               fmt::print("could not parse {}\n", sval);
               exit(0);
            }
            REAL val = result.second;
            if( islb )
            {
               lb4cols[colidx] = REAL{ val };
               lb_is_default[colidx] = false;
               col_flags[colidx].unset( ColFlag::kLbInf );
            }
            if( isub )
            {
               ub4cols[colidx] = REAL{ val };
               ub_is_default[colidx] = false;
               col_flags[colidx].unset( ColFlag::kUbInf );
            }

            if( isintegral )
               col_flags[colidx].set( ColFlag::kIntegral );

            if( col_flags[colidx].test( ColFlag::kIntegral ) )
            {
               col_flags[colidx].set( ColFlag::kIntegral );
               if( !islb && lb_is_default[colidx] )
                  lb4cols[colidx] = REAL{ 0.0 };
               if( !isub && ub_is_default[colidx] )
                  col_flags[colidx].set( ColFlag::kUbInf );
            }
         };

         std::istringstream is( strline );
         std::vector<std::string> tokens;
         std::string tmp;
         while( is >> tmp )
            tokens.push_back( tmp );
         if( tokens.size() != 4 )
            return ParseKey::kFail;
         parsename( tokens[2] );
         adddomains( tokens[3] );

      }

      return ParseKey::kFail;
   }

   template <typename REAL>
   bool
   MpsParser<REAL>::parseFile( const std::string& filename )
   {
      std::ifstream file( filename, std::ifstream::in );
      boost::iostreams::filtering_istream in;

      if( !file )
         return false;

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_ZLIB
      if( boost::algorithm::ends_with( filename, ".gz" ) )
      in.push( boost::iostreams::gzip_decompressor() );
#endif

#ifdef PAPILO_USE_BOOST_IOSTREAMS_WITH_BZIP2
      if( boost::algorithm::ends_with( filename, ".bz2" ) )
      in.push( boost::iostreams::bzip2_decompressor() );
#endif

      in.push( file );

      return parse( in );
   }

   template <typename REAL>
   REAL
   MpsParser<REAL>::pow(int base, int exponent)
   {
      REAL answer = 1;
      for(int i = 0; i < exponent; i++)
         answer *= 10;
      return answer;
   }

   template <typename REAL>
   std::pair<bool, REAL>
   MpsParser<REAL>::read_number( const String& s )
   {
      using Integral = boost::multiprecision::mpz_int;
      REAL number;
      try
      {
         std::stringstream string_stream;
         string_stream.str( s );
         string_stream >> number;
         if( !string_stream.fail() && string_stream.eof() )
            return { false, number };
      }
      catch( ... ) { }

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
                  if( num_traits<REAL>::is_rational )
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
                  if( num_traits<REAL>::is_rational )
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
               else if( ( c == 'e' || c == 'E' ) && num_traits<REAL>::is_rational )
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
      number = REAL( Rational( numerator, denominator ) );

      return { failure, number };
   }


   template <typename REAL>
   bool
   MpsParser<REAL>::parse( boost::iostreams::filtering_istream& file )
   {
      nnz = 0;
      ParseKey keyword = ParseKey::kNone;
      ParseKey keyword_old = ParseKey::kNone;

      // parsing loop
      while( keyword != ParseKey::kFail && keyword != ParseKey::kEnd &&
             !file.eof() && file.good() )
      {
         keyword_old = keyword;
         switch( keyword )
         {
            case ParseKey::kRows:
               keyword = parseRows( file, row_type );
               break;
            case ParseKey::kCols:
               keyword = parseCols( file, row_type );
               break;
            case ParseKey::kRhs:
               keyword = parseRhs( file );
               break;
            case ParseKey::kRanges:
               keyword = parseRanges( file );
               break;
            case ParseKey::kBounds:
               keyword = parseBounds( file );
               break;
            case ParseKey::kFail:
               break;
            default:
               keyword = parseDefault( file );
               break;
         }
      }

      if( keyword == ParseKey::kFail || keyword != ParseKey::kEnd )
      {
         printErrorMessage( keyword_old );
         return false;
      }

      assert( row_type.size() == unsigned( nRows ) );

//      std::ranlux24 randgen( 1 );
//
//      shuffle(randgen, colname2idx);
      nCols = colname2idx.size();
      nRows = rowname2idx.size() - 1; // subtract obj row

      return true;
   }

   template <typename REAL>
   void
   MpsParser<REAL>::shuffle( std::ranlux24& random_generator, Vec<int>& array )
   {
      int tmp;
      int i;
      int end = (int)array.size();

      int begin = 0;
      // loop backwards through all elements and always swap the current last
      // element to a random position
      while( end > begin + 1 )
      {
         end--;

         // get a random position into which the last entry should be shuffled
         boost::random::uniform_int_distribution<> distrib( begin, end );
         i = distrib( random_generator );

         // swap the last element and the random element
         tmp = array[i];
         array[i] = array[end];
         array[end] = tmp;
      }
   }

} // namespace exact

#endif /* _PARSING_MPS_PARSER_HPP_ */
