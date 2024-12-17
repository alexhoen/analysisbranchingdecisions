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

#ifndef _EXACT_INTERFACES_SOPLEX_INTERFACE_HPP_
#define _EXACT_INTERFACES_SOPLEX_INTERFACE_HPP_

#include <cassert>
#include <stdexcept>
#include <type_traits>

#include "exact/data/Problem.hpp"

#include "soplex.h"

namespace exact {

   struct SoplexResult {
      NodeStatus status;
      Vec<Rational> primal;
      Vec<Rational> reduced;
      Vec<Rational> dual;
   };

   class SoplexInterface {
   private:
      soplex::SoPlex spx;

      bool rational = true;

   public:

      SoplexInterface(bool _rational) {
         rational = _rational;
      }


      void
      write_problem(const Problem<Rational> &problem) {
         soplex::NameSet cnames { };
         soplex::NameSet vnames { };
         for( const auto &i: problem.getConstraintNames( ))
            cnames.add(i.c_str( ));
         for( const auto &i: problem.getVariableNames( ))
            vnames.add(i.c_str( ));
         if( rational )
            spx.writeFileRational("test.mps", &cnames, &vnames, nullptr, true);
         else
            spx.writeFileReal("test.mps", &cnames, &vnames, nullptr, true);
      }


      void reset_to_bounds(const VariableDomains<Rational> &domains) {
         using namespace soplex;

         if( !rational )
         {
            spx.setIntParam(SoPlexBase<double>::SIMPLIFIER, SoPlexBase<double>::SIMPLIFIER_AUTO);
            spx.clearBasis( );
            for( unsigned int i = 0; i < domains.lower_bounds.size( ); ++i )
            {
               Real lb = domains.flags[ i ].test(ColFlag::kLbInf) ? -infinity : Real(domains.lower_bounds[ i ]);
               Real ub = domains.flags[ i ].test(ColFlag::kUbInf) ? infinity : Real(domains.upper_bounds[ i ]);
               spx.changeBoundsReal(i, lb, ub);
            }
         }
         else
         {
            spx.clearBasis( );
            for( unsigned int i = 0; i < domains.lower_bounds.size( ); ++i )
            {
               soplex::Rational lb = domains.flags[ i ].test(ColFlag::kLbInf)
                                     ? -spx.realParam(SoPlex::INFTY)
                                     : soplex::Rational(domains.lower_bounds[ i ]);
               soplex::Rational ub = domains.flags[ i ].test(ColFlag::kUbInf)
                                     ? spx.realParam(SoPlex::INFTY)
                                     : soplex::Rational(domains.upper_bounds[ i ]);
               spx.changeBoundsRational(i, lb, ub);
            }
         }
      }

      void reset_and_fix_integer(const VariableDomains<Rational> &domains, const Vec<Rational> &primal) {
         using namespace soplex;

         if( !rational )
         {
            spx.setIntParam(SoPlexBase<double>::SIMPLIFIER, SoPlexBase<double>::SIMPLIFIER_AUTO);
            spx.clearLPReal( );
            for( unsigned int i = 0; i < domains.lower_bounds.size( ); ++i )
            {
               Real lb, ub;
               if( !domains.flags[ i ].test(ColFlag::kIntegral))
               {
                  lb = domains.flags[ i ].test(ColFlag::kLbInf) ? -infinity : Real(domains.lower_bounds[ i ]);
                  ub = domains.flags[ i ].test(ColFlag::kUbInf) ? infinity : Real(domains.upper_bounds[ i ]);
               }
               else
               {
                  lb = Real(floor(primal[ i ] + 0.5));
                  ub = Real(lb);
               }
               spx.changeBoundsReal(i, lb, ub);
            }
         }
         else
         {
            spx.clearLPRational( );
            for( unsigned int i = 0; i < domains.lower_bounds.size( ); ++i )
            {
               soplex::Rational lb, ub;
               if( !domains.flags[ i ].test(ColFlag::kIntegral))
               {
                  lb = domains.flags[ i ].test(ColFlag::kLbInf)
                       ? -spx.realParam(SoPlex::INFTY)
                       : soplex::Rational(domains.lower_bounds[ i ]);
                  ub = domains.flags[ i ].test(ColFlag::kUbInf)
                       ? spx.realParam(SoPlex::INFTY)
                       : soplex::Rational(domains.upper_bounds[ i ]);
               }
               else
               {
                  lb = soplex::Rational(floor(primal[ i ] + 0.5));
                  ub = lb;
               }
               spx.changeBoundsRational(i, lb, ub);
            }
         }
      }


      void
      doSetUpWithFixedIntegers(const Problem<Rational> &problem, const VariableDomains<Rational> &domains,
                               const Vec<Rational> &primal) {
         using namespace soplex;

         int ncols = problem.getNCols( );
         int nrows = problem.getNRows( );
         const Objective<Rational> &obj = problem.getObjective( );
         const auto &consMatrix = problem.getConstraintMatrix( );
         const auto &lhs_values = consMatrix.getLeftHandSides( );
         const auto &rhs_values = consMatrix.getRightHandSides( );
         const auto &rflags = problem.getRowFlags( );

         /* set the objective sense and offset */
         spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);

         if( obj.offset != 0 )
            spx.setRealParam(SoPlex::OBJ_OFFSET, Real(obj.offset));

         if( !rational )
         {
            LPRowSet rows(nrows);
            LPColSet cols(ncols);
            DSVector vec(ncols);
            for( int i = 0; i < nrows; ++i )
            {
               Real lhs = rflags[ i ].test(RowFlag::kLhsInf) ? -infinity : Real(lhs_values[ i ]);
               Real rhs = rflags[ i ].test(RowFlag::kRhsInf) ? infinity : Real(rhs_values[ i ]);

               rows.add(lhs, vec, rhs);
            }

            spx.addRowsReal(rows);

            for( int i = 0; i < ncols; ++i )
            {
               assert(!domains.flags[ i ].test(ColFlag::kInactive));

               Real lb, ub;
               if( !domains.flags[ i ].test(ColFlag::kIntegral))
               {
                  lb = domains.flags[ i ].test(ColFlag::kLbInf) ? -infinity : Real(domains.lower_bounds[ i ]);
                  ub = domains.flags[ i ].test(ColFlag::kUbInf) ? infinity : Real(domains.upper_bounds[ i ]);
               }
               else
               {
                  assert(primal[ i ].convert_to<int>( ) == primal[ i ]);
                  lb = Real(primal[ i ]);
                  ub = Real(primal[ i ]);
               }
               auto colvec = consMatrix.getColumnCoefficients(i);

               int collen = colvec.getLength( );
               const int *colrows = colvec.getIndices( );
               const Rational *colvals = colvec.getValues( );

               vec.clear( );

               if( std::is_same<Rational, Real>::value )
               {
                  vec.add(collen, colrows, ( const Real * ) colvals);
               }
               else
               {
                  for( int j = 0; j != collen; ++j )
                     vec.add(colrows[ j ], Real(colvals[ j ]));
               }

               cols.add(Real(obj.coefficients[ i ]), lb, vec, ub);
            }

            spx.addColsReal(cols);
            return;
         }

         spx.setRealParam(SoPlex::OPTTOL, 0);
         spx.setRealParam(SoPlex::FEASTOL, 0);
         spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
         spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
         spx.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);

         for( int col = 0; col < ncols; ++col )
         {

            LPColRational var { };
            soplex::Rational lb, ub;

            if( !domains.flags[ col ].test(ColFlag::kIntegral))
            {
               lb = domains.flags[ col ].test(ColFlag::kLbInf)
                    ? -spx.realParam(SoPlex::INFTY)
                    : soplex::Rational(domains.lower_bounds[ col ]);
               ub = domains.flags[ col ].test(ColFlag::kUbInf)
                    ? spx.realParam(SoPlex::INFTY)
                    : soplex::Rational(domains.upper_bounds[ col ]);
            }
            else
            {
               assert(primal[ col ].convert_to<int>( ) == primal[ col ]);
               lb = Rational(primal[ col ]);
               ub = Rational(primal[ col ]);
            }
            assert(!domains.flags[ col ].test(ColFlag::kInactive) || lb == ub);
            var.setLower(lb);
            var.setUpper(ub);
            var.setObj(soplex::Rational(obj.coefficients[ col ]));
            spx.addColRational(var);

         }

         for( int row = 0; row < nrows; ++row )
         {

            assert(!rflags[ row ].test(RowFlag::kLhsInf) || !rflags[ row ].test(RowFlag::kRhsInf));
            const auto &rowvec = consMatrix.getRowCoefficients(row);
            const auto &rowinds = rowvec.getIndices( );
            const auto &rowvals = rowvec.getValues( );
            int nrowcols = rowvec.getLength( );
            soplex::Rational lhs = rflags[ row ].test(RowFlag::kLhsInf)
                                   ? -spx.realParam(SoPlex::INFTY)
                                   : soplex::Rational(lhs_values[ row ]);
            soplex::Rational rhs = rflags[ row ].test(RowFlag::kRhsInf)
                                   ? spx.realParam(SoPlex::INFTY)
                                   : soplex::Rational(rhs_values[ row ]);
            DSVectorRational cons(nrowcols);
            for( int i = 0; i < nrowcols; ++i )
            {
               assert(rowvals[ i ] != 0);
               cons.add(rowinds[ i ], soplex::Rational(rowvals[ i ]));
            }
            spx.addRowRational(LPRowRational(lhs, cons, rhs));
         }
      }


      void
      doSetUp(const Problem<Rational> &problem, const VariableDomains<Rational> &domains) {
         using namespace soplex;

         int ncols = problem.getNCols( );
         int nrows = problem.getNRows( );
         const Objective<Rational> &obj = problem.getObjective( );
         const auto &consMatrix = problem.getConstraintMatrix( );
         const auto &lhs_values = consMatrix.getLeftHandSides( );
         const auto &rhs_values = consMatrix.getRightHandSides( );
         const auto &rflags = problem.getRowFlags( );
//         const auto &domains = problem.getVariableDomains( );

         /* set the objective sense and offset */
         spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);

         if( obj.offset != 0 )
            spx.setRealParam(SoPlex::OBJ_OFFSET, Real(obj.offset));


         spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);

         if( !rational )
         {
            LPRowSet rows(nrows);
            LPColSet cols(ncols);
            DSVector vec(ncols);
            for( int i = 0; i < nrows; ++i )
            {
               Real lhs = rflags[ i ].test(RowFlag::kLhsInf) ? -infinity : Real(lhs_values[ i ]);
               Real rhs = rflags[ i ].test(RowFlag::kRhsInf) ? infinity : Real(rhs_values[ i ]);
               rows.add(lhs, vec, rhs);
            }

            spx.addRowsReal(rows);

            for( int i = 0; i < ncols; ++i )
            {
               assert(!domains.flags[ i ].test(ColFlag::kInactive));

               Real lb = domains.flags[ i ].test(ColFlag::kLbInf) ? -infinity : Real(domains.lower_bounds[ i ]);
               Real ub = domains.flags[ i ].test(ColFlag::kUbInf) ? infinity : Real(domains.upper_bounds[ i ]);

               auto colvec = consMatrix.getColumnCoefficients(i);

               int collen = colvec.getLength( );
               const int *colrows = colvec.getIndices( );
               const Rational *colvals = colvec.getValues( );

               vec.clear( );

               for( int j = 0; j != collen; ++j )
                  vec.add(colrows[ j ], Real(colvals[ j ]));

               cols.add(Real(obj.coefficients[ i ]), lb, vec, ub);
            }

            spx.addColsReal(cols);
            return;
         }

         spx.setRealParam(SoPlex::OPTTOL, 0);
         spx.setRealParam(SoPlex::FEASTOL, 0);
         spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
         spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
         spx.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);

         for( int col = 0; col < ncols; ++col )
         {

            LPColRational var { };
            soplex::Rational lb = domains.flags[ col ].test(ColFlag::kLbInf)
                                  ? -spx.realParam(SoPlex::INFTY)
                                  : soplex::Rational(domains.lower_bounds[ col ]);
            soplex::Rational ub = domains.flags[ col ].test(ColFlag::kUbInf)
                                  ? spx.realParam(SoPlex::INFTY)
                                  : soplex::Rational(domains.upper_bounds[ col ]);
            assert(!domains.flags[ col ].test(ColFlag::kInactive) || lb == ub);
            var.setLower(lb);
            var.setUpper(ub);

            var.setObj(soplex::Rational(obj.coefficients[ col ]));
            spx.addColRational(var);

         }

         for( int row = 0; row < nrows; ++row )
         {

            assert(!rflags[ row ].test(RowFlag::kLhsInf) || !rflags[ row ].test(RowFlag::kRhsInf));
            const auto &rowvec = consMatrix.getRowCoefficients(row);
            const auto &rowinds = rowvec.getIndices( );
            const auto &rowvals = rowvec.getValues( );
            int nrowcols = rowvec.getLength( );
            soplex::Rational lhs = rflags[ row ].test(RowFlag::kLhsInf)
                                   ? -spx.realParam(SoPlex::INFTY)
                                   : soplex::Rational(lhs_values[ row ]);
            soplex::Rational rhs = rflags[ row ].test(RowFlag::kRhsInf)
                                   ? spx.realParam(SoPlex::INFTY)
                                   : soplex::Rational(rhs_values[ row ]);
            DSVectorRational cons(nrowcols);
            for( int i = 0; i < nrowcols; ++i )
            {
               assert(rowvals[ i ] != 0);
               cons.add(rowinds[ i ], soplex::Rational(rowvals[ i ]));
            }
            spx.addRowRational(LPRowRational(lhs, cons, rhs));
         }
      }

      void
      setBasis(Vec<soplex::SPxSolver::VarStatus> &rowstat, Vec<soplex::SPxSolver::VarStatus> &colstat) {
         soplex::SPxSolver::VarStatus *rows = rowstat.data( );
         soplex::SPxSolver::VarStatus *cols = colstat.data( );
         spx.setBasis(rows, cols);
      }

      void
      setBasis(soplex::SPxSolver::VarStatus *rowstat, soplex::SPxSolver::VarStatus *colstat) {
         spx.setBasis(rowstat, colstat);
      }

      SoplexResult
      solve(const Vec<RowFlags> &row_flags, const Problem<Rational> &problem) {
         using namespace soplex;


         SPxSolver::Status stat = spx.optimize( );

         if( stat == SPxSolver::Status::OPTIMAL )
         {
            Vec<exact::Rational> reducedCosts;
            Vec<exact::Rational> dual;
            Vec<exact::Rational> primal;

            if( !rational )
            {

               Vec<soplex::Real> buffer;
               int numcols = spx.numColsReal( );
               buffer.resize(numcols);

               bool success = spx.getPrimalReal(buffer.data( ), numcols);
               assert(success);

               primal.resize(numcols);
               for( int i = 0; i != numcols; ++i )
                  primal[ i ] = Rational(buffer[ i ]);

               buffer.resize(numcols);

               success = spx.getRedCostReal(buffer.data( ), numcols);
               assert(success);

               reducedCosts.resize(numcols);
               for( int i = 0; i != numcols; ++i )
                  reducedCosts[ i ] = exact::Rational(buffer[ i ]);

               int numrows = spx.numRowsReal( );
               buffer.resize(numrows);
               success = spx.getDualReal(buffer.data( ), numrows);
               assert(success);

               dual.resize(numrows);
               for( int i = 0; i != numrows; ++i )
               {
                  dual[ i ] = exact::Rational(buffer[ i ]);
                  if(( dual[ i ] > 0 && row_flags[ i ].test(RowFlag::kLhsInf)) ||
                     ( dual[ i ] < 0 && row_flags[ i ].test(RowFlag::kRhsInf)))
                  {
#ifdef DEBUG
                     fmt::print("setting dual solution value of row {} to 0 ({})\n", i, dual[ i ]);
#endif
                     dual[ i ] = 0;
                  }
               }
            }
            else
            {
               int numcols = spx.numColsRational( );

               VectorBase<soplex::Rational> buffer(numcols);
               bool success = spx.getPrimalRational(buffer);
               assert(success);

               primal.resize(numcols);
               for( int i = 0; i != numcols; ++i )
                  primal[ i ] = exact::Rational(buffer[ i ]);

               VectorBase<soplex::Rational> buffer2(numcols);

               success = spx.getRedCostRational(buffer2);
               assert(success);

               reducedCosts.resize(numcols);
               for( int i = 0; i != numcols; ++i )
                  reducedCosts[ i ] = exact::Rational(buffer2[ i ]);

               int numrows = spx.numRowsRational( );
               VectorBase<soplex::Rational> buffer3(numrows);
               success = spx.getDualRational(buffer3);
               assert(success);

               dual.resize(numrows);
               for( int i = 0; i != numrows; ++i )
                  dual[ i ] = exact::Rational(buffer3[ i ]);
            }
            return { NodeStatus::Solved, primal, reducedCosts, dual };
         }
         else if( stat == SPxSolver::Status::INFEASIBLE )
         {
            Vec<exact::Rational> var_farkas;
            Vec<exact::Rational> dual_farkas;
            if( !spx.hasDualFarkas( ))
            {
               spx.setIntParam(soplex::SoPlexBase<double>::SIMPLIFIER, soplex::SoPlexBase<double>::SIMPLIFIER_OFF);
               return solve(row_flags, problem);
            }

            if( !rational )
            {
               Vec<soplex::Real> buffer;
               int numcols = spx.numColsReal( );
               buffer.resize(numcols);


               int numrows = spx.numRowsReal( );
               buffer.resize(numrows);
               bool success = spx.getDualFarkasReal(buffer.data( ), numrows);
               assert(success);

               dual_farkas.resize(numrows);
               std::fill_n(std::back_inserter(var_farkas), numcols, 0);
               for( int i = 0; i != numrows; ++i )
               {
                  dual_farkas[ i ] = exact::Rational(buffer[ i ]);
                  if(( dual_farkas[ i ] > 0 && row_flags[ i ].test(RowFlag::kLhsInf)) ||
                     ( dual_farkas[ i ] < 0 && row_flags[ i ].test(RowFlag::kRhsInf)))
                  {
#ifdef DEBUG
                     fmt::print("setting dual solution value of row {} to 0 ({})\n", i, dual_farkas[ i ]);
#endif
                     dual_farkas[ i ] = 0;
                  }
                  if( dual_farkas[ i ] == 0 )
                     continue;
                  auto row_data = problem.getConstraintMatrix( ).getRowCoefficients(i);
                  for( int j = 0; j < row_data.getLength( ); j++ )
                     var_farkas[ row_data.getIndices( )[ j ]] -= dual_farkas[ i ] * row_data.getValues( )[ j ];

               }
            }
            else
            {
               int numcols = spx.numColsRational( );
               int numrows = spx.numRowsRational( );
               VectorBase<soplex::Rational> buffer2(numrows);
               bool success = spx.getDualFarkasRational(buffer2);
               assert(success);
               dual_farkas.resize(numrows);

               std::fill_n(std::back_inserter(var_farkas), numcols, 0);
               for( int i = 0; i != numrows; ++i )
               {
                  dual_farkas[ i ] = exact::Rational(buffer2[ i ]);
                  if( dual_farkas[ i ] == 0 )
                     continue;
                  auto row_data = problem.getConstraintMatrix( ).getRowCoefficients(i);
                  for( int j = 0; j < row_data.getLength( ); j++ )
                     var_farkas[ row_data.getIndices( )[ j ]] -= dual_farkas[ i ] * row_data.getValues( )[ j ];
               }
            }

            return { NodeStatus::Infeasible, { }, var_farkas, dual_farkas };
         }
         fmt::print("\nFAILED: SoPlex returned status {}!!\n", stat);
         return { NodeStatus::LPError, { }, { }};
      }


      void
      get_basis(soplex::SPxSolver::VarStatus *rowstat, soplex::SPxSolver::VarStatus *colstat) {
         spx.getBasis(rowstat, colstat);
      }

      std::pair<bool, int>
      factorize(Vec<Rational> &primal, Vec<Rational> &dual, Vec<Rational> &redCosts,
                Vec<soplex::SPxSolver::VarStatus> rowstat, Vec<soplex::SPxSolver::VarStatus> colstat, bool abort_on_non_optimal) {
         return factorize(primal, dual, redCosts, rowstat.data( ), colstat.data( ), abort_on_non_optimal);
      }

      std::pair<bool, int>
      factorize(Vec<Rational> &primal, Vec<Rational> &dual, Vec<Rational> &redCosts,
                soplex::SPxSolver::VarStatus *rowstat, soplex::SPxSolver::VarStatus *colstat, bool abort_on_non_optimal) {
         soplex::VectorBase<Rational> vecprimal(primal.size( ), primal.data( ));
         soplex::VectorBase<Rational> vecdual(dual.size( ), dual.data( ));
         soplex::VectorBase<Rational> vecredCosts(redCosts.size( ), redCosts.data( ));
         soplex::DataArray<soplex::SPxSolver::VarStatus> daRow { };
         soplex::DataArray<soplex::SPxSolver::VarStatus> daCol { };
         daRow.insert(0, dual.size( ), rowstat);
         daCol.insert(0, primal.size( ), colstat);
         soplex::SolRational solRational(vecprimal, vecdual, vecredCosts);
         bool stoppedTime = false;
         bool stoppedIter = false;
         bool error = false;
         bool optimal = false;

         spx.factorizeColumnRational(solRational, daRow, daCol, stoppedTime, stoppedIter, error, optimal);

         if( stoppedIter || stoppedTime || error )
            return { false, 0 };

         if( abort_on_non_optimal && !optimal)
            return { false, 0 };

         solRational.getDualSol(vecdual);
         solRational.getPrimalSol(vecprimal);
         solRational.getRedCostSol(vecredCosts);

         //TODO: optimize
         int nnz = 0;
         for( int i = 0; i < vecprimal.dim( ); i++ )
            primal[ i ] = vecprimal[ i ];
         for( int i = 0; i < vecdual.dim( ); i++ )
         {
            dual[ i ] = vecdual[ i ];
            if( dual[ i ] != 0 )
               nnz++;
         }
         for( int i = 0; i < vecredCosts.dim( ); i++ )
         {
            redCosts[ i ] = vecredCosts[ i ];
            if( redCosts[ i ] != 0 )
               nnz++;
         }

         return { true, nnz };
      };

   };


} // namespace papilo

#endif
