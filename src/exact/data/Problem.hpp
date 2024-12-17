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

#ifndef _EXACT_CORE_PROBLEM_HPP_
#define _EXACT_CORE_PROBLEM_HPP_

#include "exact/data/ConstraintMatrix.hpp"
#include "exact/data/Objective.hpp"
#include "exact/data/SingleRow.hpp"
#include "exact/data/VariableDomains.hpp"
#include "exact/io/Message.hpp"
#include "exact/misc/MultiPrecision.hpp"
#include "exact/misc/StableSum.hpp"
#include "exact/misc/String.hpp"
#include "exact/misc/Vec.hpp"
#include "exact/misc/fmt.hpp"
#ifdef PAPILO_TBB
#include "exact/misc/tbb.hpp"
#endif

namespace exact
{

/// class representing the problem consisting of the constraint matrix, the left
/// and right hand side values, the variable domains, the column bounds,
/// column integrality restrictions, and the objective function
template <typename REAL>
class Problem
{
 public:
   /// set objective function
   void
   setObjective( Vec<REAL> coefficients, REAL offset = 0.0 )
   {
      objective = Objective<REAL>{ std::move( coefficients ), offset };
   }

   /// set objective function
   void
   setObjective( Objective<REAL>&& obj )
   {
      objective = obj;
   }

   /// set (transposed) constraint matrix
   void
   setConstraintMatrix( SparseStorage<REAL> cons_matrix, Vec<REAL> lhs_values,
                        Vec<REAL> rhs_values, Vec<RowFlags> row_flags,
                        bool transposed = false )
   {
      assert( lhs_values.size() == rhs_values.size() );
      assert( lhs_values.size() == row_flags.size() );
      assert( ( transposed ? cons_matrix.getNCols()
                           : cons_matrix.getNRows() ) == (int) row_flags.size() );

      auto cons_matrix_other = cons_matrix.getTranspose();
      if( transposed )
         constraintMatrix = ConstraintMatrix<REAL>{
             std::move( cons_matrix_other ), std::move( cons_matrix ),
             std::move( lhs_values ), std::move( rhs_values ),
             std::move( row_flags ) };
      else
         constraintMatrix = ConstraintMatrix<REAL>{
             std::move( cons_matrix ), std::move( cons_matrix_other ),
             std::move( lhs_values ), std::move( rhs_values ),
             std::move( row_flags ) };
   }

   /// set constraint matrix
   void
   setConstraintMatrix( ConstraintMatrix<REAL>&& cons_matrix )
   {
      constraintMatrix = cons_matrix;
   }

   /// set domains of variables
   void
   setVariableDomains( VariableDomains<REAL>&& domains )
   {
      variableDomains = domains;

      nintegers = 0;
      ncontinuous = 0;

      for( ColFlags cf : variableDomains.flags )
      {
         if( cf.test( ColFlag::kIntegral ) )
            ++nintegers;
         else
            ++ncontinuous;
      }
   }

   /// set domains of variables
   void
   setVariableDomains( Vec<REAL> lower_bounds, Vec<REAL> upper_bounds,
                       Vec<ColFlags> col_flags )
   {
      variableDomains = VariableDomains<REAL>{ std::move( lower_bounds ),
                                               std::move( upper_bounds ),
                                               std::move( col_flags ) };
      nintegers = 0;
      ncontinuous = 0;

      for( ColFlags cf : variableDomains.flags )
      {
         if( cf.test( ColFlag::kIntegral ) )
            ++nintegers;
         else
            ++ncontinuous;
      }
   }

   /// returns number of active integral columns
   int
   getNumIntegralCols() const
   {
      return nintegers;
   }

   /// returns number of active integral columns
   int&
   getNumIntegralCols()
   {
      return nintegers;
   }

   /// returns number of active continuous columns
   int
   getNumContinuousCols() const
   {
      return ncontinuous;
   }

   /// returns number of active continuous columns
   int&
   getNumContinuousCols()
   {
      return ncontinuous;
   }

   /// set variable names
   void
   setVariableNames( Vec<String> var_names )
   {
      variableNames = std::move( var_names );
   }

   /// set constraint names
   void
   setConstraintNames( Vec<String> cons_names )
   {
      constraintNames = std::move( cons_names );
   }

   /// set problem name
   void
   setName( String name_ )
   {
      this->name = std::move( name_ );
   }

   /// get the problem matrix
   const ConstraintMatrix<REAL>&
   getConstraintMatrix() const
   {
      return constraintMatrix;
   }

   /// get the problem matrix
   ConstraintMatrix<REAL>&
   getConstraintMatrix()
   {
      return constraintMatrix;
   }

   /// get number of columns
   int
   getNCols() const
   {
      return constraintMatrix.getNCols();
   }

   /// get number of rows
   int
   getNRows() const
   {
      return constraintMatrix.getNRows();
   }

   /// get the objective function
   const Objective<REAL>&
   getObjective() const
   {
      return objective;
   }

   /// get the objective function
   Objective<REAL>&
   getObjective()
   {
      return objective;
   }

   /// get the variable domains
   const VariableDomains<REAL>&
   getVariableDomains() const
   {
      return variableDomains;
   }

   /// get the variable domains
   VariableDomains<REAL>&
   getVariableDomains()
   {
      return variableDomains;
   }

   const Vec<ColFlags>&
   getColFlags() const
   {
      return variableDomains.flags;
   }

   Vec<ColFlags>&
   getColFlags()
   {
      return variableDomains.flags;
   }

   const Vec<RowFlags>&
   getRowFlags() const
   {
      return constraintMatrix.getRowFlags();
   }

   Vec<RowFlags>&
   getRowFlags()
   {
      return constraintMatrix.getRowFlags();
   }

   /// get the variable names
   const Vec<String>&
   getVariableNames() const
   {
      return variableNames;
   }

   /// get the constraint names
   const Vec<String>&
   getConstraintNames() const
   {
      return constraintNames;
   }

   /// get the problem name
   const String&
   getName() const
   {
      return name;
   }

   /// get the (dense) vector of variable lower bounds
   const Vec<REAL>&
   getLowerBounds() const
   {
      return variableDomains.lower_bounds;
   }

   /// get the (dense) vector of variable lower bounds
   Vec<REAL>&
   getLowerBounds()
   {
      return variableDomains.lower_bounds;
   }

   /// get the (dense) vector of variable upper bounds
   const Vec<REAL>&
   getUpperBounds() const
   {
      return variableDomains.upper_bounds;
   }

   /// get the (dense) vector of variable upper bounds
   Vec<REAL>&
   getUpperBounds()
   {
      return variableDomains.upper_bounds;
   }

   /// get the (dense) vector of column sizes
   const Vec<int>&
   getColSizes() const
   {
      return constraintMatrix.getColSizes();
   }

   /// get the (dense) vector of column sizes
   Vec<int>&
   getColSizes()
   {
      return constraintMatrix.getColSizes();
   }

   /// get the (dense) vector of row sizes
   const Vec<int>&
   getRowSizes() const
   {
      return constraintMatrix.getRowSizes();
   }

   /// get the (dense) vector of row sizes
   Vec<int>&
   getRowSizes()
   {
      return constraintMatrix.getRowSizes();
   }

   /// substitute a variable in the objective using an equality constraint
   /// given by a row index
   void
   substituteVarInObj( const Num<REAL>& num, int col, int equalityrow );

   bool
   computeSolViolations( const Num<REAL>& num, const Vec<REAL>& sol,
                         REAL& boundviolation, REAL& rowviolation,
                         REAL& intviolation ) const
   {
      if( (int) sol.size() != getNCols() )
         return false;

      boundviolation = 0;
      intviolation = 0;

      for( int i = 0; i != getNCols(); ++i )
      {
         if( !variableDomains.flags[i].test( ColFlag::kLbInf ) &&
             sol[i] < variableDomains.lower_bounds[i] )
         {
            REAL thisviol = variableDomains.lower_bounds[i] - sol[i];

            if( !num.isFeasZero( thisviol ) )
               Message::debug( this,
                               "lower bound {} of column {} with solution "
                               "value {} is violated by {}\n",
                               double( variableDomains.lower_bounds[i] ), i,
                               double( sol[i] ), double( thisviol ) );

            boundviolation = num.max( boundviolation, thisviol );
         }

         if( !variableDomains.flags[i].test( ColFlag::kUbInf ) &&
             sol[i] > variableDomains.upper_bounds[i] )
         {
            REAL thisviol = sol[i] - variableDomains.upper_bounds[i];

            if( !num.isFeasZero( thisviol ) )
               Message::debug( this,
                               "upper bound {} of column {} with solution "
                               "value {} is violated by {}\n",
                               double( variableDomains.upper_bounds[i] ), i,
                               double( sol[i] ), double( thisviol ) );

            boundviolation = num.max( boundviolation, thisviol );
         }

         if( variableDomains.flags[i].test( ColFlag::kIntegral ) )
         {
            REAL thisviol = abs( num.round( sol[i] ) - sol[i] );

            if( !num.isFeasZero( thisviol ) )
               Message::debug( this,
                               "integrality of column {} with solution value "
                               "{} is violated by {}\n",
                               i, double( sol[i] ), double( thisviol ) );

            intviolation = num.max( intviolation, thisviol );
         }
      }

      rowviolation = 0;

      const Vec<RowFlags>& rflags = getRowFlags();
      const Vec<REAL>& lhs = constraintMatrix.getLeftHandSides();
      const Vec<REAL>& rhs = constraintMatrix.getRightHandSides();

      for( int i = 0; i != getNRows(); ++i )
      {
         auto rowvec = constraintMatrix.getRowCoefficients( i );
         const REAL* vals = rowvec.getValues();
         const int* inds = rowvec.getIndices();

         StableSum<REAL> activitySum;
         for( int j = 0; j != rowvec.getLength(); ++j )
            activitySum.add( sol[inds[j]] * vals[j] );

         REAL activity = activitySum.get();

         if( !rflags[i].test( RowFlag::kRhsInf )
             && num.isFeasGT( activity, rhs[i] ) )
         {
            Message::debug( this,
                            "the activity {} of constraint {}  "
                            "{} is greater than the righthandside {}\n",
                            activity, i, rhs[i] );
            rowviolation = num.max( rowviolation, activity - rhs[i] );
         }

         if( !rflags[i].test( RowFlag::kLhsInf )
             && num.isFeasLT( activity, rhs[i] ) )
         {
            Message::debug( this,
                            "the activity {} of constraint {}  "
                            "{} is greater than the lefthandside {}\n",
                            activity, i, lhs[i] );
            rowviolation = num.max( rowviolation, lhs[i] - activity );
         }
      }

      return num.isFeasZero( boundviolation ) &&
             num.isFeasZero( intviolation ) && num.isFeasZero( rowviolation );
   }

   REAL
   computeSolObjective( const Vec<REAL>& sol ) const
   {
      assert( (int) sol.size() == getNCols() );

      StableSum<REAL> obj( objective.offset );
      for( int i = 0; i < getNCols(); ++i )
         obj.add( sol[i] * objective.coefficients[i] );

      return obj.get();
   }

   /// return const reference to vector of row activities
   const Vec<RowActivity<REAL>>&
   getRowActivities() const
   {
      return rowActivities;
   }

   /// return reference to vector of row activities
   Vec<RowActivity<REAL>>&
   getRowActivities()
   {
      return rowActivities;
   }

   void
   recomputeAllActivities()
   {
      rowActivities.resize( getNRows() );

      // loop through rows once, compute initial acitvities, detect trivial
      // redundancy
#ifdef PAPILO_TBB
      tbb::parallel_for(
          tbb::blocked_range<int>( 0, getNRows() ),
          [this]( const tbb::blocked_range<int>& r ) {
             for( int row = r.begin(); row < r.end(); ++row )
#else
      for( int row = 0; row < getNRows(); ++row )
#endif
             {
                auto rowvec = constraintMatrix.getRowCoefficients( row );
                rowActivities[row] = compute_row_activity(
                    rowvec.getValues(), rowvec.getIndices(), rowvec.getLength(),
                    variableDomains.lower_bounds, variableDomains.upper_bounds,
                    variableDomains.flags );
             }
#ifdef PAPILO_TBB
          } );
#endif
   }

   const bool&
   isObjectiveInteger() const
   {
      return integer_objective;
   }

   const bool&
   isDecisionProblem() const
   {
      return decision_problem;
   }


   void initflags( ) {
      bool is_integer = true;
      bool decision = true;
      for( unsigned int i = 0; i < objective.coefficients.size( ); ++i )
      {
         const Rational &obj_coeff = objective.coefficients[ i ];
         if( obj_coeff == 0 )
            continue;
         decision = false;
         if( !variableDomains.flags[ i ].test(ColFlag::kIntegral) || ceil(obj_coeff) != obj_coeff )
         {
            is_integer = false;
            break;
         }
      }
      decision_problem = decision;
      integer_objective = is_integer;
   }

   bool isScaledObjectiveInteger(double scale) const
   {
      for( unsigned int i = 0; i < objective.coefficients.size(); ++i )
      {
         const Rational &obj_coeff = objective.coefficients[ i ];
         if( obj_coeff == 0)
            continue;
         const Rational &number = obj_coeff / Rational(scale);
         if( !variableDomains.flags[ i ].test(ColFlag::kIntegral) || ( ceil(number) != number ))
            return  false;
      }
      return  true;
   }


   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& name;
      ar& inputTolerance;
      ar& objective;

      ar& constraintMatrix;
      ar& variableDomains;
      ar& ncontinuous;
      ar& nintegers;

      ar& variableNames;
      ar& constraintNames;
      ar& rowActivities;
   }

   void delete_unbounded_rows( ) {
      Vec<int> unbounded_rows{};
      for( int row = 0; row < getNRows(); ++row )
      {
         auto &rowflag = constraintMatrix.getRowFlags( )[ row ];
         if( rowflag.test(RowFlag::kRhsInf) && rowflag.test(RowFlag::kLhsInf))
         {
            rowflag.set(RowFlag::kRedundant);
            unbounded_rows.emplace_back(row);
         }
      }
      constraintMatrix.deleteRows(unbounded_rows);
      constraintMatrix.compress(true);
   }

private:
   String name;
   REAL inputTolerance{ 0 };
   Objective<REAL> objective;
   ConstraintMatrix<REAL> constraintMatrix;
   VariableDomains<REAL> variableDomains;
   int ncontinuous;
   int nintegers;

   Vec<String> variableNames;
   Vec<String> constraintNames;

   /// minimal and maximal row activities
   Vec<RowActivity<REAL>> rowActivities;

   bool integer_objective = false;
   bool decision_problem = false;
};

template <typename REAL>
void
Problem<REAL>::substituteVarInObj( const Num<REAL>& num, int col, int row )
{
   auto& consMatrix = getConstraintMatrix();
   auto& objcoefficients = getObjective().coefficients;
   REAL freevarCoefInObj = objcoefficients[col];

   if( freevarCoefInObj == REAL{ 0 } )
      return;

   const auto equalityrow = consMatrix.getRowCoefficients( row );
   const int length = equalityrow.getLength();
   const REAL* values = equalityrow.getValues();
   const int* indices = equalityrow.getIndices();

   int consid = consMatrix.getSparseIndex( col, row );
   assert( consid >= 0 );
   assert( indices[consid] == col );
   REAL freevarCoefInCons = values[consid];

   REAL substscale = -freevarCoefInObj / freevarCoefInCons;

   objcoefficients[col] = REAL{ 0.0 };
   for( int j = 0; j < length; ++j )
   {
      if( indices[j] == col )
         continue;

      REAL newobjcoeff = objcoefficients[indices[j]] + values[j] * substscale;
      if( num.isZero( newobjcoeff ) )
         newobjcoeff = 0;

      objcoefficients[indices[j]] = newobjcoeff;
   }

   assert( consMatrix.getRowFlags()[row].test( RowFlag::kEquation ) &&
           !consMatrix.getRowFlags()[row].test( RowFlag::kRhsInf ) &&
           !consMatrix.getRowFlags()[row].test( RowFlag::kLhsInf ) &&
           consMatrix.getLeftHandSides()[row] ==
               consMatrix.getRightHandSides()[row] );
   getObjective().offset -= consMatrix.getLeftHandSides()[row] * substscale;
}

} // namespace exact

#endif
