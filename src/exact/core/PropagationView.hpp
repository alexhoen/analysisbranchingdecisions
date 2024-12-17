/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*               This file is part of the program and library                */
/*    EXACT                                                                  */
/*                                                                           */
/* Copyright (C) 2020-2024 Zuse Institute Berlin (ZIB)                       */
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

#ifndef _EXACT_CORE_PROBING_VIEW_HPP_
#define _EXACT_CORE_PROBING_VIEW_HPP_

#include "exact/data/Problem.hpp"
#include "exact/core/SingleRow.hpp"
#include <memory>
#include "exact/misc/Vec.hpp"

namespace exact
{

struct ProbingBoundChg
{
   Rational bound;
   unsigned int col : 31;
   unsigned int upper : 1;
   int probing_col : 32;

   ProbingBoundChg( bool upper_, int col_, Rational bound_, int probing_col_ )
   {
      this->upper = upper_ ? 1 : 0;
      this->col = static_cast<unsigned int>( col_ );
      this->bound = bound_;
      this->probing_col = ( probing_col_ );
   }
};

struct ProbingBoundChgReason {
   Rational bound;
   int var_index;
   bool upper;
   int row_reason;

   ProbingBoundChgReason(bool upper_, int col_, Rational bound_, int row_reason_) {
      this->upper = upper_;
      this->var_index = col_;
      this->bound = bound_;
      this->row_reason = ( row_reason_ );
   }
};


struct ProbingSubstitution
{
   Rational col2scale;
   Rational col2const;
   int col1;
   int col2;

   ProbingSubstitution( int col1_, Rational col2scale_, int col2_, Rational col2const_ )
       : col2scale( col2scale_ ), col2const( col2const_ ), col1( col1_ ),
         col2( col2_ )
   {
   }
};

class ProbingView
{

private:
   // reference to problem and numerics class
   const Problem<Rational>& problem;
   Num<Rational> num;
   Rational minintdomred;
   Rational mincontdomred;

   Vec<ProbingBoundChgReason> changed_bounds;
   Vec<int> changed_activities;
   Vec<Rational> propagated_lower_bounds;
   Vec<Rational> propagated_upper_bounds;
   Vec<ColFlags> propagated_domain_flags;
   Vec<RowActivity<Rational>> probing_activities;

   Vec<int> prop_activities;
   Vec<int> next_prop_activities;

   bool infeasible;
   int round;
   int probingCol;
   Rational probingValue;

   // results of probing and statistics
   Vec<ProbingBoundChg> boundChanges;

 public:
   ProbingView( const Problem<Rational>& problem_ )
   //TODO: get rid of this
         : problem( problem_ ), num({ } ),
           propagated_lower_bounds(problem_.getLowerBounds() ),
           propagated_upper_bounds(problem_.getUpperBounds() ),
           propagated_domain_flags(problem_.getColFlags() ),
           probing_activities( problem_.getRowActivities() )
   {
      round = -2;
      infeasible = false;
      probingCol = -1;
      probingValue = false;
      minintdomred = num.getFeasTol() * 1000;
      mincontdomred = 0.3;
   }
   void
   setMinIntDomRed( const Rational& value )
   {
      this->minintdomred = value;
   }

   void
   setMinContDomRed( const Rational& value )
   {
      this->mincontdomred = value;
   }


   void
   activityChanged( ActivityChange actchange, int rowid,
                    RowActivity<Rational>& activity );
   void
   changeLb( int col, Rational& newlb, int row );

   void
   changeUb( int col, Rational& newub, int row );

   void
   propagateDomains();

   bool
   isInfeasible() const
   {
      return infeasible;
   }


   Vec<ProbingBoundChgReason>&
   getProbingBoundChanges()
   {
      return changed_bounds;
   }

   const Vec<Rational>&
   getProbingLowerBounds() const
   {
      return propagated_lower_bounds;
   }

   const Vec<Rational>&
   getProbingUpperBounds() const
   {
      return propagated_upper_bounds;
   }

   const Vec<ColFlags>&
   getProbingDomainFlags() const
   {
      return propagated_domain_flags;
   }

   void mark_all_activity_as_changed( ) {
      for( int i = 0; i < problem.getNRows( ); i++ )
         next_prop_activities.push_back(i);
   }

   void copyBounds( Problem<Rational>& problem );
};

void
ProbingView::activityChanged( ActivityChange actchange, int rowid,
                                    RowActivity<Rational>& activity )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhs = consMatrix.getLeftHandSides();
   const auto& rhs = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();

   // mark the lastchange fields with round values starting from -2
   // and counting backward By doing this we avoid that we need to
   // alter the lastchange field after copying the original activities
   // since they are always larger or equal to -1
   if( activity.lastchange > -2 )
      changed_activities.push_back( rowid ); // activity was changed for the first time

   if( activity.lastchange != round )
      next_prop_activities.push_back( rowid );

   activity.lastchange = round;

   // check if the updated activity is reliable or if it is zero relative to
   // the initial activity
   const RowActivity<Rational>& origactivity = problem.getRowActivities()[rowid];

   bool unreliable;

   if( actchange == ActivityChange::kMin )
      unreliable = ( activity.ninfmin <= 1 && activity.min != 0 &&
                     origactivity.min != 0 &&
                     num.isZero( activity.min / origactivity.min ) );
   else
      unreliable = ( activity.ninfmax <= 1 && activity.max != 0 &&
                     origactivity.max != 0 &&
                     num.isZero( activity.max / origactivity.max ) );

   if( unreliable )
   {
      auto rowvec = problem.getConstraintMatrix().getRowCoefficients( rowid );

      activity = compute_row_activity(rowvec.getValues(), rowvec.getIndices(),
                                      rowvec.getLength(), propagated_lower_bounds,
                                      propagated_upper_bounds,
                                      propagated_domain_flags, round );
   }

   // check for infeasibility
   if( actchange == ActivityChange::kMin && activity.ninfmin == 0 &&
       !rflags[rowid].test( RowFlag::kRhsInf ) &&
       num.isFeasLT( rhs[rowid], activity.min ) &&
       num.isSafeLT( rhs[rowid], activity.min ) )
   {
      Message::debug( this,
                      "[{}:{}] probing on col {} with val {} is infeasible min "
                      "activity is {:.15}, right hand side is {:.15}, and "
                      "original max activity was {:.15}\n",
                      __FILE__, __LINE__, probingCol, probingValue,
                      double( activity.min ), double( rhs[rowid] ),
                      double( problem.getRowActivities()[rowid].min ) );
      infeasible = true;
   }

   if( actchange == ActivityChange::kMax && activity.ninfmax == 0 &&
       !rflags[rowid].test( RowFlag::kLhsInf ) &&
       num.isFeasGT( lhs[rowid], activity.max ) &&
       num.isSafeGT( lhs[rowid], activity.max ) )
   {
      Message::debug( this,
                      "[{}:{}] probing on col {} with val {} is infeasible max "
                      "activity is {:.15}, left hand side is {:.15}, and "
                      "original max activity was {:.15}\n",
                      __FILE__, __LINE__, probingCol, probingValue,
                      double( activity.max ), double( lhs[rowid] ),
                      double( problem.getRowActivities()[rowid].max ) );
      infeasible = true;
   }
}

void
ProbingView::changeLb( int col, Rational& newlb, int row )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   auto colvec = consMatrix.getColumnCoefficients( col );

   bool lbinf = propagated_domain_flags[col].test(ColFlag::kLbUseless );

   // assumptions seem to be not as tight as necessary
   if(row == -1 && ( !lbinf && propagated_lower_bounds[col] >= newlb))
         return;

   Message::debug( this, "changing probing lower bound of col {} to {}\n", col,
                   double( newlb ) );

   assert(lbinf || propagated_lower_bounds[ col ] <= newlb);
   if( lbinf )
   {
      // bound was not altered yet, store the negative (index + 1) to
      // indicate that the infinity flag was altered
      propagated_domain_flags[col].unset(ColFlag::kLbUseless );
      changed_bounds.emplace_back( false, col, newlb, row );
   }
   else if( propagated_lower_bounds[col] <= newlb
//           && !problem.getColFlags()[col].test( ColFlag::kLbUseless )
            )
      // if bound was not altered yet remember it in the index vector
      changed_bounds.emplace_back( false, col, newlb, row );

   // change the bound in the domain vector
   Rational oldlb = propagated_lower_bounds[col];
   propagated_lower_bounds[col] = newlb;

   // update the probing activities by using the column view
   update_activities_after_boundchange(
       colvec.getValues(), colvec.getIndices(), colvec.getLength(),
       BoundChange::kLower, oldlb, newlb, lbinf, probing_activities,
       [this]( ActivityChange actChange, int rowid,
               RowActivity<Rational>& activity ) {
          activityChanged( actChange, rowid, activity );
       },
       true );
}

void
ProbingView::changeUb( int col, Rational& newub, int row )
{
   const auto& consMatrix = problem.getConstraintMatrix();
   auto colvec = consMatrix.getColumnCoefficients( col );

   // bound must be tighter than current domains
   bool ubinf = propagated_domain_flags[col].test(ColFlag::kUbUseless );

   if(row == -1 && ( !ubinf && propagated_upper_bounds[col] <= newub))
      return;

   Message::debug( this, "changing probing upper bound of col {} to {}\n", col,
                   double( newub ) );
   assert(ubinf || propagated_upper_bounds[ col ] >= newub);

   if( ubinf )
   {
      // bound was not altered yet, store the negative (index + 1) to
      // indicate that the infinity flag was altered
      propagated_domain_flags[col].unset(ColFlag::kUbUseless );
      changed_bounds.emplace_back( true, col, newub, row );
   }
   else if( propagated_upper_bounds[col] >= newub
//   && !problem.getColFlags()[col].test( ColFlag::kUbUseless )
   )
      // if bound was not altered yet remember it in the index vector
      changed_bounds.emplace_back( true, col, newub, row );


   // change the bound in the domain vector
   Rational oldub = propagated_upper_bounds[col];
   propagated_upper_bounds[col] = newub;

   // update the probing activities by using the column view
   update_activities_after_boundchange(
       colvec.getValues(), colvec.getIndices(), colvec.getLength(),
       BoundChange::kUpper, oldub, newub, ubinf, probing_activities,
       [this]( ActivityChange actChange, int rowid,
               RowActivity<Rational>& activity ) {
          activityChanged( actChange, rowid, activity );
       },
       true );
}

void
ProbingView::propagateDomains()
{
   const auto& consMatrix = problem.getConstraintMatrix();
   const auto& lhs = consMatrix.getLeftHandSides();
   const auto& rhs = consMatrix.getRightHandSides();
   const auto& rflags = consMatrix.getRowFlags();

   using std::swap;

   swap( prop_activities, next_prop_activities );
   next_prop_activities.clear();

   while( !prop_activities.empty() )
   {
      int curr_round = round--;

      Message::debug( this,
                      "starting probing propagation round {} on {} rows\n",
                      -curr_round - 2, prop_activities.size() );

      for( int candrow : prop_activities )
      {
         bool propagate = false;

         if( !rflags[candrow].test( RowFlag::kRhsInf ) && probing_activities[candrow].ninfmin <= 1 )
            propagate = true;

         if( !rflags[candrow].test( RowFlag::kLhsInf ) &&
             probing_activities[candrow].ninfmax <= 1 )
            propagate = true;

         if( !propagate )
            continue;

         auto rowvec = consMatrix.getRowCoefficients( candrow );

         propagate_row(candrow,
                       rowvec.getValues(), rowvec.getIndices(), rowvec.getLength(),
                       probing_activities[candrow], lhs[candrow], rhs[candrow],
                       rflags[candrow], propagated_lower_bounds, propagated_upper_bounds,
                       propagated_domain_flags,
                       [this]( BoundChange bndChg, int colid, Rational newbound , int row ) {
                if( num.isHugeVal( newbound ) )
                   return;

                bool isint = propagated_domain_flags[colid].test(
                    ColFlag::kIntegral, ColFlag::kImplInt );

                if( bndChg == BoundChange::kLower )
                {
                   if( isint )
                      newbound = num.feasCeil( newbound );

                   if( !propagated_domain_flags[colid].test(ColFlag::kUbInf ) &&
                       newbound > propagated_upper_bounds[colid] )
                   {
                      if( num.isFeasGT(newbound,
                                       propagated_upper_bounds[colid] ) )
                      {
                         Message::debug( this,
                                         "[{}:{}] probing on col {} with "
                                         "val {} is infeasible\n",
                                         __FILE__, __LINE__, probingCol,
                                         probingValue );
                         infeasible = true;
                         return;
                      }

                      newbound = propagated_upper_bounds[colid];
                   }

                   Rational delta = newbound - propagated_lower_bounds[colid];
                   bool finiteDomain =
                       !propagated_domain_flags[colid].test(ColFlag::kUbInf );

                   Rational mindomred = isint ? minintdomred : mincontdomred;

                   if( propagated_domain_flags[colid].test(
                           ColFlag::kLbUseless ) ||
                       ( finiteDomain && delta > 0 &&
                         ( delta / ( propagated_upper_bounds[colid] -
                                     propagated_lower_bounds[colid] ) >=
                           mindomred ) ) )
                      changeLb( colid, newbound, row );
                }
                else
                {
                   assert( bndChg == BoundChange::kUpper );
                   if( isint )
                      newbound = num.feasFloor( newbound );

                   if( !propagated_domain_flags[colid].test(ColFlag::kLbInf ) &&
                       newbound < propagated_lower_bounds[colid] )
                   {
                      if( num.isFeasLT(newbound,
                                       propagated_lower_bounds[colid] ) )
                      {
                         Message::debug( this,
                                         "[{}:{}] probing on col {} with "
                                         "val {} is infeasible\n",
                                         __FILE__, __LINE__, probingCol,
                                         probingValue );
                         infeasible = true;
                         return;
                      }

                      newbound = propagated_lower_bounds[colid];
                   }

                   Rational delta = propagated_upper_bounds[colid] - newbound;
                   bool finiteDomain =
                       !propagated_domain_flags[colid].test(ColFlag::kLbInf );

                   Rational mindomred = isint ? minintdomred : mincontdomred;

                   if( propagated_domain_flags[colid].test(
                           ColFlag::kUbUseless ) ||
                       ( finiteDomain && delta > 0 &&
                         ( delta / ( propagated_upper_bounds[colid] -
                                     propagated_lower_bounds[colid] ) >=
                           mindomred ) ) )
                      changeUb( colid, newbound, row );
                }
             } );
         if( infeasible )
            return;
      }

      swap( prop_activities, next_prop_activities );
      next_prop_activities.clear();
   }
}

   void ProbingView::copyBounds( Problem<Rational>& problem ) {
      for( int i = 0; i < problem.getNCols( ); i++ )
      {
         problem.getLowerBounds()[i] = propagated_lower_bounds[i];
         problem.getUpperBounds()[i] = propagated_upper_bounds[i];
         problem.getColFlags()[i] = propagated_domain_flags[i];
      }
   }

} // namespace papilo

#endif
