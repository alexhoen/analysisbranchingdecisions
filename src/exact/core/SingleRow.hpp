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

#ifndef _EXACT_CORE_SINGLE_ROW_HPP_
#define _EXACT_CORE_SINGLE_ROW_HPP_

#include "exact/data/RowFlags.hpp"
#include "exact/data/VariableDomains.hpp"
#include "exact/misc/Flags.hpp"
#include "exact/misc/Num.hpp"
#include "exact/misc/Vec.hpp"
#include <tuple>

namespace papilo
{

enum class BoundChange
{
   kLower,
   kUpper
};

enum class ActivityChange
{
   kMin,
   kMax
};

enum class RowStatus
{
   kInfeasible,
   kRedundant,
   kRedundantLhs,
   kRedundantRhs,
   kUnknown,
};

struct RowActivity
{
   /// minimal activity of the row
   Rational min;

   /// maximal activity of the row
   Rational max;

   /// number of variables that contribute with an infinite bound to the minimal
   /// activity of this row
   int ninfmin;

   /// number of variables that contribute with an infinite bound to the maximal
   /// activity of this row
   int ninfmax;

   /// last presolving round where this activity changed
   int lastchange;

   bool
   repropagate( ActivityChange actChange, RowFlags rflags )
   {
      if( actChange == ActivityChange::kMin &&
          !rflags.test( RowFlag::kRhsInf ) && ninfmin <= 1 )
         return true;

      if( actChange == ActivityChange::kMax &&
          !rflags.test( RowFlag::kLhsInf ) && ninfmax <= 1 )
         return true;

      return false;
   }

   RowStatus
   checkStatus( const Num<Rational>& num, RowFlags rflags, const Rational& lhs,
                const Rational& rhs ) const
   {
      RowStatus status = RowStatus::kRedundant;

      if( !rflags.test( RowFlag::kLhsInf ) )
      {
         if( ninfmax == 0 && num.isFeasLT( max, lhs ) &&
             num.isSafeLT( max, lhs ) )
            return RowStatus::kInfeasible;

         if( ninfmin == 0 && num.isFeasGE( min, lhs ) )
            status = RowStatus::kRedundantLhs;
         else
            status = RowStatus::kUnknown;
      }

      if( !rflags.test( RowFlag::kRhsInf ) )
      {
         if( ninfmin == 0 && num.isFeasGT( min, rhs ) &&
             num.isSafeGT( min, rhs ) )
            return RowStatus::kInfeasible;

         if( ninfmax == 0 && num.isFeasLE( max, rhs ) )
         {
            if( status == RowStatus::kUnknown )
               status = RowStatus::kRedundantRhs;
            else
               status = RowStatus::kRedundant;
         }
         else if( status == RowStatus::kRedundant )
            status = RowStatus::kUnknown;
      }
      else if( status == RowStatus::kRedundantLhs )
         status = RowStatus::kRedundant;

      return status;
   }

   template <typename Archive>
   void
   serialize( Archive& ar, const unsigned int version )
   {
      ar& min;
      ar& max;
      ar& ninfmin;

      ar& ninfmax;
      ar& lastchange;
   }

   RowActivity() = default;
};

/// counts the locks for the given row entry
void
count_locks( const Rational& val, RowFlags rflags, int& ndownlocks, int& nuplocks )
{
   assert( val != 0 );

   if( val < 0 )
   {
      if( !rflags.test( RowFlag::kLhsInf ) )
         ++nuplocks;

      if( !rflags.test( RowFlag::kRhsInf ) )
         ++ndownlocks;
   }
   else
   {
      if( !rflags.test( RowFlag::kLhsInf ) )
         ++ndownlocks;

      if( !rflags.test( RowFlag::kRhsInf ) )
         ++nuplocks;
   }
}

RowActivity<Rational>
compute_row_activity( const Rational* rowvals, const int* colindices, int rowlen,
                      const Vec<Rational>& lower_bounds,
                      const Vec<Rational>& upper_bounds, const Vec<ColFlags>& flags,
                      int presolveround = -1 )
{
   RowActivity<Rational> activity;

   activity.min = 0.0;
   activity.max = 0.0;
   activity.ninfmin = 0;
   activity.ninfmax = 0;
   activity.lastchange = presolveround;

   for( int j = 0; j < rowlen; ++j )
   {
      int col = colindices[j];
      if( !flags[col].test( ColFlag::kUbUseless ) )
      {
         if( rowvals[j] < 0 )
            activity.min += rowvals[j] * upper_bounds[col];
         else
            activity.max += rowvals[j] * upper_bounds[col];
      }
      else
      {
         assert( flags[col].test( ColFlag::kUbUseless ) );
         if( rowvals[j] < 0 )
            ++activity.ninfmin;
         else
            ++activity.ninfmax;
      }

      if( !flags[col].test( ColFlag::kLbUseless ) )
      {
         if( rowvals[j] < 0 )
            activity.max += rowvals[j] * lower_bounds[col];
         else
            activity.min += rowvals[j] * lower_bounds[col];
      }
      else
      {
         assert( flags[col].test( ColFlag::kLbUseless ) );
         if( rowvals[j] < 0 )
            ++activity.ninfmax;
         else
            ++activity.ninfmin;
      }
   }

   return activity;
}

Rational
compute_minimal_row_activity( const Rational* rowvals, const int* colindices, int rowlen,
                      const Vec<Rational>& lower_bounds,
                      const Vec<Rational>& upper_bounds, const Vec<ColFlags>& flags)
{
   Rational min = 0.0;

   for( int j = 0; j < rowlen; ++j )
   {
      int col = colindices[j];
      if( !flags[col].test( ColFlag::kUbUseless ) &&  rowvals[j] < 0 )
            min += rowvals[j] * upper_bounds[col];
      if( !flags[col].test( ColFlag::kLbUseless ) && rowvals[j] > 0 )
            min += rowvals[j] * lower_bounds[col];
   }
   return min ;
}

Rational
compute_maximal_row_activity( const Rational* rowvals, const int* colindices, int rowlen,
                      const Vec<Rational>& lower_bounds,
                      const Vec<Rational>& upper_bounds, const Vec<ColFlags>& flags)
{
   Rational max = 0.0;

   for( int j = 0; j < rowlen; ++j )
   {
      int col = colindices[j];
      if( !flags[col].test( ColFlag::kUbUseless ) && rowvals[j] > 0 )
            max += rowvals[j] * upper_bounds[col];
      if( !flags[col].test( ColFlag::kLbUseless ) && rowvals[j] < 0 )
            max += rowvals[j] * lower_bounds[col];
   }

   return max;
}

/// update the vector of row activities after lower or upper bounds of a column
/// changed. The last argument must be callable with arguments (ActivityChange,
/// rowid, RowActivity) and is called to inform about row activities that
/// changed
ActivityChange
update_activity_after_boundchange( const Rational& colval, BoundChange type,
                                   const Rational& oldbound, const Rational& newbound,
                                   bool oldbound_inf,
                                   RowActivity<Rational>& activity )
{
   assert( oldbound_inf ||
           ( type == BoundChange::kLower && newbound != oldbound ) ||
           ( type == BoundChange::kUpper && newbound != oldbound ) );

   if( type == BoundChange::kLower )
   {
      if( colval < Rational{ 0.0 } )
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmax > 0 );
            --activity.ninfmax;

            activity.max += newbound * colval;
         }
         else
         {
            activity.max += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::kMax;
      }
      else
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmin > 0 );
            --activity.ninfmin;

            activity.min += newbound * colval;
         }
         else
         {
            activity.min += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::kMin;
      }
   }
   else
   {
      if( colval < Rational{ 0.0 } )
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmin > 0 );
            --activity.ninfmin;

            activity.min += newbound * colval;
         }
         else
         {
            activity.min += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::kMin;
      }
      else
      {
         if( oldbound_inf )
         {
            assert( activity.ninfmax > 0 );
            --activity.ninfmax;

            activity.max += newbound * colval;
         }
         else
         {
            activity.max += ( newbound - oldbound ) * colval;
         }

         return ActivityChange::kMax;
      }
   }
}

/// update the vector of row activities after removing a finite lower or upper
/// bound of a column
void
update_activities_remove_finite_bound( const int* colinds, const Rational* colvals,
                                       int collen, BoundChange type,
                                       const Rational& oldbound,
                                       Vec<RowActivity<Rational>>& activities )
{
   if( type == BoundChange::kLower )
   {
      for( int i = 0; i != collen; ++i )
      {
         const Rational& colval = colvals[i];
         RowActivity<Rational>& activity = activities[colinds[i]];

         if( colval < Rational{ 0.0 } )
         {
            activity.max -= oldbound * colval;
            ++activity.ninfmax;
         }
         else
         {
            activity.min -= oldbound * colval;
            ++activity.ninfmin;
         }
      }
   }
   else
   {
      for( int i = 0; i != collen; ++i )
      {
         const Rational& colval = colvals[i];
         RowActivity<Rational>& activity = activities[colinds[i]];

         if( colval < Rational{ 0.0 } )
         {
            activity.min -= oldbound * colval;
            ++activity.ninfmin;
         }
         else
         {
            activity.max -= oldbound * colval;
            ++activity.ninfmax;
         }
      }
   }
}

/// update the vector of row activities after lower or upper bounds of a column
/// changed. The last argument must be callable with arguments (ActivityChange,
/// rowid, RowActivity) and is called to inform about row activities that
/// changed
template <typename ACTIVITYCHANGE>
void
update_activities_after_boundchange( const Rational* colvals, const int* colrows,
                                     int collen, BoundChange type,
                                     Rational oldbound, Rational newbound,
                                     bool oldbound_inf,
                                     Vec<RowActivity<Rational>>& activities,
                                     ACTIVITYCHANGE&& activityChange,
                                     bool watchInfiniteActivities = false )
{
   assert( oldbound_inf ||
           ( type == BoundChange::kLower && newbound != oldbound ) ||
           ( type == BoundChange::kUpper && newbound != oldbound ) );

   for( int i = 0; i < collen; ++i )
   {
      RowActivity<Rational>& activity = activities[colrows[i]];

      ActivityChange actChange = update_activity_after_boundchange(
          colvals[i], type, oldbound, newbound, oldbound_inf, activity );

      if( actChange == ActivityChange::kMin &&
          ( activity.ninfmin == 0 || watchInfiniteActivities ) )
         activityChange( ActivityChange::kMin, colrows[i], activity );

      if( actChange == ActivityChange::kMax &&
          ( activity.ninfmax == 0 || watchInfiniteActivities ) )
         activityChange( ActivityChange::kMax, colrows[i], activity );
   }
}

/**
 * updates the row activity for a changed coefficient in the matrix.
 * In case that the difference between the old and new coefficient is large,
 * the activity is recalculated entirely to prevent numerical difficulties.
 * The last argument must be callable with arguments (ActivityChange,
 * RowActivity) and is called to inform about row activities that changed.
 * @tparam Rational
 * @tparam ACTIVITYCHANGE
 * @param collb
 * @param colub
 * @param cflags
 * @param oldcolcoef
 * @param newcolcoef
 * @param activity
 * @param rowLength
 * @param colindices
 * @param rowvals
 * @param domains
 * @param num
 * @param activityChange
 */
template <typename ACTIVITYCHANGE>
void
update_activity_after_coeffchange( Rational collb, Rational colub, ColFlags cflags,
                                   Rational oldcolcoef, Rational newcolcoef,
                                   RowActivity<Rational>& activity,
                                   int rowLength, const int* colindices,
                                   const Rational* rowvals,
                                   const VariableDomains<Rational>& domains,
                                   const Num<Rational> num,
                                   ACTIVITYCHANGE&& activityChange )
{
   assert( oldcolcoef != newcolcoef );

   if( oldcolcoef * newcolcoef <= 0.0 )
   { // the sign of the coefficient flipped, so the column bounds now contribute
     // to the opposite activity bound

      // remember old activity
      RowActivity<Rational> oldactivity = activity;

      if( oldcolcoef != 0.0 )
      { // if the old coefficient was not 0.0 we remove its contributions to the
        // minimum and maximum activity
         // remove old contributions of the lower bound
         if( cflags.test( ColFlag::kLbUseless ) )
         {
            if( oldcolcoef < 0.0 )
               --activity.ninfmax;
            else
               --activity.ninfmin;
         }
         else
         {
            if( oldcolcoef < 0.0 )
               activity.max -= oldcolcoef * collb;
            else
               activity.min -= oldcolcoef * collb;
         }

         // remove old contributions of the upper bound
         if( cflags.test( ColFlag::kUbUseless ) )
         {
            if( oldcolcoef < 0.0 )
               --activity.ninfmin;
            else
               --activity.ninfmax;
         }
         else
         {
            if( oldcolcoef < 0.0 )
               activity.min -= oldcolcoef * colub;
            else
               activity.max -= oldcolcoef * colub;
         }
      }

      if( newcolcoef != 0.0 )
      { // if the new coefficient is not 0.0 we add its contributions to the
        // minimum and maximum activity
         // add new contributions of the lower bound
         if( cflags.test( ColFlag::kLbUseless ) )
         {
            if( newcolcoef < 0.0 )
               ++activity.ninfmax;
            else
               ++activity.ninfmin;
         }
         else
         {
            if( newcolcoef < 0.0 )
               activity.max += newcolcoef * collb;
            else
               activity.min += newcolcoef * collb;
         }

         // addnewold contributions of the upper bound
         if( cflags.test( ColFlag::kUbUseless ) )
         {
            if( newcolcoef < 0.0 )
               ++activity.ninfmin;
            else
               ++activity.ninfmax;
         }
         else
         {
            if( newcolcoef < 0.0 )
               activity.min += newcolcoef * colub;
            else
               activity.max += newcolcoef * colub;
         }
      }

      if( ( oldactivity.ninfmin != 0 && activity.ninfmin == 0 ) ||
          ( oldactivity.ninfmin == 0 && activity.ninfmin == 0 &&
            oldactivity.min != activity.min ) )
         activityChange( ActivityChange::kMin, activity );

      if( ( oldactivity.ninfmax != 0 && activity.ninfmax == 0 ) ||
          ( oldactivity.ninfmax == 0 && activity.ninfmax == 0 &&
            oldactivity.max != activity.max ) )
         activityChange( ActivityChange::kMax, activity );
   }
   else
   { // the sign of the coefficient did not flip, so the column bounds still
     // contribute to the same activity bound
      bool isDifferenceHugeVal = num.isHugeVal( newcolcoef - oldcolcoef );
      if( !cflags.test( ColFlag::kLbUseless ) && collb != 0.0 )
      {
         if( newcolcoef < Rational{ 0.0 } )
         {
            if( isDifferenceHugeVal )
               activity.max = compute_maximal_row_activity(
                   rowvals, colindices, rowLength, domains.lower_bounds,
                   domains.upper_bounds, domains.flags );
            else
               activity.max += collb * ( newcolcoef - oldcolcoef );

            if( activity.ninfmax == 0 )
               activityChange( ActivityChange::kMax, activity );
         }
         else
         {
            if( isDifferenceHugeVal )
               activity.min = compute_minimal_row_activity(
                   rowvals, colindices, rowLength, domains.lower_bounds,
                   domains.upper_bounds, domains.flags );

            else
               activity.min += collb * ( newcolcoef - oldcolcoef );
            if( activity.ninfmin == 0 )
               activityChange( ActivityChange::kMin, activity );
         }
      }

      if( !cflags.test( ColFlag::kUbUseless ) && colub != 0.0 )
      {
         if( newcolcoef < Rational{ 0.0 } )
         {
            if( isDifferenceHugeVal )
               activity.min = compute_minimal_row_activity(
                   rowvals, colindices, rowLength, domains.lower_bounds,
                   domains.upper_bounds, domains.flags );
            else
               activity.min += colub * ( newcolcoef - oldcolcoef );
            if( activity.ninfmin == 0 )
               activityChange( ActivityChange::kMin, activity );
         }
         else
         {
            if( isDifferenceHugeVal )
               activity.max = compute_maximal_row_activity(
                   rowvals, colindices, rowLength, domains.lower_bounds,
                   domains.upper_bounds, domains.flags );
            else
               activity.max += colub * ( newcolcoef - oldcolcoef );
            if( activity.ninfmax == 0 )
               activityChange( ActivityChange::kMax, activity );
         }
      }
   }
}

/// propagate domains of variables using the given a row and its activity. The
/// last argument must be callable with arguments (BoundChange, colid, newbound, row)
/// and is called to inform about column bounds that changed.
template <typename BOUNDCHANGE>
void
propagate_row( const Num<Rational>& num, int row, const Rational* rowvals, const int* colindices, int rowlen,
               const RowActivity<Rational>& activity, Rational lhs, Rational rhs,
               RowFlags rflags, const Vec<Rational>& lower_bounds,
               const Vec<Rational>& upper_bounds, const Vec<ColFlags>& domainFlags,
               BOUNDCHANGE&& boundchange )
{

   bool adj_rhs = false;
   if( activity.ninfmin == 1 && activity.ninfmax == 0 &&
       rflags.test( RowFlag::kRhsInf ) )
   {
      adj_rhs = true;
      rhs = activity.max;
   }

   if( ( !rflags.test( RowFlag::kRhsInf ) && activity.ninfmin <= 1 ) ||
       adj_rhs )
   {
      for( int j = 0; j < rowlen; ++j )
      {
         int col = colindices[j];
         Rational lb = lower_bounds[col];
         Rational ub = upper_bounds[col];
         Rational minresact = activity.min;
         Rational val = rowvals[j];

         if( val < Rational{ 0.0 } )
         {
            if( activity.ninfmin == 1 )
            {
               if( !domainFlags[col].test( ColFlag::kUbUseless ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::kUbUseless ) );
               minresact -= val * ub;
            }
            Rational new_lb = (rhs - minresact) / val;
            if( domainFlags[col].test(ColFlag::kIntegral) )
            {
               Rational floored = num.epsFloor(new_lb);
               if( num.isGE(rhs, floored * val + minresact) )
                  new_lb = floored;
               else
                  new_lb = num.epsCeil(new_lb);
            }
            if( domainFlags[col].test( ColFlag::kLbInf ) || new_lb > lb )
               boundchange( BoundChange::kLower, col, new_lb, row );
         }
         else
         {
            if( activity.ninfmin == 1 )
            {
               if( !domainFlags[col].test( ColFlag::kLbUseless ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::kLbUseless ) );
               minresact -= val * lb;
            }
            Rational new_ub = (rhs - minresact) / val;
            if( domainFlags[col].test(ColFlag::kIntegral) )
            {
               Rational ceiled = num.epsCeil(new_ub);
               if( num.isGE(rhs, ceiled * val + minresact) )
                  new_ub = ceiled;
               else
                  new_ub = num.epsFloor(new_ub);
            }
            if( domainFlags[col].test( ColFlag::kUbInf ) || new_ub < ub )
               boundchange( BoundChange::kUpper, col, new_ub, row );
         }
      }
   }

   bool adj_lhs = false;
   if( activity.ninfmax == 1 && activity.ninfmin == 0 &&
       rflags.test( RowFlag::kLhsInf ) )
   {
      adj_lhs = true;
      lhs = activity.min;
   }

   if( ( !rflags.test( RowFlag::kLhsInf ) && activity.ninfmax <= 1 ) ||
       adj_lhs )
   {
      for( int j = 0; j < rowlen; ++j )
      {
         int col = colindices[j];
         Rational lb = lower_bounds[col];
         Rational ub = upper_bounds[col];
         Rational maxresact = activity.max;
         Rational val = rowvals[j];

         if( val < Rational{ 0.0 } )
         {
            if( activity.ninfmax == 1 )
            {
               if( !domainFlags[col].test( ColFlag::kLbUseless ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::kLbUseless ) );
               maxresact -= val * lb;
            }
            Rational new_ub = (lhs - maxresact) / val;
            if( domainFlags[col].test(ColFlag::kIntegral) )
            {
               Rational ceiled = num.epsCeil(new_ub);
               if( num.isLE(lhs, ceiled * val + maxresact) )
                  new_ub = ceiled;
               else
                  new_ub = num.epsFloor(new_ub);
            }
            if( domainFlags[col].test( ColFlag::kUbInf ) || new_ub < ub )
               boundchange( BoundChange::kUpper, col, new_ub, row );

         }
         else
         {
            if( activity.ninfmax == 1 )
            {
               if( !domainFlags[col].test( ColFlag::kUbUseless ) )
                  continue;

               j = rowlen;
            }
            else
            {
               assert( !domainFlags[col].test( ColFlag::kUbUseless ) );
               maxresact -= val * ub;
            }
            Rational new_lb = (lhs - maxresact) / val;
            if( domainFlags[col].test(ColFlag::kIntegral) )
            {
               Rational floored = num.epsFloor(new_lb);
               if( num.isLE(lhs, floored * val + maxresact) )
                  new_lb = floored;
               else
                  new_lb = num.epsCeil(new_lb);
            }
            if( domainFlags[col].test( ColFlag::kLbInf ) || new_lb > lb )
               boundchange( BoundChange::kLower, col, new_lb, row );
         }
      }
   }
}

bool
row_implies_LB( const Num<Rational>& num, Rational lhs, Rational rhs, RowFlags rflags,
                const RowActivity<Rational>& activity, Rational colcoef, Rational collb,
                Rational colub, ColFlags cflags )

{
   if( cflags.test( ColFlag::kLbInf ) )
      return true;

   Rational resact;
   Rational side;

   if( colcoef > 0.0 && !rflags.test( RowFlag::kLhsInf ) )
   {
      if( activity.ninfmax == 0 )
      {
         assert( !cflags.test( ColFlag::kUbUseless ) );
         resact = activity.max - colub * colcoef;
      }
      else if( activity.ninfmax == 1 && cflags.test( ColFlag::kUbUseless ) )
         resact = activity.max;
      else
         return false;

      side = lhs;
   }
   else if( colcoef < 0.0 && !rflags.test( RowFlag::kRhsInf ) )
   {
      if( activity.ninfmin == 0 )
      {
         assert( !cflags.test( ColFlag::kUbUseless ) );
         resact = activity.min - colub * colcoef;
      }
      else if( activity.ninfmin == 1 && cflags.test( ColFlag::kUbUseless ) )
         resact = activity.min;
      else
         return false;

      side = rhs;
   }
   else
      return false;

   return num.isFeasGE( ( side - resact ) / colcoef, collb );
}

bool
row_implies_UB( const Num<Rational>& num, Rational lhs, Rational rhs, RowFlags rflags,
                const RowActivity<Rational>& activity, Rational colcoef, Rational collb,
                Rational colub, ColFlags cflags )
{
   if( cflags.test( ColFlag::kUbInf ) )
      return true;

   Rational resact;
   Rational side;

   if( colcoef > 0.0 && !rflags.test( RowFlag::kRhsInf ) )
   {
      if( activity.ninfmin == 0 )
      {
         assert( !cflags.test( ColFlag::kLbUseless ) );
         resact = activity.min - collb * colcoef;
      }
      else if( activity.ninfmin == 1 && cflags.test( ColFlag::kLbUseless ) )
         resact = activity.min;
      else
         return false;

      side = rhs;
   }
   else if( colcoef < 0.0 && !rflags.test( RowFlag::kLhsInf ) )
   {
      if( activity.ninfmax == 0 )
      {
         assert( !cflags.test( ColFlag::kLbUseless ) );
         resact = activity.max - collb * colcoef;
      }
      else if( activity.ninfmax == 1 && cflags.test( ColFlag::kLbUseless ) )
         resact = activity.max;
      else
         return false;

      side = lhs;
   }
   else
      return false;

   return num.isFeasLE( ( side - resact ) / colcoef, colub );
}

} // namespace papilo

#endif
