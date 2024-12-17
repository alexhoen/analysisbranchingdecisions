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

#ifndef _EXACT_SOLVER_HPP_
#define _EXACT_SOLVER_HPP_

#include "exact/misc/Vec.hpp"
#include "exact/core/PropagationView.hpp"
#include "exact/data/Tree.hpp"
#include "exact/data/Statistics.hpp"
#include "exact/data/Bound.hpp"
#include "exact/data/ExactOptions.hpp"

#include "tbb/concurrent_priority_queue.h"

namespace exact {


   class ExactSolver {

      Problem<Rational> &problem;
      Tree tree;

   public:
      ExactSolver(Problem<exact::Rational> &problem_) :
            problem(problem_) {
      }

      void start( ExactOptions& options) {

         double start_time = 0;
         Timer timer {start_time};

         exact::Message msg { };
         exact::Statistics stats{};


         bool unbounded_variables = false;
         for( int i = 0; i < problem.getNCols( ); ++i )
         {
            if( problem.getVariableDomains( ).flags[ i ].test(ColFlag::kUbInf) ||
                problem.getVariableDomains( ).flags[ i ].test(ColFlag::kLbInf))
            {
               unbounded_variables = true;
               break;
            }
         }
         if(unbounded_variables)
            msg.info("Before presolving problem has unbounded variables after presolving!\n");
         else
            msg.info("Before presolving problem has no unbounded variables after presolving!\n");

         //TODO: EXACT_HAVE_PAPILO not working
//         exact::PaPILOInterface presolve { };
//         presolve.presolve(problem);

         ViprInterface vipr { problem };
         Vec<std::pair<exact::Vec<exact::Rational>, exact::Rational>> solutions { };
         ScipInterface scip { };
         Bound best_bound {};
         Bound dual_bound {};

         problem.delete_unbounded_rows();

         vipr.init_vipr_certificate(problem);

         vipr.switch_to_tmp( false );
         vipr.invert_bounds_disguised_as_neg_constraints(problem);
         bool propagation_detected_infeasibility = false;
         if( options.propagation )
            propagation_detected_infeasibility = do_propagation(vipr);
         vipr.switch_to_original();

         #ifdef DEBUG
         for( int i = 0; i < problem.getNCols(); i++ )
            msg.info("{}: {} {}\n", problem.getVariableNames()[i], problem.getLowerBounds()[i], problem.getUpperBounds()[i]);
         #endif

         bool timelimit_reached = false;
         double seconds_in_events = 0;

         unbounded_variables = false;
         for( int i = 0; i < problem.getNCols( ); ++i )
         {
            if( problem.getVariableDomains( ).flags[ i ].test(ColFlag::kUbInf) ||
                problem.getVariableDomains( ).flags[ i ].test(ColFlag::kLbInf))
            {
               unbounded_variables = true;
               break;
            }
         }

         if(unbounded_variables)
            msg.info("Problem has unbounded variables after presolving!\n");
         else
            msg.info("Problem has no unbounded variables after presolving!\n");


         if(!propagation_detected_infeasibility)
         {
            vipr.switch_to_tmp(true);
            //setup SCIP and solve

            scip.setseed(options.randomseed);
            scip.doSetUp(problem, problem.getVariableDomains( ));
            scip.enable_branch_and_bound_only( options.feastol, options.randomseed );
            EventCatcher catcher { solutions, vipr, problem, best_bound, dual_bound, stats, msg, seconds_in_events };
            catcher.register_scip(scip);
            scip.set_time_limit(10800);
            auto status = scip.solve( );
            timelimit_reached = status == SCIP_STATUS_TIMELIMIT;
            if(timelimit_reached)
               dual_bound.update_neg_infinity();
            vipr.switch_to_original( );
         }

         if( solutions.empty( ) && !dual_bound.is_initialized())
            vipr.write_rtp_infeas( );
         else
            vipr.write_rtp_range(dual_bound, best_bound);

         vipr.write_solutions(solutions, true, best_bound.get());

         if(timelimit_reached)
            vipr.write_derivations(0);
         else
         {
            vipr.write_derivations(vipr.get_derivation_number( ));
            vipr.merge_tmp_to_certificate( );
         }
//#ifdef DEBUG
//         if(!timelimit_reached)
//            vipr.analysetree(problem.getVariableNames( ));
//#endif

         if(timelimit_reached)
            msg.info("\nTime limit reached in {}s spent {:.3f} in events  with bound [-inf,{}].\n", timer.getTime( ), seconds_in_events, best_bound.get());
         else if( solutions.empty( ))
         {
            if( !dual_bound.is_initialized())
               msg.info("\nInstance solved in {}s spent {:.3f} in events with status infeasible.\n", timer.getTime( ),
                     seconds_in_events);
            else
               msg.info("\nInstance solved with no feasible solution in {}s spent {:.3f} in events with bound [{}, inf].\n", timer.getTime( ),
                        seconds_in_events, dual_bound.get());
         }
         else
         {
            if(dual_bound.is_negative_infinity())
               msg.info("\nInstance solved in {}s spent {:.3f} in events with bound [-inf,{}].\n", timer.getTime( ), seconds_in_events, best_bound.get());
            else
               msg.info("\nInstance solved in {}s spent {:.3f} in events with bound [{},{}].\n", timer.getTime( ), seconds_in_events, dual_bound.get(), best_bound.get());
            msg.info("Exact-Gap: {} ({})\n", ( best_bound.get() - dual_bound.get() ).convert_to<double>( ),
                     ( best_bound.get() - dual_bound.get() ));
         }
         stats.print(msg);
      }

   private:

      bool do_propagation(ViprInterface &vipr) {
         problem.recomputeAllActivities( );
         ProbingView global_view { problem };
         global_view.mark_all_activity_as_changed( );
         global_view.propagateDomains( );
         remove_singleton_rows(global_view);
         global_view.copyBounds(problem);
         problem.recomputeAllActivities();

         Vec<Assumption> empty_assumptions { };
         for( auto &b: global_view.getProbingBoundChanges( ))
         {
            int vipr_index = vipr.write_boundchange(b, problem, empty_assumptions, vipr.getLbMapping( ),
                                                    vipr.getUbMapping( ));
            vipr.setVarIndex(b.var_index, vipr_index, b.upper);
         }
         if( !global_view.isInfeasible( ))
            return false;

         Vec<Rational> var_values(problem.getNCols(), 0);
         Vec<Rational> row_values(problem.getNRows(), 0);
         Vec<Assumption> assumptions{};

         auto constraintMatrix = problem.getConstraintMatrix( );
         for( int row = 0; row < problem.getNRows( ); row++ )
         {
            auto activity = problem.getRowActivities( )[ row ];
            if( activity.ninfmax == 0 && !constraintMatrix.getRowFlags( )[ row ].test(RowFlag::kLhsInf) &&
                constraintMatrix.getLeftHandSides( )[ row ] > activity.max )
            {
               row_values[row] = 1;
               auto data = constraintMatrix.getRowCoefficients(row);
               for( int i = 0; i < data.getLength(); ++i )
                  var_values[data.getIndices()[i]] = data.getValues()[i];
               vipr.write_node(NodeStatus::Infeasible, assumptions, var_values,
                               row_values, activity.max, data.getLength() +1, problem);
               return true;
            }
            if( activity.ninfmin == 0 && !constraintMatrix.getRowFlags( )[ row ].test(RowFlag::kRhsInf) &&
                constraintMatrix.getRightHandSides( )[ row ] < activity.min )
            {
               row_values[row] = -1;
               auto data = constraintMatrix.getRowCoefficients(row);
               for( int i = 0; i < data.getLength(); ++i )
                  var_values[data.getIndices()[i]] = data.getValues()[i];
               vipr.write_node(NodeStatus::Infeasible, assumptions, var_values,
                               row_values, activity.min, data.getLength() + 1, problem);
               return true;
            }
         }

         assert(false);
         return true;
      }

      void remove_singleton_rows( ProbingView &global_view) {
         // variables that are [lb, inf) in singleton rows might not be propagated
         for( int row = 0; row < this->problem.getNRows( ); row++ )
            if( this->problem.getRowSizes( )[ row ] == 1 )
            {
               auto vector = this->problem.getConstraintMatrix( ).getRowCoefficients(row);
               Rational coeff = vector.getValues( )[ 0 ];
               int col = vector.getIndices( )[ 0 ];
               if( coeff > 0 )
               {
                  if( !this->problem.getRowFlags( )[ row ].test(RowFlag::kLhsInf))
                  {
                     Rational newlb = this->problem.getConstraintMatrix( ).getLeftHandSides( )[ row ] / coeff;
                     if( global_view.getProbingDomainFlags( )[ col ].test(ColFlag::kLbInf) ||
                         global_view.getProbingLowerBounds( )[ col ] < newlb )
                        global_view.changeLb(col, newlb, row);

                  }
                  if( !this->problem.getRowFlags( )[ row ].test(RowFlag::kRhsInf))
                  {
                     Rational newUb = this->problem.getConstraintMatrix( ).getRightHandSides( )[ row ] / coeff;
                     if( global_view.getProbingDomainFlags( )[ col ].test(ColFlag::kUbInf) ||
                         global_view.getProbingUpperBounds( )[ col ] > newUb )
                        global_view.changeUb(col, newUb, row);
                  }
               }
               else if( coeff < 0 )
               {
                  if( !this->problem.getRowFlags( )[ row ].test(RowFlag::kLhsInf))
                  {
                     Rational newUb = this->problem.getConstraintMatrix( ).getLeftHandSides( )[ row ] / coeff;
                     if( global_view.getProbingDomainFlags( )[ col ].test(ColFlag::kUbInf) ||
                         global_view.getProbingUpperBounds( )[ col ] > newUb )
                        global_view.changeUb(col, newUb, row);
                  }
                  if( !this->problem.getRowFlags( )[ row ].test(RowFlag::kRhsInf))
                  {
                     Rational newlb = this->problem.getConstraintMatrix( ).getRightHandSides( )[ row ] / coeff;
                     if( global_view.getProbingDomainFlags( )[ col ].test(ColFlag::kLbInf) ||
                         global_view.getProbingLowerBounds( )[ col ] < newlb )
                        global_view.changeLb(col, newlb, row);
                  }
               }
            }
      }


   };
} // namespace exact

#endif
