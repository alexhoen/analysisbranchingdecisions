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

#ifndef EXACT_INTERFACES_EVENT_CATCHER_HPP_
#define EXACT_INTERFACES_EVENT_CATCHER_HPP_

#define UNUSED(expr) do { (void)(expr); } while (0)

//#define DEBUG
#include <cassert>
#include <stdexcept>

#include "scip/scip.h"
#include "scip/sol.h"

#include "exact/data/Solution.hpp"
#include "exact/data/Bound.hpp"
#include "exact/data/Statistics.hpp"
#include "exact/misc/Vec.hpp"
#include "exact/interfaces/ViprInterface.hpp"
#include "exact/interfaces/Leaf.hpp"
#include "exact/misc/Hash.hpp"
#include "exact/external/ska/bytell_hash_map.hpp"
#include "SoplexInterface.hpp"

namespace exact {

   static constexpr const char *const EVENT_NODE_ACTION = "node_action";

   static const Rational denom = Rational(0.000000001);

   static constexpr const char *const DESC_NODE_ACTION = "event handler for best solutions found";

   struct Exact_eventdata {
      Vec<std::pair<Vec<Rational>, Rational>> &solutions;
      ViprInterface &vipr;
      const Problem<Rational> &problem;
      Bound &best_bound;
      Bound &dual_bound;
      Statistics &stats;
      Message &msg;
      double& seconds_in_events;
      SoplexInterface exact_soplex;
      SoplexInterface floating_soplex;
      Vec<long long> parent_of_cutoff_nodes;
      int counter = 0;
      bool scip_found_solution= false;
      bool scaled_obj_integer = false;

   };

   class ReturnValue {

   public:

      bool success;
      Rational rhs;
      Vec<Rational> diff;
      int nnz;
      bool scaled_before_rounding;

      ReturnValue (bool _success, Rational& _rhs, Vec<Rational>& _diff, int _nnz)
         : success(_success), rhs(_rhs), diff(_diff), nnz(_nnz), scaled_before_rounding(false)
         {

         }

      ReturnValue (bool _success, Rational& _rhs, Vec<Rational>& _diff, int _nnz, bool _scaled_before_rounding )
            : success(_success), rhs(_rhs), diff(_diff), nnz(_nnz), scaled_before_rounding(_scaled_before_rounding)
      {

      }
   };


   class EventCatcher {
   private:

      Vec<std::pair<Vec<Rational>, Rational>> &solutions;
      ViprInterface &vipr;
      const Problem<Rational> &problem;
      Bound &best_bound;
      Bound &dual_bound;
      Statistics &stats;
      Message &msg;
      double &seconds_in_events;

   public:
      EventCatcher(Vec<std::pair<Vec<Rational>, Rational>> &_solutions, ViprInterface &_vipr,
                   const Problem<Rational> &_problem, Bound &_best_bound, Bound &_dual_bound, Statistics &_stats,
                   Message &_msg, double& _seconds_in_events) :
            solutions(_solutions), vipr(_vipr), problem(_problem),
            best_bound(_best_bound), dual_bound(_dual_bound), stats(_stats), msg(_msg), seconds_in_events(_seconds_in_events) {
      }

      SCIP_RETCODE
      register_scip(ScipInterface &scip, bool found_solution = false) {
         SoplexInterface exact_soplex{true};
         SoplexInterface floating_soplex{false};
//         exact_soplex.doSetUp(problem);
         floating_soplex.doSetUp(problem, problem.getVariableDomains());

         auto *eventdata = new Exact_eventdata { solutions, vipr, problem, best_bound, dual_bound, stats, msg,
                                                 seconds_in_events, exact_soplex, floating_soplex, { }};
         eventdata->scip_found_solution = found_solution;
         auto *eventhdlrdata = ( SCIP_EVENTHDLRDATA * ) eventdata;

         SCIP_EVENTHDLR *eventhdlr;

         eventhdlr = nullptr;
         SCIP_CALL(SCIPincludeEventhdlrBasic(scip.getSCIP( ), &eventhdlr, EVENT_NODE_ACTION, DESC_NODE_ACTION,
                                             eventNodeAction, eventhdlrdata));
         assert(eventhdlr != nullptr);

         SCIP_CALL(SCIPsetEventhdlrInit(scip.getSCIP( ), eventhdlr, eventInitBestsol));
         SCIP_CALL(SCIPsetEventhdlrExit(scip.getSCIP( ), eventhdlr, eventExitBestsol));
         return SCIP_OKAY;

      }

      static


      SCIP_DECL_EVENTINIT(eventInitBestsol) {  /*lint --e{715}*/
         assert(scip != nullptr);
         assert(eventhdlr != nullptr);
         assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENT_NODE_ACTION) == 0);

         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_POORSOLFOUND, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFEASIBLE, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEINFEASIBLE, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEDELETE, eventhdlr, nullptr, nullptr));

#ifdef DEBUG
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_DOMCHANGED, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_LBCHANGED, eventhdlr, nullptr, nullptr));
         SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_UBCHANGED, eventhdlr, nullptr, nullptr));

#endif

         return SCIP_OKAY;
      }

      static
      SCIP_DECL_EVENTEXIT(eventExitBestsol) {  /*lint --e{715}*/
         assert(scip != nullptr);
         assert(eventhdlr != nullptr);
         assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENT_NODE_ACTION) == 0);

         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_POORSOLFOUND, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFEASIBLE, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEINFEASIBLE, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEDELETE, eventhdlr, nullptr, -1));

#ifdef DEBUG
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEBRANCHED, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_DOMCHANGED, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_LBCHANGED, eventhdlr, nullptr, -1));
         SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_UBCHANGED, eventhdlr, nullptr, -1));
#endif

         return SCIP_OKAY;
      }



      static
      SCIP_DECL_EVENTEXEC(eventNodeAction) {
         assert(eventhdlr != nullptr);
         assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENT_NODE_ACTION) == 0);
         assert(event != nullptr);
         assert(scip != nullptr);
         auto *data = ( Exact_eventdata * ) SCIPeventhdlrGetData(eventhdlr);
         assert(data != nullptr);

         if( SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT )
            return SCIP_OKAY;

         double zero = 0;
         Timer t(zero);
         auto msg = data->msg;
         if(data->counter == 0)
         {
            double obj_scale = SCIPgetTransObjscale(scip);
            data->scaled_obj_integer = data->problem.isScaledObjectiveInteger(obj_scale);
            msg.info("Scale {} problem has integer obj {}, scaled has integer obj {}\n", obj_scale, data->problem.isObjectiveInteger(), data->scaled_obj_integer);
         }
         switch( SCIPeventGetType(event))
         {
            case SCIP_EVENTTYPE_NODEINFEASIBLE:
            {
               msg.info("------------\n");
               if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE )
               {
                  msg.info("node infeasible ({}) :", data->counter);
                  Vec<Assumption> assumptions = get_assumptions(scip, SCIPeventGetNode(event), msg);
                  SCIP_VAR **vars;
                  int nvars;
                  auto rows = SCIPgetLPRows(scip);
                  int nrows = SCIPgetNLPRows(scip);
                  SCIP_CALL(SCIPgetOrigVarsData(scip, &vars, &nvars, nullptr, nullptr, nullptr, nullptr));

                  Vec<Rational> var_farkas;
                  Vec<Rational> dual_farkas;
                  std::fill_n(std::back_inserter(var_farkas), nvars, 0);
                  std::fill_n(std::back_inserter(dual_farkas), data->problem.getNRows( ), 0);

                  VariableDomains<Rational> updated_domains = ( data->problem.getVariableDomains( ));
                  update_domains(updated_domains, assumptions, assumptions.size( ));

                  for( int var = 0; var < nvars; var++ )
                  {
                     int index = SCIPvarGetIndex(vars[ var ]);
                     var_farkas[ index ] = Rational(-SCIPgetVarFarkasCoef(scip, vars[ var ]));
                     if(( var_farkas[ index ] < 0 && updated_domains.flags[ index ].test(ColFlag::kUbInf)) ||
                        ( var_farkas[ index ] > 0 && updated_domains.flags[ index ].test(ColFlag::kLbInf)))
                        var_farkas[ index ] = 0;
                  }
                  for( int i = 0; i < nrows; i++ )
                  {
                     int index = SCIProwGetIndex(rows[ i ]);
                     dual_farkas[ index ] = Rational(SCIProwGetDualfarkas(rows[ i ]));
                     if(( dual_farkas[ index ] < 0 && data->problem.getRowFlags( )[ index ].test(RowFlag::kRhsInf))
                        || ( dual_farkas[ index ] > 0 && data->problem.getRowFlags( )[ index ].test(RowFlag::kLhsInf)))
                     {
                        dual_farkas[ index ] = 0;
                     }
                  }
                  ReturnValue result =
                        get_farkas_value_and_evaluate_farkas_proof(data, var_farkas, dual_farkas, updated_domains);
                  if( result.success )
                  {
                     msg.info("\tverified\n");
                     data->vipr.write_node(NodeStatus::Infeasible, assumptions, var_farkas,
                                           dual_farkas, result.rhs, result.nnz, data->problem);
                     data->stats.node_infeasible_easy( );
                     break;
                  }
                  auto diff = apply_Neumaier_Shcherbina(data, updated_domains, var_farkas, result.diff, result.nnz);
                  result.rhs += diff.second;
                  if( diff.first && result.rhs > 0 )
                  {
                     msg.info("\tverified by applying Neumaier Shcherbina\n");
                     data->vipr.write_node(NodeStatus::Infeasible, assumptions, var_farkas,
                                           dual_farkas, result.rhs, result.nnz, data->problem);
                     data->stats.node_infeasible_neumaier_shcherbina( );
                     break;
                  }

                  reconstructVector(var_farkas, denom);
                  reconstructVector(dual_farkas, denom);
                  result = get_farkas_value_and_evaluate_farkas_proof(data, var_farkas, dual_farkas, updated_domains);
                  if( result.success )
                  {
                     msg.info("\tverified after reconstruction.\n");
                     data->vipr.write_node(NodeStatus::Infeasible, assumptions, var_farkas,
                                           dual_farkas, result.rhs, result.nnz, data->problem);
                     data->stats.node_infeasible_reconstruct( );
                     break;
                  }
                  SoplexInterface exact_soplex{true};
                  exact_soplex.doSetUp(data->problem, updated_domains);
//                  data->exact_soplex.reset_to_bounds(updated_domains);

                  auto soplex_result = exact_soplex.solve(data->problem.getRowFlags( ), data->problem);

                  if( soplex_result.status == NodeStatus::Infeasible)
                  {
                     result = get_farkas_value_and_evaluate_farkas_proof(data, soplex_result.reduced,
                                                                         soplex_result.dual,
                                                                         updated_domains);
                     assert( result.success );
                        msg.info("\tNode infeasible verified by exact calculation\n");
                        data->stats.node_infeasible_exact( );
                        data->vipr.write_node(NodeStatus::Infeasible, assumptions, soplex_result.reduced,
                                              soplex_result.dual, result.rhs, result.nnz, data->problem);
                        break;
                  }
                  assert(soplex_result.status == NodeStatus::Solved);
                  Rational value = compute_dual(data, soplex_result.reduced, soplex_result.dual, updated_domains);

                  data->vipr.write_node(NodeStatus::ObjLimit, assumptions, soplex_result.reduced,
                                        soplex_result.dual, value,
                                        count_nonzeros(soplex_result.dual) + count_nonzeros(soplex_result.reduced),
                                        data->problem);
                  update_lower_bound(data, value);
                  msg.info("FAILED: Instead of infeasible the node has a lower bound of {}", value);
                  data->stats.node_infeasible_error( );


               }
               else if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED )
               {
                  msg.info("node not solved ({}) :", data->counter);
                  compute_pseudoobjective(data, SCIPeventGetNode(event), scip);
                  data->stats.node_objlimit_not_solved_easy( );
                  break;
               }
               else if(( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT && !data->scip_found_solution))
               {
                  msg.info("node objlimit pseudo ({}) :", data->counter);

                  Vec<Assumption> assumptions = get_assumptions(scip, SCIPeventGetNode(event), msg);

                  VariableDomains<Rational> updated_domains = ( data->problem.getVariableDomains( ));
                  update_domains(updated_domains, assumptions, assumptions.size( ));


                  data->floating_soplex.reset_to_bounds(updated_domains);

                  auto double_result = data->floating_soplex.solve(data->problem.getRowFlags( ), data->problem);

                  if( double_result.status != NodeStatus::LPError )
                  {
                     assert(double_result.status == NodeStatus::Infeasible);
                     ReturnValue result =
                           get_farkas_value_and_evaluate_farkas_proof(data, double_result.reduced, double_result.dual,
                                                                      updated_domains);
                     if( result.success )
                     {
                        msg.info("Node (objlimit by pseudoobjective) verified\n");
                        data->vipr.write_node(NodeStatus::Infeasible, assumptions, double_result.reduced,
                                              double_result.dual, result.rhs, result.nnz, data->problem);
                        data->stats.node_objlimit_pseudo_easy( );
                        break;
                     }
                     else
                     {
                        auto diff = apply_Neumaier_Shcherbina(data, updated_domains, double_result.reduced, result.diff,
                                                              result.nnz);
                        result.rhs += diff.second;
                        if( diff.first && result.rhs > 0 )
                        {
                           data->stats.node_objlimit_pseudo_neumaier( );
                           msg.info("Node (objlimit by pseudoobjective) verified by applying polishing\n");
                           data->vipr.write_node(NodeStatus::Infeasible, assumptions, double_result.reduced,
                                                 double_result.dual, result.rhs, result.nnz, data->problem);
                           break;
                        }
                     }
                  }
                  else
                  {
                     msg.info("\tFAILED: ERROR in floating-point SoPlex\n");
                     data->stats.floating_point_lp_error();
                  }

                  SoplexInterface exact_soplex{true};
                  exact_soplex.doSetUp(data->problem, updated_domains);

//                  exact_soplex.reset_to_bounds( updated_domains);

                  auto exact_soplex_result = exact_soplex.solve(data->problem.getRowFlags( ), data->problem);

                  if( exact_soplex_result.status == NodeStatus::Infeasible)
                  {
                     auto result = get_farkas_value_and_evaluate_farkas_proof(data, exact_soplex_result.reduced,
                                                                         exact_soplex_result.dual, updated_domains);
                     if( result.success )
                     {
                        msg.info("\tverified by SoPlex exact solve\n");
                        data->vipr.write_node(NodeStatus::Infeasible, assumptions, exact_soplex_result.reduced,
                                              exact_soplex_result.dual, result.rhs, result.nnz, data->problem);
                        data->stats.node_objlimit_pseudo_exact( );
                        break;
                     }
                  }
                  else
                  {
                     assert(exact_soplex_result.status == NodeStatus::LPError);
                     msg.info("FAILED: ERROR in Exact-SoPlex\n");
                     data->stats.exact_lp_error();
                     update_lower_bound_to_neg_inf(data);
                  }
                  msg.info("FAILED: not verified");
                  data->stats.node_objlimit_pseudo_error( );
                  break;
               }
               else
               {
                  assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT);
                  msg.info("node objlimit ({}) :", data->counter);
                  double obj_scale = (SCIPgetTransObjscale(scip));
                  bool is_obj_scaled = obj_scale != 1;

                  SCIP_LPI *lpi;
                  SCIP_RETCODE returncode = SCIPgetLPI(scip, &lpi);
                  assert(returncode == SCIP_OKAY);

                  SCIP_VAR **vars;
                  int nvars;
                  auto rows = SCIPgetLPRows(scip);
                  int nrows = SCIPgetNLPRows(scip);
                  SCIP_CALL(SCIPgetOrigVarsData(scip, &vars, &nvars, nullptr, nullptr, nullptr, nullptr));
                  Vec<Assumption> assumptions = get_assumptions(scip, SCIPeventGetNode(event), msg);

                  SCIP_Real objval;
                  SCIP_Real double_redcosts[nvars];
                  SCIP_Real double_primal[nvars];
                  SCIP_Real double_dual_solution[nrows];

                  SCIPlpiGetSol(lpi, &objval, double_primal, double_dual_solution, nullptr, double_redcosts);

                  Vec<Rational> redcosts;
                  Vec<Rational> dual_solution;
                  std::fill_n(std::back_inserter(redcosts), nvars, 0);
                  std::fill_n(std::back_inserter(dual_solution), data->problem.getNRows( ), 0);

                  VariableDomains<Rational> updated_domains = data->problem.getVariableDomains( );
                  update_domains(updated_domains, assumptions, assumptions.size( ));

                  Rational obj = 0;
                  const Vec<Rational> &obj_coeff = data->problem.getObjective( ).coefficients;
                  for( int var = 0; var < nvars; var++ )
                  {
                     int index = SCIPvarGetIndex(vars[ var ]);
                     redcosts[ index ] = Rational(double_redcosts[ var ]);
                     if(( redcosts[ index ] < 0 && updated_domains.flags[ index ].test(ColFlag::kUbInf)) ||
                        ( redcosts[ index ] > 0 && updated_domains.flags[ index ].test(ColFlag::kLbInf)))
                        redcosts[ index ] = 0;
                     obj += obj_coeff[ index ] * double_primal[ var ];
                  }

                  assert(abs(objval * obj_scale - obj.convert_to<double>( )) <= 0.0000001 ||
                         ( objval * obj_scale / obj.convert_to<double>( ) <= 1.0000001 &&
                           objval * obj_scale / obj.convert_to<double>( ) >= 0.9999999 ));

                  for( int row = 0; row < nrows; row++ )
                  {
                     int row_index = SCIProwGetIndex(rows[ row ]);
                     dual_solution[ row_index ] = double_dual_solution[ row ];
                     if(( dual_solution[ row_index ] > 0 &&
                          data->problem.getRowFlags( )[ row_index ].test(RowFlag::kLhsInf)) ||
                        ( dual_solution[ row_index ] < 0 &&
                          data->problem.getRowFlags( )[ row_index ].test(RowFlag::kRhsInf)))
                     {
#ifdef DEBUG
                        msg.info("setting dual solution value of row {} to 0 ({})\n", row_index, dual_solution[ row_index ])
#endif
                        dual_solution[ row_index ] = 0;
                     }

                  }

                  auto result = get_dualbound_and_check_for_obj_limit(data, redcosts, dual_solution,
                                                                      updated_domains, scip, obj_scale);

                  if( result.success )
                  {
                     msg.info("\tverified\n");
                     ScaleInformation scaleInformation{ is_obj_scaled, obj_scale };
                     if(result.scaled_before_rounding)
                     {
                        Rational number = result.rhs/obj_scale;
                        ceil_rational(number);
                        scaleInformation = { is_obj_scaled, obj_scale, number };
                     }
                     data->vipr.write_node(NodeStatus::ObjLimit, assumptions, redcosts,
                                           dual_solution, result.rhs, result.nnz, data->problem, scaleInformation );
                     data->stats.node_objlimit_easy( );
                     break;
                  }

                  Vec<Rational> neumaier_redcosts(redcosts);

                  auto diff = apply_Neumaier_Shcherbina(data, updated_domains, neumaier_redcosts, result.diff,
                                                        result.nnz);

                  if( diff.first )
                  {
                     Rational sum = recalculate(data, neumaier_redcosts, dual_solution, updated_domains);
                     Rational scaled_obj = sum * obj_scale;
                     auto comp_result = compare(data->problem.isObjectiveInteger(), data->scaled_obj_integer, obj_scale, data->best_bound, scaled_obj);
                     if( comp_result == 1 )
                     {
                        if(data->problem.isObjectiveInteger())
                           ceil_rational(scaled_obj);
                        msg.info("\tverified by applying Neumaier Shcherbina\n");
                        data->vipr.write_node(NodeStatus::ObjLimit, assumptions, neumaier_redcosts, dual_solution,
                                              scaled_obj, result.nnz, data->problem, { is_obj_scaled, obj_scale });
                        data->stats.node_objlimit_neumaier_shcherbina( );
                        break;
                     }
                     if( comp_result == 2 )
                     {
                        if(data->problem.isObjectiveInteger())
                           ceil_rational(sum);
                        msg.info("\tverified by applying Neumaier Shcherbina\n");
                        data->vipr.write_node(NodeStatus::ObjLimit, assumptions, neumaier_redcosts, dual_solution,
                                              scaled_obj, result.nnz, data->problem,
                                              { is_obj_scaled, obj_scale, sum });
                        data->stats.node_objlimit_neumaier_shcherbina( );
                        break;
                     }
                     msg.info("\tObjLimit - Safe bounding failed {} !>= {}\n", scaled_obj, obj);
                  }


                  reconstructVector(redcosts, Rational(denom));
                  reconstructVector(dual_solution, Rational(denom));
                  result = get_dualbound_and_check_for_obj_limit(data, redcosts, dual_solution,
                                                                 updated_domains, scip, obj_scale);

                  if( result.success )
                  {
                     msg.info("\treconstruct is successful\n");
                     ScaleInformation scaleInformation{ is_obj_scaled, obj_scale };
                     if(result.scaled_before_rounding)
                     {
                        Rational number = result.rhs/obj_scale;
                        ceil_rational(number);
                        scaleInformation = { is_obj_scaled, obj_scale, number };
                     }
                     data->vipr.write_node(NodeStatus::ObjLimit, assumptions, redcosts, dual_solution,
                                           result.rhs, result.nnz, data->problem, scaleInformation);
                     data->stats.node_objlimit_reconstruct( );
                     break;
                  }

                  SoplexInterface exact_soplex{true};
                  exact_soplex.doSetUp(data->problem, updated_domains);
//                  exact_soplex.reset_to_bounds(updated_domains);
                  auto basis = convert_basis(vars, nvars, rows, data->problem.getNRows( ), scip);
                  exact_soplex.setBasis(basis.first, basis.second);

                  Vec<Rational> primal(redcosts.size( ));
                  auto success = exact_soplex.factorize(primal, dual_solution, redcosts, basis.first, basis.second, false);

                  if( success.first )
                  {
                     result = get_dualbound_and_check_for_obj_limit(data, redcosts, dual_solution,
                                                                    updated_domains, scip, obj_scale);
                     if( result.success )
                     {
                        msg.info("\tfactorize is successful\n");
                        ScaleInformation scaleInformation{ is_obj_scaled, obj_scale };
                        if(result.scaled_before_rounding)
                        {
                           Rational number = result.rhs/obj_scale;
                           ceil_rational(number);
                           scaleInformation = { is_obj_scaled, obj_scale, number };
                        }
                        data->vipr.write_node(NodeStatus::ObjLimit, assumptions, redcosts,
                                              dual_solution, result.rhs, result.nnz, data->problem, scaleInformation);
                        data->stats.node_objlimit_factorize( );
                        break;
                     }
                  }

                  auto exact_soplex_result = exact_soplex.solve(data->problem.getRowFlags( ), data->problem);

                  if( exact_soplex_result.status == NodeStatus::Infeasible )
                  {
                     result = get_farkas_value_and_evaluate_farkas_proof(data, exact_soplex_result.reduced,
                                                                         exact_soplex_result.dual, updated_domains);
                     assert(result.success);
                     data->stats.node_objlimit_exact( );
                     ScaleInformation scaleInformation{ is_obj_scaled, obj_scale };
                     if(result.scaled_before_rounding)
                     {
                        Rational number = result.rhs/obj_scale;
                        ceil_rational(number);
                        scaleInformation = { false, obj_scale, number };
                     }
                     msg.info("\tNode is infeasible and verified\n");
                     data->vipr.write_node(NodeStatus::Infeasible, assumptions,
                                           exact_soplex_result.reduced, exact_soplex_result.dual, result.rhs, result.nnz,
                                           data->problem, scaleInformation);
                     break;
                  }
                  else if( exact_soplex_result.status == NodeStatus::Solved)
                  {
                     result = get_dualbound_and_check_for_obj_limit(data, exact_soplex_result.reduced, exact_soplex_result.dual,
                                                                    updated_domains, scip, 1);

                     if( result.success )
                     {
                        data->stats.node_objlimit_exact( );
                        msg.info("\texact recalculation verified {} {}\n", result.rhs.convert_to<double>( ),
                                 data->dual_bound.get( ).convert_to<double>( ));
                        data->vipr.write_node(NodeStatus::ObjLimit, assumptions,
                                              exact_soplex_result.reduced, exact_soplex_result.dual, result.rhs, result.nnz,
                                              data->problem);
                        break;
                     }

                     data->stats.node_objlimit_error( );
                     msg.info("\thas better objective {} than current best solution {}\n", result.rhs,
                              data->dual_bound.get( ));
                     bool feasible = check_feasibility(data, exact_soplex_result.primal, false);
                     if( feasible )
                     {
                        assert(result.rhs == calculate_primal_obj(data->problem, exact_soplex_result.primal));
                        auto result_comp = compare_dual_bound_to_primal(data, exact_soplex_result.reduced, exact_soplex_result.dual,
                                                                        updated_domains, scip,
                                                                        result.rhs);
                        if( result_comp.success )
                        {
                           msg.info("\tObjLimit (Feasible solution) verified with obj {} - new solution added\n", result.rhs);
                           add_new_solution(data, exact_soplex_result.primal, result.rhs);
                           update_lower_bound(data, result.rhs);
                           data->vipr.write_node(NodeStatus::Solved, assumptions, redcosts,
                                                 dual_solution, result_comp.rhs, result_comp.nnz, data->problem);
                           break;
                        }

                     }
                     msg.info("FAILED: not verified {} {} (diff: {})\n", result.rhs,
                              data->dual_bound.get( ), result.rhs - data->dual_bound.get( ));
                     msg.info("\t\t\tnot verified {} {}\n", result.rhs.convert_to<double>(),
                              data->dual_bound.get( ).convert_to<double>());
                     data->vipr.write_node(NodeStatus::ObjLimit, assumptions, exact_soplex_result.reduced,
                                           exact_soplex_result.dual, result.rhs, result.nnz, data->problem);
                     update_lower_bound(data, result.rhs);
                  }
                  else{

                  }
               }
               break;
            }
            case SCIP_EVENTTYPE_NODEFEASIBLE:
            {
               msg.info("------------\n");
               msg.info("node feasible ({}) :", data->counter);
               SCIP_VAR **vars;
               int nvars;
               auto rows = SCIPgetLPRows(scip);
               int nrows = SCIPgetNLPRows(scip);
               SCIP_CALL(SCIPgetOrigVarsData(scip, &vars, &nvars, nullptr, nullptr, nullptr, nullptr));
               Vec<Assumption> assumptions = get_assumptions(scip, SCIPeventGetNode(event), msg);
               SCIP_SOL *bestsol = SCIPgetBestSol(scip);
               assert(bestsol != nullptr);
               Vec<Rational> primal, redcosts, dual_solution;
               std::fill_n(std::back_inserter(redcosts), nvars, 0);
               std::fill_n(std::back_inserter(primal), nvars, 0);
               std::fill_n(std::back_inserter(dual_solution), data->problem.getNRows( ), 0);

               for( int var = 0; var < nvars; var++ )
               {
                  int index = SCIPvarGetIndex(vars[ var ]);
                  redcosts[ index ] = Rational(SCIPgetVarRedcost(scip, vars[ var ]));
                  primal[ index ] = Rational { SCIPgetSolVal(scip, bestsol, vars[ var ]) };
               }
               for( int i = 0; i < nrows; i++ )
               {
                  int index = SCIProwGetIndex(rows[ i ]);
                  dual_solution[ index ] = Rational(SCIProwGetDualsol(rows[ i ]));
                  if(( dual_solution[ index ] < 0 && data->problem.getRowFlags( )[ index ].test(RowFlag::kRhsInf))
                     || ( dual_solution[ index ] > 0 && data->problem.getRowFlags( )[ index ].test(RowFlag::kLhsInf)))
                  {
                     assert(dual_solution[ index ] < 0.1);
                     dual_solution[ index ] = 0;
                  }
               }

               VariableDomains<Rational> updated_domains = ( data->problem.getVariableDomains( ));
               update_domains(updated_domains, assumptions, assumptions.size( ));


               bool feasible = check_feasibility(data, primal, true);
               if( feasible )
               {
                  Rational primal_obj = calculate_primal_obj(data->problem, primal);
                  auto result = compare_dual_bound_to_primal(data, redcosts, dual_solution, updated_domains, scip,
                                                             primal_obj);
                  if( result.success )
                  {
                     data->stats.node_feasible_easy( );
                     msg.info("\tverified with solution value {}\n", result.rhs);
                     add_new_solution(data, primal, primal_obj);
                     update_lower_bound(data, primal_obj);
                     data->vipr.write_node(NodeStatus::Solved, assumptions, redcosts,
                                           dual_solution, result.rhs, result.nnz, data->problem);
                     break;
                  }

               }


               reconstructVector(primal, Rational(denom));
               reconstructVector(redcosts, Rational(denom));
               reconstructVector(dual_solution, Rational(denom));

               Rational primal_obj = calculate_primal_obj(data->problem, primal);
               auto result = compare_dual_bound_to_primal(data, redcosts, dual_solution, updated_domains, scip,
                                                          primal_obj);

               if( check_feasibility(data, primal, true) && result.success )
               {
                  assert(check_feasibility(data, primal, false));
                  add_new_solution(data, primal, primal_obj);
                  update_lower_bound(data, primal_obj);
                  data->stats.node_feasible_reconstruction( );
                  msg.info("\tverified after reconstruction with solution value {}\n", primal_obj);
                  data->vipr.write_node(NodeStatus::Solved, assumptions, redcosts,
                                        dual_solution, result.rhs, result.nnz, data->problem);
                  break;
               }

               SoplexInterface exact_soplex{true};
               exact_soplex.doSetUp(data->problem, updated_domains);
//               exact_soplex.reset_to_bounds( updated_domains);
               auto basis = convert_basis(vars, nvars, rows, data->problem.getNRows( ), scip);
               exact_soplex.setBasis(basis.first, basis.second);

               auto success = exact_soplex.factorize(primal, dual_solution, redcosts, basis.first, basis.second, true);
               if( success.first && check_integrality(data, primal) )
               {
                  primal_obj = calculate_primal_obj(data->problem, primal);
                  assert( compare_dual_bound_to_primal(data, redcosts, dual_solution,
                                                        updated_domains, scip, primal_obj).success);
                  add_new_solution(data, primal, primal_obj);
                  update_lower_bound(data, primal_obj);
                  data->stats.node_feasible_factorization( );
                  msg.info("\tverified after factorization with solution value {}\n", primal_obj);
                  data->vipr.write_node(NodeStatus::Solved, assumptions, redcosts,
                                        dual_solution, primal_obj, success.second, data->problem);
                  break;
               }

               exact_soplex.setBasis(basis.first, basis.second);
               auto exact_soplex_result = exact_soplex.solve(data->problem.getRowFlags( ), data->problem);

               if( exact_soplex_result.status == NodeStatus::Solved)
               {
                  bool integrality = check_integrality(data, exact_soplex_result.primal);
                  if( integrality )
                  {
                     assert(check_feasibility(data, exact_soplex_result.primal, false));
                     primal_obj = calculate_primal_obj(data->problem, exact_soplex_result.primal);
                     add_new_solution(data, primal, primal_obj);
                     result = compare_dual_bound_to_primal(data, exact_soplex_result.reduced, exact_soplex_result.dual,
                                                           updated_domains, scip, primal_obj);
                     if( result.success )
                     {
                        data->stats.node_feasible_exact_solve( );
                        update_lower_bound(data, primal_obj);
                        msg.info("\tverified after exact-solve with solution value {}\n", primal_obj);
                        data->vipr.write_node(NodeStatus::Solved, assumptions, exact_soplex_result.reduced,
                                              exact_soplex_result.dual, result.rhs, result.nnz, data->problem);
                        break;
                     }

                  }
                  const Rational &dual_bound = compute_dual(data, exact_soplex_result.reduced, exact_soplex_result.dual,
                                                            updated_domains);
                  data->vipr.write_node(NodeStatus::Solved, assumptions, exact_soplex_result.reduced, exact_soplex_result.dual,
                                        dual_bound, count_nonzeros(exact_soplex_result.reduced) + count_nonzeros(exact_soplex_result.dual),
                                        data->problem);
                  if( !data->solutions.empty() && data->best_bound.is_initialized() && dual_bound > data->best_bound.get() )
                  {
                     msg.info("\t\tPrimal bound passed {}, {}\n", primal_obj.convert_to<double>(),  data->best_bound.get().convert_to<double>());
                     break;
                  }
                  if( !data->solutions.empty() && data->dual_bound.is_initialized() && dual_bound > data->dual_bound.get() )
                  {
                     data->stats.node_feasible_exact_solve( );
                     msg.info("FAILED: Dual bound exceeded {}, {}\n", primal_obj,  data->best_bound.get());
                     msg.info("\t\tDual bound exceeded {}, {}\n", primal_obj.convert_to<double>(),  data->best_bound.get().convert_to<double>());
                     break;
                  }
                  update_lower_bound(data, dual_bound);
                  msg.info("FAILED: NODE_FEASIBLE not verified with primal {} current primal {} integral {}\n", primal_obj, data->best_bound.get(), integrality);

               }
               else if( exact_soplex_result.status == NodeStatus::Infeasible)
               {
                  auto result = get_farkas_value_and_evaluate_farkas_proof(data, exact_soplex_result.reduced, exact_soplex_result.dual,
                                                                           updated_domains);
                  assert(result.success);
                  data->vipr.write_node(NodeStatus::Infeasible, assumptions, exact_soplex_result.reduced,
                                        exact_soplex_result.dual, result.rhs, result.nnz, data->problem);
                  msg.info("FAILED: NODE_FEASIBLE infeasible instead of feasible with {}\n", primal_obj);
               }
               else
               {
                  assert(exact_soplex_result.status == NodeStatus::LPError);
                  msg.info("FAILED: ERROR in Exact-SoPlex\n");
                  data->stats.exact_lp_error();
                  update_lower_bound_to_neg_inf(data);
               }
               data->stats.node_feasible_error( );

               break;
            }
            case SCIP_EVENTTYPE_POORSOLFOUND:
               assert(false);
            case SCIP_EVENTTYPE_BESTSOLFOUND:
            {
               msg.info("------------\n");
               msg.info("Best Solution found ({}) :\n", data->counter);
               data->scip_found_solution = true;
               SCIP_SOL *bestsol = SCIPgetBestSol(scip);
               assert(bestsol != nullptr);
               SCIP_VAR **vars;
               int nvars;
               SCIP_CALL(SCIPgetOrigVarsData(scip, &vars, &nvars, nullptr, nullptr, nullptr, nullptr));
               Vec<Rational> sol;
               sol.resize(nvars);
               for( int var = 0; var < nvars; var++ )
               {
                  int index = SCIPvarGetIndex(vars[ var ]);
                  sol[ index ] = Rational { SCIPgetSolVal(scip, bestsol, vars[ var ]) };
               }
               bool feasible = check_feasibility(data, sol, true);
               if( feasible )
               {
                  Rational primal_obj = calculate_primal_obj(data->problem, sol);
                  add_new_solution(data, sol, primal_obj);
                  msg.info("\tfeasible with solution value {}\n", primal_obj);
                  data->stats.solution_feasible( );
                  break;
               }
               reconstructVector(sol, Rational(denom));
               feasible = check_feasibility(data, sol, true);
               if( feasible )
               {
                  Rational primal_obj = calculate_primal_obj(data->problem, sol);
                  add_new_solution(data, sol, primal_obj);
                  msg.info("\tfeasible after reconstruction with solution value {}\n", primal_obj);
                  data->stats.solution_reconstruction_feasible( );
                  break;
               }

               SoplexInterface exact_soplex { true };

               auto updated_domains = data->problem.getVariableDomains();

               for(int i=0; i < data->problem.getNCols(); i++)
                  if(data->problem.getColFlags()[i].test(ColFlag::kIntegral))
                  {
                     Rational value = floor(sol[i]+0.5);
                     updated_domains.lower_bounds[i] = value;
                     updated_domains.upper_bounds[i] = value;
                  }

               exact_soplex.doSetUp(data->problem, updated_domains);

               auto exact_soplex_result = exact_soplex.solve(data->problem.getRowFlags( ), data->problem);

               if ( exact_soplex_result.status == NodeStatus::Solved)
               {
                  if( check_feasibility(data, exact_soplex_result.primal, false))
                  {
                     assert(check_feasibility(data, exact_soplex_result.primal, false));
                     Rational primal_obj = calculate_primal_obj(data->problem, exact_soplex_result.primal);
                     add_new_solution(data, exact_soplex_result.primal, primal_obj);
                     data->stats.solution_integer_fixing_feasible( );
                     msg.info("\tfeasible after integer-fixing with solution value {}\n", primal_obj);
                     break;
                  }
               }
               else if(exact_soplex_result.status == NodeStatus::LPError)
               {
                  msg.info("FAILED: ERROR in Exact-SoPlex\n");
                  data->stats.exact_lp_error();
               }
               add_new_solution(data, sol, calculate_primal_obj(data->problem, sol));
               data->stats.solution_error( );
               msg.info("FAILED: found new best solution with solution value {} could not be verified\n",
                        SCIPgetSolOrigObj(scip, bestsol));
               break;
            }
            case SCIP_EVENTTYPE_NODEDELETE:
            {
               auto node = SCIPeventGetNode(event);
               long long node_number = SCIPnodeGetNumber(node);
               SCIP_NODETYPE type = SCIPnodeGetType(node);
               auto parent_node = SCIPnodeGetParent(node);
               if( type == SCIP_NODETYPE_LEAF || type == SCIP_NODETYPE_SIBLING || type == SCIP_NODETYPE_CHILD )
               {
                  bool parent = type != SCIP_NODETYPE_CHILD;

                  msg.info("------------\n");
                  // TODO: use depth instead
                  double lb = SCIPnodeGetLowerbound(node);
                  double lb_parent = SCIPnodeGetLowerbound(parent_node);
                  if( !SCIPisInfinity(scip, lb) && lb != lb_parent )
                  {
                     msg.info("Leaf (scip node #{}): node deleted pseudo ({}) :", node_number, data->counter);
                     compute_pseudoobjective(data, node, scip);
                     data->stats.node_deletion_pseudo_easy( );
                     break;
                  }
                  if( parent && lb != lb_parent && lb_parent < data->dual_bound.get() )
                     parent = false;
                  NodeStatus status = parent ? NodeStatus::ObjLimitParent : NodeStatus::ObjLimit;
                  msg.info("Leaf node cutoff ({}) (scip node #{}) (Status: {}/{}):\n", data->counter, node_number, type, parent);

                  if( parent )
                  {
                     long long parent_number = SCIPnodeGetNumber(SCIPnodeGetParent(node));
                     if( std::find(data->parent_of_cutoff_nodes.begin( ), data->parent_of_cutoff_nodes.end( ),
                                   parent_number) != data->parent_of_cutoff_nodes.end( ))
                     {
                        msg.info("\tcutoff leaf already registered ");
                        const Vec<Assumption> assumptions = get_assumptions(scip, node, msg);
                        break;
                     }
                     data->parent_of_cutoff_nodes.push_back(parent_number);
                  }

                  Vec<Assumption> assumptions = get_assumptions(scip, node, msg);

                  VariableDomains<Rational> updated_domains = ( data->problem.getVariableDomains( ));
                  update_domains(updated_domains, assumptions, assumptions.size( ) - ( parent ? 1 : 0 ));

//                  if( data->solutions.empty() )
//                  {
//                     ScipInterface subscip{};
//                     subscip.doSetUp(data->problem, updated_domains);
//                     subscip.enable_branch_and_bound_only( );
//                     double time = 0;
//                     EventCatcher catcher { data->solutions, data->vipr, data->problem, data->best_bound, data->dual_bound, data->stats, data->msg, time };
//                     catcher.register_scip(subscip, true);
//                     subscip.set_time_limit(10800);
//                     data->vipr.load_assumptions(assumptions);
//                     auto status = subscip.solve( );
//                     data->vipr.reset_assumptions();
//                     break;
//                  }


                  data->floating_soplex.reset_to_bounds(updated_domains);

                  double obj_scale = SCIPgetTransObjscale(scip);
                  ScaleInformation scaleInformation;

                  auto soplex_result = data->floating_soplex.solve(data->problem.getRowFlags( ), data->problem);

                  for( unsigned int var = 0; var < soplex_result.reduced.size( ); var++ )
                     if(( soplex_result.reduced[ var ] < 0 && updated_domains.flags[ var ].test(ColFlag::kUbInf)) ||
                        ( soplex_result.reduced[ var ] > 0 && updated_domains.flags[ var ].test(ColFlag::kLbInf)))
                        soplex_result.reduced[ var ] = 0;
                  for( unsigned int row = 0; row < soplex_result.dual.size( ); row++ )
                     if(( soplex_result.dual[ row ] < 0 && data->problem.getRowFlags( )[ row ].test(RowFlag::kRhsInf))
                        ||
                        ( soplex_result.dual[ row ] > 0 && data->problem.getRowFlags( )[ row ].test(RowFlag::kLhsInf)))
                        soplex_result.dual[ row ] = 0;

                  if( soplex_result.status == NodeStatus::Solved )
                  {
                     auto result = get_dualbound_and_check_for_obj_limit(data, soplex_result.reduced,
                                                                         soplex_result.dual, updated_domains, scip, 1);

                     if( result.success )
                     {
                        data->stats.node_deletion_soplex( );
                        msg.info("\trecalculating verified by solving SoPlex\n");

                        if(result.scaled_before_rounding)
                        {
                           Rational number = result.rhs / obj_scale;
                           ceil_rational(number);
                           scaleInformation = { false, obj_scale, number };
                        }
                        data->vipr.write_node(status, assumptions,
                                              soplex_result.reduced, soplex_result.dual, result.rhs, result.nnz,
                                              data->problem, scaleInformation);
                        break;
                     }

                     auto diff = apply_Neumaier_Shcherbina(data, updated_domains, soplex_result.reduced, result.diff,
                                                           result.nnz);
                     if( diff.first )
                     {
                        result.rhs += diff.second;
                        Rational updated_obj_value = recalculate(data, soplex_result.reduced, soplex_result.dual,
                                                                 updated_domains);
                        auto comp_result = compare(data->problem.isObjectiveInteger(), data->scaled_obj_integer, obj_scale, data->best_bound, updated_obj_value);
                        if( comp_result == 1 )
                        {
                           if(data->problem.isObjectiveInteger())
                              ceil_rational(updated_obj_value);
                           msg.info("\tverified by applying Neumaier Shcherbina\n");
                           data->vipr.write_node(status, assumptions, soplex_result.reduced, soplex_result.dual,
                                                 updated_obj_value, result.nnz, data->problem );
                           data->stats.node_objlimit_neumaier_shcherbina( );
                           break;
                        }
                        if( comp_result == 2 )
                        {
                           Rational number = updated_obj_value/obj_scale;
                           ceil_rational(number);
                           msg.info("\tverified by applying Neumaier Shcherbina\n");
                           data->vipr.write_node(status, assumptions, soplex_result.reduced, soplex_result.dual,
                                                 updated_obj_value, result.nnz, data->problem,
                                                 { false, obj_scale, number });
                           data->stats.node_objlimit_neumaier_shcherbina( );
                           break;
                        }
                        msg.info("\trecalculating floating-point failed {} !>= {}\n",
                                 updated_obj_value,
                                 data->dual_bound.get().convert_to<double>( ));
                     }

                     reconstructVector(soplex_result.reduced, denom);
                     reconstructVector(soplex_result.dual, denom);

                     result = get_dualbound_and_check_for_obj_limit(data, soplex_result.reduced,
                                                                    soplex_result.dual, updated_domains, scip, 1);
                     if( result.success )
                     {
                        data->stats.node_deletion_reconstruct( );
                        msg.info("\treconstruct verified\n");
                        if(result.scaled_before_rounding)
                        {
                           Rational number = result.rhs / obj_scale;
                           ceil_rational(number);
                           scaleInformation = { false, obj_scale, number };
                        }
                        data->vipr.write_node(status, assumptions, soplex_result.reduced, soplex_result.dual,
                                              result.rhs,
                                              result.nnz, data->problem, scaleInformation);
                        break;
                     }
                  }
                  else if( soplex_result.status == NodeStatus::Infeasible )
                  {
                     NodeStatus status = parent ? NodeStatus::InfeasibleParent : NodeStatus::Infeasible;
                     msg.info("\tSoPlex: Infeasible instead of ObjLimit.\n");
                     ReturnValue result =
                           get_farkas_value_and_evaluate_farkas_proof(data,soplex_result.reduced, soplex_result.dual, updated_domains);
                     if( result.success )
                     {
                        msg.info("\tverified\n");
                        data->vipr.write_node(status, assumptions, soplex_result.reduced,
                                              soplex_result.dual, result.rhs, result.nnz, data->problem);
                        data->stats.node_objlimit_easy( );
                        break;
                     }
                     auto diff = apply_Neumaier_Shcherbina(data, updated_domains, soplex_result.reduced, result.diff, result.nnz);
                     result.rhs += diff.second;
                     if( diff.first && result.rhs > 0 )
                     {
                        msg.info("\tverified by applying Neumaier Shcherbina\n");
                        data->vipr.write_node(status, assumptions, soplex_result.reduced,
                                              soplex_result.dual, result.rhs, result.nnz, data->problem);
                        data->stats.node_objlimit_neumaier_shcherbina( );
                        break;
                     }

                     reconstructVector(soplex_result.reduced, denom);
                     reconstructVector(soplex_result.dual, denom);
                     result = get_farkas_value_and_evaluate_farkas_proof(data, soplex_result.reduced, soplex_result.dual, updated_domains);
                     if( result.success )
                     {
                        msg.info("\tverified after reconstruction.\n");
                        data->vipr.write_node(status, assumptions, soplex_result.reduced,
                                              soplex_result.dual, result.rhs, result.nnz, data->problem);
                        data->stats.node_objlimit_reconstruct( );
                        break;
                     }
                  }
                  else
                  {
                     assert(soplex_result.status == NodeStatus::LPError);
                     msg.info("\tFAILED: ERROR in floating-point SoPlex\n");
                     data->stats.floating_point_lp_error( );
                  }

                  soplex::SPxSolver::VarStatus rowstat[data->problem.getNRows()];
                  soplex::SPxSolver::VarStatus colstat[data->problem.getNCols()];
                  data->floating_soplex.get_basis(rowstat, colstat);

                  SoplexInterface exact_soplex{true};
                  exact_soplex.doSetUp(data->problem, updated_domains);

                  Vec<Rational> primal(soplex_result.reduced.size( ));
                  auto success = exact_soplex.factorize(primal, soplex_result.dual, soplex_result.reduced, rowstat, colstat, false);

                  if( success.first )
                  {
                     auto result = get_dualbound_and_check_for_obj_limit(data, soplex_result.reduced, soplex_result.dual,
                                                                    updated_domains, scip, 1);
                     if( result.success )
                     {
                        msg.info("\tfactorize is successful\n");
                        if(result.scaled_before_rounding)
                        {
                           Rational number = result.rhs / obj_scale;
                           ceil_rational(number);
                           scaleInformation = { false, obj_scale, number };
                        }
                        data->vipr.write_node(status, assumptions, soplex_result.reduced,
                                              soplex_result.dual, result.rhs, result.nnz, data->problem, scaleInformation);
                        data->stats.node_deletion_factorization( );
                        break;
                     }
                  }

                  exact_soplex.setBasis(rowstat, colstat);

                  auto soplex_exact_result = exact_soplex.solve(data->problem.getRowFlags( ), data->problem);

                  if( soplex_exact_result.status == NodeStatus::Solved)
                  {
                     auto result = get_dualbound_and_check_for_obj_limit(data, soplex_exact_result.reduced,
                                                                    soplex_exact_result.dual, updated_domains, scip, 1);

                     if( result.success )
                     {
                        data->stats.node_deletion_exact( );
                        msg.info("\tSoPlex-Exact: verified\n");
                        if(result.scaled_before_rounding)
                        {
                           Rational number = result.rhs / obj_scale;
                           ceil_rational(number);
                           scaleInformation = { false, obj_scale, number };
                        }
                        data->vipr.write_node(status, assumptions, soplex_exact_result.reduced,
                                              soplex_exact_result.dual,
                                              result.rhs, result.nnz, data->problem, scaleInformation);
                        break;
                     }
                     data->stats.node_deletion_error( );
                     msg.info("FAILED: Best bound {} exceeds bound {}\n", data->dual_bound.get(), result.rhs);
                     msg.info("\t\tBest bound {} exceeds bound {}\n", data->dual_bound.get().convert_to<double>(), result.rhs.convert_to<double>());
                     update_lower_bound(data, result.rhs);
                     data->vipr.write_node(status, assumptions, soplex_exact_result.reduced, soplex_exact_result.dual,
                                           result.rhs, result.nnz, data->problem);
                     break;
                  }
                  else if ( soplex_exact_result.status == NodeStatus::Infeasible)
                  {
                     status = parent ? NodeStatus::InfeasibleParent : NodeStatus::Infeasible;
                     ReturnValue result = get_farkas_value_and_evaluate_farkas_proof(data, soplex_exact_result.reduced,
                                                                                     soplex_exact_result.dual, updated_domains);
                     assert(soplex_exact_result.status == NodeStatus::Infeasible);
                     assert(result.success);
                     data->stats.node_deletion_exact( );
                     msg.info("\tSoPlex-Exact: Infeasible instead of ObjLimit.\n");
                     data->vipr.write_node(status, assumptions, soplex_exact_result.reduced, soplex_exact_result.dual,
                                           result.rhs, result.nnz, data->problem);
                  }
                  else
                  {
                     assert(soplex_exact_result.status == NodeStatus::LPError);
                     msg.info("FAILED: ERROR in Exact-SoPlex\n");
                     data->stats.exact_lp_error();
                     update_lower_bound_to_neg_inf(data);
                  }

               }
               break;
            }
#ifdef DEBUG
               case SCIP_EVENTTYPE_NODEFOCUSED:
               {
                  data->counter++;
                  auto node = SCIPeventGetNode(event);
                  msg.info("Focused nodes {} (SCIP number: {})\n" , data->counter, SCIPnodeGetNumber(node));
                  const Vec<Assumption> assumptions = get_assumptions(scip, SCIPeventGetNode(event), msg);
                  break;
               }
            case SCIP_EVENTTYPE_NODEBRANCHED:
            {
               auto node = SCIPeventGetNode(event);
               msg.info("Branched nodes {} (SCIP number: {})\n" , data->counter, SCIPnodeGetNumber(node));
               const Vec<Assumption> assumptions = get_assumptions(scip, SCIPeventGetNode(event), msg);
               break;
            }
#endif
            default:
               assert(false);
         }

         data->seconds_in_events += t.getTime();
         data->counter++;
         return SCIP_OKAY;
      }


      static
      void update_domains(VariableDomains<Rational> &domains, const Vec<Assumption> &assumptions, unsigned int size) {
         for( unsigned int i = 0; i < size; i++ )
         {
            Assumption a = assumptions[ i ];
            if( a.is_lower )
            {
               domains.lower_bounds[ a.var_index ] = a.value;
               domains.flags[ a.var_index ].unset(ColFlag::kLbInf);
            }
            else
            {
               domains.upper_bounds[ a.var_index ] = a.value;
               domains.flags[ a.var_index ].unset(ColFlag::kUbInf);
            }
         }
      }

      static std::pair<bool, Rational>
      apply_Neumaier_Shcherbina(const Exact_eventdata *data, const VariableDomains<Rational> &updated_domains,
                                Vec<Rational> &var_values, Vec<Rational> &diff, int &nnz) {
         Rational modification_rhs = 0;
         for( unsigned int i = 0; i < var_values.size( ); i++ )
         {
            if( diff[ i ] != 0 && var_values[ i ] == 0 )
            {
               if(( -diff[ i ] < 0 && updated_domains.flags[ i ].test(ColFlag::kUbInf)) ||
                  ( -diff[ i ] > 0 && updated_domains.flags[ i ].test(ColFlag::kLbInf)))
                  return { false, 0 };
               nnz++;
            }
            var_values[ i ] -= diff[ i ];
            if( diff[ i ] != 0 )
            {
               if(( -diff[ i ] < 0 && updated_domains.flags[ i ].test(ColFlag::kUbInf)) ||
                  ( -diff[ i ] > 0 && updated_domains.flags[ i ].test(ColFlag::kLbInf)))
                  return { false, 0 };
               if( var_values[ i ] < 0 )
                  modification_rhs -= ( diff[ i ] * updated_domains.upper_bounds[ i ] );
               else
                  modification_rhs -= ( diff[ i ] * updated_domains.lower_bounds[ i ] );
               if( var_values[ i ] == 0 )
                  nnz--;
            }
         }
         return { true, modification_rhs };
      }

   private:

      static void
      compute_pseudoobjective(Exact_eventdata *data, SCIP_NODE *node, SCIP *scip) {
         Vec<Assumption> assumptions = get_assumptions(scip, node, data->msg);
         VariableDomains<Rational> updated_domains = ( data->problem.getVariableDomains( ));
         update_domains(updated_domains, assumptions, assumptions.size( ));
         Rational bound = 0;
         int nnz = 0;
         Vec<Rational> zeros;
         std::fill_n(std::back_inserter(zeros), data->problem.getNRows( ), 0);
         for( int i = 0; i < data->problem.getNCols( ); i++ )
         {
            if( data->problem.getObjective( ).coefficients[ i ] == 0 )
               continue;
            nnz++;
            if( data->problem.getObjective( ).coefficients[ i ] > 0 )
            {
               bound += data->problem.getObjective( ).coefficients[ i ] *
                        updated_domains.lower_bounds[ i ];
            }
            else if( data->problem.getObjective( ).coefficients[ i ] < 0 )
               bound += data->problem.getObjective( ).coefficients[ i ] *
                        updated_domains.upper_bounds[ i ];
         }
         if(data->problem.isObjectiveInteger())
            ceil_rational(bound);
         assert(data->best_bound.get() <= bound);
         data->vipr.write_node(NodeStatus::ObjLimit, assumptions,
                               data->problem.getObjective( ).coefficients, zeros, bound, nnz, data->problem);
      }


      static std::pair<Vec<soplex::SPxSolver::VarStatus>, Vec<soplex::SPxSolver::VarStatus>>
      convert_basis(SCIP_VAR **vars, int nvars, SCIP_ROW **rows, int nrows, SCIP *scip) {
         int nlprows = SCIPgetNLPRows(scip);
         assert(nlprows <= nrows);
         int cstat[nvars];
         int rstat[nlprows];
         SCIP_LPI *lpi;
         auto retcode = SCIPgetLPI(scip, &lpi);
         assert(retcode == SCIP_OKAY);
         retcode = SCIPlpiGetBase(lpi, cstat, rstat);
         assert(retcode == SCIP_OKAY);

         retcode = SCIPlpiGetBase(lpi, cstat, rstat);
         assert(retcode == SCIP_OKAY);
         Vec<soplex::SPxSolver::VarStatus> soplex_colstat { };
         Vec<soplex::SPxSolver::VarStatus> soplex_rowstat { };
         soplex_colstat.resize(nvars, soplex::SPxSolver::BASIC);
         soplex_rowstat.resize(nrows, soplex::SPxSolver::BASIC);

         for( int var = 0; var < nvars; ++var )
         {
            int index = SCIPvarGetIndex(vars[ var ]);
            switch( cstat[ var ] ) /*lint !e613*/
            {
               case SCIP_BASESTAT_LOWER:
                  soplex_colstat[ index ] = soplex::SPxSolver::ON_LOWER;
                  break;
               case SCIP_BASESTAT_BASIC:
                  soplex_colstat[ index ] = soplex::SPxSolver::BASIC;
                  break;
               case SCIP_BASESTAT_UPPER:
                  soplex_colstat[ index ] = soplex::SPxSolver::ON_UPPER;
                  break;
               case SCIP_BASESTAT_ZERO:
                  soplex_colstat[ index ] = soplex::SPxSolver::ZERO;
                  break;
               default:
                  assert(false);
            }
         }

         for( int row = 0; row < nlprows; ++row )
         {
            int index = SCIProwGetIndex(rows[ row ]);
            switch( rstat[ row ] ) /*lint !e613*/
            {
               case SCIP_BASESTAT_LOWER:
                  soplex_rowstat[ index ] = soplex::SPxSolver::ON_LOWER;
                  break;
               case SCIP_BASESTAT_BASIC:
                  soplex_rowstat[ index ] = soplex::SPxSolver::BASIC;
                  break;
               case SCIP_BASESTAT_UPPER:
                  soplex_rowstat[ index ] = soplex::SPxSolver::ON_UPPER;
                  break;
               case SCIP_BASESTAT_ZERO:
                  assert(false);
               default:
                  assert(false);

            }
         }

         return { soplex_rowstat, soplex_colstat };
      }


      static Rational
      calculate_primal_obj(const Problem<Rational> &problem, const Vec<Rational> &primal) {
         Rational sum = 0;
         for( int i = 0; i < problem.getNCols( ); i++ )
            sum += primal[ i ] * problem.getObjective( ).coefficients[ i ];
         return sum;
      };

      static ReturnValue
      get_farkas_value_and_evaluate_farkas_proof(
            Exact_eventdata *data, const Vec<Rational> &var_farkas, const Vec<Rational> &dual_farkas,
            const VariableDomains<Rational> &domains) {

         const Vec<Rational> &rhs = data->problem.getConstraintMatrix( ).getRightHandSides( );
         const Vec<Rational> &lhs = data->problem.getConstraintMatrix( ).getLeftHandSides( );

         int nnz_farkas_coefficients = 0;
         Vec<Rational> sum;
         sum.resize(data->problem.getNCols( ), 0);
         Rational sum_rhs = 0;

         for( int i = 0; i < data->problem.getNCols( ); i++ )
         {
            sum[ i ] = var_farkas[ i ];
            if( var_farkas[ i ] < 0 )
            {
               assert(!domains.flags[ i ].test(ColFlag::kUbInf));
               sum_rhs += var_farkas[ i ] * domains.upper_bounds[ i ];
               nnz_farkas_coefficients++;
            }
            else if( var_farkas[ i ] > 0 )
            {
               assert(!domains.flags[ i ].test(ColFlag::kLbInf));
               sum_rhs += var_farkas[ i ] * domains.lower_bounds[ i ];
               nnz_farkas_coefficients++;
            }
         }
         for( int i = 0; i < data->problem.getNRows( ); i++ )
         {
            if( dual_farkas[ i ] > 0 )
            {
               assert(!data->problem.getRowFlags( )[ i ].test(RowFlag::kLhsInf));
               sum_rhs += dual_farkas[ i ] * lhs[ i ];
               nnz_farkas_coefficients++;
            }
            else if( dual_farkas[ i ] < 0 )
            {
               assert(!data->problem.getRowFlags( )[ i ].test(RowFlag::kRhsInf));
               sum_rhs += dual_farkas[ i ] * rhs[ i ];
               nnz_farkas_coefficients++;
            }
            if( dual_farkas[ i ] != 0 )
            {
               auto row_data = data->problem.getConstraintMatrix( ).getRowCoefficients(i);
               for( int j = 0; j < row_data.getLength( ); j++ )
                  sum[ row_data.getIndices( )[ j ]] += dual_farkas[ i ] * row_data.getValues( )[ j ];
            }
         }
         if( sum_rhs <= 0 )
            return ReturnValue { false, sum_rhs, sum, nnz_farkas_coefficients };
         for( int i = 0; i < data->problem.getNCols( ); i++ )
            if( sum[ i ] != 0 )
               return ReturnValue { false, sum_rhs, sum, nnz_farkas_coefficients };
         return ReturnValue { true, sum_rhs, sum, nnz_farkas_coefficients };
      }

      static bool
      check_integrality(Exact_eventdata *data, Vec<Rational> &primal) {
         const Problem<Rational> &problem = data->problem;
         for( int i = 0; i < problem.getNCols( ); i++ )
            if( problem.getColFlags( )[ i ].test(ColFlag::kIntegral) && primal[ i ].convert_to<int>( ) != primal[ i ] )
               return false;
         return true;
      }

      static bool
      check_feasibility(Exact_eventdata *data, Vec<Rational> &primal, bool round) {
         const Problem<Rational> &problem = data->problem;
         for( int i = 0; i < problem.getNCols( ); i++ )
         {
            if( problem.getColFlags( )[ i ].test(ColFlag::kIntegral))
            {
               if( primal[ i ].convert_to<int>( ) != primal[ i ] )
               {
                  if( round )
                     primal[ i ] = Rational(( primal[ i ] + 0.5 ).convert_to<int>( ));
                  else
                     return false;
               }
            }
            if( !problem.getColFlags( )[ i ].test(ColFlag::kUbInf) && problem.getUpperBounds( )[ i ] < primal[ i ] )
            {
               if( !round )
                  return false;
               primal[i] = problem.getUpperBounds( )[ i ];
            }
            if( !problem.getColFlags( )[ i ].test(ColFlag::kLbInf) && problem.getLowerBounds( )[ i ] > primal[ i ] )
            {
               if( !round )
                  return false;
               primal[i] = problem.getLowerBounds( )[ i ];
            }
         }
         for( int row = 0; row < problem.getNRows( ); row++ )
         {
            Rational sum = 0;
            auto row_vector = problem.getConstraintMatrix( ).getRowCoefficients(row);
            for( int j = 0; j < row_vector.getLength( ); ++j )
               sum += primal[ row_vector.getIndices( )[ j ]] * row_vector.getValues( )[ j ];
            bool violated = ( !problem.getRowFlags( )[ row ].test(RowFlag::kLhsInf) &&
                              sum < problem.getConstraintMatrix( ).getLeftHandSides( )[ row ] )
                            || ( !problem.getRowFlags( )[ row ].test(RowFlag::kRhsInf) &&
                                 sum > problem.getConstraintMatrix( ).getRightHandSides( )[ row ] );
           if( violated )
               return false;
         }
         return true;
      }


      static int
      count_nonzeros(const Vec<Rational> &vector) {
         int nnz = 0;
         for(const auto & index : vector)
            if( index != 0 )
               nnz++;
         return nnz;
      }

      static ReturnValue
      compare_dual_bound_to_primal(
            Exact_eventdata *data, const Vec<Rational> &redcosts, const Vec<Rational> &dual_solution,
            const VariableDomains<Rational> &updated_domains, SCIP *scip, Rational &dual_bound) {
         const Vec<Rational> &upper_bounds = updated_domains.upper_bounds;
         const Vec<Rational> &lower_bounds = updated_domains.lower_bounds;
         const Vec<Rational> &rhs = data->problem.getConstraintMatrix( ).getRightHandSides( );
         const Vec<Rational> &lhs = data->problem.getConstraintMatrix( ).getLeftHandSides( );


         int nnz = 0;
         Vec<Rational> sum;
         std::fill_n(std::back_inserter(sum), data->problem.getNCols( ), 0);
         Rational sum_rhs = 0;

         for( int index = 0; index < data->problem.getNCols( ); index++ )
         {
            sum[ index ] = redcosts[ index ];
            if( redcosts[ index ] < 0 )
            {
               assert(!updated_domains.flags[ index ].test(ColFlag::kUbInf));
               sum_rhs += redcosts[ index ] * upper_bounds[ index ];
               nnz++;
            }
            else if( redcosts[ index ] > 0 )
            {
               assert(!updated_domains.flags[ index ].test(ColFlag::kLbInf));
               sum_rhs += redcosts[ index ] * lower_bounds[ index ];
               nnz++;
            }
         }

         for( int index = 0; index < data->problem.getNRows( ); index++ )
         {
            if( dual_solution[ index ] > 0 )
            {
               assert(!data->problem.getRowFlags( )[ index ].test(RowFlag::kLhsInf));
               sum_rhs += dual_solution[ index ] * lhs[ index ];
               nnz++;
            }
            else if( dual_solution[ index ] < 0 )
            {
               assert(!data->problem.getRowFlags( )[ index ].test(RowFlag::kRhsInf));
               sum_rhs += dual_solution[ index ] * rhs[ index ];
               nnz++;
            }
            if( dual_solution[ index ] != 0 )
            {
               auto row_data = data->problem.getConstraintMatrix( ).getRowCoefficients(index);
               for( int j = 0; j < row_data.getLength( ); j++ )
                  sum[ row_data.getIndices( )[ j ]] += dual_solution[ index ] * row_data.getValues( )[ j ];
            }
         }

         for( int i = 0; i < data->problem.getNCols( ); i++ )
            sum[ i ] -= data->problem.getObjective( ).coefficients[ i ];

         Rational diff = sum_rhs - dual_bound;
         if( sum_rhs != dual_bound )
            return ReturnValue { false, sum_rhs, sum, nnz };
         for( int i = 0; i < data->problem.getNCols( ); i++ )
            if( sum[ i ] != 0 )
               return ReturnValue { false, sum_rhs, sum, nnz };
         return ReturnValue { true, sum_rhs, sum, nnz };
      }


      static Rational
      compute_dual(Exact_eventdata *data, const Vec<Rational> &redcosts, const Vec<Rational> &dual_solution,
                   const VariableDomains<Rational> &updated_domains) {
         Rational result = 0;
         const Vec<Rational> &rhs = data->problem.getConstraintMatrix( ).getRightHandSides( );
         const Vec<Rational> &lhs = data->problem.getConstraintMatrix( ).getLeftHandSides( );
         for( int index = 0; index < data->problem.getNCols( ); index++ )
         {
            if( redcosts[ index ] < 0 )
               result += redcosts[ index ] * updated_domains.upper_bounds[ index ];
            else if( redcosts[ index ] > 0 )
               result += redcosts[ index ] * updated_domains.lower_bounds[ index ];
         }

         for( int row = 0; row < data->problem.getNRows( ); row++ )
         {
            if( dual_solution[ row ] > 0 )
               result += dual_solution[ row ] * lhs[ row ];
            else if( dual_solution[ row ] < 0 )
               result += dual_solution[ row ] * rhs[ row ];
         }

         if( data->problem.isObjectiveInteger( ) )
            ceil_rational(result);

         return result;
      }

      //if scaling is 1 the objective does not have to be scaled anymore
      static ReturnValue
      get_dualbound_and_check_for_obj_limit(
            Exact_eventdata *data, const Vec<Rational> &redcosts, const Vec<Rational> &dual_solution,
            const VariableDomains<Rational> &updated_domains, SCIP *scip, Rational scaling) {
         assert(scaling != 0);
         const Vec<Rational> &rhs = data->problem.getConstraintMatrix( ).getRightHandSides( );
         const Vec<Rational> &lhs = data->problem.getConstraintMatrix( ).getLeftHandSides( );

         Rational scale = Rational(SCIPgetTransObjscale(scip));
         assert(scaling == scale || scaling == 1);
         int nnz = 0;
         Vec<Rational> sum;
         sum.resize(data->problem.getNCols( ), 0);
         Rational sum_rhs = 0;

         for( int index = 0; index < data->problem.getNCols( ); index++ )
         {
            sum[ index ] = redcosts[ index ];
            if( redcosts[ index ] < 0 )
            {
               sum_rhs += redcosts[ index ] * updated_domains.upper_bounds[ index ];
               nnz++;
            }
            else if( redcosts[ index ] > 0 )
            {
               sum_rhs += redcosts[ index ] * updated_domains.lower_bounds[ index ];
               nnz++;
            }
         }

         for( int row = 0; row < data->problem.getNRows( ); row++ )
         {
            if( dual_solution[ row ] > 0 )
            {
               assert(!data->problem.getRowFlags( )[ row ].test(RowFlag::kLhsInf));
               sum_rhs += dual_solution[ row ] * lhs[ row ];
               nnz++;
            }
            else if( dual_solution[ row ] < 0 )
            {
               assert(!data->problem.getRowFlags( )[ row ].test(RowFlag::kRhsInf));
               sum_rhs += dual_solution[ row ] * rhs[ row ];
               nnz++;
            }
            if( dual_solution[ row ] != 0 )
            {
               auto row_data = data->problem.getConstraintMatrix( ).getRowCoefficients(row);
               for( int j = 0; j < row_data.getLength( ); j++ )
                  sum[ row_data.getIndices( )[ j ]] += dual_solution[ row ] * row_data.getValues( )[ j ];
            }
         }

         for( int i = 0; i < data->problem.getNCols( ); i++ )
            sum[ i ] -= data->problem.getObjective( ).coefficients[ i ] / scaling;

         sum_rhs = sum_rhs * scaling;
         auto result = compare(data->problem.isObjectiveInteger(), data->scaled_obj_integer, scale, data->best_bound, sum_rhs);
         // fail
         if( result == 0 )
            return ReturnValue { false, sum_rhs, sum, nnz };
         bool rounded_before_scaling = result == 2;
         for( int i = 0; i < data->problem.getNCols( ); i++ )
            if( sum[ i ] != 0 )
               return ReturnValue { false, sum_rhs, sum, nnz, rounded_before_scaling };
         if( !rounded_before_scaling && data->problem.isObjectiveInteger())
            ceil_rational(sum_rhs);
         return ReturnValue { true, sum_rhs, sum, nnz, rounded_before_scaling };
      }

      static int
      compare(bool integer_obj, bool scaled_integer_obj, Rational obj_scale, Bound &dual_bound,
              Rational &current_scaled_obj) {
         bool is_obj_scaled = obj_scale != 1;
         if( !scaled_integer_obj )
         {
            assert(!integer_obj);
            if( current_scaled_obj >= dual_bound.get() )
               return  1;
            return  0;
         }
         if( !is_obj_scaled && integer_obj )
         {
            ceil_rational(current_scaled_obj);
            if( current_scaled_obj >= dual_bound.get() )
               return 1;
            return 0;
         }
         if( scaled_integer_obj && !integer_obj )
         {
            assert(is_obj_scaled);
            Rational result = current_scaled_obj / obj_scale;
            ceil_rational(result);
            result = result * obj_scale;
            if( result >= dual_bound.get() )
               return 2;
            return 0;
         }
         assert(scaled_integer_obj && integer_obj);
         assert(is_obj_scaled);
         Rational result = current_scaled_obj;
         ceil_rational(result);
         if( result >= dual_bound.get() )
            return 1;
         result = current_scaled_obj / obj_scale;
         ceil_rational(result);
         result = result * obj_scale;
         if( result >= dual_bound.get() )
            return 2 ;
         return 0;
      }


      static void update_lower_bound(Exact_eventdata *data, Rational val)
      {
         if( !data->dual_bound.is_initialized() || val < data->dual_bound.get())
            data->dual_bound.update(val);
      }

      static void update_lower_bound_to_neg_inf(Exact_eventdata *data)
      {
         data->dual_bound.update_neg_infinity();
      }


      static void add_new_solution(Exact_eventdata *data, const Vec<Rational> &primal, Rational val)
      {
         if( !data->best_bound.is_initialized() || val < data->best_bound.get() )
            data->best_bound.update( val );
         data->solutions.emplace_back(primal, val);
         update_lower_bound(data, val);
      }

      static Rational
      recalculate(Exact_eventdata *data, const Vec<Rational> &redcosts, const Vec<Rational> &dual_solution,
                  const VariableDomains<Rational> &domains) {

         const Vec<Rational> &rhs = data->problem.getConstraintMatrix( ).getRightHandSides( );
         const Vec<Rational> &lhs = data->problem.getConstraintMatrix( ).getLeftHandSides( );
         const Vec<Rational> &ub = domains.upper_bounds;
         const Vec<Rational> &lb = domains.lower_bounds;

         Rational result = 0;

         for( int i = 0; i < data->problem.getNCols( ); i++ )
            if( redcosts[ i ] < 0 )
               result += redcosts[ i ] * ub[ i ];
            else if( redcosts[ i ] > 0 )
               result += redcosts[ i ] * lb[ i ];

         for( int i = 0; i < data->problem.getNRows( ); i++ )
            if( dual_solution[ i ] > 0 )
               result += dual_solution[ i ] * lhs[ i ];
            else if( dual_solution[ i ] < 0 )
               result += dual_solution[ i ] * rhs[ i ];

         return result;
      }

      static void ceil_rational(Rational &value) {

         value = ceil(value);
         //Alternative from VIPR
         // mpz_t z;
         // auto num = numerator(value);
         // auto den = denominator(value);
         // mpz_init (z);
         // mpz_fdiv_q(z, num.backend().data(), den.backend().data() ); // Divide numerator by denominator and floor the result
         // Rational value3 = Rational(z + 1);
         // mpz_clear (z);


         // Alternative not working in Boost <1.79
         // if( value != Rational(value.convert_to<int>( )))
         //    value = value < 0 ? Rational(value.convert_to<int>( )) : Rational(value.convert_to<int>( ) + 1);
      }



      /** compare SCIPprintNodeRootPath*/
      static Vec<Assumption> get_assumptions(SCIP *scip, SCIP_NODE *node, Message& msg) {
         SCIP_VAR **branchvars;
         SCIP_Real *branchbounds;
         SCIP_BOUNDTYPE *boundtypes;
         int *nodeswitches;
         int nbranchvars;
         int nnodes;
         int branchvarssize = SCIPnodeGetDepth(node);
         int nodeswitchsize = branchvarssize;

         /* memory allocation */
         SCIPallocBufferArray(scip, &branchvars, branchvarssize);
         SCIPallocBufferArray(scip, &branchbounds, branchvarssize);
         SCIPallocBufferArray(scip, &boundtypes, branchvarssize);
         SCIPallocBufferArray(scip, &nodeswitches, nodeswitchsize);

         SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize,
                                          nodeswitches, &nnodes, nodeswitchsize);

         /* if the arrays were too small, we have to reallocate them and recall SCIPnodeGetAncestorBranchingPath */
         if( nbranchvars > branchvarssize || nnodes > nodeswitchsize )
         {
            branchvarssize = nbranchvars;
            nodeswitchsize = nnodes;

            /* memory reallocation */
            SCIPreallocBufferArray(scip, &branchvars, branchvarssize);
            SCIPreallocBufferArray(scip, &branchbounds, branchvarssize);
            SCIPreallocBufferArray(scip, &boundtypes, branchvarssize);
            SCIPreallocBufferArray(scip, &nodeswitches, nodeswitchsize);
            SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize,
                                             nodeswitches, &nnodes, nodeswitchsize);
            assert(nbranchvars == branchvarssize);
         }


         msg.info("\t");
         Vec<Assumption> assumptions;
         /* we only want to create output, if branching were performed */
         if( nbranchvars >= 1 )
         {

            /* print all nodes, starting from the root, which is last in the arrays */
            for( int j = nnodes - 1; j >= 0; --j )
            {
               int end;
               if( j == nnodes - 1 )
                  end = nbranchvars;
               else
                  end = nodeswitches[ j + 1 ];


               for( int i = nodeswitches[ j ]; i < end; ++i )
               {
                  SCIP_Real scalar = 0;
                  SCIP_Real constant = 0;
                  SCIPvarGetOrigvarSum(&branchvars[ i ], &scalar, &constant);
                  assumptions.emplace_back(SCIPvarGetIndex(( branchvars[ i ] )), branchbounds[ i ],
                                           boundtypes[ i ] == SCIP_BOUNDTYPE_LOWER);

                  msg.info("{} {} {} -", SCIPvarGetName(branchvars[ i ]),
                             boundtypes[ i ] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", branchbounds[ i ]);


               }
            }

         }
         msg.info("\n");
         /* free all local memory */
         SCIPfreeBufferArray(scip, &nodeswitches);
         SCIPfreeBufferArray(scip, &boundtypes);
         SCIPfreeBufferArray(scip, &branchbounds);
         SCIPfreeBufferArray(scip, &branchvars);
         return assumptions;
      }


      using Integer = boost::multiprecision::number<boost::multiprecision::gmp_int, boost::multiprecision::et_off>;

      /** reconstruct a rational vector */
      static bool reconstructVector(Vec<Rational> &input, const Rational &denomBoundSquared) {
         Vec<Integer> xnum(input.size( )); /* numerator of input vector */
         Integer denom = 1; /* common denominator of input vector */
         int dim = input.size( );

         /* find common denominator */
//         if(indexSet == nullptr)
//         {
         for( int i = 0; i < dim; i++ )
            lcm(denom, denom, denominator(input[ i ]));

         for( int i = 0; i < dim; i++ )
         {
            xnum[ i ] = denom * Integer(numerator(input[ i ]));
            xnum[ i ] = xnum[ i ] / Integer(denominator(input[ i ]));
         }

         /* reconstruct */
         return Reconstruct(input, xnum.data( ), denom, dim, denomBoundSquared);

      }

      /** this reconstruction routine will set x equal to the mpq vector where each component is the best rational
       *  approximation of xnum / denom with where the GCD of denominators of x is at most Dbound; it will return true on
       *  success and false if more accuracy is required: specifically if componentwise rational reconstruction does not
       *  produce such a vector
       */
      static int Reconstruct(Vec<Rational> &resvec, Integer *xnum, Integer denom, int dim,
                             const Rational &denomBoundSquared) {
         bool rval = true;
         int done = 0;

         /* denominator must be positive */
         assert(denom > 0);
         assert(denomBoundSquared > 0);

         Integer temp = 0;
         Integer td = 0;
         Integer tn = 0;
         Integer Dbound = 0;
         Integer gcd = 1;

         Dbound = numerator(denomBoundSquared) /
                  boost::multiprecision::numerator(denomBoundSquared);

         Dbound = ( Integer ) sqrt(Dbound);

         /* if Dbound is below 2^24 increase it to this value, this avoids changing input vectors that have low denominator
          * because they are floating point representable
          */
         if( Dbound < 16777216 )
            Dbound = 16777216;

         /* The following represent a_i, the cont frac representation and p_i/q_i, the convergents */
         Integer a0 = 0;
         Integer ai = 0;

         /* here we use p[2]=pk, p[1]=pk-1,p[0]=pk-2 and same for q */
         Integer p[3];
         Integer q[3];

//         for( int c = 0; ( indexSet == nullptr && c < dim ) || ( indexSet != nullptr && c < indexSet->size( )); c++ )
         for( int c = 0; c < dim; c++ )
         {
            int j = c;

            assert(j >= 0);
            assert(j < dim);

            /* if xnum =0 , then just leave x[j] as zero */
            if( xnum[ j ] != 0 )
            {
               /* setup n and d for computing a_i the cont. frac. rep */
               tn = xnum[ j ];
               td = denom;

               /* divide tn and td by gcd */
               mpz_gcd(temp.backend( ).data( ), tn.backend( ).data( ), td.backend( ).data( ));
               tn = tn / temp;
               td = td / temp;

               if( td <= Dbound )
                  resvec[ j ] = Rational(tn, td);

               else
               {

                  temp = 1;

                  divide_qr(tn, td, a0, temp);

                  tn = td;
                  td = temp;

                  divide_qr(tn, td, ai, temp);
                  tn = td;
                  td = temp;

                  p[ 1 ] = a0;
                  p[ 2 ] = 1;
                  p[ 2 ] += a0 * ai;

                  q[ 1 ] = 1;
                  q[ 2 ] = ai;

                  done = 0;

                  /* if q is already big, skip loop */
                  if( q[ 2 ] > Dbound )
                     done = 1;


                  while( !done && td != 0 )
                  {
                     /* update everything: compute next ai, then update convergents */

                     /* update ai */
                     divide_qr(tn, td, ai, temp);
                     tn = td;
                     td = temp;

                     /* shift p,q */
                     q[ 0 ] = q[ 1 ];
                     q[ 1 ] = q[ 2 ];
                     p[ 0 ] = p[ 1 ];
                     p[ 1 ] = p[ 2 ];

                     /* compute next p,q */
                     p[ 2 ] = p[ 0 ];
                     p[ 2 ] += p[ 1 ] * ai;
                     q[ 2 ] = q[ 0 ];
                     q[ 2 ] += q[ 1 ] * ai;

                     if( q[ 2 ] > Dbound )
                        done = 1;

                  }

                  assert(q[ 1 ] != 0);

                  /* Assign the values */
                  if( q[ 1 ] >= 0 )
                     resvec[ j ] = Rational(p[ 1 ], q[ 1 ]);
                  else
                     resvec[ j ] = Rational(-p[ 1 ], -q[ 1 ]);

                  calc_gcd(temp, gcd, denominator(resvec[ j ]));

                  gcd *= temp;

                  if( gcd > Dbound )
                  {
                     rval = false;
                     break;
                  }
               }
            }
         }

         return rval;
      }

      static void lcm(Integer &result, Integer a, Integer b) {
         mpz_lcm(result.backend( ).data( ), a.backend( ).data( ), b.backend( ).data( ));
      }

      static void calc_gcd(Integer &result, Integer a, Integer b) {
         mpz_gcd(result.backend( ).data( ), a.backend( ).data( ), b.backend( ).data( ));
      }


   };

} // namespace exact

#endif
