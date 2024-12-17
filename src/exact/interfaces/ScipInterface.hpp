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

#ifndef EXACT_INTERFACES_SCIP_INTERFACE_HPP_
#define EXACT_INTERFACES_SCIP_INTERFACE_HPP_

#define UNUSED(expr) do { (void)(expr); } while (0)

#include "exact/misc/Vec.hpp"
#include <cassert>
#include <stdexcept>

#include "exact/data/Problem.hpp"
#include "scip/cons_linear.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_paramset.h"
#include "Leaf.hpp"

namespace exact {


   class ScipInterface {
   private:
      SCIP *scip;
      Vec<SCIP_VAR *> variables;


   public:
      ScipInterface( ) : scip(nullptr) {
         if( SCIPcreate(&scip) != SCIP_OKAY )
            throw std::runtime_error("could not create SCIP");
      }

      const Vec<SCIP_Var *>
      getVariables( ) {
         return variables;
      }

      SCIP *
      getSCIP( ) {
         return scip;
      }

      SCIP_RETCODE
      doSetUp(const Problem<exact::Rational> &problem, const VariableDomains<exact::Rational>& domains) {
         auto return_code = SCIPcreate(&scip);
         SCIP_CALL_ABORT(SCIPincludeDefaultPlugins(scip));
         assert(return_code == SCIP_OKAY);
         int ncols = problem.getNCols( );
         int nrows = problem.getNRows( );
         const Vec<String> &varNames = problem.getVariableNames( );
         const Vec<String> &consNames = problem.getConstraintNames( );
         const Objective<exact::Rational> &obj = problem.getObjective( );
         const auto &consMatrix = problem.getConstraintMatrix( );
         const auto &lhs_values = consMatrix.getLeftHandSides( );
         const auto &rhs_values = consMatrix.getRightHandSides( );
         const auto &rflags = problem.getRowFlags( );

         SCIP_CALL(SCIPcreateProbBasic(scip, problem.getName( ).c_str( )));

         variables.resize(problem.getNCols( ));

         for( int i = 0; i < ncols; ++i )
         {
            SCIP_VAR *var;
            assert(!domains.flags[ i ].test(ColFlag::kInactive));

            SCIP_Real lb = domains.flags[ i ].test(ColFlag::kLbInf)
                           ? -SCIPinfinity(scip)
                           : domains.lower_bounds[ i ].convert_to<double>();
            SCIP_Real ub = domains.flags[ i ].test(ColFlag::kUbInf)
                           ? SCIPinfinity(scip)
                           : domains.upper_bounds[ i ].convert_to<double>();
            SCIP_VARTYPE type;
            if( domains.flags[ i ].test(ColFlag::kIntegral))
            {
               if( lb == exact::Rational { 0 } && ub == exact::Rational { 1 } )
                  type = SCIP_VARTYPE_BINARY;
               else
                  type = SCIP_VARTYPE_INTEGER;
            }
            else if( domains.flags[ i ].test(ColFlag::kImplInt))
               type = SCIP_VARTYPE_IMPLINT;
            else
               type = SCIP_VARTYPE_CONTINUOUS;

            SCIP_CALL(SCIPcreateVarBasic(
                  scip, &var, varNames[ i ].c_str( ), lb, ub,
                  SCIP_Real(obj.coefficients[ i ]), type));
            SCIP_CALL(SCIPaddVar(scip, var));
            variables[ i ] = var;

            SCIP_CALL(SCIPreleaseVar(scip, &var));
         }

         Vec<SCIP_VAR *> consvars;
         Vec<SCIP_Real> consvals;
         consvars.resize(problem.getNCols( ));
         consvals.resize(problem.getNCols( ));

         for( int i = 0; i < nrows; ++i )
         {
            SCIP_CONS *cons;

            auto rowvec = consMatrix.getRowCoefficients(i);
            const exact::Rational *vals = rowvec.getValues( );
            const int *inds = rowvec.getIndices( );
            SCIP_Real lhs = rflags[ i ].test(RowFlag::kLhsInf)
                            ? -SCIPinfinity(scip)
                            : SCIP_Real(lhs_values[ i ]);
            SCIP_Real rhs = rflags[ i ].test(RowFlag::kRhsInf)
                            ? SCIPinfinity(scip)
                            : SCIP_Real(rhs_values[ i ]);

            for( int k = 0; k < rowvec.getLength( ); ++k )
            {
               consvars[ k ] = variables[ inds[ k ]];
               consvals[ k ] = SCIP_Real(vals[ k ]);
            }

            SCIP_CALL(SCIPcreateConsBasicLinear(
                  scip, &cons, consNames[ i ].c_str( ), rowvec.getLength( ),
                  consvars.data( ), consvals.data( ), lhs, rhs));
            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));

         }

         if( obj.offset != exact::Rational { 0 } )
            SCIP_CALL(SCIPaddOrigObjoffset(scip, SCIP_Real(obj.offset)));

         return SCIP_OKAY;
      }


      SCIP_Status
      solve( ) {
         SCIP_RETCODE returncode = SCIPsolve(scip);
         assert(returncode == SCIP_OKAY);
         return SCIPgetStatus(scip);
      };

      void
      set_time_limit(int time_in_seconds)
      {
         SCIP_RETCODE status = SCIPsetRealParam(scip, "limits/time", time_in_seconds);
         assert(status == SCIP_OKAY);
      }

      void
      setseed(int seed) {
         auto status = SCIPsetIntParam(scip, "randomization/permutationseed", seed);
//         status = SCIPsetBoolParam(scip, "randomization/permutevars", TRUE);

      }
      void
      enable_branch_and_bound_only( double feastol, int seed ) {

         SCIP_RETCODE status = SCIPsetBoolParam(scip, "conflict/enable", false);
         assert(status == SCIP_OKAY);

         status = SCIPsetRealParam(scip, "numerics/feastol", feastol);
         assert(status == SCIP_OKAY);
//         status = SCIPsetIntParam(scip, "randomization/permutationseed", seed);
//         assert(status == SCIP_OKAY);
         status = SCIPsetBoolParam(scip, "branching/fullstrong/probingbounds", false);
         assert(status == SCIP_OKAY);
         status = SCIPsetBoolParam(scip, "branching/multaggr/probingbounds", false);
         assert(status == SCIP_OKAY);
         status = SCIPsetBoolParam(scip, "branching/relpscost/probingbounds", false);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "branching/fullstrong/maxproprounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "branching/relpscost/maxproprounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "branching/relpscost/maxbdchgs", 0);
         assert(status == SCIP_OKAY);



         status = SCIPsetIntParam(scip, "presolving/maxrounds", 0);
         assert(status == SCIP_OKAY);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/maxrestarts", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/maxroundsroot", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/nonlinear/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/nonlinear/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/linear/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/linear/maxprerounds", 0);
         assert(status == SCIP_OKAY);
//         status = SCIPsetIntParam(scip, "constraints/linear-exact/sepafreq", -1);
//         assert(status == SCIP_OKAY);
//         status = SCIPsetIntParam(scip, "constraints/linear-exact/maxprerounds", 0);
//         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/and/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/and/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/bounddisjunction/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/cardinality/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/cardinality/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/conjunction/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/cumulative/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/cumulative/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/disjunction/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/indicator/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/indicator/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/knapsack/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/knapsack/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/linking/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/linking/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/logicor/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/logicor/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/or/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/or/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/orbisack/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/orbisack/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/orbitope/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/pseudoboolean/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/setppc/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/setppc/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/SOS1/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/SOS1/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/SOS2/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/SOS2/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/superindicator/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/symresack/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/symresack/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/varbound/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/varbound/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/xor/sepafreq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/xor/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "constraints/components/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/domcol/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/dualcomp/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/gateextraction/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/implics/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/inttobinary/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/milp/maxrounds", 0);
//         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/trivial/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/sparsify/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "presolving/dualsparsify/maxrounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/adaptivediving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/clique/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/completesol/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/conflictdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/crossover/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/distributiondiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/farkasdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/feaspump/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/fracdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/gins/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/guideddiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/indicator/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/indicatordiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/intshifting/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/locks/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/lpface/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/alns/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/nlpdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/multistart/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/mpec/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/ofins/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/oneopt/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/padm/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/pscostdiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/randrounding/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/rens/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/reoptsols/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/rins/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/rootsoldiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/rounding/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/shifting/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/simplerounding/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/subnlp/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/trivial/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/trivialnegation/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/trysol/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/undercover/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/vbounds/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/veclendiving/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "heuristics/zirounding/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/dualfix/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/genvbounds/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/obbt/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/nlobbt/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/pseudoobj/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/redcost/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/rootredcost/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/symmetry/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "propagating/vbounds/maxprerounds", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "misc/usesymmetry", 0);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/clique/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/flowcover/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/cmir/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/knapsackcover/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/aggregation/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/disjunctive/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/gomory/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/strongcg/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/gomorymi/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/impliedbounds/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/mcf/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/minor/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/mixing/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/rapidlearning/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/rlt/freq", -1);
         assert(status == SCIP_OKAY);
         status = SCIPsetIntParam(scip, "separating/zerohalf/freq", -1);
         assert(status == SCIP_OKAY);

      }

      ~ScipInterface( ) {
         if( scip != nullptr )
         {

            SCIP_RETCODE retcode = SCIPfree(&scip);
            UNUSED(retcode);
            assert(retcode == SCIP_OKAY);
         }
      }
   };

} // namespace exact

#endif
