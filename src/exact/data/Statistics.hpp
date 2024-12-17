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

#ifndef _EXACT_CORE_STATS_HPP_
#define _EXACT_CORE_STATS_HPP_

namespace exact {

   class Statistics {


      int best_solution[4] = { };
      int node_feasible[5] = { };
      int node_deletion_pseudo[2] = { };
      int node_deletion[6] = { };
      int node_infeasible[5] = { };
      int node_objlimit[6] = { };
      int node_objlimit_not_solved[2] = { };
      int node_objlimit_pseudo[4] = { };
      int lperror[2] = { };

   public:

      Statistics( ) { }

      void solution_feasible( ) {
         best_solution[ 0 ]++;
      }

      void solution_reconstruction_feasible( ) {
         best_solution[ 1 ]++;
      }

      void solution_integer_fixing_feasible( ) {
         best_solution[ 2 ]++;
      }

      void solution_error( ) {
         best_solution[ 3 ]++;
      }

      void node_feasible_easy( ) {
         node_feasible[ 0 ]++;
      }

      void node_feasible_reconstruction( ) {
         node_feasible[ 1 ]++;
      }

      void node_feasible_factorization( ) {
         node_feasible[ 2 ]++;
      }

      void node_feasible_exact_solve( ) {
         node_feasible[ 3 ]++;
      }

      void node_feasible_error( ) {
         node_feasible[ 4 ]++;
      }

      void node_deletion_pseudo_easy( ) {
         node_deletion_pseudo[ 0 ]++;
      }

      void node_deletion_pseudo_error( ) {
         node_deletion_pseudo[ 1 ]++;
      }

      void node_deletion_soplex( ) {
         node_deletion[ 0 ]++;
      }

      void node_deletion_neumaier_shcherbina( ) {
         node_deletion[ 1 ]++;
      }

      void node_deletion_reconstruct( ) {
         node_deletion[ 2 ]++;
      }

      void node_deletion_factorization( ) {
         node_deletion[ 3 ]++;
      }

      void node_deletion_exact( ) {
         node_deletion[ 4 ]++;
      }

      void node_deletion_error( ) {
         node_deletion[ 5 ]++;
      }

      void node_infeasible_easy( ) {
         node_infeasible[ 0 ]++;
      }

      void node_infeasible_neumaier_shcherbina( ) {
         node_infeasible[ 1 ]++;
      }

      void node_infeasible_reconstruct( ) {
         node_infeasible[ 2 ]++;
      }

      void node_infeasible_exact( ) {
         node_infeasible[ 3 ]++;
      }

      void node_infeasible_error( ) {
         node_infeasible[ 4 ]++;
      }

      void node_objlimit_easy( ) {
         node_objlimit[ 0 ]++;
      }

      void node_objlimit_neumaier_shcherbina( ) {
         node_objlimit[ 1 ]++;
      }

      void node_objlimit_reconstruct( ) {
         node_objlimit[ 2 ]++;
      }

      void node_objlimit_factorize( ) {
         node_objlimit[ 3 ]++;
      }

      void node_objlimit_exact( ) {
         node_objlimit[ 4 ]++;
      }

      void node_objlimit_error( ) {
         node_objlimit[ 5 ]++;
      }

      void node_objlimit_not_solved_easy( ) {
         node_objlimit_not_solved[ 0 ]++;
      }

      void node_objlimit_not_solved_error( ) {
         node_objlimit_not_solved[ 1 ]++;
      }

      void node_objlimit_pseudo_easy( ) {
         node_objlimit_pseudo[ 0 ]++;
      }

      void node_objlimit_pseudo_neumaier( ) {
         node_objlimit_pseudo[ 1 ]++;
      }

      void node_objlimit_pseudo_exact( ) {
         node_objlimit_pseudo[ 2 ]++;
      }

      void node_objlimit_pseudo_error( ) {
         node_objlimit_pseudo[ 3 ]++;
      }

      void floating_point_lp_error()
      {
         lperror[0]++;
      }

      void exact_lp_error()
      {
         lperror[1]++;
      }

      void print(Message &msg) {
         msg.info("\nStatistics:\n");
         msg.info("Best solution:   {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n", best_solution[ 0 ], "feasible",
                  best_solution[ 1 ], "reconstr", best_solution[ 2 ], "fix-ints", best_solution[ 3 ], "error");
         msg.info("Node feasible:   {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n", node_feasible[ 0 ],
                  "feasible", node_feasible[ 1 ], "reconstr", node_feasible[ 2 ], "factor", node_feasible[ 3 ], "exact",
                  node_feasible[ 4 ], "error");
         msg.info("Node deletion:   {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n", node_deletion[ 0 ], "soplex",
                  node_deletion[ 1 ], "NeuShb", node_deletion[ 2 ], "reconstr", node_deletion[ 3 ], "factor", node_deletion[ 4 ],
                  "exact", node_deletion[ 5 ], "error");
         msg.info("Node del-pseudo: {:>8} {:>10} {:>8} {:>8} \n", node_deletion_pseudo[ 0 ], "feasible",
                  node_deletion_pseudo[ 1 ], "error");
         msg.info("Node infeasible: {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n", node_infeasible[ 0 ],
                  "infeasible", node_infeasible[ 1 ], "NeuShb", node_infeasible[ 2 ], "reconstr", node_infeasible[ 3 ],
                  "exact", node_infeasible[ 4 ], "error");
         msg.info("Node objlimit:   {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n",
                  node_objlimit[ 0 ], "feasible", node_objlimit[ 1 ], "NeuShb", node_objlimit[ 2 ], "reconstr",
                  node_objlimit[ 3 ], "factor", node_objlimit[ 4 ], "exact", node_objlimit[ 5 ], "error");
         msg.info("Node not solved: {:>8} {:>10} {:>8} {:>8} \n", node_objlimit_not_solved[ 0 ], "pseudo",
                  node_objlimit_not_solved[ 1 ], "error");
         msg.info("Node objlim-ps:  {:>8} {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} \n", node_objlimit_pseudo[ 0 ], "easy",
                  node_objlimit_pseudo[ 1 ], "NeuShb", node_objlimit_pseudo[ 2 ], "exact",
                  node_objlimit_pseudo[ 3 ], "error");
         msg.info("\nLP Errors: {:>8} {:>10} {:>8} {:>8}", lperror[0], "floating", lperror[1], "exact");

         msg.info("\nTODO:\n");
         msg.info("- set objective limit for SoPlex\n");

      }
   };

} // namespace exact

#endif