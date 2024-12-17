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

#ifndef VIPR_INTERFACES_SCIP_INTERFACE_HPP_
#define VIPR_INTERFACES_SCIP_INTERFACE_HPP_


#include "exact/misc/Vec.hpp"
#include <cassert>
#include <stdexcept>

#include "exact/data/Problem.hpp"
#include "exact/data/Bound.hpp"
#include "exact/interfaces/Leaf.hpp"
#include "exact/interfaces/ScaleInformation.hpp"
#include "exact/data/Tree.hpp"

namespace exact {

   class ViprInterface {

   private:
      std::ofstream proof_out;
      int counter;
      Vec<int> lhs_mapping;
      Vec<int> rhs_mapping;
      Vec<int> ub_mapping;
      Vec<int> lb_mapping;
      Tree tree;

      HashMap<int, Vec<std::pair<int, Rational>>> assumptions_mapping_ub;
      HashMap<int, Vec<std::pair<int, Rational>>> assumptions_mapping_lb;

      String viprpath;
      String temppath;
      int constraints_for_problem = 0;
      const int UNKNOWN = -1;

      Vec<Assumption> prev_assumptions;

   public:

      ViprInterface(const Problem<exact::Rational> &problem) :
            counter(0), tree({ }), prev_assumptions({}) {
         const auto &problem_name = problem.getName( );
         int length = problem_name.length( );
         int ending = 4;
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_ZLIB
         if( problem_name.substr(length - 3) == ".gz" )
            ending = 7;
#endif
#ifdef EXACT_USE_BOOST_IOSTREAMS_WITH_BZIP2
         if( problem_name.substr(length - 4) == ".bz2" )
            ending = 8;
#endif
         viprpath = problem_name.substr(0, length - ending) + ".vipr";
         temppath = problem_name.substr(0, length - ending) + "-temp.vipr";
         proof_out = std::ofstream(viprpath);
         counter = 0;
      }

      void
      load_assumptions(Vec<Assumption> &_assumptions)
      {
         prev_assumptions = _assumptions;
      }

      void
      reset_assumptions()
      {
         prev_assumptions.clear();
      }

      void
      init_vipr_certificate(const Problem<exact::Rational> &problem) {
         int nvars = problem.getNCols( );
         int nboundconss = 0;
         int non_zero_obj = 0;
         int nintvars = 0;
         auto domains = problem.getVariableDomains( );

         lb_mapping.reserve(nvars);
         ub_mapping.reserve(nvars);

         for( int var = 0; var < nvars; var++ )
         {
            if( !domains.flags[ var ].test(ColFlag::kLbInf))
               nboundconss++;
            if( !domains.flags[ var ].test(ColFlag::kUbInf))
               nboundconss++;
            if( domains.flags[ var ].test(ColFlag::kIntegral))
               nintvars++;
            if( problem.getObjective( ).coefficients[ var ] != 0 )
               non_zero_obj++;
         }

         proof_out << "VER 1.0\n";
         proof_out << "VAR " << nvars << "\n";
         for( int var = 0; var < nvars; var++ )
         {
            const char *varname = problem.getVariableNames( )[ var ].c_str( );
            //TODO: SCIPvarSetCertificateIndex(vars[ var ], var);
            //      SCIPvarSetCertificateIndex(SCIPvarGetTransVar(vars[ var ]), var);
            if( strstr(varname, " ") != NULL || strstr(varname, "\t") != NULL || strstr(varname, "\n") != NULL
                || strstr(varname, "\v") != NULL || strstr(varname, "\f") != NULL || strstr(varname, "\r") != NULL )
            {
               SCIPerrorMessage(
                     "Variable name <%s> cannot be printed to certificate file because it contains whitespace.\n",
                     varname);
               return;
            }
            proof_out << varname << "\n";
         }
         proof_out << "INT " << nintvars << "\n";
         for( int var = 0; var < nvars; var++ )
            if( domains.flags[ var ].test(ColFlag::kIntegral))
               proof_out << var << "\n";

         proof_out << "OBJ min\n" << non_zero_obj;

         for( int var = 0; var < nvars; var++ )
            if( problem.getObjective( ).coefficients[ var ] != 0 )
               proof_out << " " << var << " " << problem.getObjective( ).coefficients[ var ];
         proof_out << "\n";
         int ncertcons = 0;
         for( int cons = 0; cons < problem.getNRows( ); cons++ )
         {
            if( !problem.getRowFlags( )[ cons ].test(RowFlag::kLhsInf))
               ncertcons++;
            if( !problem.getRowFlags( )[ cons ].test(RowFlag::kRhsInf))
               ncertcons++;
         }

         proof_out << "CON " << ncertcons + nboundconss << " " << nboundconss << "\n";

         for( int var = 0; var < nvars; var++ )
         {
            if( !domains.flags[ var ].test(ColFlag::kLbInf))
            {
               proof_out << "B" << counter << " G "
                         << domains.lower_bounds[ var ] << " " << 1 << " " << var << " " << 1 << "\n";
               lb_mapping.push_back(counter);
               counter++;
            }
            else
               lb_mapping.push_back(UNKNOWN);
            if( !domains.flags[ var ].test(ColFlag::kUbInf))
            {
               proof_out << "B" << counter << " L "
                         << domains.upper_bounds[ var ] << " " << 1 << " " << var << " " << 1 << "\n";
               ub_mapping.push_back(counter);
               counter++;
            }
            else
               ub_mapping.push_back(UNKNOWN);
         }

         rhs_mapping.reserve(problem.getNRows( ));
         lhs_mapping.reserve(problem.getNRows( ));

         for( int cons = 0; cons < problem.getNRows( ); cons++ )
         {
            auto cons_data = problem.getConstraintMatrix( ).getRowCoefficients(cons);
            const int *indices = cons_data.getIndices( );
            const Rational *values = cons_data.getValues( );
            if( !problem.getRowFlags( )[ cons ].test(RowFlag::kLhsInf))
            {

               proof_out << "C" << counter << " G " << problem.getConstraintMatrix( ).getLeftHandSides( )[ cons ]
                         << " " << cons_data.getLength( );
               for( int i = 0; i < cons_data.getLength( ); i++ )
                  proof_out << " " << indices[ i ] << " " << values[ i ];
               proof_out << "\n";
               lhs_mapping.push_back(counter);
               counter++;
            }
            else
               lhs_mapping.push_back(UNKNOWN);
            if( !problem.getRowFlags( )[ cons ].test(RowFlag::kRhsInf))
            {
               proof_out << "C" << counter << " L " << problem.getConstraintMatrix( ).getRightHandSides( )[ cons ]
                         << " " << cons_data.getLength( );
               for( int i = 0; i < cons_data.getLength( ); i++ )
                  proof_out << " " << indices[ i ] << " " << values[ i ];
               proof_out << "\n";
               rhs_mapping.push_back(counter);
               counter++;
            }
            else
               rhs_mapping.push_back(counter);
         }
         constraints_for_problem = counter;
         proof_out.flush();
      }

      void write_rtp_range(Bound &lower, Bound &upper) {
         proof_out << "RTP range ";
         if(lower.is_negative_infinity())
            proof_out << "-inf" << " ";
         else
            proof_out << lower.get( ) << " ";
         if( upper.is_initialized( ))
            proof_out << upper.get( ) << "\n";
         else
            proof_out << "inf" << "\n";
         proof_out.flush( );
      }

      void write_rtp_infeas( ) {
         proof_out << "RTP infeas \n";
         proof_out.flush();
      }

      void write_assumption(Assumption &assumption) {
         proof_out << "A" << counter << " ";
         if( assumption.is_lower )
            proof_out << "G" << " ";
         else
            proof_out << "L" << " ";
         proof_out << assumption.value << " " << 1 << " ";
         if( assumption.is_lower )
            proof_out << assumption.var_index;
         else
            proof_out << assumption.var_index;
         proof_out << " 1 { asm } -1 \n";

         assert(assumption.vipr_index == -1);
         assumption.vipr_index = counter;
         counter++;
      }

      void
      write_solutions(Vec<std::pair<Vec<Rational>, Rational>> &solutions, bool optimum_found, const Rational& optimal_value) {
         proof_out << "SOL " << solutions.size( ) << "\n";

         for( const auto &solution: solutions )
         {
            int nonzeros = 0;
            for( const auto &val: solution.first )
               if( val != 0 )
                  nonzeros++;
            if( optimum_found && optimal_value == solution.second )
               proof_out << "best " << nonzeros;
            else
               proof_out << "feas " << nonzeros;
            for( unsigned int i = 0; i < solution.first.size( ); i++ )
               if( solution.first[ i ] != 0 )
                  proof_out << " " << i << " " << solution.first[ i ];
            proof_out << "\n";
         }
         proof_out.flush();
      }

      void write_derivations(unsigned int number) {
         proof_out << "DER " << number << "\n";
      }

      int write_boundchange(ProbingBoundChgReason &reason, const Problem<Rational> &problem,
                            Vec<Assumption> &assumptions, Vec<int> &local_lb_mapping,
                            Vec<int> &local_ub_mapping) {
         assert(reason.row_reason != -1);
         auto row = problem.getConstraintMatrix( ).getRowCoefficients(reason.row_reason);
         Rational value = 0;
         int col_index_in_row = -1;
         for( int i = 0; i < row.getLength( ); i++ )
            if( reason.var_index == row.getIndices( )[ i ] )
            {
               col_index_in_row = row.getIndices( )[ i ];
               value = row.getValues( )[ i ];
               break;
            }
         assert(value != 0);

         if( reason.upper )
         {
            proof_out << "PU" << counter << " L " << reason.bound << " 1 " << reason.var_index << " 1 " << " { ";
            if( problem.getColFlags( )[ reason.var_index ].test(ColFlag::kIntegral))
               proof_out << "rnd ";
            else
               proof_out << "lin ";
            proof_out << row.getLength( ) << " ";
            assert(( value > 0 && rhs_mapping[ reason.row_reason ] != -1 ) ||
                   ( value < 0 && lhs_mapping[ reason.row_reason ] != -1 ));
            if( value > 0 )
               proof_out << rhs_mapping[ reason.row_reason ] << " " << 1 / value << " ";
            else
               proof_out << lhs_mapping[ reason.row_reason ] << " " << 1 / value << " ";

            for( int i = 0; i < row.getLength( ); i++ )
            {
               if( col_index_in_row == row.getIndices( )[ i ] )
                  continue;
               int vipr_index;
               if(( value > 0 && row.getValues( )[ i ] > 0 ) ||
                  ( value < 0 && row.getValues( )[ i ] < 0 ))
                  vipr_index = local_lb_mapping[ row.getIndices( )[ i ]];
               else
                  vipr_index = local_ub_mapping[ row.getIndices( )[ i ]];
               assert(vipr_index != -1);
               proof_out << vipr_index << " " << -row.getValues( )[ i ] / value << " ";
            }

            proof_out << " } -1 \n";
         }
         else
         {
            proof_out << "PL" << counter << " G " << reason.bound << " 1 " << reason.var_index << " 1 " << " { ";
            if( problem.getColFlags( )[ reason.var_index ].test(ColFlag::kIntegral))
               proof_out << "rnd ";
            else
               proof_out << "lin ";
            proof_out << row.getLength( ) << " ";
            assert(( value < 0 && rhs_mapping[ reason.row_reason ] != -1 ) ||
                   ( value > 0 && lhs_mapping[ reason.row_reason ] != -1 ));
            if( value > 0 )
               proof_out << lhs_mapping[ reason.row_reason ] << " " << 1 / value << " ";
            else
               proof_out << rhs_mapping[ reason.row_reason ] << " " << 1 / value << " ";
            for( int i = 0; i < row.getLength( ); i++ )
            {
               if( col_index_in_row == row.getIndices( )[ i ] )
                  continue;
               int vipr_index;
               if(( value > 0 && row.getValues( )[ i ] > 0 ) ||
                  ( value < 0 && row.getValues( )[ i ] < 0 ))
                  vipr_index = local_ub_mapping[ row.getIndices( )[ i ]];
               else
                  vipr_index = local_lb_mapping[ row.getIndices( )[ i ]];
               proof_out << vipr_index << " " << -row.getValues( )[ i ] / value << " ";
            }

            proof_out << " } -1 \n";
         }
         proof_out.flush();
         counter++;
         return counter - 1;
      }

      void write_node(NodeStatus status, Vec<Assumption> &assumptions,
                      const Vec<Rational> &var_values, const Vec<Rational> &row_values, const Rational& rhs,
                      int nonzeros, const Problem<Rational> &problem ) {
         write_node(status, assumptions, var_values, row_values, rhs, nonzeros, problem, {});
      }

      void write_node(NodeStatus status, Vec<Assumption> &assumptions,
                      const Vec<Rational> &var_values, const Vec<Rational> &row_values, const Rational& rhs,
                      int nonzeros, const Problem<Rational> &problem, ScaleInformation scaling ) {
         Vec<int> local_ub_mapping { ub_mapping };
         Vec<int> local_lb_mapping { lb_mapping };

         unsigned int current_node = 0;
         for( unsigned int i = 0; i < prev_assumptions.size(); i++ )
         {
            Assumption& assumption = prev_assumptions[i];
            unsigned int childnode = tree.get_child(current_node, assumption.is_lower);
            if( childnode == Tree::UNDEFINED )
            {
               write_assumption(assumption);
               childnode = tree.add_child(current_node, assumption);
            }
            else
            {
               //TODO: maybe move this to the actual solver since this might not be necessary
               if( tree.is_cutoff(childnode))
                  return;
               assumption.vipr_index = tree.get_assumption_index(childnode);
            }
            if( assumption.is_lower )
               local_lb_mapping[ assumption.var_index ] = assumption.vipr_index;
            else
               local_ub_mapping[ assumption.var_index ] = assumption.vipr_index;
            current_node = childnode;
         }

         for( unsigned int i = 0; i < assumptions.size(); i++ )
         {
            // skip last not if cutoff and add it later if necessary
            if( i == assumptions.size( ) - 1
                && ( NodeStatus::ObjLimitParent == status || NodeStatus::InfeasibleParent == status ))
               break;
            Assumption& assumption = assumptions[i];
            unsigned int childnode = tree.get_child(current_node, assumption.is_lower);
            if( childnode == Tree::UNDEFINED )
            {
               write_assumption(assumption);
               childnode = tree.add_child(current_node, assumption);
            }
            else
            {
               //TODO: maybe move this to the actual solver since this might not be necessary
               if( tree.is_cutoff(childnode))
                  return;
               assumption.vipr_index = tree.get_assumption_index(childnode);
            }
            if( assumption.is_lower )
               local_lb_mapping[ assumption.var_index ] = assumption.vipr_index;
            else
               local_ub_mapping[ assumption.var_index ] = assumption.vipr_index;
            current_node = childnode;
         }


//         fmt::print("\tVIPR number {}\n", counter);
         if( status == NodeStatus::Infeasible || status == NodeStatus::InfeasibleParent )
            proof_out << "InfNode";
         else if( status == NodeStatus::Solved )
            proof_out << "DB";
         else if( status == NodeStatus::ObjLimit )
            proof_out << "ObjLimit";
         else if( status == NodeStatus::ObjLimitParent )
            proof_out << "LeafCutoff";
         else
            assert(false);
         proof_out << counter << " G " << rhs;

         if( status == NodeStatus::Infeasible || status == NodeStatus::InfeasibleParent )
            proof_out << " 0 ";
         else
         {
            assert(status == NodeStatus::Solved || status == NodeStatus::ObjLimit ||
                   status == NodeStatus::ObjLimitParent);
            proof_out << " OBJ ";
         }

         if( problem.isObjectiveInteger( ))
            proof_out << "{ rnd ";
         else
            proof_out << "{ lin ";

         proof_out << nonzeros;

         for( unsigned int i = 0; i < lb_mapping.size( ); i++ )
         {
            if( var_values[ i ] != 0 )
            {
               int vipr_bound_index;
               if( var_values[ i ] < 0 )
                  vipr_bound_index = local_ub_mapping[ i ];
               else
                  vipr_bound_index = local_lb_mapping[ i ];
               assert(vipr_bound_index >= 0);
               proof_out << "  " << vipr_bound_index << " ";
               if(scaling.is_scaled)
                  proof_out << var_values[ i ] * scaling.scale;
               else
                  proof_out << var_values[ i ];

            }
         }
         for( unsigned int i = 0; i < rhs_mapping.size( ); i++ )
         {
            if( row_values[ i ] != 0 )
            {
               int vipr_bound_index;
               if( row_values[ i ] < 0 )
                  vipr_bound_index = rhs_mapping[ i ];
               else
                  vipr_bound_index = lhs_mapping[ i ];
               assert(vipr_bound_index >= 0);
               proof_out << "  " << vipr_bound_index << " ";
               if(scaling.is_scaled)
                  proof_out << row_values[ i ] * scaling.scale;
               else
                  proof_out << row_values[ i ];
            }
         }
         proof_out << " } -1\n";
         if( scaling.derive_strongest_bound )
         {
            derive_scaled_bound(problem, scaling.scale, scaling.bound);
         }
         tree.set_constraint_index(current_node, counter);
         if( status == NodeStatus::ObjLimit || status == NodeStatus::Solved )
         {
            if( !scaling.derive_strongest_bound)
               tree.set_dual_bound(current_node, rhs);
            else
            {
               const Rational &bound = scaling.scale * scaling.bound;
               tree.set_dual_bound(current_node, bound);
            }
            perform_unsplit(current_node);
            counter++;
         }
         else if( status == NodeStatus::Infeasible )
         {
            perform_unsplit(current_node);
            counter++;
         }
         else if ( status == NodeStatus::InfeasibleParent)
         {
            unsigned int parent = current_node;
            int sibling = tree.is_parent_resolved_and_get_sibling(current_node);
            /* nodes of the parent does not exist or the sibling is not fully resolved*/
            if( (tree.get_child(parent, true) == Tree::UNDEFINED && tree.get_child(parent, false) == Tree::UNDEFINED) ||
                  tree.get_vipr_index(sibling) == Tree::UNDEFINED )
            {
               tree.set_constraint_index(parent, counter);
               perform_unsplit(parent);
               counter++;
               // "delete" sibling
               if( tree.get_vipr_index(sibling) == Tree::UNDEFINED )
                  tree.resolve_node(parent);
            }
            else
            {
               int db_vipr_id = counter;
               counter++;
               Assumption &assumption = assumptions[ assumptions.size( ) - 1 ];
               write_assumption(assumption);
               current_node = tree.add_child(current_node, assumption);
               tree.set_constraint_index(current_node, db_vipr_id);
               assert(tree.get_vipr_index(sibling) != 0);
               counter--;
               perform_unsplit(current_node);
               counter++;
            }
         }
         else
         {
            assert(status == NodeStatus::ObjLimitParent );
            unsigned int parent = current_node;
            int sibling = tree.is_parent_resolved_and_get_sibling(current_node);
            /* nodes of the parent does not exist or the sibling is not fully resolved*/
            if( (tree.get_child(parent, true) == Tree::UNDEFINED && tree.get_child(parent, false) == Tree::UNDEFINED) ||
                  tree.get_vipr_index(sibling)== Tree::UNDEFINED )
            {
               if( !scaling.derive_strongest_bound)
                  tree.set_dual_bound(current_node, rhs);
               else
               {
                  const Rational &bound = scaling.scale * scaling.bound;
                  tree.set_dual_bound(current_node, bound);
               }
               tree.set_constraint_index(parent, counter);
               perform_unsplit(parent);
               counter++;
               // "delete" sibling
               if( tree.get_vipr_index(sibling) == Tree::UNDEFINED )
                  tree.resolve_node(parent);
            }
            else
            {
               int db_vipr_id = counter;
               counter++;
               Assumption &assumption = assumptions[ assumptions.size( ) - 1 ];
               write_assumption(assumption);
               current_node = tree.add_child(current_node, assumption);
               if( !scaling.derive_strongest_bound)
                  tree.set_dual_bound(current_node, rhs);
               else
               {
                  const Rational &bound = scaling.scale * scaling.bound;
                  tree.set_dual_bound(current_node, bound);
               }
               tree.set_constraint_index(current_node, db_vipr_id);
               assert(tree.get_vipr_index(sibling)!= 0);
               counter--;
               perform_unsplit(current_node);
               counter++;
            }
         }
         proof_out.flush();

      }

      void derive_scaled_bound(const Problem <Rational> &problem, Rational& scale_factor, Rational& bound) {
         int nnz = 0;
         for( const auto &item: problem.getObjective( ).coefficients )
            if(item != 0)
               nnz++;
         proof_out << "Round" << ( counter + 1 ) << " G " << bound << " " << nnz << " ";
         for( int i = 0; i < problem.getObjective( ).coefficients.size( ); i++ )
            if( problem.getObjective( ).coefficients[ i ] != 0 )
               proof_out << " " << i << " " << problem.getObjective( ).coefficients[ i ]/scale_factor << " ";
         proof_out << "{ rnd 1 " << counter << " " << 1 / scale_factor << " } -1\n";
         counter ++;
         proof_out << "Scale" << ( counter + 1 ) << " G " << bound * scale_factor << " OBJ ";
         proof_out << "{ lin 1 " << counter << " " << scale_factor << " } -1\n";
         counter++;
      }

      void setVarIndex(int var_index, int vipr_index, bool upper) {
         if( upper )
            ub_mapping[ var_index ] = vipr_index;
         else
            lb_mapping[ var_index ] = vipr_index;
      }

      Vec<int> &getUbMapping( ) { return ub_mapping; }

      Vec<int> &getLbMapping( ) { return lb_mapping; }

      void analysetree(const Vec<String> &varnames) {
         tree.analyse(varnames);
      }

      void switch_to_tmp( bool append ) {
         proof_out.flush( );
         if( append )
            proof_out = std::ofstream(temppath, std::ios_base::app);
         else
            proof_out = std::ofstream(temppath);

      }

      void switch_to_original( ) {
         proof_out.flush();
         proof_out = std::ofstream(viprpath, std::ios_base::app);
      }

      void merge_tmp_to_certificate( ) {
         proof_out.flush();
         std::ofstream of_a(viprpath, std::ios_base::binary | std::ios_base::app);
         std::ifstream if_b(temppath, std::ios_base::binary);

         of_a.seekp(0, std::ios_base::end);
         of_a << if_b.rdbuf();
         }

      int get_derivation_number(){
         return counter - constraints_for_problem;
      }

      void invert_bounds_disguised_as_neg_constraints(const Problem<exact::Rational> &problem) {
         auto domains = problem.getVariableDomains( );


         for( int cons = 0; cons < problem.getNRows( ); cons++ )
         {
            if( problem.getRowSizes( )[ cons ] != 1 )
               continue;
            auto cons_data = problem.getConstraintMatrix( ).getRowCoefficients(cons);
            const int colindex = cons_data.getIndices( )[0];
            const Rational value = cons_data.getValues( )[0];
            assert(cons_data.getLength( ) == 1);
            if( !problem.getRowFlags( )[ cons ].test(RowFlag::kLhsInf))
            {
               Rational new_value = problem.getConstraintMatrix( ).getLeftHandSides( )[ cons ] / value;
               if( value > 0 &&
                   ( domains.flags[ colindex ].test(ColFlag::kLbInf) ||
                     new_value > domains.lower_bounds[ colindex ] ))
                  lb_mapping[ colindex ] = counter;
               else if( value < 0 && ( domains.flags[ colindex ].test(ColFlag::kUbInf) ||
                                       new_value < domains.upper_bounds[ colindex ] ))
               {
                  proof_out << "NU" << counter << " L " << new_value << " 1 " << colindex << " 1 " << " { lin 1 " << lhs_mapping[cons] << " " << 1/value << "} -1\n";
                  ub_mapping[ colindex ] = counter;
                  counter ++;
               }
            }
            if( !problem.getRowFlags( )[ cons ].test(RowFlag::kRhsInf))
            {
               Rational new_value = problem.getConstraintMatrix( ).getRightHandSides( )[ cons ] / value;
               if( value < 0 &&
                   ( domains.flags[ colindex ].test(ColFlag::kLbInf) ||
                     new_value > domains.lower_bounds[ colindex ] ))
               {
                  proof_out << "NL" << counter << " G " << new_value << " 1 " << colindex << " 1 " << " { lin 1 " << rhs_mapping[cons] << " " << 1/value << "} -1\n";
                  lb_mapping[ colindex ] = counter;
                  counter ++;
               }
               else if( value > 0 && ( domains.flags[ colindex ].test(ColFlag::kUbInf) ||
                                       problem.getConstraintMatrix( ).getRightHandSides( )[ cons ] / value<
                                             domains.upper_bounds[ colindex ] ))
                  ub_mapping[ colindex ] = rhs_mapping[cons];
            }
         }
      }

   private:

      void perform_unsplit(unsigned int node) {
         unsigned int sibling = tree.is_parent_resolved_and_get_sibling(node);
         //TODO: root node
         if( sibling == Tree::UNDEFINED )
            return;
         if( node == 0 )
            return;
         counter++;
         assert(tree.get_assumption_index(sibling) != Tree::UNDEFINED);
         assert(tree.get_assumption_index(node) != Tree::UNDEFINED);
         assert(tree.get_vipr_index(sibling) != Tree::UNDEFINED);
         assert(tree.get_vipr_index(node) != Tree::UNDEFINED);

         proof_out << "UNS" << counter;
         if( tree.is_infeasible(sibling) && tree.is_infeasible(node))
            proof_out << " G 1 0";
         else
         {
            Rational dual_bound;
            if( tree.is_infeasible(sibling) && !tree.is_infeasible(node))
               dual_bound = tree.get_dual_bound(node);
            else if( !tree.is_infeasible(sibling) && tree.is_infeasible(node))
               dual_bound = tree.get_dual_bound(sibling);
            else
            {
               Rational d1 = tree.get_dual_bound(node);
               Rational d2 = tree.get_dual_bound(sibling);
               if( d1 < d2 )
                  dual_bound = d1;
               else
                  dual_bound = d2;
            }
            proof_out << " G " << dual_bound << " OBJ ";
            tree.set_dual_bound(tree.get_parent(node), dual_bound);
         }

         proof_out << " { uns ";
         proof_out << tree.get_vipr_index(sibling) << " " << tree.get_assumption_index(sibling) << " ";
         proof_out << tree.get_vipr_index(node) << " " << tree.get_assumption_index(node) << " } -1\n";
         tree.set_constraint_index(tree.get_parent(node), counter);
         perform_unsplit(tree.get_parent(node));
      }


   };

} // namespace exact

#endif
