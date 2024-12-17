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

#ifndef _EXACT_DATA_TREE_HPP_
#define _EXACT_DATA_TREE_HPP_

#include "exact/misc/Vec.hpp"
#include "exact/core/PropagationView.hpp"

namespace exact {


   struct TreeNode {
      unsigned int parent;
      unsigned int child1;
      unsigned int child2;
      int vipr_constraint;
      int assumption_vipr_index;
      bool infeasible;
      bool cutoff;
      Rational dual_bound;

      /** debugging*/
      Rational value;
      int varindex;

   public:

      TreeNode(unsigned int _parent) :
            parent(_parent), child1(0), child2(0), vipr_constraint(0), assumption_vipr_index(0), infeasible(true), cutoff(false),
            dual_bound(0),
            value(0), varindex(0) {
      }
   };

   class Tree {


      Vec<TreeNode> nodes { };

   public:

      const static int UNDEFINED = 0;


      Tree( ) {
         nodes.emplace_back(0);
      }



      unsigned int add_child(unsigned int current_node, const Assumption &assumption) {

         unsigned int pointer = current_node;
         assert(nodes[ current_node ].child1 == 0 || nodes[ current_node ].child2 == 0);
         unsigned int pos_child = nodes.size( );
         add_node(current_node);
         if( assumption.is_lower )
            nodes[ pointer ].child2 = pos_child;
         else
            nodes[ pointer ].child1 = pos_child;

         nodes[ pos_child ].value = assumption.value;
         nodes[ pos_child ].varindex = assumption.var_index;
         nodes[ pos_child ].assumption_vipr_index = assumption.vipr_index;
         return pos_child;
      }

      unsigned int get_child(unsigned int parent, bool is_lower) {
         if( !is_lower )
            return nodes[ parent ].child1;
         return nodes[ parent ].child2;
      }

      void set_assumption_index(unsigned int node, int vipr_index) {
         nodes[ node ].assumption_vipr_index = vipr_index;
      }

      bool is_infeasible(unsigned int node) {
         return nodes[ node ].infeasible;
      }

      Rational get_dual_bound(unsigned int node) {
         return nodes[ node ].dual_bound;
      }

      bool is_cutoff(unsigned int node) {
         return nodes[ node ].cutoff;
      }

      int get_assumption_index(unsigned int node) {
         return nodes[ node ].assumption_vipr_index;
      }

      void set_constraint_index(unsigned int node, int vipr_index) {
         nodes[ node ].vipr_constraint = vipr_index;
      }

      unsigned int get_parent(unsigned int node) {
         return nodes[ node ].parent;
      }

      unsigned int get_vipr_index(unsigned int node) {
         return nodes[ node ].vipr_constraint;
      }

      unsigned int is_parent_resolved_and_get_sibling(unsigned int child) {

         unsigned parent = nodes[ child ].parent;
         unsigned child1 = nodes[ parent ].child1;
         unsigned child2 = nodes[ parent ].child2;
         if( child1 != UNDEFINED && child2 != UNDEFINED
             && nodes[ child1 ].vipr_constraint != UNDEFINED && nodes[ child2 ].vipr_constraint != UNDEFINED )
         {
            return child1 == child ? child2 : child1;
         }
         return UNDEFINED;

      }

      void set_dual_bound(unsigned int node, const Rational &dual_bound) {
         nodes[ node ].infeasible = false;
         nodes[ node ].dual_bound = dual_bound;
      }
      
      void resolve_node( unsigned int node)
      {
         unsigned int child1 = nodes[ node ].child1;
         unsigned int child2 = nodes[ node ].child2;

         if( child1 != UNDEFINED)
         {
            nodes[ child1 ].cutoff = true;
            resolve_node(child1);
         }
         if( child2 != UNDEFINED)
         {
            nodes[ child2 ].cutoff = true;
            resolve_node(child2);
         }
      }

      void analyse(const Vec<String> &varnames) {
         for( unsigned  int node = 0; node < nodes.size( ); node++ )
         {
            if( nodes[ node ].vipr_constraint != 0 || nodes[node].cutoff )
               continue;

            if( nodes[ node ].child1 == 0 )
            {
               fmt::print("Child1 missing with parent {} \n", nodes[ node ].parent);
               Vec<unsigned int> parents { };
               unsigned int current_node = node;
               while( true )
               {
                  parents.push_back(current_node);
                  if( nodes[ current_node ].parent == 0 )
                  {
                     parents.push_back(0);
                     break;
                  }
                  current_node = nodes[ current_node ].parent;
               }
               unsigned int i = parents.size( ) - 2;
               while( true )
               {

                  fmt::print("{}", varnames[ nodes[ parents[ i ]].varindex ]);
                  bool lower = nodes[parents[i+1]].child2 == parents[i];
                  fmt::print(" {}= {} -", lower? '>': '<', nodes[ parents[ i ]].value);
                  if( i == 0 )
                     break;
                  i--;
               }
               fmt::print("\n");
            }
            else if( nodes[ node ].child2 == 0 )
            {
               fmt::print("Child2 missing with parent {} \n", nodes[ node ].parent);
               Vec<unsigned int> parents { };
               unsigned int current_node = node;
               while( true )
               {
                  parents.push_back(current_node);
                  if( nodes[ current_node ].parent == 0 )
                  {
                     parents.push_back(0);
                     break;
                  }
                  current_node = nodes[ current_node ].parent;
               }
               unsigned int i = parents.size( ) - 2;
               while( true )
               {
                  fmt::print("{}", varnames[ nodes[ parents[ i ]].varindex ]);
                  bool lower = nodes[parents[i+1]].child2 == parents[i];
                  fmt::print(" {}= {} -", lower? '>': '<', nodes[ parents[ i ]].value);
                  if( i == 0 )
                     break;
                  i--;
               }
               fmt::print("\n");
            }
            else if( nodes[ node ].vipr_constraint == 0 ){
               fmt::print("node {} missing constraint\n", node);
               Vec<unsigned int> parents { };
               unsigned int current_node = node;
               while( true )
               {
                  parents.push_back(current_node);
                  if( nodes[ current_node ].parent == 0 )
                  {
                     parents.push_back(0);
                     break;
                  }
                  current_node = nodes[ current_node ].parent;
               }
               unsigned int i = parents.size( ) - 2;
               while( true )
               {
                  fmt::print("{}", varnames[ nodes[ parents[ i ]].varindex ]);
                  bool lower = nodes[parents[i+1]].child2 == parents[i];
                  fmt::print(" {}= {} -", lower? '>': '<', nodes[ parents[ i ]].value);
                  if( i == 0 )
                     break;
                  i--;
               }
               fmt::print("\n");
            }
         }
      }

      void delete_child(unsigned int parent) {
         //TODO: get rid of them in the list
         if( nodes[parent].child1 != Tree::UNDEFINED)
            nodes[nodes[parent].child1].parent = Tree::UNDEFINED;
         if( nodes[parent].child2 != Tree::UNDEFINED)
            nodes[nodes[parent].child2].parent = Tree::UNDEFINED;
         nodes[parent].child2 = Tree::UNDEFINED;
         nodes[parent].child1 = Tree::UNDEFINED;

      }

   private:
      void add_node(unsigned int parent) {
         nodes.emplace_back(parent);
      }

   };



} // namespace exact

#endif
