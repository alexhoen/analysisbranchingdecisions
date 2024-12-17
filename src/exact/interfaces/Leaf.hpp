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

#ifndef _EXACT_CORE_NODE_INFEASIBLE_HPP_
#define _EXACT_CORE_NODE_INFEASIBLE_HPP_

#include "exact/misc/Vec.hpp"
#include "exact/core/PropagationView.hpp"

namespace exact {


   struct Assumption {
   public:

      int var_index;
      Rational value;
      bool is_lower;
      int vipr_index;

      Assumption(int _index, Rational _value, bool _is_lower)
            : var_index(_index), value(_value), is_lower(_is_lower), vipr_index(-1) {
      }
   };

   struct Propagation {
      int var;
      bool is_upper;
      Rational value;
   };

   enum class NodeStatus {
      Infeasible,
      InfeasibleParent,
      Solved,
      ObjLimit,
      ObjLimitParent,
      LPError,
   };

   struct Leaf{
      NodeStatus status;
      Vec<Assumption> assumptions;
      Vec<Propagation> propagations;
      Vec<Rational> var_values;
      Vec<Rational> row_values;
      Rational rhs;
      int nonzeros;
      Vec<ProbingBoundChgReason> boundchanges;
      bool certified;

   public:

      //TOOD: use address
      Leaf(NodeStatus _status, Vec<Assumption> _assumptions, Vec<Propagation> _propagations, Vec<Rational> _var_values,
           Vec<Rational> _row_values, Rational _rhs, int _nonzeros)
      : status(_status), assumptions(_assumptions), propagations(_propagations), var_values(_var_values), row_values(_row_values), rhs(_rhs),
        nonzeros(_nonzeros), boundchanges { }, certified(_propagations.empty( )) {
      }

      void setBoundchanges(const Vec<ProbingBoundChgReason> &_boundchanges) {
         boundchanges = _boundchanges;
         for( const auto& propagation: propagations )
         {
            bool certified_propagation = false;
            for( const auto& boundchange: boundchanges )
               if( boundchange.var_index == propagation.var )
                  if( propagation.is_upper == boundchange.upper && propagation.value == boundchange.bound )
                  {
                     certified_propagation = true;
                     break;
                  }
            if( !certified_propagation )
               return;
         }
         certified = true;
      }

      void
      setCertified( bool _certified){certified = _certified;}

   };

} // namespace exact

#endif
