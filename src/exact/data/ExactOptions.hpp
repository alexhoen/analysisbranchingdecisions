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

#ifndef _EXACT_CORE_PRESOLVE_OPTIONS_HPP_
#define _EXACT_CORE_PRESOLVE_OPTIONS_HPP_

#include "exact/misc/ParameterSet.hpp"
#include <type_traits>

namespace exact
{

struct ExactOptions
{

   int threads = 0;

   bool propagation = true;

   double feastol = 1e-6;

   unsigned int randomseed = 0;

//   double tlim = std::numeric_limits<double>::max();

   void
   addParameters( ParameterSet& paramSet )
   {

      paramSet.addParameter( "solver.seed", "random seed value for MIP Solver",
                             randomseed );
      paramSet.addParameter( "propagation", "perform propagation",
                             propagation );
//      paramSet.addParameter( "presolve.tlim", "time limit for presolve", tlim,
//                             0.0 );
      paramSet.addParameter( "presolve.threads",
                             "maximal number of threads to use (0: automatic)",
                             threads, 0 );
      paramSet.addParameter( "solver.feastol",
                             "feasible tolerance of the MIP solver",
                             feastol, 0, 1e-1 );
   }


};

} // namespace exact

#endif
