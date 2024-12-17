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

#ifndef _EXACT_DATA_BOUND_HPP_
#define _EXACT_DATA_BOUND_HPP_

#include "exact/misc/Flags.hpp"

namespace exact
{

class Bound
{

private:
   Rational bound;
   bool initialized;
   bool neg_inf;

public:

   Bound(): bound(100000), initialized(false), neg_inf(false){};

   void
   update(Rational& new_bound)
   {
      initialized = true;
      bound = new_bound;
   }

   void
   update_neg_infinity()
   {
      neg_inf = true;
   }

   //TODO: move comparison to own function
   Rational
   get()
   {
      return bound;
   }

   bool
   is_initialized()
   {
      return initialized;
   }

   bool
   is_negative_infinity()
   {
      return neg_inf;
   }


};

} // namespace exact

#endif
