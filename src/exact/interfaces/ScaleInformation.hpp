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

#ifndef _EXACT_CORE_SCALE_INFORMATION_HPP_
#define _EXACT_CORE_SCALE_INFORMATION_HPP_


namespace exact {


   struct ScaleInformation {
   public:

      bool is_scaled;
      Rational scale;
      bool derive_strongest_bound;
      Rational bound;


      ScaleInformation(bool _is_scaled, double _scale)
            : is_scaled(_is_scaled), scale(_scale), derive_strongest_bound(false), bound(0) {
      }


      ScaleInformation(bool _is_scaled, double _scale, Rational& _bound)
            : is_scaled(_is_scaled), scale(_scale), derive_strongest_bound(true), bound(_bound) {
      }

      ScaleInformation()
            : is_scaled(false), scale(0), derive_strongest_bound(false), bound(0) {
      }
   };

} // namespace exact

#endif
