/*! \file polymorphic.hpp
    \brief Support for pointers to polymorphic base classes
    \ingroup OtherTypes */
/*
  Copyright (c) 2014, Randolph Voorhies, Shane Grant
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither the name of cereal nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL RANDOLPH VOORHIES OR SHANE GRANT BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef CEREAL_TYPES_POLYMORPHIC_FWD_HPP_
#define CEREAL_TYPES_POLYMORPHIC_FWD_HPP_

//! Forces dynamic initialization of polymorphic support in a
//! previously registered source file
/*! @sa CEREAL_REGISTER_DYNAMIC_INIT

    See CEREAL_REGISTER_DYNAMIC_INIT for detailed explanation
    of how this macro should be used.  The name used should
    match that for CEREAL_REGISTER_DYNAMIC_INIT. */
#define CEREAL_FORCE_DYNAMIC_INIT(LibName)              \
  namespace cereal {                                    \
  namespace detail {                                    \
    void dynamic_init_dummy_##LibName();                \
  } /* end detail */                                    \
  namespace {                                           \
    void dynamic_init_##LibName()                       \
    {                                                   \
      ::cereal::detail::dynamic_init_dummy_##LibName(); \
    }                                                   \
  } } /* end namespaces */

#endif
