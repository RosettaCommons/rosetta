#ifndef INCLUDED_ObjexxFCL_ObserverMediator_hh
#define INCLUDED_ObjexxFCL_ObserverMediator_hh


// ObserverMediator: Observer Mediator Functions
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>


namespace ObjexxFCL {
namespace internal {
namespace ObserverMediator {


// Types
typedef  SetWrapper< Observer * >  Observers;


/// @brief Notify Observers About Change in a Subject
void
notify( Subject const & s );


/// @brief Acyclic After Adding a Subject-Observer Relation?
bool
acyclic( Subject const & s, Observer & o );


/// @brief Accumulate a Subject's Observers into Accumulated Observers and Recurse: Return Acyclicity
bool
accumulate( Subject const & s_root, Subject const & s, Observers & accum_observers );


} // namespace ObserverMediator
} // namespace internal
} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ObserverMediator_HH
