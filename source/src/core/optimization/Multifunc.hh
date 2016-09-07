// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/Multifunc.hh
/// @brief  Multifunction interface class
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_Multifunc_hh
#define INCLUDED_core_optimization_Multifunc_hh

// Unit headers
#include <core/optimization/Multifunc.fwd.hh>

// Package headers
#include <core/optimization/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {


/// @brief Multifunction interface class
class Multifunc : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;

protected: // Creation


	/// @brief Default constructor
	inline
	Multifunc() : parent()
	{}


	/// @brief Copy constructor
	inline
	Multifunc( Multifunc const & ) : parent()
	{}


public: // Creation


	/// @brief Destructor

	~Multifunc() override
	= default;


protected: // Assignment


	/// @brief Copy assignment
	inline
	Multifunc &
	operator =( Multifunc const & )
	{
		return *this;
	}


public:


	virtual
	Real
	operator ()( Multivec const & phipsi ) const = 0;


	virtual
	void
	dfunc( Multivec const & phipsi, Multivec & dE_dphipsi ) const = 0;

	/// @brief Christophe added the following to allow premature end of minimization
	/// If you want to abort the minimizer under specific circonstances
	/// overload this function and return true if you want to stop, false if you want to continue.
	/// FOR THE MOMENT, ONLY IN DFPMIN!
	virtual
	bool
	abort_min(Multivec const & ) const{
		return false; //By default, we don't abort.
	}
	//End Christophe modifications


	/// @brief Error state reached -- derivative does not match gradient
	/// Derived classes have the oportunity to now output and or analyze the two
	/// vars assignments vars, vars+delta where the derivatives are incorrect.
	virtual
	void
	dump( Multivec const & /*vars*/, Multivec const & /*vars2*/ ) const
	{}


}; // Multifunc


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_Multifunc_HH
