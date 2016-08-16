// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   utility/signals/LinkUnit.hh
/// @brief  class specifying actual link data between observers and subjects
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_utility_signals_LinkUnit_hh
#define INCLUDED_utility_signals_LinkUnit_hh

// unit headers
#include <utility/signals/LinkUnit.fwd.hh>

// package headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/assert.hh>

namespace utility {
namespace signals {


/// @brief struct specifying actual link data between observers and subjects
/// @remarks Used internally by SignalHub.  There is no copy constructor and
///  copy assignment, a LinkUnit is meant to be created only once for each connection.
struct LinkUnit : public utility::pointer::ReferenceCount {


private:


	typedef utility::pointer::ReferenceCount Super;


	/// @brief disallow default constructor
	LinkUnit();


	/// @brief disallow copy constructor
	LinkUnit( LinkUnit const & rval );


	/// @brief disallow copy assignment
	LinkUnit & operator =( LinkUnit const & rval );


public:


	/// @brief function pointer constructor
	/// @param f pointer to function object with proper signature
	inline
	LinkUnit( void * f ) :
		Super(),
		fn( f ),
		valid( true ),
		blocked( false )
	{}


	/// @brief default destructor
	inline
	~LinkUnit() {}


	/// @brief return a reference to the function object with requested cast
	template< typename Function >
	inline
	Function & fref() {
		return ( *static_cast< Function * >( fn ) );
	}


	/// @brief send a signal using the member function
	/// @tparam the boost function type
	/// @tparam the Signal to send, must match the necessary function signature
	template< typename Function, typename Signal >
	inline
	void send( Signal s ) {
		if ( valid && !blocked ) {
			debug_assert( fn );
			fref< Function >()( s );
		}
	}


	/// @brief the observer's member function
	/// @details dependent upon SignalHub to manage memory
	void * fn;


	/// @brief connection still valid?
	bool valid;


	/// @brief connection locally blocked?
	bool blocked;


};


/// @brief struct useful as combined predicate + deallocation of function object
///  contained inside a LinkUnit
/// @details if LinkUnit is no longer valid, deallocates 'fn'
template< typename Function >
struct IsLinkUnitInvalid {
	bool
	operator ()( LinkUnitOP & lu ) const {
		if ( !lu->valid && lu->fn ) {
			delete static_cast< Function * >( lu->fn );
			lu->fn = NULL;
		}
		return !lu->valid;
	}
};


} // namespace signals
} // namespace utility


#endif /* INCLUDED_utility_signals_LinkUnit_HH */
