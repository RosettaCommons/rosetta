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

/// @file   utility/signals/Link.hh
/// @brief  class storing the link between observers and subjects
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_utility_signals_Link_hh
#define INCLUDED_utility_signals_Link_hh


// unit headers
#include <utility/signals/Link.fwd.hh>

// package headers
#include <utility/signals/LinkUnit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/assert.hh>


namespace utility {
namespace signals {


/// @brief interface wrapper around utility::signals::LinkUnit
/// @details Can be used to manage the link between Subject
///  and Observer.  Signals can be temporarily blocked and the
///  connection can be invalidated (disconnected) by the user.
///  If a Subject is destroyed or a manual disconnect is called
///  via the SignalHub, the Link will automatically become invalidated.
/// @remarks Returned by utility::signals::SignalHub upon connection.
class Link {


private: // typedefs


	typedef utility::pointer::shared_ptr< LinkUnit > LinkUnitOP;


public: // construct/destruct


	/// @brief default constructor
	inline
	Link() {}


	/// @brief LinkUnit constructor
	inline
	Link( LinkUnitOP unit ) :
		unit_( unit ) // must store the same LinkUnit
	{}


	/// @brief copy constructor
	inline
	Link( Link const & rval ) :
		unit_( rval.unit_ ) // must store the same LinkUnit
	{}


	/// @brief default destructor
	inline
	~Link() {}


public: // copy assignment


	/// @brief copy assignment
	Link & operator =( Link const & rval ) {
		if ( this != &rval ) {
			unit_ = rval.unit_; // must store the same LinkUnit
		}

		return *this;
	}


public: // operators


	/// @brief equality comparison
	bool
	operator ==( Link const & rval ) {
		return unit_.get() == rval.unit_.get();
	}


public: // interface


	/// @brief link empty?
	inline
	bool empty() const {
		return ( unit_ == NULL );
	}


	/// @brief link still valid?
	/// @return false if link invalid or empty, otherwise true
	inline
	bool valid() const {
		if ( unit_ ) {
			return unit_->valid;
		}
		return false;
	}


	/// @brief cut the connection, safe to call even when link is empty
	inline
	void invalidate() {
		if ( unit_ ) {
			unit_->valid = false;
		}
	}


	/// @brief link blocked?
	/// @brief true if link is blocked or empty, otherwise false
	inline
	bool blocked() const {
		if ( unit_ ) {
			return unit_->blocked;
		}
		return true;
	}


	/// @brief locally block the link
	inline
	void block() {
		if ( unit_ ) {
			unit_->blocked = true;
		}
	}


	/// @brief locally unblock the link
	inline
	void unblock() {
		if ( unit_ ) {
			unit_->blocked = false;
		}
	}


private: // data


	LinkUnitOP unit_;

};


} // namespace signals
} // namespace utility

#endif /* INCLUDED_utility_signals_Link_HH */
