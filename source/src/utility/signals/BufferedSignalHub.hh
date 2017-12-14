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

/// @file   utility/signals/BufferedSignalHub.hh
/// @brief  SignalHub capable of buffering while blocking.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_utility_signals_BufferedSignalHub_hh
#define INCLUDED_utility_signals_BufferedSignalHub_hh

// type headers
#include <platform/types.hh>

// unit headers
/// #include <core/conformation/signals/GeneralEvent.hh> -- Illegal and unneccessary
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/assert.hh>
/// #include <core/pose/signals/DestructionEvent.hh> -- Illegal and unneccessary
// boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <vector>


namespace utility {
namespace signals {


/// @brief SignalHub capable of buffering while blocking.
/// @details This class is meant for use in implementations of the
///  observer pattern that pass a single Signal object from the subject
///  to its observers.  Signals must be copyable to allow for buffering capability.
///  An observer connects by binding one of its member functions to the
///  BufferedSignalHub.  The member function must be a unary function
///  taking the desired Signal as its argument.  Binding functions that
///  pass references and const Signals will be interpreted correctly even if
///  only 'Signal' is given as the template type.  For example, when using a
///  <tt> BufferedSignalHub< void, Signal > </tt> it's possible to bind the following:
///  <ul>
///       <li> <tt> void f( Signal ) </tt>
///       <li> <tt> void f( Signal const ) </tt>
///       <li> <tt> void f( Signal const & ) </tt>
///  </ul>
///  However, using something like <tt> Signal & </tt> or <tt> Signal const </tt>
///  as the template type will enforce references and/or const-ness.
/// @warning BufferedSignalHub objects are not copyable.
template< typename ReturnType, typename Signal >
class BufferedSignalHub : public SignalHub< ReturnType, Signal > {


private: // typedefs


	typedef SignalHub< ReturnType, Signal > Super;


public: // typedefs


	typedef typename Super::Size Size;
	typedef std::vector< Signal > Buffer;


public: // construct/destruct


	/// @brief default constructor
	inline
	BufferedSignalHub() :
		Super(),
		buffering_( false )
	{}


	/// @brief default destructor
	inline
	virtual
	~BufferedSignalHub() = default;


private: // disallow copy construction and assignment


	/// @brief disallow copy constructor
	BufferedSignalHub( BufferedSignalHub const & rval );


	/// @brief disallow copy assignment
	BufferedSignalHub & operator =( BufferedSignalHub const & rval );


public: // signal management


	/// @brief block signals and buffer them for release upon unblocking
	inline
	void buffer() {
		Super::block();
		buffering_ = true;
	}


	/// @brief allow signals to be sent and release any signals that were buffered
	inline
	virtual
	void unblock() {
		Super::unblock();
		buffering_ = false;
		release_buffer();
	}


	/// @brief are signals being buffered?
	inline
	bool buffering() const {
		return buffering_;
	}


	/// @brief clear the buffer
	inline
	void clear_buffer() {
		buffer_.clear();
	}


	/// @brief number of signals left in the buffer
	inline
	Size buffer_size() const {
		return buffer_.size();
	}


protected: // methods


	/// @brief release the signals in the buffer
	inline
	void release_buffer() {
		for ( typename Buffer::iterator i = buffer_.begin(), ie = buffer_.end(); i != ie; ++i ) {
			this->operator ()( *i );
		}

		clear_buffer();
	}


	/// @brief is signal allowed to be passed?
	/// @details If buffering is enabled, signal will be stored in buffer.
	/// @return true if hub is not blocked, otherwise false
	inline
	virtual
	bool signal_allowed( Signal const s ) {
		if ( buffering_ ) {
			buffer_.push_back( s );
			return false;
		}
		return !Super::blocked();
	}


private: // data


	/// @brief signal buffer, only used when hub is blocked and buffering
	Buffer buffer_;


	/// @brief flag indicating blocked signals should be buffered
	bool buffering_;


};


} // namespace signals
} // namespace utility


#endif // INCLUDED_utility_signals_BufferedSignalHub_HH
