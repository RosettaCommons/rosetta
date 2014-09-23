// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/signals/SignalHub.hh
/// @brief  primary class for function based implementation of observer pattern
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_utility_signals_SignalHub_hh
#define INCLUDED_utility_signals_SignalHub_hh

// type headers
#include <platform/types.hh>

// unit headers
#include <utility/signals/SignalHub.fwd.hh>

// package headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.hh>

// boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <algorithm>
#include <vector>


namespace utility {
namespace signals {


/// @brief  primary class for function based implementation of observer pattern
/// @details This class is meant for use in implementations of the
///  observer pattern that pass a single Signal object from the subject
///  to its observers.  An observer connects by binding one of its member
///  functions to the SignalHub.  The member function must be a unary
///  function taking the desired Signal as its argument.  Binding functions that
///  pass references and const Signals will be interpreted correctly even if
///  only 'Signal' is given as the template type.  For example, when using a
///  <tt> SignalHub< void, Signal > </tt> it's possible to bind the following:
///  <ul>
///       <li> <tt> void f( Signal ) </tt>
///       <li> <tt> void f( Signal const ) </tt>
///       <li> <tt> void f( Signal const & ) </tt>
///  </ul>
///  However, using something like <tt> Signal & </tt> or <tt> Signal const </tt>
///  as the template type will enforce references and/or const-ness.
/// @warning SignalHub objects are not copyable.
template< typename ReturnType, typename Signal >
class SignalHub {


public: // typedefs


	typedef platform::Size Size;


protected: // typedefs


//	typedef boost::function< ReturnType ( Signal ) > Function; // preferred boost::function syntax
	typedef boost::function1< ReturnType, Signal > Function; // portable boost::function syntax
	typedef std::vector< LinkUnitOP > LinkUnits;
	typedef std::vector< Link > Links;


public: // construct/destruct


	/// @brief default constructor
	inline
	SignalHub() :
		blocked_( false ),
		currently_sending_( false )
	{}


	/// @brief default destructor
	inline
	virtual
	~SignalHub() {
		invalidate_all( links_ );
	}


private: // disallow copy construction and assignment


	/// @brief disallow copy constructor
	SignalHub( SignalHub const & rval );


	/// @brief disallow copy assignment
	SignalHub & operator =( SignalHub const & rval );


public: // operators


	/// @brief send signal to all observers connected to the hub
	inline
	void operator ()( Signal s ) {
		if ( signal_allowed( s ) ) {
			send_signal( s );
			after_send();
		}
	}


public: // methods


	/// @brief number of connections to the hub
	inline
	Size size() const {
		remove_invalid( links_ );
		return links_.size();
	}


public: // signal management


	/// @brief block signals from being sent
	inline
	void block() {
		blocked_ = true;
	}


	/// @brief allow signals to be sent
	/// @remarks override this to allow for custom behavior upon unblocking,
	///  such as buffer release
	inline
	virtual
	void unblock() {
		blocked_ = false;
	}


	/// @brief are signals blocked?
	inline
	bool blocked() const {
		return blocked_;
	}


public: // connection management


	/// @brief clear all connections
	inline
	void clear() {
		invalidate_all( links_ );
		links_.clear();
	}


	/// @brief connect an observer's member function to this hub
	/// @tparam MemFn unary member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @return Link that can manage the connection
	template< typename MemFn, typename Ptr >
	inline
	Link connect( MemFn fn, Ptr ptr ) {
		return connect( fn, ptr, links_ );
	}


	/// @brief disconnect an observer's member function from this hub
	/// @tparam MemFn unary member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool disconnect( MemFn fn, Ptr ptr ) {
		return disconnect( fn, ptr, links_ );
	}


	/// @brief transfer linkunits from source to this SignalHub
	/// @param src The source SignalHub to use.
	/// @return the number of links transferred
	/// @warning When using this function, make sure any observers not using
	///  Link connection management become aware of any changes in subject.
	inline
	Size receive_linkunits_from( SignalHub const & src ) {
		Size const n_to_transfer = src.size(); // also removes invalid links from rval
		remove_invalid( links_ );

		links_.insert( links_.end(), src.links_.begin(), src.links_.end() );
		src.links_.clear();

		return n_to_transfer;
	}


	/// @brief return a copy of the list of links currently being managed by
	///  this SignalHub
	inline
	Links current_links() const {
		Links l;
		l.assign( links_.begin(), links_.end() );
		return l;
	}


protected: // virtual methods for sending signals via operator()


	/// @brief is signal allowed to be passed?
	/// @return true if hub is not blocked, otherwise false
	/// @remarks override this function to get custom signal passing behavior,
	///  such as filtering or buffering
	inline
	virtual
	bool signal_allowed( Signal const ) {
		return !blocked_;
	}


	/// @brief send signal to all observers
	/// @remarks override this function to get custom signal sending behavior
	inline
	virtual
	void send_signal( Signal s ) {
		send( s, links_ );
	}


	/// @brief initiate an action after send (no-op in SignalHub)
	/// @return true
	/// @remarks override this function to get custom behavior after signal
	///  has been sent to all observers
	inline
	virtual
	void after_send() {}


protected: // re-usable link container methods


	/// @brief connect an observer's member function to the given LinkUnits
	/// @tparam MemFn unary member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @param links function object gets added to this set of LinkUnits
	/// @return Link that can manage the connection
	template< typename MemFn, typename Ptr >
	inline
	Link connect( MemFn fn, Ptr ptr, LinkUnits & links ) const {
		using namespace boost;

		remove_invalid( links ); // search and destroy first

		LinkUnits::iterator i = find_connection( bind( fn, ptr, _1 ), links );
		if ( i != links.end() ) {
			return *i; // implicit LinkUnitOP -> Link creation
		}

		links.push_back( utility::pointer::shared_ptr<struct utility::signals::LinkUnit>( new LinkUnit( new Function( bind( fn, ptr, _1 ) ) ) ) );
		return links.back(); // implicit LinkUnitOP -> Link creation
	}


	/// @brief disconnect an observer's member function from the given LinkUnits
	/// @tparam MemFn unary member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @param links function object gets removed from this set of LinkUnits
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool disconnect( MemFn fn, Ptr ptr, LinkUnits & links ) const {
		using namespace boost;

		remove_invalid( links ); // search and destroy first

		LinkUnits::iterator i = find_connection( bind( fn, ptr, _1 ), links );
		if ( i == links.end() ) {
			return false;
		}

		// Only mark and deallocate -- invalid link will automatically be
		// swept up at a later time.  Erasing here will cause an error
		// if an observer needs to disconnect itself while SignalHub is
		// still iterating through its links.
		(*i)->valid = false; // mark unit as false
		deallocate( **i );
		return true;
	}


	/// @brief find a specific connection to this hub
	/// @tparam Functor unary functor storable in Function wrapper and equality comparable with Function wrapper, typically result of boost::bind()
	/// @return iterator if found, otherwise links.end()
	template< typename Functor >
	inline
	LinkUnits::iterator find_connection(
		Functor f,
		LinkUnits & links
	) const
	{
		// Function wrapper objects are not equality comparable due to
		// implementation ambiguities, so instead we run through
		// and compare all functions against the functor.
		for ( LinkUnits::iterator i = links.begin(), ie = links.end(); i != ie ; ++i ) {
			if ( (*i)->valid && (*i)->fref< Function >() == f ) {
				return i;
			}
		}

		return links.end();
	}


	/// @brief send signal to all links taking into account validity and local blocking
	inline
	void send( Signal s, LinkUnits & links ) const {
		if ( links.empty() ) {
			return; // nothing to do
		}

		// first check for and remove all invalid links
		remove_invalid( links );

		// mark beginning of send, prevent reallocation of links_ list
		currently_sending_ = true;

		// the only links left are valid links
		for ( LinkUnits::iterator i = links.begin(), ie = links.end(); i != ie ; ++i ) {
			(*i)->send< Function >( s );
		}

		// mark end of send, allow reallocation of links_ list
		currently_sending_ = false;
	}


	/// @brief invalidate and destroy internals of all links
	inline
	void invalidate_all( LinkUnits & links ) const {
		for ( LinkUnits::iterator i = links.begin(), ie = links.end(); i != ie; ++i ) {
			(*i)->valid = false;
			deallocate( **i );
		}
	}


	/// @brief remove and destroy internals of invalid links
	inline
	void remove_invalid( LinkUnits & links ) const {
		if ( !currently_sending_ ) {
			LinkUnits::iterator i = std::remove_if( links.begin(), links.end(), IsLinkUnitInvalid< Function >() );
			links.erase( i, links.end() );
		}
	}


protected: // LinkUnit memory management


	/// @brief deallocate function stored in a LinkUnit
	void deallocate( LinkUnit & lu ) const {
		if ( lu.fn ) {
			delete static_cast< Function * >( lu.fn );
			lu.fn = NULL;
		}
	}


protected: // access to private data


	/// @brief returns the list of LinkUnits
	inline
	LinkUnits & linkunits() {
		return links_;
	}


private: // data


	/// @brief links specifying connections to observers
	mutable LinkUnits links_;


	/// @brief flag blocking signals from being sent
	bool blocked_;


	/// @brief flag indicating SignalHub is currently in the process of sending
	/// signals; e.g. true will prevent remove_invalid() from causing any
	/// reallocation of the links_ list by erasing elements
	mutable bool currently_sending_;


};


} // namespace signals
} // namespace utility


#endif // INCLUDED_utility_signals_SignalHub_HH
