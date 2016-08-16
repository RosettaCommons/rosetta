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

/// @file   utility/signals/TokenHub.hh
/// @brief  BufferedSignalHub capable of passing observers
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_utility_signals_TokenHub_hh
#define INCLUDED_utility_signals_TokenHub_hh

// type headers
#include <platform/types.hh>

// unit headers
#include <utility/signals/TokenHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>

// boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <iostream>
#include <vector>


namespace utility {
namespace signals {


/// @brief BufferedSignalHub capable of passing observers
/// @details If an observer's member function is connected with the 'tokenize'
///  flag, then this function is allowed to permanently pass from a source TokenHub to
///  a destination TokenHub via a call to <tt>receive_tokens_from()</tt>, which will
///  erase the observer function from the source hub.  The observer member function,
///  in effect, becomes a single token which is intended to be passed from Subject
///  to Subject upon copy (see remarks in <tt>receive_tokens_from()</tt> for different
///  behavior in copy construct vs copy assign).  Invalidating the Link will stop the
///  token from being passed from Subject to Subject.  Use the token ability with care, as it
///  depends on the programmer to consistently do the right thing and make sure Observer
///  and Subject(s) do not loose sync with each other.  Observers that communicate via a
///  tokenized member function are recommended to follow the following guidelines:
///  <ul>
///    <li> the Observer should be created and never destroyed as long as any Subject objects are in use
///    <li> the Observer should not store a pointer to its Subject
///  </ul>
/// @warning Only use this class when necessary, as usage of the token ability
///  can introduce too much ambiguity into deciding whether an Observer or Subject still
///  exists.  Most users should think about using a BufferedSignalHub instead and managing
///  the connections manually.
template< typename ReturnType, typename Signal >
class TokenHub : public BufferedSignalHub< ReturnType, Signal > {


private: // typedefs


	typedef BufferedSignalHub< ReturnType, Signal > Super;
	typedef typename Super::LinkUnits LinkUnits;


public: // typedefs


	typedef typename Super::Size Size;


public: // construct/destruct


	/// @brief default constructor
	inline
	TokenHub() :
		Super()
	{}


	/// @brief default destructor
	inline
	virtual
	~TokenHub() {
		this->invalidate_all( tokens_ );
	}


private: // disallow copy construction and assignment


	/// @brief disallow copy constructor
	TokenHub( TokenHub const & rval );


	/// @brief disallow copy assignment
	TokenHub & operator =( TokenHub const & rval );


public: // token copy


	/// @brief receive tokenized connections from source and remove them from source
	/// @remarks Use this within either copy constructor or copy assignment of Subject.
	///  Use within a copy constructor will give behavior that tends to keep the token
	///  from descending too far down in the call graph, while use within copy
	///  assignment will give behavior that can drill deep down into the call graph and
	///  has the potential of getting trapped at the lowest level.
	/// @return the number of tokens transferred
	inline
	Size
	receive_tokens_from( TokenHub const & src ) {
		this->remove_invalid( src.tokens_ );
		this->remove_invalid( tokens_ );

		Size const n_to_transfer = src.tokens_.size();
		tokens_.insert( tokens_.end(), src.tokens_.begin(), src.tokens_.end() );
		src.tokens_.clear();

		return n_to_transfer;
	}


public: // connection management


	/// @brief number of tokenized connections
	inline
	Size
	n_tokens() const {
		this->remove_invalid( tokens_ );
		return tokens_.size();
	}


	/// @brief clear only tokenized connections, will invalidate all tokens
	inline
	void
	clear_tokens() {
		invalidate_all( tokens_ );
		tokens_.clear();
	}


	/// @brief connect and tokenize an observer's member function
	/// @tparam MemFn unary member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @return Link that can manage the connection
	template< typename MemFn, typename Ptr >
	inline
	Link connect_tokenize( MemFn fn, Ptr ptr ) {
		return this->connect( fn, ptr, tokens_ );
	}


	/// @brief disconnect an observer's member function
	/// @details will search both regular and tokenized functions
	/// @tparam MemFn unary member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool disconnect( MemFn fn, Ptr ptr ) {
		bool const flag_regular = Super::disconnect( fn, ptr );
		bool const flag_tokenized = Super::disconnect( fn, ptr, tokens_ );
		return flag_regular || flag_tokenized;
	}


	/// @brief disconnect an observer's tokenized member function
	/// @param fn pointer to observer's member function
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool disconnect_tokenize( MemFn fn, Ptr ptr ) {
		return disconnect( fn, ptr, tokens_ );
	}


protected: // virtual methods


	/// @brief send signal to both regular and tokenized connections
	inline
	virtual
	void send_signal( Signal s ) {
		this->send( s, Super::linkunits() );
		this->send( s, tokens_ );
	}


private: // data


	/// @brief tokenized observer function objects
	mutable LinkUnits tokens_;


};


} // namespace signals
} // namespace utility


#endif // INCLUDED_utility_signals_TokenHub
