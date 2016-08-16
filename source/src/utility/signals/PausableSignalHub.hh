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

/// @file   utility/signals/PausableSignalHub.hh
/// @brief  BufferedSignalHub capable of pausing and waiting for stdin.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_utility_signals_PausableSignalHub_hh
#define INCLUDED_utility_signals_PausableSignalHub_hh

// type headers
#include <platform/types.hh>

// unit headers
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/assert.hh>

// boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <iostream>
#include <vector>


namespace utility {
namespace signals {


/// @brief BufferedSignalHub capable of pausing and waiting for stdin.
/// @details When pause is set, a Signal that's passed in to the Hub will
///  first be sent to all observers, and then afterwards the Hub will
///  attempt to read from stdin, waiting for the user to press enter.
///  Please see BufferedSignalHub documentation for other usage.
/// @note If <tt> #define BOINC </tt> pause behavior is non-operative.
/// @warning PausableSignalHub objects are not copyable.
template< typename ReturnType, typename Signal >
class PausableSignalHub : public BufferedSignalHub< ReturnType, Signal > {


private: // typedefs


	typedef BufferedSignalHub< ReturnType, Signal > Super;


public: // typedefs


	typedef typename Super::Size Size;


public: // construct/destruct


	/// @brief default constructor
	inline
	PausableSignalHub() :
		Super(),
		pausing_( false )
	{}


	/// @brief default destructor
	inline
	virtual
	~PausableSignalHub() {}


private: // disallow copy construction and assignment


	/// @brief disallow copy constructor
	PausableSignalHub( PausableSignalHub const & rval );


	/// @brief disallow copy assignment
	PausableSignalHub & operator =( PausableSignalHub const & rval );


public: // signal management


	/// @brief "pause" by waiting for stdin after a signal is sent
	inline
	void pause() {
		pausing_ = true;
	}


	/// @brief do not wait for stdin after a signal is sent
	inline
	void unpause() {
		pausing_ = false;
	}


	/// @brief will a pause occur after sending a signal?
	inline
	bool pausing() const {
		return pausing_;
	}


protected: // methods


	/// @brief wait for stdin after sending a signal
	/// @details user presses enter to continue
	/// @return true if successful, false otherwise
	inline
	virtual
	void after_send() {
#ifndef BOINC
		if ( pausing_ ) {
			char hold;
			std::cin.clear();
			std::cout << "Press enter to continue...\n";
			std::cin >> hold;
		}
#endif
	}


private: // data


	/// @brief flag indicating blocked signals should be buffered
	bool pausing_;


};


} // namespace signals
} // namespace utility


#endif // INCLUDED_utility_signals_PausableSignalHub
