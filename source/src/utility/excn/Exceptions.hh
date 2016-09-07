// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/excn/Exceptions.hh
/// @brief  common derived classes for thrown exceptions
/// @author Oliver Lange


#ifndef INCLUDED_utility_excn_Exceptions_HH
#define INCLUDED_utility_excn_Exceptions_HH


// Unit Headers
#include <utility>
#include <utility/excn/Exceptions.fwd.hh>

#include <utility/excn/EXCN_Base.hh>

// Package Headers
#include <string>
#include <ostream>

namespace utility {
namespace excn {

/* *********************************************************************************************************
************************************************************************************************************
*********************                                                                  *********************
*********************                        W  A  R  N  I  N  G                       *********************
*********************                                                                  *********************
************************************************************************************************************
************************************************************************************************************
Please wait until this note is gone before you start using this interface
so far it is very experimental and should remain fluid.
We will have something definite soon. bug us per email if you need to know details
Oliver <olange@u.washington.edu>
Matthew O'Meara <mattjomeara@gmail.com>


generally:
include files will be found in
/<namespace>/Exceptions.hh
for specialized Exceptions e.g. a EXCN_InvalidFoldTree

all-purpose exceptions are all bundled together in this header.
if this gets to big we will have extra forward declarations in Exceptions.fwd.hh

************************************************************************************************************
************************************************************************************************************
*/

class EXCN_Exception : public EXCN_Base {
public:

};

class EXCN_Msg_Exception : public EXCN_Exception {
public:
	EXCN_Msg_Exception( std::string  msg ) : msg_(std::move( msg )) {};
	void show( std::ostream& ) const override;
	std::string const msg() const override { return msg_; };
	virtual void add_msg( std::string const& str ) {
		msg_ = msg_+"\n"+str;
	}
protected:
	EXCN_Msg_Exception() {};
private:
	std::string msg_;
};

class EXCN_IO : public virtual EXCN_Msg_Exception {
protected:
	EXCN_IO() {};
};

class EXCN_BadInput : public EXCN_IO {
public:
	EXCN_BadInput( std::string const& msg ) : EXCN_Msg_Exception( msg ) {};
protected:
	EXCN_BadInput() {};
private:
};

class EXCN_FileNotFound : public EXCN_IO {
public:
	EXCN_FileNotFound( std::string const& file ) :
		EXCN_Msg_Exception( "unable to open file " + file ), file_( file ) {};
private:
	std::string file_;
};

class EXCN_RangeError : public EXCN_Msg_Exception {
public:
	EXCN_RangeError( std::string const& msg ) :
		EXCN_Msg_Exception( msg ) {};
private:
};

class EXCN_KeyError : public EXCN_Msg_Exception {
public:
	EXCN_KeyError(std::string const & msg) :
		EXCN_Msg_Exception( msg ) {};
private:
};

class EXCN_NullPointer: public EXCN_RangeError {
public:
	EXCN_NullPointer( std::string const& msg ) :
		EXCN_RangeError( msg ) {};
private:
};

class EXCN_RosettaScriptsOption : public EXCN_Msg_Exception {
public:
	EXCN_RosettaScriptsOption( std::string const& msg ) :
		EXCN_Msg_Exception( msg ) {};
private:
};

class EXCN_JD2Failure : public EXCN_Msg_Exception {
public:
	EXCN_JD2Failure( std::string const& msg ) :
		EXCN_Msg_Exception( msg ) {};
private:
};

}
}

#endif
