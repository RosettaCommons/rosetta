// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/JumpSampleData.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_JumpSampleData_hh
#define INCLUDED_protocols_abinitio_abscript_JumpSampleData_hh

// Unit Headers
#include <protocols/abinitio/abscript/JumpSampleData.fwd.hh>
#include <protocols/abinitio/abscript/JumpSampleDataCreator.hh>

// Package headers

// Project headers
#include <protocols/jumping/JumpSample.hh>

// Utility Headers
#include <basic/datacache/WriteableCacheableData.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript {

class JumpSampleData : public basic::datacache::WriteableCacheableData {
	typedef basic::datacache::WriteableCacheableData Parent;
	typedef basic::datacache::WriteableCacheableDataOP ParentOP;

public:
	JumpSampleData( std::istream &in );

	JumpSampleData( std::string const & moverkey,
		jumping::JumpSample const & );

	void write( std::ostream &out ) const;

	basic::datacache::CacheableDataOP clone() const;

	std::string datatype() const;

	jumping::JumpSample const& jump_sample() const{
		return jump_sample_;
	}

	std::string const& moverkey() const {
		return moverkey_;
	}

private:
	jumping::JumpSample jump_sample_;
	std::string moverkey_;
};

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_JumpSampleData_hh
