// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Residue.fwd.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_ResidueMatcher_hh
#define INCLUDED_core_conformation_ResidueMatcher_hh


// Project headers
#include <core/conformation/ResidueMatcher.fwd.hh>
#include <core/conformation/Residue.fwd.hh>


// Utility headers
//#include <utility/vector1.fwd.hh>
//#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers


namespace core {
namespace conformation {

class ResidueMatcher : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ResidueMatcher() override;
	virtual
	bool
	operator()( Residue const & rsd1, Residue const & rsd2 ) const = 0;

};

class WatsonCrickResidueMatcher : public ResidueMatcher {
public:

	bool
	operator()( Residue const & rsd1, Residue const & rsd2 ) const override;
};

class ExactResidueMatcher : public ResidueMatcher {
public:

	bool
	operator()( Residue const & rsd1, Residue const & rsd2 ) const override;
};


} // namespace conformation
} // namespace core


#endif // INCLUDED_core_conformation_Residue_FWD_HH
