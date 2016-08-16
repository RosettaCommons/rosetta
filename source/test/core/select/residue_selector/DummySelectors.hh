// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/DummySelectors.hh
/// @brief  Two simple derived ResidueSelectors to aid in unit testing
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_test_core_pack_task_residue_selector_DummySelectors_HH
#define INCLUDED_test_core_pack_task_residue_selector_DummySelectors_HH

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>

using namespace core::select::residue_selector;

class OddResidueSelector : public ResidueSelector {
public:
	OddResidueSelector() {}
	
	OddResidueSelector(OddResidueSelector const &)
	{}
	
	ResidueSelectorOP clone() const { return ResidueSelectorOP( new OddResidueSelector(*this) ); }

	virtual
	ResidueSubset
	apply( core::pose::Pose const & pose ) const
	{
		ResidueSubset subset( pose.total_residue(), false );
		for ( core::Size ii = 1; ii <= subset.size(); ii += 2 ) {
			subset[ ii ] = true;
		}
		return subset;
	}

	virtual std::string get_name() const { return "Odd"; }
};

class XModYResidueSelector : public ResidueSelector {
public:
	XModYResidueSelector( core::Size x, core::Size y ) :
		x_( x ),
		y_( y )
	{}

	XModYResidueSelector(XModYResidueSelector const &src):
		x_( src.x_ ),
		y_( src.y_ )
	{}
	
	ResidueSelectorOP clone() const { return ResidueSelectorOP( new XModYResidueSelector(*this) ); }


	virtual
	ResidueSubset apply( core::pose::Pose const & pose ) const
	{
		ResidueSubset subset( pose.total_residue(), false );
		std::fill( subset.begin(), subset.end(), false );
		for ( core::Size ii = x_; ii <= subset.size(); ii += y_ ) {
			subset[ ii ] = true;
		}
		return subset;
	}
	virtual std::string get_name() const { return "XModY"; }

private:
	core::Size x_, y_;
};

#endif
