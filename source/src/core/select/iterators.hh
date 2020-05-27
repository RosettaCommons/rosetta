// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/iterators.hh
/// @brief  Pose iterator utilities
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_select_iterators_hh
#define INCLUDED_core_select_iterators_hh

#include <utility/iterators.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <utility/vector1.hh>

namespace core {
namespace select {

///@details USAGE: for( core::Size const resid : resids( pose ) )
inline
utility::SimpleRange1
resids( core::pose::Pose const & pose ){
	return utility::enumerate1( pose.size() );
}


class Enumerate1WithSelector {
public:
	Enumerate1WithSelector(
		core::Size r,
		core::select::residue_selector::ResidueSubsetOP const & selected_residues
	) :
		res_( r ),
		selected_residues_( selected_residues )
	{
		advance_if_needed();
	}

	bool operator == ( Enumerate1WithSelector const & o ) const { return res_ == o.res_; }
	bool operator != ( Enumerate1WithSelector const & o ) const { return res_ != o.res_; }

	core::Size const & operator * () const { return res_; }
	core::Size & operator * () { return res_; }

	void operator++(){
		++res_;
		advance_if_needed();
	}
private:
	void advance_if_needed(){
		//advance to the next selected residue if the current residue is not selected
		//stop if we go too far
		while ( res_ <= selected_residues_->size() && ! (*selected_residues_)[ res_ ] ) {
			++res_;
		}
	}

	core::Size res_ = 0;
	core::select::residue_selector::ResidueSubsetOP selected_residues_;
};


//This method is not optimal but the overhead is expected to be minimial compared to the residue selector overhead.
class PreSelectedResidRange {
public:
	PreSelectedResidRange(
		core::pose::Pose const & pose,
		utility::vector1< select::residue_selector::ResidueSelectorCOP > const & selectors
	) :
		nres_( pose.size() )
	{
		selected_residues_ = utility::pointer::make_shared< core::select::residue_selector::ResidueSubset >( pose.size(), true );

		//Perform AND logic with all selectors
		for ( select::residue_selector::ResidueSelectorCOP const & selector : selectors ) {
			utility::vector1< bool > const selection = selector->apply( pose );
			debug_assert( selection.size() == selected_residues_->size() );
			for ( core::Size const ii : indices1( *selected_residues_ ) ) {
				(*selected_residues_)[ ii ] = (*selected_residues_)[ ii ] && selection[ ii ];
			}
		}
	}

	Enumerate1WithSelector begin() {
		return Enumerate1WithSelector( 1, selected_residues_ );
	}

	Enumerate1WithSelector end() {
		return Enumerate1WithSelector( nres_ + 1, selected_residues_ );
	}


private:
	core::Size nres_ = 0;
	core::select::residue_selector::ResidueSubsetOP selected_residues_;
};


///@details USAGES: for( core::Size const resid : resids( pose, {s1, s2} ) )
inline
PreSelectedResidRange
resids(
	core::pose::Pose const & pose,
	utility::vector1< select::residue_selector::ResidueSelectorCOP > const & selectors
){
	return PreSelectedResidRange( pose, selectors );
}



} // select
} // core

#endif // INCLUDED_core_pose_util_HH
