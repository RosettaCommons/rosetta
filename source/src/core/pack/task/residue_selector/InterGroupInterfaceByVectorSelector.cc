// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/InterGroupByVectorSelector.cc
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/pack/task/residue_selector/InterGroupInterfaceByVectorSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorFactory.hh>
#include <core/pack/task/operation/util/interface_vector_calculate.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <set>
#include <string>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

// derived from base class
InterGroupInterfaceByVectorSelector::InterGroupInterfaceByVectorSelector() :
	cb_dist_cut_( 11.0 ),
	nearby_atom_cut_( 5.5 ),
	vector_angle_cut_( 75.0 ),
	vector_dist_cut_( 9.0 )
{}

InterGroupInterfaceByVectorSelector::~InterGroupInterfaceByVectorSelector() {}

ResidueSubset
InterGroupInterfaceByVectorSelector::apply( core::pose::Pose const & pose ) const
{
	std::set< Size > local_set1, local_set2;

	// collect the residues in set1
	if ( group1_selector_ ) {
		set_from_residue_selector( pose, *group1_selector_, local_set1 );
	} else if ( ! group1_resstring_.empty() ) {
		set_from_residue_list_string( pose, group1_resstring_, local_set1 );
	}

	// collect the residues in set2
	if ( group2_selector_ ) {
		set_from_residue_selector( pose, *group2_selector_, local_set2 );
	} else if ( ! group2_resstring_.empty() ) {
		set_from_residue_list_string( pose, group2_resstring_, local_set2 );
	}

	std::set< Size > const & set1( group1_selector_ || ! group1_resstring_.empty() ? local_set1 : group1_set_ );
	std::set< Size > const & set2( group2_selector_ || ! group2_resstring_.empty() ? local_set2 : group2_set_ );

	ResidueSubset subset = operation::util::calc_interacting_vector(
		pose, set1, set2, cb_dist_cut_, nearby_atom_cut_, vector_angle_cut_, vector_dist_cut_ );
	return subset;
}

/// @brief Initialize this object from a Tag.
/// @throws utility::excn::EXCN_Msg_Exception if neither a grp1_selector nor a grp1_residues option is provided,
/// or if neither a grp2_selector nor a grp2_residues option is provided.
void
InterGroupInterfaceByVectorSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	cb_dist_cut( tag->getOption< Real >( "cb_dist_cut", 11.0 )); // should be bigger than vector_dist_cut
	nearby_atom_cut( tag->getOption< Real >( "nearby_atom_cut", 5.5 ) );
	vector_angle_cut( tag->getOption< Real >( "vector_angle_cut", 75.0 ) );
	vector_dist_cut( tag->getOption< Real >( "vector_dist_cut", 9.0 ) );

	// add selectors from tags if any are present
	utility::vector0< utility::tag::TagCOP > const & subtags = tag->getTags();

	if( subtags.size() == 2 ) {
		ResidueSelectorCOP rs1 = ResidueSelectorFactory::get_instance()->new_residue_selector(
					subtags[0]->getName(),
					subtags[0],
					datamap
				);
		ResidueSelectorCOP rs2 = ResidueSelectorFactory::get_instance()->new_residue_selector(
					subtags[1]->getName(),
					subtags[1],
					datamap
				);
		 group1_selector( rs1 );
		 group2_selector( rs2 );
	} else if( subtags.size() == 0 ) { // all needs to be parsed from options
		std::string grp1_selector_name, grp2_selector_name;
		std::string grp1resstring, grp2resstring;
		if ( tag->hasOption( "grp1_selector" ) ) {
			grp1_selector_name = tag->getOption< std::string >( "grp1_selector" );
		} else if ( tag->hasOption( "grp1_residues" ) ) {
			grp1resstring = tag->getOption< std::string >( "grp1_residues" );
		} else {
			std::string error_message = "InterGroupInterfaceByVectorSelector::parse_my_tag requires either grp1_selector or grp1_residues to be specified\n";
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}

		if ( tag->hasOption( "grp2_selector" )) {
			grp2_selector_name = tag->getOption< std::string >( "grp2_selector" );
		} else if ( tag->hasOption( "grp2_residues" ) ) {
			grp2resstring = tag->getOption< std::string >( "grp2_residues" );
		} else {
			std::string error_message = "InterGroupInterfaceByVectorSelector::parse_my_tag requires either grp2_selector or grp2_residues to be specified\n";
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}

		ResidueSelectorCOP grp1_sel_op, grp2_sel_op;
		if ( ! grp1_selector_name.empty() ) {
			try {
				grp1_sel_op = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", grp1_selector_name );
			} catch ( utility::excn::EXCN_Msg_Exception e ) {
				std::string error_message = "Failed to find ResidueSelector named '" + grp1_selector_name + "' from the Datamap from InterGroupInterfaceByVectorSelector::parse_my_tag\n" + e.msg();
				throw utility::excn::EXCN_Msg_Exception( error_message );
			}
		}

		if ( ! grp2_selector_name.empty() ) {
			try {
				grp2_sel_op = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", grp2_selector_name );
			} catch ( utility::excn::EXCN_Msg_Exception e ) {
				std::string error_message = "Failed to find ResidueSelector named '" + grp2_selector_name + "' from the Datamap from InterGroupInterfaceByVectorSelector::parse_my_tag\n" + e.msg();
				throw utility::excn::EXCN_Msg_Exception( error_message );
			}
		}

		if ( grp1_sel_op ) { group1_selector( grp1_sel_op ); }
		else { group1_resstring( grp1resstring ); }

		if ( grp2_sel_op ) { group2_selector( grp2_sel_op ); }
		else { group2_resstring( grp2resstring ); }
	} else { // weird number of subtags
		throw utility::excn::EXCN_Msg_Exception( "InterGroupInterfaceByVectorSelector takes either two or zero subtags to specify residue groups!\n" );
	}
}

ResidueSelectorCOP
InterGroupInterfaceByVectorSelector::group1_selector() const { return group1_selector_; }

/// @details this has the side-effect of clearing out the group1_set_, since the ResidueSelector
/// and std:set schemes are meant to be mutually exclusive.
void
InterGroupInterfaceByVectorSelector::group1_selector( ResidueSelectorCOP setting ) {
	group1_selector_ = setting;
	group1_set_.clear();
	group1_resstring_ = "";
}

std::set< Size > const &
InterGroupInterfaceByVectorSelector::group1_set() const { return group1_set_; }

/// @brief This has the side-effect of resetting the group1_selector_ pointer to 0
/// since the ResidueSelector and explicit std::set schemes are mutually exclusive
void
InterGroupInterfaceByVectorSelector::group1_set( std::set< Size > const & setting )
{
	group1_set_ = setting;
	group1_selector_.reset();
	group1_resstring_ = "";
}

std::string const & InterGroupInterfaceByVectorSelector::group1_resstring() const { return group1_resstring_; }

void InterGroupInterfaceByVectorSelector::group1_resstring( std::string const & setting ) {
	group1_resstring_ = setting;
	group1_set_.clear();
	group1_selector_.reset();
}


ResidueSelectorCOP
InterGroupInterfaceByVectorSelector::group2_selector() const { return group2_selector_; }

/// @details this has the side-effect of clearing out the group2_set_, since the ResidueSelector
/// and std:set schemes are meant to be mutually exclusive.
void
InterGroupInterfaceByVectorSelector::group2_selector( ResidueSelectorCOP setting )
{
	group2_selector_ = setting;
	group2_set_.clear();
	group2_resstring_ = "";
}

std::set< Size > const &
InterGroupInterfaceByVectorSelector::group2_set() const { return group2_set_; }

/// @brief This has the side-effect of resetting the group2_selector_ pointer to 0
/// since the ResidueSelector and explicit std::set schemes are mutually exclusive.
void
InterGroupInterfaceByVectorSelector::group2_set( std::set< Size > const & setting )
{
	group2_set_ = setting;
	group2_selector_.reset();
	group2_resstring_ = "";
}

std::string const & InterGroupInterfaceByVectorSelector::group2_resstring() const { return group2_resstring_; }

void InterGroupInterfaceByVectorSelector::group2_resstring( std::string const & setting ) {
	group2_resstring_ = setting;
	group2_set_.clear();
	group2_selector_.reset();
}

Real InterGroupInterfaceByVectorSelector::cb_dist_cut() const { return cb_dist_cut_; }
Real InterGroupInterfaceByVectorSelector::nearby_atom_cut() const { return nearby_atom_cut_; }
Real InterGroupInterfaceByVectorSelector::vector_angle_cut() const { return vector_angle_cut_; }
Real InterGroupInterfaceByVectorSelector::vector_dist_cut() const { return vector_dist_cut_; }

void InterGroupInterfaceByVectorSelector::cb_dist_cut( Real setting ) { cb_dist_cut_ = setting; }
void InterGroupInterfaceByVectorSelector::nearby_atom_cut( Real setting ) { nearby_atom_cut_ = setting; }
void InterGroupInterfaceByVectorSelector::vector_angle_cut( Real setting ) { vector_angle_cut_ = setting; }
void InterGroupInterfaceByVectorSelector::vector_dist_cut( Real setting ) { vector_dist_cut_ = setting; }


std::string InterGroupInterfaceByVectorSelector::get_name() const {
	return InterGroupInterfaceByVectorSelector::class_name();
}

std::string InterGroupInterfaceByVectorSelector::class_name() {
	return "InterfaceByVector";
}

void
InterGroupInterfaceByVectorSelector::set_from_residue_selector(
	core::pose::Pose const & pose,
	ResidueSelector const & selector,
	std::set< Size > & subset
) const
{
	subset.clear();
	ResidueSubset ressubset = selector.apply( pose );
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ressubset[ ii ] ) {
			subset.insert( ii );
		}
	}
}

void
InterGroupInterfaceByVectorSelector::set_from_residue_list_string(
	core::pose::Pose const & pose,
	std::string const & res_list_string,
	std::set< Size > & subset
) const
{
	subset = pose::get_resnum_list( res_list_string, pose );
}

ResidueSelectorOP
InterGroupInterfaceByVectorSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new InterGroupInterfaceByVectorSelector );
}

std::string
InterGroupInterfaceByVectorSelectorCreator::keyname() const {
	return InterGroupInterfaceByVectorSelector::class_name();
}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
