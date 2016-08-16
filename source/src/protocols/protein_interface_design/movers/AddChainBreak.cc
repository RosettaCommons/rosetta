// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/AddChainBreak.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <protocols/protein_interface_design/movers/AddChainBreakCreator.hh>
#include <core/pose/util.hh>
// Package headers
#include <core/chemical/VariantType.hh>
// Project headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <boost/foreach.hpp>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.AddChainBreak" );

std::string
AddChainBreakCreator::keyname() const
{
	return AddChainBreakCreator::mover_name();
}

protocols::moves::MoverOP
AddChainBreakCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddChainBreak );
}

std::string
AddChainBreakCreator::mover_name()
{
	return "AddChainBreak";
}

AddChainBreak::AddChainBreak() :
	protocols::moves::Mover( AddChainBreakCreator::mover_name() ),
	resnum_( "" ),
	change_foldtree_( true ),
	find_automatically_( false ),
	automatic_distance_cutoff_( 2.5 ),
	remove_( false )
{}

AddChainBreak::~AddChainBreak() {}

protocols::moves::MoverOP
AddChainBreak::clone() const {
	return (protocols::moves::MoverOP( new AddChainBreak( *this ) ) );
}

void
AddChainBreak::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	/// resnum & pdb_num are now equivalent
	if ( tag->hasOption( "resnum" ) ) {
		resnum_ = tag->getOption< std::string > ("resnum" );
	} else if ( tag->hasOption( "pdb_num" ) ) {
		resnum_ = tag->getOption< std::string > ("pdb_num" );
	}
	if ( tag->hasOption( "find_automatically" ) ) {
		find_automatically( tag->getOption< bool >( "find_automatically" ) );
		automatic_distance_cutoff( tag->getOption< core::Real >( "distance_cutoff", 2.5 ));
	}
	remove( tag->getOption< bool >( "remove", false ) );
	change_foldtree( tag->getOption< bool >( "change_foldtree", true ) );
	TR<<"resnum: "<<resnum_<<" change foldtree "<<change_foldtree()<<" find cutpoints automatically "<<find_automatically()<<" remove: "<<remove()<<std::endl;
}//end parse my tag

void
AddChainBreak::apply( core::pose::Pose & pose )
{
	using namespace core::chemical;
	using namespace pose;
	core::kinematics::FoldTree f( pose.fold_tree() );

	if ( resnum() != "" ) {
		core::Size const resn( core::pose::parse_resnum( resnum(), pose ) );

		if ( change_foldtree() ) {
			f.new_jump( resn, resn+1, resn );
		}
		if ( pose.residue(resn  ).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(pose,resn);
		if ( pose.residue(resn+1).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(pose,resn+1);
		if ( remove() ) {
			remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, resn );
			remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, resn +1);
		} else {
			add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, resn );
			add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, resn +1);
		}
	}
	utility::vector1< core::Size > cuts;
	cuts.clear();
	if ( find_automatically() ) {
		for ( core::Size i = 1; i < pose.total_residue(); ++i ) {
			core::Real const distance( pose.residue( i ).xyz( "C" ).distance( pose.residue( i + 1 ).xyz( "N" ) ) );
			if ( distance >= automatic_distance_cutoff() ) {
				cuts.push_back( i );
				TR<<"Detecting cut at "<<i<<" with distance "<<distance<<std::endl;
			}
		}
	}
	BOOST_FOREACH ( core::Size const res, cuts ) {
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, res );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, res +1);
		if ( change_foldtree() ) {
			f.new_jump( res, res+1, res );
		}
	}
	if ( change_foldtree() ) {
		pose.fold_tree( f );
		TR<<"New fold tree: "<<pose.fold_tree()<<std::endl;
	}
}

std::string
AddChainBreak::get_name() const {
	return AddChainBreakCreator::mover_name();
}

} //movers
} //protein_interface_design
} //protocols

