// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.AddChainBreak" );

std::string
AddChainBreakCreator::keyname() const
{
	return AddChainBreakCreator::mover_name();
}

protocols::moves::MoverOP
AddChainBreakCreator::create_mover() const {
	return new AddChainBreak;
}

std::string
AddChainBreakCreator::mover_name()
{
	return "AddChainBreak";
}

AddChainBreak::AddChainBreak() :
	protocols::moves::Mover( AddChainBreakCreator::mover_name() ),
	resnum_( "" )
{}

AddChainBreak::~AddChainBreak() {}

protocols::moves::MoverOP
AddChainBreak::clone() const {
	return (protocols::moves::MoverOP( new AddChainBreak( *this ) ) );
}

void
AddChainBreak::parse_my_tag( TagPtr const tag, DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	/// resnum & pdb_num are now equivalent
	if( tag->hasOption( "resnum" ) )
		resnum_ = tag->getOption< std::string > ("resnum" );
	else if( tag->hasOption( "pdb_num" ) ){
		resnum_ = tag->getOption< std::string > ("pdb_num" );
	}
	TR<<"resnum: "<<resnum_<<std::endl;
}//end parse my tag

void
AddChainBreak::apply( core::pose::Pose & pose )
{
	core::Size const resn( protocols::rosetta_scripts::parse_resnum( resnum(), pose ) );

	using namespace core::chemical;
	using namespace pose;
	core::kinematics::FoldTree f( pose.fold_tree() );
	f.new_jump( resn, resn+1, resn );
	add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, resn );
	add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, resn +1);
	pose.fold_tree( f );
	TR<<"New fold tree: "<<pose.fold_tree()<<std::endl;
}

std::string
AddChainBreak::get_name() const {
	return AddChainBreakCreator::mover_name();
}

} //movers
} //protein_interface_design
} //protocols

