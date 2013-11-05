// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/movers/ShoveResidueMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/ShoveResidueMover.hh>
#include <protocols/protein_interface_design/movers/ShoveResidueMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/protein_interface_design/util.hh>
#include <utility/string_util.hh>

#include <core/graph/Graph.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.ShoveResidueMover" );

std::string
ShoveResidueMoverCreator::keyname() const
{
	return ShoveResidueMoverCreator::mover_name();
}

protocols::moves::MoverOP
ShoveResidueMoverCreator::create_mover() const {
	return new ShoveResidueMover;
}

std::string
ShoveResidueMoverCreator::mover_name()
{
	return "ShoveResidueMover";
}


ShoveResidueMover::ShoveResidueMover() :
	protocols::moves::Mover( ShoveResidueMoverCreator::mover_name() ),
	resnum_( 0 )
{}

ShoveResidueMover::ShoveResidueMover( Size resnum ) :
	protocols::moves::Mover( ShoveResidueMoverCreator::mover_name() ),
	resnum_( resnum )
{}

void
ShoveResidueMover::apply ( pose::Pose & pose )
{
	//using namespace rotamer_set;
	using namespace core::scoring;
	using namespace core::pack::task;
	using namespace core::pack::rotamer_set;
	foreach( core::Size const resid, shove_residues_ ) {
		if(pose.residue(resid).name3() == "GLY")
		{
			//The SHOVE_BB patch does not properly work with Glycine.  Looking at it, it's not clear that it even
			//Makes sense to use this mover with Glycine.  Given this, it's probably best to exit here.
			//If this isn't true, a special case should be added to the patch.  I added this if statement
			//Because running this mover on a glycine causes select_atoms_to_orient() fails during repacking
			// - Sam DeLuca

			utility_exit_with_message("ERROR: ShoveResidueMover does not currently work properly with glycines, and residue " + utility::to_string(resid) + " is a glycine.");

		}
		if ( remove_shove_variant_ ) {
			core::pose::remove_variant_type_from_pose_residue( pose, "SHOVE_BB", resid );
		} else {
			core::pose::add_variant_type_to_pose_residue( pose, "SHOVE_BB", resid );
		}
	}
}

std::string
ShoveResidueMover::get_name() const {
	return ShoveResidueMoverCreator::mover_name();
}

void
ShoveResidueMover::parse_my_tag( TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & pose)
{
	resnum_ = core::pose::get_resnum( tag, pose );
	remove_shove_variant_ = tag->getOption<bool>( "remove_shove_variant", false );
	if( tag->hasOption( "shove" ) ){
  	std::string const shove_val( tag->getOption< std::string >( "shove" ) );
  	utility::vector1< std::string > const shove_keys( utility::string_split( shove_val, ',' ) );
  	foreach( std::string const key, shove_keys ){
			core::Size const resnum( core::pose::parse_resnum( key, pose ) );
			shove_residues_.push_back( resnum );
			TR<<"Using shove atomtype for "<< key <<'\n';
		}
	} else {
		shove_residues_.push_back( resnum_ );
	}
}

} //movers
} //protein_interface_design
} //protocols

