// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 6/26/2009
/// @edit Yang Hsia <yhsia@uw.edu> added extra option for mutating a residue to itself (to remove TERcards for RotamerLinks)

// Unit headers
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/MutateResidueCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
//parsing
#include <utility>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


// Utility Headers

// Unit Headers

// C++ headers

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace core::chemical;
using namespace std;

using core::pose::Pose;
using core::conformation::Residue;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MutateResidue" );

std::string
MutateResidueCreator::keyname() const
{
	return MutateResidueCreator::mover_name();
}

protocols::moves::MoverOP
MutateResidueCreator::create_mover() const {
	return protocols::moves::MoverOP( new MutateResidue );
}

std::string
MutateResidueCreator::mover_name()
{
	return "MutateResidue";
}

/// @brief default ctor
MutateResidue::MutateResidue() :
	parent(),
	target_(""),
	res_name_(""),
	preserve_atom_coords_(false),
	mutate_self_(false)
{}

/// @brief copy ctor
MutateResidue::MutateResidue(MutateResidue const& dm) :
	//utility::pointer::ReferenceCount(),
	parent( dm ),
	target_(dm.target_),
	res_name_(dm.res_name_),
	preserve_atom_coords_(dm.preserve_atom_coords_),
	mutate_self_(dm.mutate_self_)
{}

/// @brief Mutate a single residue to a new amino acid
/// @param target The residue index to mutate
/// @param new_res The name of the replacement residue

MutateResidue::MutateResidue( core::Size const target, string new_res ) :
	parent(),
	target_(""),
	res_name_(std::move(new_res)),
	preserve_atom_coords_(false),
	mutate_self_(false)
{
	set_target( target );
}

MutateResidue::MutateResidue( core::Size const target, int const new_res ) :
	parent(),
	target_(""),
	res_name_( name_from_aa( aa_from_oneletter_code( new_res ) ) ),
	preserve_atom_coords_(false),
	mutate_self_(false)
{
	set_target( target );
}

MutateResidue::MutateResidue( core::Size const target, core::chemical::AA const aa) :
	parent(),
	target_(""),
	res_name_( name_from_aa( aa )),
	preserve_atom_coords_(false),
	mutate_self_(false)
{
	set_target( target );
}

/**
* @brief Reinitialize this protocols::moves::Mover with parameters from the specified tags.
* @details Parameters recognized:
*  - target_pdb_num or target_res_num. A single target residue to form disulfides to
*  - target_pdb_nums or target_res_nums. A list of possible target residues
*/
void MutateResidue::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & /*pose*/)
{

	// Set target to the residue specified by "target_pdb_num" or "target_res_num":
	if ( !tag->hasOption("target") ) {
		TR.Error << "Error: no 'target' parameter specified." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("");
	}
	set_target( tag->getOption<string>("target") );
	//set_target( core::pose::parse_resnum( tag->getOption<string>("target"), pose ) );

	//Set the identity of the new residue:
	if ( !tag->hasOption("new_res") ) {
		TR.Error << "Error: no 'new_res' parameter specified." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("");
	}
	set_res_name( tag->getOption<string>("new_res") );

	//set if you want to mutate the residue to itself, default false.
	mutate_self_ = tag->getOption< bool >( "mutate_self", false );

	//Set whether the mover should try to preserve atom XYZ coordinates or not.  (Default false).
	set_preserve_atom_coords( tag->getOption<bool>("preserve_atom_coords", false) );

	return;
}

/// @brief Set this mover's target residue index, based on the Rosetta indexing.
///
void MutateResidue::set_target(core::Size const target_in)
{
	std::stringstream target_in_string;
	target_in_string << target_in;
	target_ = target_in_string.str();
}

void MutateResidue::set_res_name( core::chemical::AA const & aa){
	res_name_ = name_from_aa( aa );
}

void MutateResidue::apply( Pose & pose ) {

	// Converting the target string to target residue index must be done at apply time,
	// since the string might refer to a residue in a reference pose.
	core::Size const rosetta_target( core::pose::parse_resnum( target(), pose, true /*"true" must be specified to check for refpose numbering when parsing the string*/ ) );

	if ( rosetta_target < 1 ) {
		// Do nothing for 0
		return;
	}
	if ( rosetta_target > pose.total_residue() ) {
		TR.Error << "Error: Residue "<< rosetta_target <<" is out of bounds." << std::endl;
		utility_exit();
	}

	if ( mutate_self_ ) {
		TR << "Setting target residue: " << rosetta_target << " to self (" << pose.residue( rosetta_target ).name3() << ")" << std::endl;
		set_res_name( pose.residue( rosetta_target ).name3() ); //sets res_name to the residue of target
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Mutating residue " << rosetta_target << " from "
			<< pose.residue( rosetta_target ).name3() << " to " << res_name_ <<" ." << std::endl;
	}

	chemical::ResidueTypeSetCOP restype_set( pose.residue( rosetta_target ).residue_type_set() );

	// Create the new residue and replace it
	conformation::ResidueOP new_res = conformation::ResidueFactory::create_residue(
		restype_set->name_map(res_name_), pose.residue( rosetta_target ),
		pose.conformation());
	// Make sure we retain as much info from the previous res as possible
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( rosetta_target ),
		*new_res, pose.conformation(), !preserve_atom_coords() );
	pose.replace_residue( rosetta_target, *new_res, false );

}

std::string
MutateResidue::get_name() const {
	return MutateResidueCreator::mover_name();
}






} // moves
} // protocols
