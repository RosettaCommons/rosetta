// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 6/26/2009

// Unit headers
#include <protocols/moves/MutateResidue.hh>
#include <protocols/moves/MutateResidueCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
//parsing
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
//Auto Headers
#include <core/kinematics/Jump.hh>


// Utility Headers

// Unit Headers

// C++ headers

namespace protocols {
namespace moves {

using namespace core;
using namespace core::chemical;
using namespace std;

using core::pose::Pose;
using core::conformation::Residue;

static basic::Tracer TR( "protocols.moves.MutateResidue" );

std::string
MutateResidueCreator::keyname() const
{
	return MutateResidueCreator::mover_name();
}

protocols::moves::MoverOP
MutateResidueCreator::create_mover() const {
	return new MutateResidue;
}

std::string
MutateResidueCreator::mover_name()
{
	return "MutateResidue";
}

///@brief default ctor
MutateResidue::MutateResidue() :
	parent(),
	rb_jump_(1),
	target_(0)
{}

///@brief copy ctor
MutateResidue::MutateResidue(MutateResidue const& dm) :
	//utility::pointer::ReferenceCount(),
	parent( dm ),
	rb_jump_(dm.rb_jump_),
	target_(dm.target_),
	res_name_(dm.res_name_)
{}

///@brief Mutate a single residue to a new amino acid
///@param target The residue index to mutate
///@param new_res The name of the replacement residue
MutateResidue::MutateResidue( Size const target, string const new_res ) :
	parent(),
	rb_jump_(1),
	target_(target),
	res_name_(new_res)
{}

MutateResidue::MutateResidue( Size const target, char const new_res ) :
	parent(),
	rb_jump_(1),
	target_(target),
	res_name_( name_from_aa( aa_from_oneletter_code( new_res ) ) )
{}


void MutateResidue::apply( Pose & pose ) {

	if( target_ < 1 ) {
		// Do nothing for 0
		return;
	}
	if( target_ > pose.total_residue() ) {
		TR.Error << "Error: Residue "<<target_<<" is out of bounds." << std::endl;
		utility_exit();
	}

	TR.Debug << "Mutating residue " << target_ << " from "
		<< pose.residue(target_).name3() << " to " << res_name_ <<" ." << std::endl;

	chemical::ResidueTypeSet const& restype_set( pose.residue(target_).residue_type_set() );

	// Create the new residue and replace it
	conformation::ResidueOP new_res = conformation::ResidueFactory::create_residue(
		restype_set.name_map(res_name_), pose.residue(target_),
		pose.conformation());
	// Make sure we retain as much info from the previous res as possible
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue(target_),
		*new_res, pose.conformation() );
	pose.replace_residue(target_, *new_res, false );

}

std::string
MutateResidue::get_name() const {
	return MutateResidueCreator::mover_name();
}


/**
 * @brief Reinitialize this Mover with parameters from the specified tags.
 * @details Parameters recognized:
 *  - target_pdb_num or target_res_num. A single target residue to form disulfides to
 *  - target_pdb_nums or target_res_nums. A list of possible target residues
 */
void MutateResidue::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & pose)
{

	// Set target to the residue specified by "target_pdb_num" or "target_res_num"
	if( !tag->hasOption("target") ){
		TR.Error << "Error: no 'target' parameter specified." << std::endl;
		utility_exit();
	}
	target_ = protocols::rosetta_scripts::parse_resnum(
		tag->getOption<string>("target"),pose);

	if( !tag->hasOption("new_res") ){
		TR.Error << "Error: no 'new_res' parameter specified." << std::endl;
		utility_exit();
	}
	res_name_ = tag->getOption<string>("new_res");
}


} // moves
} // protocols
