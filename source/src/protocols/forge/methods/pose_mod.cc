// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/methods/pose_mod.cc
/// @brief methods for Pose modifications
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/methods/pose_mod.hh>

// project headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/OptCysHG.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

// numeric headers

// C++ headers
#include <cmath>
#include <list>
#include <map>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace methods {


static basic::Tracer TR( "protocols.forge.methods.pose_mod" );


/// @brief add cutpoint variants at a specific position
/// @param[in,out] Pose to modify
/// @param[in] position at which to add cutpoint variants
/// @return true if cutpoint variants added, false if position not a topological cutpoint
///  or cutpoint variants already exist
bool
add_cutpoint_variants(
	core::pose::Pose & pose,
	core::Size const pos
)
{ // This function is implemented only because I don't want to
	// waste time changing all the calls for the the existing function
	// in protocols::loops code right now, which, for some reason, takes
	// a Loop object instead of a residue position...
	using core::Size;
	using core::pose::add_variant_type_to_pose_residue;
	using core::chemical::CUTPOINT_LOWER;
	using core::chemical::CUTPOINT_UPPER;

	bool op_performed = false;

	if ( !pose.fold_tree().is_cutpoint( pos ) || pos == 0 || pos >= pose.size() ) {
		return op_performed;
	}

	if ( !pose.residue( pos ).has_variant_type( CUTPOINT_LOWER ) ) {
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, pos );
		op_performed = true;
	}

	Size const next_pos = pos + 1;
	if ( !pose.residue( next_pos ).has_variant_type( CUTPOINT_UPPER ) ) {
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, next_pos );
		op_performed = true;
	}

	if ( op_performed ) {
		pose.conformation().declare_chemical_bond( pos, pose.residue( pos ).atom_name( pose.residue( pos ).upper_connect_atom() ),
			pos + 1, pose.residue( pos + 1 ).atom_name( pose.residue( pos + 1 ).lower_connect_atom() ) );
	}

	return op_performed;
}


/// @brief remove cutpoint variants at a specific position
/// @param[in,out] Pose to modify
/// @param[in] position at which to remove cutpoint variants
/// @return true if cutpoint variants removed, false if no cutpoint variants found
///  or position not a topological cutpoint
bool
remove_cutpoint_variants(
	core::pose::Pose & pose,
	core::Size const pos
)
{
	using core::Size;
	using core::pose::remove_variant_type_from_pose_residue;
	using core::chemical::CUTPOINT_LOWER;
	using core::chemical::CUTPOINT_UPPER;

	if ( !pose.fold_tree().is_cutpoint( pos ) || pos == 0 || pos >= pose.size() ) {
		return false;
	}

	bool op_done = false;

	if ( pose.residue( pos ).has_variant_type( CUTPOINT_LOWER ) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, pos );
		op_done = true;
	}

	Size const next_pos = pos + 1;
	if ( pose.residue( next_pos ).has_variant_type( CUTPOINT_UPPER ) ) {
		core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, next_pos );
		op_done = true;
	}

	// should we remove chemical bond?

	return op_done;
}


/// @brief restore residues (i.e. sidechains)
/// @param[in] old2new map indicating residues to be transferred and
///  the mapping from archive_pose position -> pose position
/// @param[in] archive_pose the original Pose to take residues from
/// @param[out] pose the altered Pose that needs residue restoration
void
restore_residues(
	std::map< core::Size, core::Size > const & old2new,
	core::pose::Pose & archive_pose,
	core::pose::Pose & pose
)
{
	using core::pack::task::TaskFactory;
	using core::pack::task::TaskFactoryOP;
	using core::pack::task::operation::OptCysHG;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;

	using core::pack::pack_rotamers;

	typedef std::map< core::Size, core::Size > Old2New;


	for ( Old2New::const_iterator i = old2new.begin(), ie = old2new.end(); i != ie; ++i ) {

		if ( i->second > pose.size() || i->first > archive_pose.size() ) {
			//TR << "mappign beyond pose length.  only possible in repeat generation!!" << std::endl;
			continue;
		}

		/*
		// check variant type
		core::chemical::ResidueType const & rsd_type (pose.residue(i->second).type());
		core::chemical::ResidueType const & archive_type ( archive_pose.residue(i->first).type());

		if (! rsd_type.variants_match( archive_type )){
		utility::vector1<core::chemical::VariantType> const & variant_types ( rsd_type.variant_types() );
		utility::vector1<core::chemical::VariantType> missing_variant_types;
		for ( utility::vector1<core::chemical::VariantType>::const_iterator it = variant_types.begin(), it_end = variant_types.end(); it != it_end; ++it) {
		if (!archive_type.has_variant_type( *it )) missing_variant_types.push_back(*it);
		}
		for (utility::vector1<core::chemical::VariantType>::const_iterator it = missing_variant_types.begin(), it_end=missing_variant_types.end(); it != it_end; ++it) {
		core::pose::add_variant_type_to_pose_residue( archive_pose, *it, i->first);
		}
		}
		*/
		pose.replace_residue( i->second, archive_pose.residue( i->first ), true );
	}

	//archive_pose.dump_pdb("ARCposeInprogress.pdb");
	//pose.dump_pdb("inProgress.pdb");

	// safety
	pose.energies().clear();

	// (re-)detect disulfides, will convert CYD to CYS if disulfide bond is lost
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	if ( option[ in::detect_disulf ].user() ?
			option[ in::detect_disulf ]() : // detect_disulf true
			pose.is_fullatom() // detect_disulf default but fa pose
			) {
		pose.conformation().detect_disulfides();
	}

	// fix HG of CYS to relieve clashes of any newly converted CYS
	ScoreFunctionOP sfx = core::scoring::get_score_function();
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( core::pack::task::operation::TaskOperationOP( new OptCysHG() ) );
	pack_rotamers( pose, *sfx, tf->create_task_and_apply_taskoperations( pose ) );

	// safety
	pose.energies().clear();
}


/// @brief restore residues (i.e. sidechains)
/// @param[in] archive_pose the original Pose to take residues from
/// @param[out] pose the altered Pose that needs residue restoration
/// @remarks length between two poses must be equal
void
restore_residues(
	core::pose::Pose & archive_pose,
	core::pose::Pose & pose
)
{
	using core::Size;

	debug_assert( archive_pose.size() == pose.size() );

	std::map< Size, Size > old2new;
	for ( Size i = 1, ie = archive_pose.size(); i <= ie; ++i ) {
		old2new[ i ] = i;
	}

	restore_residues( old2new, archive_pose, pose );
}


} // methods
} // forge
} // protocols
