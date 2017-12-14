// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/align/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.align.util" );

using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace align {

///////////////////////////////////////////////////////////////////////
// By using the StepWisePoseAligner,
// this is really careful with changes in virtual atoms, etc.
// And the cross-checks have been very useful at finding bugs in stepwise code.
///////////////////////////////////////////////////////////////////////
core::Real
get_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
	utility::vector1< core::Size > const & calc_rms_res,
	bool const check_align_at_superimpose_res,
	bool const check_switch ) {

	StepWisePoseAligner pose_aligner( pose2 );
	pose_aligner.set_user_defined_calc_rms_res( calc_rms_res );
	if ( check_align_at_superimpose_res ) {
		pose_aligner.set_root_partition_res( modeler::figure_out_root_partition_res( pose2, calc_rms_res ) );
	}
	pose_aligner.initialize( pose1 );
	// HOW did it accept check_align_at_superimpose_res instead of a pose at first??
	core::Real const rmsd = pose_aligner.get_rmsd_no_superimpose( pose1, check_align_at_superimpose_res );
	//TR << "Immediately before problem call to get_rmsd_no_superimpose, check_align is " << check_align_at_superimpose_res << std::endl;
	//core::Real const rmsd = pose_aligner.get_rmsd_no_superimpose( pose1, pose2, check_align_at_superimpose_res );

	if ( check_switch ) {
		core::Real const rmsd_switch = get_rmsd( pose2, pose1, calc_rms_res,
			check_align_at_superimpose_res, false /* check switch */);
		runtime_assert( std::abs( rmsd  - rmsd_switch ) < 1.0e-3 );
	}

	return rmsd;
}

///////////////////////////////////////////////////////////////////////
core::Real
get_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
	bool const check_align_at_superimpose_res,
	bool const check_switch ) {
	utility::vector1< core::Size > blank_calc_rms_res;
	return get_rmsd( pose1, pose2, blank_calc_rms_res,
		check_align_at_superimpose_res, check_switch );
}


////////////////////////////////////////////////////////////////////////////////
void
align_pose_and_add_rmsd_constraints( pose::Pose & pose,
	pose::PoseCOP align_pose,
	utility::vector1< Size > const & moving_res_list,
	Real const rmsd_screen ) {

	if ( align_pose == nullptr ) return;

	utility::vector1< Size > root_partition_res = figure_out_root_partition_res( pose, moving_res_list );
	if ( root_partition_res.size() == 0 ) root_partition_res.push_back( pose.fold_tree().root() );

	// can later generalize to use 'reference_pose', not necessarily align_pose.
	modeler::align::StepWisePoseAligner pose_aligner( *align_pose );
	pose_aligner.set_root_partition_res( root_partition_res );
	Pose pose_save = pose;
	pose_aligner.apply( pose );

	Real const rms_in_superimpose_atoms = core::scoring::rms_at_corresponding_atoms_no_super( pose_save, *align_pose, pose_aligner.superimpose_atom_id_map() );
	TR.Debug << "SUPERIMPOSE RMSD: " << rms_in_superimpose_atoms  << std::endl;
	if ( rms_in_superimpose_atoms < 1.0e-5 ) {
		TR.Debug << "Not reorienting because rmsd_over_alignment_atoms is very small." << std::endl;
		pose = pose_save; // to avoid floating point deviations.
	}

	// Don't do this if we have new_align_pdb!!
	using namespace basic::options;
	if ( !option[ basic::options::OptionKeys::stepwise::new_align_pdb ].user() ) {
		if ( rmsd_screen > 0.0 ) pose_aligner.create_coordinate_constraints( pose, rmsd_screen );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate rmsd for pose and any 'other poses', and add up in quadrature.
Real
superimpose_with_stepwise_aligner( pose::Pose & pose, pose::Pose const & align_pose,
	bool const superimpose_over_all_instantiated /* = false */ ){
	modeler::align::StepWisePoseAligner pose_aligner( align_pose );
	pose_aligner.set_superimpose_over_all_instantiated( superimpose_over_all_instantiated );
	return pose_aligner.get_rmsd_over_all_poses( pose ); // will include "other_poses"
}

///////////////////////////////////////////////////////////////////////////////////////////
// Following functions (superimpose_pose, creat_alignment_id_map) use legacy code
// for choosing which atoms to superimpose on -- but are called by InputStreamWithResidueInfo
// and a couple other classes that should probably ALL BE DEPRECATED.
//
// Almost all (perhaps all) choices in atoms for superimposition are nicely formalized in
//  StepwisePoseAligner -- switch over completely...
//
//   -- rhiju, 2014
///////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Superimpose mod_pose onto ref_pose using the mapping of residues from
/// mod_pose to ref_pose given by res_map. Simple wrapper around superimpose_pose using IDs.
Real
superimpose_pose_legacy(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	std::map< Size, Size > const & res_map )
{
	id::AtomID_Map< id::AtomID > atom_ID_map = create_alignment_id_map_legacy( mod_pose, ref_pose, res_map );
	if ( atom_ID_map.all_default() ) return 0.0;
	return scoring::superimpose_pose( mod_pose, ref_pose, atom_ID_map );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
id::AtomID_Map< id::AtomID >
create_alignment_id_map_legacy( pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & superimpose_res ){

	std::map< core::Size, core::Size > res_map;

	for ( Size seq_num = 1; seq_num <= mod_pose.size(); ++seq_num ) {
		if ( !superimpose_res.has_value(seq_num) ) continue;
		res_map[ seq_num ] = seq_num;
	}

	return create_alignment_id_map_legacy( mod_pose, ref_pose, res_map );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
id::AtomID_Map< id::AtomID >
create_alignment_id_map_legacy( pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::Size, core::Size > res_map ){

	using namespace chemical;
	using namespace protocols::stepwise::modeler::protein;
	using namespace protocols::stepwise::modeler::rna;
	using namespace core::id;

	AtomID_Map< AtomID > atom_ID_map;
	pose::initialize_atomid_map( atom_ID_map, mod_pose, AtomID::BOGUS_ATOM_ID() );

	for ( Size seq_num = 1; seq_num <= mod_pose.size(); ++seq_num ) {
		if ( mod_pose.residue_type( seq_num ).is_RNA() && res_map.find( seq_num ) != res_map.end() && res_map[ seq_num ] > 0 ) {
			// Parin please update this function!!! Can't we just superimpose over C4'?
			setup_suite_atom_id_map( mod_pose, ref_pose, seq_num,  res_map[ seq_num ], atom_ID_map);
		} else if ( mod_pose.residue( seq_num ).is_protein() ) { // superimpose over CA.
			setup_protein_backbone_atom_id_map( mod_pose, ref_pose, seq_num, res_map[ seq_num ], atom_ID_map); // This will superimpose over N, C-alpha, C
		}
	}

	return atom_ID_map;
}


} //align
} //modeler
} //stepwise
} //protocols
