// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/sap/util.cc
/// @brief  Utility functions for SapScore
/// @author Brian Coventry (bcov@uw.edu)

#include <core/pack/guidance_scoreterms/sap/util.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.hh>

#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/init_id_map.hh>
#include <core/scoring/packing/surf_vol.hh>
#include <core/scoring/sasa.hh>
#include <utility/pointer/memory.hh>

#include <basic/Tracer.hh>


#include <chrono>
#include <thread>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

static basic::Tracer TR( "core.pack.guidance_scoreterms.sap.util" );

pack::rotamer_set::RotamerSetsOP
rotamer_sets_from_pose(
	pose::Pose const & pose,
	utility::vector1< core::conformation::ResidueCOP > & res_vector
) {
	res_vector.clear();
	res_vector.reserve( pose.size() );

	// Blank packer task that says we're going to pack every position
	pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose );
	task->restrict_to_repacking();
	pack::rotamer_set::RotamerSetsOP rotsets = utility::pointer::make_shared< pack::rotamer_set::RotamerSets >();
	rotsets->set_task( task );

	for ( core::Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {

		pack::rotamer_set::RotamerSetOP rotset = utility::pointer::make_shared<pack::rotamer_set::RotamerSet_>();
		rotset->set_resid( seqpos );
		rotset->add_rotamer( pose.residue(seqpos) );
		rotsets->set_explicit_rotamers( rotsets->resid_2_moltenres( seqpos ), rotset );

		res_vector.push_back( pose.residue(seqpos).get_self_ptr());
	}
	rotsets->update_offset_data();

	return rotsets;
}


SapConstraintHelperOP
common_setup(
	SapConstraintOptionsOP const & options,
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel,
	utility::vector1< core::conformation::ResidueCOP > & res_vector
) {

	runtime_assert( score_sel );
	options->score_selector( score_sel );
	options->sap_calculate_selector( sap_calculate_sel );
	options->sasa_selector( sasa_sel );

	pack::rotamer_set::RotamerSetsOP rotsets = rotamer_sets_from_pose( pose, res_vector );

	SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( options );
	helper->init_with_pose( pose, *rotsets );

	return helper;
}



Real
calculate_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel
) {
	SapConstraintOptionsOP options = utility::pointer::make_shared<SapConstraintOptions>();
	utility::vector1< core::conformation::ResidueCOP > res_vector;

	SapConstraintHelperOP helper = common_setup( options, pose, score_sel, sap_calculate_sel, sasa_sel, res_vector );

	helper->calculate_energy( res_vector, 0 );

	Real score = helper->set_accurate_sasa_and_recalc( pose );

	return score;
}

core::id::AtomID_Map<Real>
calculate_per_atom_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel
) {
	SapConstraintOptionsOP options = utility::pointer::make_shared<SapConstraintOptions>();
	utility::vector1< core::conformation::ResidueCOP > res_vector;

	SapConstraintHelperOP helper = common_setup( options, pose, score_sel, sap_calculate_sel, sasa_sel, res_vector );

	helper->calculate_energy( res_vector, 0 );

	helper->set_accurate_sasa_and_recalc( pose );

	return helper->get_per_atom_sap( pose );
}

utility::vector1<Real>
calculate_per_res_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel
) {

	core::id::AtomID_Map<Real> atom_saps = calculate_per_atom_sap( pose, score_sel, sap_calculate_sel, sasa_sel );

	utility::vector1<Real> res_saps;
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		core::conformation::Residue const & res = pose.residue(seqpos);

		Real res_sum = 0;
		for ( Size iatom = 1; iatom <= res.natoms(); iatom++ ) {
			if ( atom_saps( seqpos, iatom ) > 0 ) res_sum += atom_saps( seqpos, iatom );
		}
		res_saps.push_back( res_sum );
	}

	return res_saps;
}

Real
calculate_slow_approx_sap(
	pose::Pose const & pose,
	select::residue_selector::ResidueSelectorCOP const & score_sel,
	select::residue_selector::ResidueSelectorCOP const & sap_calculate_sel,
	select::residue_selector::ResidueSelectorCOP const & sasa_sel
) {
	SapConstraintOptionsOP options = utility::pointer::make_shared<SapConstraintOptions>();
	utility::vector1< core::conformation::ResidueCOP > res_vector;

	SapConstraintHelperOP helper = common_setup( options, pose, score_sel, sap_calculate_sel, sasa_sel, res_vector );

	Real score = helper->calculate_energy( res_vector, 0 );

	return score;
}

core::id::AtomID_Map<Real>
sap_atom_sasa(
	core::pose::Pose const & pose,
	select::residue_selector::ResidueSubset const & sasa_sub
) {
	core::id::AtomID_Mask atom_mask;
	core::pose::initialize_atomid_map( atom_mask, pose, false );
	for ( Size seqpos = 1; seqpos <= pose.size(); seqpos++ ) {
		if ( ! sasa_sub[seqpos] ) continue;
		atom_mask.fill_with( seqpos, true );
	}

	core::id::AtomID_Map<Real> atom_sasa;
	if ( core::scoring::packing::surf_vol_available() ) {
		core::scoring::packing::SurfVol surf_vol = core::scoring::packing::get_surf_vol( pose, atom_mask, SAP_PROBE_SIZE );
		atom_sasa = surf_vol.surf;

	} else {
		core::pose::initialize_atomid_map( atom_sasa, pose, Real(0) );
		utility::vector1<Real> rsd_sasa( pose.size() );
		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, SAP_PROBE_SIZE, false, atom_mask );
	}
	return atom_sasa;
}


} //sap
} //guidance_scoreterms
} //pack
} //core
