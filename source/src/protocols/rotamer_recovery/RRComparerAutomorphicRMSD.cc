// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRComparerAutomorphicRMSD.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// and apps/pilot/chrisk/rotamer_repack.cc (Chris King)


// Unit Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRComparerAutomorphicRMSD.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <sstream>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>


using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::Real;
using core::chemical::num_canonical_aas;
using core::chemical::ResidueType;
using core::chemical::ResidueTypeSetCOP;
using core::conformation::Residue;
using core::conformation::ResidueOP;
using core::conformation::ResidueFactory;
using core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms;
using core::pose::Pose;
using core::scoring::automorphic_rmsd;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static THREAD_LOCAL Tracer TR("protocol.moves.RRComparerAutomorphicRMSD");

RRComparerAutomorphicRMSD::RRComparerAutomorphicRMSD() :
	include_backbone_atoms_( false ),
	recovery_threshold_( .05 )
{}

RRComparerAutomorphicRMSD::RRComparerAutomorphicRMSD( RRComparerAutomorphicRMSD const & src ) :
	RRComparer(),
	include_backbone_atoms_( src.include_backbone_atoms_ ),
	recovery_threshold_( src.recovery_threshold_ )
{}

RRComparerAutomorphicRMSD::~RRComparerAutomorphicRMSD() = default;

bool
RRComparerAutomorphicRMSD::measure_rotamer_recovery(
	Pose const & pose1,
	Pose const & pose2,
	Residue const & res1,
	Residue const & res2,
	Real & score,
	bool & recovered
) {

	if ( res1.aa() != res2.aa()  || res1.nheavyatoms() != res2.nheavyatoms() ) {
		TR << "Cannot measure rotamer recovery of residue " << res1.seqpos() << " because" << endl;
		TR << "\nresidue 1 has type '" << res1.type().name() << "'" << endl;
		TR << "\nresidue 2 has type '" << res2.type().name() << "'" << endl;
		TR << "\nMake sure the protocol to generate the conformations did not 'design' the sequence identity too." << endl;
		score = -1; recovered = false; return false;
	}

	// TODO: Can this restriction be relaxed? What about using 'is_polymer()'?
	if ( res1.aa() > num_canonical_aas ) {
		TR.Warning << "trying to compare rotamer bins for non-canonical amino acid '" << res1.name() << "'" << endl;
		score = -1; recovered = false; return false;
	}

	if ( get_include_backbone_atoms() ) {
		score = automorphic_rmsd( res1, res2, false /*superimpose*/ );
		recovered = (score <= get_recovery_threshold() );
	} else {
		// We assume here that res1 came from pose1 and res2 comes from pose2 (or at the least their typesets are compatable).
		ResidueTypeSetCOP res1_set( pose1.residue_type_set_for_pose( res1.type().mode() ) );
		ResidueType const & working_res1_type(
			res1_set->get_residue_type_with_variant_added( res1.type(), core::chemical::VIRTUAL_BB ) );
		ResidueOP working_res1 = ResidueFactory::create_residue( working_res1_type );
		copy_residue_coordinates_and_rebuild_missing_atoms(
			res1, *working_res1, pose1.conformation() );

		ResidueTypeSetCOP res2_set( pose2.residue_type_set_for_pose( res2.type().mode() ) );
		ResidueType const & working_res2_type(
			res2_set->get_residue_type_with_variant_added( res2.type(), core::chemical::VIRTUAL_BB ) );
		ResidueOP working_res2 = ResidueFactory::create_residue( working_res2_type );
		copy_residue_coordinates_and_rebuild_missing_atoms(
			res2, *working_res2, pose2.conformation() );

		score = automorphic_rmsd(*working_res1, *working_res2, false /*superimpose*/);
		recovered = (score <= get_recovery_threshold() );
	}
	return true;
}

string
RRComparerAutomorphicRMSD::get_name() const {
	return "RRComparerAutomorphicRMSD";
}

string
RRComparerAutomorphicRMSD::get_parameters() const {
	stringstream ret;
	ret << "include_backbone_atoms:" << get_include_backbone_atoms()
		<< ",recovery_threshold:" << get_recovery_threshold();
	return ret.str();
}

void
RRComparerAutomorphicRMSD::set_include_backbone_atoms(
	bool const include_backbone_atoms
) {
	include_backbone_atoms_ = include_backbone_atoms;
}

bool
RRComparerAutomorphicRMSD::get_include_backbone_atoms() const{
	return include_backbone_atoms_;
}

void
RRComparerAutomorphicRMSD::set_recovery_threshold(
	Real const recovery_threshold
) {
	recovery_threshold_ = recovery_threshold;
}

Real
RRComparerAutomorphicRMSD::get_recovery_threshold() const {
	return recovery_threshold_;
}

void
RRComparerAutomorphicRMSD::set_absolute_threshold(
	Real const absolute_threshold
) {
	absolute_threshold_ = absolute_threshold;
}

Real
RRComparerAutomorphicRMSD::get_absolute_threshold() const {
	return absolute_threshold_;
}

} // rotamer_recovery
} // protocols
