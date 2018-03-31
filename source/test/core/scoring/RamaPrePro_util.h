// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/RamaPrePro.cxxtest.hh
/// @brief  Unit tests for the RamaPrePro energy.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef test_core_scoring_RamaPrePro_util_h
#define test_core_scoring_RamaPrePro_util_h

// Test headers
#include <cxxtest/TestSuite.h>

// Project Headers
#include <core/scoring/RamaPrePro.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

// Protocol Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <sstream>

static basic::Tracer TR_util("RamaPrePro_util");

inline
void do_ramaprepro_test(
	std::string const &seq,
	bool const /*add_nmethyl*/,
	bool const flip_chirality
) {
	core::pose::PoseOP pose( new core::pose::Pose );
	core::pose::make_pose_from_sequence(*pose, seq, "fa_standard", false);

	if ( flip_chirality ) {
		protocols::cyclic_peptide::FlipChiralityMover flipper;
		flipper.apply( *pose );
	}

	core::scoring::RamaPrePro const & rama( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );
	if( pose->residue(2).mainchain_torsions().size() == 3 ) {
		TR_util << "\nPHI\tPSI\n";
	} else if ( pose->residue(2).mainchain_torsions().size() == 5 ) {
		TR_util << "\nPHI\tTHETA\tPSI\tMU";
	}
	core::Size leftcount(0);
	for ( core::Size i=1; i<=1000; ++i ) {
		utility::vector1 < core::Real > tors;
		rama.random_mainchain_torsions( pose->conformation(), pose->residue_type_ptr(2), pose->residue_type_ptr(3), tors);
		for(core::Size i(1), imax(tors.size()); i<=imax; ++i ) {
			TR_util << tors[i];
			if(i<imax) TR_util << "\t";
		}
		TR_util << "\n";
		if ( tors[1] <= 0 ) ++leftcount;
	}
	TR_util << "LEFT: " << leftcount << "\tRIGHT: " << 1000-leftcount << std::endl;
	if( basic::options::option[ basic::options::OptionKeys::score::symmetric_gly_tables ]() && pose->residue_type(2).is_achiral_backbone() ) {
		TS_ASSERT( leftcount < 600 && leftcount > 400 );
	} else {
		if ( flip_chirality ) {
			TS_ASSERT( leftcount < 500 );
		} else {
			TS_ASSERT( leftcount > 500 );
		}
	}
}

#endif

