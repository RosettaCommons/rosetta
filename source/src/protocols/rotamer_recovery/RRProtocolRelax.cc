// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRProtocolRelax.cc
/// @brief  Perform the rotamer recovery test using the FastRelax protocol
/// @author Patrick Conway (ptconway@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolRelax.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

// C++ Headers
#include <string>

//Auto Headers
#include <utility/vector1.hh>

using std::string;
using std::stringstream;
using core::Size;
using core::pack::pack_rotamers;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::get_score_function;
using core::pack::task::PackerTask;
using basic::Tracer;
using core::kinematics::MoveMap;
using core::kinematics::MoveMapOP;
using core::id::DOF_ID;
using core::id::AtomID;
using protocols::relax::RelaxProtocolBaseOP;
using protocols::relax::generate_relax_from_cmd;

namespace protocols {
namespace rotamer_recovery {

static THREAD_LOCAL Tracer TR("protocol.moves.RRProtocolRelax");

RRProtocolRelax::RRProtocolRelax() :
	nonideal_(false),
	cartesian_(false)
{}

RRProtocolRelax::RRProtocolRelax( RRProtocolRelax const & ) = default;

RRProtocolRelax::~RRProtocolRelax() = default;

string
RRProtocolRelax::get_name() const {
	return "RRProtocolRelax";
}

string
RRProtocolRelax::get_parameters() const {
	stringstream params;
	params
		<< "nonideal:" << get_nonideal()
		<< ",cartesian:" << get_cartesian();
	return params.str();
}

void
RRProtocolRelax::set_nonideal(
	bool setting
) {
	nonideal_ = setting;
}

bool
RRProtocolRelax::get_nonideal() const {
	return nonideal_;
}

void
RRProtocolRelax::set_cartesian(
	bool setting
) {
	cartesian_ = setting;
}

bool
RRProtocolRelax::get_cartesian() const {
	return cartesian_;
}


/// @details apply FastRelax and measure rotamer recovery for each residue
void
RRProtocolRelax::run(
	RRComparerOP comparer,
	RRReporterOP reporter,
	Pose const & pose,
	ScoreFunction const &,
	PackerTask const & packer_task
) {
	// Assume score_function.setup_for_scoring(pose) has already been called.

	Pose working_pose = pose; // deep copy
	RelaxProtocolBaseOP relax = generate_relax_from_cmd();
	MoveMapOP movemap( new MoveMap() );

	if ( cartesian_ ) {
		relax->cartesian(true);
	}

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( !packer_task.pack_residue(ii) ) continue;

		// set movemap to allow sidechain DOFs only for residue ii
		core::conformation::Residue res = pose.residue(ii);
		movemap->clear();
		movemap->set_chi( ii, true );
		if ( nonideal_ ) {
			for ( Size jj = 1; jj <= res.natoms(); ++jj ) {
				if ( res.atom_is_backbone(jj) ) continue;
				if ( res.is_virtual(jj) ) continue;
				movemap->set( DOF_ID( AtomID( jj, ii ), core::id::THETA ), true );
				movemap->set( DOF_ID( AtomID( jj, ii ), core::id::D ), true );
			}
		}
		relax->set_movemap( movemap );

		// relax!
		relax->apply(working_pose);

		// measure recovery
		measure_rotamer_recovery(
			comparer, reporter,
			pose, working_pose,
			pose.residue(ii), working_pose.residue(ii) );

		// reset
		working_pose.replace_residue( ii, pose.residue( ii ), false );  // restore the original conformation.
	}
}

} // rotamer_recovery
} // protocols
