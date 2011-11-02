// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SuperimposeMover.cc
///
/// @brief
/// @author Ingemar Andre

// unit headers
#include <protocols/moves/SuperimposeMover.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

// options
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


namespace protocols {
namespace moves {

static basic::Tracer TR("protocols.moves.SuperimposeMover");

SuperimposeMover::SuperimposeMover() :
	Mover("SuperimposeMover")
{}

SuperimposeMover::SuperimposeMover( Pose const & pose ) :
  Mover("SuperimposeMover"),
	ref_pose_(pose)
	{}

SuperimposeMover::~SuperimposeMover() {}

void
SuperimposeMover::set_reference_pose( Pose const & pose ) {
	ref_pose_ = pose;
}

void
SuperimposeMover::apply( Pose & pose ) {
	using namespace basic::options;

//	if ( option[ OptionKeys::in::file::native ].user() ) {
//		Pose native;
		// For now only align proteins with the same number of residues
	if ( ref_pose_.total_residue() == pose.total_residue() ) {
		core::Real rms = core::scoring::calpha_superimpose_pose( pose, ref_pose_ );
		rms = core::scoring::CA_rmsd( pose, ref_pose_ );
		TR << "superimpose to native. Rms to native: " << rms << std::endl;
	}
}

std::string
SuperimposeMover::get_name() const {
	return "SuperimposeMover";
}


} // moves
} // protocols
