// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRProtocol.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolRTMin.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRComparerElecDensDiff.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/pack/rtmin.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

// C++ Headers
#include <string>
#include <sstream>

//Auto Headers
#include <utility/vector1.hh>
using std::string;
using std::stringstream;
using core::Size;
using core::pack::RTMin;
using core::pack::task::PackerTask;
using core::pack::task::PackerTaskOP;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::getScoreFunction;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.moves.RRProtocolRTMin");

RRProtocolRTMin::RRProtocolRTMin() :
	nonideal_(false),
	cartesian_(false)
{}

RRProtocolRTMin::RRProtocolRTMin(RRProtocolRTMin const & src) : RRProtocol(src),
	nonideal_(src.nonideal_),
	cartesian_(src.cartesian_)
{}

RRProtocolRTMin::~RRProtocolRTMin() {}

string
RRProtocolRTMin::get_name() const {
	return "RRProtocolRTMin";
}

string
RRProtocolRTMin::get_parameters() const {
	stringstream params;
	params
		<< "nonideal:" << get_nonideal()
		<< ",cartesian:" << get_cartesian();
	return params.str();
}

void
RRProtocolRTMin::set_nonideal(
	bool setting
) {
	nonideal_ = setting;
}

bool
RRProtocolRTMin::get_nonideal() const {
	return nonideal_;
}

void
RRProtocolRTMin::set_cartesian(
	bool setting
) {
	cartesian_ = setting;
}

bool
RRProtocolRTMin::get_cartesian() const {
	return cartesian_;
}

/// @details For each residue, minimize it, and measure the rotamer
/// compared to where it started
void
RRProtocolRTMin::run(
  RRComparerOP comparer,
  RRReporterOP reporter,
  Pose const & pose,
	ScoreFunction const & score_function,
  PackerTask const & packer_task
) {

	using core::chemical::aa_unk;

	// Assume score_function.setup_for_scoring(pose) has already been called.

	PackerTaskOP one_res_task( packer_task.clone() );
	core::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, score_function, one_res_task );

	RTMin rtmin;
	rtmin.set_nonideal(nonideal_);
	rtmin.set_cartesian(cartesian_);

	// I don't know if rtmin looks at more than pack_residue(..)
	one_res_task->temporarily_fix_everything();
	Pose working_pose = pose;  // deep copy

	// For each residue in the packer task,
	// rtmin residue -> and measure recovery
	for( Size ii = 1; ii <= pose.total_residue(); ++ii ){
		if ( !packer_task.pack_residue(ii) ) continue;

		one_res_task->temporarily_set_pack_residue( ii, true );
    if ( ! packer_task.include_current( ii ) ) {
      // if we're not asking for the input sidechains, then don't use them -- replace the input sidechain with a rotamer that will be sampled inside
      // rotamer trials anyways
			core::pack::rotamer_set::RotamerSetFactory rsf;
			core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( pose.residue( ii ) ));
      rotset->set_resid( ii );
      rotset->build_rotamers( pose, score_function, *one_res_task, packer_neighbor_graph );
      if ( rotset->num_rotamers() > 0 ) {
        working_pose.replace_residue( ii, *rotset->rotamer(1), false );
      }
    }

		rtmin.rtmin( working_pose, score_function, one_res_task );
		measure_rotamer_recovery(
			comparer, reporter,
			pose, working_pose,
			pose.residue(ii), working_pose.residue(ii) );
		working_pose.replace_residue( ii, pose.residue( ii ), false ); // restore the original conformation.
		one_res_task->temporarily_set_pack_residue( ii, false );
	}
}

} // rotamer_recovery
} // protocols
