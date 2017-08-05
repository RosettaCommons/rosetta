// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRProtocolReferenceStructure.cc
/// @brief  Preform the rotamer recovery against a reference structure
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolReferenceStructure.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <sstream>

using std::string;
using std::stringstream;
using core::Size;
using core::pose::Pose;
using core::pose::PoseCOP;
using core::pack::task::PackerTask;
using core::scoring::ScoreFunction;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static THREAD_LOCAL Tracer TR("protocol.rotamer_recovery.RRProtocolReferenceStructure");

RRProtocolReferenceStructure::RRProtocolReferenceStructure() :
	reference_pose_(/* NULL */)
{}

RRProtocolReferenceStructure::RRProtocolReferenceStructure(
	PoseCOP reference_pose
) :
	reference_pose_(std::move(reference_pose))
{}

RRProtocolReferenceStructure::RRProtocolReferenceStructure( RRProtocolReferenceStructure const & src) :
	RRProtocol(),
	reference_pose_(src.reference_pose_)
{}

RRProtocolReferenceStructure::~RRProtocolReferenceStructure() = default;

string
RRProtocolReferenceStructure::get_name() const {
	return "RRProtocolReferenceStructure";
}

string
RRProtocolReferenceStructure::get_parameters() const {
	return "";
}

void
RRProtocolReferenceStructure::reference_structure(
	PoseCOP reference_pose){
	reference_pose_ = reference_pose;
}

/// @details  measure rotamer recovery for each residue
void
RRProtocolReferenceStructure::run(
	RRComparerOP comparer,
	RRReporterOP reporter,
	Pose const & pose,
	ScoreFunction const &,
	PackerTask const & packer_task
) {
	// Assume score_function.setup_for_scoring(pose) has already been called.

	using namespace basic::resource_manager;

	if ( !reference_pose_ ) {
		if ( ResourceManager::get_instance()->
				has_resource_with_description("native") ) {
			reference_pose_ = get_resource< Pose >("native");
		} else {
			stringstream err_msg;
			err_msg
				<< "Attempting to run the Rotamer Recovery against a Reference Structure, "
				<< "but no pose with resource decription could be found.";
			utility_exit_with_message(err_msg.str());
		}
	}


	if ( pose.size() != reference_pose_->size() ) {
		stringstream err_msg;
		err_msg
			<< "Attempting to run the Rotamer Recovery against Reference Structure protocol, "
			<< "but the saved structure has a different number of residues.";
		utility_exit_with_message(err_msg.str());
	}

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( !packer_task.pack_residue(ii) ) continue;
		measure_rotamer_recovery(
			comparer, reporter,
			pose, *reference_pose_,
			pose.residue(ii), reference_pose_->residue(ii) );
	}
}

} // rotamer_recovery
} // protocols
