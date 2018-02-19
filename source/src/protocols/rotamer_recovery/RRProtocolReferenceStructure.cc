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
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

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

static Tracer TR("protocol.rotamer_recovery.RRProtocolReferenceStructure");

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

void
RRProtocolReferenceStructure::parse_attributes_from_tag(
	basic::datacache::DataMap const & datamap,
	utility::tag::TagCOP tag
) {

	if ( ! tag->hasOption( "reference_pdb" ) && ! tag->hasOption( "reference_pose" ) ) {
		std::ostringstream oss;
		oss << "ERROR in parsing attributes from the RRProtocolReferenceStructure: one of either the \"reference_pdb\" or \"reference_pose\" attributes must be given\n";
		CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	if ( tag->hasOption( "reference_pdb" ) && tag->hasOption( "reference_pose" ) ) {
		std::ostringstream oss;
		oss << "ERROR in parsing attributes from the RRProtocolReferenceStructure: only one of the \"reference_pdb\" or \"reference_pose\" attributes may be given\n";
		CREATE_EXCEPTION( utility::excn::Exception, oss.str() );
	}

	if ( tag->hasOption( "reference_pdb" ) ) {
		// load in a PDB
		std::string fname = tag->getOption< std::string >( "reference_pdb" );
		TR << "Loading reference Pose from PDB: " << fname << std::endl;
		core::pose::PoseOP pose = core::import_pose::pose_from_file( fname );
		reference_structure( pose );
	}

	if ( tag->hasOption( "reference_pose" ) ) {
		std::string posename = tag->getOption< std::string >( "reference_pose" );
		TR << "Loading reference Pose from resource pose named: " << posename << std::endl;
		reference_structure( datamap.get_ptr< core::pose::Pose >( "poses", posename ));
	}

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

	if ( !reference_pose_ ) {
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			TR << "Setting reference Pose from the command-line \"native\" flag" << std::endl;
			reference_structure( core::import_pose::pose_from_file(
				basic::options::option[ basic::options::OptionKeys::in::file::native ] ));
		} else {
			stringstream err_msg;
			err_msg
				<< "Attempting to run the Rotamer Recovery against a Reference Structure, "
				<< "but native has not been set.";
			utility_exit_with_message(err_msg.str());
		}
	}


	if ( pose.size() != reference_pose_->size() ) {
		// if waters are added to a pose, size may be longer than reference_pose
		TR << "Warning:: -- attempting to compare reference pose of size " << reference_pose_->size() << " to pose of size " << pose.size() << std::endl;
		//  stringstream err_msg;
		//  err_msg
		//   << "Attempting to run the Rotamer Recovery against Reference Structure protocol, "
		//   << "but the saved structure has a different number of residues.";
		//  utility_exit_with_message(err_msg.str());
	}

	for ( Size ii = 1; ii <= reference_pose_->size(); ++ii ) {
		if ( !packer_task.pack_residue(ii) || reference_pose_->residue(ii).is_water() ) continue;
		measure_rotamer_recovery(
			comparer, reporter,
			*reference_pose_, pose,
			reference_pose_->residue(ii), pose.residue(ii) );
	}
}

void RRProtocolReferenceStructure::append_attributes(
	utility::tag::AttributeList & attlist
){
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;
	attlist
		+ Attr( "reference_pdb", xs_string, "For use with the RRProtocolReferenceStructure. The PDB formatted file that should be compared against. Mutually exclusive with the 'reference_pose' attribute, but at least one of the two must be provided." )
		+ Attr( "reference_pose", xs_string, "For use with the RRProtocolReferenceStructure. The Pose held in the DataMap that should be compared against. This Pose should be loaded into the DataMap using the ResourceManager and the PoseFromPoseResourceMover. Mutually exclusive with the 'reference_pose' attribute, but at least one of the two must be provided." );
}


} // rotamer_recovery
} // protocols
