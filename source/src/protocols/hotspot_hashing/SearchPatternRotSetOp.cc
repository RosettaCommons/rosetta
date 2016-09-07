// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/hotspot_hashing/SearchPatternRotSetOp.cc
/// @brief  Creates rigid body variants from search pattern during repacking.
/// @author Alex Ford (fordas@uw.edu)

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/RigidBodyMoveRotSetOps.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <basic/Tracer.hh>

//Project headers
#include <utility>
#include <utility/exit.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/hotspot_hashing/StubGenerator.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>
#include <protocols/hotspot_hashing/SearchPatternRotSetOp.hh>

#include <utility/vector1.hh>
#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

namespace protocols {
namespace hotspot_hashing {

static THREAD_LOCAL basic::Tracer TR( "protocols.hotspot_hashing.SearchPatternRotSetOp" );

SearchPatternRotSetOp::SearchPatternRotSetOp( SearchPatternOP pattern ) :
	protocols::toolbox::rotamer_set_operations::RigidBodyMoveBaseRSO(),
	search_stubs_(pattern->Searchpoints())
{}

SearchPatternRotSetOp::SearchPatternRotSetOp( SearchPatternRotSetOp const & other ) :
	protocols::toolbox::rotamer_set_operations::RigidBodyMoveBaseRSO(),
	search_stubs_(other.search_stubs_)
{}

core::pack::rotamer_set::RotamerSetOperationOP
SearchPatternRotSetOp::clone() const
{
	return core::pack::rotamer_set::RotamerSetOperationOP( new SearchPatternRotSetOp(*this) );
}


utility::vector1< core::conformation::ResidueCOP > SearchPatternRotSetOp::get_rigid_body_confs(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & /*ptask*/,
	core::Size residue_index)
{
	utility::vector1< core::conformation::ResidueCOP > result_residues;


	// Create residue copy in centroid frame
	core::conformation::ResidueOP source_residue( new core::conformation::Residue(pose.residue(residue_index)) );
	TR.Debug << "Setting rigid body confs for residue: " << residue_index << " " << source_residue->name() << std::endl;
	core::kinematics::Stub source_frame = StubGenerator::residueStubCentroidFrame(source_residue);
	TR.Trace << "Generated target residue centroid frame: ";
	stub_to_points(TR.Trace, source_frame) << std::endl;

	StubGenerator::moveFromStubFrame(source_residue, source_frame);

	for ( core::Size i = 1; i <= search_stubs_.size(); ++i ) {
		// Create residue copies translating the search frame into the centroid frame

		core::conformation::ResidueOP search_residue( new core::conformation::Residue(*source_residue) );

		StubGenerator::moveIntoStubFrame(search_residue, search_stubs_[i]);
		StubGenerator::moveIntoStubFrame(search_residue, source_frame);

		if ( TR.Trace.visible() ) {
			core::kinematics::Stub search_frame = StubGenerator::residueStubCentroidFrame(search_residue);
			TR.Trace << "Generated alt rb residue: " << search_residue->name() << std::endl;
			TR.Trace << "Generated alt rb residue residue centroid frame: ";
			stub_to_points(TR.Trace, search_frame) << std::endl;
		}

		result_residues.push_back(core::conformation::ResidueCOP(search_residue));
	}

	TR.Debug << "Rigid body confs count: " << result_residues.size() << std::endl;

	return result_residues;
}

core::Real SearchPatternRotSetOp::increase_packer_residue_radius(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskCOP,
	core::Size residue_index
)
{
	// Calculate the centroid frame and centroid location in the centroid frame
	core::conformation::ResidueCOP source_residue( core::conformation::ResidueOP( new core::conformation::Residue(pose.residue(residue_index)) ) );
	TR.Debug << "Increasing packer residue radius for residue: " << residue_index << " " << source_residue->name() << std::endl;

	core::kinematics::Stub centroid_frame = StubGenerator::residueStubCentroidFrame(source_residue);
	Vector centroid_nbr_location = centroid_frame.global2local(source_residue->xyz(source_residue->nbr_atom()));

	// Calculate centroid frame coords of nbr atom location translated into all search frames
	core::Real max_sq_dist(0.0);
	for ( core::Size i = 1; i <= search_stubs_.size(); ++i ) {
		core::Vector search_nbr_location = search_stubs_[i].local2global(centroid_nbr_location);
		core::Real sq_dist = (centroid_nbr_location - search_nbr_location).length_squared();

		if ( sq_dist > max_sq_dist ) {
			max_sq_dist = sq_dist;
		}
	}

	core::Real max_distance = std::sqrt( max_sq_dist );

	// Return largest nbr_atom location displacement
	TR.Debug << "Increased packer residue radius size: " << max_distance << std::endl;
	return max_distance;
}

AddSearchPatternRotSetOp::AddSearchPatternRotSetOp(
	core::Size target_residue,
	SearchPatternOP pattern) :
	target_residue_(target_residue),
	pattern_(std::move(pattern))
{}

AddSearchPatternRotSetOp::~AddSearchPatternRotSetOp()= default;

core::pack::task::operation::TaskOperationOP
AddSearchPatternRotSetOp::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new AddSearchPatternRotSetOp(*this) );
}

/// @details not valid now, extend with search pattern parsing
void
AddSearchPatternRotSetOp::parse_tag( TagCOP /*tag*/ , DataMap & )
{
	utility_exit_with_message("Can not parse tag for add search pattern rot set op.");
}

void
AddSearchPatternRotSetOp::apply(
	core::pose::Pose const & /*pose*/,
	core::pack::task::PackerTask & task
) const
{
	if ( !task.being_packed( target_residue_) ) {
		return;
	}

	SearchPatternRotSetOpOP rotsetop( new SearchPatternRotSetOp(pattern_) );
	task.nonconst_residue_task( target_residue_ ).append_rotamerset_operation( rotsetop );
}

}
}
