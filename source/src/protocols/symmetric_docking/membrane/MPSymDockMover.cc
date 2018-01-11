// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/symmetric_docking/membrane/MPSymDockMover.cc
///
/// @brief      Membrane symmetric docking protocol
/// @details    Given an asymmetric starting pose, setup the pose for symmetry,
///             add a membrane representation for the full symmetric complex,
///             position the pose at the center of mass of all transmembrane spans,
///             and proceed with the standard symmetric docking protocol.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 2/9/15

// Unit Headers
#include <protocols/symmetric_docking/membrane/MPSymDockMover.hh>
#include <protocols/symmetric_docking/membrane/MPSymDockMoverCreator.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/symmetry/SymmetricAddMembraneMover.hh>
#include <protocols/membrane/symmetry/util.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>

#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <core/pose/symmetry/util.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.symmetric_docking.membrane.MPSymDockMover" );

namespace protocols {
namespace symmetric_docking {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Defualt constructor for the membrane protein symmetric
/// docking protocol
MPSymDockMover::MPSymDockMover() :
	protocols::moves::Mover() {}

/// @brief Destructor
MPSymDockMover::~MPSymDockMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPSymDockMover::clone() const {
	return ( protocols::moves::MoverOP( new MPSymDockMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPSymDockMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPSymDockMover );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MPSymDockMover::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
}

/// @brief Create a new copy of this mover
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MPSymDockMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MPSymDockMover );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP MPSymDockMoverCreator::keyname() const {
// XRW TEMP  return MPSymDockMover::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP MPSymDockMover::mover_name() {
// XRW TEMP  return "MPSymDockMover";
// XRW TEMP }

//////////////////////
/// Mover Methods  ///
//////////////////////

/// @brief Return the name of this mover (MPSymDockMover)
// XRW TEMP std::string
// XRW TEMP MPSymDockMover::get_name() const {
// XRW TEMP  return "MPSymDockMover";
// XRW TEMP }

/// @brief Apply Method: Symmetric Docking in the membrane
/// @details Setup pose for symmetry, add the membrane components to the total pose,
/// and then perform symmetric docking
void
MPSymDockMover::apply( Pose & pose ) {

	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace core::scoring;
	using namespace protocols::simple_moves::symmetry;
	using namespace protocols::membrane::symmetry;
	using namespace protocols::symmetric_docking;

	// Check that the pose is neither symmetric nor a membrane protein
	// Protocol is not that advanced yet - build a clean new conf every time
	if ( is_symmetric( pose ) || pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot setup a new symmetric membrane protein if the conformation is already symmetric and/or a membrane pose" );
	}

	// Setup the pose for symmetry based on inputs
	SetupForSymmetryMoverOP setup_for_symm = SetupForSymmetryMoverOP( new SetupForSymmetryMover() );
	setup_for_symm->apply( pose );

	SymmetricAddMembraneMoverOP add_symm_memb = SymmetricAddMembraneMoverOP( new SymmetricAddMembraneMover() );
	add_symm_memb->apply( pose );

	// Position the membrane at the center of mass of transmembrane spans
	// of the full symmetric complex
	position_membrane_at_topology_com( pose );

	// Configure Symmetric MP Docking. Rosetta will symmetrize the
	// score functions automatically
	ScoreFunctionOP sfxn_high = ScoreFunctionFactory::create_score_function( "mpframework_symdock_fa_2015" );
	ScoreFunctionOP sfxn_low = ScoreFunctionFactory::create_score_function( "mpframework_symdock_cen_2015" );

	// Setup repulsives based slide criteria (better for membranes than
	// initial contact scoring)
	auto & symm_conf ( dynamic_cast< SymmetricConformation & > ( pose.conformation()) );
	symm_conf.Symmetry_Info()->get_slide_info().set_SlideCriteriaType( FA_REP_SCORE );

	// Setup a docking protocol, don't apply filters during docking runs
	// for now
	SymDockProtocolOP symdock( new SymDockProtocol( true, false, false, sfxn_low, sfxn_high ) );
	symdock->apply( pose );

	// Done!
}

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Position Pose at TopologyCOM
/// @details Position the initial pose at the center of mass of
/// predicted transmembrane spanning regions of the whole complex
void
MPSymDockMover::position_membrane_at_topology_com( core::pose::Pose & pose ) {

	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
	using namespace protocols::membrane::symmetry;

	// Symmetrize the spanning topology
	SpanningTopologyOP symmetrized_topology = symmetrize_spans( pose, *pose.membrane_info()->spanning_topology() );
	symmetrized_topology->show();
	Embedding embeddings = Embedding( *symmetrized_topology, pose );
	EmbeddingDefOP final_embed = embeddings.total_embed();

	// Set Initial Membrane position to be the symmetrized position
	SetMembranePositionMoverOP set_initial_position( new SetMembranePositionMover( final_embed->center(), final_embed->normal() ) );
	set_initial_position->apply( pose );

}

std::string MPSymDockMover::get_name() const {
	return mover_name();
}

std::string MPSymDockMover::mover_name() {
	return "MPSymDockMover";
}

void MPSymDockMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// Empty parse_my_tag so this is intentional
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Given an asymmetric starting pose, setup the pose for symmetry, add a membrane representation for the full symmetric complex, position the pose at the center of mass of all transmembrane spans, and proceed with the standard symmetric docking protocol.", attlist );
}

std::string MPSymDockMoverCreator::keyname() const {
	return MPSymDockMover::mover_name();
}

protocols::moves::MoverOP
MPSymDockMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MPSymDockMover );
}

void MPSymDockMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MPSymDockMover::provide_xml_schema( xsd );
}


} // membrane
} // symmetric_docking
} // protocols

