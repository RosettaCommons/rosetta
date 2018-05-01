// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author  Ingemar Andre

// Unit headers
#include <protocols/minimization_packing/symmetry/SymRotamerTrialsMoverCreator.hh>
#include <protocols/minimization_packing/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <utility>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace minimization_packing {
namespace symmetry {

SymRotamerTrialsMover::~SymRotamerTrialsMover() = default;

std::string SymRotamerTrialsMover::get_name() const {
	return mover_name();
}

std::string SymRotamerTrialsMover::mover_name() {
	return "SymRotamerTrialsMover";
}

void SymRotamerTrialsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_rotamer_trials_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "This mover goes through each repackable/redesignable position in the pose, taking every permitted rotamer in turn, and evaluating the energy. Each position is then updated to the lowest energy rotamer. It does not consider coordinated changes at multiple residues, and may need several invocations to reach convergence.\n"
		"NOTE: This Mover is provided for historical support only. The regular RotamerTrialsMover should handle symmetry transparently and is preferred." )
		.write_complex_type_to_schema( xsd );

	//SymRotamersTrial description: The symmetric versions of pack rotamers and rotamer trials movers (they take the same tags as asymmetric versions)

}

std::string SymRotamerTrialsMoverCreator::keyname() const {
	return SymRotamerTrialsMover::mover_name();
}

protocols::moves::MoverOP
SymRotamerTrialsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SymRotamerTrialsMover );
}

void SymRotamerTrialsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymRotamerTrialsMover::provide_xml_schema( xsd );
}


/////////////////////////
// default constructor

SymEnergyCutRotamerTrialsMover::~SymEnergyCutRotamerTrialsMover() = default;

std::string
SymEnergyCutRotamerTrialsMover::get_name() const {
	return "SymEnergyCutRotamerTrialsMover";
}


} // symmetry
} // moves
} // protocols
