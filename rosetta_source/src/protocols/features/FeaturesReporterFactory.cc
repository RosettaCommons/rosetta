// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/FeaturesReporterFactory.cc
/// @brief  Create FeaturesReporters
/// @author Matthew O'Meara (mattjomeara@gmail.com)

/// @detail Right now this is a big switch statement, but it might be
/// worth refactoring this into a load-time factory registration
/// system like the score terms, movers etc.

// Unit Headers
#include <protocols/features/FeaturesReporterFactory.hh>

// Project Headers
#include <protocols/features/AtomAtomPairFeatures.hh>
#include <protocols/features/GeometricSolvationFeatures.hh>
#include <protocols/features/HBondFeatures.hh>
#include <protocols/features/HBondParameterFeatures.hh>
#include <protocols/features/OrbitalsFeatures.hh>
#include <protocols/features/PairFeatures.hh>
#include <protocols/features/PoseCommentsFeatures.hh>
#include <protocols/features/PoseConformationFeatures.hh>
#include <protocols/features/ProteinBackboneTorsionAngleFeatures.hh>
#include <protocols/features/ProteinBackboneAtomAtomPairFeatures.hh>
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/RadiusOfGyrationFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/ResidueTypesFeatures.hh>
#include <protocols/features/RotamerRecoveryFeatures.hh>
#include <protocols/features/RotamerBoltzmannWeightFeatures.hh>
#include <protocols/features/ResidueBurialFeatures.hh>
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/StructureScoresFeatures.hh>

// C++ Headers
#include <sstream>

namespace protocols {
namespace features {

	using std::endl;
	using std::string;
	using std::stringstream;
	using core::scoring::ScoreFunctionOP;

	FeaturesReporterFactory::FeaturesReporterFactory() {}

	FeaturesReporterFactory::FeaturesReporterFactory(
		const FeaturesReporterFactory &
	) {}

	FeaturesReporterFactory::~FeaturesReporterFactory() {}

	FeaturesReporterOP
	FeaturesReporterFactory::create_features_reporter(
		string const & name,
		ScoreFunctionOP scfxn) {

		if(name == "AtomAtomPairFeatures"){
			return new AtomAtomPairFeatures();
		} else if(name == "GeometricSolvationFeatures") {
			return new GeometricSolvationFeatures();
		} else if(name == "HBondFeatures") {
			return new HBondFeatures(scfxn);
		} else if(name == "HBondParameterFeatures") {
			return new HBondParameterFeatures();
		} else if (name == "OrbitalsFeatures") {
			return new OrbitalsFeatures();
		} else if (name == "PairFeatures") {
			return new PairFeatures();
		} else if (name == "PoseCommentsFeatures") {
			return new PoseCommentsFeatures();
		} else if (name == "PoseConformationFeatures") {
			return new PoseConformationFeatures();
		} else if (name == "ProteinBackboneTorsionAngleFeatures") {
			return new ProteinBackboneTorsionAngleFeatures();
		} else if (name == "ProteinBackboneAtomAtomPairFeatures") {
			return new ProteinBackboneAtomAtomPairFeatures();
		} else if (name == "ProteinResidueConformationFeatures") {
			return new ProteinResidueConformationFeatures();
		} else if (name == "ProtocolFeatures") {
			return new ProtocolFeatures();
		} else if (name == "RadiusOfGyrationFeatures") {
			return new RadiusOfGyrationFeatures();
		} else if (name == "ResidueFeatures") {
			return new ResidueFeatures();
		} else if (name == "ResidueTypesFeatures") {
			return new ResidueTypesFeatures();
		} else if (name == "RotamerBoltzmannWeightFeatures") {
			return new RotamerBoltzmannWeightFeatures(scfxn);
		} else if (name == "RotamerRecoveryFeatures") {
			return new RotamerRecoveryFeatures(scfxn);
		} else if (name == "ResidueBurialFeatures") {
			return new ResidueBurialFeatures();
		} else if (name == "ResidueSecondaryStructureFeatures") {
			return new ResidueSecondaryStructureFeatures();
		} else if (name == "StructureFeatures") {
			return new StructureFeatures();
		} else if (name == "StructureScoresFeatures") {
			return new StructureScoresFeatures(scfxn);
		} else {
			stringstream error_message;
			error_message << "Attempting to create unrecognized FeaturesReporter '" << name << "'." << endl;
			utility_exit_with_message(error_message.str());
		}
		return NULL;
	}

} // namespace
} // namespace
