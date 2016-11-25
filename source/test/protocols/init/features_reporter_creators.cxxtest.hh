// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/StandardJobQueen.cxxtest.hh
/// @brief  test suite for the StandardJobQueen
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// basic headers
#include <basic/options/option.hh>

// Utility headers
#include <utility/string_util.hh>

// C++ headers
#include <sstream>

//  creator headers
#include <protocols/features/AtomInResidueAtomInResiduePairFeaturesCreator.hh>
#include <protocols/features/AtomAtomPairFeaturesCreator.hh>
#include <protocols/features/AtomTypesFeaturesCreator.hh>
#include <protocols/features/BetaTurnDetectionFeaturesCreator.hh>
#include <protocols/features/ChargeChargeFeaturesCreator.hh>
#include <protocols/features/DdGFeaturesCreator.hh>
#include <protocols/features/GeometricSolvationFeaturesCreator.hh>
#include <protocols/features/HBondFeaturesCreator.hh>
#include <protocols/features/HBondParameterFeaturesCreator.hh>
#include <protocols/features/HelixCapFeaturesCreator.hh>
#include <protocols/features/JobDataFeaturesCreator.hh>
#include <protocols/features/LoopAnchorFeaturesCreator.hh>
#include <protocols/features/ModelFeaturesCreator.hh>
#include <protocols/features/OrbitalsFeaturesCreator.hh>
#include <protocols/features/PairFeaturesCreator.hh>
#include <protocols/features/PdbDataFeaturesCreator.hh>
#include <protocols/features/PoseCommentsFeaturesCreator.hh>
#include <protocols/features/PoseConformationFeaturesCreator.hh>
#include <protocols/features/ProteinBackboneTorsionAngleFeaturesCreator.hh>
#include <protocols/features/ProteinBackboneAtomAtomPairFeaturesCreator.hh>
#include <protocols/features/ProteinBondGeometryFeaturesCreator.hh>
#include <protocols/features/ProteinResidueConformationFeaturesCreator.hh>
#include <protocols/features/ProteinRMSDFeaturesCreator.hh>
#include <protocols/features/ProteinRMSDNoSuperpositionFeaturesCreator.hh>
//#include <protocols/features/ProtocolFeaturesCreator.hh>
#include <protocols/features/RadiusOfGyrationFeaturesCreator.hh>
#include <protocols/features/ResidueBurialFeaturesCreator.hh>
#include <protocols/features/ResidueConformationFeaturesCreator.hh>
#include <protocols/features/ResidueFeaturesCreator.hh>
#include <protocols/features/ResidueScoresFeaturesCreator.hh>
#include <protocols/features/ResidueTotalScoresFeaturesCreator.hh>
#include <protocols/features/ResidueSecondaryStructureFeaturesCreator.hh>
#include <protocols/features/ResidueTypesFeaturesCreator.hh>
#include <protocols/features/RotamerFeaturesCreator.hh>
#include <protocols/features/RotamerBoltzmannWeightFeaturesCreator.hh>
#include <protocols/features/RotamerRecoveryFeaturesCreator.hh>
#include <protocols/features/RuntimeFeaturesCreator.hh>
#include <protocols/features/SaltBridgeFeaturesCreator.hh>
#include <protocols/features/ScoreTypeFeaturesCreator.hh>
#include <protocols/features/SecondaryStructureSegmentFeaturesCreator.hh>
#include <protocols/features/SmotifFeaturesCreator.hh>
#include <protocols/features/StructureFeaturesCreator.hh>
#include <protocols/features/StructureScoresFeaturesCreator.hh>
#include <protocols/features/ScreeningFeaturesCreator.hh>
#include <protocols/features/ScoreFunctionFeaturesCreator.hh>
#include <protocols/features/TaskOperationFeaturesCreator.hh>
#include <protocols/features/TotalScoreFeaturesCreator.hh>
#include <protocols/features/TrajectoryMapFeaturesCreator.hh>
#include <protocols/features/UnrecognizedAtomFeaturesCreator.hh>
//#include <protocols/features/BatchFeaturesCreator.hh>
#include <protocols/features/helixAssembly/HelixBundleFeaturesCreator.hh>
#include <protocols/features/helixAssembly/ConcurrencyTestCreator.hh>
#include <protocols/features/strand_assembly/SandwichFeaturesCreator.hh>
#include <protocols/features/strand_assembly/StrandBundleFeaturesCreator.hh>
//#include <protocols/features/strand_assembly/ReportAADistributionFromDBCreator.hh>
#include <protocols/features/ResidueGridScoresFeaturesCreator.hh>
#include <protocols/features/WaterFeaturesCreator.hh>
#include <protocols/features/InterfaceFeaturesCreator.hh>
#include <protocols/antibody/clusters/CDRClusterFeaturesCreator.hh>
#include <protocols/antibody/AntibodyFeaturesCreator.hh>


class BackwardsProtocolsFeaturesReporterCreatorTests : public CxxTest::TestSuite
{
public:

	void write_name_test( std::string const & creator_name, std::string const & features_reporter_name ) {
		std::cout << "{ " << creator_name << "  creator; TS_ASSERT_EQUALS( creator->type_name(), \"" << features_reporter_name << "\" ); } " << std::endl;
	}

	void dont_test_find_all_test()
	{
		{ protocols::features::AtomAtomPairFeaturesCreator creator; write_name_test( "protocols::features::AtomAtomPairFeaturesCreator", creator.type_name() );}
		{ protocols::features::AtomTypesFeaturesCreator creator; write_name_test( "protocols::features::AtomTypesFeaturesCreator", creator.type_name() );}
		{ protocols::features::AtomInResidueAtomInResiduePairFeaturesCreator creator; write_name_test( "protocols::features::AtomInResidueAtomInResiduePairFeaturesCreator", creator.type_name() );}
		{ protocols::features::BetaTurnDetectionFeaturesCreator creator; write_name_test( "protocols::features::BetaTurnDetectionFeaturesCreator", creator.type_name() );}
		{ protocols::features::ChargeChargeFeaturesCreator creator; write_name_test( "protocols::features::ChargeChargeFeaturesCreator", creator.type_name() );}
		{ protocols::features::DdGFeaturesCreator creator; write_name_test( "protocols::features::DdGFeaturesCreator", creator.type_name() );}
		{ protocols::features::GeometricSolvationFeaturesCreator creator; write_name_test( "protocols::features::GeometricSolvationFeaturesCreator", creator.type_name() );}
		{ protocols::features::HBondFeaturesCreator creator; write_name_test( "protocols::features::HBondFeaturesCreator", creator.type_name() );}
		{ protocols::features::HBondParameterFeaturesCreator creator; write_name_test( "protocols::features::HBondParameterFeaturesCreator", creator.type_name() );}
		{ protocols::features::HelixCapFeaturesCreator creator; write_name_test( "protocols::features::HelixCapFeaturesCreator", creator.type_name() );}
		{ protocols::features::JobDataFeaturesCreator creator; write_name_test( "protocols::features::JobDataFeaturesCreator", creator.type_name() );}
		{ protocols::features::LoopAnchorFeaturesCreator creator; write_name_test( "protocols::features::LoopAnchorFeaturesCreator", creator.type_name() );}
		{ protocols::features::ModelFeaturesCreator creator; write_name_test( "protocols::features::ModelFeaturesCreator", creator.type_name() );}
		{ protocols::features::OrbitalsFeaturesCreator creator; write_name_test( "protocols::features::OrbitalsFeaturesCreator", creator.type_name() );}
		{ protocols::features::PairFeaturesCreator creator; write_name_test( "protocols::features::PairFeaturesCreator", creator.type_name() );}
		{ protocols::features::PdbDataFeaturesCreator creator; write_name_test( "protocols::features::PdbDataFeaturesCreator", creator.type_name() );}
		{ protocols::features::PoseCommentsFeaturesCreator creator; write_name_test( "protocols::features::PoseCommentsFeaturesCreator", creator.type_name() );}
		{ protocols::features::PoseConformationFeaturesCreator creator; write_name_test( "protocols::features::PoseConformationFeaturesCreator", creator.type_name() );}
		{ protocols::features::ProteinBackboneAtomAtomPairFeaturesCreator creator; write_name_test( "protocols::features::ProteinBackboneAtomAtomPairFeaturesCreator", creator.type_name() );}
		{ protocols::features::ProteinBackboneTorsionAngleFeaturesCreator creator; write_name_test( "protocols::features::ProteinBackboneTorsionAngleFeaturesCreator", creator.type_name() );}
		{ protocols::features::ProteinBondGeometryFeaturesCreator creator; write_name_test( "protocols::features::ProteinBondGeometryFeaturesCreator", creator.type_name() );}
		{ protocols::features::ProteinResidueConformationFeaturesCreator creator; write_name_test( "protocols::features::ProteinResidueConformationFeaturesCreator", creator.type_name() );}
		{ protocols::features::ProteinRMSDFeaturesCreator creator; write_name_test( "protocols::features::ProteinRMSDFeaturesCreator", creator.type_name() );}
		{ protocols::features::ProteinRMSDNoSuperpositionFeaturesCreator creator; write_name_test( "protocols::features::ProteinRMSDNoSuperpositionFeaturesCreator", creator.type_name() );}
		{ protocols::features::RadiusOfGyrationFeaturesCreator creator; write_name_test( "protocols::features::RadiusOfGyrationFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueBurialFeaturesCreator creator; write_name_test( "protocols::features::ResidueBurialFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueConformationFeaturesCreator creator; write_name_test( "protocols::features::ResidueConformationFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueFeaturesCreator creator; write_name_test( "protocols::features::ResidueFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueScoresFeaturesCreator creator; write_name_test( "protocols::features::ResidueScoresFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueSecondaryStructureFeaturesCreator creator; write_name_test( "protocols::features::ResidueSecondaryStructureFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueTotalScoresFeaturesCreator creator; write_name_test( "protocols::features::ResidueTotalScoresFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueTypesFeaturesCreator creator; write_name_test( "protocols::features::ResidueTypesFeaturesCreator", creator.type_name() );}
		{ protocols::features::RotamerBoltzmannWeightFeaturesCreator creator; write_name_test( "protocols::features::RotamerBoltzmannWeightFeaturesCreator", creator.type_name() );}
		{ protocols::features::RotamerFeaturesCreator creator; write_name_test( "protocols::features::RotamerFeaturesCreator", creator.type_name() );}
		{ protocols::features::RotamerRecoveryFeaturesCreator creator; write_name_test( "protocols::features::RotamerRecoveryFeaturesCreator", creator.type_name() );}
		{ protocols::features::RuntimeFeaturesCreator creator; write_name_test( "protocols::features::RuntimeFeaturesCreator", creator.type_name() );}
		{ protocols::features::SaltBridgeFeaturesCreator creator; write_name_test( "protocols::features::SaltBridgeFeaturesCreator", creator.type_name() );}
		{ protocols::features::ScreeningFeaturesCreator creator; write_name_test( "protocols::features::ScreeningFeaturesCreator", creator.type_name() );}
		{ protocols::features::ScoreTypeFeaturesCreator creator; write_name_test( "protocols::features::ScoreTypeFeaturesCreator", creator.type_name() );}
		{ protocols::features::SecondaryStructureSegmentFeaturesCreator creator; write_name_test( "protocols::features::SecondaryStructureSegmentFeaturesCreator", creator.type_name() );}
		{ protocols::features::SmotifFeaturesCreator creator; write_name_test( "protocols::features::SmotifFeaturesCreator", creator.type_name() );}
		{ protocols::features::StructureFeaturesCreator creator; write_name_test( "protocols::features::StructureFeaturesCreator", creator.type_name() );}
		{ protocols::features::StructureScoresFeaturesCreator creator; write_name_test( "protocols::features::StructureScoresFeaturesCreator", creator.type_name() );}
		{ protocols::features::ScoreFunctionFeaturesCreator creator; write_name_test( "protocols::features::ScoreFunctionFeaturesCreator", creator.type_name() );}
		{ protocols::features::TaskOperationFeaturesCreator creator; write_name_test( "protocols::features::TaskOperationFeaturesCreator", creator.type_name() );}
		{ protocols::features::TotalScoreFeaturesCreator creator; write_name_test( "protocols::features::TotalScoreFeaturesCreator", creator.type_name() );}
		{ protocols::features::TrajectoryMapFeaturesCreator creator; write_name_test( "protocols::features::TrajectoryMapFeaturesCreator", creator.type_name() );}
		{ protocols::features::UnrecognizedAtomFeaturesCreator creator; write_name_test( "protocols::features::UnrecognizedAtomFeaturesCreator", creator.type_name() );}
		{ protocols::features::helixAssembly::HelixBundleFeaturesCreator creator; write_name_test( "protocols::features::helixAssembly::HelixBundleFeaturesCreator", creator.type_name() );}
		{ protocols::features::helixAssembly::ConcurrencyTestCreator creator; write_name_test( "protocols::features::helixAssembly::ConcurrencyTestCreator", creator.type_name() );}
		{ protocols::features::strand_assembly::SandwichFeaturesCreator creator; write_name_test( "protocols::features::strand_assembly::SandwichFeaturesCreator", creator.type_name() );}
		{ protocols::features::strand_assembly::StrandBundleFeaturesCreator creator; write_name_test( "protocols::features::strand_assembly::StrandBundleFeaturesCreator", creator.type_name() );}
		{ protocols::features::ResidueGridScoresFeaturesCreator creator; write_name_test( "protocols::features::ResidueGridScoresFeaturesCreator", creator.type_name() );}
		{ protocols::features::WaterFeaturesCreator creator; write_name_test( "protocols::features::WaterFeaturesCreator", creator.type_name() );}
		{ protocols::features::InterfaceFeaturesCreator creator; write_name_test( "protocols::features::InterfaceFeaturesCreator", creator.type_name() );}
		{ protocols::antibody::clusters::CDRClusterFeaturesCreator creator; write_name_test( "protocols::antibody::clusters::CDRClusterFeaturesCreator", creator.type_name() );}
		{ protocols::antibody::AntibodyFeaturesCreator creator; write_name_test( "protocols::antibody::AntibodyFeaturesCreator", creator.type_name() );}
	}
	void test_all_creators_generate_correct_name_strings()
	{
		{ protocols::features::AtomAtomPairFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "AtomAtomPairFeatures" ); }
		{ protocols::features::AtomTypesFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "AtomTypesFeatures" ); }
		{ protocols::features::AtomInResidueAtomInResiduePairFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "AtomInResidueAtomInResiduePairFeatures" ); }
		{ protocols::features::BetaTurnDetectionFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "BetaTurnDetectionFeatures" ); }
		{ protocols::features::ChargeChargeFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ChargeChargeFeatures" ); }
		{ protocols::features::DdGFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "DdGFeatures" ); }
		{ protocols::features::GeometricSolvationFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "GeometricSolvationFeatures" ); }
		{ protocols::features::HBondFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "HBondFeatures" ); }
		{ protocols::features::HBondParameterFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "HBondParameterFeatures" ); }
		{ protocols::features::HelixCapFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "HelixCapFeatures" ); }
		{ protocols::features::JobDataFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "JobDataFeatures" ); }
		{ protocols::features::LoopAnchorFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "LoopAnchorFeatures" ); }
		{ protocols::features::ModelFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ModelFeatures" ); }
		{ protocols::features::OrbitalsFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "OrbitalsFeatures" ); }
		{ protocols::features::PairFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "PairFeatures" ); }
		{ protocols::features::PdbDataFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "PdbDataFeatures" ); }
		{ protocols::features::PoseCommentsFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "PoseCommentsFeatures" ); }
		{ protocols::features::PoseConformationFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "PoseConformationFeatures" ); }
		{ protocols::features::ProteinBackboneAtomAtomPairFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ProteinBackboneAtomAtomPairFeatures" ); }
		{ protocols::features::ProteinBackboneTorsionAngleFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ProteinBackboneTorsionAngleFeatures" ); }
		{ protocols::features::ProteinBondGeometryFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ProteinBondGeometryFeatures" ); }
		{ protocols::features::ProteinResidueConformationFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ProteinResidueConformationFeatures" ); }
		{ protocols::features::ProteinRMSDFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ProteinRMSDFeatures" ); }
		{ protocols::features::ProteinRMSDNoSuperpositionFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ProteinRMSDNoSuperpositionFeatures" ); }
		{ protocols::features::RadiusOfGyrationFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "RadiusOfGyrationFeatures" ); }
		{ protocols::features::ResidueBurialFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueBurialFeatures" ); }
		{ protocols::features::ResidueConformationFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueConformationFeatures" ); }
		{ protocols::features::ResidueFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueFeatures" ); }
		{ protocols::features::ResidueScoresFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueScoresFeatures" ); }
		{ protocols::features::ResidueSecondaryStructureFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueSecondaryStructureFeatures" ); }
		{ protocols::features::ResidueTotalScoresFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueTotalScoresFeatures" ); }
		{ protocols::features::ResidueTypesFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueTypesFeatures" ); }
		{ protocols::features::RotamerBoltzmannWeightFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "RotamerBoltzmannWeightFeatures" ); }
		{ protocols::features::RotamerFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "RotamerFeatures" ); }
		{ protocols::features::RotamerRecoveryFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "RotamerRecoveryFeatures" ); }
		{ protocols::features::RuntimeFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "RuntimeFeatures" ); }
		{ protocols::features::SaltBridgeFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "SaltBridgeFeatures" ); }
		{ protocols::features::ScreeningFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ScreeningFeatures" ); }
		{ protocols::features::ScoreTypeFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ScoreTypeFeatures" ); }
		{ protocols::features::SecondaryStructureSegmentFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "SecondaryStructureSegmentFeatures" ); }
		{ protocols::features::SmotifFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "SmotifFeatures" ); }
		{ protocols::features::StructureFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "StructureFeatures" ); }
		{ protocols::features::StructureScoresFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "StructureScoresFeatures" ); }
		{ protocols::features::ScoreFunctionFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ScoreFunctionFeatures" ); }
		{ protocols::features::TaskOperationFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "TaskOperationFeatures" ); }
		{ protocols::features::TotalScoreFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "TotalScoreFeatures" ); }
		{ protocols::features::TrajectoryMapFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "TrajectoryMapFeatures" ); }
		{ protocols::features::UnrecognizedAtomFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "UnrecognizedAtomFeatures" ); }
		{ protocols::features::helixAssembly::HelixBundleFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "HelixBundleFeatures" ); }
		{ protocols::features::helixAssembly::ConcurrencyTestCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ConcurrencyTest" ); }
		{ protocols::features::strand_assembly::SandwichFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "SandwichFeatures" ); }
		{ protocols::features::strand_assembly::StrandBundleFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "StrandBundleFeatures" ); }
		{ protocols::features::ResidueGridScoresFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "ResidueGridScoresFeatures" ); }
		{ protocols::features::WaterFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "WaterFeatures" ); }
		{ protocols::features::InterfaceFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "InterfaceFeatures" ); }
		{ protocols::antibody::clusters::CDRClusterFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "CDRClusterFeatures" ); }
		{ protocols::antibody::AntibodyFeaturesCreator creator; TS_ASSERT_EQUALS( creator.type_name(), "AntibodyFeatures" ); }

	}

};
