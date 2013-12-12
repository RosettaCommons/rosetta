// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol headers
#include <protocols/features/ReportToDB.hh>
#include <protocols/features/PoseConformationFeatures.hh>

// Utility headers
#include <devel/init.hh>

#include <protocols/features/ReportToDB.hh>
#include <protocols/features/FeaturesReporter.hh>

#include <protocols/features/ScoreTypeFeatures.hh>
#include <protocols/features/StructureScoresFeatures.hh>
#include <protocols/features/PdbDataFeatures.hh>
#include <protocols/features/PoseCommentsFeatures.hh>
#include <protocols/features/PoseConformationFeatures.hh>
#include <protocols/features/ProteinResidueConformationFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/ResidueConformationFeatures.hh>
#include <protocols/features/JobDataFeatures.hh>
#include <protocols/features/DatabaseFilters.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/util.hh>

using namespace std;
using core::pose::Pose;
using core::import_pose::pose_from_pdb;
using protocols::features::ReportToDB;
using protocols::features::ReportToDBOP;
using protocols::features::PoseConformationFeatures;

// This program is just supposed to read in a pdb and spit it back out again.

int main(int argc, char** argv) {

	devel::init(argc, argv);

	Pose pose;
	pose_from_pdb(pose, "structures/linear/4.symmetry.pdb");

	//ReportToDBOP reporter = new ReportToDB;
	//reporter->add_features_reporter(new PoseConformationFeatures);

	protocols::features::ReportToDBOP reporter = new protocols::features::ReportToDB();
	//reporter_->set_relevant_residues_mode( protocols::features::RelevantResiduesMode::Implicit );
	reporter->set_batch_description( "Rosetta: FragMixFeatures" );

	// TODO: better checking for batch names. Check and see if db_job_outputter exists already?
	reporter->set_batch_name( "FragMix" );

	reporter->add_features_reporter( new protocols::features::PdbDataFeatures() );
	reporter->add_features_reporter( new protocols::features::ScoreTypeFeatures() );
	reporter->add_features_reporter( new protocols::features::PoseConformationFeatures() );
	reporter->add_features_reporter( new protocols::features::PoseCommentsFeatures() );
	reporter->add_features_reporter( new protocols::features::ResidueFeatures() );
	reporter->add_features_reporter( new protocols::features::ProteinResidueConformationFeatures() );
	reporter->add_features_reporter( new protocols::features::ResidueConformationFeatures() );
	reporter->add_features_reporter( new protocols::features::JobDataFeatures() );

	reporter->apply(pose);
}

