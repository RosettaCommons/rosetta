// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RotamerRecoveryFeatures.cc
/// @brief  report rotamer recover features and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/RotamerRecoveryFeatures.hh>

// Package Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>
#include <protocols/rotamer_recovery/RRReporterSQLite.hh>
#include <protocols/rotamer_recovery/RRProtocolMover.hh>


// Project Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <sstream>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>

//Auto Headers
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>



namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using basic::Tracer;
using core::Size;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunction;
using core::scoring::getScoreFunction;
using core::pose::Pose;
using core::pack::task::PackerTaskOP;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::InitializeFromCommandline;
using core::pack::task::operation::RestrictToRepacking;
using protocols::filters::Filters_map;
using protocols::moves::DataMap;
using protocols::moves::MoverOP;
using protocols::moves::Movers_map;
using protocols::rosetta_scripts::parse_mover;
using protocols::rotamer_recovery::RotamerRecovery;
using protocols::rotamer_recovery::RotamerRecoveryFactory;
using protocols::rotamer_recovery::RRProtocolMover;
using protocols::rotamer_recovery::RRReporterSQLite;
using protocols::rotamer_recovery::RRReporterSQLiteOP;
using utility::sql_database::sessionOP;
using utility::tag::TagPtr;
using utility::vector1;
using cppdb::statement;

static Tracer TR("protocols.features.RotamerRecoveryFeatures");

RotamerRecoveryFeatures::RotamerRecoveryFeatures() :
	scfxn_(getScoreFunction()),
	protocol_(),
	comparer_()
{}

RotamerRecoveryFeatures::RotamerRecoveryFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(scfxn),
	protocol_(),
	comparer_()
{}

RotamerRecoveryFeatures::RotamerRecoveryFeatures(
	RotamerRecoveryFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_),
	protocol_(),
	comparer_()
{}

RotamerRecoveryFeatures::~RotamerRecoveryFeatures() {}

string
RotamerRecoveryFeatures::type_name() const { return "RotamerRecoveryFeatures"; }

string
RotamerRecoveryFeatures::schema() const {
	return RRReporterSQLite::schema(
		RRReporterSQLite::OutputLevel::features );
}

utility::vector1<std::string>
RotamerRecoveryFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}


void
RotamerRecoveryFeatures::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & movers,
	Pose const & /*pose*/
) {
	if(tag->hasOption("scorefxn")){
		string scorefxn_name = tag->getOption<string>("scorefxn");
		scfxn_ = data.get<ScoreFunction*>("scorefxns", scorefxn_name);
	} else {
		stringstream error_msg;
		error_msg
			<< "The " << type_name() << " reporter requires a 'scorefxn' tag:" << endl
			<< endl
			<< "    <feature name=" << type_name() <<" scorefxn=(name_of_score_function) />" << endl;
		utility_exit_with_message(error_msg.str());
	}

	RotamerRecoveryFactory * factory(RotamerRecoveryFactory::get_instance());

	if(tag->hasOption("mover") || tag->hasOption("mover_name")){
		MoverOP mover = parse_mover(tag->hasOption("mover") ?
			tag->getOption<string>("mover") : tag->getOption<string>("mover_name"), movers);
		protocol_ = new RRProtocolMover(mover);
	} else {
		string const & protocol_name(tag->getOption<string>("protocol", "RRProtocolMinPack"));
		protocol_ = factory->get_rotamer_recovery_protocol(protocol_name);
	}

	string const & comparer_name(tag->getOption<string>("comparer", "RRComparerAutomorphicRMSD"));
	comparer_ = factory->get_rotamer_recovery_comparer(comparer_name);
}



Size
RotamerRecoveryFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	RRReporterSQLiteOP reporter(
		new RRReporterSQLite(db_session, RRReporterSQLite::OutputLevel::features));

	reporter->set_struct_id1(struct_id);
	RotamerRecovery rotamer_recovery(protocol_, comparer_, reporter);

	// I would like to assert that this has been called but I don't know how.
	//scfxn.setup_for_scoring(pose);

	TaskFactory task_factory;
	PackerTaskOP packer_task(task_factory.create_packer_task(pose));
	packer_task->restrict_to_repacking();
	packer_task->restrict_to_residues(relevant_residues);

	rotamer_recovery.run(pose, *scfxn_, *packer_task);

	return 0;
}


} // namesapce
} // namespace
