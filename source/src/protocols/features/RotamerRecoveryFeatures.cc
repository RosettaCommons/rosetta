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
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/features/ReportToDB.hh>
#include <protocols/rotamer_recovery/RotamerRecovery.hh>
#include <protocols/rotamer_recovery/RotamerRecoveryFactory.hh>
#include <protocols/rotamer_recovery/RRReporterSQLite.hh>
#include <protocols/rotamer_recovery/RRProtocolMover.hh>
#include <protocols/rotamer_recovery/RRProtocolReferenceStructure.hh>
#include <protocols/rotamer_recovery/RRProtocolRTMin.hh>
#include <protocols/rotamer_recovery/RRProtocolRelax.hh>


// Project Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <string>
#include <sstream>

#include <utility/vector0.hh>

//Auto Headers
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <utility/excn/Exceptions.hh>
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
using core::scoring::get_score_function;
using core::pose::Pose;
using core::pose::PoseCOP;
using core::pack::task::PackerTaskOP;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::InitializeFromCommandline;
using core::pack::task::operation::RestrictToRepacking;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::MoverOP;
using protocols::moves::Movers_map;
using protocols::rosetta_scripts::parse_mover;
using protocols::rosetta_scripts::saved_reference_pose;
using protocols::rosetta_scripts::parse_task_operations;
using protocols::rosetta_scripts::parse_score_function;
using protocols::rotamer_recovery::RotamerRecovery;
using protocols::rotamer_recovery::RotamerRecoveryFactory;
using protocols::rotamer_recovery::RRProtocolMover;
using protocols::rotamer_recovery::RRProtocolReferenceStructure;
using protocols::rotamer_recovery::RRReporterSQLite;
using protocols::rotamer_recovery::RRReporterSQLiteOP;
using utility::sql_database::sessionOP;
using utility::tag::TagCOP;
using utility::vector1;
using cppdb::statement;

static Tracer TR("protocols.features.RotamerRecoveryFeatures");

RotamerRecoveryFeatures::RotamerRecoveryFeatures() :
	scfxn_(get_score_function()),
	reporter_(protocols::rotamer_recovery::RRReporterSQLiteOP( new RRReporterSQLite() ) ),
	protocol_(),
	comparer_(),
	task_factory_()
{
	reporter_->set_output_level(rotamer_recovery::OL_features);
}

RotamerRecoveryFeatures::RotamerRecoveryFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(scfxn),
	reporter_(protocols::rotamer_recovery::RRReporterSQLiteOP( new RRReporterSQLite() ) ),
	protocol_(),
	comparer_(),
	task_factory_()
{
	reporter_->set_output_level(rotamer_recovery::OL_features);
}

RotamerRecoveryFeatures::RotamerRecoveryFeatures(
	RotamerRecoveryFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_),
	reporter_(src.reporter_),
	protocol_(src.protocol_),
	comparer_(src.comparer_),
	task_factory_(src.task_factory_)
{}

RotamerRecoveryFeatures::~RotamerRecoveryFeatures() {}

string
RotamerRecoveryFeatures::type_name() const { return "RotamerRecoveryFeatures"; }

vector1<string>
RotamerRecoveryFeatures::features_reporter_dependencies() const {
	vector1<string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
RotamerRecoveryFeatures::initialize_task_factory(
	core::pack::task::TaskFactoryOP task_factory
) {
	task_factory_ = task_factory;
}

void
RotamerRecoveryFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	Movers_map const & movers,
	Pose const & /*pose*/
) {

	scfxn_ = parse_score_function( tag, data );

	RotamerRecoveryFactory * factory(RotamerRecoveryFactory::get_instance());

	if(tag->hasOption("mover") || tag->hasOption("mover_name")){
		if(tag->hasOption("reference_name")){
			throw utility::excn::EXCN_RosettaScriptsOption(
				"Both 'mover_name' and 'reference_name' were supplied. "
				"Please specify at most one to indicate which protocols should "
				"be used to run RotamerRecovery.");
		}

		MoverOP mover = parse_mover(tag->hasOption("mover") ?
			tag->getOption<string>("mover") :
			tag->getOption<string>("mover_name"), movers);
		protocol_ = protocols::rotamer_recovery::RRProtocolOP( new RRProtocolMover(mover) );
	} else if(tag->hasOption("reference_name")){
		if(tag->hasOption("mover")){
			throw utility::excn::EXCN_RosettaScriptsOption(
				"Both 'mover' and 'reference_name' were supplied. "
				"Please specify at most one to indicate which protocols "
				"should be used to run RotamerRecovery.");
		}

		if(tag->hasOption("mover_name")){
			throw utility::excn::EXCN_RosettaScriptsOption(
				"Both 'mover_name' and 'reference_name' were supplied. "
				"Please specify at most one to indicate which protocols "
				"should be used to run RotamerRecovery.");
		}

		if(tag->hasOption("protocol") &&
			!tag->getOption<string>("protocol").compare(
				"RRProtocolReferenceStructure")){
			throw utility::excn::EXCN_RosettaScriptsOption(
				"Specifying 'reference_name' is only compatible with the "
				"'RRProtocolReferenceStructure' protocol.");
		}

		// Use with SavePoseMover
		// WARNING! reference_pose is not initialized until apply time
		PoseCOP reference_pose(saved_reference_pose(tag, data));
		protocol_ = protocols::rotamer_recovery::RRProtocolOP( new RRProtocolReferenceStructure(reference_pose) );
	} else {
		string const & protocol_name(tag->getOption<string>(
				"protocol", "RRProtocolMinPack"));
		if(!protocol_name.compare("RRProtocolMover")){
			throw utility::excn::EXCN_RosettaScriptsOption(
				"Please specify 'mover_name' with the 'RRProtocolMover' "
				"rotamer recovery protocol.");
		} else if(!protocol_name.compare("RRProtocolReferenceStructure")){
			throw utility::excn::EXCN_RosettaScriptsOption(
				"Please specify 'reference_name' with 'RRProtocolReferenceStructure' "
				"rotamer recovery protocol.");
		} else {
			protocol_ = factory->get_rotamer_recovery_protocol(protocol_name);
		}
	}

	//mjo if there are many options to be passed to the components,
	//consider passing the tag to the components themselves to do their
	//own tag parsing
	if(protocol_->get_name() == "RRProtocolRTMin"){
		if(tag->hasOption("nonideal")){
			static_cast<rotamer_recovery::RRProtocolRTMin &>(*protocol_).set_nonideal(tag->getOption<bool>("nonideal"));
		}
		if(tag->hasOption("cartesian")){
			static_cast<rotamer_recovery::RRProtocolRTMin &>(*protocol_).set_cartesian(tag->getOption<bool>("cartesian"));
		}
	}

	if(protocol_->get_name() == "RRProtocolRelax"){
		if(tag->hasOption("nonideal")){
			static_cast<rotamer_recovery::RRProtocolRelax &>(*protocol_).set_nonideal(tag->getOption<bool>("nonideal"));
		}
		if(tag->hasOption("cartesian")){
			static_cast<rotamer_recovery::RRProtocolRelax &>(*protocol_).set_cartesian(tag->getOption<bool>("cartesian"));
		}
	}

	string const & comparer_name(tag->getOption<string>(
			"comparer", "RRComparerAutomorphicRMSD"));
	comparer_ = factory->get_rotamer_recovery_comparer(comparer_name);

	if(tag->hasOption("recovery_threshold")){
		comparer_->set_recovery_threshold(tag->getOption<Real>("recovery_threshold"));
	}

	task_factory_ = parse_task_operations(tag, data);

	if(tag->hasOption("predicted_features_reporter")){
		string report_to_db_name(
			tag->getOption<string>("predicted_features_reporter"));
		MoverOP report_to_db(parse_mover(report_to_db_name, movers));

		if(report_to_db->get_name() != "ReportToDB"){
			stringstream errmsg;
			errmsg
				<< "For tag '" << tag->getName() << "',"
				<< " the predicted_features_reporter item should refer to"
				<< " a previously defined mover of type ReportToDB."
				<< " Instead type specified mover '" << report_to_db_name << "'"
				<< " is a mover of type '" << report_to_db->get_name() << "'."
				<< std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption(errmsg.str());
		}

		ReportToDBOP report_to_db_ptr(
			utility::pointer::static_pointer_cast<ReportToDB>(report_to_db));

		reporter_->set_predicted_report_to_db(report_to_db_ptr);

	}
}


void
RotamerRecoveryFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	reporter_->write_schema_to_db(db_session);
}


Size
RotamerRecoveryFeatures::report_features(
	Pose const & pose_in,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	reporter_->db_session(db_session);


	reporter_->set_struct_id1(struct_id);
	RotamerRecovery rotamer_recovery(protocol_, comparer_, reporter_);

	Pose pose=pose_in;
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose, scfxn_);
	(*scfxn_)(pose);

	if(task_factory_ == 0){
		task_factory_ = core::pack::task::TaskFactoryOP( new TaskFactory() );
	}

	PackerTaskOP packer_task(task_factory_->create_task_and_apply_taskoperations(pose));
	packer_task->restrict_to_repacking();
	packer_task->restrict_to_residues(relevant_residues);
	packer_task->initialize_from_command_line();

	rotamer_recovery.run(pose, *scfxn_, *packer_task);

	return 0;
}


} // namesapce
} // namespace
