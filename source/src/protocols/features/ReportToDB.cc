// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ReportToDB.cc
///
/// @brief  report all data to a database
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/features/ReportToDB.hh>
#include <string>

// Setup Mover
#include <protocols/features/ReportToDBCreator.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.fwd.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/chemical/ResidueType.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/ProteinRMSDFeatures.fwd.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/BatchFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Numeric Headers
#include <numeric>

// Boost Headers
#include <boost/foreach.hpp>


// C++ Headers
#include <sstream>
#include <algorithm>
#include <cctype>

namespace protocols {
namespace features {

std::string
ReportToDBCreator::keyname() const
{
	return ReportToDBCreator::mover_name();
}

moves::MoverOP
ReportToDBCreator::create_mover() const {
	return moves::MoverOP( new ReportToDB );
}

std::string
ReportToDBCreator::mover_name()
{
	return "ReportToDB";
}

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script

using basic::T;
using basic::Tracer;
using basic::Error;
using basic::Warning;
using basic::datacache::CacheableString;
using basic::database::safely_prepare_statement;
using basic::database::safely_write_to_database;
using basic::database::safely_read_from_database;
using basic::database::get_db_session;
using basic::database::set_cache_size;
using core::Size;
using core::pack::task::PackerTaskCOP;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::symmetry::is_symmetric;
using cppdb::cppdb_error;
using cppdb::statement;
using cppdb::result;
using protocols::jd2::JobDistributor;
using protocols::moves::MoverOP;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using protocols::rosetta_scripts::parse_task_operations;
using protocols::rosetta_scripts::parse_score_function;
using basic::database::parse_database_connection;
using std::string;
using std::endl;
using std::accumulate;
using std::stringstream;
using utility::file::FileName;
using utility::vector0;
using utility::vector1;
using utility::tag::TagCOP;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::session;
using utility::sql_database::sessionOP;

static Tracer TR("protocols.features.ReportToDB");

ReportToDB::ReportToDB():
	Mover("ReportToDB"),
	db_session_(),
	batch_name_("features"),
	batch_description_("Rosetta: Unknown Protocol"),
	use_transactions_(true),
	cache_size_(2000),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	custom_structure_tag_(false),
	structure_tag_(""),
	custom_structure_input_tag_(false),
	structure_input_tag_(""),
	last_struct_id_(0),
	task_factory_(core::pack::task::TaskFactoryOP( new TaskFactory() )),
	relevant_residues_mode_(RelevantResiduesMode::Exclusive),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(string const & type):
	Mover(type),
	db_session_(),
	batch_name_("features"),
	batch_description_("Rosetta: Unknown Protocol"),
	use_transactions_(true),
	cache_size_(2000),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	custom_structure_tag_(false),
	structure_tag_(""),
	custom_structure_input_tag_(false),
	structure_input_tag_(""),
	last_struct_id_(0),
	task_factory_(core::pack::task::TaskFactoryOP( new TaskFactory() )),
	relevant_residues_mode_(RelevantResiduesMode::Exclusive),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(
	sessionOP db_session,
	string const & batch_name,
	string const & batch_description,
	bool use_transactions,
	Size cache_size) :
	Mover("ReportToDB"),
	db_session_(db_session),
	batch_name_(batch_name),
	batch_description_(batch_description),
	use_transactions_(use_transactions),
	cache_size_(cache_size),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	custom_structure_tag_(false),
	structure_tag_(""),
	custom_structure_input_tag_(false),
	structure_input_tag_(""),
	last_struct_id_(0),
	task_factory_(core::pack::task::TaskFactoryOP( new TaskFactory() )),
	relevant_residues_mode_(RelevantResiduesMode::Exclusive),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	if ( batch_name == "" ) {
		utility::excn::EXCN_BadInput("Failed to create ReportToDB instance because the batch name must not be ''.");
	}

	initialize_reporters();
}

ReportToDB::ReportToDB(
	string const & name,
	sessionOP db_session,
	string const & batch_name,
	string const & batch_description,
	bool use_transactions,
	Size cache_size
) :
	Mover(name),
	db_session_(db_session),
	batch_name_(batch_name),
	batch_description_(batch_description),
	use_transactions_(use_transactions),
	cache_size_(cache_size),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	custom_structure_tag_(false),
	structure_tag_(""),
	custom_structure_input_tag_(false),
	structure_input_tag_(""),
	last_struct_id_(0),
	task_factory_(core::pack::task::TaskFactoryOP( new TaskFactory() )),
	relevant_residues_mode_(RelevantResiduesMode::Exclusive),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	if ( batch_name == "" ) {
		utility::excn::EXCN_BadInput("Failed to create ReportToDB instance because the batch name must not be ''.");
	}

	initialize_reporters();
}

ReportToDB::ReportToDB( ReportToDB const & src):
	Mover(src),
	db_session_(src.db_session_),
	batch_name_(src.batch_name_),
	batch_description_(src.batch_description_),
	use_transactions_(src.use_transactions_),
	cache_size_(src.cache_size_),
	remove_xray_virt_(src.remove_xray_virt_),
	protocol_id_(src.protocol_id_),
	batch_id_(src.batch_id_),
	custom_structure_tag_(src.custom_structure_tag_),
	structure_tag_(src.structure_tag_),
	custom_structure_input_tag_(src.custom_structure_input_tag_),
	structure_input_tag_(src.structure_input_tag_),
	last_struct_id_(src.last_struct_id_),
	task_factory_(src.task_factory_),
	relevant_residues_mode_(src.relevant_residues_mode_),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	protocol_features_(src.protocol_features_),
	batch_features_(src.batch_features_),
	structure_features_(src.structure_features_),
	features_reporters_(src.features_reporters_),
	initialized(src.initialized)
{
}

ReportToDB::~ReportToDB(){}

void
ReportToDB::register_options() const {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// This mover is equiped to work with the Rosetta Scripts interface
	option.add_relevant( parser::protocol );

	//TODO call relevant_options on FeaturesMover objects
}

MoverOP ReportToDB::fresh_instance() const { return MoverOP( new ReportToDB ); }

MoverOP ReportToDB::clone() const
{
	return MoverOP( new ReportToDB( *this ) );
}

void
ReportToDB::set_batch_name(
	std::string const & name
) {
	if ( name == "" ) {
		utility::excn::EXCN_BadInput("Setting the batch name for a ReporToDB instance to '' is not allowed.");
	}

	batch_name_ = name;
}

std::string
ReportToDB::get_batch_name() const {
	return batch_name_;
}

void
ReportToDB::set_batch_description(
	std::string const & batch_description
) {
	batch_description_ = batch_description;
}

std::string
ReportToDB::get_batch_description() const {
	return batch_description_;
}

void
ReportToDB::set_relevant_residues_task_factory(
	TaskFactoryOP task_factory
) {
	task_factory_ = task_factory;
}

TaskFactoryOP
ReportToDB::get_relevant_residues_task_factory() const {
	return task_factory_;
}

void
ReportToDB::set_relevant_residues(
	vector1<bool> const & relevant_residues
) {
	relevant_residues_ = relevant_residues;
}

vector1<bool>
ReportToDB::get_relevant_residues() const {
	return relevant_residues_;
}

void
ReportToDB::set_relevant_residues_mode(
	RelevantResiduesMode::T setting
) {
	relevant_residues_mode_ = setting;
}

RelevantResiduesMode::T
ReportToDB::get_relevant_residues_mode() const {
	return relevant_residues_mode_;
}

void
ReportToDB::set_structure_tag(
	std::string const & setting
) {
	structure_tag_ = setting;
	custom_structure_tag_ = true;
}

std::string
ReportToDB::get_structure_tag() const {
	return structure_tag_;
}

void
ReportToDB::set_structure_input_tag(
	std::string const & setting
) {
	structure_input_tag_ = setting;
	custom_structure_input_tag_ = true;
}

std::string
ReportToDB::get_structure_input_tag() const {
	return structure_input_tag_;
}

void
ReportToDB::add_features_reporter(
	FeaturesReporterOP features_reporter
) {
	check_features_reporter_dependencies(features_reporter);
	// TODO IMPLMENT THIS:
	//check_multiple_features_reporter_definitions(features_reporter);

	features_reporters_.push_back(features_reporter);
}

void
ReportToDB::parse_batch_description_tag_item(
	TagCOP const tag
){
	batch_description_ = tag->getOption<string>("batch_description", "");
}

void
ReportToDB::parse_batch_name_tag_item(TagCOP const tag){
	if ( tag->hasOption("name") ) {
		batch_name_ = tag->getOption<string>("name");
		utility::replace_in(batch_name_, ' ', "_");
	} else {
		TR << "Field 'name' required for use of ReportToDB in Rosetta Scripts and it will be used as the name field in the batches table." << endl;
	}
}

void
ReportToDB::parse_protocol_id_tag_item(
	TagCOP const tag){

	if ( tag->hasOption("protocol_id") ) {
		//#ifdef USEMPI
		//  int mpi_rank(0);
		//  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		//  protocol_id_ = tag->getOption<Size>("protocol_id") + mpi_rank;
		//#else
		protocol_id_ = tag->getOption<Size>("protocol_id");
		//#endif
	}
#ifdef USEMPI
	else {
				protocol_id_ = 0;
	}
#endif

}

void
ReportToDB::parse_batch_id_tag_item(
	TagCOP const tag){

	if ( tag->hasOption("batch_id") ) {
		//#ifdef USEMPI
		//  int mpi_rank(0);
		//  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		//  batch_id_ = tag->getOption<Size>("batch_id") + mpi_rank;
		//#else
		batch_id_ = tag->getOption<Size>("batch_id");
		//#endif
	}
#ifdef USEMPI
	else {
				batch_id_ = 0;
	}
#endif

}

void
ReportToDB::parse_use_transactions_tag_item(
	TagCOP const tag) {
	if ( tag->hasOption("use_transactions") ) {
		use_transactions_ = tag->getOption<bool>("use_transactions");
	}
}

void
ReportToDB::parse_cache_size_tag_item(
	TagCOP const tag) {
	if ( tag->hasOption("cache_size") ) {
		cache_size_ = tag->getOption<bool>("cache_size");
	}
}

void
ReportToDB::parse_remove_xray_virt_tag_item(
	TagCOP const tag) {
	if ( tag->hasOption("remove_xray_virt") ) {
		remove_xray_virt_ = tag->getOption<bool>("remove_xray_virt");
	}
}

void
ReportToDB::parse_relevant_residues_mode_tag_item(
	TagCOP const tag) {
	string rel_res_mode = tag->getOption<string>(
		"relevant_residues_mode", "explicit");
	std::transform(
		rel_res_mode.begin(), rel_res_mode.end(), rel_res_mode.begin(),
		::toupper);

	if ( rel_res_mode == "EXPLICIT" ) {
		relevant_residues_mode_ = RelevantResiduesMode::Exclusive;
	} else if ( rel_res_mode == "IMPLICIT" ) {
		relevant_residues_mode_ = RelevantResiduesMode::Inclusive;
	} else {
		throw utility::excn::EXCN_RosettaScriptsOption
			( "Bad value for relevant_residues_mode: '" + rel_res_mode + "'. It must be either 'EXPLICIT' or 'IMPLICIT' (case insensitive). This indicates which features should be reported given the relevant residue specification (determined by the packable residues in the given task operation:\n\tEXCLUSIVE: All residues in a feature must be specified as 'relevant'. (DEFAULT)\n\tINCLUSIVE: At least one residue in the the feature must be specified as 'relevant' to be reported.");
	}
}

/// Allow ReportToDB to be called from RosettaScripts
/// See
void
ReportToDB::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose )
{
	if ( tag->hasOption("db") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The 'db' tag has been deprecated. Please use 'database_name' instead.");
	}

	if ( tag->hasOption("db_mode") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The 'db_mode' tag has been deprecated. Please use 'database_mode' instead.");
	}

	if ( tag->hasOption("separate_db_per_mpi_process") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The 'parse_separate_db_per_mpi_process' tag has been deprecated. Please use 'database_parse_separate_db_per_mpi_process' instead.");
	}

	if ( tag->hasOption("sample_source") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The 'sample_source' tag has been deprecated. Please use 'batch_description' instead.");
	}

	// Name of output features database:
	// EXAMPLE: db=features_<batch_description>.db3
	// REQUIRED
	if ( tag->hasOption("resource_description") ) {
		std::string resource_description = tag->getOption<string>("resource_description");
		if ( ! basic::resource_manager::ResourceManager::get_instance()->has_resource_with_description( resource_description ) ) {
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a resource_description of " + resource_description +
				" for ReportToDB, but the ResourceManager doesn't have a resource with that description" );
		}
		db_session_ = basic::resource_manager::get_resource< utility::sql_database::session >( resource_description );
	} else {
		db_session_ = parse_database_connection(tag);
	}

	// Name of report to db mover. A new batch will be created for each uniquely named
	// ReportToDb mover
	// EXAMPLE: name="initial_feature_extraction"
	// RECOMMENDED
	parse_batch_name_tag_item(tag);

	// Description of features database
	// EXAMPLE: batch_description="This is a description of the sample source."
	// RECOMMENDED
	parse_batch_description_tag_item(tag);

	// Manually control the protocol_id associated with this ReportToDB tag
	// EXAMPLE: protocol_id=6
	// OPTIONAL default is to autoincrement the protocol_id in the protocols table
	parse_protocol_id_tag_item(tag);

	// Manually control the batch_id associated with this ReportToDB tag
	// EXAMPLE: batch_id=6
	// OPTIONAL default is to autoincrement the batch_id in the batches table
	parse_batch_id_tag_item(tag);

	// Use transactions to group database i/o to be more efficient. Turning them off can help debugging.
	// EXAMPLE: use_transactions=true
	// DEFAULT: TRUE
	parse_use_transactions_tag_item(tag);

	// Specify the maximum number 1k pages to keep in memory before writing to disk
	// EXAMPLE: cache_size=1000000  // this uses ~ 1GB of memory
	// DEFAULT: 2000
	parse_cache_size_tag_item(tag);

	// Remove virtual residue attached during xray refine process
	// EXAMPLE: remove_xray_virt=true
	// DEFAULT: FALSE
	parse_remove_xray_virt_tag_item(tag);

	// Determine what features are reported given the relevant residues
	// EXAMPLE: relevant_residues_mode=implicit
	// DEFAULT: explicit
	parse_relevant_residues_mode_tag_item(tag);

	//This is probably not necessary, but do it for completeness.
	structure_features_->set_relevant_residues_mode(relevant_residues_mode_);

	task_factory_ = parse_task_operations(tag, data);

	vector0< TagCOP >::const_iterator begin=tag->getTags().begin();
	vector0< TagCOP >::const_iterator end=tag->getTags().end();

	for ( ; begin != end; ++begin ) {
		TagCOP feature_tag= *begin;
		// BOOST_FOREACH(TagCOP const & feature_tag, tag->getTags()){

		if ( feature_tag->getName() != "feature" ) {
			TR.Error << "Please include only tags with name 'feature' as subtags of ReportToDB" << endl;
			TR.Error << "Tag with name '" << feature_tag->getName() << "' is invalid" << endl;
			throw utility::excn::EXCN_RosettaScriptsOption("");
		}

		FeaturesReporterOP features_reporter(
			features_reporter_factory_->get_features_reporter(
			feature_tag, data, filters, movers, pose));
		features_reporter->set_relevant_residues_mode(relevant_residues_mode_);

		add_features_reporter(features_reporter);
	}
}


void
ReportToDB::check_features_reporter_dependencies(
	FeaturesReporterOP test_features_reporter
) const {

	BOOST_FOREACH ( string const dependency,
			test_features_reporter->features_reporter_dependencies() ) {

		// These are defined by default
		if ( dependency == "ProtocolFeatures" || dependency == "BatchFeatures" || dependency == "StructureFeatures" ) {
			continue;
		}

		bool exists(false);
		BOOST_FOREACH ( FeaturesReporterOP features_reporter, features_reporters_ ) {
			if ( features_reporter->type_name() == dependency ) {
				exists = true;
				break;
			}
		}
		if ( !exists ) {
			stringstream error_msg;
			error_msg
				<< "For batch '" << batch_name_ << "'," << endl
				<< "the dependencies for the '" << test_features_reporter->type_name() << "'"
				<< " reporter are not satisfied because the '" << dependency << "' has not been defined yet." << endl
				<< "These are the FeaturesReporters that have been defined:" << endl
				<< "\tProtocolFeatures (included by default)" << endl
				<< "\tStructureFeatures (included by default)" << endl;
			BOOST_FOREACH ( FeaturesReporterOP features_reporter, features_reporters_ ) {
				error_msg
					<< "\t" << features_reporter->type_name() << endl;
			}

			// Theoretically this could be thrown in a non-RosettaScripts context,
			// but that would be bad coding on the C++ developer's part
			throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
		}
	}
}

void
ReportToDB::initialize_reporters()
{
	// the protocols, batches, and structure features are special
	protocol_features_ = protocols::features::ProtocolFeaturesOP( new ProtocolFeatures() );
	batch_features_ = protocols::features::BatchFeaturesOP( new BatchFeatures() );
	structure_features_ = protocols::features::StructureFeaturesOP( new StructureFeatures() );
}

void
ReportToDB::initialize_database(){

	if ( basic::options::option[basic::options::OptionKeys::out::database_protocol_id].user() ) {
		if ( protocol_id_ != 0 ) {
			utility_exit_with_message("You have specified a protocol id in both the ReportToDB tag and the options system. Please use only one.");
		}
		protocol_id_ = basic::options::option[basic::options::OptionKeys::out::database_protocol_id];
	}

	// KAB - check if db_session_ has been set.
	// Necessary if parse_my_tag not called because mover is being called
	// outside of RosettaScripts
	if ( !db_session_ ) {
		db_session_ = basic::database::get_db_session();
	}

	if ( !initialized ) {
		//  if(use_transactions_) db_session_->begin();

		protocol_features_->write_schema_to_db(db_session_, protocol_id_);
		batch_features_->write_schema_to_db(db_session_, batch_id_);
		structure_features_->write_schema_to_db(db_session_);

		//deferred until after batch_id has been set
		//  write_linking_tables();

		BOOST_FOREACH ( FeaturesReporterOP const & reporter, features_reporters_ ) {
			reporter->write_schema_to_db(db_session_);
		}

		//  if(use_transactions_) db_session_->commit();

		initialized = true;
	}
}

void
ReportToDB::initialize_pose(
	Pose & pose
) const {
	pose.update_residue_neighbors(); // As some of the sub-features may need neighbor information
	if ( remove_xray_virt_ ) {
		TR << "Removing virtual residue left behind by xray refinement" << endl;
		while ( pose.residue( pose.total_residue() ).aa() == core::chemical::aa_vrt )
				pose.conformation().delete_residue_slow( pose.total_residue() );
	}
}

vector1< bool >
ReportToDB::initialize_relevant_residues(
	Pose const & pose
) {
	if ( relevant_residues_.size() && task_factory_ ) {
		TR.Warning
			<< "Both relevant_residues and task_factory has been specified;"
			<< " using the relevant residues." << std::endl;
	}
	if ( relevant_residues_.size() ) {
		if ( relevant_residues_.size() != pose.total_residue() ) {
			TR.Warning
				<< "The size of relevant_residues is: " << relevant_residues_.size()
				<< " while the pose has '" << pose.total_residue() << "' residues,"
				<< " verify that this mismatch intended."
				<< " Resetting the relevant_residues to be all residues in the pose."
				<< std::endl;
			return vector1< bool >(pose.total_residue(), true);
		} else {
			return relevant_residues_;
		}
	} else {
		PackerTaskCOP task(task_factory_->create_task_and_apply_taskoperations(pose));
		return task->repacking_residues();
	}
}

void
ReportToDB::apply( Pose& pose ){

	initialize_pose(pose);
	vector1< bool > relevant_residues(initialize_relevant_residues(pose));

	TR
		<< "Reporting features for "
		<< accumulate(relevant_residues.begin(), relevant_residues.end(), 0)
		<< " of the " << pose.total_residue()
		<< " total residues in the pose "
		<< JobDistributor::get_instance()->current_output_name()
		<< " for batch '" << batch_name_ << "'." << endl;


	initialize_database();

	set_cache_size(db_session_, cache_size_);

	try{
		//If no protcol and batch id were set by the user, let them be autogenerated
		if ( protocol_id_==0 && batch_id_==0 ) {
			std::pair<Size, Size> ids = get_protocol_and_batch_id(batch_name_, batch_description_, features_reporters_, db_session_);
			protocol_id_ = ids.first;
			batch_id_ = ids.second;
		} else {
			set_protocol_and_batch_id(protocol_id_, batch_id_, batch_name_, batch_description_, features_reporters_, db_session_);
		}
	}
catch(cppdb::cppdb_error & except)
{
	TR.Warning << "There was an error writing protocol and batch id. You are likely using MPI mode "
		<< "with a user-specified protocol_id and/or batch_id. This workflow is for advanced users "
		<< "only and is not recommended. The cppdb error is: " << except.what() << std::endl;
}
//Write linking tables after we have a valid batch_id
// write_linking_tables();

	ensure_structure_tags_are_ready();
	StructureID struct_id = report_structure_features();
	report_features(pose, struct_id, relevant_residues);

	last_struct_id_ = struct_id;
}

void
ReportToDB::ensure_structure_tags_are_ready()
{
	if ( (!custom_structure_tag_) || (structure_tag_ == "") ) {
		structure_tag_ = JobDistributor::get_instance()->current_output_name();
	}

	if ( (!custom_structure_input_tag_) || (structure_input_tag_ == "") ) {
		structure_input_tag_ = JobDistributor::get_instance()->current_job()->input_tag();
	}
}

StructureID
ReportToDB::report_structure_features() const {
	StructureID struct_id;
	try {
		if ( use_transactions_ ) {
			db_session_->begin_transaction();
		}

		struct_id = structure_features_->report_features(
			batch_id_, db_session_, structure_tag_, structure_input_tag_
		);

		if ( use_transactions_ ) {
			db_session_->commit_transaction();
		}
	} catch (cppdb_error & error){
		db_session_->rollback();
		stringstream err_msg;
		err_msg
			<< "Failed to report structure features for:" << endl
			<< "\tprotocol_id: '" << protocol_id_ << "'" << endl
			<< "\tbatch name: '" << batch_name_ << "'" <<  endl
			<< "\tbatch description: '" << batch_description_ << "'" <<  endl
			<< "\tbatch_id: '" << batch_id_ << "'" << endl
			<< "Error Message:" << endl << error.what() << endl;
		utility_exit_with_message(err_msg.str());
	} catch (utility::excn::EXCN_Base & error){
		if ( use_transactions_ ) {
			db_session_->rollback();
		}
		stringstream err_msg;
		err_msg
			<< "Failed to report structure features for:" << endl
			<< "\tprotocol_id: '" << protocol_id_ << "'" << endl
			<< "\tbatch name: '" << batch_name_ << "'" <<  endl
			<< "\tbatch description: '" << batch_description_ << "'" <<  endl
			<< "\tbatch_id: '" << batch_id_ << "'" << endl
			<< "Error Message:" << endl << error << endl;
		utility_exit_with_message(err_msg.str());
	}
	return struct_id;
}

void
ReportToDB::report_features(
	Pose const & pose,
	StructureID const struct_id,
	utility::vector1<bool> const & relevant_residues
) const {

	for ( Size i=1; i <= features_reporters_.size(); ++i ) {
		string report_name = features_reporters_[i]->type_name();

		TR << "Reporting " << report_name << std::endl;

		try {
			if ( use_transactions_ ) {
				db_session_->begin_transaction();
				features_reporters_[i]->report_features(
					pose, relevant_residues, struct_id, db_session_);
				db_session_->commit_transaction();
			} else {
				features_reporters_[i]->report_features(
					pose, relevant_residues, struct_id, db_session_);
			}
		} catch (cppdb_error & error){
			stringstream err_msg;
			err_msg
				<< "Failed to report features for the "
				<< "'" << report_name << "' reporter with:" << endl
				<< "with:" << endl
				<< "\tprotocol_id: '" << protocol_id_ << "' " << endl
				<< "\tbatch name: '" << batch_name_ << "' " << endl
				<< "\tbatch description: '" << batch_description_ << "' " << endl
				<< "\tbatch_id: '" << batch_id_ << "'" << endl
				<< "\tstruct_id: '" << struct_id << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}
	}
}

StructureID
ReportToDB::get_last_struct_id() const {
	return last_struct_id_;
}


} // namespace
} // namespace
