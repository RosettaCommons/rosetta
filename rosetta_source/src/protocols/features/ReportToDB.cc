// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

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


// Platform Headers
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.fwd.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/chemical/ResidueType.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
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
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Numeric Headers
#include <numeric>

// Boost Headers
#include <boost/foreach.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>


// C++ Headers
#include <sstream>


namespace protocols{
namespace features{

std::string
ReportToDBCreator::keyname() const
{
	return ReportToDBCreator::mover_name();
}

moves::MoverOP
ReportToDBCreator::create_mover() const {
	return new ReportToDB;
}

std::string
ReportToDBCreator::mover_name()
{
	return "ReportToDB";
}

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script
#define foreach BOOST_FOREACH

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
using core::pose::Pose;
using core::pose::PoseOP;
using core::pose::symmetry::is_symmetric;
using core::scoring::getScoreFunction;
using core::scoring::ScoreFunctionFactory;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using core::scoring::STANDARD_WTS;
using cppdb::cppdb_error;
using cppdb::statement;
using cppdb::result;
using protocols::features::FeaturesReporterOP;
using protocols::features::ProteinRMSDFeatures;
using protocols::features::ProtocolFeatures;
using protocols::features::StructureFeatures;
using protocols::features::FeaturesReporterFactory;
using protocols::jd2::JobDistributor;
using protocols::moves::MoverOP;
using protocols::moves::DataMap;
using protocols::moves::Movers_map;
using protocols::rosetta_scripts::parse_task_operations;
using protocols::rosetta_scripts::parse_score_function;
using basic::database::parse_database_connection;
using std::string;
using std::endl;
using std::accumulate;
using std::stringstream;
using boost::uuids::uuid;
using utility::file::FileName;
using utility::vector0;
using utility::vector1;
using utility::tag::TagPtr;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::session;
using utility::sql_database::sessionOP;

static Tracer TR("protocols.features.ReportToDB");

ReportToDB::ReportToDB():
	Mover("ReportToDB"),
	db_session_(),
	sample_source_("Rosetta: Unknown Protocol"),
	use_transactions_(true),
	cache_size_(2000),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(string const & name):
	Mover(name),
	db_session_(),
	sample_source_("Rosetta: Unknown Protocol"),
	use_transactions_(true),
	cache_size_(2000),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(
	string const & name,
	sessionOP db_session,
	string const & sample_source,
	bool use_transactions,
	Size cache_size) :
	Mover(name),
	db_session_(db_session),
	sample_source_(sample_source),
	use_transactions_(use_transactions),
	cache_size_(cache_size),
	remove_xray_virt_(false),
	protocol_id_(0),
	batch_id_(0),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB( ReportToDB const & src):
	Mover(src),
	db_session_(src.db_session_),
	sample_source_(src.sample_source_),
	name_(src.name_),
	use_transactions_(src.use_transactions_),
	cache_size_(src.cache_size_),
	remove_xray_virt_(src.remove_xray_virt_),
	protocol_id_(src.protocol_id_),
	batch_id_(src.batch_id_),
	task_factory_(src.task_factory_),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	protocol_features_(src.protocol_features_),
	batch_features_(src.batch_features_),
	structure_features_(src.structure_features_),
	features_reporters_(src.features_reporters_),
	initialized(src.initialized)
{
	TR << "ReportToDB copy ctor called" << std::endl;
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

MoverOP ReportToDB::fresh_instance() const { return new ReportToDB; }

MoverOP ReportToDB::clone() const
{
	return new ReportToDB( *this );
}

void
ReportToDB::parse_sample_source_tag_item(
	TagPtr const tag){
	if( tag->hasOption("sample_source") ){
		sample_source_ = tag->getOption<string>("sample_source");
	} else {
		TR << "Field 'sample_source' required for use of ReportToDB in Rosetta Scripts." << endl;
		TR << "The sample_source should describe where the samples came from. To access the description run \"sqlite3 'select description from sample_source' fname.db3\"" << endl;
		TR << "For example: Top4400 natives from Richardson Lab. Reduce placed hydrogens with -correct flag." << endl;
	}
}

void
ReportToDB::parse_name_tag_item(TagPtr const tag){
		if( tag->hasOption("name") ){
				name_=tag->getOption<string>("name");
		} else {
				TR << "Field 'name' required for use of ReportToDB in Rosetta Scripts." << endl;
		}
}

void
ReportToDB::parse_protocol_id_tag_item(
	TagPtr const tag){

	if(tag->hasOption("protocol_id")){
#ifdef USEMPI
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		protocol_id_ = tag->getOption<Size>("protocol_id") + mpi_rank;
#else
		protocol_id_ = tag->getOption<Size>("protocol_id");
#endif
	}
#ifdef USEMPI
	else {
				protocol_id_ = 0;
	}
#endif

}

void
ReportToDB::parse_use_transactions_tag_item(
	TagPtr const tag) {
	if(tag->hasOption("use_transactions")){
		use_transactions_ = tag->getOption<bool>("use_transactions");
	}
}

void
ReportToDB::parse_cache_size_tag_item(
	TagPtr const tag) {
	if(tag->hasOption("cache_size")){
		cache_size_ = tag->getOption<bool>("cache_size");
	}
}

void
ReportToDB::parse_remove_xray_virt_tag_item(
									  TagPtr const tag) {
	if(tag->hasOption("remove_xray_virt")){
		remove_xray_virt_ = tag->getOption<bool>("remove_xray_virt");
	}
}

/// Allow ReportToDB to be called from RosettaScripts
/// See
void
ReportToDB::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose )
{
	if(tag->hasOption("db")){
		throw utility::excn::EXCN_RosettaScriptsOption("The 'db' tag has been deprecated. Please use 'database_name' instead.");
	}

	if(tag->hasOption("db_mode")){
		throw utility::excn::EXCN_RosettaScriptsOption("The 'database_mode' tag has been deprecated. Please use 'database_mode' instead.");
	}

	if(tag->hasOption("separate_db_per_mpi_process")){
		throw utility::excn::EXCN_RosettaScriptsOption("The 'parse_separate_db_per_mpi_process' tag has been deprecated. Please use 'database_parse_separate_db_per_mpi_process' instead.");
	}

	// Name of output features database:
	// EXAMPLE: db=features_<sample_source>.db3
	// REQUIRED
	if(tag->hasOption("resource_description")){
		std::string resource_description = tag->getOption<string>("resource_description");
		if ( ! basic::resource_manager::ResourceManager::get_instance()->has_resource_with_description( resource_description ) )
		{
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a resource_description of " + resource_description +
					" for ReportToDB, but the ResourceManager doesn't have a resource with that description" );
		}
		db_session_ = basic::resource_manager::get_resource< utility::sql_database::session >( resource_description );
	}
	else{
		db_session_ = parse_database_connection(tag);
	}

	// Description of features database
	// EXAMPLE: sample_source="This is a description of the sample source."
	// RECOMMENDED
	parse_sample_source_tag_item(tag);

		// Name of report to db mover. A new batch will be created for each uniquely named
		// ReportToDb mover
	// EXAMPLE: name="initial_feature_extraction"
	// RECOMMENDED
	parse_name_tag_item(tag);

	// Manually control the id of associated with this protocol
	// EXAMPLE: protocol_id=6
	// OPTIONAL default is to autoincrement the protocol_id in the protocols table
	parse_protocol_id_tag_item(tag);

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

	task_factory_ = parse_task_operations(tag, data);

	vector0< TagPtr >::const_iterator begin=tag->getTags().begin();
	vector0< TagPtr >::const_iterator end=tag->getTags().end();

	for(; begin != end; ++begin){
		TagPtr feature_tag= *begin;
		//	foreach(TagPtr const & feature_tag, tag->getTags()){

		if(feature_tag->getName() != "feature"){
			TR.Error << "Please include only tags with name 'feature' as subtags of ReportToDB" << endl;
			TR.Error << "Tag with name '" << feature_tag->getName() << "' is invalid" << endl;
			throw utility::excn::EXCN_RosettaScriptsOption("");
		}

		FeaturesReporterOP features_reporter(
			features_reporter_factory_->get_features_reporter(
				feature_tag, data, filters, movers, pose));

		check_features_reporter_dependencies(features_reporter);

		// TODO IMPLMENT THIS:
		//check_multiple_features_reporter_definitions(features_reporter);

		features_reporters_.push_back(features_reporter);

	}

}


void
ReportToDB::check_features_reporter_dependencies(
	FeaturesReporterOP test_features_reporter
) const {

	foreach(string const dependency,
		test_features_reporter->features_reporter_dependencies()){

		// These are defined by default
		if(dependency == "ProtocolFeatures" || dependency == "BatchFeatures" || dependency == "StructureFeatures"){
			continue;
		}

		bool exists(false);
		foreach(FeaturesReporterOP features_reporter, features_reporters_){
			if(features_reporter->type_name() == dependency){
				exists = true;
				break;
			}
		}
		if(!exists){
			stringstream error_msg;
			error_msg
				<< "For batch '" << name_ << "'," << endl
				<< "the dependencies for the '" << test_features_reporter->type_name() << "'"
				<< " reporter are not satisfied because the '" << dependency << "' has not been defined yet." << endl
				<< "These are the FeaturesReporters that have been defined:" << endl
				<< "\tProtocolFeatures (included by default)" << endl
				<< "\tStructureFeatures (included by default)" << endl;
			foreach(FeaturesReporterOP features_reporter, features_reporters_){
				error_msg
					<< "\t" << features_reporter->type_name() << endl;
			}
			utility_exit_with_message(error_msg.str());
		}
	}
}

void
ReportToDB::initialize_reporters()
{
	// the protocols, batches, and structure features are special
	protocol_features_ = new ProtocolFeatures();
	batch_features_ = new BatchFeatures();
	structure_features_ = new StructureFeatures();
}

void
ReportToDB::initialize_database(){

	if (!initialized){
		if(use_transactions_) db_session_->begin();

		protocol_features_->write_schema_to_db(db_session_);
		batch_features_->write_schema_to_db(db_session_);
		structure_features_->write_schema_to_db(db_session_);

		write_linking_tables();

		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			reporter->write_schema_to_db(db_session_);
		}

		if(use_transactions_) db_session_->commit();

		initialized = true;
	}
}

vector1< bool >
ReportToDB::initialize_pose(
	Pose & pose
) const {

	if (remove_xray_virt_) {
		TR << "Removing virtual residue left behind by xray refinement" << endl;
		while (pose.residue( pose.total_residue() ).aa() == core::chemical::aa_vrt )
			pose.conformation().delete_residue_slow( pose.total_residue() );
	}
	
	PackerTaskCOP task(task_factory_->create_task_and_apply_taskoperations(pose));
	vector1< bool > relevant_residues(task->repacking_residues());
	
	TR
		<< "Reporting features for "
		<< accumulate(relevant_residues.begin(), relevant_residues.end(), 0)
		<< " of the " << pose.total_residue()
		<< " total residues in the pose "
		<< JobDistributor::get_instance()->current_output_name()
		<< " for batch '" << name_ << "'." << endl;

	return relevant_residues;
}

/// @detail The 'features_reporters' table lists the type_names of the
/// all defined features reporters. The 'batch_reports' table link the
/// features reporters with each batch defined in the 'batches' table.
void
ReportToDB::write_features_reporters_table() const {
	using namespace basic::database::schema_generator;

	Schema features_reporters(
		"features_reporters",
		PrimaryKey( Column("report_name", new DbTextKey())));

	features_reporters.write(db_session_);

	//Only report features that aren't already in the database
	string select_string =
		"SELECT *\n"
		"FROM\n"
		"	features_reporters\n"
		"WHERE\n"
		"	report_name = ?;";
	statement select_stmt(safely_prepare_statement(select_string, db_session_));

	string insert_string = "INSERT INTO features_reporters (report_name) VALUES (?);";
	statement insert_stmt(safely_prepare_statement(insert_string, db_session_));

	foreach(FeaturesReporterOP const & reporter, features_reporters_){
		string const report_name(reporter->type_name());
		select_stmt.bind(1,report_name);

		result res(safely_read_from_database(select_stmt));
		if(!res.next()) {
			insert_stmt.bind(1, report_name);
			safely_write_to_database(insert_stmt);
		}
	}
}

void
ReportToDB::write_batch_reports_table() const {
	using namespace basic::database::schema_generator;

	Schema batch_reports("batch_reports");
	Column report_name("report_name", new DbTextKey());
	Column batch_id("batch_id", new DbInteger());

	batch_reports.add_foreign_key(
		ForeignKey(batch_id, "batches", "batch_id", true /*defer*/));
	batch_reports.add_foreign_key(
		ForeignKey(report_name, "features_reporters", "report_name", true /*defer*/));

	vector1<Column> batch_reports_unique;
	batch_reports_unique.push_back(batch_id);
	batch_reports_unique.push_back(report_name);
	batch_reports.add_constraint( new UniqueConstraint(batch_reports_unique) );

	batch_reports.write(db_session_);
}

void
ReportToDB::write_linking_tables() const {

	try{
		write_features_reporters_table();
	} catch(cppdb_error error){
		stringstream err_msg;
		err_msg
			<< "The ReportToDB Mover failed to write the 'features_reporters' table "
			<< "to the database for batch '" << name_ << "'." << endl
			<< "Error Message:" << endl << error.what() << endl;
		utility_exit_with_message(err_msg.str());
	}

	try{
		write_batch_reports_table();
	} catch(cppdb_error error){
		stringstream err_msg;
		err_msg
			<< "The ReportToDB Mover failed to write the 'batch_reports' table "
			<< "to the database." << endl
			<< "Error Message:" << endl << error.what() << endl;
		utility_exit_with_message(err_msg.str());
	}
}

void
ReportToDB::apply( Pose& pose ){

	vector1<bool> relevant_residues(initialize_pose(pose));

	initialize_database();

	if(use_transactions_) db_session_->begin();

	set_cache_size(db_session_, cache_size_);

	std::pair<Size, Size> ids = get_protocol_and_batch_id(name_, db_session_);
	protocol_id_ = ids.first;
	batch_id_ = ids.second;

	uuid struct_id = report_structure_features(relevant_residues);

	report_features(pose, struct_id, relevant_residues);

	if(use_transactions_) db_session_->commit();
}

uuid
ReportToDB::report_structure_features(
	vector1<bool> const & relevant_residues
) const {
	uuid struct_id;
	try {
		struct_id = structure_features_->report_features(
			relevant_residues, batch_id_, db_session_);
	} catch (cppdb_error error){
		stringstream err_msg;
		err_msg
			<< "Failed to report structure features for:" << endl
			<< "\tprotocol_id: '" << protocol_id_ << "'" << endl
			<< "\tbatch name: '" << name_ << "'" <<  endl
			<< "\tbatch_id: '" << batch_id_ << "'" << endl
			<< "Error Message:" << endl << error.what() << endl;
		utility_exit_with_message(err_msg.str());
	}
	return struct_id;
}

void
ReportToDB::report_features(
	Pose const & pose,
	uuid const struct_id,
	utility::vector1<bool> const & relevant_residues
) const {

//	string batch_reports_string =
//		"INSERT INTO batch_reports (batch_id, report_name) VALUES (?,?);";
//	statement batch_reports_stmt(
//		safely_prepare_statement(batch_reports_string, db_session));

	for(Size i=1; i <= features_reporters_.size(); ++i){
		string report_name = features_reporters_[i]->type_name();

		TR << "Reporting " << report_name << std::endl;

		try {
			features_reporters_[i]->report_features(
				pose, relevant_residues, struct_id, db_session_);
		} catch (cppdb_error error){
			stringstream err_msg;
			err_msg
				<< "Failed to report features for the "
				<< "'" << report_name << "' reporter with:" << endl
				<< "with:" << endl
				<< "\tprotocol_id: '" << protocol_id_ << "' " << endl
				<< "\tbatch name: '" << name_ << "' " << endl
				<< "\tbatch_id: '" << batch_id_ << "'" << endl
				<< "\tstruct_id: '" << struct_id << "'" << endl
				<< "Error Message:" << endl << error.what() << endl;
			utility_exit_with_message(err_msg.str());
		}

		//Need to check for preexisting entry to avoid constraint failure caused by having multiple structures in a batch. Alternatively, we could add struct_id to batch_reports table
		//        batch_reports_stmt.bind(1, batch_id_);
		//        batch_reports_stmt.bind(2, report_name);
		//        safely_write_to_database(batch_reports_stmt);
	}

}

} // namespace
} // namespace
