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
/// @author Matthew O'Meara


#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/features/ReportToDB.hh>
#include <string>

// Setup Mover
#include <protocols/features/ReportToDBCreator.hh>

#include <basic/database/sql_utils.hh>

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/features/FeaturesReporterFactory.hh>
#include <protocols/features/ProteinRMSDFeatures.fwd.hh>
#include <protocols/features/ProtocolFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <numeric>
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.fwd.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/foreach.hpp>


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
using core::Size;
using core::pack::task::PackerTaskCOP;
using core::pack::task::TaskFactory;
using core::pose::Pose;
using core::pose::PoseOP;
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
using std::string;
using std::endl;
using std::accumulate;
using std::stringstream;
using utility::file::FileName;
using utility::vector0;
using utility::vector1;
using utility::tag::TagPtr;
using utility::sql_database::DatabaseSessionManager;
using utility::sql_database::session;

static Tracer TR("protocols.features.ReportToDB");

Size ReportToDB::struct_id_ = 0;
Size ReportToDB::protocol_id_ = 0;
bool ReportToDB::autoincrement_struct_id_ = true;
bool ReportToDB::protocol_table_initialized_ = false;


ReportToDB::ReportToDB():
	Mover("ReportToDB"),
	database_fname_("FeatureStatistics.db3"),
	database_mode_("sqlite3"),
	sample_source_("Rosetta: Unknown Protocol"),
	scfxn_(getScoreFunction()),
	use_transactions_(true),
	cache_size_(2000),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(string const & name):
	Mover(name),
	database_fname_("FeatureStatistics.db3"),
	database_mode_("sqlite3"),
	sample_source_("Rosetta: Unknown Protocol"),
	scfxn_( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) ),
	cache_size_(2000),
	use_transactions_(true),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB(
	string const & name,
	string const & database_fname,
	string const & sample_source,
	ScoreFunctionOP scfxn,
	bool use_transactions,
	Size cache_size) :
	Mover(name),
	database_fname_(database_fname),
	database_mode_("sqlite3"),
	sample_source_(sample_source),
	scfxn_(scfxn),
	use_transactions_(use_transactions),
	cache_size_(cache_size),
	task_factory_(new TaskFactory()),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	features_reporters_(),
	initialized( false )
{
	initialize_reporters();
}

ReportToDB::ReportToDB( ReportToDB const & src):
	Mover(src),
	database_fname_(src.database_fname_),
	database_mode_(src.database_mode_),
	sample_source_(src.sample_source_),
	scfxn_(new ScoreFunction(* src.scfxn_)),
	use_transactions_(src.use_transactions_),
	cache_size_(src.cache_size_),
	task_factory_(src.task_factory_),
	features_reporter_factory_(FeaturesReporterFactory::get_instance()),
	protocol_features_(src.protocol_features_),
	structure_features_(src.structure_features_),
	features_reporters_(src.features_reporters_),
	initialized(src.initialized)
{}

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
ReportToDB::parse_db_tag_item(
	TagPtr const tag){

	if( tag->hasOption("db") ){
		database_fname_ = tag->getOption<string>("db");
	} else {
		TR << "Field 'db' required for use of ReportToDB mover in Rosetta Scripts." << endl;
	}
}

void
ReportToDB::parse_sample_source_tag_item(
	TagPtr const tag){
	if( tag->hasOption("sample_source") ){
		sample_source_ = tag->getOption<string>("sample_source");
	} else {
		TR << "Field 'sample_source' required for use of ReportToDB in Rosetta Scripts." << endl;
		TR << "The sample_sourOBce should describe where the samples came from. To access the description run \"sqlite3 'select description from sample_source' fname.db3\"" << endl;
		TR << "For example: Top4400 natives from Richardson Lab. Reduce placed hydrogens with -correct flag." << endl;
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
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		protocol_id_ = mpi_rank;
	}
#endif

}

void
ReportToDB::parse_first_struct_id_tag_item(
	TagPtr const tag) {

	if(tag->hasOption("first_struct_id")){
		Size first_struct_id = tag->getOption<Size>("first_struct_id");
		// Initialize struct_id_ to first_struct_id only if this is the
		// first time ReportToDB has been executed. It is sufficient to
		// test autoincrement_struct_id_.
		if (autoincrement_struct_id_){
#ifdef USEMPI
			int mpi_rank(0);
			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
			struct_id_ = first_struct_id + mpi_rank - 1; // minus 1 to setup invariant
#else
			struct_id_ = first_struct_id - 1; // minus 1 to setup invariant
#endif
			autoincrement_struct_id_ = false;
		}
	}
#ifdef USEMPI
	else if (autoincrement_struct_id_){
		int mpi_rank(0);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		struct_id_ = mpi_rank - 1;
		autoincrement_struct_id_ = false;
	}
#endif

}

void
ReportToDB::parse_db_mode_tag_item(
	TagPtr const tag) {
	if(tag->hasOption("db_mode")){
		database_mode_ = tag->getOption<std::string>("db_mode");
	} else {
		database_mode_ = "sqlite3";
	}
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

	// Name of output features database:
	// EXAMPLE: db=features_<sample_source>.db3
	// REQUIRED
	parse_db_tag_item(tag);

	// Description of features database
	// EXAMPLE: sample_source="This is a description of the sample source."
	// RECOMMENDED
	parse_sample_source_tag_item(tag);

	// Manually control the id of associated with this protocol
	// EXAMPLE: protocol_id=6
	// OPTIONAL
	parse_protocol_id_tag_item(tag);

	// Manually control the first structure id
	// EXAMPLE: first_struct_id=100
	// OPTIONAL
	parse_first_struct_id_tag_item(tag);

	// The database backend to use ('sqlite3', 'mysql' etc)
	// EXAMPLE: db_mode=sqlite3
	parse_db_mode_tag_item(tag);

	// Use transactions to group database i/o to be more efficient. Turning them off OBcan help debugging.
	// EXAMPLE: use_transactions=true
	// DEFAULT: TRUE
	parse_use_transactions_tag_item(tag);

	// Specify the maximum number 1k pages to keep in memory before writing to disk
	// EXAMPLE: cache_size=1000000  // this uses ~ 1GB of memory
	// DEFAULT: 2000
	parse_cache_size_tag_item(tag);

	task_factory_ = parse_task_operations(tag, data);

	vector0< TagPtr >::const_iterator begin=tag->getTags().begin();
	vector0< TagPtr >::const_iterator end=tag->getTags().end();

	for(; begin != end; ++begin){
		TagPtr feature_tag= *begin;
		//	foreach(TagPtr const & feature_tag, tag->getTags()){

		if(feature_tag->getName() != "feature"){
			TR.Error << "Please include only tags with name 'feature' as subtags of ReportToDB" << endl;
			TR.Error << "Tag with name '" << feature_tag->getName() << "' is invalid" << endl;
			utility_exit();
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
		if(dependency == "ProtocolFeatures" || dependency == "StructureFeatures"){
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
				<< "The dependencies for the '" << test_features_reporter->type_name() << "'"
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
	// the protocol and structure features are special
	protocol_features_ = new ProtocolFeatures();
	structure_features_ = new StructureFeatures();
}

void
ReportToDB::initialize_database(
	utility::sql_database::sessionOP db_session
){
	if (!initialized){
		protocol_features_->write_schema_to_db(db_session);
		structure_features_->write_schema_to_db(db_session);
		foreach( FeaturesReporterOP const & reporter, features_reporters_ ){
			reporter->write_schema_to_db(db_session);
		}
		initialized = true;
	}

}

void
ReportToDB::apply( Pose& pose ){

	utility::sql_database::sessionOP db_session(basic::database::get_db_session(database_fname_, database_mode_, false, true));

	// Make sure energy objects are initialized
	(*scfxn_)(pose);

	PackerTaskCOP task(task_factory_->create_task_and_apply_taskoperations(pose));
	vector1< bool > relevant_residues(task->repacking_residues());

	TR << "Reporting features for "
		<< accumulate(relevant_residues.begin(), relevant_residues.end(), 0)
		<< " of the " << pose.total_residue()
		<< " total residues in the pose." << endl;

	if(use_transactions_) db_session->begin();
	initialize_database(db_session);
	if(use_transactions_) db_session->commit();

	if(use_transactions_) db_session->begin();

	stringstream stmt_ss; stmt_ss << "PRAGMA cache_size = " << cache_size_ << ";";
	statement stmt = (*db_session) << stmt_ss.str(); stmt.exec();

	if(!protocol_table_initialized_){
		protocol_id_ = protocol_features_->report_features(
			protocol_id_, db_session);
		protocol_table_initialized_ = true;

		statement stmt = (*db_session) <<
			"CREATE TABLE IF NOT EXISTS sample_source ( description TEXT );";
		stmt.exec();
		stmt = (*db_session)
			<< "INSERT INTO sample_source VALUES (?);"
			<< sample_source_;
		stmt.exec();
	}

	if(autoincrement_struct_id_){
		struct_id_ = structure_features_->report_features(
		  pose, relevant_residues, protocol_id_, db_session);
	} else {
#ifdef USEMPI
		int mpi_size(0);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
		struct_id_+=mpi_size;
#else
		struct_id_++;
#endif
		structure_features_->report_features(
			pose, relevant_residues, struct_id_, protocol_id_, db_session);
	}

	for(Size i=1; i <= features_reporters_.size(); ++i){
		features_reporters_[i]->report_features(
			pose, relevant_residues, struct_id_, db_session);
	}

	if(use_transactions_) db_session->commit();
}

} // namespace
} // namespace
