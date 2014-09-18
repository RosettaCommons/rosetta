// Headers {{{1

// Core headers
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/svn_version.hh>

// Protocol headers
#include <protocols/canonical_sampling/DbTrajectoryRecorder.hh>
#include <protocols/canonical_sampling/FixedTemperatureController.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/MpiParallelTempering.hh>
#include <protocols/canonical_sampling/ProgressBarObserver.hh>
#include <protocols/canonical_sampling/MultiTempTrialCounter.hh>
#include <protocols/canonical_sampling/TemperatureController.hh>
#include <protocols/kinematic_closure/BalancedKicMover.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/perturbers/WalkingPerturber.hh>
#include <protocols/kinematic_closure/perturbers/WalkingBondAnglePerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/FixedOffsetPivots.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialCounter.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

// Utility headers
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/sql_utils.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>

// MPI headers
#ifdef USEMPI
#include <mpi.h>
#include <protocols/jd2/util.hh>
#endif

// External headers
#include <boost/foreach.hpp>
#include <cppdb/frontend.h>

// C++ headers
#include <ctime>
#include <string>


// Global Names {{{1

using namespace std;
using namespace basic::options;

using core::Size;
using core::Real;
using core::pose::Pose;
using core::import_pose::pose_from_pdb;
using core::pack::task::TaskFactory;
using core::pack::task::TaskFactoryOP;
using core::pack::task::operation::RestrictToRepacking;
using core::pack::task::operation::PreventRepacking;
using core::pack::task::operation::PreventRepackingOP;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;

using protocols::canonical_sampling::DbTrajectoryRecorder;
using protocols::canonical_sampling::DbTrajectoryRecorderOP;
using protocols::canonical_sampling::FixedTemperatureController;
using protocols::canonical_sampling::MetropolisHastingsMover;
using protocols::canonical_sampling::MetropolisHastingsMoverOP;
using protocols::canonical_sampling::MpiParallelTempering;
using protocols::canonical_sampling::ProgressBarObserver;
using protocols::canonical_sampling::MultiTempTrialCounter;
using protocols::canonical_sampling::MultiTempTrialCounterOP;
using protocols::canonical_sampling::MultiTempTrialCounterCOP;
using protocols::canonical_sampling::TemperatureControllerCOP;
using protocols::kinematic_closure::BalancedKicMover;
using protocols::kinematic_closure::BalancedKicMoverOP;
using protocols::kinematic_closure::perturbers::RamaPerturber;
using protocols::kinematic_closure::perturbers::WalkingPerturber;
using protocols::kinematic_closure::perturbers::WalkingBondAnglePerturber;
using protocols::kinematic_closure::pivot_pickers::FixedOffsetPivots;
using protocols::loops::Loop;
using protocols::moves::Mover;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarlo;
using protocols::moves::MonteCarloOP;
using protocols::simple_moves::sidechain_moves::SidechainMover;
using protocols::simple_moves::sidechain_moves::SidechainMoverOP;

using utility::tools::make_vector1;
typedef utility::vector1<Size> IndexList;
static thread_local basic::Tracer tr( "apps.native_ensemble" );

// Options {{{1

OPT_1GRP_KEY(IntegerVector, native_ensemble, loop)
OPT_1GRP_KEY(IntegerVector, native_ensemble, weights)
OPT_1GRP_KEY(String, native_ensemble, algorithm)
OPT_1GRP_KEY(Integer, native_ensemble, iterations)
OPT_1GRP_KEY(Real, native_ensemble, temperature)
OPT_1GRP_KEY(Real, native_ensemble, frequency)
OPT_1GRP_KEY(String, native_ensemble, message)
OPT_1GRP_KEY(Boolean, native_ensemble, quiet)

// }}}1

class NativeEnsemble; // {{{1

typedef utility::pointer::owning_ptr<NativeEnsemble> NativeEnsembleOP;
typedef utility::pointer::owning_ptr<NativeEnsemble const> NativeEnsembleCOP;

class NativeEnsemble : public Mover {

public:
	NativeEnsemble()
		: loop_(0, 0),
		  algorithm_("rama"),
			iterations_(1000),
			temperature_(1),
			backbone_weight_(1),
			sidechain_weight_(1),
			dump_frequency_(1),
			quiet_(false) {}

	void apply(Pose & pose);
	string get_name() const { return "Native Ensemble"; }
	MoverOP fresh_instance() const { return new NativeEnsemble(*this); }

	void set_loop(Loop loop) { loop_ = loop; }
	void set_loop(IndexList values) { set_loop(Loop(values[1], values[2])); }
	void set_algorithm(string value) { algorithm_ = value; }
	void set_iterations(Size value) { iterations_ = value; }
	void set_temperature(Real value) { temperature_ = value; }
	void set_weights(Real bb, Real sc) {
		backbone_weight_ = bb;
		sidechain_weight_ = sc;
	}
	void set_weights(IndexList values) { set_weights(values[1], values[2]); }
	void set_dump_frequency(Size value) { dump_frequency_ = value; }
	void set_input_file(string value) { input_file_ = value; }
	void set_description(string value) { description_ = value; }
	void set_quiet(bool value) { quiet_ = value; }

public:
	static void register_options();

private:
	void show_header() const;
	void show_progress(Size current_iteration) const;
	void show_footer() const;

	void write_schema_to_db() const;
	Size write_job_to_db() const;
	void write_stop_time_to_db(Size id) const;
	void write_temperatures_to_db(Size id, TemperatureControllerCOP tc) const;
	void write_stats_to_db(Size id, MultiTempTrialCounterCOP counter) const;

	bool ok_to_use_db() const;

private:
	Loop loop_;
	string algorithm_;
	Size iterations_;
	Real temperature_;
	Real backbone_weight_;
	Real sidechain_weight_;
	Size dump_frequency_;
	string input_file_;
	string description_;
	bool quiet_;
};

void NativeEnsemble::apply(Pose & pose) { // {{{1

	// Setup the monte carlo run.

	MetropolisHastingsMoverOP canonical_mc = new MetropolisHastingsMover;
	ScoreFunctionOP score_function = core::scoring::get_score_function();
	MonteCarloOP monte_carlo = new MonteCarlo(pose, *score_function, 0);

	canonical_mc->set_ntrials(iterations_);
	canonical_mc->set_monte_carlo(monte_carlo);
	canonical_mc->set_tempering(
			new FixedTemperatureController(temperature_));

	if (loop_.start() <= 1) loop_ = Loop(2, pose.total_residue() - 1);
	pose.update_residue_neighbors();

	// Setup the standard backbone and sidechain moves.

	SidechainMoverOP sidechain_mover = new SidechainMover;
	BalancedKicMoverOP backbone_mover = new BalancedKicMover;

	TaskFactoryOP task_factory = new TaskFactory;
	PreventRepackingOP prevent_repacking = new PreventRepacking();

	utility::vector1<bool> is_near_loop = 
		protocols::loops::select_loop_residues(pose, loop_, true, 10.0);

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if (not is_near_loop[i]) prevent_repacking->include_residue(i);
	}

	task_factory->push_back(new RestrictToRepacking);
	task_factory->push_back(prevent_repacking);

	sidechain_mover->set_task_factory(task_factory);
	backbone_mover->set_loop(loop_);

	canonical_mc->add_mover(sidechain_mover, sidechain_weight_);
	canonical_mc->add_mover(backbone_mover, backbone_weight_);

	// Customize the backbone and tempering moves.

	if (algorithm_ == "rama") {
		backbone_mover->add_perturber(new RamaPerturber);
	}
	else if (algorithm_ == "walking") {
		backbone_mover->add_perturber(new WalkingPerturber(5));
	}
	else if (algorithm_ == "rama/para-temp") {
		canonical_mc->set_tempering(new MpiParallelTempering);
	}
	else if (algorithm_ == "walking/para-temp") {
		canonical_mc->set_tempering(new MpiParallelTempering);
		backbone_mover->add_perturber(new WalkingPerturber(5));
	}
	else if (algorithm_ == "rama/bond-angle") {
		backbone_mover->add_perturber(new RamaPerturber);
		backbone_mover->add_perturber(new WalkingBondAnglePerturber);
	}
	else if (algorithm_ == "walking/bond-angle") {
		backbone_mover->add_perturber(new WalkingPerturber(5));
		backbone_mover->add_perturber(new WalkingBondAnglePerturber);
	}
	else {
		utility_exit_with_message("Unknown algorithm: " + algorithm_);
	}

	TemperatureControllerCOP tempering = canonical_mc->tempering();
	MultiTempTrialCounterOP trial_counter = new MultiTempTrialCounter(tempering);
	monte_carlo->set_counter(trial_counter);

	// Setup output.

	write_schema_to_db();
	Size id = write_job_to_db();
	DbTrajectoryRecorderOP db_output = new DbTrajectoryRecorder;
	db_output->set_job_id(id);
	db_output->set_temp_level(1);
	db_output->stride(dump_frequency_);
	db_output->cache_limit(1);
	canonical_mc->add_observer(db_output);
	canonical_mc->add_observer(new ProgressBarObserver);
	show_header();

	// Run the monte carlo loop.

	canonical_mc->apply(pose);

	// Finalize output and clean up.

	write_stop_time_to_db(id);
	write_temperatures_to_db(id, canonical_mc->tempering());
	write_stats_to_db(id, trial_counter);
	show_footer();
}

void NativeEnsemble::register_options() { // {{{1
	MetropolisHastingsMover::register_options();
	MpiParallelTempering::register_options();
	SidechainMover::register_options();
	BalancedKicMover::register_options();
}

void NativeEnsemble::show_header() const { // {{{1
	tr << "Residues:    " << loop_.start() << "/" << loop_.stop() << endl;
	tr << "Algorithm:   " << algorithm_ << endl;
	tr << "Iterations:  " << iterations_ << endl;
	tr << "Temperature: " << temperature_ << endl;
	tr << "BB:SC Ratio: " << backbone_weight_ << ":" << sidechain_weight_ << endl;
}

void NativeEnsemble::show_progress(Size current_iteration) const { // {{{1
	if (quiet_ == true) return;
	cerr << "\r[" << current_iteration << "/" << iterations_ << "]";
}

void NativeEnsemble::show_footer() const { // {{{1
	cerr << endl;
}

void NativeEnsemble::write_schema_to_db() const { // {{{1
	using utility::sql_database::sessionOP;
	using namespace basic::database::schema_generator;

	if (not ok_to_use_db()) return;

	sessionOP db_session = basic::database::get_db_session();

	// Define the jobs table.
	// Strictly speaking, the loop_begin and loop_end fields should be in their 
	// own table.  As it is, this schema cannot support multiple loop regions.
	{
		Column id("id", new DbBigInt, false, true);
		Column command("command", new DbText);
		Column revision("revision", new DbText);
		Column algorithm("algorithm", new DbText);
		Column start_time("start_time", new DbText);
		Column stop_time("stop_time", new DbText);
		Column input_file("input_file", new DbText);
		Column loop_begin("loop_begin", new DbInteger);
		Column loop_end("loop_end", new DbInteger);
		Column iterations("iterations", new DbInteger);
		Column frames("frames", new DbInteger);
		Column description("description", new DbText);

		Schema jobs("jobs", PrimaryKey(id));

		jobs.add_column(id);
		jobs.add_column(command);
		jobs.add_column(revision);
		jobs.add_column(algorithm);
		jobs.add_column(start_time);
		jobs.add_column(stop_time);
		jobs.add_column(input_file);
		jobs.add_column(loop_begin);
		jobs.add_column(loop_end);
		jobs.add_column(iterations);
		jobs.add_column(frames);
		jobs.add_column(description);
		jobs.write(db_session);
	}

	// Define the moves table.
	{
		Column job_id("job_id", new DbBigInt, false);
		Column type("type", new DbText(50), false);
		Column temp_level("temp_level", new DbInteger, false);
		Column num_trials("num_trials", new DbBigInt, false);
		Column num_accepted("num_accepted", new DbBigInt, false);

		Columns key_columns = make_vector1(job_id, type, temp_level);
		Schema moves("moves", PrimaryKey(key_columns));

		moves.add_column(num_trials);
		moves.add_column(num_accepted);
		moves.write(db_session);
	}

	// Define the temperatures table.
	{
		Column job_id("job_id", new DbBigInt, false);
		Column level("level", new DbInteger, false);
		Column temperature("temperature", new DbReal, false);

		Columns key_columns = make_vector1(job_id, level);
		Schema temperatures("temperatures", PrimaryKey(key_columns));

		temperatures.add_column(temperature);
		temperatures.write(db_session);
	}
}

Size NativeEnsemble::write_job_to_db() const { // {{{1
	using utility::sql_database::sessionOP;
	using basic::database::safely_prepare_statement;
	using basic::database::safely_write_to_database;

	Size job_id = 0;

	if (ok_to_use_db()) {
		sessionOP db_session = basic::database::get_db_session();
		db_session->begin_transaction();

		// Create a new entry in the jobs table, fill it in with all the 
		// information that's available right now, and return the job id.

		string const insert_string =
			"INSERT INTO jobs ("
				"command, revision, algorithm, start_time, input_file, "
				"loop_begin, loop_end, iterations, frames, description) "
			"VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);";

		string const command = option.get_argv();
		string const revision = core::minirosetta_svn_version();
		time_t now = time(0); string const start_time = ctime(&now);

		cppdb::statement insert_statement =
			safely_prepare_statement(insert_string, db_session);

		insert_statement.bind( 1, command);
		insert_statement.bind( 2, revision);
		insert_statement.bind( 3, algorithm_);
		insert_statement.bind( 4, start_time);
		insert_statement.bind( 5, input_file_);
		insert_statement.bind( 6, loop_.start());
		insert_statement.bind( 7, loop_.stop());
		insert_statement.bind( 8, iterations_);
		insert_statement.bind( 9, iterations_ / dump_frequency_);
		insert_statement.bind(10, description_);

		safely_write_to_database(insert_statement);

		job_id = insert_statement.sequence_last("jobs_id_seq");
		db_session->force_commit_transaction();
	}

	// If MPI is being used, broadcast the job id from the root node.  It is 
	// assumed that the root node also inserted the job into the database and 
	// received a valid id.  This is enforced by ok_to_use_db().

#ifdef USEMPI
	MPI_Bcast(&job_id, 1, MPI_INT, 0, protocols::jd2::current_mpi_comm());
#endif

	return job_id;
}

void NativeEnsemble::write_stop_time_to_db(Size id) const { // {{{1
	using utility::sql_database::sessionOP;
	using basic::database::safely_prepare_statement;
	using basic::database::safely_write_to_database;

	if (not ok_to_use_db()) return;

	// Record the stop time in the jobs table.

	sessionOP db_session = basic::database::get_db_session();

	string const update_string =
		"UPDATE jobs SET stop_time=? WHERE id=?;";

	time_t now = time(0);
	string const stop_time = ctime(&now);

	cppdb::statement update_statement =
		safely_prepare_statement(update_string, db_session);

	update_statement.bind(1, stop_time);
	update_statement.bind(2, id);

	safely_write_to_database(update_statement);
}

void NativeEnsemble::write_temperatures_to_db( // {{{1
		Size id, TemperatureControllerCOP temp_controller) const {

	using utility::sql_database::sessionOP;
	using basic::database::safely_prepare_statement;
	using basic::database::safely_write_to_database;

	if (not ok_to_use_db()) return;

	sessionOP db_session = basic::database::get_db_session();

	string const insert_string = 
		"INSERT INTO temperatures (job_id, level, temperature) "
		"VALUES (?, ?, ?);";

	for (Size level = 1; level <= temp_controller->n_temp_levels(); level++) {
		cppdb::statement insert_statement =
			safely_prepare_statement(insert_string, db_session);

		insert_statement.bind(1, id);
		insert_statement.bind(2, level);
		insert_statement.bind(3, temp_controller->temperature(level));

		safely_write_to_database(insert_statement);
	}
}

void NativeEnsemble::write_stats_to_db( // {{{1
		Size id, MultiTempTrialCounterCOP counter) const {

	using utility::sql_database::sessionOP;
	using basic::database::safely_prepare_statement;
	using basic::database::safely_write_to_database;
	using protocols::moves::TrialCounter;

	if (not ok_to_use_db()) return;

	// Create a new entry in the moves table with information about which moves 
	// were used in this simulation and how effective they each were.

	sessionOP db_session = basic::database::get_db_session();

	string const insert_string =
		"INSERT INTO moves (job_id, type, temp_level, num_trials, num_accepted) "
		"VALUES (?, ?, ?, ?, ?);";

	for (Size temp_level = 1;
	     temp_level <= counter->num_temp_levels();
	     temp_level++) {

		TrialCounter const & temp_counter = counter->temp_level(temp_level);

		BOOST_FOREACH(string tag, temp_counter.tags()) {
			Size num_trials = temp_counter.trial(tag);
			Size num_accepted = temp_counter.accepted(tag);

			cppdb::statement insert_statement =
				safely_prepare_statement(insert_string, db_session);

			insert_statement.bind(1, id);
			insert_statement.bind(2, tag);
			insert_statement.bind(3, temp_level);
			insert_statement.bind(4, num_trials);
			insert_statement.bind(5, num_accepted);

			safely_write_to_database(insert_statement);
		}
	}
}

bool NativeEnsemble::ok_to_use_db() const { // {{{1
	// This script should only ever access the database from one thread at a 
	// time.  The reason is that this everything this script writes to the 
	// database is metadata about a job, and is doesn't make sense to write that 
	// information from each child node.  This function simply returns true in 
	// only one thread per job.  If MPI is being used, it will return true in the 
	// root node.  Otherwise it will always return true.

#ifdef USEMPI
	int rank; MPI_Comm_rank(protocols::jd2::current_mpi_comm(), &rank);
	return rank == 0;
#else
	return true;
#endif
}
// }}}1

int main(int argc, char * argv[]) { // {{{1
	IndexList default_loop = make_vector1(0, 0);
	IndexList default_weights = make_vector1(1, 9);

	NEW_OPT(native_ensemble::loop, "Residues to sample", default_loop);
	NEW_OPT(native_ensemble::weights, "BB-to-SC move ratio", default_weights);
	NEW_OPT(native_ensemble::algorithm, "Backbone sampling algorithm", "rama");
	NEW_OPT(native_ensemble::iterations, "Monte Carlo iterations", 1000);
	NEW_OPT(native_ensemble::temperature, "Monte Carlo temperature", 1);
	NEW_OPT(native_ensemble::frequency, "Pose output frequency", 1);
	NEW_OPT(native_ensemble::message, "Description of the job", "");
	NEW_OPT(native_ensemble::quiet, "Hide progress updates", false);

	NativeEnsemble::register_options();

	devel::init(argc, argv);

	NativeEnsembleOP app = new NativeEnsemble;

	if (not option[OptionKeys::in::file::s].active()) {
		utility_exit_with_message("Specify an input PDB using the '-s' flag.");
	}
	if (not option[OptionKeys::inout::dbms::database_name].user()) {
		option[OptionKeys::inout::dbms::database_name].value("sandbox.db");
	}

	Pose pose;
	pose_from_pdb(pose, option[OptionKeys::in::file::s][1].name());

	if (option[OptionKeys::native_ensemble::loop].active()) {
		app->set_loop(option[OptionKeys::native_ensemble::loop]());
	}

	option[OptionKeys::jd2::no_output].value(false);

	app->set_weights(option[OptionKeys::native_ensemble::weights]());
	app->set_algorithm(option[OptionKeys::native_ensemble::algorithm]());
	app->set_iterations(option[OptionKeys::native_ensemble::iterations]());
	app->set_temperature(option[OptionKeys::native_ensemble::temperature]());
	app->set_dump_frequency(option[OptionKeys::native_ensemble::frequency]());
	app->set_input_file(option[OptionKeys::in::file::s][1].name());
	app->set_description(option[OptionKeys::native_ensemble::message]());
	app->set_quiet(option[OptionKeys::native_ensemble::quiet]());

	app->apply(pose);

#ifdef USEMPI
	MPI_Finalize();
#endif
}
// }}}1
