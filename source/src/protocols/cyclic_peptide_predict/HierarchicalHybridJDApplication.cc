// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.cc
/// @brief Application-level code for the simple_cycpep_predict app, MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.  This
/// version uses MPI for cross-communication with parallel processes.
/// On 4 Aug 2017, work commenced to make this appliction compatible with a job distribution scheme in which the last level of the hierarchy splits
/// jobs over many threads (hierarchical process-level parallelism plus thread-level parallelism).
/// On 29 Oct 2018, this code was moved to the HierarchicalHybridJDApplication base class, from which both the SimpleCycpepPredictApplication_MPI and
/// HelicalBundlePredictApplication_MPI classes derive.
/// On 16 June 2020, this code was updated to replace "emperor"/"master"/"slave", which some people found objectionable.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

// Unit Headers
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_RMSDToBestSummary.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_SASASummary.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_PNearToArbitraryStateSummary.hh>
#include <protocols/cyclic_peptide_predict/PNearCalculator.hh>

// Basic Headers
#include <basic/Tracer.hh>
#include <basic/random/init_random_generator.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/PDBInfo.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/util.hh>
#include <core/chemical/AA.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/simple_metrics/metrics/SasaMetric.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/numbers.hh>
#include <utility/pointer/memory.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>
#include <basic/options/keys/helical_bundle_predict.OptionKeys.gen.hh>

//numeric headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <stdio.h>
#include <chrono>

#ifdef MULTI_THREADED
#include <thread>
#endif

#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__)
typedef std::chrono::system_clock HHJDA_CLOCK;
#else
typedef std::chrono::high_resolution_clock HHJDA_CLOCK;
#endif

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Tracer initialization constructor
/// @details Base class must be initialized to point to derived class tracers.
HierarchicalHybridJDApplication::HierarchicalHybridJDApplication(
	basic::Tracer & tracer,
	basic::Tracer & summary_tracer
) :
	derivedTR_(tracer),
	derivedTR_summary_(summary_tracer),
	MPI_rank_(0),
	MPI_n_procs_(0),
	worker_job_count_(0),
	scorefxn_(nullptr),
	hierarchy_level_(0),
	total_hierarchy_levels_(0),
	procs_per_hierarchy_level_(),
	batchsize_(0),
	my_children_(),
	my_parent_(0),
	sequence_(""),
	native_(nullptr),
	sort_type_(SORT_BY_ENERGIES),
	select_highest_(false),
	output_fraction_(1.0),
	output_filename_(""),
	lambda_(0.5),
	kbt_(1.0),
	compute_rmsd_to_lowest_(false),
	compute_pnear_to_this_fract_(0.0),
	compute_sasa_metrics_(false)
#ifdef MULTI_THREADED
	,
	threads_per_worker_proc_(1),
	results_mutex_(),
	joblist_mutex_(),
	use_const_random_seed_(false),
	random_seed_offset_(0),
	random_seed_(11111111),
	rgtype_("")
#endif
	//TODO -- Initialize vars here as they are added.
{

	scorefxn_ = core::scoring::get_score_function(); //Reads from file.
#ifdef MULTI_THREADED
	init_random_info_from_options();
#endif
}

/// @brief Constructor with options
///
HierarchicalHybridJDApplication::HierarchicalHybridJDApplication(
	basic::Tracer & tracer,
	basic::Tracer & summary_tracer,
	int const MPI_rank,
	int const MPI_n_procs,
	core::scoring::ScoreFunctionCOP sfxn_in,
	core::Size total_hierarchy_levels,
	utility::vector1 < core::Size > const & procs_per_hierarchy_level,
	utility::vector1< core::Size > const &batchsize_per_level,
	std::string const &sort_type,
	bool const select_highest,
	core::Real const output_fraction,
	std::string const &output_filename,
	core::Real const lambda,
	core::Real const kbt,
	bool const compute_rmsd_to_lowest,
	core::Real const compute_pnear_to_this_fract,
	bool const compute_sasa_metrics,
#ifdef MULTI_THREADED
	core::Size const threads_per_worker_proc
#else
	core::Size const /*threads_per_worker_proc*/
#endif
) :
	derivedTR_(tracer),
	derivedTR_summary_(summary_tracer),
	MPI_rank_( MPI_rank ),
	MPI_n_procs_( MPI_n_procs ),
	worker_job_count_( 0 ),
	scorefxn_( sfxn_in != nullptr ? sfxn_in->clone() : nullptr ),
	hierarchy_level_( 0 ),
	total_hierarchy_levels_( total_hierarchy_levels ),
	procs_per_hierarchy_level_(),
	batchsize_(0),
	my_children_(),
	my_parent_(0),
	sequence_(""),
	native_(nullptr),
	sort_type_(SORT_BY_ENERGIES),
	select_highest_(select_highest),
	output_fraction_(1.0),
	output_filename_( output_filename ),
	lambda_(lambda),
	kbt_(kbt),
	compute_rmsd_to_lowest_( compute_rmsd_to_lowest ),
	compute_pnear_to_this_fract_(compute_pnear_to_this_fract),
	compute_sasa_metrics_( compute_sasa_metrics )
#ifdef MULTI_THREADED
	,
	threads_per_worker_proc_(threads_per_worker_proc),
	results_mutex_(),
	joblist_mutex_(),
	use_const_random_seed_(false),
	random_seed_offset_(0),
	random_seed_(11111111),
	rgtype_("")
#endif
{
#ifdef MULTI_THREADED
	init_random_info_from_options();
#endif
	set_sort_type( sort_type );
	set_output_fraction( output_fraction );
	set_procs_per_hierarchy_level( procs_per_hierarchy_level );
	assign_level_children_and_parent();
	set_batchsize( batchsize_per_level );
	std::string const errmsg( "Error in constructor for HierarchicalHybridJDApplication class:  " );
	runtime_assert_string_msg( output_filename_ != "", errmsg + "The output filename cannot be empty." );
	runtime_assert_string_msg( compute_pnear_to_this_fract_ == 0.0 || !compute_rmsd_to_lowest_, errmsg  + "The -compute_pnear_to_this_fract option cannot be used with the -compute_rmsd_to_lowest option." );
}

/// @brief Explicit virtual destructor.
///
HierarchicalHybridJDApplication::~HierarchicalHybridJDApplication() {}


/// @brief Explicit copy constructor.
///
HierarchicalHybridJDApplication::HierarchicalHybridJDApplication(
	HierarchicalHybridJDApplication const &src
) :
	VirtualBase( src ),
	derivedTR_(src.derivedTR_),
	derivedTR_summary_(src.derivedTR_summary_),
	MPI_rank_( src.MPI_rank_ ),
	MPI_n_procs_( src.MPI_n_procs_ ),
	worker_job_count_( src.worker_job_count_ ),
	scorefxn_(), //Cloned below
	hierarchy_level_( 0 ),
	total_hierarchy_levels_( src.total_hierarchy_levels_ ),
	procs_per_hierarchy_level_(), //Assigned below
	batchsize_(src.batchsize_),
	my_children_(), //Assigned below
	my_parent_(0),
	sequence_(src.sequence_),
	native_(), //Cloned below.
	sort_type_(src.sort_type_),
	select_highest_(src.select_highest_),
	output_fraction_(src.output_fraction_),
	output_filename_(src.output_filename_),
	lambda_(src.lambda_),
	kbt_(src.kbt_),
	compute_rmsd_to_lowest_( src.compute_rmsd_to_lowest_ ),
	compute_pnear_to_this_fract_( src.compute_pnear_to_this_fract_ ),
	compute_sasa_metrics_( src.compute_sasa_metrics_ )
#ifdef MULTI_THREADED
	,
	threads_per_worker_proc_( src.threads_per_worker_proc_ ),
	results_mutex_(), //Don't copy mutexes
	joblist_mutex_(), //Don't copy mutexes.
	use_const_random_seed_(src.use_const_random_seed_),
	random_seed_offset_(src.random_seed_offset_),
	random_seed_(src.random_seed_),
	rgtype_(src.rgtype_)
#endif
	//TODO -- copy variables here as they are added.
{
	set_procs_per_hierarchy_level( src.procs_per_hierarchy_level_ );
	assign_level_children_and_parent();
	if(src.native_) native_ = src.native_->clone();
	if(src.scorefxn_) scorefxn_ = src.scorefxn_->clone();

	runtime_assert( hierarchy_level_ == src.hierarchy_level_ );
	runtime_assert( my_children_ == src.my_children_ );
	runtime_assert( my_parent_ == src.my_parent_ );
	runtime_assert( output_fraction_ >= 0.0 && output_fraction_ <= 1.0 );
	runtime_assert_string_msg( output_filename_ != "", "Error in copy constructor for HierarchicalHybridJDApplication class: The output filename cannot be empty." );
}

#ifdef MULTI_THREADED
/// @brief Initialize private member variables related to thread random seeds from the options
/// system.  Does nothing if this isn't a multi-threaded compilation.
void
HierarchicalHybridJDApplication::init_random_info_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ run::constant_seed ].active() ) use_const_random_seed_ = option[ run::constant_seed ]();
	if ( option[ run::seed_offset ].active() ) random_seed_offset_ = option[ run::seed_offset ]();
	if ( option[ run::jran ].active() )  random_seed_ = option[ run::jran ]();
	rgtype_ = option[ run::rng ]();
}
#endif


/// @brief Actually run the application.
/// @details On worker nodes, this creates an instance of the relevant application and runs that.  Nodes higher
/// in the communications hierarchy are just involved in sending out jobs and pulling in results.
/// @note The run() function is nonconst due to some setup steps that it performs.  It then calls const run functions
/// for director, manager, and worker nodes.
void
HierarchicalHybridJDApplication::run() {
	debug_assert(scorefxn_ != nullptr); //Must be assigned before running.

	get_sequence(); //Director reads sequence from disk and broadcasts it to everyone else.
	get_native(); //Director reads native pose from disk and broadcasts it to everyone else.
	get_protocol_specific_settings(); //Get settings specific for the application in question.  For example, for simple_cycpep_predict, if we're doing design, then the director reads design settings from disk and broadcasts them to everyone else.

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Proc " << MPI_rank_ << " reports that it is configured and ready to start." << std::endl;

	if( i_am_director() ) {
		run_director();
	} else if ( i_am_worker() ) {
		//core::Size k(5); //For GDB debugging in MPI mode.
		//while(k!=4) {  } //For GDB debugging in MPI mode.
		run_worker();
	} else {
		run_intermediate_manager();
	}
	return;
}

/// ------------- Public Methods ---------------------

/// @brief Set the sort type, by string.
///
void
HierarchicalHybridJDApplication::set_sort_type(
	std::string const &sort_type
) {
	runtime_assert_string_msg( sort_type == "energy" ||  sort_type == "rmsd" || sort_type == "hbonds" || sort_type == "cis_peptide_bonds",
		"Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::set_sort_type(): The sort type " + sort_type + " is unknown." );

	if( sort_type == "energy" ) {
		set_sort_type( SORT_BY_ENERGIES );
	} else if ( sort_type == "rmsd" ) {
		set_sort_type( SORT_BY_RMSD );
	} else if (sort_type == "hbonds" ) {
		set_sort_type( SORT_BY_HBONDS );
	} else if (sort_type == "cis_peptide_bonds" ) {
		set_sort_type( SORT_BY_CIS_PEPTIDE_BONDS );
	}

	return;
}

/// @brief Set the sort type, by enum.
///
void
HierarchicalHybridJDApplication::set_sort_type(
	HIERARCHICAL_HYBRID_JD_MPI_SORT_TYPE const sort_type
) {
	sort_type_ = sort_type;
}

/// @brief Set the ouput fraction.
/// @details Checks that this is between 0 and 1.
void
HierarchicalHybridJDApplication::set_output_fraction(
	core::Real const &val
) {
	runtime_assert_string_msg( val <= 1.0,
		"Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::set_output_fraction(): The output fraction must be less than or equal to one.");
	runtime_assert_string_msg( val >= 0.0,
		"Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::set_output_fraction(): The output fraction must be greater than or equal to zero.");
	output_fraction_ = val;
}

/// ------------- Methods ----------------------------

/// @brief Reads a FASTA file and returns a string of space-separated full names for the sequence.
/// @details TRIGGERS A READ FROM DISK.
std::string
HierarchicalHybridJDApplication::sequence_from_fasta(
	std::string const &fasta_file
) const {
	utility::vector1< std::string > const fasta_file_contents( core::sequence::read_fasta_file_str( fasta_file ) );
	runtime_assert_string_msg( fasta_file_contents.size() == 1, "Error in HierarchicalHybridJDApplication::sequence_from_fasta():  Found " + std::to_string( fasta_file_contents.size() ) + " sequences in file \"" + fasta_file + "\".  Please provide a file with exactly one sequence." );
	std::string const & fasta_seq( fasta_file_contents[1] );

	// Accumulator for full sequence
	std::stringstream ss;

	// Loop through the string.  Note that strings are zero-based
	for( core::Size i(0); i<fasta_seq.size(); ++i ) {
		if( i == ' ' || i == '\n' || i == '\t' ) continue; //Skip whitespace (oughtn't be necessary, but just in case).
		if( fasta_seq[i] != 'X' && fasta_seq[i] != 'Z' && core::chemical::oneletter_code_specifies_aa( fasta_seq[i]) ) {
			if( !ss.str().empty() ) {
				ss << " ";
			}
			ss << core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( fasta_seq[i] ) );
		} else {
			utility_exit_with_message( "Could not parse FASTA sequence \"" + fasta_seq + "\".  For noncanonical structure prediction, please provide a sequence file containing whitespace-separated basenames using the \"-cyclic_peptide:sequence_file\" flag, instead of a FASTA file." );
		}
	}

	derivedTR_ << "Read sequence \"" << ss.str() << "\" from FASTA file " << fasta_file << "." << std::endl;

	return ss.str();
}

/// @brief Compute SASA metrics for this pose and bundle them in a SASASummary.
/// @details Although derived classes MAY override this, they needn't.  The default behaviour
/// is just to use SASA metrics, which should be pretty general.
HierarchicalHybridJD_SASASummaryOP
HierarchicalHybridJDApplication::generate_sasa_summary(
	core::pose::Pose const &pose,
	HierarchicalHybridJD_JobResultsSummaryCOP jobsummary
) const {
	using namespace core::simple_metrics::metrics;
	using namespace core::scoring::sasa;

	SasaMetric sasa_metric( nullptr, SasaMethodHPMode::ALL_SASA );
	SasaMetric polar_metric( nullptr, SasaMethodHPMode::POLAR_SASA );
	SasaMetric hydrophobic_metric( nullptr, SasaMethodHPMode::HYDROPHOBIC_SASA );

	return utility::pointer::make_shared< HierarchicalHybridJD_SASASummary >(
		jobsummary->originating_node_MPI_rank(),
		jobsummary->jobindex_on_originating_node(),
		sasa_metric.calculate(pose),
		polar_metric.calculate(pose),
		hydrophobic_metric.calculate(pose)
	);
}

/// @brief Check the current time and determine whether it's past the timeout time.
///
bool
HierarchicalHybridJDApplication::halting_condition(
	clock_t const start_time,
	core::Size const timeout
) const {
	return ( static_cast< core::Real >( clock() - start_time ) / static_cast<core::Real>( CLOCKS_PER_SEC ) > static_cast<core::Real>(timeout) );
}


/// @brief Set the number of processes at each level of the communications hierarchy.
/// @details The total_hierarchy_levels, MPI_rank, and MPI_n_procs variables must be set first.  This does
/// some checks to make sure that data_in[1] is 1, that data_in[n] >= data_in[n-1], and that the sum of entries
/// equals the total number of processes.  The total number of hierarchy levels must also match total_hierarchy_levels.
void
HierarchicalHybridJDApplication::set_procs_per_hierarchy_level(
	utility::vector1 < core::Size > const &data_in
) {
	runtime_assert( data_in.size() == total_hierarchy_levels_ );
	runtime_assert( data_in[1] == 1 );
	core::Size sum(0);
	for(core::Size i=1, imax=data_in.size(); i<=imax; ++i) {
		sum+=data_in[i];
		if(i>1) runtime_assert( data_in[i] >= data_in[i-1] );
	}
	runtime_assert(sum == static_cast<core::Size>(MPI_n_procs_));

	//All checks passed at this point.
	procs_per_hierarchy_level_ = data_in;
}

	/// @brief Sets the batch size for this proc.
	/// @details The total_hierarchy_levels, MPI_rank, MPI_n_procs, my_children_, my_parent_, and hierarchy_level_ variables
	/// must be set first.  This does some checks to ensure that the batch sizes of subsequent levels are smaller and previous
	/// levels are larger.
void
HierarchicalHybridJDApplication::set_batchsize(
	utility::vector1< core::Size > const &size_in
) {
	runtime_assert( size_in.size() == total_hierarchy_levels_ - 1 );
	for(core::Size i=1; i<total_hierarchy_levels_; ++i) {
		runtime_assert( size_in[i] > 0 );
		if(i > 1) runtime_assert( size_in[i] <= size_in[i-1] );
	}

	if( hierarchy_level_ == total_hierarchy_levels_ ) {
		batchsize_ = 0;
	} else {
		batchsize_ = size_in[hierarchy_level_];
	}
	derivedTR_.Debug << "Proc " << MPI_rank_ << " set batchsize to " << batchsize_ << "." << std::endl;
}

/// @brief Figure out which processes are my children, and which is my parent.  Also, figure out the level in the
/// communications hierarchy that I'm in.
/// @details This must be done AFTER the MPI_rank_, MPI_n_procs_, total_hierarchy_levels_, and procs_per_hierarchy_level_,
/// variables are set.
void
HierarchicalHybridJDApplication::assign_level_children_and_parent() {
	{ //Scope for checks.
		//Let's check some stuff first:
		runtime_assert(MPI_n_procs_ > 0);
		runtime_assert( total_hierarchy_levels_ > 1 && procs_per_hierarchy_level_.size() == total_hierarchy_levels_ );
		runtime_assert(procs_per_hierarchy_level_[1] == 1);
		core::Size accumulator(0);
		for(core::Size i=1; i<=total_hierarchy_levels_; ++i) {
			accumulator += procs_per_hierarchy_level_[i];
			if(i>1) runtime_assert(procs_per_hierarchy_level_[i] >= procs_per_hierarchy_level_[i-1]);
		}
		runtime_assert(accumulator==static_cast<core::Size>(MPI_n_procs_));
	} //End of scope for checks.

	//OK, now I'm reasonably sure that this class has been set up properly.  Now let's actually do the assignments:
	hierarchy_level_ = 0;
	my_parent_ = 0;
	my_children_.clear();

	//Let's figure out which level we're in, first.  Levels are assigned to processes in sequential order.:
	core::Size curindex(0);
	core::Size curlevel(1);
	utility::vector1 < std::pair < core::Size, core::Size > > level_starts_and_ends;
	for(core::Size i=1; i<=total_hierarchy_levels_; ++i) {
		level_starts_and_ends.push_back( std::pair<core::Size,core::Size>( curindex, curindex + procs_per_hierarchy_level_[i] - 1 ) );

		if(level_starts_and_ends[level_starts_and_ends.size()].first <= static_cast<core::Size>(MPI_rank_) && level_starts_and_ends[level_starts_and_ends.size()].second >= static_cast<core::Size>(MPI_rank_)) {
			hierarchy_level_ = curlevel;
		}
		curindex += procs_per_hierarchy_level_[i];
		++curlevel;
	}

	//Next, let's figure out our parent.
	core::Size const rank_in_level( static_cast<core::Size>(MPI_rank_) - level_starts_and_ends[hierarchy_level_].first + 1); //Declared out of if statement scope for re-use later.
	if( hierarchy_level_ > 1 ) { //Directors have no parent.
		core::Size const parent_level(hierarchy_level_ - 1);
		core::Size const children_per_parent( procs_per_hierarchy_level_[hierarchy_level_] / procs_per_hierarchy_level_[parent_level] );
		core::Size const remainder( procs_per_hierarchy_level_[hierarchy_level_] % procs_per_hierarchy_level_[parent_level] );
		if( rank_in_level <= procs_per_hierarchy_level_[hierarchy_level_] - remainder  ) { //If we're not part of the remainder, divide into even chunks and assign sequentially.
			my_parent_ = ( (rank_in_level - 1) / children_per_parent ) + level_starts_and_ends[parent_level].first;
		} else { // If we're part of the remainder, assign us to parents round-robin style.
			my_parent_ = rank_in_level - children_per_parent*procs_per_hierarchy_level_[parent_level] - 1 + level_starts_and_ends[parent_level].first;
		}
	}

	//Finally, let's figure out our children.
	my_children_.clear();
	if( hierarchy_level_ < total_hierarchy_levels_ ) { //Workers have no children (i.e. nothing below them in the hierarchy -- they do the work themselves).
		core::Size const child_level( hierarchy_level_ + 1 );
		core::Size const children_per_parent( procs_per_hierarchy_level_[child_level] / procs_per_hierarchy_level_[hierarchy_level_] );
		core::Size const remainder( procs_per_hierarchy_level_[child_level] % procs_per_hierarchy_level_[hierarchy_level_] );
		core::Size const first_child ( (rank_in_level-1)*children_per_parent + level_starts_and_ends[child_level].first );
		core::Size const last_child( first_child + children_per_parent - 1 );
		core::Size const extra_child( remainder >= rank_in_level ? children_per_parent*procs_per_hierarchy_level_[hierarchy_level_] + level_starts_and_ends[child_level].first - 1 + rank_in_level : 0 );
		my_children_.reserve( children_per_parent + (extra_child > 0 ? 1 : 0) );
		for(core::Size i=first_child; i<=last_child; ++i ) {
			my_children_.push_back( static_cast<int>(i) );
		}
		if(extra_child > 0 ) my_children_.push_back(static_cast<int>(extra_child));
	}

	if ( derivedTR_.Debug.visible() ) {
		derivedTR_.Debug << "Proc " << MPI_rank_ << " was assigned to communication level " << hierarchy_level_ << ", with parent " << my_parent_;
		if( my_children_.size() > 0 ) {
			derivedTR_.Debug << " and children ";
			for(core::Size i=1, imax=my_children_.size(); i<=imax; ++i) {
				derivedTR_.Debug << my_children_[i];
				if(i<imax) derivedTR_.Debug << ", ";
			}
		} else {
			derivedTR_.Debug << " and no children";
		}
		derivedTR_.Debug << "." << std::endl;
		derivedTR_.Debug.flush();
	}

}

/// @brief Get the amino acid sequence of the peptide we're going to predict; set the sequence_ private member variable.
/// @details The director reads this from disk and broadcasts it to all other nodes.  This function can be called from any node;
/// it figures out which behaviour it should be performing.
void
HierarchicalHybridJDApplication::get_sequence() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( i_am_director() ) {
		runtime_assert_string_msg(
			option[basic::options::OptionKeys::helical_bundle_predict::sequence_file].user() || option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user() || option[basic::options::OptionKeys::in::file::fasta].user(),
			"Error in get_sequence() in HierarchicalHybridJDApplication: No sequence file was specified.  A sequence file must be provided either with \"-cyclic_peptide:sequence_file <filename>\" or with \"-in:file:fasta <filename>\"."
		);
		if( option[basic::options::OptionKeys::helical_bundle_predict::sequence_file].user() ) {
			runtime_assert_string_msg( !option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user() && !option[basic::options::OptionKeys::in::file::fasta].user(), "Error in get_sequence() in HierarchicalHybridJDApplication: Only one of \"-helical_bundle_predict:sequence_file\" or \"-cyclic_peptide:sequence_file\" or \"-in:file:fasta\" should be provided." );
			sequence_=utility::file_contents( option[basic::options::OptionKeys::helical_bundle_predict::sequence_file]() );
			utility::trim( sequence_, " \n\t" );
		} else if( option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user() ) {
			runtime_assert_string_msg( !option[basic::options::OptionKeys::in::file::fasta].user(), "Error in get_sequence() in HierarchicalHybridJDApplication: Only one of \"-helical_bundle_predict:sequence_file\" or \"-cyclic_peptide:sequence_file\" or \"-in:file:fasta\" should be provided." );
			sequence_=utility::file_contents( option[basic::options::OptionKeys::cyclic_peptide::sequence_file]() );
			utility::trim( sequence_, " \n\t" );
		} else if( option[basic::options::OptionKeys::in::file::fasta].user() ) {
			runtime_assert_string_msg( option[basic::options::OptionKeys::in::file::fasta]().size() == 1, "Error in get_sequence() in HierarchicalHybridJDApplication: Only one FASTA file may be provided with the \"-in:file:fasta\" flag." );
			sequence_ = sequence_from_fasta( option[basic::options::OptionKeys::in::file::fasta]()[1] );
		}
	}

	int stringlen( i_am_director() ? sequence_.size() + 1 : 0);
	MPI_Bcast( &stringlen, 1, MPI_INT, 0, MPI_COMM_WORLD); //Broadcast the length of the sequence string.

	utility::vector0< char > charseq( stringlen );
	if( i_am_director() ) sprintf( charseq.data(), "%s", sequence_.c_str() );

	MPI_Bcast( charseq.data(), stringlen, MPI_CHAR, 0, MPI_COMM_WORLD );
	if( i_am_director() ) {
		derivedTR_ << "Director read sequence from disk and and sent \"" << sequence_ << "\" to all other nodes." << std::endl;
		derivedTR_.flush();
	} else {
		sequence_ = std::string(charseq.data());
		derivedTR_.Debug << "Received sequence \"" << sequence_ << "\" in broadcast from director." << std::endl;
		derivedTR_.Debug.flush();
	}
}

/// @brief Get the native structure of the peptide we're going to predict (if the user has specified one with the -in:file:native flag).
/// @details The director reads this from disk and broadcasts it to all other nodes.  This function should be called from all nodes;
/// it figures out which behaviour it should be performing.
void
HierarchicalHybridJDApplication::get_native() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( !option[basic::options::OptionKeys::in::file::native].user() ) return; // Do nothing if no native state has been specified.

	if( i_am_director() ) { //The director reads the native state from disk.
		std::string natfile( option[basic::options::OptionKeys::in::file::native]() );
		derivedTR_ << "Director reading native pose " << natfile << " from disk." << std::endl;
		core::pose::PoseOP native( utility::pointer::make_shared< core::pose::Pose >() );
		core::import_pose::pose_from_file(*native, natfile, core::import_pose::PDB_file);
		core::io::silent::SilentFileOptions opts;
		opts.in_fullatom(true);
		core::io::silent::SilentStructOP native_ss( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( "binary", opts ) );
		native_ss->fill_struct(*native, "native");

		broadcast_silent_struct_from_this_node( native_ss );

		native_ = native; //Store the pose (though maybe the director should discard it at this point to free memory).
		derivedTR_ << "Read " << native_->size() << "-residue pose from disk and broadcasted it to all other nodes." << std::endl;
	} else { //Other nodes do this
		native_ = receive_broadcast_silent_struct_and_build_pose();
		if(derivedTR_.Debug.visible()) {
			derivedTR_.Debug << "Process " << MPI_rank_ << " received a " << native_->size() << "-residue pose from the director's broadcast, with sequence \"";
			for(core::Size i=1, imax=native_->size(); i<=imax; ++i) {
				derivedTR_.Debug << native_->residue(i).name3();
				if(i<imax) derivedTR_.Debug << " ";
			}
			derivedTR_.Debug << "\"." << std::endl;
			derivedTR_.Debug.flush();
		}
	}
}

/// @brief Given a map of indicies to lists of residue names, broadcast it to all MPI ranks.
void
HierarchicalHybridJDApplication::broadcast_res_list(
	std::map< core::Size, utility::vector1 < std::string > > &res_list
) const {
	if( i_am_director() ) { //The director converts the map into a string and sends it.
		std::stringstream reslist_string;
		bool first_map_entry(true);
		for (std::map<core::Size, utility::vector1 < std::string > >::const_iterator it=res_list.begin(); it!=res_list.end(); ++it) {
			if(first_map_entry) { first_map_entry=false; }
			else { reslist_string << " "; }
			reslist_string << it->first << " " << it->second.size() << " ";
			for(core::Size i=1, imax=it->second.size(); i<=imax; ++i) {
				reslist_string << it->second[i];
				if( i < imax ) reslist_string << " ";
			}
		}
		std::string reslist_str( reslist_string.str() );
		broadcast_string_from_director( reslist_str );
		if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Reslist string broadcast from director:\n" << reslist_str << std::endl;

	} else {  //Everyone else receives the string and converts it into a map.
		res_list.clear();
		std::string conversion_string( "" );
		broadcast_string_from_director( conversion_string );
		if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Reslist string received by proc " << MPI_rank_ << ":\n" << conversion_string << std::endl;
		std::stringstream received_string( conversion_string );

		while( !received_string.eof() ) {
			core::Size index, entries;
			received_string >> index >> entries;
			res_list[index].reserve(entries);
			if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Parsing " << entries << " entries." << std::endl;
			for(core::Size i=1; i<=entries; ++i ) {
				std::string curname;
				received_string >> curname;
				res_list[index].push_back(curname);
				if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Added " << curname << "." << std::endl;
			}
		}
		if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Finished splitting reslist on proc " << MPI_rank_ << "." << std::endl;
	}
}


/// @brief The director sends a message to everyone letting them know it's time to start.
/// @details Following the go signal, workers send requests for jobs upward.
void
HierarchicalHybridJDApplication::go_signal() const {
	char mybuf('G'); //A dummy piece of information to send.
	MPI_Bcast( &mybuf, 1, MPI_CHAR, 0, MPI_COMM_WORLD );
}

/// @brief The director sends a message to everyone letting them know it's time to stop.
/// @details Following the stop signal, HierarchicalHybridJDApplication::run() terminates.
void
HierarchicalHybridJDApplication::stop_signal() const {
	char mybuf('S'); //A dummy piece of information to send.
	MPI_Bcast( &mybuf, 1, MPI_CHAR, 0, MPI_COMM_WORLD );
}

/// @brief Any node can wait for a node, above or below in the hierarchy, to send it some sort of request.
/// @details Only messags with tag GENERAL_REQUEST.
/// @param[out] requesting_node The node from which the request came.
/// @param[out] message The type of request received.
void
HierarchicalHybridJDApplication::wait_for_request(
	int &requesting_node,
	HIERARCHICAL_MPI_COMMUNICATION_TYPE &message
) const {
	int buf(0);
	MPI_Status status;
	MPI_Recv( &buf, 1, MPI_INT, MPI_ANY_SOURCE, static_cast<int>(GENERAL_REQUEST), MPI_COMM_WORLD, &status );
	requesting_node = status.MPI_SOURCE;
	message = static_cast<HIERARCHICAL_MPI_COMMUNICATION_TYPE>(buf);
}

/// @brief Send a signal to stop job distribution to a set of intermediate managers.
///
void
HierarchicalHybridJDApplication::send_halt_signal(
	utility::vector1 < int > const &ranks_to_target
) const {
	int buf(0);
	for(core::Size i=1, imax=ranks_to_target.size(); i<=imax; ++i) {
		MPI_Send( &buf, 1, MPI_INT, ranks_to_target[i], static_cast<int>(HALT_SIGNAL), MPI_COMM_WORLD );
		if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Proc " << MPI_rank_ << " sent halt signal to proc " << ranks_to_target[i] << "." << std::endl;
	}
}


/// @brief Any node can send some sort of request to a specific node in the hierarchy (parent or child).
/// @details Sends messags with tag GENERAL_REQUEST.
/// @param[in] target_node The node to which we're sending the request.
/// @param[in] message The type of request to send.
void
HierarchicalHybridJDApplication::send_request(
	int const target_node,
	HIERARCHICAL_MPI_COMMUNICATION_TYPE const message
) const {
	int buf( static_cast<int>(message) );
	MPI_Send( &buf, 1, MPI_INT, target_node, static_cast<int>(GENERAL_REQUEST), MPI_COMM_WORLD );
}

/// @brief Send an integer to a child indicating that the child is now holding this number of jobs.
/// @details Sends messags with tag NEW_JOBS_DOWNWARD.  Sends zero if no jobs remain to send.
/// @param[in] target_node The rank of the node to which we're sending the request.
/// @param[in,out] structures_remaining_to_send The number of structures remaining in the send buffer.  Decremented by this operation to a minimum of zero.
/// @param[in] number_to_send How many jobs should I send down?
void
HierarchicalHybridJDApplication::send_new_jobs_downward(
	int const target_node,
	core::Size &structures_remaining_to_send,
	core::Size const number_to_send
) const {
	unsigned long buf( static_cast<unsigned long>( ( structures_remaining_to_send < number_to_send ? structures_remaining_to_send : number_to_send) ) );

	MPI_Send( &buf, 1, MPI_UNSIGNED_LONG, target_node, static_cast<int>(NEW_JOBS_DOWNWARD), MPI_COMM_WORLD );

	if( structures_remaining_to_send > number_to_send ) { structures_remaining_to_send -= number_to_send; }
	else { structures_remaining_to_send = 0; }
}

/// @brief Receive an integer from my parent indicating that I should add N jobs to my set to carry out or to pass to my children.
/// @details Receives message with tag NEW_JOBS_DOWNARD.  Receives zero if no jobs remain in parent to send.
/// @param[in,out] njobs The number of jobs held on this node that are to be done.  Incremented by this function with however many are received from above.
void
HierarchicalHybridJDApplication::receive_njobs_from_above(
	core::Size &njobs
) const {
	unsigned long buf(0);
	MPI_Status status;
	MPI_Recv( &buf, 1, MPI_UNSIGNED_LONG, my_parent_, static_cast<int>(NEW_JOBS_DOWNWARD), MPI_COMM_WORLD, &status);
	njobs += static_cast<core::Size>(buf);
}

/// @brief Non-originating nodes must call this when the originating node calls broadcast_silent_struct_from_this_node.
/// @details This will build a pose and return an owning pointer to it.
/// @note If skip_pose_build is true, we don't bother to convert the structure into a pose; we just participate in the
/// broadcast.  In that case, this function returns nullptr.
core::pose::PoseCOP
HierarchicalHybridJDApplication::receive_broadcast_silent_struct_and_build_pose(
	int node_to_receive_from /*=0*/,
	bool const skip_pose_build /*=false*/
) const {
	int strlen(0);
	MPI_Bcast( &strlen, 1, MPI_INT, node_to_receive_from, MPI_COMM_WORLD); //Receive the length of the string.
	utility::vector0< char > inchar( strlen );
	MPI_Bcast( inchar.data(), strlen, MPI_CHAR, node_to_receive_from, MPI_COMM_WORLD); //Receive the string.

	if( skip_pose_build ) return nullptr;

	std::string instring( inchar.data() );

	std::istringstream in_ss( instring );
	utility::vector1 < std::string > tagvector;
	core::io::silent::SilentFileOptions opts;
	opts.in_fullatom(true);
	core::io::silent::SilentFileData silentfiledata( opts );
	silentfiledata.read_stream( in_ss, tagvector, true, "suppress_bitflip" );
	core::pose::PoseOP mypose( utility::pointer::make_shared< core::pose::Pose >());
	core::chemical::ResidueTypeSet const & restypeset( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ) ) );
	silentfiledata[silentfiledata.tags()[1]]->fill_pose(*mypose, restypeset);

	return core::pose::PoseCOP(mypose);
}

/// @brief Convert a silent struct into a character string and send it to a node.
/// @details Intended to be used with receive_silent_struct() to allow another node to receive the broadcast.
/// Message is tagged with SILENT_STRUCT_TRANSMISSION
void
HierarchicalHybridJDApplication::send_silent_structs(
	utility::vector1 < core::io::silent::SilentStructOP > const &ss_vect,
	int const target_node
) const {
	core::io::silent::SilentFileOptions opts;
	opts.in_fullatom(true);
	core::io::silent::SilentFileData silentfiledata( opts );
	std::stringbuf sb;
	std::ostream outstream(&sb);

	for(core::Size i=1, imax=ss_vect.size(); i<=imax; ++i) {
		silentfiledata._write_silent_struct( *(ss_vect[i]), outstream );
	}
	std::string outstring( sb.str() );

	int strlen( outstring.length() + 1 ); //The plus one is for the /0 terminal character
	utility::vector0< char > outchar( strlen );
	sprintf( outchar.data(), "%s", outstring.c_str() );

	MPI_Send( &strlen, 1, MPI_INT, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the length of the string.
	MPI_Send( outchar.data(), strlen, MPI_CHAR, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the string.
}

/// @brief Receive a transmitted set of poses (as silent strings).
/// @details Appends received information to results_string.
void
HierarchicalHybridJDApplication::receive_pose_batch_as_string(
	int const originating_node,
	std::string &results_string
) const {
	MPI_Status status;
	int strlen( 0 );
	MPI_Recv( &strlen, 1, MPI_INT, originating_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD, &status); //Send the length of the string.

	utility::vector0<char> charvect( strlen );
	MPI_Recv( charvect.data(), strlen, MPI_CHAR, originating_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD, &status); //Send the length of the string.
	std::string const received_string( charvect.data() );
	results_string += received_string;
}

/// @brief Receive the number of jobs attempted by all of the nodes below a child node, and add this to the total.
///
void
HierarchicalHybridJDApplication::receive_jobs_attempted_count(
	core::Size &total_jobs_attempted,
	int const originating_node
) const {
	MPI_Status status;
	unsigned long sizebuf(0);
	MPI_Recv( &sizebuf, 1, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(JOBS_ATTEMPTED_COUNT_UPWARD), MPI_COMM_WORLD, &status );
	total_jobs_attempted += static_cast<core::Size>(sizebuf);
}

/// @brief Send the number of jobs attempted by this nodes or all of the nodes below this node.
///
void
HierarchicalHybridJDApplication::send_jobs_attempted_count(
	core::Size const total_jobs_attempted,
	int const target_node
) const {
	unsigned long sizebuf( static_cast<unsigned long>(total_jobs_attempted) );
	MPI_Send( &sizebuf, 1, MPI_UNSIGNED_LONG, target_node, static_cast<int>(JOBS_ATTEMPTED_COUNT_UPWARD), MPI_COMM_WORLD);
}


/// @brief Recieve a sorted list of job summaries, and merge them with an existing sorted list to make a combined sorted list.
/// @details To be used in conjunction with send_job_summaries().  Sending and receiving procs must send messages to synchronize, first.
void
HierarchicalHybridJDApplication::receive_and_sort_job_summaries(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > &original_summary_list,
	int const originating_node,
	bool const append_to_handler_list
) const {
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > summary_list;
	MPI_Status status;

	unsigned long sizebuf(0);
	MPI_Recv( &sizebuf, 1, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //First, receive the size of the list.
	summary_list.reserve( sizebuf );

	if(sizebuf == 0) return;

	//Buffers for data that we'll send:
	utility::vector0 < int > procbuf( sizebuf );
	utility::vector0 < unsigned long > jobindexbuf( sizebuf );
	utility::vector0 < double > energiesbuf( sizebuf );
	utility::vector0 < double > rmsdbuf( sizebuf );
	utility::vector0 < unsigned long > hbondsbuf( sizebuf );
	utility::vector0 < unsigned long > cispepbondbuf( sizebuf );

	MPI_Recv( procbuf.data(), sizebuf, MPI_INT, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process indices.
	MPI_Recv( jobindexbuf.data(), sizebuf, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process job indices.
	MPI_Recv( energiesbuf.data(), sizebuf, MPI_DOUBLE, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the pose energies.
	MPI_Recv( rmsdbuf.data(), sizebuf, MPI_DOUBLE, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the RMSD values.
	MPI_Recv( hbondsbuf.data(), sizebuf, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the hydrogen bond counts.
	MPI_Recv( cispepbondbuf.data(), sizebuf, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the cis-peptide bond counts.

	for(core::Size i=1, imax=static_cast<core::Size>(sizebuf); i<=imax; ++i) {
		summary_list.push_back( utility::pointer::make_shared< HierarchicalHybridJD_JobResultsSummary >() );

		summary_list[i]->set_originating_node_MPI_rank( procbuf[i-1] );
		summary_list[i]->set_jobindex_on_originating_node( static_cast<core::Size>(jobindexbuf[i-1]) );
		summary_list[i]->set_pose_energy( static_cast<core::Real>(energiesbuf[i-1]) );
		summary_list[i]->set_rmsd( static_cast<core::Real>(rmsdbuf[i-1]) );
		summary_list[i]->set_hbonds( static_cast<core::Size>(hbondsbuf[i-1]) );
		summary_list[i]->set_cis_peptide_bonds( static_cast<core::Size>(cispepbondbuf[i-1]) );
	}

	debug_assert(summary_list.size() == static_cast<core::Size>(sizebuf));

	//Receive the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1, imax=summary_list.size(); i<=imax; ++i) {
		int nodelist_size(0);
		MPI_Recv( &nodelist_size, 1, MPI_INT, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );

		if(nodelist_size != 0) {
			utility::vector0< int > nodelist( nodelist_size );
			MPI_Recv( nodelist.data(), nodelist_size, MPI_INT, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );
			for( core::Size j=1; j<=static_cast<core::Size>(nodelist_size); ++j ) {
				summary_list[i]->add_MPI_rank_handling_message( nodelist[j-1] );
			}
		}
		if( append_to_handler_list ) {
			summary_list[i]->add_MPI_rank_handling_message( originating_node ); //Add the originating node to the list of MPI ranks that handled this message.
		}

	}

	mergesort_jobsummaries_list( original_summary_list, summary_list, sort_type_ );
}

/// @brief Send a list of job summaries.
/// @details To be used in conjunction with receive_and_sort_job_summaries().  Sending and receiving procs must send messages to synchronize, first.
void
HierarchicalHybridJDApplication::send_job_summaries(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_list,
	int const target_node
) const {
	unsigned long sizebuf( summary_list.size() );
	MPI_Send( &sizebuf, 1, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //First, send the size of the list.

	if(sizebuf == 0) return;

	//Buffers for data that we'll send:
	utility::vector0 < int > procbuf( sizebuf );
	utility::vector0 < unsigned long > jobindexbuf( sizebuf );
	utility::vector0 < double > energiesbuf( sizebuf );
	utility::vector0 < double > rmsdbuf( sizebuf );
	utility::vector0 < unsigned long > hbondsbuf( sizebuf );
	utility::vector0 < unsigned long > cispepbondbuf( sizebuf );

	for(core::Size i=1, imax=summary_list.size(); i<=imax; ++i) {
		procbuf[i-1] = summary_list[i]->originating_node_MPI_rank();
		jobindexbuf[i-1] = static_cast<unsigned long>( summary_list[i]->jobindex_on_originating_node() );
		energiesbuf[i-1] = static_cast<double>( summary_list[i]->pose_energy() );
		rmsdbuf[i-1] = static_cast<double>( summary_list[i]->rmsd() );
		hbondsbuf[i-1] = static_cast<unsigned long>( summary_list[i]->hbonds() );
		cispepbondbuf[i-1] = static_cast<unsigned long>( summary_list[i]->cis_peptide_bonds() );
	}

	MPI_Send( procbuf.data(), sizebuf, MPI_INT, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process indices.
	MPI_Send( jobindexbuf.data(), sizebuf, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process job indices.
	MPI_Send( energiesbuf.data(), sizebuf, MPI_DOUBLE, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the pose energies.
	MPI_Send( rmsdbuf.data(), sizebuf, MPI_DOUBLE, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the RMSD values.
	MPI_Send( hbondsbuf.data(), sizebuf, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the hydrogen bond counts.
	MPI_Send( cispepbondbuf.data(), sizebuf, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the cis-peptide bond counts.

	//Send the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1, imax=summary_list.size(); i<=imax; ++i) {
		int nodelist_size(summary_list[i]->MPI_ranks_handling_message().size());
		MPI_Send( &nodelist_size, 1, MPI_INT, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
		if(nodelist_size == 0) continue;

		utility::vector0< int > nodelist( summary_list[i]->MPI_ranks_handling_message().size() );
		for( core::Size j=1, jmax=summary_list[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j ) {
			nodelist[j-1] = summary_list[i]->MPI_ranks_handling_message()[j];
		}
		MPI_Send( nodelist.data(), summary_list[i]->MPI_ranks_handling_message().size(), MPI_INT, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
	}
}

/// @brief Given a short list of job summaries, split the list by the index of the node that I'd have to send the request to, and send requests for full poses down the hierarchy.
/// @details Throws an error if any of the jobs in the list cannot be reached by sending a request by way of my_children_.  To be used with receive_pose_requests_from_above().
/// @param[in] summary_shortlist The list of job summaries to split up and send downard.
/// @param[out] children_receiving_nonzero_requests The number of children receiving a request for at least one pose.
void
HierarchicalHybridJDApplication::request_poses_from_below(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_shortlist,
	core::Size &children_receiving_nonzero_requests
) const {
	utility::vector1< utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > > split_shortlists; //Storage for splitting the shortlist into one for each child node.
	core::Size const nchildren( my_children_.size() );
	split_shortlists.resize( nchildren );

	for(core::Size i=1, imax=summary_shortlist.size(); i<=imax; ++i) {
		core::Size child_index(0); //0 will be the failure signal.
		for(core::Size j=1; j<=nchildren; ++j) {
			core::Size const rank_to_check( summary_shortlist[i]->MPI_ranks_handling_message().size() - hierarchy_level_ + 1 );
			if( summary_shortlist[i]->MPI_ranks_handling_message()[ rank_to_check ] == my_children_[j] ) {
				child_index=j;
				break;
			}
		}
		//if( child_index == 0 ) TR.Error << "Node " << MPI_rank_ << " was unable to assign a child given the jobs list.  My children are: " << my_children_ << ".  My hierarchy level is " << hierarchy_level_ << "." << std::endl; //FOR DEBUGGING ONLY
		runtime_assert_string_msg( child_index > 0, "Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::request_poses_from_below(): Function was passed a set of jobs that included jobs not handled by child nodes." );
		split_shortlists[child_index].push_back( summary_shortlist[i] );
	}

	children_receiving_nonzero_requests=0;
	for(core::Size j=1; j<=nchildren; ++j) {
		if(split_shortlists[j].size() > 0) ++children_receiving_nonzero_requests;
		send_job_summaries( split_shortlists[j], my_children_[j] );
	}
}

/// @brief Recieve a short list of job summaries from above.
/// @param[out] summary_shortlist The job summaries list, cleared and populated by this function.
void
HierarchicalHybridJDApplication::receive_pose_requests_from_above(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > &summary_shortlist
) const {
	summary_shortlist.clear();
	receive_and_sort_job_summaries( summary_shortlist, my_parent_, false /*No need to list all the nodes that we pass through this time*/ );
}

/// @brief The director is sending out a request that the worker that produced the top pose broadcast
/// that pose to all other nodes.  All nodes participate in the broadcast.  At the end of this operation,
/// all nodes know:
/// - The node index that produced the top structure.
/// - The index of the top structure on that node index.
/// - The top structure (the pose), as a binary silent structure.
/// They can then compare the top structure to all of their structures and compute an RMSD for each.
/// @note If broadcast_no_poses_found is true, the director announces to all processes that there is no best
/// pose with which to compare.
bool
HierarchicalHybridJDApplication::top_pose_bcast(
	HierarchicalHybridJD_JobResultsSummaryOP top_summary /*=nullptr*/,
	core::io::silent::SilentStructOP top_pose_silentstruct /*=nullptr*/,
	bool const broadcast_no_poses_found/*=false*/
) const {
	// 1.  The director broadcasts to everyone else the node index and job index for the best sample:
	if( i_am_director() ) {
		runtime_assert( broadcast_no_poses_found || top_summary != nullptr );
	}
	int there_were_no_poses( broadcast_no_poses_found ? 1 : 0 );
	MPI_Bcast( &there_were_no_poses, 1, MPI_INT, 0 /*Director broadcasts*/, MPI_COMM_WORLD );
	if( there_were_no_poses==1 ) return true;

	int node_that_produced_best( top_summary == nullptr ? 0 : static_cast< int >( top_summary->originating_node_MPI_rank() ) );
	unsigned long jobindex_on_node_that_produced_best( top_summary == nullptr ? 0 : static_cast< unsigned long >( top_summary->jobindex_on_originating_node() ) );
	MPI_Bcast( &node_that_produced_best, 1, MPI_INT, 0 /*Director broadcasts*/, MPI_COMM_WORLD  );
	MPI_Bcast( &jobindex_on_node_that_produced_best, 1, MPI_UNSIGNED_LONG, 0 /*Director broadcasts*/, MPI_COMM_WORLD  );

	// 2.  The node that produced the top structure transmits it to all others.
	std::string silentstring;
	if( MPI_rank_ == node_that_produced_best ) {
		runtime_assert( top_pose_silentstruct != nullptr );
		core::io::silent::SilentFileOptions opts;
		opts.in_fullatom(true);
		core::io::silent::SilentFileData silentfiledata( opts );
		std::stringbuf sb;
		std::ostream outstream(&sb);

		silentfiledata._write_silent_struct( (*top_pose_silentstruct), outstream );
		silentstring = sb.str();
	}
	broadcast_string_from_node( silentstring, node_that_produced_best );

	// 3.  All other worker nodes reconstruct the silent structure.
	if( MPI_rank_ != node_that_produced_best && top_pose_silentstruct != nullptr ) {
		std::istringstream in_ss( silentstring );
		utility::vector1 < std::string > tagvector;
		core::io::silent::SilentFileOptions opts;
		opts.in_fullatom(true);
		core::io::silent::SilentFileData silentfiledata2( opts );
		silentfiledata2.read_stream( in_ss, tagvector, true, "suppress_bitflip" );
		core::pose::PoseOP mypose( utility::pointer::make_shared< core::pose::Pose >());
		core::io::silent::BinarySilentStruct & top_pose_binarysilentstruct( dynamic_cast< core::io::silent::BinarySilentStruct & >( *top_pose_silentstruct ) );
		top_pose_binarysilentstruct = dynamic_cast< core::io::silent::BinarySilentStruct const & >( *( silentfiledata2[silentfiledata2.tags()[1]] ) );
		//(*top_pose_silentstruct) = ( *( silentfiledata2[silentfiledata2.tags()[1]] ) );
	}

	return false;
}

/// @brief Given a string on the director node, send it to all nodes.
/// @details The "mystring" string is the input on the director node, and the output on all other nodes.
void
HierarchicalHybridJDApplication::broadcast_string_from_director(
	std::string &mystring
) const {
	broadcast_string_from_node( mystring, 0 );
}

/// @brief Given a string on a given node, send it to all nodes.
/// @details The "mystring" string is the input on the originating node, and the output on all other nodes.
/// @note Protected, not private.
void
HierarchicalHybridJDApplication::broadcast_string_from_node(
	std::string &mystring,
	int const originating_node_index
) const {
	if( MPI_rank_ == originating_node_index ) {
		//if(TR.Debug.visible()) derivedTR_.Debug << "Broadcasting the following string from the director:\n" << mystring << std::endl;
		unsigned long charsize( mystring.size() + 1 );
		utility::vector0< char > bcastchar( charsize ); // Plus one for null terminator.
		for(core::Size i=0; i < static_cast<core::Size>(charsize); ++i) bcastchar[i] = mystring.c_str()[i]; //Copy the string to the char array for broadcast.
		MPI_Bcast( &charsize, 1, MPI_UNSIGNED_LONG, originating_node_index, MPI_COMM_WORLD ); //Broadcast the length of the string.
		MPI_Bcast( bcastchar.data(), charsize, MPI_CHAR, originating_node_index, MPI_COMM_WORLD ); //Broadcast the string.
		//if(TR.Debug.visible()) derivedTR_.Debug << "String broadcast from director complete." << std::endl;
	} else {
		unsigned long charsize(0);
		MPI_Bcast( &charsize, 1, MPI_UNSIGNED_LONG, originating_node_index, MPI_COMM_WORLD ); //Receive the length of the string.
		utility::vector0 < char > bcastchar( charsize );
		MPI_Bcast( bcastchar.data(), charsize, MPI_CHAR, originating_node_index, MPI_COMM_WORLD ); //Receive the string.
		mystring.clear();
		mystring = bcastchar.data();
	}
}

/// @brief Given a vector of RMSDs to best summaries, send them up the hierarchy.
/// @details Made for use with receive_and_sort_all_rmsd_to_best_summaries().
void
HierarchicalHybridJDApplication::send_rmsds_to_best_summaries_upward(
	utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > const &rmsds_to_best_summaries,
	int const receiving_node
) const {
	send_request( receiving_node, OFFER_NEW_RMSD_TO_BEST_SUMMARY_BATCH_UPWARD );
	unsigned long sizebuf( rmsds_to_best_summaries.size() );
	MPI_Send( &sizebuf, 1, MPI_UNSIGNED_LONG, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //First, send the size of the list.

	if(sizebuf == 0) return;

	//Buffers for data that we'll send:
	utility::vector0<int> procbuf( sizebuf );
	utility::vector0<unsigned long> jobindexbuf( sizebuf );
	utility::vector0<double> rmsdbuf( sizebuf );

	for(core::Size i=1; i<=static_cast<core::Size>(sizebuf); ++i) {
		procbuf[i-1] = rmsds_to_best_summaries[i]->originating_node_MPI_rank();
		jobindexbuf[i-1] = static_cast<unsigned long>( rmsds_to_best_summaries[i]->jobindex_on_originating_node() );
		rmsdbuf[i-1] = static_cast<double>( rmsds_to_best_summaries[i]->rmsd_to_best() );
	}

	MPI_Send( procbuf.data(), sizebuf, MPI_INT, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process indices.
	MPI_Send( jobindexbuf.data(), sizebuf, MPI_UNSIGNED_LONG, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process job indices.
	MPI_Send( rmsdbuf.data(), sizebuf, MPI_DOUBLE, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the RMSD values.

	derivedTR_.Debug << "Sent " << sizebuf << " RMSD-to-best summaries to node " << receiving_node << "." << std::endl;

	//Send the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1; i<=static_cast<core::Size>(sizebuf); ++i) {
		int nodelist_size(rmsds_to_best_summaries[i]->MPI_ranks_handling_message().size());
		MPI_Send( &nodelist_size, 1, MPI_INT, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
		if(nodelist_size == 0) continue;

		utility::vector0< int > nodelist( rmsds_to_best_summaries[i]->MPI_ranks_handling_message().size() );
		for( core::Size j=1, jmax=rmsds_to_best_summaries[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j ) {
			nodelist[j-1] = rmsds_to_best_summaries[i]->MPI_ranks_handling_message()[j];
		}
		MPI_Send( nodelist.data(), rmsds_to_best_summaries[i]->MPI_ranks_handling_message().size(), MPI_INT, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
	}
}

/// @brief From one of the nodes below me in the next level down, receive one vector of RMSD-to-best summaries.
/// @details Made for use with send_rmsds_to_best_summaries_upward().  THe new_rmsd_to_best_summaries vector is cleared
/// and populated by this operation.
void
HierarchicalHybridJDApplication::receive_rmsd_to_best_summaries_from_below(
	utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > &new_rmsd_to_best_summaries,
	int const requesting_node,
	bool const append_to_handler_list
) const {
	runtime_assert( new_rmsd_to_best_summaries.empty() );

	MPI_Status status;
	unsigned long sizebuf(0);
	MPI_Recv( &sizebuf, 1, MPI_UNSIGNED_LONG, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //First, receive the size of the list.
	new_rmsd_to_best_summaries.reserve( sizebuf );

	if(sizebuf == 0) return;

	//Buffers for data that we'll send:
	utility::vector0<int> procbuf( sizebuf );
	utility::vector0<unsigned long> jobindexbuf( sizebuf );
	utility::vector0<double> rmsdbuf( sizebuf );

	MPI_Recv( procbuf.data(), sizebuf, MPI_INT, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process indices.
	MPI_Recv( jobindexbuf.data(), sizebuf, MPI_UNSIGNED_LONG, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process job indices.
	MPI_Recv( rmsdbuf.data(), sizebuf, MPI_DOUBLE, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the RMSD values.

	for(core::Size i=1, imax=static_cast<core::Size>(sizebuf); i<=imax; ++i) {
		new_rmsd_to_best_summaries.push_back(
			utility::pointer::make_shared< HierarchicalHybridJD_RMSDToBestSummary >( procbuf[i-1], jobindexbuf[i-1], rmsdbuf[i-1] )
		);
	}

	runtime_assert(new_rmsd_to_best_summaries.size() == static_cast<core::Size>(sizebuf));

	derivedTR_.Debug << "Received " << sizebuf << " RMSD-to-best summaries from node " << requesting_node << "." << std::endl;

	//Receive the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1, imax=new_rmsd_to_best_summaries.size(); i<=imax; ++i) {
		int nodelist_size(0);
		MPI_Recv( &nodelist_size, 1, MPI_INT, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );

		if(nodelist_size != 0) {
			utility::vector0<int> nodelist( nodelist_size );
			MPI_Recv( nodelist.data(), nodelist_size, MPI_INT, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );
			for( core::Size j=1; j<=static_cast<core::Size>(nodelist_size); ++j ) {
				new_rmsd_to_best_summaries[i]->add_MPI_rank_handling_message( nodelist[j-1] );
			}
		}
		if( append_to_handler_list ) {
			new_rmsd_to_best_summaries[i]->add_MPI_rank_handling_message( requesting_node ); //Add the originating node to the list of MPI ranks that handled this message.
		}

	}
}

/// @brief From all of the nodes below me in the next level down, receive RMSD-to-best summaries.  Then sort them
/// based on the already-received job summaries.
/// @details Made for use with send_rmsds_to_best_summaries_upward().  THe rmsds_to_best_summaries vector is cleared
/// and populated by this operation.
void
HierarchicalHybridJDApplication::receive_and_sort_all_rmsd_to_best_summaries (
	utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > &rmsds_to_best_summaries,
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &job_summary_list,
	bool const append_to_handler_list
) const {
	runtime_assert( rmsds_to_best_summaries.empty() );
	if( job_summary_list.size() == 0 ) return;
	rmsds_to_best_summaries.resize( job_summary_list.size(), nullptr );

	core::Size n_summaries_received(0);
	core::Size const n_children( my_children_.size() );

	do {
		int requesting_node(0);
		HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_below(NULL_MESSAGE);

		wait_for_request( requesting_node, message_from_below );
		runtime_assert_string_msg( message_from_below == OFFER_NEW_RMSD_TO_BEST_SUMMARY_BATCH_UPWARD, "Error in HierarchicalHybridJDApplication::receive_and_sort_all_rmsd_to_best_summaries(): Children nodes of node " + std::to_string( MPI_rank_ ) + " should be transmitting RMSD-to-best summaries at this point, but I received something else from node " + std::to_string( requesting_node ) + "!" );

		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > new_rmsd_to_best_summaries;
		receive_rmsd_to_best_summaries_from_below( new_rmsd_to_best_summaries, requesting_node, append_to_handler_list );
		++n_summaries_received;

		// Sort the new list into the old:
		core::Size counter(1);
		for( core::Size i(1), imax( job_summary_list.size() ); i<=imax; ++i ) { //Loop through the job summary list, placing each summary.
			//Note that the received list is already sorted in the same order as the job_summary_list, except that the job_summary_list is longer (has more elements).
			if( counter > new_rmsd_to_best_summaries.size() ) break; //Stop when all have been placed.

			if( job_summary_list[i]->originating_node_MPI_rank() == new_rmsd_to_best_summaries[counter]->originating_node_MPI_rank() ) {
				runtime_assert( job_summary_list[i]->jobindex_on_originating_node() == new_rmsd_to_best_summaries[counter]->jobindex_on_originating_node() );
				runtime_assert( rmsds_to_best_summaries[i] == nullptr );
				rmsds_to_best_summaries[i] = new_rmsd_to_best_summaries[counter];
				++counter;
			}
		}

		//Sanity checks
		runtime_assert( counter == new_rmsd_to_best_summaries.size() + 1 );

	} while( n_summaries_received < n_children );  //Keep looping until all children have sent summaries.

	//Another sanity check:
	for( core::Size i(1), imax(rmsds_to_best_summaries.size()); i<=imax; ++i ) {
		runtime_assert( rmsds_to_best_summaries[i] != nullptr );
	}
}

/// @brief Given a vector of SASA summaries, send them to a target node.
void
HierarchicalHybridJDApplication::send_sasa_summaries_upward(
	utility::vector1< HierarchicalHybridJD_SASASummaryOP > const & sasa_summaries,
	int const receiving_node
) const {
	send_request( receiving_node, OFFER_NEW_SASA_SUMMARY_BATCH_UPWARD );
	unsigned long sizebuf( sasa_summaries.size() );
	MPI_Send( &sizebuf, 1, MPI_UNSIGNED_LONG, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //First, send the size of the list.

	if(sizebuf == 0) return;

	//Buffers for data that we'll send:
	utility::vector0< int > procbuf( sizebuf );
	utility::vector0< unsigned long > jobindexbuf( sizebuf );
	utility::vector0< double > sasabuf( sizebuf );
	utility::vector0< double > polarsasabuf( sizebuf );
	utility::vector0< double > hydrophobicsasabuf( sizebuf );

	for(core::Size i=1; i<=static_cast<core::Size>(sizebuf); ++i) {
		procbuf[i-1] = sasa_summaries[i]->originating_node_MPI_rank();
		jobindexbuf[i-1] = static_cast<unsigned long>( sasa_summaries[i]->jobindex_on_originating_node() );
		sasabuf[i-1] = static_cast<double>( sasa_summaries[i]->sasa() );
		polarsasabuf[i-1] = static_cast<double>( sasa_summaries[i]->polar_sasa() );
		hydrophobicsasabuf[i-1] = static_cast<double>( sasa_summaries[i]->hydrophobic_sasa() );
	}

	MPI_Send( procbuf.data(), sizebuf, MPI_INT, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process indices.
	MPI_Send( jobindexbuf.data(), sizebuf, MPI_UNSIGNED_LONG, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process job indices.
	MPI_Send( sasabuf.data(), sizebuf, MPI_DOUBLE, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the SASA values.
	MPI_Send( polarsasabuf.data(), sizebuf, MPI_DOUBLE, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the polar SASA values.
	MPI_Send( hydrophobicsasabuf.data(), sizebuf, MPI_DOUBLE, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the hydrophobic SASA values.

	//Send the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1; i<=static_cast<core::Size>(sizebuf); ++i) {
		int nodelist_size(sasa_summaries[i]->MPI_ranks_handling_message().size());
		MPI_Send( &nodelist_size, 1, MPI_INT, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
		if(nodelist_size == 0) continue;

		utility::vector0< int > nodelist( sasa_summaries[i]->MPI_ranks_handling_message().size() );
		for( core::Size j=1, jmax=sasa_summaries[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j ) {
			nodelist[j-1] = sasa_summaries[i]->MPI_ranks_handling_message()[j];
		}
		MPI_Send( nodelist.data(), sasa_summaries[i]->MPI_ranks_handling_message().size(), MPI_INT, receiving_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
	}
}

/// @brief A child node has sent a batch of SASA summaries.  Receive them.
/// @details This function complements send_sasa_summaries_upward().
void
HierarchicalHybridJDApplication::receive_sasa_summaries_from_below(
	utility::vector1< HierarchicalHybridJD_SASASummaryOP > & sasa_summaries,
	int const requesting_node,
	bool const append_to_handler_list
) const {
	runtime_assert( sasa_summaries.empty() );

	MPI_Status status;
	unsigned long sizebuf(0);
	MPI_Recv( &sizebuf, 1, MPI_UNSIGNED_LONG, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //First, receive the size of the list.
	sasa_summaries.reserve( sizebuf );

	if(sizebuf == 0) return;

	//Buffers for data that we'll send:
	utility::vector0< int > procbuf( sizebuf );
	utility::vector0< unsigned long > jobindexbuf( sizebuf );
	utility::vector0< double > sasabuf( sizebuf );
	utility::vector0< double > polarsasabuf( sizebuf );
	utility::vector0< double > hydrophobicsasabuf( sizebuf );

	MPI_Recv( procbuf.data(), sizebuf, MPI_INT, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process indices.
	MPI_Recv( jobindexbuf.data(), sizebuf, MPI_UNSIGNED_LONG, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process job indices.
	MPI_Recv( sasabuf.data(), sizebuf, MPI_DOUBLE, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the SASA values.
	MPI_Recv( polarsasabuf.data(), sizebuf, MPI_DOUBLE, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the polar SASA values.
	MPI_Recv( hydrophobicsasabuf.data(), sizebuf, MPI_DOUBLE, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the hydrophobic SASA values.

	for(core::Size i=1, imax=static_cast<core::Size>(sizebuf); i<=imax; ++i) {
		sasa_summaries.push_back(
			utility::pointer::make_shared< HierarchicalHybridJD_SASASummary >( procbuf[i-1], jobindexbuf[i-1], sasabuf[i-1], polarsasabuf[i-1], hydrophobicsasabuf[i-1] )
		);
	}

	runtime_assert(sasa_summaries.size() == static_cast<core::Size>(sizebuf));

	//Receive the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1, imax=sasa_summaries.size(); i<=imax; ++i) {
		int nodelist_size(0);
		MPI_Recv( &nodelist_size, 1, MPI_INT, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );

		if(nodelist_size != 0) {
			utility::vector0< int > nodelist( nodelist_size );
			MPI_Recv( nodelist.data(), nodelist_size, MPI_INT, requesting_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );
			for( core::Size j=1; j<=static_cast<core::Size>(nodelist_size); ++j ) {
				sasa_summaries[i]->add_MPI_rank_handling_message( nodelist[j-1] );
			}
		}
		if( append_to_handler_list ) {
			sasa_summaries[i]->add_MPI_rank_handling_message( requesting_node ); //Add the originating node to the list of MPI ranks that handled this message.
		}

	}
}

/// @brief Recieve all SASA summaries from all children, and sort them into the order that jobsummaries is in.
/// @details This is intended for use with send_sasa_summaries_upward().
void
HierarchicalHybridJDApplication::receive_and_sort_all_sasa_summaries(
	utility::vector1< HierarchicalHybridJD_SASASummaryOP > & sasa_summaries,
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const & job_summary_list,
	bool const append_to_handler_list
) const {
	runtime_assert( sasa_summaries.empty() );
	if( job_summary_list.size() == 0 ) return;
	sasa_summaries.resize( job_summary_list.size(), nullptr );

	core::Size n_summaries_received(0);
	core::Size const n_children( my_children_.size() );

	do {
		int requesting_node(0);
		HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_below(NULL_MESSAGE);

		wait_for_request( requesting_node, message_from_below );
		runtime_assert_string_msg( message_from_below == OFFER_NEW_SASA_SUMMARY_BATCH_UPWARD, "Error in HierarchicalHybridJDApplication::receive_and_sort_all_sasa_summaries(): Children nodes of node " + std::to_string( MPI_rank_ ) + " should be transmitting SASA summaries at this point, but I received something else from node " + std::to_string( requesting_node ) + "!" );

		utility::vector1< HierarchicalHybridJD_SASASummaryOP > new_sasa_summaries;
		receive_sasa_summaries_from_below( new_sasa_summaries, requesting_node, append_to_handler_list );
		++n_summaries_received;

		// Sort the new list into the old:
		core::Size counter(1);
		for( core::Size i(1), imax( job_summary_list.size() ); i<=imax; ++i ) { //Loop through the job summary list, placing each summary.
			//Note that the received list is already sorted in the same order as the job_summary_list, except that the job_summary_list is longer (has more elements).
			if( counter > new_sasa_summaries.size() ) break; //Stop when all have been placed.

			if( job_summary_list[i]->originating_node_MPI_rank() == new_sasa_summaries[counter]->originating_node_MPI_rank() ) {
				runtime_assert( job_summary_list[i]->jobindex_on_originating_node() == new_sasa_summaries[counter]->jobindex_on_originating_node() );
				runtime_assert( sasa_summaries[i] == nullptr );
				sasa_summaries[i] = new_sasa_summaries[counter];
				++counter;
			}
		}

		//Sanity checks
		runtime_assert( counter == new_sasa_summaries.size() + 1 );

	} while( n_summaries_received < n_children );  //Keep looping until all children have sent summaries.

	//Another sanity check:
	for( core::Size i(1), imax(sasa_summaries.size()); i<=imax; ++i ) {
		runtime_assert( sasa_summaries[i] != nullptr );
	}
}

/// ------------- Director Methods --------------------

/// @brief Is this an director (root) node?
/// @details The director is responsible for sending out and retrieving all jobs, and for all file output.
bool
HierarchicalHybridJDApplication::i_am_director() const {
	return MPI_rank_ == 0;
}

/// @brief The jobs done by the director during parallel execution.
/// @details The director is responsible for sending out and retrieving all jobs, and for all file output.
void
HierarchicalHybridJDApplication::run_director() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size const nstruct( option[out::nstruct].user() ? option[out::nstruct]() : 1 ); //Number of structures to generate

	bool const halt_after_timeout( option[basic::options::OptionKeys::cyclic_peptide::MPI_stop_after_time].user() ); //Should we stop distributing jobs after a certain amount of time has elapsed?
	core::Size const timeout( halt_after_timeout ? option[basic::options::OptionKeys::cyclic_peptide::MPI_stop_after_time]() : 0 );
	bool halt_signal_fired(false);
	bool halt_signal_sent_downward(false);

	core::Size structures_remaining_to_send( nstruct );

	core::Size n_children_done(0);
	core::Size n_summaries_received(0);
	core::Size const n_children( my_children_.size() );
	core::Size total_jobs_attempted_by_children(0);
	core::Size children_reporting_job_counts(0);

	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > results_summary_from_below;
	utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > rmsds_to_best_pose; //Only computed if compute_rmsd_to_lowest_ is true.
	utility::vector1< HierarchicalHybridJD_SASASummaryOP > sasa_summaries; //Only computed if compute_sasa_metrics_ is true.

	utility::vector1< HierarchicalHybridJD_PNearToArbitraryStateSummaryCOP > pnears_to_lowest_fract; //Only computed if compute_pnear_to_this_fract_ is true.

	clock_t start_time( clock() ); //The start of the run, for timing information.
	go_signal(); //Send signal to everyone that it's time to start.

	do {
		int requesting_node(0);
		HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_below(NULL_MESSAGE);

		wait_for_request( requesting_node, message_from_below );
		switch( message_from_below ) {
			case OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD:
				receive_jobs_attempted_count( total_jobs_attempted_by_children, requesting_node );
				++children_reporting_job_counts;
			break;
			case REQUEST_NEW_JOBS_BATCH_UPWARD:
				if( halt_after_timeout && halting_condition(start_time, timeout) ) {
					halt_signal_fired = true;
					structures_remaining_to_send = 0;
				}
				send_new_jobs_downward( requesting_node, structures_remaining_to_send, batchsize_ );
				if( halt_signal_fired && !halt_signal_sent_downward && hierarchy_level_ < total_hierarchy_levels_ - 1 ) {
					send_halt_signal( my_children_ );
					halt_signal_sent_downward = true;
				}
				break;
			case OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD:
				receive_and_sort_job_summaries( results_summary_from_below, requesting_node, true);
				++n_summaries_received;
				if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") has received " << n_summaries_received << " set(s) of job summaries, representing " << results_summary_from_below.size() << " job(s), the last from node " << requesting_node << "." << std::endl;
				break;
			case GIVE_COMPLETION_SIGNAL_UPWARD:
				++n_children_done;
				//TODO -- cover other communications here as other options are added.
				break;
			default:
				utility_exit_with_message( "Unknown signal (" + std::to_string(message_from_below) + ") received from below!" );
				break;
		}
	} while( n_summaries_received < n_children || n_children_done < n_children || children_reporting_job_counts < n_children );  //Keep looping until all children have reported that they're done and have sent summaries.

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") has received completion signals, job attempt counts, and job summaries from all children." << std::endl;

	if( compute_rmsd_to_lowest_ ) {
		if( !results_summary_from_below.empty() ) {
			top_pose_bcast( results_summary_from_below[1], nullptr, false );
			if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") is waiting for RMSD-to-best summaries from nodes below." << std::endl;
			receive_and_sort_all_rmsd_to_best_summaries( rmsds_to_best_pose, results_summary_from_below, true );
			if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") has received " << rmsds_to_best_pose.size() << " RMSD-to-best summaries from nodes below." << std::endl;
		} else {
			top_pose_bcast( nullptr, nullptr, true );
			derivedTR_.Warning << "Warning: no structures were generated, so no RMSDs to best pose are expected." << std::endl;
		}
	}

	// If we're computing the PNears to the lowest-energy fraction of structures, do so here:
	if( compute_pnear_to_this_fract_ > 0.0 ) {
		director_compute_pnear_to_lowest_fract( compute_pnear_to_this_fract_, results_summary_from_below, pnears_to_lowest_fract );
	}

	if( compute_sasa_metrics_ ) {
		if( !results_summary_from_below.empty() ) {
			if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") is requesting SASA summaries from nodes below." << std::endl;
			director_send_request_for_sasa_summaries_downward();
			receive_and_sort_all_sasa_summaries( sasa_summaries, results_summary_from_below, true );
			if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") has received " << sasa_summaries.size() << " SASA summaries from nodes below." << std::endl;
		} else {
			derivedTR_.Warning << "Warning: no structures were generated, so no SASA summaries are expected." << std::endl;
		}
	}

	director_write_summaries_to_tracer( results_summary_from_below, rmsds_to_best_pose, sasa_summaries, pnears_to_lowest_fract );

	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > summary_shortlist;
	director_select_best_summaries( summary_shortlist, results_summary_from_below, output_fraction_, select_highest_ );  //This populates summary_shortlist.

	core::Size children_receiving_nonzero_requests( 0 ); //How many of my children are receiving a request for at least one pose?
	request_poses_from_below( summary_shortlist, children_receiving_nonzero_requests ); //This makes new lists of summaries, each of which it sends to the appropriate node for transmission down the hierarchy.

	core::Size children_that_sent_poses(0);
	std::string results_collected_from_below("");
	if( children_receiving_nonzero_requests > 0) {
		do {
			int requesting_node(0);
			HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_above_or_below(NULL_MESSAGE);
			wait_for_request( requesting_node, message_from_above_or_below );

			switch( message_from_above_or_below ) {
				case OFFER_NEW_POSE_BATCH_UPWARD:
					receive_pose_batch_as_string( requesting_node, results_collected_from_below );
					++children_that_sent_poses;
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Director (proc " << MPI_rank_ << ") received structures from node " << requesting_node << "." << std::endl;
					break;
				default:
					break;
			}
		} while( children_that_sent_poses < children_receiving_nonzero_requests );
	}

	if(derivedTR_.visible()) derivedTR_ << "Director is writing final structures to disk." << std::endl;
	director_write_to_disk( results_collected_from_below );

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Director (proc " << MPI_rank_ << ") reached end of run_director() function." << std::endl;
	stop_signal(); //Wait for signal to everyone that it's time to stop.
	clock_t end_time( clock() );

	if(derivedTR_summary_.visible()) {
		derivedTR_summary_ << "The simple_cycpep_predict application completed " << total_jobs_attempted_by_children << " of " << nstruct
			<< " jobs in " << static_cast< core::Real >( end_time - start_time ) / static_cast<core::Real>( CLOCKS_PER_SEC )
			<< " seconds.  " << results_summary_from_below.size() << " jobs returned structures, the top "
			<< summary_shortlist.size() << " of which were written to disk." << std::endl;
		derivedTR_summary_.flush();
	}

} //run_director()

/// @brief Convert a silent struct into a character string and broadcast it to all nodes.
/// @details Intended to be used with receive_broadcast_silent_struct_and_build_pose() to allow all other nodes to receive the broadcast.
void
HierarchicalHybridJDApplication::broadcast_silent_struct_from_this_node(
	core::io::silent::SilentStructOP ss
) const {
	core::io::silent::SilentFileOptions opts;
	opts.in_fullatom(true);
	core::io::silent::SilentFileData silentfiledata( opts );
	std::stringbuf sb;
	std::ostream outstream(&sb);

	silentfiledata._write_silent_struct( (*ss), outstream );
	std::string outstring( sb.str() );

	int strlen( outstring.length() + 1 ); //The plus one is for the /0 terminal character
	utility::vector0< char > outchar( strlen );
	sprintf( outchar.data(), "%s", outstring.c_str() );

	MPI_Bcast( &strlen, 1, MPI_INT, MPI_rank_, MPI_COMM_WORLD); //Broadcast the length of the string.
	MPI_Bcast( outchar.data(), strlen, MPI_CHAR, MPI_rank_, MPI_COMM_WORLD); //Broadcast the string.
}

/// @brief Write out a summary of the jobs completed (node, job index on node, total energy, rmsd, handler path) to the summary tracer.
/// @details The RMSD to best pose vector will only be populated if the -compute_rmsd_to_lowest option is used.
void
HierarchicalHybridJDApplication::director_write_summaries_to_tracer(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_list,
	utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > const &rmsds_to_best_pose,
	utility::vector1< HierarchicalHybridJD_SASASummaryOP > const &sasa_summaries,
	utility::vector1< HierarchicalHybridJD_PNearToArbitraryStateSummaryCOP > const & pnear_to_lowest_fract_summaries
) const {
	if( !derivedTR_summary_.visible() ) return; //Do nothing if the tracer is off.
	bool const write_rmsds_to_best( !rmsds_to_best_pose.empty() );
	bool const write_sasa_summaries( !sasa_summaries.empty() );

	if( write_rmsds_to_best ) {
		runtime_assert( rmsds_to_best_pose.size() == summary_list.size() );
	}
	if( write_sasa_summaries ) {
		runtime_assert( sasa_summaries.size() == summary_list.size() );
	}

	derivedTR_summary_ << "Summary for " << summary_list.size() << " job(s) returning results:\n";
	bool const original_print_channel_name_setting( derivedTR_summary_.get_local_print_channel_name() );
	derivedTR_summary_.set_local_print_channel_name( false ); //Disable channel name printing from this point forward.

	derivedTR_summary_ << "MPI_worker_node\tJobindex_on_node\tRMSD" << (write_rmsds_to_best ? "\tRMSD_to_best" : "" ) << "\tEnergy\tHbonds\tCisPepBonds" << ( write_sasa_summaries ? "\tSASA\tPolarSASA\tApolarSASA\tPolarSASAFract\tApolarSASAFract" : "") << "\tNode_path_to_director\n";

	core::Real numerator_sasa(0), partfxn_sasa(0);
	core::Real numerator_polar_sasa(0), numerator_hydrophobic_sasa(0);
	core::Real numerator_fraction_polar_sasa(0), numerator_fraction_hydrophobic_sasa(0);

	PNearCalculator pnear_native( lambda(), kbt() );
	PNearCalculator pnear_lowest( lambda(), kbt() );

	for( core::Size i=1, imax=summary_list.size(); i<=imax; ++i ) {
		derivedTR_summary_ << summary_list[i]->originating_node_MPI_rank() << "\t" << summary_list[i]->jobindex_on_originating_node() << "\t"
			<< summary_list[i]->rmsd() << "\t" << ( write_rmsds_to_best ? std::to_string( rmsds_to_best_pose[i]->rmsd_to_best() ) + "\t" : "" )
			<< summary_list[i]->pose_energy() << "\t" << summary_list[i]->hbonds() << "\t" << summary_list[i]->cis_peptide_bonds() << "\t";
		if( write_sasa_summaries ) {
			derivedTR_summary_ << sasa_summaries[i]->sasa() << "\t" << sasa_summaries[i]->polar_sasa() << "\t"
				<< sasa_summaries[i]->hydrophobic_sasa() << "\t" << sasa_summaries[i]->fraction_polar_sasa()
				<< "\t" << sasa_summaries[i]->fraction_hydrophobic_sasa() << "\t";
		}
		for(core::Size j=1, jmax=summary_list[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j) {
			derivedTR_summary_ << summary_list[i]->MPI_ranks_handling_message()[j];
			if( j<jmax ) derivedTR_summary_ << ",";
		}
		derivedTR_summary_ << "\n";
		if( i % 128 == 0 ) derivedTR_summary_.flush(); //Flush every hundred and twenty-eight lines.

		//Calculations for PNear and other ensemble metrics:
		if( native_ != nullptr || write_rmsds_to_best || write_sasa_summaries ) {
			if(native_ != nullptr) {
				pnear_native.add_data_point( summary_list[i]->pose_energy(), summary_list[i]->rmsd() );
			}
			if( write_rmsds_to_best ) {
				pnear_lowest.add_data_point( summary_list[i]->pose_energy(), rmsds_to_best_pose[i]->rmsd_to_best() );
			}
			if( write_sasa_summaries ) {
				core::Real const Pcurrent(  std::exp( -1.0*summary_list[i]->pose_energy()/kbt() ) );
				numerator_sasa += ( sasa_summaries[i]->sasa() * Pcurrent );
				numerator_polar_sasa += ( sasa_summaries[i]->polar_sasa() * Pcurrent );
				numerator_hydrophobic_sasa += ( sasa_summaries[i]->hydrophobic_sasa() * Pcurrent );
				numerator_fraction_polar_sasa += ( sasa_summaries[i]->fraction_polar_sasa() * Pcurrent );
				numerator_fraction_hydrophobic_sasa += ( sasa_summaries[i]->fraction_hydrophobic_sasa() * Pcurrent );

				partfxn_sasa += Pcurrent;
			}
		}
	}
	derivedTR_summary_ << "End summary for " << summary_list.size() << " job(s) returning results.\n";

	if( pnear_to_lowest_fract_summaries.size() > 0 ) {
		derivedTR_summary_ << "\nSummary of PNear values for lowest-energy " << pnear_to_lowest_fract_summaries.size() << " job(s) returning results:\n";
		derivedTR_summary_ << "MPI_worker_node\tJobindex_on_node\tPNear\tKeq\t-kbt*ln(Keq)\n";
		for( core::Size i(1), imax(pnear_to_lowest_fract_summaries.size()); i<=imax; ++i ) {
			HierarchicalHybridJD_PNearToArbitraryStateSummaryCOP const & current( pnear_to_lowest_fract_summaries[i] );
			derivedTR_summary_ << current->originating_mpi_node() << "\t";
			derivedTR_summary_ << current->jobindex_on_originating_node() << "\t";
			derivedTR_summary_ << current->pnear() << "\t";
			derivedTR_summary_ << current->Keq() << "\t";
			derivedTR_summary_ << current->deltaG_folding() << "\n";
			if( i % 128 == 0 ) derivedTR_summary_.flush();
		}
		derivedTR_summary_ << "End summary of PNear values for lowest-energy " << pnear_to_lowest_fract_summaries.size() << " job(s) returning results.\n";
	}

	if(native_ && summary_list.size() > 0 ) {
		core::Real PNear(0.0), Keq(0.0), minus_kbt_ln_Keq(0.0);
		pnear_native.compute_pnear_and_dgfolding(PNear,Keq,minus_kbt_ln_Keq);

		derivedTR_summary_ << "\nPNear:\t" << PNear << "\n";
		derivedTR_summary_ << "\nKeq:\t" << Keq << "\n";
		derivedTR_summary_ << "-kB*T*ln(Keq):\t" << minus_kbt_ln_Keq << "\n";
		derivedTR_summary_ << "lambda:\t" << lambda() << "\n";
		derivedTR_summary_ << "kB*T:\t" << kbt() << "\n";
	}

	if( write_rmsds_to_best ) {
		core::Real PNearBest(0.0), KeqBest(0.0), minus_kbt_ln_KeqBest(0.0);
		pnear_lowest.compute_pnear_and_dgfolding( PNearBest, KeqBest, minus_kbt_ln_KeqBest );

		derivedTR_summary_ << "\nPNearLowest:\t" << PNearBest << "\n";
		derivedTR_summary_ << "\nKeqLowest:\t" << KeqBest << "\n";
		derivedTR_summary_ << "-kB*T*ln(KeqLowest):\t" << minus_kbt_ln_KeqBest << "\n";
		derivedTR_summary_ << "lambda:\t" << lambda() << "\n";
		derivedTR_summary_ << "kB*T:\t" << kbt() << "\n";
	}

	if( write_sasa_summaries && partfxn_sasa > 1e-14 ) {
		derivedTR_summary_ << "BoltzAvg(SASA):\t" << numerator_sasa/partfxn_sasa << "\n";
		derivedTR_summary_ << "BoltzAvg(PolarSASA):\t" << numerator_polar_sasa/partfxn_sasa << "\n";
		derivedTR_summary_ << "BoltzAvg(ApolarSASA):\t" << numerator_hydrophobic_sasa/partfxn_sasa << "\n";
		derivedTR_summary_ << "BoltzAvg(PolarSASAFract):\t" << numerator_fraction_polar_sasa/partfxn_sasa << "\n";
		derivedTR_summary_ << "BoltzAvg(ApolarSASAFract):\t" << numerator_fraction_hydrophobic_sasa/partfxn_sasa << "\n";
	}

	derivedTR_summary_ << std::endl;
	derivedTR_summary_.set_local_print_channel_name( original_print_channel_name_setting );
	derivedTR_summary_.flush();
}

/// @brief Based on the sorted list of summaries, populate a short list of jobs, the results of which will be collected from below for output to disk.
/// @param[out] summary_shortlist The short list of job summaries populated by this function.
/// @param[in] summary_full_sorted_list The full list of job summaries collected from below, sorted.
/// @param[in] output_fraction The fraction of total jobs to collect.
/// @param[in] select_highest Should we select from the top of the summary list (lowest values) or from the bottom (highest)?
void
HierarchicalHybridJDApplication::director_select_best_summaries(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > &summary_shortlist,
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_full_sorted_list,
	core::Real const &output_fraction,
	bool const select_highest
) const {
	if(summary_full_sorted_list.size() == 0 ) {
		summary_shortlist.clear();
		return; //Do nothing if we have no jobs.
	}

	//If we're outputting everything or nothing, then this function has an easy job:
	if( output_fraction == 0.0) {
		summary_shortlist.clear();
		return;
	}	else if( output_fraction == 1.0 ) {
		summary_shortlist = summary_full_sorted_list;
		return;
	}
	summary_shortlist.clear();

	core::Size jobs_to_write( static_cast <core::Size> ( utility::round( output_fraction * static_cast<core::Real>( summary_full_sorted_list.size() ) ) ) );
	if( jobs_to_write == 0 && output_fraction > 0.0 && summary_full_sorted_list.size() > 0 ) jobs_to_write = 1; //Write at least one job.
	summary_shortlist.reserve(jobs_to_write);

	core::Size curjob( select_highest ? summary_full_sorted_list.size() - jobs_to_write + 1 : 1 );
	for(core::Size i=1; i<=jobs_to_write; ++i) {
		summary_shortlist.push_back( summary_full_sorted_list[curjob] );
		++curjob;
	}

	if( derivedTR_.Debug.visible() ) {
		derivedTR_.Debug << "Shortlisted jobs:\n";
		derivedTR_.Debug << "MPI_worker_node\tJobindex_on_node\tRMSD\tEnergy\tHbonds\tCisPepBonds\tNode_path_to_director\n";

		for(core::Size i=1, imax=summary_shortlist.size(); i<=imax; ++i) {
			derivedTR_.Debug << summary_shortlist[i]->originating_node_MPI_rank() << "\t" << summary_shortlist[i]->jobindex_on_originating_node() << "\t" << summary_shortlist[i]->rmsd() << "\t" << summary_shortlist[i]->pose_energy() << "\t" << summary_shortlist[i]->hbonds() << "\t" << summary_shortlist[i]->cis_peptide_bonds() << "\t";
			for(core::Size j=1, jmax=summary_shortlist[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j) {
				derivedTR_.Debug << summary_shortlist[i]->MPI_ranks_handling_message()[j];
				if( j<jmax ) derivedTR_.Debug << ",";
			}
			derivedTR_.Debug << "\n";
			derivedTR_.Debug.flush();
		}
		derivedTR_.Debug << std::endl;
		derivedTR_.Debug.flush();
	}
}

/// @brief Write all the collected results from below to disk.
/// @details Assumes silent output.
void
HierarchicalHybridJDApplication::director_write_to_disk(
	std::string const &output
) const {
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Writing to " << output_filename_ << "." << std::endl;
	utility::io::ozstream file(output_filename_);
	runtime_assert_string_msg( file,
		"Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::director_write_to_disk(): Couldn't open " + output_filename_ + " for writing!" );
	file << output << std::endl;
	file.close();
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Wrote to " << output_filename_ << "." << std::endl;
}

/// @brief The director is asking for SASA metrics to be sent up the hierarchy.
void
HierarchicalHybridJDApplication::director_send_request_for_sasa_summaries_downward() const {
	for( core::Size i(1), imax( my_children_.size() ); i<=imax; ++i ) {
		send_request( my_children_[i], REQUEST_SASA_SUMMARIES_DOWNWARD );
	}
}

/// @brief The director is asking for data with which to compute PNear values from below.
/// @details The steps are:
///     - Sends a request for each state in turn from the director node to all worker nodes.
///     - Each worker node broadcasts that state to all other worker nodes.
///     - All worker nodes compute RMSD to that state.
///     - Director collects RMSDs up the hierarchy and carries out PNear calculation.
///     - Repeat for each relevant state.
/// @note The pnears_to_lowest_fract vector is cleared and populated by this operation.
void
HierarchicalHybridJDApplication::director_compute_pnear_to_lowest_fract(
	core::Real const & fraction,
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const & results_summary_sorted,
	utility::vector1< HierarchicalHybridJD_PNearToArbitraryStateSummaryCOP > & pnears_to_lowest_fract
) const {
	std::chrono::time_point<HHJDA_CLOCK> const starttime( HHJDA_CLOCK::now() );

	runtime_assert( fraction > 0.0 && fraction <= 1.0 ); //Should be guaranteed true.
	pnears_to_lowest_fract.clear();

	//If there were no samples, we're done:
	if( results_summary_sorted.empty() ) {
		if( derivedTR_.Debug.visible() ) derivedTR_.Debug << "Director node (" << MPI_rank_ << ") is signalling that we have no samples, and can skip the PNear calculation to the lowest N% of samples." << std::endl;
		for( core::Size i : my_children_ ) {
			send_request( i, SKIP_PNEAR_TO_LOWEST_FRACT_DOWNWARD );
		}
		//Block until everything reaches here.
		MPI_Barrier( MPI_COMM_WORLD );
		return;
	}

	// If we reach here, there were samples.
	for( int i : my_children_ ) {
		if( derivedTR_.Debug.visible() ) derivedTR_.Debug << "Director node (" << MPI_rank_ << ") is signalling that we are starting the PNear calculation to the lowest N% of samples." << std::endl;
		send_request( i, BEGIN_PNEAR_TO_LOWEST_FRACT_DOWNWARD );
	}

	// Number of PNears to compute:
	core::Size const num_to_compute(
		std::min(
			std::max(
				core::Size(1),
				static_cast< core::Size >( std::round( fraction * static_cast< core::Real >( results_summary_sorted.size() ) ) )
			),
			results_summary_sorted.size()
		)
	);
	pnears_to_lowest_fract.reserve( num_to_compute );
	for( core::Size j(1); j<= num_to_compute; ++j ) {
		if( derivedTR_.Debug.visible() ) derivedTR_.Debug << "Director node (" << MPI_rank_ << ") is requesting PNear data for sample " << j << "." << std::endl;

		for( int i : my_children_ ) {
			send_request( i, REQUEST_PNEAR_TO_PARTICULAR_SAMPLE_DOWNWARD );
			//Send node index and job index downward:
			unsigned long sizebuf[2] = { static_cast< core::Size>(results_summary_sorted[j]->originating_node_MPI_rank()), results_summary_sorted[j]->jobindex_on_originating_node() };
			MPI_Send( &sizebuf, 2, MPI_UNSIGNED_LONG, i, static_cast<int>(JOB_NODE_AND_INDEX_DOWNWARD), MPI_COMM_WORLD );
		}

		//Participate of broadcast of relevant structure:
		receive_broadcast_silent_struct_and_build_pose( static_cast<int>(results_summary_sorted[j]->originating_node_MPI_rank()), true );

		//Receive results summary from below:
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > rmsds_to_best_summaries;
		receive_and_sort_all_rmsd_to_best_summaries( rmsds_to_best_summaries, results_summary_sorted, true );

		//Compute PNear and store:
		PNearCalculator pnear_calc( lambda(), kbt() );
		runtime_assert( rmsds_to_best_summaries.size() == results_summary_sorted.size() );
		for( core::Size i(1), imax(rmsds_to_best_summaries.size()); i<=imax; ++i ) {
			pnear_calc.add_data_point( results_summary_sorted[i]->pose_energy(), rmsds_to_best_summaries[i]->rmsd_to_best() );
		}
		core::Real pnear, kbt, dgfolding;
		pnear_calc.compute_pnear_and_dgfolding( pnear, kbt, dgfolding );
		pnears_to_lowest_fract.push_back(
			utility::pointer::make_shared< HierarchicalHybridJD_PNearToArbitraryStateSummary >(
				pnear, kbt, dgfolding,
				static_cast< core::Size>(results_summary_sorted[j]->originating_node_MPI_rank()),
				static_cast<core::Size>(results_summary_sorted[j]->jobindex_on_originating_node()
		) ) );
	}

	// And we're done.
	for( int i : my_children_ ) {
		if( derivedTR_.Debug.visible() ) derivedTR_.Debug << "Director node (" << MPI_rank_ << ") is signalling that we have completed the PNear calculation to the lowest N% of samples." << std::endl;
		send_request( i, END_PNEAR_TO_LOWEST_FRACT_DOWNWARD );
	}

	//Block until everything reaches here.
	MPI_Barrier( MPI_COMM_WORLD );

	//Report end time.
	if( derivedTR_.visible() ) {
		std::chrono::time_point<HHJDA_CLOCK> const endtime( HHJDA_CLOCK::now() );
		derivedTR_ << "Completed PNear computation to lowest " << fraction * 100 << "% of samples in " << std::chrono::duration_cast< std::chrono::microseconds >( endtime-starttime ).count() << " microseconds." << std::endl;
	}
}

/// ------------- Intermediate Manager Methods --------

/// @brief Is this an intermediate manager node?
/// @details The managers are responsible for distributing jobs to other managers lower in the hierarchy, and/or to workers, and for
/// collecting results from managers/workers lower in the hierarchy and sending them up to the director.
bool
HierarchicalHybridJDApplication::i_am_intermediate_manager() const {
	return !( i_am_director() || i_am_worker() );
}

/// @brief The jobs done by the intermediate managers during parallel execution.
/// @details The managers are responsible for distributing jobs to other managers lower in the hierarchy, and/or to workers, and for
/// collecting results from managers/workers lower in the hierarchy and sending them up to the director.
void
HierarchicalHybridJDApplication::run_intermediate_manager() const {
	go_signal(); //Wait for signal to everyone that it's time to start.

	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > jobsummaries; //Job summaries received from children.  Sorted list.
	utility::vector1 < HierarchicalHybridJD_RMSDToBestSummaryOP > rmsd_to_best_summaries; //RMSD-to-best summaries received from children.  Sorted list.  Only used if we're calculating RMSD to best; remains an empty vector otherwise.
	utility::vector1 < HierarchicalHybridJD_SASASummaryOP > sasa_summaries; //SASA summaries received from children.  Sorted list.  Only used if we're calculating SASA metrics; otherwise this remains an empty vector.

	bool halt_signal_received(false);
	core::Size const n_children( my_children_.size() ); //How many children do I have?
	core::Size n_children_done(0), n_summaries_received(0) /*, n_posesets_received(0)*/; //How many of my children have reported that they've finished all of their jobs?  How many have sent me summaries of what they've done?  How many have sent me their poses?

	core::Size n_jobs(0); //Jobs that I'm holding, waiting to assign to lower-level nodes.
	core::Size total_jobs_attempted_by_children(0); //How many jobs have my children attempted?
	core::Size children_reporting_job_counts(0); //How many children have reported job counts

	do {
		int requesting_node(0);
		HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_above_or_below(NULL_MESSAGE);
		wait_for_request( requesting_node, message_from_above_or_below );

		switch( message_from_above_or_below ) {
			case OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD:
				receive_jobs_attempted_count( total_jobs_attempted_by_children, requesting_node );
				++children_reporting_job_counts;
				if( children_reporting_job_counts == n_children ) {
					//Offer the total number of jobs attempted upward:
					send_request( my_parent_, OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD);
					send_jobs_attempted_count( total_jobs_attempted_by_children,	my_parent_);
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager process " << MPI_rank_ << " reported to its parent that its children attempted " << total_jobs_attempted_by_children << " jobs." << std::endl;
				}
			break;
			case HALT_SIGNAL:
				if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Proc " << MPI_rank_ << " received halt signal from proc " << requesting_node << "." << std::endl;
				halt_signal_received = true;
				n_jobs = 0; //Clear all my jobs; don't send any more downwards.
				if( hierarchy_level_ < total_hierarchy_levels_ - 1 ) {
					send_halt_signal( my_children_ );
				}
			break;
			case REQUEST_NEW_JOBS_BATCH_UPWARD:
				if( n_jobs < batchsize_ && !halt_signal_received ) {
					send_request( my_parent_, REQUEST_NEW_JOBS_BATCH_UPWARD );
					receive_njobs_from_above( n_jobs );
				}
				send_new_jobs_downward( requesting_node, n_jobs, batchsize_ );
				break;
			case GIVE_COMPLETION_SIGNAL_UPWARD:
				//Afll jobs in one worker node are done.  Check off this child as done.  When all children are done, send the same signal upward.
				++n_children_done;
				if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager process " << MPI_rank_ << " reports that child node " << requesting_node <<  " has completed." << std::endl;
				if( n_children_done == n_children ) {
					send_request( my_parent_, GIVE_COMPLETION_SIGNAL_UPWARD);
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager process " << MPI_rank_ << " reports that all child nodes have completed." << std::endl;
				}
				break;
			case OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD:
				//A worker is offering a summary of the jobs that it completed.  Signal for it to send the summary upward, and sort the
				//list into the list that we're holding.
				if(derivedTR_.Debug.visible()) {
					derivedTR_.Debug << "Intermediate manager process " << MPI_rank_ << " has been offered new job summaries from node " << requesting_node << ".  Receiving..." << std::endl;
					derivedTR_.Debug.flush();
				}
				receive_and_sort_job_summaries( jobsummaries, requesting_node, true );
				++n_summaries_received;
				if(derivedTR_.Debug.visible()) {
					derivedTR_.Debug << "Intermediate manager process " << MPI_rank_ << " has received a total of " << n_summaries_received << " summary set(s) for a total of " << jobsummaries.size() << " job(s), the last from node " << requesting_node << "." << std::endl;
					derivedTR_.Debug.flush();
				}
				//Once the summary of all worker jobs is in hand, offer and then send the full list upward.
				if(n_summaries_received == n_children) {
					send_request( my_parent_, OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD);
					send_job_summaries( jobsummaries, my_parent_);
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager process " << MPI_rank_ << " sent job summaries to node " << my_parent_ << "." << std::endl;
				}
				break;
			case OFFER_NEW_RMSD_TO_BEST_SUMMARY_BATCH_UPWARD:
				utility_exit_with_message( "Error!  Intermediate manager " + std::to_string(MPI_rank_) + " received an offer of RMSD-to-best summaries from node " + std::to_string( requesting_node ) + ", but we're still collecting job summaries!  This should not be possible.");
			case BEGIN_PNEAR_TO_LOWEST_FRACT_DOWNWARD:
			case SKIP_PNEAR_TO_LOWEST_FRACT_DOWNWARD:
			case REQUEST_PNEAR_TO_PARTICULAR_SAMPLE_DOWNWARD:
			case END_PNEAR_TO_LOWEST_FRACT_DOWNWARD:
				utility_exit_with_message( "Error!  Intermediate manager " + std::to_string(MPI_rank_) + " received a message related to PNear computation to the lowest energy N% of structures, but we're not yet at that stage of the protocol.  We're still collecting job summaries!  This should not be possible." );
				//TODO Handle other cases here as they are added.
			default:
				utility_exit_with_message( "Unknown signal (" + std::to_string( message_from_above_or_below ) + ") received!" );
				break;
		}
	} while( n_children_done < n_children || n_summaries_received < n_children || children_reporting_job_counts < n_children );

	// Relay information about RMSDs to best pose, if we're computing RMSD to best:
	if( compute_rmsd_to_lowest_ ) {
		if( !top_pose_bcast() ) {
			receive_and_sort_all_rmsd_to_best_summaries( rmsd_to_best_summaries, jobsummaries, true );
			send_rmsds_to_best_summaries_upward( rmsd_to_best_summaries, my_parent_ );
		}
	}

	// If we're computing the PNears to the lowest-energy fraction of structures, do so here:
	if( compute_pnear_to_this_fract_ ) {
		intermediate_manager_compute_pnear_to_lowest_fract(jobsummaries);
	}

	// Relay information about SASA, if we're computing SASA metrics:
	if( compute_sasa_metrics_ ) {
		intermediate_manager_relay_request_for_sasa_summaries_downward();
		receive_and_sort_all_sasa_summaries( sasa_summaries, jobsummaries, true );
		send_sasa_summaries_upward( sasa_summaries, my_parent_ );
	}

	//Transmit requests for full poses down the hierarchy:
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > requested_jobs;
	receive_pose_requests_from_above( requested_jobs );
	core::Size children_receiving_nonzero_requests( 0 ); //How many of my children are receiving a request for at least one pose?
	request_poses_from_below( requested_jobs, children_receiving_nonzero_requests );

	core::Size children_that_sent_poses(0);
	std::string results_collected_from_below("");
	if( children_receiving_nonzero_requests > 0) {
		do {
			int requesting_node(0);
			HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_above_or_below(NULL_MESSAGE);
			wait_for_request( requesting_node, message_from_above_or_below );

			switch( message_from_above_or_below ) {
				case OFFER_NEW_POSE_BATCH_UPWARD:
					receive_pose_batch_as_string( requesting_node, results_collected_from_below );
					++children_that_sent_poses;
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager proc " << MPI_rank_ << " received structures from node " << requesting_node << ".  " << children_that_sent_poses << " of " << children_receiving_nonzero_requests << " have sent structures." << std::endl;
					break;
				default:
					break;
			}
		} while( children_that_sent_poses < children_receiving_nonzero_requests );
	}

	if(children_receiving_nonzero_requests > 0 ) {
		send_request( my_parent_, OFFER_NEW_POSE_BATCH_UPWARD);
		intermediate_manager_send_poses_as_string_upward( results_collected_from_below, my_parent_ );
	}
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager proc " << MPI_rank_ << " sent structures to node " << my_parent_ << "." << std::endl;

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate manager proc " << MPI_rank_ << " reached end of run_intermediate_manager() function." << std::endl;
	stop_signal(); //Wait for signal to everyone that it's time to stop.
} //run_intermediate_manager()

/// @brief Relay the jobs received from below, held as a concatenated string, up the hierarchy.
/// @details Transmission to be received with receive_pose_batch_as_string().
void
HierarchicalHybridJDApplication::intermediate_manager_send_poses_as_string_upward(
	std::string const &results,
	int const target_node
) const {
	int strlen( results.length() + 1 ); //The +1 is for the null terminator ('/0').
	MPI_Send( &strlen, 1, MPI_INT, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the length of the string.

	//Irritating.  MPI_Send wants a non-const pointer, so I have to copy my const string.
	utility::vector0< char > str_to_send( strlen );
	for(core::Size i=0; i<static_cast<core::Size>(strlen); ++i) { str_to_send[i] = results.c_str()[i]; }
	MPI_Send( str_to_send.data(), strlen, MPI_CHAR, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the string.
}

/// @brief Receive a request for SASA summaries from above, and send it to all children.
/// @details This function expects that the only possible message that can be received
/// at this point is the request for SASA summaries!
void
HierarchicalHybridJDApplication::intermediate_manager_relay_request_for_sasa_summaries_downward() const {
	HIERARCHICAL_MPI_COMMUNICATION_TYPE request(NULL_MESSAGE);
	int requesting_node(0);
	wait_for_request( requesting_node, request );
	runtime_assert( request == REQUEST_SASA_SUMMARIES_DOWNWARD && requesting_node == static_cast<int>( my_parent_ ) );
	for( core::Size i(1), imax(my_children_.size()); i<=imax; ++i ) {
		send_request( my_children_[i], REQUEST_SASA_SUMMARIES_DOWNWARD );
	}
}

/// @brief Send requests for PNear data for the lowest-energy N% of structures down
/// the hierarchy, and facilitate results going up the hierarchy.
void
HierarchicalHybridJDApplication::intermediate_manager_compute_pnear_to_lowest_fract(
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > const & jobsummaries
) const {
	HIERARCHICAL_MPI_COMMUNICATION_TYPE request( NULL_MESSAGE );
	int requesting_node(0);

	//First, we get a request either to skip this phase or to perform this phase.
	wait_for_request( requesting_node, request );
	if( request == SKIP_PNEAR_TO_LOWEST_FRACT_DOWNWARD ) {
		for( int i : my_children_ ) {
			send_request( i, SKIP_PNEAR_TO_LOWEST_FRACT_DOWNWARD );
		}
		//Block until everything reaches here.
		MPI_Barrier( MPI_COMM_WORLD );
		return;
	} else if ( request == BEGIN_PNEAR_TO_LOWEST_FRACT_DOWNWARD ) {
		for( int i : my_children_ ) {
			send_request( i, BEGIN_PNEAR_TO_LOWEST_FRACT_DOWNWARD );
		}
	} else {
		utility_exit_with_message( "Intermediate manager node " + std::to_string( MPI_rank_ ) + " recieved signal " + std::to_string( static_cast<int>(request) ) + ".  This was unexpected." );
	}

	//If we reach here, there were samples.  Loop until we get a message indicating that we've done them all.
	do {
		// Receive and send a request for a sample OR a requeset to end.
		wait_for_request( requesting_node, request );
		if( request == END_PNEAR_TO_LOWEST_FRACT_DOWNWARD ) {
			for( int i: my_children_ ) {
				send_request( i, END_PNEAR_TO_LOWEST_FRACT_DOWNWARD );
			}
			break;
		}
		runtime_assert_string_msg( request == REQUEST_PNEAR_TO_PARTICULAR_SAMPLE_DOWNWARD, "Intermediate manager node " + std::to_string( MPI_rank_ ) + " recieved signal " + std::to_string( static_cast<int>(request) ) + ".  This was unexpected." );
		for( int i: my_children_ ) {
			send_request( i, REQUEST_PNEAR_TO_PARTICULAR_SAMPLE_DOWNWARD );
		}

		//Receive and retransmit node index and job index downard.
		unsigned long sizebuf[2] = { 0, 0 };
		MPI_Status mpi_status;
		MPI_Recv( &sizebuf, 2, MPI_UNSIGNED_LONG, requesting_node, static_cast<int>(JOB_NODE_AND_INDEX_DOWNWARD), MPI_COMM_WORLD, &mpi_status );
		for( int i : my_children_ ) {
			MPI_Send( &sizebuf, 2, MPI_UNSIGNED_LONG, i, static_cast<int>(JOB_NODE_AND_INDEX_DOWNWARD), MPI_COMM_WORLD );
		}

		//Participate in broadcast of relevant structure:
		receive_broadcast_silent_struct_and_build_pose( static_cast<int>(sizebuf[0]), true );

		//Receive results summary from below:
		utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > rmsds_to_best_summaries;
		receive_and_sort_all_rmsd_to_best_summaries( rmsds_to_best_summaries, jobsummaries, true );
		send_rmsds_to_best_summaries_upward( rmsds_to_best_summaries, my_parent_ );

	} while(true);

	//Block until everything reaches here.
	MPI_Barrier( MPI_COMM_WORLD );
}


/// ------------- Worker Methods ----------------------

/// @brief Is this a worker (or worker) node?
/// @details The workers receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
bool
HierarchicalHybridJDApplication::i_am_worker() const {
	return hierarchy_level_ == total_hierarchy_levels_;
}

/// @brief The jobs done by the workers during parallel execution.
/// @details The workers receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
/// @note In multi-threaded mode, if threads_per_worker_process_ is greater than 1, then the workers launch threads
/// to do the work.  Only the manager thread does MPI calls.
void
HierarchicalHybridJDApplication::run_worker() const {
	go_signal(); //Wait for signal to everyone that it's time to start.

	//Temporary lines for GDB debugging with MPI.
	//core::Size i(0);
	//while(i==0) {
	//}

#ifdef MULTI_THREADED
	core::Size batchcount(0); //Counter for total number of batches sent out to worker threads
#endif

	core::Size total_jobs_attempted(0);
	core::Size njobs_from_above(0);
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > jobsummaries;
	utility::vector1 < core::io::silent::SilentStructOP > all_output; //All of the output poses, in silent structure format..

	do {
		send_request( my_parent_, REQUEST_NEW_JOBS_BATCH_UPWARD );
		receive_njobs_from_above( njobs_from_above );
		if (derivedTR_.Debug.visible()) derivedTR_.Debug << "Worker process " << MPI_rank_ << " was assigned " << njobs_from_above << " job(s)." << std::endl;
		if( njobs_from_above != 0 ) {
			total_jobs_attempted += njobs_from_above;
#ifdef MULTI_THREADED
			core::Size const jobs_in_this_batch( njobs_from_above );
			++batchcount; //We're sending out another batch.
			// In multi-threaded mode, if we're launching more than one thread per worker process, do that here.
			if( threads_per_worker_proc_ > 1 ) {
				utility::vector1< std::thread > worker_threads;
				for(core::Size i(1); i<threads_per_worker_proc_; ++i) { //Note that one thread (this one) is the MPI communication thread.  There are N-1 worker threads.
					core::pose::PoseOP native_copy;
					if( native_ != nullptr ) {
						native_copy = utility::pointer::make_shared< core::pose::Pose >();
						native_copy->detached_copy(*native_);
					}

					worker_threads.push_back( std::thread( std::bind(
					   &protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::worker_carry_out_njobs_in_thread,
					   this,
					   i,
					   &njobs_from_above,
					   jobs_in_this_batch,
					   &jobsummaries,
					   &all_output,
					   scorefxn_->clone(),
					   native_copy,
					   sequence_,
					   batchcount
					) ) ); //Note that, in the C++11 standard, objects created with push_back() are moved rather than copied into the container, so this is correct.
				}
				for(core::Size i(1); i<threads_per_worker_proc_; ++i) {
					worker_threads[i].join();
				}

				{ //Scope for lock.
					//Update the indices of jobs in the jobsummaries list.  Note that this is done between threads being open, and requires no mutex locking.
					//I'm going to lock anyways to be extra safe.
					std::lock_guard < std::mutex > lock( results_mutex_ );
					for(core::Size i(1), imax(jobsummaries.size()); i<=imax; ++i) {
						jobsummaries[i]->set_jobindex_on_originating_node(i);
					}
				}

			} else {
				worker_carry_out_njobs( njobs_from_above, jobsummaries, all_output );
			}
#else
			worker_carry_out_njobs( njobs_from_above, jobsummaries, all_output );
#endif
		} else {
			break;
		}
	} while( true );

	send_request( my_parent_, GIVE_COMPLETION_SIGNAL_UPWARD);
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Worker process " << MPI_rank_ << " reporting completion." << std::endl;
	send_request( my_parent_, OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD);
	sort_jobsummaries_list( jobsummaries, sort_type_ ); //Workers sort the job summaries lists, so that merge sorts can be used at all steps by above layers.
	send_job_summaries( jobsummaries, my_parent_ );
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Worker process " << MPI_rank_ << " sent job summaries for " << jobsummaries.size() << " job(s) to node " << my_parent_ << "." << std::endl;

	//Offer the total number of jobs attempted upward:
	send_request( my_parent_, OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD);
	send_jobs_attempted_count( total_jobs_attempted, my_parent_);
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Worker process " << MPI_rank_ << " reported to its parent that it attempted " << total_jobs_attempted << " jobs." << std::endl;

	//If we're computing RMSD to top pose, receive the top pose and compute:
	if( compute_rmsd_to_lowest_ ) {
		core::io::silent::SilentFileOptions opts;
		opts.in_fullatom(true);
		core::io::silent::SilentStructOP top_pose_silentstruct( all_output.empty() ? utility::pointer::make_shared< core::io::silent::BinarySilentStruct >(opts) : all_output[ jobsummaries[1]->jobindex_on_originating_node() ]->clone() );
		if( !top_pose_bcast( nullptr, top_pose_silentstruct ) ) {
			utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > rmsds_to_best_summaries;
			worker_compute_sorted_rmsds_to_best( *top_pose_silentstruct, rmsds_to_best_summaries, jobsummaries, all_output );
			send_rmsds_to_best_summaries_upward( rmsds_to_best_summaries, my_parent_ );
		}
	}

	// If we're computing the PNears to the lowest-energy fraction of structures, do so here:
	if( compute_pnear_to_this_fract_ ) {
		worker_compute_pnear_to_lowest_fract( jobsummaries, all_output );
	}

	//If we're computing SASA metrics, receive the go signal and proceed:
	if( compute_sasa_metrics_ ) {
		utility::vector1< HierarchicalHybridJD_SASASummaryOP > sasa_summaries;
		worker_receive_request_for_sasa_summaries();
		worker_generate_and_sort_sasa_summaries( sasa_summaries, jobsummaries, all_output );
		send_sasa_summaries_upward( sasa_summaries, my_parent_ );
	}

	//Receive requests for full poses from up the hierarchy:
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > requested_jobs;
	receive_pose_requests_from_above( requested_jobs );
	if(derivedTR_.Debug.visible() && requested_jobs.size() > 0) { //Debug output only: list the jobs that this worker node was instructed to send up the hierarchy.
		derivedTR_.Debug << "Worker process " << MPI_rank_ << " received requests for the following job(s): ";
		for(core::Size i=1, imax=requested_jobs.size(); i<=imax; ++i) {
			derivedTR_.Debug << requested_jobs[i]->originating_node_MPI_rank() << "-" << requested_jobs[i]->jobindex_on_originating_node();
			if(i<imax) derivedTR_.Debug << ", ";
		}
		derivedTR_.Debug << "." << std::endl;
	}

	if(requested_jobs.size() > 0) {	//Only those workers that received nonzero requests bother to send anything.
		send_request( my_parent_, OFFER_NEW_POSE_BATCH_UPWARD);
		worker_send_poses_upward( requested_jobs, all_output );
	}

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Worker proc " << MPI_rank_ << " reached end of run_worker() function." << std::endl;
	derivedTR_.flush();

	stop_signal(); //Wait for signal to everyone that it's time to stop.
} //run_worker()

/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
///
void
HierarchicalHybridJDApplication::worker_carry_out_njobs(
	core::Size &njobs_from_above,
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > &all_output
) const {
	if (derivedTR_.Debug.visible()) {
		derivedTR_.Debug << "Starting " << njobs_from_above << " job(s) on worker node " << MPI_rank_ << "." << std::endl;
	}

	try {
		derived_worker_carry_out_n_jobs( njobs_from_above, jobsummaries, all_output, scorefxn_, native_, sequence_ );
		worker_job_count_ += njobs_from_above;
	} catch ( utility::excn::Exception &excn ) {
		derivedTR_.Error << "Exception in SimpleCycpepPredictApplication caught:" << std::endl;
		excn.show( derivedTR_.Error );
		derivedTR_.Error << "\nRecovering from error and continuing to next job." << std::endl;
		derivedTR_.Error.flush();
	}

	njobs_from_above=0;
} //worker_carry_out_njobs()

#ifdef MULTI_THREADED

/// @brief Decrement the job counter.  Return true if the job counter was greater than zero.
/// @details Does this with proper locking to prevent threads from stepping on one another.
bool
HierarchicalHybridJDApplication::worker_decrement_jobcount_multithreaded(
	core::Size * available_job_count,
	core::Size &already_completed_job_count,
	core::Size const jobs_in_this_batch,
	core::Size const /*thread_index*/
) const {
	std::lock_guard< std::mutex > lock( joblist_mutex_ ); //Lock the mutex
	if( (*available_job_count) == 0 ) return false;
	already_completed_job_count = worker_job_count_ + ( jobs_in_this_batch - (*available_job_count) );
	--(*available_job_count);
	//TR << "Thread " << thread_index << " decrementing available jobs.  " << already_completed_job_count << " already completed.  " << (*available_job_count) << " jobs now remain." << std::endl;
	return true;
}

/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
/// @details This is the multi-threaded version, which locks the job count to decrement it, and locks the jobsummaries and all_output vectors to add to them.
void
HierarchicalHybridJDApplication::worker_carry_out_njobs_in_thread(
	core::Size const thread_index,
	core::Size * njobs_from_above,
	core::Size const jobs_in_this_batch,
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > * jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > * all_output,
	core::scoring::ScoreFunctionOP sfxn,
	core::pose::PoseOP native,
	std::string const sequence, //Deliberately copied
	core::Size const batch_index //Number of batches that have been sent out on this proc.
) const {
	derivedTR_.Debug << "Launching worker thread " << thread_index << " from worker process " << MPI_rank_ << std::endl;

	// Initialize the random number generators for this thread:
	int seed(1111111);
	if( !use_const_random_seed_ ) {
		seed = time(nullptr);
	}
	seed += (thread_index - 1)*MPI_n_procs_ + MPI_rank_ + (threads_per_worker_proc_ * MPI_n_procs_ * (batch_index - 1) ); //Unique seed for each proc, thread, and job batch.
	basic::random::init_random_generators( seed, rgtype_ );
	derivedTR_.Debug << "Using random seed " << seed << " for worker thread " << thread_index << " from worker process " << MPI_rank_ << std::endl;

	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > local_jobsummaries;
	utility::vector1< core::io::silent::SilentStructOP > local_all_output;
	core::Size already_completed_job_count(0);

	while( worker_decrement_jobcount_multithreaded(njobs_from_above, already_completed_job_count, jobs_in_this_batch, thread_index) ) {
		try {
			//Create and initialize the relevant application:
			derived_worker_carry_out_n_jobs(1, local_jobsummaries, local_all_output, sfxn, native, sequence);
			{ //Scope for lock
				std::lock_guard< std::mutex > lock( joblist_mutex_ ); //Lock the mutex to access worker_job_count_.
				worker_job_count_ += 1;
			} //Unlock here
		} catch ( utility::excn::Exception &excn ) {
			derivedTR_.Error << "Exception in SimpleCycpepPredictApplication caught:" << std::endl;
			excn.show( derivedTR_.Error );
			derivedTR_.Error << "\nRecovering from error and continuing to next job." << std::endl;
			derivedTR_.Error.flush();
		}
	} //Loop over available jobs.

	{ //Lock result mutex for this scope.
		std::lock_guard< std::mutex > lock( results_mutex_ );
		debug_assert( local_jobsummaries.size() == local_all_output.size() );
		for(core::Size i(1), imax(local_jobsummaries.size()); i<=imax; ++i) {
			(*jobsummaries).push_back( local_jobsummaries[i] );
			(*all_output).push_back( local_all_output[i] );
			//TR << "Thread " << thread_index << " adding summary for " << (local_jobsummaries[i])->originating_node_MPI_rank() << "-" << (local_jobsummaries[i])->jobindex_on_originating_node() << "." << std::endl;
		}
	} //Unlock results mutex.

	derivedTR_.Debug << "Terminating worker thread " << thread_index << " from worker process " << MPI_rank_ << std::endl;

} //worker_carry_out_njobs_in_thread()
#endif

/// @brief If we're computing the RMSDs to the very best pose, do so.
/// @details This function clears and populates the rmsds_to_best_summaries vector, ensuring that RMSDs to the
/// pose represented by top_pose_silentstruct are computed in the order that matches jobsummaries.
void
HierarchicalHybridJDApplication::worker_compute_sorted_rmsds_to_best(
	core::io::silent::SilentStruct const &top_pose_silentstruct,
	utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > & rmsds_to_best_summaries,
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &jobsummaries,
	utility::vector1< core::io::silent::SilentStructOP > const &poses_from_this_worker
) const {
	runtime_assert( rmsds_to_best_summaries.empty() );
	runtime_assert( jobsummaries.size() == poses_from_this_worker.size() );
	if( jobsummaries.size() == 0 ) return;

	rmsds_to_best_summaries.reserve( jobsummaries.size() );

	core::chemical::ResidueTypeSet const & restypeset( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ) ) );

	// Rebuild the top pose from the binary silent structure:
	core::pose::PoseOP top_pose( utility::pointer::make_shared< core::pose::Pose >() );
	top_pose_silentstruct.fill_pose( *top_pose, restypeset );

	for( core::Size isummary(1), isummarymax(jobsummaries.size()); isummary<=isummarymax; ++isummary ) { //Loop through all job summaries.
		core::Size const cur_job_index( jobsummaries[isummary]->jobindex_on_originating_node() );
		core::pose::PoseOP worker_pose( utility::pointer::make_shared< core::pose::Pose >() );
		poses_from_this_worker[ cur_job_index ]->fill_pose( *worker_pose, restypeset );

		rmsds_to_best_summaries.push_back(
			utility::pointer::make_shared< HierarchicalHybridJD_RMSDToBestSummary >(
				MPI_rank_, cur_job_index,
				derived_worker_compute_rmsd( *worker_pose, *top_pose, sequence_ )
			)
		);
		if( derivedTR_.Debug.visible() ) {
			derivedTR_.Debug << "Computed RMSD of " << rmsds_to_best_summaries[rmsds_to_best_summaries.size()]->rmsd_to_best() << " to best structure for job " << rmsds_to_best_summaries[rmsds_to_best_summaries.size()]->jobindex_on_originating_node() << " on node " << rmsds_to_best_summaries[rmsds_to_best_summaries.size()]->originating_node_MPI_rank() << "." << std::endl;
		}
	}
}

/// @brief Wait for requests from the director for RMSDs to a given structure, participate in a broadcast of
/// that structure, and compute RMSDs for all output to that structure to send up the hierarchy.
void
HierarchicalHybridJDApplication::worker_compute_pnear_to_lowest_fract(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryCOP > const & jobsummaries,
	utility::vector1 < core::io::silent::SilentStructCOP > const & all_output 
) const {
	// Receive instruction from up the hierarchy to either skip this step or begin this step:
	HIERARCHICAL_MPI_COMMUNICATION_TYPE message( NULL_MESSAGE );
	int requesting_node(0);
	wait_for_request( requesting_node, message );
	runtime_assert( requesting_node >= 0 && static_cast<core::Size>(requesting_node) == my_parent_ );
	if( message == SKIP_PNEAR_TO_LOWEST_FRACT_DOWNWARD ) {
		//Block until everything reaches here.
		MPI_Barrier( MPI_COMM_WORLD );
		return;
	}
	runtime_assert( message == BEGIN_PNEAR_TO_LOWEST_FRACT_DOWNWARD ); //Should be true.

	core::chemical::ResidueTypeSet const & restypeset( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ) ) );

	//If we reach here, we are go to do the comparisons.  Loop until we're done.
	do {
		wait_for_request( requesting_node, message );
		runtime_assert( requesting_node >= 0 && static_cast<core::Size>(requesting_node) == my_parent_ );
		if( message == REQUEST_PNEAR_TO_PARTICULAR_SAMPLE_DOWNWARD ) {
			unsigned long sizebuf[2] = {0, 0};
			MPI_Status status;
			MPI_Recv( &sizebuf, 2, MPI_UNSIGNED_LONG, my_parent_, static_cast<int>(JOB_NODE_AND_INDEX_DOWNWARD), MPI_COMM_WORLD, &status ); //Recive node and job index.

			core::pose::PoseCOP comparison_pose;
			if( static_cast<core::Size>(MPI_rank_) == sizebuf[0] ) {
				// This node should send out the structure.
				runtime_assert( sizebuf[1] <= all_output.size() );
				core::io::silent::SilentStructOP ss_copy( all_output[ sizebuf[1] ]->clone() ); //Cloned because for some reason the functions that this calls have to take a nonconst silent struct, and I don't want the silent struct changed.
				broadcast_silent_struct_from_this_node( ss_copy );
				core::pose::PoseOP comparison_pose_nonconst( utility::pointer::make_shared< core::pose::Pose >() );
				ss_copy->fill_pose( *comparison_pose_nonconst, restypeset );
				comparison_pose = comparison_pose_nonconst;
			} else {
				// This node should receive the structure.
				comparison_pose = receive_broadcast_silent_struct_and_build_pose( static_cast<int>(sizebuf[0]), false );
			}

			// Compute RMSD to all other structures.
			utility::vector1< HierarchicalHybridJD_RMSDToBestSummaryOP > rmsd_summaries;
			for( core::Size isummary(1), isummarymax(jobsummaries.size()); isummary<=isummarymax; ++isummary ) { //Loop through all job summaries.
				core::Size const cur_job_index( jobsummaries[isummary]->jobindex_on_originating_node() );
				core::pose::PoseOP worker_pose( utility::pointer::make_shared< core::pose::Pose >() );
				all_output[ cur_job_index ]->fill_pose( *worker_pose, restypeset );

				rmsd_summaries.push_back(
					utility::pointer::make_shared< HierarchicalHybridJD_RMSDToBestSummary >(
						MPI_rank_, cur_job_index,
						derived_worker_compute_rmsd( *worker_pose, *comparison_pose, sequence_ )
					)
				);
				if( derivedTR_.Debug.visible() ) {
					derivedTR_.Debug << "Computed RMSD of " << rmsd_summaries[rmsd_summaries.size()]->rmsd_to_best() << " to best structure for job " << rmsd_summaries[rmsd_summaries.size()]->jobindex_on_originating_node() << " on node " << rmsd_summaries[rmsd_summaries.size()]->originating_node_MPI_rank() << "." << std::endl;
				}
			}

			// Send RMSDs up the hierarchy.
			send_rmsds_to_best_summaries_upward( rmsd_summaries, my_parent_ );

		} else if ( message == END_PNEAR_TO_LOWEST_FRACT_DOWNWARD ) {
			//We're done, so break out of this loop.
			break;
		} else {
			utility_exit_with_message( "Program error in HierarchicalHybridJDApplication::worker_compute_pnear_to_lowest_fract(): Received signal " + std::to_string( message ) + ", but this is meaningless in this context.  Please consult a developer." );
		}
	} while( true );

	//Block until everything reaches here.
	MPI_Barrier( MPI_COMM_WORLD );

}

/// @brief Given a list of jobs that have been requested from above, send the corresponding poses up the hierarchy.
/// @details Throws an error if any jbo was completed on a different node than this worker.
void
HierarchicalHybridJDApplication::worker_send_poses_upward(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &requested_jobs,
	utility::vector1 < core::io::silent::SilentStructOP > const &all_output
) const {

	core::Size const n_requests( requested_jobs.size() );

	utility::vector1 < core::io::silent::SilentStructOP > relevant_output;
	relevant_output.reserve( n_requests );

	for(core::Size i=1; i<=n_requests; ++i) {
		runtime_assert_string_msg( requested_jobs[i]->originating_node_MPI_rank() == MPI_rank_, "Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::worker_send_poses_upward(): A worker node was asked for a job that it didn't do!" );
		core::Size const requested_index( requested_jobs[i]->jobindex_on_originating_node() );
		runtime_assert_string_msg( requested_index > 0 && requested_index <= all_output.size(), "Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::worker_send_poses_upward(): A worker node was asked for an out-of-range job index." );
		relevant_output.push_back( all_output[ requested_index ] );
	}

	send_silent_structs( relevant_output, my_parent_ );
	if(derivedTR_.Debug.visible()) {
		derivedTR_.Debug << "Worker proc " << MPI_rank_ << " sent " << n_requests << " structures to node " << my_parent_ << "." << std::endl;
	}

}

/// @details This must be the ONLY type of request that this worker can receive at this time!
void
HierarchicalHybridJDApplication::worker_receive_request_for_sasa_summaries() const {
	HIERARCHICAL_MPI_COMMUNICATION_TYPE request(NULL_MESSAGE);
	int requesting_node(0);
	wait_for_request( requesting_node, request );
	runtime_assert( request == REQUEST_SASA_SUMMARIES_DOWNWARD && requesting_node == static_cast<int>(my_parent_) );
}

/// @brief Generate SASA metrics, and sort these for transmission up the hierarchy.
/// @details Sort order matches the order of jobsummaries.
void
HierarchicalHybridJDApplication::worker_generate_and_sort_sasa_summaries(
	utility::vector1< HierarchicalHybridJD_SASASummaryOP > & sasa_summaries,
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &jobsummaries,
	utility::vector1< core::io::silent::SilentStructOP > const &poses_from_this_worker
) const {
	runtime_assert( sasa_summaries.empty() );
	if( jobsummaries.size() == 0 ) return;
	sasa_summaries.reserve( jobsummaries.size() );
	core::chemical::ResidueTypeSet const & restypeset( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t ) ) );
	for( core::Size i(1), imax(jobsummaries.size()); i<=imax; ++i ) {
		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		poses_from_this_worker[ jobsummaries[i]->jobindex_on_originating_node() ]->fill_pose( *pose, restypeset );
		sasa_summaries.push_back( generate_sasa_summary( *pose, jobsummaries[i] ) );
	}
}


} //cyclic_peptide_predict
} //protocols

#endif //USEMPI
