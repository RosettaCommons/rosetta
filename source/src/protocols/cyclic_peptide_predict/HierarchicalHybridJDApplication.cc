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
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

// Unit Headers
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>

// Package Headers
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <core/init/init.hh>
#include <numeric/conversions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/numbers.hh>
#include <utility/pointer/memory.hh>
#include <basic/random/init_random_generator.hh>
#include <core/sequence/util.hh>
#include <core/chemical/AA.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>

//numeric headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <stdio.h>

#ifdef MULTI_THREADED
#include <thread>
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
	slave_job_count_(0),
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
	kbt_(1.0)
#ifdef MULTI_THREADED
	,
	threads_per_slave_proc_(1),
	results_mutex_(),
	joblist_mutex_(),
	use_const_random_seed_(false),
	random_seed_offset_(0),
	random_seed_(11111111),
	rgtype_("")
#endif
	//TODO -- Initialize vars here.
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
	core::Real const &output_fraction,
	std::string const &output_filename,
	core::Real const &lambda,
	core::Real const &kbt,
#ifdef MULTI_THREADED
	core::Size const threads_per_slave_proc
#else
	core::Size const /*threads_per_slave_proc*/
#endif
) :
	derivedTR_(tracer),
	derivedTR_summary_(summary_tracer),
	MPI_rank_( MPI_rank ),
	MPI_n_procs_( MPI_n_procs ),
	slave_job_count_( 0 ),
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
	kbt_(kbt)
#ifdef MULTI_THREADED
	,
	threads_per_slave_proc_(threads_per_slave_proc),
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
	runtime_assert_string_msg( output_filename_ != "", "Error in constructor for HierarchicalHybridJDApplication class: The output filename cannot be empty." );
}

/// @brief Explicit virtual destructor.
///
HierarchicalHybridJDApplication::~HierarchicalHybridJDApplication() {}


/// @brief Explicit copy constructor.
///
HierarchicalHybridJDApplication::HierarchicalHybridJDApplication(
	HierarchicalHybridJDApplication const &src
) :
	derivedTR_(src.derivedTR_),
	derivedTR_summary_(src.derivedTR_summary_),
	MPI_rank_( src.MPI_rank_ ),
	MPI_n_procs_( src.MPI_n_procs_ ),
	slave_job_count_( src.slave_job_count_ ),
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
	kbt_(src.kbt_)
#ifdef MULTI_THREADED
	,
	threads_per_slave_proc_( src.threads_per_slave_proc_ ),
	results_mutex_(), //Don't copy mutexes
	joblist_mutex_(), //Don't copy mutexes.
	use_const_random_seed_(src.use_const_random_seed_),
	random_seed_offset_(src.random_seed_offset_),
	random_seed_(src.random_seed_),
	rgtype_(src.rgtype_)
#endif
	//TODO -- copy variables here.
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
/// @details On slave nodes, this creates an instance of the relevant application and runs that.  Nodes higher
/// in the communications hierarchy are just involved in sending out jobs and pulling in results.
/// @note The run() function is nonconst due to some setup steps that it performs.  It then calls const run functions
/// for emperor, master, and slave nodes.
void
HierarchicalHybridJDApplication::run() {
	debug_assert(scorefxn_ != nullptr); //Must be assigned before running.

	get_sequence(); //Emperor reads sequence from disk and broadcasts it to everyone else.
	get_native(); //Emperor reads native pose from disk and broadcasts it to everyone else.
	get_protocol_specific_settings(); //Get settings specific for the application in question.  For example, for simple_cycpep_predict, if we're doing design, then the emperor reads design settings from disk and broadcasts them to everyone else.

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Proc " << MPI_rank_ << " reports that it is configured and ready to start." << std::endl;

	if( i_am_emperor() ) {
		run_emperor();
	} else if ( i_am_slave() ) {
		//core::Size k(5); //DELETE ME -- for GDB debugging in MPI mode.
		//while(k!=4) {  } //DELETE ME -- for GDB debugging in MPI mode.
		run_slave();
	} else {
		run_intermediate_master();
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
	//TR << "RANK IN LEVEL = " << rank_in_level << std::endl; //DELETE ME
	if( hierarchy_level_ > 1 ) { //Emperors have no parent.
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
	if( hierarchy_level_ < total_hierarchy_levels_ ) { //Slaves have no children.  (Goodness -- that's a terrible-sounding statement, isn't it?  I mean that there are no nodes to which worker nodes assign work; they do the work themselves.)
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
/// @details The emperor reads this from disk and broadcasts it to all other nodes.  This function can be called from any node;
/// it figures out which behaviour it should be performing.
/// @note This function necessarily uses new and delete[].
void
HierarchicalHybridJDApplication::get_sequence() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( i_am_emperor() ) {
		runtime_assert_string_msg(
			option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user() || option[basic::options::OptionKeys::in::file::fasta].user(),
			"Error in get_sequence() in HierarchicalHybridJDApplication: No sequence file was specified.  A sequence file must be provided either with \"-cyclic_peptide:sequence_file <filename>\" or with \"-in:file:fasta <filename>\"."
		);
		if( option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user() ) {
			runtime_assert_string_msg( !option[basic::options::OptionKeys::in::file::fasta].user(), "Error in get_sequence() in HierarchicalHybridJDApplication: Only one of \"-cyclic_peptide:sequence_file\" or \"-in:file:fasta\" should be provided." );
			sequence_=utility::file_contents( option[basic::options::OptionKeys::cyclic_peptide::sequence_file]() );
			utility::trim( sequence_, " \n\t" );
		} else if( option[basic::options::OptionKeys::in::file::fasta].user() ) {
			runtime_assert_string_msg( option[basic::options::OptionKeys::in::file::fasta]().size() == 1, "Error in get_sequence() in HierarchicalHybridJDApplication: Only one FASTA file may be provided with the \"-in:file:fasta\" flag." );
			sequence_ = sequence_from_fasta( option[basic::options::OptionKeys::in::file::fasta]()[1] );
		}
	}

	int stringlen( i_am_emperor() ? sequence_.size() + 1 : 0);
	MPI_Bcast( &stringlen, 1, MPI_INT, 0, MPI_COMM_WORLD); //Brodcast the length of the sequence string.

	char * charseq( new char[stringlen] );
	if( i_am_emperor() ) sprintf( charseq, "%s", sequence_.c_str() );

	MPI_Bcast( charseq, stringlen, MPI_CHAR, 0, MPI_COMM_WORLD );
	if( i_am_emperor() ) {
		derivedTR_ << "Emperor read sequence from disk and and sent \"" << sequence_ << "\" to all other nodes." << std::endl;
		derivedTR_.flush();
	} else {
		sequence_ = std::string(charseq);
		derivedTR_.Debug << "Received sequence \"" << sequence_ << "\" in broadcast from emperor." << std::endl;
		derivedTR_.Debug.flush();
	}
	delete[] charseq;
}

/// @brief Get the native structure of the peptide we're going to predict (if the user has specified one with the -in:file:native flag).
/// @details The emperor reads this from disk and broadcasts it to all other nodes.  This function should be called from all nodes;
/// it figures out which behaviour it should be performing.
void
HierarchicalHybridJDApplication::get_native() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( !option[basic::options::OptionKeys::in::file::native].user() ) return; // Do nothing if no native state has been specified.

	if( i_am_emperor() ) { //The emperor reads the native state from disk.
		std::string natfile( option[basic::options::OptionKeys::in::file::native]() );
		derivedTR_ << "Emperor reading native pose " << natfile << " from disk." << std::endl;
		core::pose::PoseOP native( new core::pose::Pose );
		core::import_pose::pose_from_file(*native, natfile, core::import_pose::PDB_file);
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentStructOP native_ss( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( "binary", opts ) );
		native_ss->fill_struct(*native, "native");

		emperor_broadcast_silent_struct( native_ss );

		native_ = native; //Store the pose (though maybe the emperor should discard it at this point to free memory).
		derivedTR_ << "Read " << native_->size() << "-residue pose from disk and broadcasted it to all other nodes." << std::endl;
	} else { //Other nodes do this
		native_ = receive_broadcast_silent_struct_and_build_pose();
		if(derivedTR_.Debug.visible()) {
			derivedTR_.Debug << "Process " << MPI_rank_ << " received a " << native_->size() << "-residue pose from the emperor's broadcast, with sequence \"";
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
/// @note This necessarily uses new and delete for the data sent via MPI_Bcast.
void
HierarchicalHybridJDApplication::broadcast_res_list(
	std::map< core::Size, utility::vector1 < std::string > > &res_list
) const {
	if( i_am_emperor() ) { //The emperor converts the map into a string and sends it.
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
		broadcast_string_from_emperor( reslist_str );
		if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Reslist string broadcast from emperor:\n" << reslist_str << std::endl;

	} else {  //Everyone else receives the string and converts it into a map.
		res_list.clear();
		std::string conversion_string( "" );
		broadcast_string_from_emperor( conversion_string );
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


/// @brief The emperor sends a message to everyone letting them know it's time to start.
/// @details Following the go signal, slaves send requests for jobs upward.
void
HierarchicalHybridJDApplication::go_signal() const {
	char mybuf('G'); //A dummy piece of information to send.
	MPI_Bcast( &mybuf, 1, MPI_CHAR, 0, MPI_COMM_WORLD );
}

/// @brief The emperor sends a message to everyone letting them know it's time to stop.
/// @details Following the stop signal, HierarchicalHybridJDApplication::run() terminates.
void
HierarchicalHybridJDApplication::stop_signal() const {
	char mybuf('S'); //A dummy piece of information to send.
	MPI_Bcast( &mybuf, 1, MPI_CHAR, 0, MPI_COMM_WORLD );
}

/// @brief Any non-slave node can wait for a node, above or below in the hierarchy, to send it some sort of request.
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

/// @brief Send a signal to stop job distribution to a set of intermediate masters.
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

/// @brief Non-emperor nodes must call this when the emperor calls emperor_broadcast_silent_struct.
/// @details This will build a pose and return an owning pointer to it.
/// @note This function necessarily uses new and delete[].
core::pose::PoseCOP
HierarchicalHybridJDApplication::receive_broadcast_silent_struct_and_build_pose() const {
	int strlen(0);
	MPI_Bcast( &strlen, 1, MPI_INT, 0, MPI_COMM_WORLD); //Receive the length of the string.
	char * inchar( new char[strlen] );
	MPI_Bcast( inchar, strlen, MPI_CHAR, 0, MPI_COMM_WORLD); //Receive the string.

	std::string instring( inchar );

	std::istringstream in_ss( instring );
	utility::vector1 < std::string > tagvector;
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData silentfiledata( opts );
	silentfiledata.read_stream( in_ss, tagvector, true, "suppress_bitflip" );
	core::pose::PoseOP mypose( new core::pose::Pose );
	silentfiledata[silentfiledata.tags()[1]]->fill_pose(*mypose);

	delete[] inchar;

	return core::pose::PoseCOP(mypose);
}

/// @brief Convert a silent struct into a character string and send it to a node.
/// @details Intended to be used with receive_silent_struct() to allow another node to receive the broadcast.
/// Message is tagged with SILENT_STRUCT_TRANSMISSION
/// @note This function necessarily uses new and delete[].
void
HierarchicalHybridJDApplication::send_silent_structs(
	utility::vector1 < core::io::silent::SilentStructOP > const &ss_vect,
	int const target_node
) const {
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData silentfiledata( opts );
	std::stringbuf sb;
	std::ostream outstream(&sb);

	for(core::Size i=1, imax=ss_vect.size(); i<=imax; ++i) {
		silentfiledata._write_silent_struct( *(ss_vect[i]), outstream );
	}
	std::string outstring( sb.str() );

	int strlen( outstring.length() + 1 ); //The plus one is for the /0 terminal character
	char * outchar( new char[ strlen ] );
	sprintf( outchar, "%s", outstring.c_str() );

	MPI_Send( &strlen, 1, MPI_INT, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the length of the string.
	MPI_Send( outchar, strlen, MPI_CHAR, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the string.

	delete[] outchar;
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

	char * charvect( new char[ strlen ] );
	MPI_Recv( charvect, strlen, MPI_CHAR, originating_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD, &status); //Send the length of the string.
	std::string const received_string( charvect );
	results_string += received_string;
	delete[] charvect;
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
	int * procbuf( new int[sizebuf] );
	unsigned long * jobindexbuf( new unsigned long[sizebuf] );
	double * energiesbuf( new double[sizebuf] );
	double * rmsdbuf( new double[sizebuf] );
	unsigned long * hbondsbuf( new unsigned long[sizebuf] );
	unsigned long * cispepbondbuf( new unsigned long[sizebuf] );

	MPI_Recv( procbuf, sizebuf, MPI_INT, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process indices.
	MPI_Recv( jobindexbuf, sizebuf, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the originating process job indices.
	MPI_Recv( energiesbuf, sizebuf, MPI_DOUBLE, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the pose energies.
	MPI_Recv( rmsdbuf, sizebuf, MPI_DOUBLE, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the RMSD values.
	MPI_Recv( hbondsbuf, sizebuf, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the hydrogen bond counts.
	MPI_Recv( cispepbondbuf, sizebuf, MPI_UNSIGNED_LONG, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status ); //Receive the cis-peptide bond counts.

	for(core::Size i=1, imax=static_cast<core::Size>(sizebuf); i<=imax; ++i) {
		summary_list.push_back( HierarchicalHybridJD_JobResultsSummaryOP( utility::pointer::make_shared< HierarchicalHybridJD_JobResultsSummary >() ) );

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
			int * nodelist( new int[ nodelist_size ] );
			MPI_Recv( nodelist, nodelist_size, MPI_INT, originating_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD, &status );
			for( core::Size j=1; j<=static_cast<core::Size>(nodelist_size); ++j ) {
				summary_list[i]->add_MPI_rank_handling_message( nodelist[j-1] );
			}
			delete[] nodelist;
		}
		if( append_to_handler_list ) {
			summary_list[i]->add_MPI_rank_handling_message( originating_node ); //Add the originating node to the list of MPI ranks that handled this message.
		}

	}

	mergesort_jobsummaries_list( original_summary_list, summary_list, sort_type_ );

	delete[] procbuf;
	delete[] jobindexbuf;
	delete[] energiesbuf;
	delete[] rmsdbuf;
	delete[] hbondsbuf;
	delete[] cispepbondbuf;
}

/// @brief Recieve a list of job summaries.
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
	int * procbuf( new int[sizebuf] );
	unsigned long * jobindexbuf( new unsigned long[sizebuf] );
	double * energiesbuf( new double[sizebuf] );
	double * rmsdbuf( new double[sizebuf] );
	unsigned long * hbondsbuf( new unsigned long[sizebuf] );
	unsigned long * cispepbondbuf( new unsigned long[sizebuf] );

	for(core::Size i=1, imax=summary_list.size(); i<=imax; ++i) {
		procbuf[i-1] = summary_list[i]->originating_node_MPI_rank();
		jobindexbuf[i-1] = static_cast<unsigned long>( summary_list[i]->jobindex_on_originating_node() );
		energiesbuf[i-1] = static_cast<double>( summary_list[i]->pose_energy() );
		rmsdbuf[i-1] = static_cast<double>( summary_list[i]->rmsd() );
		hbondsbuf[i-1] = static_cast<unsigned long>( summary_list[i]->hbonds() );
		cispepbondbuf[i-1] = static_cast<unsigned long>( summary_list[i]->cis_peptide_bonds() );
	}

	MPI_Send( procbuf, sizebuf, MPI_INT, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process indices.
	MPI_Send( jobindexbuf, sizebuf, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the originating process job indices.
	MPI_Send( energiesbuf, sizebuf, MPI_DOUBLE, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the pose energies.
	MPI_Send( rmsdbuf, sizebuf, MPI_DOUBLE, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the RMSD values.
	MPI_Send( hbondsbuf, sizebuf, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the hydrogen bond counts.
	MPI_Send( cispepbondbuf, sizebuf, MPI_UNSIGNED_LONG, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD ); //Send the cis-peptide bond counts.

	//Send the message history lists (MPI_ranks_handling_message()):
	for(core::Size i=1, imax=summary_list.size(); i<=imax; ++i) {
		int nodelist_size(summary_list[i]->MPI_ranks_handling_message().size());
		MPI_Send( &nodelist_size, 1, MPI_INT, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
		if(nodelist_size == 0) continue;

		int * nodelist( new int[ summary_list[i]->MPI_ranks_handling_message().size() ] );
		for( core::Size j=1, jmax=summary_list[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j ) {
			nodelist[j-1] = summary_list[i]->MPI_ranks_handling_message()[j];
		}
		MPI_Send( nodelist, summary_list[i]->MPI_ranks_handling_message().size(), MPI_INT, target_node, static_cast<int>(RESULTS_SUMMARY_UPWARD), MPI_COMM_WORLD );
		delete[] nodelist;
	}

	delete[] procbuf;
	delete[] jobindexbuf;
	delete[] energiesbuf;
	delete[] rmsdbuf;
	delete[] hbondsbuf;
	delete[] cispepbondbuf;
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

/// @brief Given a string on the emperor node, send it to all nodes.
/// @details The "mystring" string is the input on the emperor node, and the output on all other nodes.
void
HierarchicalHybridJDApplication::broadcast_string_from_emperor(
	std::string &mystring
) const {
	if( i_am_emperor() ) {
		//if(TR.Debug.visible()) derivedTR_.Debug << "Broadcasting the following string from the emperor:\n" << mystring << std::endl;
		unsigned long charsize( mystring.size() + 1 );
		char * bcastchar( new char[charsize] ); // Plus one for null terminator.
		for(core::Size i=0; i < static_cast<core::Size>(charsize); ++i) bcastchar[i] = mystring.c_str()[i]; //Copy the string to the char array for broadcast.
		MPI_Bcast( &charsize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD ); //Broadcast the length of the string.
		MPI_Bcast( bcastchar, charsize, MPI_CHAR, 0, MPI_COMM_WORLD ); //Broadcast the string.
		delete[] bcastchar;
		//if(TR.Debug.visible()) derivedTR_.Debug << "String broadcast from emperor complete." << std::endl;
	} else {
		unsigned long charsize(0);
		MPI_Bcast( &charsize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD ); //Receive the length of the string.
		char * bcastchar( new char[charsize] );
		MPI_Bcast( bcastchar, charsize, MPI_CHAR, 0, MPI_COMM_WORLD ); //Receive the string.
		mystring.clear();
		mystring = bcastchar;
		delete[] bcastchar;
	}
}


/// ------------- Emperor Methods --------------------

/// @brief Is this an emperor (root) node?
/// @details The emperor is responsible for sending out and retrieving all jobs, and for all file output.
bool
HierarchicalHybridJDApplication::i_am_emperor() const {
	return MPI_rank_ == 0;
}

/// @brief The jobs done by the emperor during parallel execution.
/// @details The emperor is responsible for sending out and retrieving all jobs, and for all file output.
void
HierarchicalHybridJDApplication::run_emperor() const {
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

	clock_t start_time( clock() ); //The start of the run, for timing information.
	go_signal(); //Send signal to everyone that it's time to start.

	do {
		int requesting_node(0);
		HIERARCHICAL_MPI_COMMUNICATION_TYPE message_from_below(NULL_MESSAGE);

		wait_for_request( requesting_node, message_from_below );
		switch( message_from_below ) {
			case 	OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD:
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
			//TODO -- cover other communications here.
				break;
			default:
				//TODO
				break;
		}
	} while( n_summaries_received < n_children || n_children_done < n_children || children_reporting_job_counts < n_children );  //Keep looping until all children have reported that they're done and have sent summaries.

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Empreror node (" << MPI_rank_ << ") has received completion signals, job attempt counts, and job summaries from all children." << std::endl;

	emperor_write_summaries_to_tracer( results_summary_from_below );

	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > summary_shortlist;
	emperor_select_best_summaries( summary_shortlist, results_summary_from_below, output_fraction_, select_highest_ );  //This populates summary_shortlist.

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
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Emperor (proc " << MPI_rank_ << ") received structures from node " << requesting_node << "." << std::endl;
					break;
				default:
					break;
			}
		} while( children_that_sent_poses < children_receiving_nonzero_requests );
	}

	if(derivedTR_.visible()) derivedTR_ << "Emperor is writing final structures to disk." << std::endl;
	emperor_write_to_disk( results_collected_from_below );

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Emperor (proc " << MPI_rank_ << ") reached end of run_emperor() function." << std::endl;
	stop_signal(); //Wait for signal to everyone that it's time to stop.
	clock_t end_time( clock() );

	if(derivedTR_summary_.visible()) {
		derivedTR_summary_ << "The simple_cycpep_predict application completed " << total_jobs_attempted_by_children << " of " << nstruct
			<< " jobs in " << static_cast< core::Real >( end_time - start_time ) / static_cast<core::Real>( CLOCKS_PER_SEC )
			<< " seconds.  " << results_summary_from_below.size() << " jobs returned structures, the top "
			<< summary_shortlist.size() << " of which were written to disk." << std::endl;
		derivedTR_summary_.flush();
	}

} //run_emperor()

/// @brief Convert a silent struct into a character string and broadcast it to all nodes.
/// @details Intended to be used with receive_broadcast_silent_struct_and_build_pose() to allow all other nodes to receive the broadcast.
/// @note This function necessarily uses new and delete[].
void
HierarchicalHybridJDApplication::emperor_broadcast_silent_struct(
	core::io::silent::SilentStructOP ss
) const {
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData silentfiledata( opts );
	std::stringbuf sb;
	std::ostream outstream(&sb);

	silentfiledata._write_silent_struct( (*ss), outstream );
	std::string outstring( sb.str() );

	int strlen( outstring.length() + 1 ); //The plus one is for the /0 terminal character
	char * outchar( new char[ strlen ] );
	sprintf( outchar, "%s", outstring.c_str() );

	MPI_Bcast( &strlen, 1, MPI_INT, 0, MPI_COMM_WORLD); //Brodcast the length of the string.
	MPI_Bcast( outchar, strlen, MPI_CHAR, 0, MPI_COMM_WORLD); //Brodcast the string.

	delete[] outchar;
}

/// @brief Write out a summary of the jobs completed (node, job index on node, total energy, rmsd, handler path) to the summary tracer.
///
void
HierarchicalHybridJDApplication::emperor_write_summaries_to_tracer(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &summary_list
) const {
	if( !derivedTR_summary_.visible() ) return; //Do nothing if the tracer is off.

	derivedTR_summary_ << "Summary for " << summary_list.size() << " job(s) returning results:\n";
	bool const original_print_channel_name_setting( derivedTR_summary_.get_local_print_channel_name() );
	derivedTR_summary_.set_local_print_channel_name( false ); //Disable channel name printing from this point forward.

	derivedTR_summary_ << "MPI_slave_node\tJobindex_on_node\tRMSD\tEnergy\tHbonds\tCisPepBonds\tNode_path_to_emperor\n";

	core::Real numerator(0), denominator(0); //For calculating PNear.

	for( core::Size i=1, imax=summary_list.size(); i<=imax; ++i ) {
		derivedTR_summary_ << summary_list[i]->originating_node_MPI_rank() << "\t" << summary_list[i]->jobindex_on_originating_node() << "\t" << summary_list[i]->rmsd() << "\t" << summary_list[i]->pose_energy() << "\t" << summary_list[i]->hbonds() << "\t" << summary_list[i]->cis_peptide_bonds() << "\t";
		for(core::Size j=1, jmax=summary_list[i]->MPI_ranks_handling_message().size(); j<=jmax; ++j) {
			derivedTR_summary_ << summary_list[i]->MPI_ranks_handling_message()[j];
			if( j<jmax ) derivedTR_summary_ << ",";
		}
		derivedTR_summary_ << "\n";
		if( i % 128 == 0 ) derivedTR_summary_.flush(); //Flush every hundred and twenty-eight lines.

		//Calculations for PNear:
		if(native_) {
			core::Real const Pcurrent( std::exp( -1.0*summary_list[i]->pose_energy()/kbt() ) );
			denominator += Pcurrent;
			numerator += std::exp( -1.0 * std::pow(-summary_list[i]->rmsd() / lambda() , 2.0 ) ) * Pcurrent;
		}
	}

	if(native_ && summary_list.size() > 0 && denominator > 1e-14) {
		core::Real const PNear( numerator/denominator );
		derivedTR_summary_ << "\nPNear:\t" << PNear << "\n";
		derivedTR_summary_ << "-kB*T*ln(PNear):\t" << -1.0*kbt()*std::log( PNear ) << "\n";
		derivedTR_summary_ << "lambda:\t" << lambda() << "\n";
		derivedTR_summary_ << "kB*T:\t" << kbt() << "\n";
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
HierarchicalHybridJDApplication::emperor_select_best_summaries(
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
		derivedTR_.Debug << "MPI_slave_node\tJobindex_on_node\tRMSD\tEnergy\tHbonds\tCisPepBonds\tNode_path_to_emperor\n";

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
HierarchicalHybridJDApplication::emperor_write_to_disk(
	std::string const &output
) const {
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Writing to " << output_filename_ << "." << std::endl;
	utility::io::ozstream file(output_filename_);
	runtime_assert_string_msg( file,
		"Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::emperor_write_to_disk(): Couldn't open " + output_filename_ + " for writing!" );
	file << output << std::endl;
	file.close();
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Wrote to " << output_filename_ << "." << std::endl;
}


/// ------------- Intermediate Master Methods --------

/// @brief Is this an intermediate master node?
/// @details The masters are responsible for distributing jobs to other masters lower in the hierarchy, and/or to slaves, and for
/// collecting results from masters/slaves lower in the hierarchy and sending them up to the emperor.
bool
HierarchicalHybridJDApplication::i_am_intermediate_master() const {
	return !( i_am_emperor() || i_am_slave() );
}

/// @brief The jobs done by the intermediate masters during parallel execution.
/// @details The masters are responsible for distributing jobs to other masters lower in the hierarchy, and/or to slaves, and for
/// collecting results from masters/slaves lower in the hierarchy and sending them up to the emperor.
void
HierarchicalHybridJDApplication::run_intermediate_master() const {
	go_signal(); //Wait for signal to everyone that it's time to start.

	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > jobsummaries; //Job summaries received from children.  Sorted list.

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
			case 	OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD:
				receive_jobs_attempted_count( total_jobs_attempted_by_children, requesting_node );
				++children_reporting_job_counts;
				if( children_reporting_job_counts == n_children ) {
					//Offer the total number of jobs attempted upward:
					send_request( my_parent_, OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD);
					send_jobs_attempted_count( total_jobs_attempted_by_children,	my_parent_);
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master process " << MPI_rank_ << " reported to its parent that its children attempted " << total_jobs_attempted_by_children << " jobs." << std::endl;
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
				//Afll jobs in one slave node are done.  Check off this child as done.  When all children are done, send the same signal upward.
				++n_children_done;
				if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master process " << MPI_rank_ << " reports that child node " << requesting_node <<  " has completed." << std::endl;
				if( n_children_done == n_children ) {
					send_request( my_parent_, GIVE_COMPLETION_SIGNAL_UPWARD);
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master process " << MPI_rank_ << " reports that all child nodes have completed." << std::endl;
				}
				break;
			case OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD:
				//A slave is offering a summary of the jobs that it completed.  Signal for it to send the summary upward, and sort the list into the list that we're holding.
				if(derivedTR_.Debug.visible()) {
					derivedTR_.Debug << "Intermediate master process " << MPI_rank_ << " has been offered new job summaries from node " << requesting_node << ".  Receiving..." << std::endl;
					derivedTR_.Debug.flush();
				}
				receive_and_sort_job_summaries( jobsummaries, requesting_node, true );
				++n_summaries_received;
				if(derivedTR_.Debug.visible()) {
					derivedTR_.Debug << "Intermediate master process " << MPI_rank_ << " has received a total of " << n_summaries_received << " summary set(s) for a total of " << jobsummaries.size() << " job(s), the last from node " << requesting_node << "." << std::endl;
					derivedTR_.Debug.flush();
				}
				//Once the summary of all slave jobs is in hand, offer and then send the full list upward.
				if(n_summaries_received == n_children) {
					send_request( my_parent_, OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD);
					send_job_summaries( jobsummaries, my_parent_);
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master process " << MPI_rank_ << " sent job summaries to node " << my_parent_ << "." << std::endl;
				}
				break;
				//TODO other cases
			default:
				//TODO
				break;
		}
	} while( n_children_done < n_children || n_summaries_received < n_children || children_reporting_job_counts < n_children );

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
					if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master proc " << MPI_rank_ << " received structures from node " << requesting_node << ".  " << children_that_sent_poses << " of " << children_receiving_nonzero_requests << " have sent structures." << std::endl;
					break;
				default:
					break;
			}
		} while( children_that_sent_poses < children_receiving_nonzero_requests );
	}

	if(children_receiving_nonzero_requests > 0 ) {
		send_request( my_parent_, OFFER_NEW_POSE_BATCH_UPWARD);
		intermediate_master_send_poses_as_string_upward( results_collected_from_below, my_parent_ );
	}
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master proc " << MPI_rank_ << " sent structures to node " << my_parent_ << "." << std::endl;

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Intermediate master proc " << MPI_rank_ << " reached end of run_intermediate_master() function." << std::endl;
	stop_signal(); //Wait for signal to everyone that it's time to stop.
} //run_intermediate_master()

/// @brief Relay the jobs received from below, held as a concatenated string, up the hierarchy.
/// @details Transmission to be received with receive_pose_batch_as_string().
void
HierarchicalHybridJDApplication::intermediate_master_send_poses_as_string_upward(
	std::string const &results,
	int const target_node
) const {
	int strlen( results.length() + 1 ); //The +1 is for the null terminator ('/0').
	MPI_Send( &strlen, 1, MPI_INT, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the length of the string.

	//Irritating.  MPI_Send wants a non-const pointer, so I have to copy my const string.
	char * str_to_send( new char[strlen] );
	for(core::Size i=0; i<static_cast<core::Size>(strlen); ++i) { str_to_send[i] = results.c_str()[i]; }
	MPI_Send( str_to_send, strlen, MPI_CHAR, target_node, static_cast<int>(SILENT_STRUCT_TRANSMISSION), MPI_COMM_WORLD); //Send the string.
	delete [] str_to_send;
}


/// ------------- Slave Methods ----------------------

/// @brief Is this a slave (or worker) node?
/// @details The slaves receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
bool
HierarchicalHybridJDApplication::i_am_slave() const {
	return hierarchy_level_ == total_hierarchy_levels_;
}

/// @brief The jobs done by the slaves during parallel execution.
/// @details The slaves receive jobs from higher in the hierarchy, do them, and send results back up the hierarchy.
/// @note In multi-threaded mode, if threads_per_slave_process_ is greater than 1, then the slaves launch threads
/// to do the work.  Only the master thread does MPI calls.
void
HierarchicalHybridJDApplication::run_slave() const {
	go_signal(); //Wait for signal to everyone that it's time to start.

	//Temporary lines for GDB debugging with MPI.  DELETE LATER:
	//core::Size i(0);
	//while(i==0) {
	//}

#ifdef MULTI_THREADED
	core::Size batchcount(0); //Counter for total number of batches sent out to slave threads
#endif

	core::Size total_jobs_attempted(0);
	core::Size njobs_from_above(0);
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > jobsummaries;
	utility::vector1 < core::io::silent::SilentStructOP > all_output; //All of the output poses, in silent structure format..

	do {
		send_request( my_parent_, REQUEST_NEW_JOBS_BATCH_UPWARD );
		receive_njobs_from_above( njobs_from_above );
		if (derivedTR_.Debug.visible()) derivedTR_.Debug << "Slave process " << MPI_rank_ << " was assigned " << njobs_from_above << " job(s)." << std::endl;
		if( njobs_from_above != 0 ) {
			total_jobs_attempted += njobs_from_above;
#ifdef MULTI_THREADED
			core::Size const jobs_in_this_batch( njobs_from_above );
			++batchcount; //We're sending out another batch.
			// In multi-threaded mode, if we're launching more than one thread per slave process, do that here.
			if( threads_per_slave_proc_ > 1 ) {
				utility::vector1< std::thread > worker_threads;
				for(core::Size i(1); i<threads_per_slave_proc_; ++i) { //Note that one thread (this one) is the MPI communication thread.  There are N-1 worker threads.
					core::pose::PoseOP native_copy;
					if( native_ != nullptr ) {
						native_copy = core::pose::PoseOP( new core::pose::Pose );
						native_copy->detached_copy(*native_);
					}

					worker_threads.push_back( std::thread( std::bind(
					   &protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::slave_carry_out_njobs_in_thread,
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
				for(core::Size i(1); i<threads_per_slave_proc_; ++i) {
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
				slave_carry_out_njobs( njobs_from_above, jobsummaries, all_output );
			}
#else
			slave_carry_out_njobs( njobs_from_above, jobsummaries, all_output );
#endif
		} else {
			break;
		}
	} while( true );

	send_request( my_parent_, GIVE_COMPLETION_SIGNAL_UPWARD);
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Slave process " << MPI_rank_ << " reporting completion." << std::endl;
	send_request( my_parent_, OFFER_NEW_JOBRESULTSSUMMARY_BATCH_UPWARD);
	sort_jobsummaries_list( jobsummaries, sort_type_ ); //Slaves sort the job summaries lists, so that merge sorts can be used at all steps by above layers.
	send_job_summaries( jobsummaries, my_parent_ );
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Slave process " << MPI_rank_ << " sent job summaries for " << jobsummaries.size() << " job(s) to node " << my_parent_ << "." << std::endl;

	//Offer the total number of jobs attempted upward:
	send_request( my_parent_, OFFER_NEW_JOBS_ATTEMPTED_COUNT_UPWARD);
	send_jobs_attempted_count( total_jobs_attempted, my_parent_);
	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Slave process " << MPI_rank_ << " reported to its parent that it attempted " << total_jobs_attempted << " jobs." << std::endl;

	//Receive requests for full poses from up the hierarchy:
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > requested_jobs;
	receive_pose_requests_from_above( requested_jobs );
	if(derivedTR_.Debug.visible() && requested_jobs.size() > 0) { //Debug output only: list the jobs that this slave node was instructed to send up the hierarchy.
		derivedTR_.Debug << "Slave process " << MPI_rank_ << " received requests for the following job(s): ";
		for(core::Size i=1, imax=requested_jobs.size(); i<=imax; ++i) {
			derivedTR_.Debug << requested_jobs[i]->originating_node_MPI_rank() << "-" << requested_jobs[i]->jobindex_on_originating_node();
			if(i<imax) derivedTR_.Debug << ", ";
		}
		derivedTR_.Debug << "." << std::endl;
	}

	if(requested_jobs.size() > 0) {	//Only those slaves that received nonzero requests bother to send anything.
		send_request( my_parent_, OFFER_NEW_POSE_BATCH_UPWARD);
		slave_send_poses_upward( requested_jobs, all_output );
	}

	if(derivedTR_.Debug.visible()) derivedTR_.Debug << "Slave proc " << MPI_rank_ << " reached end of run_slave() function." << std::endl;
	derivedTR_.flush();

	stop_signal(); //Wait for signal to everyone that it's time to stop.
} //run_slave()

/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
///
void
HierarchicalHybridJDApplication::slave_carry_out_njobs(
	core::Size &njobs_from_above,
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > &all_output
) const {
	if (derivedTR_.Debug.visible()) {
		derivedTR_.Debug << "Starting " << njobs_from_above << " job(s) on slave node " << MPI_rank_ << "." << std::endl;
	}

	try {
		derived_slave_carry_out_n_jobs( njobs_from_above, jobsummaries, all_output, scorefxn_, native_, sequence_ );
		slave_job_count_ += njobs_from_above;
	} catch ( utility::excn::Exception &excn ) {
		derivedTR_.Error << "Exception in SimpleCycpepPredictApplication caught:" << std::endl;
		excn.show( derivedTR_.Error );
		derivedTR_.Error << "\nRecovering from error and continuing to next job." << std::endl;
		derivedTR_.Error.flush();
	}

	njobs_from_above=0;
} //slave_carry_out_njobs()

#ifdef MULTI_THREADED

/// @brief Decrement the job counter.  Return true if the job counter was greater than zero.
/// @details Does this with proper locking to prevent threads from stepping on one another.
bool
HierarchicalHybridJDApplication::slave_decrement_jobcount_multithreaded(
	core::Size * available_job_count,
	core::Size &already_completed_job_count,
	core::Size const jobs_in_this_batch,
	core::Size const /*thread_index*/
) const {
	std::lock_guard< std::mutex > lock( joblist_mutex_ ); //Lock the mutex
	if( (*available_job_count) == 0 ) return false;
	already_completed_job_count = slave_job_count_ + ( jobs_in_this_batch - (*available_job_count) );
	--(*available_job_count);
	//TR << "Thread " << thread_index << " decrementing available jobs.  " << already_completed_job_count << " already completed.  " << (*available_job_count) << " jobs now remain." << std::endl;
	return true;
}

/// @brief Actually carry out the jobs.  This is where SimpleCycpepPredictApplication is created and invoked.
/// @details This is the multi-threaded version, which locks the job count to decrement it, and locks the jobsummaries and all_output vectors to add to them.
void
HierarchicalHybridJDApplication::slave_carry_out_njobs_in_thread(
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
	derivedTR_.Debug << "Launching worker thread " << thread_index << " from slave process " << MPI_rank_ << std::endl;

	// Initialize the random number generators for this thread:
	int seed(1111111);
	if( !use_const_random_seed_ ) {
		seed = time(nullptr);
	}
	seed += (thread_index - 1)*MPI_n_procs_ + MPI_rank_ + (threads_per_slave_proc_ * MPI_n_procs_ * (batch_index - 1) ); //Unique seed for each proc, thread, and job batch.
	basic::random::init_random_generators( seed, rgtype_ );
	derivedTR_.Debug << "Using random seed " << seed << " for worker thread " << thread_index << " from slave process " << MPI_rank_ << std::endl;

	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > local_jobsummaries;
	utility::vector1< core::io::silent::SilentStructOP > local_all_output;
	core::Size already_completed_job_count(0);

	while( slave_decrement_jobcount_multithreaded(njobs_from_above, already_completed_job_count, jobs_in_this_batch, thread_index) ) {
		try {
			//Create and initialize the relevant application:
			derived_slave_carry_out_n_jobs(1, local_jobsummaries, local_all_output, sfxn, native, sequence);
			{ //Scope for lock
				std::lock_guard< std::mutex > lock( joblist_mutex_ ); //Lock the mutex to access slave_job_count_.
				slave_job_count_ += 1;
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

	derivedTR_.Debug << "Terminating worker thread " << thread_index << " from slave process " << MPI_rank_ << std::endl;

} //slave_carry_out_njobs_in_thread()
#endif

/// @brief Given a list of jobs that have been requested from above, send the corresponding poses up the hierarchy.
/// @details Throws an error if any jbo was completed on a different node than this slave.
void
HierarchicalHybridJDApplication::slave_send_poses_upward(
	utility::vector1< HierarchicalHybridJD_JobResultsSummaryOP > const &requested_jobs,
	utility::vector1 < core::io::silent::SilentStructOP > const &all_output
) const {

	core::Size const n_requests( requested_jobs.size() );

	utility::vector1 < core::io::silent::SilentStructOP > relevant_output;
	relevant_output.reserve( n_requests );

	for(core::Size i=1; i<=n_requests; ++i) {
		runtime_assert_string_msg( requested_jobs[i]->originating_node_MPI_rank() == MPI_rank_, "Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::slave_send_poses_upward(): A slave node was asked for a job that it didn't do!" );
		core::Size const requested_index( requested_jobs[i]->jobindex_on_originating_node() );
		runtime_assert_string_msg( requested_index > 0 && requested_index <= all_output.size(), "Error in protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication::slave_send_poses_upward(): A slave node was asked for an out-of-range job index." );
		relevant_output.push_back( all_output[ requested_index ] );
	}

	send_silent_structs( relevant_output, my_parent_ );
	if(derivedTR_.Debug.visible()) {
		derivedTR_.Debug << "Slave proc " << MPI_rank_ << " sent " << n_requests << " structures to node " << my_parent_ << "." << std::endl;
	}

}


} //cyclic_peptide_predict
} //protocols

#endif //USEMPI
