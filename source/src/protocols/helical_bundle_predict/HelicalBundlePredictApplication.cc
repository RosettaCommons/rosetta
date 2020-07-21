// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HelicalBundlePredictApplication.cc
/// @brief The meat-and-potatoes for the helical_bundle_predict application, used to predict structures of helical bundles
/// made from canonical or noncanonical building-blocks.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project includes:
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication.hh>

// Core includes:
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

// Protocols includes:
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/helical_bundle_predict/HBP_MoveGenerator.hh>
#include <protocols/helical_bundle_predict/HBP_HelixCoilMoveGenerator.hh>
#include <protocols/helical_bundle_predict/HBP_FinalFullatomRefinementMoveGenerator.hh>
#include <protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.hh>
#include <protocols/helical_bundle_predict/HBP_SigmoidalTemperatureScheduleGenerator.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

// Basic includes:
#include <basic/Tracer.hh>
#include <basic/options/keys/helical_bundle_predict.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

// Utility includes:
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// Numeric includes:
#include <numeric/random/random.hh>

// Uncomment the following to debug the workings of this app.
//#define DUMP_HBP_DEBUG_OUTPUT

#ifdef DUMP_HBP_DEBUG_OUTPUT
#include <core/io/pdb/pdb_writer.hh>
#endif

static basic::Tracer TR( "protocols.helical_bundle_predict.HelicalBundlePredictApplication" );


namespace protocols {
namespace helical_bundle_predict {

///////////////////// HelicalBundlePredictApplicationOptions class /////////////////////

/// @brief Constructor
/// @details Triggers read from options system!
HelicalBundlePredictApplicationOptions::HelicalBundlePredictApplicationOptions() :
	utility::VirtualBase()
{
	initialize_from_options();
}

/// @details Destructor
HelicalBundlePredictApplicationOptions::~HelicalBundlePredictApplicationOptions() = default;

/// @brief Create a copy and return smart pointer to copy.
HelicalBundlePredictApplicationOptionsOP
HelicalBundlePredictApplicationOptions::clone() const {
	return utility::pointer::make_shared< HelicalBundlePredictApplicationOptions >( *this );
}

///////////////////////////////////// PUBLIC FUNCTIONS /////////////////////////////////////

/// @brief Indicate which commandline flags are relevant (i.e. which should be listed with the --help flag).
/// @details This is a static function that must be called BEFORE devel_init().
void
HelicalBundlePredictApplicationOptions::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( basic::options::OptionKeys::in::file::fasta );
	option.add_relevant( basic::options::OptionKeys::in::file::native );
	option.add_relevant( basic::options::OptionKeys::out::nstruct );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::helix_assignment_file );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::num_steps_per_simulated_annealing_round_centroid );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::num_simulated_annealing_rounds_centroid );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::centroid_max_temperature );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::centroid_min_temperature );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::do_final_fullatom_refinement );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::fast_relax_rounds );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::find_disulfides );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::ignore_native_residues_in_rmsd );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::ignore_prediction_residues_in_rmsd );
	option.add_relevant( basic::options::OptionKeys::helical_bundle_predict::sequence_file );

	//Options only used in MPI mode:
#ifdef USEMPI //Options that are only needed in the MPI version:
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_auto_2level_distribution         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_processes_by_level               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_batchsize_by_level               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_sort_by                          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_choose_highest                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_output_fraction                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_stop_after_time                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_pnear_lambda                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_pnear_kbt                        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::compute_rmsd_to_lowest               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::compute_pnear_to_this_fract          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::compute_ensemble_sasa_metrics        );
#ifdef MULTI_THREADED //Options that are only needed in the MPI+threads version:
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::threads_per_worker                    );
#endif //ifdef MULTI_THREADED
#endif //ifdef USEMPI

}

/// @brief Read from the options system to initialize this object.
void
HelicalBundlePredictApplicationOptions::initialize_from_options() {
	static const std::string errmsg( "Error in HelicalBundlePredictApplicationOptions::initialize_from_options(): " );
	runtime_assert_string_msg( basic::options::option[basic::options::OptionKeys::helical_bundle_predict::num_simulated_annealing_rounds_centroid]() >= 0, errmsg + "The number of simulated annealing rounds in centroid mode must be greater than or equal to zero." );
	runtime_assert_string_msg( basic::options::option[basic::options::OptionKeys::helical_bundle_predict::num_steps_per_simulated_annealing_round_centroid]() > 0, errmsg + "The number of steps per simulated annealing rounds in centroid mode must be greater than zero." );

	if ( basic::options::option[basic::options::OptionKeys::in::file::native].user() ) {
		native_file_ = basic::options::option[basic::options::OptionKeys::in::file::native]();
	} else {
		native_file_ = "";
	}
	num_simulated_annealing_rounds_centroid_ = static_cast< core::Size>( basic::options::option[basic::options::OptionKeys::helical_bundle_predict::num_simulated_annealing_rounds_centroid]() );
	num_steps_per_simulated_annealing_round_centroid_ = static_cast<core::Size>( basic::options::option[basic::options::OptionKeys::helical_bundle_predict::num_steps_per_simulated_annealing_round_centroid]() );
	centroid_max_temperature_ = basic::options::option[basic::options::OptionKeys::helical_bundle_predict::centroid_max_temperature]();
	centroid_min_temperature_ = basic::options::option[basic::options::OptionKeys::helical_bundle_predict::centroid_min_temperature]();
	signed long nstruct_int( basic::options::option[ basic::options::OptionKeys::out::nstruct ]() );
	nstruct_ = (nstruct_int < 1 ? 1 : static_cast<core::Size>(nstruct_int));
	do_fullatom_refinement_ = basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::do_final_fullatom_refinement ]();
	signed long const fastrelax_rounds_int( basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::fast_relax_rounds ]() );
	fullatom_fast_relax_rounds_ = ( fastrelax_rounds_int < 1 ? 1 : static_cast<core::Size>(fastrelax_rounds_int) );
	fullatom_find_disulfides_ = basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::find_disulfides ]();

	if ( basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::ignore_native_residues_in_rmsd ].user() ) {
		set_rmsd_residues_to_ignore_native( basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::ignore_native_residues_in_rmsd ]() );
	}

	if ( basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::ignore_prediction_residues_in_rmsd ].user() ) {
		set_rmsd_residues_to_ignore_prediction( basic::options::option[ basic::options::OptionKeys::helical_bundle_predict::ignore_prediction_residues_in_rmsd ]() );
	}
}

/// @brief Set the file containing the FASTA sequence.
void
HelicalBundlePredictApplicationOptions::set_fasta_file(
	std::string const & file_in
) {
	runtime_assert_string_msg(
		!file_in.empty(),
		"Error in HelicalBundlePredictApplicationOptions::set_fasta_file(): The input filename cannot be empty!"
	);
	runtime_assert_string_msg(
		fasta_file_.empty() && sequence_file_.empty(),
		"Error in HelicalBundlePredictApplicationOptions::set_fasta_file(): A FASTA or sequence file was already set!"
	);
	fasta_file_ = file_in;
}

/// @brief Set the file containing the full-basename sequence.
void
HelicalBundlePredictApplicationOptions::set_sequence_file(
	std::string const & file_in
) {
	runtime_assert_string_msg(
		!file_in.empty(),
		"Error in HelicalBundlePredictApplicationOptions::set_sequence_file(): The input filename cannot be empty!"
	);
	runtime_assert_string_msg(
		fasta_file_.empty() && sequence_file_.empty(),
		"Error in HelicalBundlePredictApplicationOptions::set_sequence_file(): A FASTA or sequence file was already set!"
	);
	sequence_file_ = file_in;
}

/// @brief Set the file containing the helix assignments.
void
HelicalBundlePredictApplicationOptions::set_helix_assignment_file(
	std::string const & file_in
) {
	helix_assignment_file_ = file_in;
}

/// @brief Set the contents of the FASTA file.
void
HelicalBundlePredictApplicationOptions::set_fasta_file_contents(
	std::string const & contents_in
) {
	runtime_assert_string_msg(
		fasta_file_contents_.empty() && sequence_file_contents_.empty(),
		"Error in HelicalBundlePredictApplicationOptions::set_fasta_file_contents(): A FASTA or sequence file has already been read in."
	);
	fasta_file_contents_ = contents_in;
	clean_fasta_file_contents();
}

/// @brief Set the contents of the seqeunce file.
void
HelicalBundlePredictApplicationOptions::set_sequence_file_contents(
	std::string const & contents_in
) {
	runtime_assert_string_msg(
		fasta_file_contents_.empty() && sequence_file_contents_.empty(),
		"Error in HelicalBundlePredictApplicationOptions::set_sequence_file_contents(): A FASTA or sequence file has already been read in."
	);
	sequence_file_contents_ = contents_in;
}

/// @brief Set the contents of the helix assignment file.
void
HelicalBundlePredictApplicationOptions::set_helix_assignment_file_contents(
	std::string const & contents_in
) {
	helix_assignment_file_contents_ = contents_in;
}

/// @brief Get the residues to ignore in the native pose when setting up the alignment for RMSD.
/// @details Throws errors if any are zero or negative.
void
HelicalBundlePredictApplicationOptions::set_rmsd_residues_to_ignore_native(
	utility::vector1< signed long > const &input
) {
	core::Size const vectsize( input.size() );
	rmsd_residues_to_ignore_native_.resize( vectsize );
	for ( core::Size i(1); i<=vectsize; ++i ) {
		runtime_assert_string_msg( input[i] > 0, "Error in HelicalBundlePredictApplicationOptions::set_rmsd_residues_to_ignore_native(): Residue indices to ignore in the native pose for RMSDs must be strictly positive.  Could not parse " + std::to_string( input[i] ) + "!" );
		rmsd_residues_to_ignore_native_[i] = static_cast<core::Size>(input[i]);
	}
}

/// @brief Get the residues to ignore in the generated poses when setting up the alignment for RMSD.
/// @details Throws errors if any are zero or negative.
void
HelicalBundlePredictApplicationOptions::set_rmsd_residues_to_ignore_prediction(
	utility::vector1< signed long > const &input
) {
	core::Size const vectsize( input.size() );
	rmsd_residues_to_ignore_prediction_.resize( vectsize );
	for ( core::Size i(1); i<=vectsize; ++i ) {
		runtime_assert_string_msg( input[i] > 0, "Error in HelicalBundlePredictApplicationOptions::set_rmsd_residues_to_ignore_prediction(): Residue indices to ignore in the generated poses for RMSDs must be strictly positive.  Could not parse " + std::to_string( input[i] ) + "!" );
		rmsd_residues_to_ignore_prediction_[i] = static_cast<core::Size>(input[i]);
	}
}

/// @brief Given input filenames, read the files.
/// @details INVOLVES READS FROM DISK!  WARNING!
void
HelicalBundlePredictApplicationOptions::read_inputs() {
	if ( !fasta_file_.empty() ) {
		read_fasta();
	} else {
		read_sequence_file();
	}
	read_helix_assignments();
}

///////////////////////////////////// PRIVATE FUNCTIONS /////////////////////////////////////

/// @brief Read a FASTA file from disk.
void
HelicalBundlePredictApplicationOptions::read_fasta() {
	runtime_assert_string_msg( !fasta_file_.empty(), "Error in protocols::helical_bundle_predict::HelicalBundlePredictApplication::read_fasta(): An empty FASTA filename was provided." );
	set_fasta_file_contents( utility::file_contents( fasta_file_ ) );
}

/// @brief Read a sequence file from disk.
void
HelicalBundlePredictApplicationOptions::read_sequence_file() {
	runtime_assert_string_msg( !sequence_file_.empty(), "Error in protocols::helical_bundle_predict::HelicalBundlePredictApplication::read_sequence_file(): An empty sequence file name was provided." );
	set_sequence_file_contents( utility::file_contents( sequence_file_ ) );
}

/// @brief Given a set of characters, find the first instance of any of them in a string
/// and return the (zero-based) index of that character.
core::Size
HelicalBundlePredictApplicationOptions::findchar(
	std::string const & curstring,
	utility::vector1< char > const & chars
) const {
	debug_assert( chars.size() != 0 );
	core::Size returnval(0);
	bool first( true );

	for ( char const curchar : chars ) {
		core::Size const curpos( curstring.find( curchar ) );
		if ( first || curpos < returnval ) {
			first = false;
			returnval = curpos;
		}
	}
	return returnval;
}

/// @brief Given FASTA file contents, remove comment lines.
void
HelicalBundlePredictApplicationOptions::clean_fasta_file_contents() {
	if ( fasta_file_contents_.empty() ) return; //Do nothing if no contents.
	utility::vector1< std::string > const lines( utility::split_by_newlines( fasta_file_contents_ ) );
	std::stringstream outstream;
	for ( std::string const & line : lines ) {
		std::string const linecut( line.substr( 0, findchar( line, { ';', '>' } ) ) );
		if ( !linecut.empty() ) {
			outstream << linecut << "\n";
		}
	}
	fasta_file_contents_ = outstream.str();
	TR << "Trimmed FASTA file contents to:\n" << fasta_file_contents_ << std::endl;
}

/// @brief Read a helix assignemnt file from disk.
void
HelicalBundlePredictApplicationOptions::read_helix_assignments() {
	runtime_assert_string_msg( !helix_assignment_file_.empty(), "Error in protocols::helical_bundle_predict::HelicalBundlePredictApplication::read_helix_assignemnts(): An empty helix assignment file filename was provided." );
	helix_assignment_file_contents_ = utility::file_contents( helix_assignment_file_ );
}

///////////////////// HelicalBundlePredictApplication class /////////////////////

/// @brief Options constructor.
/// @details Stores input options directly; doesn't clone.
/// @note Triggers read from disk to set up move generator!
HelicalBundlePredictApplication::HelicalBundlePredictApplication(
	HelicalBundlePredictApplicationOptionsCOP options_in
) :
	options_(options_in),
	pose_(nullptr),
	centroid_move_generator_( create_centroid_move_generator() ),
	centroid_temperature_generator_( utility::pointer::make_shared< HBP_SigmoidalTemperatureScheduleGenerator >( options_in->centroid_max_temperature(), options_in->centroid_min_temperature() ) ),
	centroid_sfxn_( create_centroid_scorefunction() ),
	fullatom_sfxn_( create_fullatom_scorefunction() ),
	final_fullatom_refinement_move_generator_( create_final_fullatom_refinement_move_generator() ),
	nstruct_(options_in->nstruct())
#ifdef USEMPI
	,
	my_rank_( 0 ),
	already_completed_job_count_(0),
	jobsummaries_( nullptr ),
	all_output_(nullptr)
#endif
{
	set_up_centroid_move_generator();
	runtime_assert(options_ != nullptr);
}

/// @brief Options + move generator constructor.
/// @details Stores input options directly; doesn't clone.
/// @note Avoids read from disk by using a move generator that was already set up.  The input move
/// generator and the input scorefunctions are used directly (not cloned).
HelicalBundlePredictApplication::HelicalBundlePredictApplication(
	HelicalBundlePredictApplicationOptionsCOP options_in,
	HBP_MoveGeneratorOP centroid_move_generator_in,
	core::scoring::ScoreFunctionOP centroid_sfxn_in,
	core::scoring::ScoreFunctionOP fullatom_sfxn_in
) :
	options_(options_in),
	native_pose_(nullptr),
	pose_(nullptr),
	centroid_move_generator_(centroid_move_generator_in),
	centroid_temperature_generator_( utility::pointer::make_shared< HBP_SigmoidalTemperatureScheduleGenerator >( options_in->centroid_max_temperature(), options_in->centroid_min_temperature() ) ),
	centroid_sfxn_( centroid_sfxn_in ),
	fullatom_sfxn_( fullatom_sfxn_in ),
	final_fullatom_refinement_move_generator_( create_final_fullatom_refinement_move_generator() ),
	nstruct_(options_in->nstruct())
#ifdef USEMPI
	,
	my_rank_(0),
	already_completed_job_count_(0),
	jobsummaries_( nullptr ),
	all_output_(nullptr)
#endif
{
	static std::string const errmsg("Error in constructor for HelicalBundlePredictApplication: ");
	runtime_assert_string_msg( centroid_move_generator_ != nullptr, errmsg + "A null pointer was provided in lieu of a centroid move generator." );
	runtime_assert_string_msg( centroid_sfxn_ != nullptr, errmsg + "A null pointer was provided in lieu of a centroid scoring function.");
	runtime_assert_string_msg( fullatom_sfxn_ != nullptr, errmsg + "A null pointer was provided in lieu of a fullatom scoring function.");
}

/// @brief Copy constructor.
HelicalBundlePredictApplication::HelicalBundlePredictApplication(
	HelicalBundlePredictApplication const &src
) :
	VirtualBase( src ),
	options_(src.options_->clone()),
	native_pose_(nullptr),
	pose_( src.pose_ == nullptr ? nullptr : src.pose_->clone() ),
	centroid_move_generator_( src.centroid_move_generator_ == nullptr ? nullptr : src.centroid_move_generator_->clone() ),
	centroid_temperature_generator_( src.centroid_temperature_generator_ == nullptr ? nullptr : src.centroid_temperature_generator_ ),
	centroid_sfxn_( src.centroid_sfxn_ == nullptr ? nullptr : src.centroid_sfxn_->clone() ),
	fullatom_sfxn_( src.fullatom_sfxn_ == nullptr  ? nullptr : src.fullatom_sfxn_->clone() ),
	final_fullatom_refinement_move_generator_( src.final_fullatom_refinement_move_generator_ == nullptr ? nullptr : src.final_fullatom_refinement_move_generator_->clone() ),
	nstruct_( src.nstruct_ )
#ifdef USEMPI
	,
	my_rank_( src.my_rank_ ),
	already_completed_job_count_( src.already_completed_job_count_ ),
	jobsummaries_( src.jobsummaries_ ),
	all_output_( src.all_output_ )
#endif
{
	static std::string const errmsg("Error in copy constructor for HelicalBundlePredictApplication: ");
	runtime_assert_string_msg( centroid_move_generator_ != nullptr, errmsg + "A null pointer was provided in lieu of a centroid move generator." );
	runtime_assert_string_msg( centroid_sfxn_ != nullptr, errmsg + "A null pointer was provided in lieu of a centroid scoring function.");
	runtime_assert_string_msg( fullatom_sfxn_ != nullptr, errmsg + "A null pointer was provided in lieu of a fullatom scoring function.");
}

/// @brief Destructor
HelicalBundlePredictApplication::~HelicalBundlePredictApplication(){}

/// @brief Clone operator: copy this object and return a smart pointer to the copy.
HelicalBundlePredictApplicationOP
HelicalBundlePredictApplication::clone() const {
	return HelicalBundlePredictApplicationOP( utility::pointer::make_shared< HelicalBundlePredictApplication >( *this ) );
}


///////////////////////////////////// PUBLIC FUNCTIONS /////////////////////////////////////

/// @brief Actually run the application and produce output.
void
HelicalBundlePredictApplication::run() {
	TR << "Starting helical_bundle_predict application." << std::endl;

	if ( !options_->native_file().empty() && native_pose_ == nullptr ) load_native_pose_from_disk();

	for ( core::Size irepeat(1); irepeat <= nstruct_; ++irepeat ) {
		if ( !options_->fasta_file_contents().empty() ) {
			pose_ = make_pose_from_fasta_contents();
		} else {
			pose_ = make_pose_from_sequence_file_contents();
		}

		check_ignore_residues_reasonable( options_->rmsd_residues_to_ignore_native(), native_pose_ );
		check_ignore_residues_reasonable( options_->rmsd_residues_to_ignore_prediction(), pose_ );

#ifdef DUMP_HBP_DEBUG_OUTPUT
		pose_->dump_pdb( "START.pdb" ); //Delete me.
#endif

		do_simulated_annealing(
			pose_,
			options_->num_simulated_annealing_rounds_centroid(),
			options_->num_steps_per_simulated_annealing_round_centroid(),
			centroid_move_generator_,
			centroid_temperature_generator_,
			*centroid_sfxn_
		);

		// Do fullatom refinement, if that's what we're supposed to do.
		if ( options_->do_fullatom_refinement() ) {
			//First, switch to fullatom mode:
			core::util::switch_to_residue_type_set( *pose_, core::chemical::FULL_ATOM_t, false, false, false );

			//Next, do the actual fullatom refinement.  This is one step, not a Monte Carlo trajectory:
			do_final_fullatom_refinement( pose_ );
		}

		core::Real rmsd(0.0);
		if ( native_pose_ != nullptr ) {
			rmsd = align_to_native_pose( *pose_ );
			TR << "Energy is " << pose_->energies().total_energy() << "; RMSD to native is " << rmsd << "." << std::endl;
		}

#ifdef USEMPI
		if( jobsummaries_ != nullptr && all_output_ != nullptr ) {
			output_to_silent_list( *pose_, *jobsummaries_, *all_output_, irepeat, native_pose_ != nullptr, rmsd );
		} else {
			utility_exit(); //TODO proper output.
		}
#else //if !USEMPI

#ifdef DUMP_HBP_DEBUG_OUTPUT
		pose_->dump_pdb( "END.pdb" );
#endif

		if ( silent_output_ ) {
			char outfile[1024];
			sprintf(outfile, "%s%s.%s", outfile_prefix_.c_str(), outfile_suffix_.c_str(), outfile_extension_.c_str() );

			core::io::silent::SilentFileOptions opts;
			core::io::silent::SilentStructOP ss( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts ) );
			char tag[1024];
			sprintf(tag, "%s%06lu%s", outfile_prefix_.c_str(), static_cast<unsigned long>(irepeat), outfile_suffix_.c_str() );
			ss->fill_struct( *pose_, std::string(tag) );
			if ( native_pose_ != nullptr ) ss->add_energy( "RMSD", rmsd ); //Add the RMSD to the energy to be written out in the silent file.
			core::io::silent::SilentFileDataOP silent_file ( utility::pointer::make_shared< core::io::silent::SilentFileData >( opts ) );
			silent_file->set_filename( std::string( outfile ) );
			silent_file->write_silent_struct( *ss, std::string( outfile ) );
		} else {
			char outfile[1024];
			sprintf(outfile, "%s%06lu%s.%s", outfile_prefix_.c_str(), irepeat, outfile_suffix_.c_str(), outfile_extension_.c_str() );
			pose_->dump_file( std::string(outfile) );
		}

#endif //USEMPI

	}

	TR << "Terminating helical_bundle_predict application." << std::endl;
}

/// @brief Set the native pose.
/// @details Does not clone the input; sets owning pointer directly.
void
HelicalBundlePredictApplication::set_native(
	core::pose::PoseCOP native
) {
	native_pose_ = native;
}

/// @brief Given a pose, align it to the native pose.
/// @details Throws an error if there's a mismatch between the pose lengths or mainchain atom counts.
/// @returns RMSD to native.
core::Real
HelicalBundlePredictApplication::align_to_native_pose(
	core::pose::Pose & pose
) const {
	static const std::string errmsg( "Error in HelicalBundlePredictApplication::align_to_native_pose(): " );
	core::Size const nres( pose.total_residue() );

	runtime_assert_string_msg( native_pose_ != nullptr && native_pose_->total_atoms() > 0, errmsg + "The native pose is null or empty!" );
	core::Size const nres_native( native_pose_->total_residue() );

	core::id::AtomID_Map< core::id::AtomID > amap;

	core::pose::initialize_atomid_map(amap, pose, core::id::AtomID::BOGUS_ATOM_ID());

	//Generate list of residues to use in native.
	utility::vector1< core::Size > alignment_residues_native;
	alignment_residues_native.reserve( nres_native );
	utility::vector1< core::Size > const & ignore_res_native( options_->rmsd_residues_to_ignore_native() );
	for ( core::Size ir(1); ir<=nres_native; ++ir ) {
		if ( ignore_res_native.has_value(ir) ) continue;
		if ( !(native_pose_->residue_type(ir).is_polymer()) ) continue;
		alignment_residues_native.push_back( ir );
	}

	//Generate list of residues to use in generated.
	utility::vector1< core::Size > alignment_residues_prediction;
	alignment_residues_prediction.reserve( nres );
	utility::vector1< core::Size > const & ignore_res_prediction( options_->rmsd_residues_to_ignore_prediction() );
	for ( core::Size ir(1); ir<=nres; ++ir ) {
		if ( ignore_res_prediction.has_value(ir) ) continue;
		if ( !(pose.residue_type(ir).is_polymer()) ) continue;
		alignment_residues_prediction.push_back( ir );
	}

	core::Size const n_alignment_residues( alignment_residues_native.size() );
	runtime_assert_string_msg( n_alignment_residues == alignment_residues_prediction.size(), errmsg + "The number of alignment residues in the native pose (" + std::to_string(n_alignment_residues) + ") does not match the number in the alignment pose (" + std::to_string(alignment_residues_prediction.size()) + ")." );
	runtime_assert_string_msg( n_alignment_residues >= 3, errmsg + "There were only " + std::to_string(n_alignment_residues) + " alignment residues.  At least 3 are needed." );

	for ( core::Size ir=1; ir<=n_alignment_residues; ++ir ) {
		core::Size const native_resindex( alignment_residues_native[ir] );
		core::Size const prediction_resindex( alignment_residues_prediction[ir] );

		core::chemical::ResidueType const & restype_native( native_pose_->residue_type(native_resindex) );
		debug_assert ( restype_native.is_polymer() ); //Should be guaranteed true at this point.

		core::chemical::ResidueType const & restype_pose( pose.residue_type(prediction_resindex) );
		core::Size const n_mainchain_atoms_native( restype_native.mainchain_atoms().size() );
		runtime_assert_string_msg( n_mainchain_atoms_native == restype_pose.mainchain_atoms().size(), errmsg + "The number of mainchain atoms in native pose residue " + restype_native.name3() + std::to_string(native_resindex) + " (" + std::to_string(n_mainchain_atoms_native) + " atoms) does not match the number of mainchain atoms in alignment pose residue " + restype_pose.name3() + std::to_string(prediction_resindex) + " (" + std::to_string( restype_pose.mainchain_atoms().size() ) + " atoms)." );

		TR << "Aligning native residue " << native_resindex << " to alignment pose residue " << prediction_resindex << "." << std::endl;

		for ( core::Size ia(1); ia<=n_mainchain_atoms_native; ++ia ) { //Loop through all mainchain heavyatoms (including atoms coming off mainchain that are not sidechain atoms, like peptide "O").
			core::Size const native_atindex( restype_native.mainchain_atoms()[ia] );
			core::Size const prediction_atindex( restype_pose.mainchain_atoms()[ia] );

			if ( restype_native.atom_is_hydrogen( native_atindex ) ) continue;
			runtime_assert( !restype_pose.atom_is_hydrogen(prediction_atindex) );

			amap[ core::id::AtomID(prediction_atindex,prediction_resindex) ] = core::id::AtomID(native_atindex,native_resindex);
		}
	}

	return core::scoring::superimpose_pose( pose, *native_pose_, amap ); //Superimpose the pose and return the RMSD.
}

/// @brief Create a new move generator.  Static, so this can be called from other classes (e.g.
/// HelicalBundlePredictApplication_MPI.)
/// @details Triggers read from disk!
HBP_MoveGeneratorOP
HelicalBundlePredictApplication::create_centroid_move_generator() {
	return utility::pointer::make_shared< HBP_HelixCoilMoveGenerator >();
}

/// @brief Create the move generator used for the final fullatom refinement step.
/// @details Not static, since it depends on the options_ object.
/// @note In its current form, this should not trigger a read from disk.
HBP_MoveGeneratorOP
HelicalBundlePredictApplication::create_final_fullatom_refinement_move_generator() const {
	return utility::pointer::make_shared< HBP_FinalFullatomRefinementMoveGenerator >( options_->fullatom_fast_relax_rounds(), options_->fullatom_find_disulfides(), fullatom_sfxn_ );
}

/// @brief Create the scorefunction used during centroid mode.
/// @details Reads from disk!  Do not use repeatedly!  Store the result rather than regenerating it!
core::scoring::ScoreFunctionOP
HelicalBundlePredictApplication::create_centroid_scorefunction() {
	return core::scoring::ScoreFunctionFactory::create_score_function( "score3.wts" );
}

/// @brief Create the scorefunction used during fullatom mode.
/// @details Reads from disk!  Do not use repeatedly!  Store the result rather than regenerating it!
/// @note This returns whatever the current default fullatom scorefunction is, currently.  This function is here to make it easy
/// to hard-code a specialized scorefunction in the future if necessary.
core::scoring::ScoreFunctionOP
HelicalBundlePredictApplication::create_fullatom_scorefunction() {
	return core::scoring::get_score_function( true ); //Gets whatever the current default fullatom scorefunction is.
}

#ifdef USEMPI

/// @brief Set the job summary list and full structure list to which this protocol should write output.
/// @details This is a rare instance in which using raw pointers is appropriate.  DO NOT EMULATE THIS UNLESS YOU KNOW WHAT YOU'RE DOING.
void
HelicalBundlePredictApplication::set_output(
	utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > * jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > * all_output
) {
	jobsummaries_ = jobsummaries;
	all_output_ = all_output;
}

/// @brief Write output, not to disk, but to a list of binary silent structures in memory.
/// @details This also writes a summary of the job (basically, energy + RMSD + jobindex) to a job summary list.
void
HelicalBundlePredictApplication::output_to_silent_list(
	core::pose::Pose const &pose,
	utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > & jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > & all_output,
	core::Size const repeat_index,
	bool const include_rmsd,
	core::Real const & rmsd
) const {
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentStructOP ss( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts ) );
	char tag[512];

	if ( my_rank_ > 0 ) {
		sprintf(tag, "result_proc%04lu_%04lu", static_cast<unsigned long>(my_rank_), static_cast<unsigned long>( repeat_index + already_completed_job_count_ ) );
	} else {
		sprintf(tag, "result_%04lu", static_cast<unsigned long>(repeat_index) );
	}
	ss->fill_struct( pose, std::string(tag) );
	if ( include_rmsd ) ss->add_energy( "RMSD", rmsd ); //Add the RMSD to the energy to be written out in the silent file.

	all_output.push_back(ss);
	core::Size curjob( jobsummaries.size() + 1 );
	jobsummaries.push_back(
		utility::pointer::make_shared< protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummary >(
			my_rank_, curjob, pose.energies().total_energy(), (include_rmsd ? rmsd : 0), 0, 0
		)
	);
}

#else // Output setters for non-MPI build:

/// @brief Set the output format.  (Automatically sets extension.)
/// @details This one sets silent output to false.
void
HelicalBundlePredictApplication::set_output_format(
	core::import_pose::FileType const type
) {
	runtime_assert( type < core::import_pose::end_of_filetype_list && static_cast< core::Size >( type ) > 0 );
	runtime_assert_string_msg( type != core::import_pose::SRLZ_file, "Error in HelicalBundlePredictApplication::set_output_format(): SRLZ file output is not supported." );
	output_filetype_ = type;
	outfile_extension_ = core::import_pose::extension_from_filetype(type);
	silent_output_ = false;
}

/// @brief Indicate that we're using silent output.  This overrides set_ouput_format(), and
/// sets the extension automatically.
void
HelicalBundlePredictApplication::set_silent_output() {
	silent_output_ = true;
	outfile_extension_ = "silent";
}

/// @brief Set the prefix and suffix for output.
void
HelicalBundlePredictApplication::set_output_prefix_and_suffix(
	std::string const & prefix,
	std::string const & suffix
) {
	outfile_prefix_ = prefix;
	outfile_suffix_ = suffix;
}

#endif //USEMPI

///////////////////////////////////// PRIVATE FUNCTIONS /////////////////////////////////////

/// @brief Read the native pose in from disk, based on the native_file_ set in the options_ object.
/// @details READS FROM DISK!  This can be avoided by passing in a pointer to the native pose.
void
HelicalBundlePredictApplication::load_native_pose_from_disk() {
	std::string const & native_file( options_->native_file() );
	if ( native_file.empty() ) return; //Do nothing if there's no native file.
	TR << "Reading native pose from " << native_file << "." << std::endl;
	native_pose_ = core::import_pose::pose_from_file(native_file);
	runtime_assert_string_msg( native_pose_->total_residue() > 0, "Error in HelicalBundlePredictApplication::load_native_pose_from_disk(): Failed to read file " + native_file + ".  Generated pose is empty." );
	TR << "Successfully read " << native_pose_->total_residue() << "-residue pose from " << native_file << "." << std::endl;
}

/// @brief Given the helix assignment file's contents, set up the HBP_HelixCoilMoveGenerator used in centroid mode..
/// @details The helix_assignment_file_content_ variable must be populated first!
void
HelicalBundlePredictApplication::set_up_centroid_move_generator() {
	static std::string const errmsg( "Error in HelicalBundlePredictApplicationOptions::set_up_move_generator(): " );
	runtime_assert_string_msg( options_->helix_assignment_file_contents() != "", errmsg + "The helix assignment file contents are empty." );
	runtime_assert_string_msg( centroid_move_generator_ != nullptr, errmsg + "The centroid-mode helix/coil transition move generator was not properly created." );
	HBP_HelixCoilMoveGeneratorOP move_generator( utility::pointer::dynamic_pointer_cast< HBP_HelixCoilMoveGenerator >(centroid_move_generator_) );
	runtime_assert_string_msg( move_generator != nullptr, errmsg + "The centroid-mode helix/coil transition move generator is of the wrong type.  This should not be possible.  Please contact Vikram K. Mulligan (vmulligan@flatironinstitute.org)." );

	move_generator->set_up_user_helix_assignments( options_->helix_assignment_file_contents() );
}

/// @brief Construct a pose from the contents of a FASTA file.  The conformation is set to linear at this point.
core::pose::PoseOP
HelicalBundlePredictApplication::make_pose_from_fasta_contents() const {
	static std::string const errmsg( "Error in protocols::helical_bundle_predict::HelicalBundlePredictApplication::make_pose_from_fasta_contents(): " );
	runtime_assert_string_msg( !options_->fasta_file_contents().empty(), errmsg + "An empty string was provided for the FASTA file contents." );

	core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
	core::pose::make_pose_from_sequence( *pose, utility::strip(options_->fasta_file_contents(), " \n\t"), core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ) );

	for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
		for ( core::Size itors(1), itorsmax( pose->residue(ir).mainchain_torsions().size() ); itors<=itorsmax; ++itors ) {
			pose->set_torsion( core::id::TorsionID( ir, core::id::BB , itors ), 180.0 );
		}
	}

	pose->update_residue_neighbors();

	return pose;
}

/// @brief Construct a pose from the contents of a sequence file.  The conformation is set to linear at this point.
core::pose::PoseOP
HelicalBundlePredictApplication::make_pose_from_sequence_file_contents() const {
	static std::string const errmsg( "Error in protocols::helical_bundle_predict::HelicalBundlePredictApplication::make_pose_from_sequence_file_contents(): " );
	runtime_assert_string_msg( !options_->sequence_file_contents().empty(), errmsg + "An empty string was provided for the sequence file contents." );

	core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );

	{ //Build the pose with the peptide stub mover.
		utility::vector1< std::string > resnames( utility::split_whitespace( options_->sequence_file_contents() ) );
		runtime_assert( resnames.size() > 0 ); //Should be true.
		protocols::cyclic_peptide::PeptideStubMover stubmover;
		stubmover.set_reset_mode(true);
		stubmover.add_residue( "Append", resnames[1], 1, true, "", 0, 0, nullptr, "" );
		if ( resnames.size() > 1 ) {
			for ( core::Size i(2), imax(resnames.size()); i<=imax; ++i ) {
				stubmover.add_residue( "Append", resnames[i], i-1, false, "", 0, 0, nullptr, "" );
			}
		}
		stubmover.apply(*pose);
	}

	//Add terminal types:
	core::pose::add_lower_terminus_type_to_pose_residue( *pose, 1 );
	core::pose::add_upper_terminus_type_to_pose_residue( *pose, pose->total_residue() );

	//Convert to centroid:
	core::util::switch_to_residue_type_set( *pose, core::chemical::CENTROID_t, false, false, false );

	for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
		for ( core::Size itors(1), itorsmax( pose->residue(ir).mainchain_torsions().size() ); itors<=itorsmax; ++itors ) {
			pose->set_torsion( core::id::TorsionID( ir, core::id::BB , itors ), 180.0 );
		}
	}

	pose->update_residue_neighbors();

	return pose;
}

/// @brief Do a single simulated annealing round (a Monte Carlo trajectory with temperature ramping from high to low, once).
void
HelicalBundlePredictApplication::do_a_simulated_annealing_round(
	core::pose::PoseOP & pose,
	core::Size num_steps,
	HBP_MoveGeneratorCOP move_generator,
	HBP_TemperatureScheduleGeneratorCOP temperature_generator,
	core::scoring::ScoreFunction const &sfxn
) const {
	protocols::rosetta_scripts::ParsedProtocolOP curmove(nullptr);
	core::Real old_energy( sfxn(*pose) );
	core::Real new_energy/*( old_energy )*/, lowest_energy( old_energy );
	//Note: I'm commenting out the initialization on new_energy to keep the clang analysis test happy, above,
	//but developers should be aware that this could create uninitialized errors down the road if anything
	//changes that causes new_energy to be read before anything is assigned to it, below.

	core::pose::PoseOP lowestE_pose( pose->clone() );

	for ( core::Size i(1); i<=num_steps; ++i ) {
		core::pose::PoseOP trial_pose( pose->clone() );

		curmove = move_generator->generate_monte_carlo_move( i, num_steps, *trial_pose );
		bool failed ( curmove == nullptr );
		if ( !failed ) { curmove->apply( *trial_pose ); }
		else {
			TR << "The current helix parameters do not generate a sensible helix.  Rejecting current move." << std::endl;
		}
		if ( failed || curmove->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
			--i; //Decrement the move count so that we try this again.
			continue; //Last move failed during application; reject it without counting the move.
		}
		core::Real const temperature( temperature_generator->calculate_current_temperature( i, num_steps ) );
		new_energy = sfxn(*trial_pose);
		if ( apply_metropolis_criterion( old_energy, new_energy, temperature ) ) { //Do we accept the move?
			//We do!
			move_generator->mark_move_accepted();
			pose = trial_pose;
			old_energy = new_energy;
			if ( new_energy < lowest_energy ) {
				lowestE_pose = trial_pose->clone();
				lowest_energy = new_energy;
			}
		} else {
			//We don't.
			move_generator->mark_move_rejected();
		}
#ifdef DUMP_HBP_DEBUG_OUTPUT
				core::io::pdb::add_to_multimodel_pdb( *lowestE_pose, "DEBUG_LOWESTE.pdb", ""  );
				core::io::pdb::add_to_multimodel_pdb( *trial_pose, "DEBUG_ATTEMPTED.pdb", ""  );
				core::io::pdb::add_to_multimodel_pdb( *pose, "DEBUG_ACCEPTED.pdb", ""  );
#endif
	}

	//Return lowest energy encountered.
	pose = lowestE_pose;
}

/// @brief Entry point into the simulated annealing.
/// @details Does multiple rounds of simulated annealing.
void
HelicalBundlePredictApplication::do_simulated_annealing(
	core::pose::PoseOP & pose,
	core::Size const num_rounds,
	core::Size const num_steps_per_round,
	HBP_MoveGeneratorOP move_generator,
	HBP_TemperatureScheduleGeneratorOP temperature_generator,
	core::scoring::ScoreFunction const &sfxn
) const {
	core::pose::PoseOP lowestE_pose( pose->clone() );
	core::Real lowest_energy( sfxn(*pose)  );
	core::Real current_energy/*( lowest_energy )*/;
	//Note: I'm commenting out the initialization on current_energy to keep the clang analysis test happy, above,
	//but developers should be aware that this could create uninitialized errors down the road if anything
	//changes that causes current_energy to be read before anything is assigned to it, below.

	temperature_generator->set_max_rounds( num_rounds );
	move_generator->set_max_rounds( num_rounds );

	for ( core::Size i(1); i<=num_rounds; ++i ) {
		temperature_generator->set_current_round( i );
		move_generator->set_current_round( i );
		do_a_simulated_annealing_round( pose, num_steps_per_round, move_generator, temperature_generator, sfxn);
		current_energy = sfxn(*pose);
		if ( current_energy < lowest_energy ) {
			lowestE_pose = pose->clone();
			lowest_energy = current_energy;
		}
	}

	// We will return the lowest energy pose:
	pose = lowestE_pose;
}

/// @brief Carry out the final full-atom refinement steps.
/// @details This is a single step, not a Monte Carlo trajectory.
void
HelicalBundlePredictApplication::do_final_fullatom_refinement(
	core::pose::PoseOP & pose
) const {
	final_fullatom_refinement_move_generator_->generate_monte_carlo_move( 1, 1, *pose )->apply( *pose );
}

/// @brief Given the old energy, the new energy, and the temperature, apply the Metropolis criterion.
bool
HelicalBundlePredictApplication::apply_metropolis_criterion(
	core::Real const &old_energy,
	core::Real const &new_energy,
	core::Real const &temperature
) const {
	if ( new_energy < old_energy ) return true;
	core::Real const prob_accept( std::exp( -(new_energy-old_energy)/temperature ) );
	return ( numeric::random::uniform() < prob_accept );
}

/// @brief Check that the list of residues to ignore in calculating RMSD is reasonable.
/// @details Throws if residues are outside the size of the posem or if residues are
/// provided but there's no pose.
void
HelicalBundlePredictApplication::check_ignore_residues_reasonable(
	utility::vector1< core::Size > const & ignore_residues,
	core::pose::PoseCOP const & pose
) const {
	std::string const errmsg( "Error in HelicalBundlePredictApplication::check_ignore_residues_reasonable(): " );
	if ( pose == nullptr ) {
		runtime_assert_string_msg( ignore_residues.empty(), errmsg + "Residues to ignore in RMSD were provided, but no native pose was provided!" );
		return;
	}
	//If we reach here, there's a pose.
	core::Size const nres( pose->total_residue() );
	for ( core::Size const i : ignore_residues ) {
		runtime_assert_string_msg( i <= nres, errmsg + "A residue to ignore was provided that's outside the residue range of the pose!  The pose is " + std::to_string(nres) + " residues long, but residue " + std::to_string(i) + " was listed as a residue to ignore." );
	}
}

} //helical_bundle_predict
} //protocols






