// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.cc
/// @brief Wrapper for SimpleCycpepPredictApplication that allows the app to use hierarchical MPI/pthreads based job distribution.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

// Project headers:
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>

// Basic headers:
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
#include <basic/random/init_random_generator.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI" );
static basic::Tracer TR_summary( "protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication_MPI_summary" );

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Default constructor.
SimpleCycpepPredictApplication_MPI::SimpleCycpepPredictApplication_MPI():
	HierarchicalHybridJDApplication( TR, TR_summary ),
	allowed_canonicals_(),
	allowed_noncanonicals_(),
	L_alpha_comp_file_exists_(false),
	D_alpha_comp_file_exists_(false),
	L_beta_comp_file_exists_(false),
	D_beta_comp_file_exists_(false),
	comp_file_contents_L_alpha_(""),
	comp_file_contents_D_alpha_(""),
	comp_file_contents_L_beta_(""),
	comp_file_contents_D_beta_(""),
	abba_bins_("")
{}

/// @brief Constructor with options
SimpleCycpepPredictApplication_MPI::SimpleCycpepPredictApplication_MPI(
	int const MPI_rank,
	int const MPI_n_procs,
	core::scoring::ScoreFunctionCOP sfxn_in,
	core::Size const total_hierarchy_levels,
	utility::vector1 < core::Size > const & procs_per_hierarchy_level,
	utility::vector1 < core::Size > const &batchsize_per_level,
	std::string const &sort_type,
	bool const select_highest,
	core::Real const &output_fraction,
	std::string const &output_filename,
	core::Real const &lambda,
	core::Real const &kbt,
	bool const compute_rmsd_to_lowest,
	core::Real const & compute_pnear_to_lowest_fract,
	bool const compute_sasa_metrics,
	core::Size const threads_per_worker_proc //Only used in multi-threaded build.
) :
	HierarchicalHybridJDApplication(
		TR, TR_summary, MPI_rank, MPI_n_procs, sfxn_in, total_hierarchy_levels,
		procs_per_hierarchy_level, batchsize_per_level, sort_type, select_highest,
		output_fraction, output_filename, lambda, kbt, compute_rmsd_to_lowest,
		compute_pnear_to_lowest_fract, compute_sasa_metrics, threads_per_worker_proc
	),
	allowed_canonicals_(),
	allowed_noncanonicals_(),
	L_alpha_comp_file_exists_(false),
	D_alpha_comp_file_exists_(false),
	L_beta_comp_file_exists_(false),
	D_beta_comp_file_exists_(false),
	comp_file_contents_L_alpha_(""),
	comp_file_contents_D_alpha_(""),
	comp_file_contents_L_beta_(""),
	comp_file_contents_D_beta_(""),
	abba_bins_("")
{}

/// @brief Copy constructor.
SimpleCycpepPredictApplication_MPI::SimpleCycpepPredictApplication_MPI( SimpleCycpepPredictApplication_MPI const &src ) :
	HierarchicalHybridJDApplication( src ),
	allowed_canonicals_(src.allowed_canonicals_),
	allowed_noncanonicals_(src.allowed_noncanonicals_),
	L_alpha_comp_file_exists_(src.L_alpha_comp_file_exists_),
	D_alpha_comp_file_exists_(src.D_alpha_comp_file_exists_),
	L_beta_comp_file_exists_(src.L_beta_comp_file_exists_),
	D_beta_comp_file_exists_(src.D_beta_comp_file_exists_),
	comp_file_contents_L_alpha_(src.comp_file_contents_L_alpha_),
	comp_file_contents_D_alpha_(src.comp_file_contents_D_alpha_),
	comp_file_contents_L_beta_(src.comp_file_contents_L_beta_),
	comp_file_contents_D_beta_(src.comp_file_contents_D_beta_),
	abba_bins_(src.abba_bins_)
{}

/// @brief Destructor.
SimpleCycpepPredictApplication_MPI::~SimpleCycpepPredictApplication_MPI(){}

/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
HierarchicalHybridJDApplicationOP
SimpleCycpepPredictApplication_MPI::clone() const {
	return utility::pointer::make_shared< SimpleCycpepPredictApplication_MPI >(*this);
}

//////////////////////////////////////////////////////////////// PROTECTED FUNCTIONS ////////////////////////////////////////////////////////////////

/// @brief Get the protocol-specific settings.
/// @details The director reads these from disk and broadcasts them to all other nodes.  This function should be called from all nodes;
/// it figures out which behaviour it should be performing.
/// @note Pure virtual in base class; implemented here.
void
SimpleCycpepPredictApplication_MPI::get_protocol_specific_settings() {
using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( !option[basic::options::OptionKeys::cyclic_peptide::design_peptide]() ) return; //Do nothing if we're not designing.

	if( i_am_director() ) {
		//The director reads from disk and broadcasts to everyone else.
		if( option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position].user() ) {
			read_peptide_design_file( option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position](), allowed_canonicals_, allowed_noncanonicals_ );
		}
		read_file_into_string( abba_bins_, "protocol_data/generalizedKIC/bin_params/ABBA.bin_params", true /*from database*/);
		if(TR.Debug.visible()) TR.Debug << "Director read protocol_data/generalizedKIC/bin_params/ABBA.bin_params from disk." << std::endl;
		if( option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file].user() ) {
			read_file_into_string( comp_file_contents_L_alpha_, option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file](), false /*not from database*/);
			L_alpha_comp_file_exists_=true;
			if(TR.Debug.visible()) TR.Debug << "Director read " <<  option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file]() << " from disk." << std::endl;
		}
		if( option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file].user() ) {
			read_file_into_string( comp_file_contents_D_alpha_, option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file](), false /*not from database*/);
			D_alpha_comp_file_exists_=true;
			if(TR.Debug.visible()) TR.Debug << "Director read " <<  option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file]() << " from disk." << std::endl;
		}
		if( option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file].user() ) {
			read_file_into_string( comp_file_contents_L_beta_, option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file](), false /*not from database*/);
			L_beta_comp_file_exists_=true;
			if(TR.Debug.visible()) TR.Debug << "Director read " <<  option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file]() << " from disk." << std::endl;
		}
		if( option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file].user() ) {
			read_file_into_string( comp_file_contents_D_beta_, option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file](), false /*not from database*/);
			D_beta_comp_file_exists_=true;
			if(TR.Debug.visible()) TR.Debug << "Director read " <<  option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file]() << " from disk." << std::endl;
		}
	}

	//Everyone participates in the broadcast (sending and/or receiving):
	if( option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position].user() ) {
		broadcast_res_list( allowed_canonicals_ );
		broadcast_res_list( allowed_noncanonicals_ );
	}
	broadcast_string_from_director( comp_file_contents_L_alpha_ );
	broadcast_string_from_director( comp_file_contents_D_alpha_ );
	broadcast_string_from_director( comp_file_contents_L_beta_ );
	broadcast_string_from_director( comp_file_contents_D_beta_ );
	broadcast_string_from_director( abba_bins_ );

	if( !i_am_director() ) {
		if( option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file].user() ) L_alpha_comp_file_exists_=true;
		if( option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file].user() ) D_alpha_comp_file_exists_=true;
		if( option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file].user() ) L_beta_comp_file_exists_=true;
		if( option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file].user() ) D_beta_comp_file_exists_=true;
	}

	if(TR.Debug.visible() ) {
		if( option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position].user() ) {
			TR.Debug << "\nProc " << MPI_rank() << " allowed canonicals:\n";
			for(std::map<core::Size, utility::vector1 < std::string > >::const_iterator it=allowed_canonicals_.begin(); it!=allowed_canonicals_.end(); ++it) {
				TR.Debug << "Proc" << MPI_rank() << "\t" << it->first << "\t" << it->second.size() << "\t";
				for(core::Size i=1, imax=it->second.size(); i<=imax; ++i) { TR.Debug << it->second[i]; if( i<imax ) TR.Debug << ","; }
				TR.Debug << "\n";
			}
			TR.Debug << "\nProc " << MPI_rank() << " allowed noncanonicals:\n";
			for(std::map<core::Size, utility::vector1 < std::string > >::const_iterator it=allowed_noncanonicals_.begin(); it!=allowed_noncanonicals_.end(); ++it) {
				TR.Debug << "Proc" << MPI_rank() << "\t" << it->first << "\t" << it->second.size() << "\t";
				for(core::Size i=1, imax=it->second.size(); i<=imax; ++i) { TR.Debug << it->second[i]; if( i<imax ) TR.Debug << ","; }
				TR.Debug << "\n";
			}
		}
		//TR.Debug << "\nProc " << MPI_rank_ << " comp_file_contents_L_alpha_ (" << ( L_alpha_comp_file_exists_ ? "provided" : "not provided" )  << "):\n" << comp_file_contents_L_alpha_ << "\n";
		//TR.Debug << "\nProc " << MPI_rank_ << " comp_file_contents_D_alpha_ (" << ( D_alpha_comp_file_exists_ ? "provided" : "not provided" )  << "):\n" << comp_file_contents_D_alpha_ << "\n";
		//TR.Debug << "\nProc " << MPI_rank_ << " comp_file_contents_L_beta_ (" << ( L_beta_comp_file_exists_ ? "provided" : "not provided" )  << "):\n" << comp_file_contents_L_beta_ << "\n";
		//TR.Debug << "\nProc " << MPI_rank_ << " comp_file_contents_D_beta_ (" << ( D_beta_comp_file_exists_ ? "provided" : "not provided" )  << "):\n" << comp_file_contents_D_beta_ << "\n";
		//TR.Debug << "\nProc " << MPI_rank_ << " abba_bins_:\n" << abba_bins_ << "\n";
		TR.Debug << std::endl;
		TR.Debug.flush();
	}
}

/// @brief Create an instance of the SimpleCycpepPredictApplication app class, and carry out N jobs on a single process.
/// @details This code is called in a single thread in multi-threaded mode, and is used in the single-threaded version too.
/// @note This implements a function that is pure virutal function in the base class.
void
SimpleCycpepPredictApplication_MPI::derived_worker_carry_out_n_jobs(
	core::Size const njobs_from_above,
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > &all_output,
	core::scoring::ScoreFunctionOP sfxn,
	core::pose::PoseCOP native,
	std::string const &sequence
) const {
	//Create and initialize the predictor:
	SimpleCycpepPredictApplicationOP predict_app( utility::pointer::make_shared< SimpleCycpepPredictApplication >(false /*prevent file read*/) );
	if( native != nullptr ) predict_app->set_native( native );
	predict_app->set_scorefxn( sfxn );
	predict_app->set_nstruct( njobs_from_above );
	predict_app->set_sequence( sequence );
	predict_app->set_silentstructure_outputlist( &all_output, &jobsummaries );
	predict_app->set_suppress_checkpoints( true );
	predict_app->set_my_rank( MPI_rank() );
	predict_app->set_already_completed_job_count( worker_job_count() );
	predict_app->set_allowed_residues_by_position( allowed_canonicals_, allowed_noncanonicals_ );
	if( L_alpha_comp_file_exists_ ) predict_app->set_L_alpha_compfile_contents( comp_file_contents_L_alpha_ );
	if( D_alpha_comp_file_exists_ ) predict_app->set_D_alpha_compfile_contents( comp_file_contents_D_alpha_ );
	if( L_beta_comp_file_exists_ ) predict_app->set_L_beta_compfile_contents( comp_file_contents_L_beta_ );
	if( D_beta_comp_file_exists_ ) predict_app->set_D_beta_compfile_contents( comp_file_contents_D_beta_ );
	predict_app->set_abba_bins_binfile_contents( abba_bins_);

	predict_app->run();
}

/// @brief Compute the RMSD between a pose and a reference pose.
/// @details Must be implemented by derived classes, since this might be done differently for
/// different classes of molecule.
core::Real
SimpleCycpepPredictApplication_MPI::derived_worker_compute_rmsd(
	core::pose::Pose const & pose,
	core::pose::Pose const & reference_pose,
	std::string const &sequence
) const {
	//Create and initialize the predictor:
	SimpleCycpepPredictApplicationOP predict_app( utility::pointer::make_shared< SimpleCycpepPredictApplication >(false /*prevent file read*/) );
	predict_app->set_sequence( sequence );
	predict_app->set_suppress_checkpoints( true );
	predict_app->set_my_rank( MPI_rank() );

	return predict_app->align_and_calculate_rmsd( *(pose.clone()), reference_pose, true );
}

} //protocols
} //cyclic_peptide_predict

#endif //USEMPI
