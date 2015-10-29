// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.cc
/// @brief Application-level code for the simple_cycpep_predict app.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#include "boinc_zip.h"
#endif // BOINC

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>

// Package Headers
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <numeric/conversions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>
#include <utility/file/file_sys_util.hh>

//Constraints
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>

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

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication" );

/// @brief Register the set of options that this application uses (for the help menu).
///
void
protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::sequence_file                );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::cyclic_permutations          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_rama_filter              );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rama_cutoff                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::high_hbond_weight_multiplier );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::min_final_hbonds             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::hbond_energy_cutoff          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::count_sc_hbonds              );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::checkpoint_job_identifier    );

	return;
}


namespace protocols {
namespace cyclic_peptide_predict {


/// @brief Constructor
///
SimpleCycpepPredictApplication::SimpleCycpepPredictApplication() :
	silent_out_(false),
	out_filename_("S_"),
	out_scorefilename_("default.sc"),
	sequence_file_(""),
	genkic_closure_attempts_(1),
	genkic_min_solution_count_(1),
	cyclic_permutations_(true),
	use_rama_filter_(true),
	rama_cutoff_(0.3),
	high_hbond_weight_multiplier_(10),
	min_genkic_hbonds_(3.0),
	min_final_hbonds_(0.0),
	hbond_energy_cutoff_(-0.25),
	fast_relax_rounds_(3),
	count_sc_hbonds_(false),
	native_exists_(false),
	native_filename_(""),
	nstruct_(1),
	checkpoint_job_identifier_(""),
	checkpoint_filename_("checkpoint.txt")
{
	initialize_from_options();
}


/// @brief Explicit virtual destructor.
///
SimpleCycpepPredictApplication::~SimpleCycpepPredictApplication() {}


/// @brief Explicit copy constructor.
///
SimpleCycpepPredictApplication::SimpleCycpepPredictApplication( SimpleCycpepPredictApplication const &src ) :
	silent_out_(src.silent_out_),
	out_filename_(src.out_filename_),
	out_scorefilename_(src.out_scorefilename_),
	sequence_file_(src.sequence_file_),
	genkic_closure_attempts_(src.genkic_closure_attempts_),
	genkic_min_solution_count_(src.genkic_min_solution_count_),
	cyclic_permutations_(src.cyclic_permutations_),
	use_rama_filter_(src.use_rama_filter_),
	rama_cutoff_(src.rama_cutoff_),
	high_hbond_weight_multiplier_(src.high_hbond_weight_multiplier_),
	min_genkic_hbonds_(src.min_genkic_hbonds_),
	min_final_hbonds_(src.min_final_hbonds_),
	hbond_energy_cutoff_(src.hbond_energy_cutoff_),
	fast_relax_rounds_(src.fast_relax_rounds_),
	count_sc_hbonds_(src.count_sc_hbonds_),
	native_exists_(src.native_exists_),
	native_filename_(src.native_filename_),
	nstruct_(src.nstruct_),
	checkpoint_job_identifier_(src.checkpoint_job_identifier_),
	checkpoint_filename_(src.checkpoint_filename_)
	//TODO -- copy variables here.
{}

/// @brief Initialize the application.
/// @details Initializes using the option system.
void
SimpleCycpepPredictApplication::initialize_from_options(
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//Initial checks:
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user(), "Error in simple_cycpep_predict app: the user MUST provide a sequence file using the \"-cyclic_peptide:sequence_file\" flag." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() >= 0, "Error in simple_cycpep_predict app: the number of GeneralizedKIC closure attempts (\"-cyclic_peptide:genkic_closure_attempts\" flag) cannot be negative.  (Note also that setting this to zero is risky, since GenKIC will continue to seek solutions until the minimum number of solutions is reached.)" );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() >= 0, "Error in simple_cycpep_predict app: the minimum number of GenKIC solutions (\"-cyclic_peptide:genkic_min_solution_count\" flag) cannot be negative.  (Note also that setting this to zero means no minimum.)" );
	runtime_assert_string_msg( !(option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() == 0 && option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() == 0), "Error in simple_cycpep_predict app: both the \"-cyclic_peptide:genkic_closure_attempts\" and \"-cyclic_peptide:genkic_min_solution_count\" flags were set to zero.  This would result in GenKIC looping infinitely." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds during GenKIC steps (\"-cyclic_peptide:min_genkic_hbonds\" flag) can be zero, but cannot be negative." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::min_final_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds after relaxation steps (\"-cyclic_peptide:min_final_hbonds\" flag) can be zero, but cannot be negative." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds]() > 0, "Error in simple_cycpep_predict app: the number of FastRelax rounds (\"-cyclic_peptide:fast_relax_rounds\" flag) must be greater than zero." );
	runtime_assert_string_msg( !( option[out::file::silent].user() && option[out::file::o].user() ), "Error in simple_cycpep_predict app: either silent file output (\"-out:file:silent\" flag) or PDB output (\"-out:file:o\") output may be used, but not both." );

	//Copy options to private member variables:
	sequence_file_ = option[basic::options::OptionKeys::cyclic_peptide::sequence_file]();
	genkic_closure_attempts_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() );
	genkic_min_solution_count_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() );
	cyclic_permutations_ = option[basic::options::OptionKeys::cyclic_peptide::cyclic_permutations]();
	use_rama_filter_ = option[basic::options::OptionKeys::cyclic_peptide::use_rama_filter]();
	rama_cutoff_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::rama_cutoff]() );
	high_hbond_weight_multiplier_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::high_hbond_weight_multiplier]() );
	min_genkic_hbonds_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds]() );
	min_final_hbonds_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::min_final_hbonds]() );
	hbond_energy_cutoff_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::hbond_energy_cutoff]() );
	fast_relax_rounds_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds]() );
	count_sc_hbonds_ = option[basic::options::OptionKeys::cyclic_peptide::count_sc_hbonds]();

	//Get the native, if it exists:
	if ( option[in::file::native].user() ) {
		native_exists_ = true;
		native_filename_ = option[in::file::native]();
	} else {
		native_exists_ = false;
		native_filename_ = "";
	}

	//Set up output file names:
	out_filename_ = "S_";
	out_scorefilename_ = "default.sc";
	if ( option[out::file::silent].user() ) {
		out_filename_=option[out::file::silent]();
		silent_out_=true;
		option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary"); //Force binary file output.
	} else if ( option[out::file::o].user() ) {
		out_filename_=option[out::file::o]();
		silent_out_=false;
	}
	if ( option[out::file::scorefile].user() ) { out_scorefilename_=option[out::file::scorefile](); }

	//Figure out number of structures to try to generate:
	if ( option[out::nstruct].user() ) {
		runtime_assert_string_msg( option[out::nstruct]() > 0, "Error in simple_cycpep_predict app: the \"-out:nstruct\" flag's value cannot be less than 1." );
		nstruct_ = static_cast<core::Size>(option[out::nstruct]());
	} else { nstruct_ = 1; }

	checkpoint_job_identifier_ = option[basic::options::OptionKeys::cyclic_peptide::checkpoint_job_identifier]();

	return;
}

/// @brief Actually run the application.
/// @details The initialize_from_options() function must be called before calling this.  (Called by default constructor.)
void
SimpleCycpepPredictApplication::run() const {

	//Get the scorefunction:
	core::scoring::ScoreFunctionOP sfxn_default( core::scoring::get_score_function() );
	//Create a scorefunction variant with upweighted backbone hbond terms:
	core::scoring::ScoreFunctionOP sfxn_highhbond( sfxn_default->clone() );
	sfxn_highhbond->set_weight( core::scoring::hbond_lr_bb, high_hbond_weight_multiplier_ * sfxn_default->get_weight(core::scoring::hbond_lr_bb) ); //Upweight the long-range backbone hbonds
	sfxn_highhbond->set_weight( core::scoring::hbond_sr_bb, high_hbond_weight_multiplier_ * sfxn_default->get_weight(core::scoring::hbond_sr_bb) ); //Upweight the short-range backbone hbonds
	//Create variants of the above two scorefunctions with constraint weights turned on:
	core::scoring::ScoreFunctionOP sfxn_default_cst( sfxn_default->clone() );
	core::scoring::ScoreFunctionOP sfxn_highhbond_cst( sfxn_highhbond->clone() );
	if ( sfxn_default->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn_default->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn_default->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }

	//Get the sequence that we're considering:
	utility::vector1 < std::string > resnames;
	read_sequence( sequence_file_, resnames );

	//Get the native sequence that we will compare to.
	core::pose::PoseOP native_pose;
	if ( native_exists_ ) {
		native_pose=core::pose::PoseOP(new core::pose::Pose);
		TR << "Importing native structure from " << native_filename_ << "." << std::endl;
		import_and_set_up_native ( native_filename_, native_pose, resnames.size() );
#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( *native_pose );
#endif
	} else {
		TR << "No native structure specified by the user.  No RMSD values will be calculated." << std::endl;
	}

	//Set up a filter for total number of hbonds:
	protocols::filters::CombinedFilterOP total_hbond( new protocols::filters::CombinedFilter );
	set_up_hbond_filter( total_hbond, resnames.size(), sfxn_default, static_cast<core::Real>( min_genkic_hbonds_ ) );

	//Get the checkpoint information:
	core::Size success_count(0);
	core::Size curstruct(0);
	initialize_checkpointing( curstruct, success_count );

	//EVERYTHING ABOVE THIS POINT IS DONE ONCE PER PROGRAM EXECUTION.
	++curstruct;
	for ( core::Size irepeat=curstruct, irepeat_max=nstruct_; irepeat<=irepeat_max; ++irepeat ) { //Loop nstruct times
#ifdef BOINC
		{ //Increment the model count for BOINC.
			protocols::boinc::BoincSharedMemory* shmem = protocols::boinc::Boinc::get_shmem();
			shmem->model_count = shmem->model_count + 1;
		}
#endif

		//Cyclic permutation of sequence.
		core::Size cyclic_offset(0);
		utility::vector1 < std::string > resnames_copy;
		if ( cyclic_permutations_ ) {
			cyclic_offset = do_cyclic_permutation( resnames, resnames_copy );
		}
		runtime_assert(cyclic_offset < resnames_copy.size() ); //Should be true.

		//Create the pose:
		core::pose::PoseOP pose( new core::pose::Pose );
		build_polymer(pose, resnames_copy);

		//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
		protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
		set_up_termini_mover( termini, pose );
		termini->apply(*pose);

		//Add cyclic constraints:
		add_cyclic_constraints(pose);

		//Set all omega values to 180 and randomize mainchain torsions:
		set_mainchain_torsions(pose);

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( *pose );
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( *pose );
#endif

		//Do the kinematic closure:
		bool const success( genkic_close(pose, sfxn_highhbond_cst, total_hbond) );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( *pose );
#endif

		if ( !success ) {
			TR << "Closure failed.";
			if ( irepeat < irepeat_max ) {
				TR << "  Continuing to next job." << std::endl;
				checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
				//Increment total jobs and check whether it's time to quit.
				if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
			} else {
				TR << std::endl;
			}
			TR.flush();
			continue;
		}

		//If we reach here, then closure was successful.  Time to relax the pose.

		TR << "Closure successful." << std::endl;
		protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax(sfxn_default_cst, 1) );
		(*sfxn_default_cst)(*pose);
		core::Real cur_energy( pose->energies().total_energy() );
		for ( core::Size i=1, imax=fast_relax_rounds_; i<=imax; ++i ) {
			core::pose::PoseOP pose_copy( pose->clone() );
			TR << "Applying final FastRelax, round " << i << "." << std::endl;
			frlx->apply( *pose_copy );
			(*sfxn_default_cst)(*pose_copy);
			if ( pose_copy->energies().total_energy() < cur_energy ) {
				cur_energy = pose_copy->energies().total_energy();
				(*pose) = (*pose_copy);
			}
#ifdef BOINC_GRAPHICS
			// attach boinc graphics pose observer
			protocols::boinc::Boinc::update_graphics_current( *pose );
#endif
		}

		//Undo the cyclic permutation in anticipation of re-aligning to the native:
		depermute( pose, cyclic_offset );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
#endif

		core::Real native_rmsd(0.0);

		//Score the pose before output:
		(*sfxn_default)(*pose);

		//Re-filter based on number of Hbonds (using option[min_final_hbonds]()):
		core::Real const final_hbonds( total_hbond->compute( *pose ) );
		if ( final_hbonds > -1.0*min_final_hbonds_ ) {
			TR << "Final hbond count is " << -1.0*final_hbonds << ", which is less than the minimum.  Failing job." << std::endl;
			TR.flush();
			checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
			//Increment total jobs and check whether it's time to quit.
			if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
			continue;
		}

		++success_count; //Increment the count of number of successes.

		if ( native_pose ) {
			native_rmsd = align_and_calculate_rmsd(pose, native_pose);
		}
		TR << "Result\tRMSD\tEnergy\tHbonds" << std::endl;
		TR << irepeat << "\t";
		if ( native_pose ) { TR << native_rmsd; }
		else { TR << "--"; }
		TR << "\t" << pose->energies().total_energy() << "\t" << -1.0*final_hbonds << std::endl;

		if ( silent_out_ ) {
			core::io::silent::SilentFileDataOP silent_file (new io::silent::SilentFileData );
			silent_file->set_filename( out_filename_ );
			core::io::silent::SilentStructOP ss( io::silent::SilentStructFactory::get_instance()->get_silent_struct_out() );
			char tag[512];
			sprintf(tag, "result_%04lu", static_cast<unsigned long>(irepeat) );
			ss->fill_struct( *pose, std::string(tag) );
			if ( native_pose ) ss->add_energy( "RMSD", native_rmsd ); //Add the RMSD to the energy to be written out in the silent file.
			ss->add_energy( "HBOND_COUNT", -1.0*final_hbonds ); //Add the hbond count to be written out in the silent file.
#ifdef BOINC_GRAPHICS
			protocols::boinc::Boinc::update_graphics_current( *pose );
			protocols::boinc::Boinc::update_graphics_current_ghost( *pose );
			protocols::boinc::Boinc::update_graphics_last_accepted( *pose, pose->energies().total_energy() );
			protocols::boinc::Boinc::update_graphics_low_energy( *pose, pose->energies().total_energy() );
#endif
			silent_file->write_silent_struct( *ss, out_filename_ );
		} else { //if pdb output
			char outstring[512];
			sprintf(outstring, "%s%04lu.pdb", out_filename_.c_str(), static_cast<unsigned long>(irepeat) );
			pose->dump_scored_pdb( std::string(outstring), *sfxn_default );
		}

		TR.flush();
		checkpoint( irepeat, success_count ); //This job has been attempted and has succeeded; don't repeat it.

#ifdef BOINC
		//Increment total jobs and check whether it's time to quit.
		if (protocols::boinc::Boinc::worker_is_finished( irepeat )) break;
#endif
	} //Looping through nstruct

	TR << nstruct_ << " jobs attempted.  " << success_count << " jobs returned solutions." << std::endl;
	TR.flush();

	end_checkpointing(); //Delete the checkpoint file at this point, since all jobs have completed.
	return;
}


/// @brief Actually build the geometry that we'll be working with.
///
void
SimpleCycpepPredictApplication::build_polymer(
	core::pose::PoseOP pose,
	utility::vector1<std::string> const &restypes
) const {
	using namespace protocols::cyclic_peptide;
	core::Size const nres( restypes.size() );
	runtime_assert(restypes.size() >=4 );

	TR << "Building sequence ";
	for ( core::Size i=1; i<=nres; ++i ) {
		TR << restypes[i];
		if ( i<nres ) TR << " ";
	}
	TR << "." << std::endl;

	PeptideStubMover stubmover;

	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	for ( core::Size i=1; i<=nres; ++i ) {
		stubmover.add_residue( "Append", restypes[i], 0, false, "", 1, 0, "" );
	}

	stubmover.apply(*pose);

	TR << "Build successful." << std::endl;

	return;
}


/// @brief Read a sequence (as a series of full names, separated by whitespace) and store
/// it in a string vector.
void
SimpleCycpepPredictApplication::read_sequence (
	std::string const &seqfile,
	utility::vector1 < std::string > &resnames
) const {
	using namespace utility::io;
	resnames.clear();

	izstream infile;
	infile.open( seqfile );
	runtime_assert_string_msg( infile.good(), "Error in read_sequence() in app simple_cycpep_predict:  Unable to open sequence file for read!" );

	std::string curline(""); //Buffer for current line.
	utility::vector1< std::string > lines; //Storing all lines

	//Read the file:
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.
		lines.push_back( curline );
	}
	infile.close();

	//Parse the lines:
	for ( core::Size i=1, imax=lines.size(); i<=imax; ++i ) { //Loop through all lines
		std::istringstream curline(lines[i]);
		std::string oneword("");
		while ( !curline.eof() ) {
			curline >> oneword;
			resnames.push_back( oneword );
		}
	}

	if ( TR.visible() ) {
		TR << "Parsed the following sequence:" << std::endl;
		for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
			TR << resnames[i];
			if ( i<imax ) TR << ", ";
		}
		TR << "." << std::endl;
	}

	runtime_assert_string_msg( resnames.size() >= 4, "Error in simple_cycpcp_predict app read_sequence() function!  The minimum number of residues for a cyclic peptide is 4.  (GenKIC requires three residues, plus a fourth to serve as an anchor)." );

	return;
}


/// @brief Set up the DeclareBond mover used to connect the termini.
///
void
SimpleCycpepPredictApplication::set_up_termini_mover (
	protocols::cyclic_peptide::DeclareBondOP termini,
	core::pose::PoseCOP pose,
	bool const native
) const {
	core::Size const nres(pose->n_residue());

	runtime_assert_string_msg(pose->residue(1).has_lower_connect(), "Error in simple_cycpep_predict app set_up_termini_mover() function: residue 1 does not have a LOWER_CONNECT.");
	runtime_assert_string_msg(pose->residue(nres).has_upper_connect(), "Error in simple_cycpep_predict app set_up_termini_mover() function: the final residue does not have an UPPER_CONNECT.");
	std::string firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string lastatom( pose->residue(nres).atom_name( pose->residue(nres).upper_connect_atom() ) );

	if ( native ) {
		TR << "Setting up terminal bond for the native pose between residue 1, atom " << firstatom << " and residue " << nres << ", atom " << lastatom << "." << std::endl;
	} else {
		TR << "Setting up terminal bond between residue 1, atom " << firstatom << " and residue " << nres << ", atom " << lastatom << "." << std::endl;
	}

	termini->set( nres, lastatom, 1, firstatom, false, false, 0, 0, false  );

	return;
}


/// @brief Takes a vector of residue names, chooses a random number for cyclic offset, and
/// does a cyclic permutation.
/// @details Returns the offset and stores the new string vector in resnames_copy.
core::Size
SimpleCycpepPredictApplication::do_cyclic_permutation (
	utility::vector1 <std::string> const &resnames,
	utility::vector1 <std::string> &resnames_copy
) const {
	core::Size const nname( resnames.size() );//Number of residue names
	core::Size const offset( static_cast<core::Size>(numeric::random::rg().random_range(0,resnames.size()-1)) );

	resnames_copy.clear();
	resnames_copy.resize(nname, "");
	core::Size counter(0);
	for ( core::Size i=offset+1; i<=nname; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}
	for ( core::Size i=1; i<=offset; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}

	TR << "Circularly shifted residue list by " << offset << ".  New list is: ";
	for ( core::Size i=1; i<=nname; ++i ) {
		TR << resnames_copy[i];
		if ( i<nname ) TR << ", ";
	}
	TR << std::endl;
	TR.flush();

	return offset;
}


/// @brief Imports the native pose and sets up a terminial peptide bond.
///
void
SimpleCycpepPredictApplication::import_and_set_up_native (
	std::string const &native_file,
	core::pose::PoseOP native_pose,
	core::Size const expected_residue_count
) const {
	core::import_pose::pose_from_pdb(*native_pose, native_file);
	TR << "Improrting native structure from " << native_file << "." << std::endl;
	runtime_assert_string_msg( native_pose->n_residue() == expected_residue_count, "Error in simple_cycpep_predict app!  The imported native pose has a different number of residues than the sequence provided." );

	TR << "Stripping termini from native structure." << std::endl;
	core::pose::remove_lower_terminus_type_from_pose_residue(*native_pose, 1);
	core::pose::remove_upper_terminus_type_from_pose_residue(*native_pose, expected_residue_count);

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( termini, native_pose, true );
	termini->apply(*native_pose);

	return;
}


/// @brief Function to add cyclic constraints to a pose.
///
void
SimpleCycpepPredictApplication::add_cyclic_constraints (
	core::pose::PoseOP pose
) const {
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	TR << "Setting up cyclic constraints." << std::endl;

	core::Size const nres(pose->n_residue());

	//The four atoms defining the peptide bond:
	AtomID const atom_a( pose->residue(nres).type().icoor(pose->residue(nres).upper_connect_atom()).stub_atom1().atomno(), nres );
	AtomID const atom_b( pose->residue(nres).upper_connect_atom(), nres );
	AtomID const atom_c( pose->residue(1).lower_connect_atom(), 1 );
	core::Size atom_d_index(0);
	for ( core::Size i=1, imax=pose->residue(1).n_mainchain_atoms(); i<=imax; ++i ) { //Find the atom index of the first mainchain atom with the lower_connect atom as a parent.
		if ( i == atom_c.atomno() ) continue;
		if ( pose->residue(1).type().icoor(i).stub_atom1().atomno() == atom_c.atomno() ) {
			atom_d_index=i;
			break;
		}
	}
	AtomID const atom_d( atom_d_index, 1 );

	TR << "The following four atoms define the terminal bond:" << std::endl;
	TR << "1.\tRes=" << atom_a.rsd() << "\tAtom=" << pose->residue(atom_a.rsd()).atom_name(atom_a.atomno()) << std::endl;
	TR << "2.\tRes=" << atom_b.rsd() << "\tAtom=" << pose->residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
	TR << "3.\tRes=" << atom_c.rsd() << "\tAtom=" << pose->residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
	TR << "4.\tRes=" << atom_d.rsd() << "\tAtom=" << pose->residue(atom_d.rsd()).atom_name(atom_d.atomno()) << std::endl;

	{//Peptide bond length constraint:
		FuncOP harmfunc1( new HarmonicFunc( SimpleCycpepPredictApplication_PEPBOND_LENGTH, 0.01) );
		ConstraintCOP distconst1( new AtomPairConstraint ( atom_b, atom_c, harmfunc1 ) );
		pose->add_constraint (distconst1);
	}

	{ //Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1( new CircularHarmonicFunc( numeric::constants::d::pi, 0.02) );
		ConstraintCOP dihedconst1( new DihedralConstraint ( atom_a, atom_b, atom_c, atom_d, circharmfunc1) );
		pose->add_constraint (dihedconst1);
	}

	{ //Peptide bond angle constraints:
		FuncOP circharmfunc2a( new CircularHarmonicFunc( SimpleCycpepPredictApplication_PEPBOND_C_ANGLE, 0.02) );
		FuncOP circharmfunc2b( new CircularHarmonicFunc( SimpleCycpepPredictApplication_PEPBOND_N_ANGLE, 0.02) );
		ConstraintCOP angleconst1( new AngleConstraint ( atom_a, atom_b, atom_c, circharmfunc2a) );
		ConstraintCOP angleconst2( new AngleConstraint ( atom_b, atom_c, atom_d, circharmfunc2b) );
		pose->add_constraint (angleconst1);
		pose->add_constraint (angleconst2);
	}

	TR << "Finished setting up constraints." << std::endl;

	return;
}


/// @brief Sets all omega values to 180, and randomizes mainchain torsions.
/// @details For alpha-amino acids, mainchain torsions are randomized by the Ramachandran plot.
/// For other residue types, just randomizes mainchain torsions other than peptide bonds.
void
SimpleCycpepPredictApplication::set_mainchain_torsions (
	core::pose::PoseOP pose
) const {
	TR << "Randomizing mainchain torsions." << std::endl;
	core::Size const nres(pose->n_residue());
	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues
		if ( pose->residue(i).type().is_alpha_aa() ) {
			core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran(); //Get the Rama scoring function
			core::Real phi(0.0), psi(0.0);
			//TR << "aa" << i << "=" << pose->residue_type(i).aa() << std::endl; //DELETE ME
			rama.random_phipsi_from_rama( pose->residue_type(i).aa(), phi, psi); //TODO -- use backbone_aa
			pose->set_phi(i,phi);
			pose->set_psi(i,psi);
			if ( i!=nres ) pose->set_omega(i, 180.0);
		} else { //If this is not an alpha-amino acid:
			for ( core::Size j=1, jmax=pose->residue(i).mainchain_torsions().size(); j<=jmax; ++j ) { //Loop through all mainchain torsions.
				if ( i==nres && j==jmax ) continue; //Skip the last mainchain torsion (not a DOF).
				core::Real setting(180.0);
				if ( j!=jmax ) {
					setting = numeric::random::rg().uniform()*360.0 - 180.0;
				}
				pose->set_torsion( core::id::TorsionID(i, core::id::BB, j), setting );
			}
		}
	}
	return;
}

/// @brief Set up the filters for the mainchain hydrogen bonds that will
/// be used to discard solutions with too little mainchain hydrogen bonding.
void
SimpleCycpepPredictApplication::set_up_hbond_filter(
	protocols::filters::CombinedFilterOP total_hbond,
	core::Size const nres,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const &min_hbonds
) const {
	total_hbond->set_threshold( -1.0 * min_hbonds );
	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues and add hbond counters
		protocols::protein_interface_design::filters::HbondsToResidueFilterOP hbondfilt( new protocols::protein_interface_design::filters::HbondsToResidueFilter );
		hbondfilt->set_resnum(i);
		hbondfilt->set_sidechain( count_sc_hbonds_ );
		hbondfilt->set_energy_cutoff( hbond_energy_cutoff_ );
		hbondfilt->set_partners(0);
		hbondfilt->set_scorefxn( sfxn );
		total_hbond->add_filter( hbondfilt, -0.5, false );
	}
	return;
}


/// @brief Use GeneralizedKIC to close the pose.
///
bool
SimpleCycpepPredictApplication::genkic_close(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP sfxn_highhbond,
	protocols::filters::CombinedFilterOP total_hbond
) const {
	using namespace protocols::generalized_kinematic_closure;
	using namespace protocols::generalized_kinematic_closure::selector;

	TR << "Performing GeneralizedKIC closure of loop." << std::endl;

	//Number of residues in the pose:
	core::Size const nres( pose->n_residue() );
	runtime_assert( nres >= 4 ); //Already checked at sequence load time, so should be true, but let's make sure.

	//Randomly pick one of the middle residues to be the anchor residue:
	core::Size const anchor_res( numeric::random::rg().random_range(2, nres-1) );
	core::Size const first_loop_res( anchor_res + 1 );
	core::Size const last_loop_res( anchor_res - 1 );

	//Randomly pick a residue to be the middle pivot residue.  Can't be first in loop, last in loop, or anchor res.
	core::Size middle_loop_res( numeric::random::rg().random_range(1, nres-3 ) );
	if ( middle_loop_res == last_loop_res ) { middle_loop_res += 3; }
	else if ( middle_loop_res == anchor_res ) { middle_loop_res +=2; }
	else if ( middle_loop_res == first_loop_res ) { middle_loop_res +=1; }
	if ( middle_loop_res > nres ) { middle_loop_res -= nres; }

	//Create the pre-selection mover and set options.
	protocols::rosetta_scripts::ParsedProtocolOP pp( new protocols::rosetta_scripts::ParsedProtocol );
	pp->add_mover_filter_pair( NULL, "Total_Hbonds", total_hbond );
	protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax(sfxn_highhbond, fast_relax_rounds_) );
	pp->add_mover_filter_pair( frlx, "High_Hbond_FastRelax", NULL );

	//Create the mover and set options:
	GeneralizedKICOP genkic( new GeneralizedKIC );
	genkic->set_selector_type( lowest_energy_selector );
	genkic->set_closure_attempts( genkic_closure_attempts_ );
	genkic->set_min_solution_count( genkic_min_solution_count_ );
	genkic->set_selector_scorefunction( sfxn_highhbond );
	genkic->set_preselection_mover(pp);

	//If we're using BOINC graphics, let the GenKIC mover update the graphics with a "ghost" of the current
	//conformation being sampled:
#ifdef BOINC_GRAPHICS
	genkic->set_attach_boinc_ghost_observer(true);
#endif


	//Define the loop residues:
	for ( core::Size i=first_loop_res; i<=nres; ++i ) { genkic->add_loop_residue(i); }
	for ( core::Size i=1; i<=last_loop_res; ++i ) { genkic->add_loop_residue(i); }

	//Set pivots:
	std::string at1(""), at2(""), at3("");
	if ( pose->residue(first_loop_res).type().is_alpha_aa() ) { at1="CA"; }
	else if ( pose->residue(first_loop_res).type().is_beta_aa() ) { at1="CM"; }
	else if ( pose->residue(first_loop_res).type().is_gamma_aa() ) { at1="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop start.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	if ( pose->residue(middle_loop_res).type().is_alpha_aa() ) { at2="CA"; }
	else if ( pose->residue(middle_loop_res).type().is_beta_aa() ) { at2="CM"; }
	else if ( pose->residue(middle_loop_res).type().is_gamma_aa() ) { at2="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	if ( pose->residue(last_loop_res).type().is_alpha_aa() ) { at3="CA"; }
	else if ( pose->residue(last_loop_res).type().is_beta_aa() ) { at3="CM"; }
	else if ( pose->residue(last_loop_res).type().is_gamma_aa() ) { at3="C3"; }
	else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids." ); }
	genkic->set_pivot_atoms( first_loop_res, at1, middle_loop_res, at2, last_loop_res, at3 );

	//Close the bond:
	std::string const firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string const lastatom( pose->residue(nres).atom_name( pose->residue(nres).upper_connect_atom() ) );
	genkic->close_bond( nres, lastatom, 1, firstatom, 0, "", 0, "", SimpleCycpepPredictApplication_PEPBOND_LENGTH, SimpleCycpepPredictApplication_PEPBOND_C_ANGLE/numeric::constants::d::pi*180.0, SimpleCycpepPredictApplication_PEPBOND_N_ANGLE/numeric::constants::d::pi*180.0, 180.0, false, false );

	//Add perturbers:
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( i==anchor_res ) continue; //Can't perturb the anchor residue.
		if ( pose->residue(i).type().is_alpha_aa() ) {
			genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_alpha_backbone_by_rama );
			genkic->add_residue_to_perturber_residue_list(i);
		} else {
			//TODO Randomize mainchain torsions here for beta- and gamma-amino acids.
			utility_exit_with_message( "Handling of beta- and gamma-amino acids in setup of the genKIC perturber in the simple_cycpep_predict app has not yet been written.  TODO." );
		}
	}

	//Add bump check filter:
	genkic->add_filter( protocols::generalized_kinematic_closure::filter::loop_bump_check );

	//Add rama check filters:
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( i!=first_loop_res && i!=middle_loop_res && i!=last_loop_res ) continue; //Just filter the pivots.
		if ( pose->residue(i).type().is_alpha_aa() ) {
			genkic->add_filter( protocols::generalized_kinematic_closure::filter::alpha_aa_rama_check );
			genkic->set_filter_resnum(i);
			genkic->set_filter_rama_cutoff_energy( rama_cutoff_ );
			if ( i==first_loop_res ) genkic->set_filter_attach_boinc_ghost_observer(true);
		}
	}

	//Apply the mover:
	genkic->apply( *pose );

	return genkic->last_run_successful();
}


/// @brief Given a pose that has undergone an N-residue cyclic permutation, restore
/// the original pose, without the permutation.
void
SimpleCycpepPredictApplication::depermute (
	core::pose::PoseOP pose,
	core::Size const offset
) const {

	/// 1 2 3 4 5 6 7 8
	/// 2 3 4 5 6 7 8 1
	/// 3 4 5 6 7 8 1 2
	/// 4 5 6 7 8 1 2 3

	if ( offset==0 ) return; //Do nothing if the pose was not offset.

	core::Size const nres(pose->n_residue());
	debug_assert(nres > offset);
	core::Size const old_first_res_index( nres-offset+1 );

	//TR << "nres=" << nres << " offset=" << offset << " old_first_res_index=" << old_first_res_index << std::endl; //DELETE ME

	core::pose::PoseOP newpose( new core::pose::Pose );

	for ( core::Size ir=old_first_res_index; ir<=nres; ++ir ) {
		if ( ir == old_first_res_index ) {
			newpose->append_residue_by_jump( *(pose->residue(ir).clone()), 0, "", "", true );
		} else {
			newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
		}
	}

	for ( core::Size ir=1; ir<old_first_res_index; ++ir ) {
		newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
	}

	//I don't bother to set up cyclic constraints, since we won't be doing any more minimization after calling this function.

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::cyclic_peptide::DeclareBondOP termini( new protocols::cyclic_peptide::DeclareBond );
	set_up_termini_mover( termini, newpose );
	termini->apply(*newpose);

	(*pose) = (*newpose);

	return;
}


/// @brief Align pose to native_pose, and return the RMSD between the two poses.
/// @details Assumes that the pose has already been de-permuted (i.e. the native and the pose line up).
core::Real
SimpleCycpepPredictApplication::align_and_calculate_rmsd(
	core::pose::PoseOP pose,
	core::pose::PoseCOP native_pose
) const {
	core::Size const nres( pose->n_residue() );
	debug_assert( native_pose->n_residue() == nres ); //Should be true.

	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, *pose, core::id::BOGUS_ATOM_ID);
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		for ( core::Size ia=1, iamax=native_pose->residue(ir).type().first_sidechain_atom(); ia<iamax; ++ia ) { //Loop through all mainchain heavyatoms (including atoms coming off mainchain that are not sidechain atoms, like peptide "O").
			if ( native_pose->residue(ir).type().atom_is_hydrogen(ia) ) continue;
			amap[ core::id::AtomID(ia,ir) ] = core::id::AtomID(ia,ir);
			//TR << "Adding ia=" << ia << " ir=" << ir << " to map." << std::endl; //DELETE ME
		}
	}
	return core::scoring::superimpose_pose( *pose, *native_pose, amap ); //Superimpose the pose and return the RMSD.
}

/// @brief Create a new checkpoint file.
///
void
SimpleCycpepPredictApplication::new_checkpoint_file() const {
	using namespace utility::io;

	ozstream outfile;
	outfile.open( checkpoint_filename_ );
	runtime_assert(outfile.good());
	outfile << checkpoint_job_identifier_ << std::endl;
	outfile << "LAST\t0\tSUCCESS\t0" << std::endl;
	outfile.flush();
	outfile.close();

	erase_random_seed_info();
	store_random_seed_info();

	return;
}


/// @brief Initialize checkpointing for this run.
/// @details  This function does several things.  First, it checks for an existing checkpoint
/// file.  If one exists, it checks whether the unique job name in the file matches the current
/// job.  If it does, then this job has already been attempted, and we're somewhere in the middle
/// of it.  The function reads the last attempt number and success count from the checkpoint
/// file, and returns these values.  Otherwise, it creates a new checkpoint file with the current
/// job name and returns (0,0).  If checkpointing is disabled, this function does nothing, and
/// returns (0,0).
/// @param[out] lastjob The index of the last job run.  Set to zero if checkpointing is disabled
/// or if we're creating a new checkpoint file (first job run).
/// @param[out] successes The number of successes so far.  Set to zero if checkpointing is
/// disabled or if we're creating a new checkpoint file (first job run).
void
SimpleCycpepPredictApplication::initialize_checkpointing(
	core::Size &lastjob,
	core::Size &successes
) const {
	using namespace utility::io;

	//If we're not using checkpointing, return 0,0:
	if ( checkpoint_job_identifier_ == "" ) {
		lastjob=0;
		successes=0;
		return;
	}

	//Check for a checkpoint file:
	izstream infile;
	infile.open( checkpoint_filename_ );
	if ( !infile.good() ) {
		//If the checkpoint file doesn't exist/isn't readable, then we're starting a new job and need a new checkpoint file.
		infile.close();
		new_checkpoint_file();
		lastjob=0;
		successes=0;
		return;
	}

	//If we've reached this point, then the checkpoint file IS good, and we need to read the job name and the last job ID/success count:
	std::string curline;
	infile.getline(curline);
	if ( infile.eof() || curline=="" || curline!=checkpoint_job_identifier_ ) {
		//If the checkpoint file isn't readable or is for a different job, then we're starting a new job and need a new checkpoint file.
		infile.close();
		new_checkpoint_file();
		lastjob=0;
		successes=0;
		return;
	}

	//Loop through the checkpoint file and read the job lines:
	while ( !infile.eof() ) {
		infile.getline(curline);
		std::istringstream ss(curline);
		ss >> curline;
		ss >> lastjob;
		ss >> curline;
		ss >> successes;
	}
	infile.close();

	get_random_seed_info();

	if ( TR.Debug.visible() ) {
		TR.Debug << "Initialized job to " << lastjob << ", successes to " << successes << "." << std::endl;
		TR.Debug.flush();
	}

	return;
}

/// @brief Add a checkpoint to the checkpoint file.
/// @details  The checkpoint file must already exist.  Does nothing if checkpointing is disabled.
/// @param[in] curjob The index of the current job just run, for writing to the checkpoint file.
/// @param[in] successes The number of successes so far, for writing to the checkpoint file.
void
SimpleCycpepPredictApplication::checkpoint(
	core::Size const curjob,
	core::Size const successes
) const {
	using namespace utility::io;

	//Do nothing if we're not using checkpointing:
	if ( checkpoint_job_identifier_ == "" ) {
		return;
	}

	ozstream outfile;
	outfile.open_append( checkpoint_filename_ );
	runtime_assert(outfile.good());
	outfile << "LAST\t" << curjob << "\tSUCCESS\t" << successes << std::endl;
	outfile.flush();
	outfile.close();

	store_random_seed_info();

#ifdef BOINC_GRAPHICS
	protocols::boinc::Boinc::update_pct_complete();
#endif

	return;

}

/// @brief End checkpointing and delete the checkpoint file.
/// @details Does nothing if checkpointing is disabled.
void
SimpleCycpepPredictApplication::end_checkpointing() const {
	using namespace utility::io;

	//Do nothing if we're not using checkpointing:
	if ( checkpoint_job_identifier_ == "" ) {
		return;
	}
	runtime_assert( remove( checkpoint_filename_.c_str() ) == 0 );
	erase_random_seed_info();
	return;
}

/// @brief Restore the state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::get_random_seed_info() const {
#ifdef BOINC
	boinc_begin_critical_section();
#endif
	if ( utility::file::file_exists("rng.state.gz") ) {
		utility::io::izstream izs("rng.state.gz");
		numeric::random::rg().restoreState(izs);
		izs.close();
	}
#ifdef BOINC
	boinc_end_critical_section();
#endif
	return;
}

/// @brief Store the state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::store_random_seed_info() const {
#ifdef BOINC
	boinc_begin_critical_section();
#endif
	utility::io::ozstream ozs("rng.state.gz");
	numeric::random::rg().saveState(ozs);
	ozs.close();
#ifdef BOINC
	boinc_end_critical_section();
#endif
	return;
}

/// @brief Erase the stored state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::erase_random_seed_info() const {
	if ( utility::file::file_exists("rng.state.gz") ) {
		utility::file::file_delete("rng.state.gz");
	}
	return;
}

} //cyclic_peptide_predict
} //protocols
