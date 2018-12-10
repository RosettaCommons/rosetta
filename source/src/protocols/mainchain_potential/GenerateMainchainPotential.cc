// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mainchain_potential/GenerateMainchainPotential.cc
/// @brief A generator for mainchain potentials.  Inputs are a noncanonical residue type with an already-generated sidechain potential;
/// outputs are a potential file suitable for use by the RamaPrePro scoreterm.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/mainchain_potential/GenerateMainchainPotential.hh>
#include <protocols/mainchain_potential/GenerateMainchainPotentialOptions.hh>

// Core headers:
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.hh>
#include <core/chemical/mainchain_potential/util.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/TorsionID.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>

// Protocols headers:
#include <protocols/minimization_packing/MinMover.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>

// C++ headers:
#include <iostream>
#include <fstream>

static basic::Tracer TR( "protocols.mainchain_potential.GenerateMainchainPotential" );


namespace protocols {
namespace mainchain_potential {

/// @brief Default constructor.
/// @details Optionally takes a GenerageMainchainPotentialOptions const-owning pointer.  If null, this object initializes itself
/// to default values.
GenerateMainchainPotential::GenerateMainchainPotential( GenerateMainchainPotentialOptionsCOP options /*= nullptr*/ ):
	utility::pointer::ReferenceCount(),
	options_( options == nullptr ? utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) : options->clone() ),
	last_generated_scoretable_(nullptr),
	last_generated_scoretables_by_scoreterm_()
{}

/// @brief Copy constructor.
GenerateMainchainPotential::GenerateMainchainPotential(
	GenerateMainchainPotential const & src
) :
	utility::pointer::ReferenceCount(src),
	options_( src.options_ != nullptr ? options_->clone() : utility::pointer::make_shared< GenerateMainchainPotentialOptions >(false) ),
	last_generated_scoretable_( src.last_generated_scoretable_ == nullptr ? nullptr : src.last_generated_scoretable_->clone() ),
	last_generated_scoretables_by_scoreterm_()
{
	for ( std::map< core::scoring::ScoreType, core::chemical::mainchain_potential::MainchainScoreTableOP >::const_iterator it( src.last_generated_scoretables_by_scoreterm_.begin() ); it != src.last_generated_scoretables_by_scoreterm_.end(); ++it ) {
		last_generated_scoretables_by_scoreterm_[it->first] = it->second->clone();
	}
}

/// @brief Destructor.
GenerateMainchainPotential::~GenerateMainchainPotential() {
	// "CHOOSE.  CHOOSE THE FORM OF THE DESTRUCTOR!"
}

/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
GenerateMainchainPotentialOP
GenerateMainchainPotential::clone() const {
	return GenerateMainchainPotentialOP( utility::pointer::make_shared< GenerateMainchainPotential >( *this ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// PUBLIC FUNCTIONS //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Entry point into protocol execution.
void
GenerateMainchainPotential::run() {
	if ( last_generated_scoretable_ != nullptr && TR.Warning.visible() ) {
		TR.Warning << TR.Red << "Warning!  Overwriting previously-generated scoring table!" << TR.Reset << std::endl;
	}

	// Create the one-residue pose, with suitable terminal patches:
	core::pose::PoseOP pose( generate_pose() );

	// Create the scoring function that we'll use:
	core::scoring::ScoreFunctionOP sfxn( generate_sfxn() );

	// A place to store the data that we generate:
	core::chemical::mainchain_potential::MainchainScoreTableOP newtable( utility::pointer::make_shared< core::chemical::mainchain_potential::MainchainScoreTable >() );

	// Actually generate the potential:
	generate_mainchain_potential( pose, sfxn, newtable, last_generated_scoretables_by_scoreterm_ );

	// Store the result:
	last_generated_scoretable_ = newtable;
}

/// @brief Write the last generated mainchain potential to disk.  The output filename is given in the
/// options_ object.
/// @details If we're writing individual scoretables for individual scoreterms, that happens here, too.
void
GenerateMainchainPotential::write_last_generated_to_disk() const {
	static const std::string errmsg( "Error in GenerateMainchainPotential::write_last_generated_to_disk(): " );
	runtime_assert_string_msg( last_generated_scoretable_ != nullptr, errmsg + "This function cannot be called before a scoretable is generated!" );
	std::stringstream outstream;
	outstream << "# Mainchain potential for residue type " << options_->residue_name() << ".\n";
	outstream << "# File created with the Rosetta make_mainchain_potential application written by Vikram K. Mulligan (vmulligan@flatironinstitute.org).\n";
	last_generated_scoretable_->write_mainchain_scoretable_to_stream(outstream);

	TR << "Generated the following mainchain potential:\n" << outstream.str() << std::endl; //TODO -- switch to debug output.
	TR.flush();

	//Do the actual file output:
	std::string const & output_filename( options_->output_filename() );
	{
		TR << "Writing generated mainchain potential to " << output_filename << "." << std::endl;
		std::ofstream outfile;
		outfile.open( output_filename );
		runtime_assert_string_msg( outfile.good(), errmsg + "Could not write to " + output_filename + "." );
		outfile << outstream.str();
		outfile.close();
		TR << "Successfully wrote generated mainchain potential to " << output_filename << "." << std::endl;
	}

	//If we're writing individual scoretables for individual scoreterms, do that here too.
	if ( options_->write_potentials_for_individual_scoreterms() ) {
		std::string const outfile_base( utility::file::file_basename( output_filename ) );
		std::string const outfile_extension( utility::file::file_extension( output_filename ) );
		for ( std::map< core::scoring::ScoreType, core::chemical::mainchain_potential::MainchainScoreTableOP >::const_iterator it( last_generated_scoretables_by_scoreterm_.begin() ); it!=last_generated_scoretables_by_scoreterm_.end(); ++it ) {
			core::scoring::ScoreType const curtype( it->first );
			std::string const curname( core::scoring::name_from_score_type( curtype ) );
			std::string const curfilename( outfile_base + "_" + curname + "." + outfile_extension );
			TR << "Writing mainchain potential for scoreterm \"" << curname << "\" to file " + curfilename + "." << std::endl;
			std::stringstream outstream2;
			debug_assert( last_generated_scoretables_by_scoreterm_.count(curtype) != 0 );
			last_generated_scoretables_by_scoreterm_.at(curtype)->write_mainchain_scoretable_to_stream(outstream2);

			std::ofstream outfile;
			outfile.open( curfilename );
			runtime_assert_string_msg( outfile.good(), errmsg + "Could not write to " + curfilename + "." );
			outfile << outstream2.str();
			outfile.close();
			TR << "Successfully wrote wrote potential for \"" + curname + "\" scoreterm to " << curfilename << "." << std::endl;
		}
	}

	TR.flush();
}

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// PRIVATE FUNCTIONS /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Get the scorefunction that we'll be using.
core::scoring::ScoreFunctionOP
GenerateMainchainPotential::generate_sfxn() const {
	if ( !(options_->get_scorefxn_filename().empty()) ) {
		TR << "The make_mainchain_potential application is using the user_specified scorefunction defined in \"" << options_->get_scorefxn_filename() << "\"." << std::endl;
		return core::scoring::ScoreFunctionFactory::create_score_function( options_->get_scorefxn_filename() );
	} else { //If the filename is empty, then we use the default.
		TR << "The make_mainchain_potential application is using the default scorefunction defined in \"mainchain_potential_generation.wts\"." << std::endl;
		return core::scoring::ScoreFunctionFactory::create_score_function("mainchain_potential_generation.wts");
	}
}

/// @brief Generate the one-residue pose, based on options.
core::pose::PoseOP
GenerateMainchainPotential::generate_pose() const {
	static const std::string errmsg("Error in protocols::mainchain_potential::GenerateMainchainPotential::generate_pose(): ");

	std::string const & resname( options_->residue_name() );
	runtime_assert_string_msg( !resname.empty(), errmsg + "The residue name cannot be empty!" );

	core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
	core::pose::make_pose_from_sequence( *pose,"X[" + resname + "]", core::chemical::FA_STANDARD );
	debug_assert( pose->total_residue() == 1 ); //Should be a 1-res pose at this point.
	core::pose::remove_lower_terminus_type_from_pose_residue( *pose, 1 );
	core::pose::remove_upper_terminus_type_from_pose_residue( *pose, 1 );

	core::chemical::ResidueTypeCOP restype( pose->residue_type_ptr(1) );
	if ( protein_patches_can_apply( restype ) ) {
		apply_protein_patches(*pose);
	} else {
		utility_exit_with_message( errmsg + "The ResidueType " + restype->base_name() + " is incompatible with any of the known patches for potential generation.  Failing." );
	}

	return pose;
}

/// @brief Can this residue type take protein terminal patches?
bool
GenerateMainchainPotential::protein_patches_can_apply(
	core::chemical::ResidueTypeCOP restype
) const {
	utility::vector1< core::chemical::VariantType > vartypes(2);
	vartypes[1] = core::chemical::ACETYLATED_NTERMINUS_VARIANT;
	vartypes[2] = core::chemical::METHYLATED_CTERMINUS_VARIANT;
	return patches_can_apply( restype, vartypes );
}

/// @brief Can this residue type take specified terminal patches?
bool
GenerateMainchainPotential::patches_can_apply(
	core::chemical::ResidueTypeCOP restype,
	utility::vector1< core::chemical::VariantType > const &vartypes
) const {
	return core::chemical::ResidueTypeFinder( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) ).base_type( restype ).variants( vartypes ).get_representative_type() != nullptr;
}

/// @brief Given a one-residue pose, apply acetylated N-terminus and aminomethylated C-terminus patches.
/// @details Applies aminomethylated or aminodimethylated C-terminus patches depending on settings in options.
void
GenerateMainchainPotential::apply_protein_patches(
	core::pose::Pose & pose
) const {
	if ( options_->make_pre_proline_potential() ) {
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::DIMETHYLATED_CTERMINUS_VARIANT, 1 );
		runtime_assert( pose.residue_type(1).has_variant_type( core::chemical::DIMETHYLATED_CTERMINUS_VARIANT ) );
	} else {
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::METHYLATED_CTERMINUS_VARIANT, 1 );
		runtime_assert( pose.residue_type(1).has_variant_type( core::chemical::METHYLATED_CTERMINUS_VARIANT ) );
	}
	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::ACETYLATED_NTERMINUS_VARIANT, 1 );
	runtime_assert( pose.residue_type(1).has_variant_type( core::chemical::ACETYLATED_NTERMINUS_VARIANT ) );

	pose.update_residue_neighbors();
}

/// @brief Given a pose, cycle through mainchain dihedrals and generate the mainchain potential.
void
GenerateMainchainPotential::generate_mainchain_potential(
	core::pose::PoseCOP masterpose,
	core::scoring::ScoreFunctionOP sfxn,
	core::chemical::mainchain_potential::MainchainScoreTableOP newtable,
	std::map< core::scoring::ScoreType, core::chemical::mainchain_potential::MainchainScoreTableOP >  & last_generated_scoretables_by_scoreterm
) const {
	debug_assert(masterpose->total_residue() == 1 ); //Should be true when this function is called.

	core::chemical::ResidueTypeCOP restype( masterpose->residue_type_ptr(1) );
	utility::vector1< core::Size > const &dimensions( options_->dimensions() );
	utility::vector1< core::Size > const &mainchain_torsions_covered( options_->mainchain_torsions_covered() );
	core::Real const kbt( options_->kbt() );
	bool const symmetrize_output( options_->symmetrize_output() );
	newtable->initialize_for_de_novo_generation( *restype, dimensions, mainchain_torsions_covered, symmetrize_output );

	// Initialize storage for individual scoreterms if we're writing those:
	last_generated_scoretables_by_scoreterm.clear();
	std::map< core::scoring::ScoreType, core::Real > curenergy_by_scoreterm;
	utility::vector1< core::scoring::ScoreType > scoreterms_in_use;
	bool const write_potentials_for_individual_scoreterms( options_->write_potentials_for_individual_scoreterms() );
	if ( write_potentials_for_individual_scoreterms ) {
		TR << "Allocating storage for individual scoreterms." << std::endl;
		for ( core::Size i(1); i<core::scoring::n_score_types; ++i ) {
			if ( sfxn->has_zero_weight( static_cast< core::scoring::ScoreType >(i) ) ) continue; //Do nothing for zero-weighted terms.
			last_generated_scoretables_by_scoreterm[static_cast< core::scoring::ScoreType >(i)] = utility::pointer::make_shared< core::chemical::mainchain_potential::MainchainScoreTable >();
			last_generated_scoretables_by_scoreterm[static_cast< core::scoring::ScoreType >(i)]->initialize_for_de_novo_generation(*restype, dimensions, mainchain_torsions_covered, symmetrize_output);
			curenergy_by_scoreterm[static_cast< core::scoring::ScoreType >(i)] = 0.0;
			scoreterms_in_use.push_back(static_cast< core::scoring::ScoreType >(i));
		}
	}

	utility::vector1< core::Size > curindex( dimensions.size(), 0 );

	core::pack::task::PackerTaskOP task( utility::pointer::make_shared< core::pack::task::PackerTask_ >( *masterpose ) );
	task->initialize_extra_rotamer_flags_from_command_line();
	task->restrict_to_repacking();
	task->set_bump_check(false);

	do {
		//Copy this pose to manipulate it.
		core::pose::PoseOP pose( masterpose->clone() );

		utility::vector1< core::Real > curtorsions;
		newtable->get_mainchain_torsions_from_coords(curindex, curtorsions);
		debug_assert( curtorsions.size() == dimensions.size() ); //Should be true
		//Set mainchain torsions:
		TR << "Setting mainchain torsions to ";
		for ( core::Size i(1), imax(curtorsions.size()); i<=imax; ++i ) {
			pose->set_torsion( core::id::TorsionID( 1, core::id::BB, ( mainchain_torsions_covered.size() == 0 ? i : mainchain_torsions_covered[i] ) ), curtorsions[i] );
			if ( TR.visible() ) {
				TR << curtorsions[i];
				if ( i<imax ) TR << ", ";
			}
		}
		TR << "." << std::endl;
		pose->update_residue_neighbors();

		//Generate rotamers:
		TR << "Generating rotamers." << std::endl;
		utility::graph::GraphOP neighbour_graph( core::pack::create_packer_graph( *pose, *sfxn, task, pose->total_residue(), utility::vector1< core::Real >({50.0}) ) );
		core::pack::rotamer_set::RotamerSetsOP rotsets( utility::pointer::make_shared< core::pack::rotamer_set::RotamerSets >() );
		rotsets->set_task(task);
		rotsets->initialize_pose_for_rotsets_creation( *pose );
		rotsets->build_rotamers( *pose, *sfxn, neighbour_graph );
		debug_assert( rotsets->nmoltenres() == 1 );
		debug_assert( rotsets->nrotamers_for_moltenres(1) > 0);
		core::Size const nrot( rotsets->rotamer_set_for_residue(1)->num_rotamers() );
		TR << "Rotamer generation complete.  " << nrot << " rotamers built." << std::endl;

		//Compute Boltzmann-weighted average energy over all rotamers:
		core::Real partition_fxn(0.0), avg_energy(0.0);
		if ( write_potentials_for_individual_scoreterms ) {
			for ( core::Size i(1), imax(scoreterms_in_use.size()); i<=imax; ++i ) {
				curenergy_by_scoreterm[ scoreterms_in_use[i] ] = 0;
			}
		}
		bool const do_minimization( options_->do_minimization() );
		core::kinematics::MoveMapOP movemap( do_minimization ? utility::pointer::make_shared<core::kinematics::MoveMap>() : nullptr );
		if ( do_minimization ) {
			movemap->set_jump(false);
			movemap->set_bb(false);
			movemap->set_chi(true);
		}
		protocols::minimization_packing::MinMoverOP minmover( do_minimization ? utility::pointer::make_shared< protocols::minimization_packing::MinMover >( movemap, sfxn, options_->minimization_type(), options_->minimization_threshold(), false, false, false ) : nullptr );

		for ( core::Size i(1); i<=nrot; ++i ) { //Loop through all rotamers
			//I'll make a temporary pose, to avoid propagating state accidentally.  There's some overhead to this, but it's better than the mistakes that might arise from not doing this:
			core::pose::PoseOP temppose( pose->clone() );

			core::conformation::Residue const & rotamer( *(rotsets->rotamer_set_for_residue(1)->rotamer(i)) );

			//Set pose to rotamer conformation
			for ( core::Size ichi(1), ichimax( rotamer.nchi() ); ichi<=ichimax; ++ichi ) {
				temppose->set_chi( ichi, 1, rotamer.chi(ichi) );
			}
			temppose->update_residue_neighbors();

			//Do some energy-minimization
			if ( do_minimization ) {
				TR << "Minimizing rotamer " << i << "." << std::endl;
				minmover->apply(*temppose);
			}

			//Compute the energy:
			core::Real const curenergy( (*sfxn)(*temppose) );
			TR << "Energy for rotamer " << i << " is " << curenergy << "." << std::endl;

			//Add to the average energy and partition function computation:
			core::Real const exp_curenergy( std::exp( -curenergy/kbt ) );
			partition_fxn += exp_curenergy;
			avg_energy += exp_curenergy * curenergy;

			//Add to individual scoreterm energies
			if ( write_potentials_for_individual_scoreterms ) {
				for ( core::Size ii(1), iimax(scoreterms_in_use.size()); ii<=iimax; ++ii ) {
					curenergy_by_scoreterm[ scoreterms_in_use[ii] ] += temppose->energies().residue_total_energies(1)[ scoreterms_in_use[ii] ];
				}
			}
		}

		//Normalize with partition function:
		avg_energy /= partition_fxn;
		if ( write_potentials_for_individual_scoreterms ) {
			for ( core::Size i(1), imax(scoreterms_in_use.size()); i<=imax; ++i ) {
				curenergy_by_scoreterm[ scoreterms_in_use[i] ] /= partition_fxn;
				last_generated_scoretables_by_scoreterm[ scoreterms_in_use[i] ]->set_energy( curindex, curenergy_by_scoreterm[ scoreterms_in_use[i] ] );
			}
		}

		if ( TR.visible() ) {
			TR << "Boltzmann-weighted average energy over all rotamers for mainchain torsions ";
			for ( core::Size j(1), jmax(curtorsions.size()); j<=jmax; ++j ) {
				TR << curtorsions[j];
				if ( j<jmax ) TR << ", ";
			}
			TR << " is " << avg_energy << "." << std::endl;
		}

		newtable->set_energy( curindex, avg_energy );

	} while( newtable->increment_energy_coords( curindex ) );

	newtable->finalize_de_novo_scoretable_from_energies( kbt, true );

	if ( write_potentials_for_individual_scoreterms ) {
		for ( core::Size i(1), imax(scoreterms_in_use.size()); i<=imax; ++i ) {
			last_generated_scoretables_by_scoreterm[ scoreterms_in_use[i] ]->finalize_de_novo_scoretable_from_energies(kbt, false /*Do not normalize!*/);
		}
	}
}


} //protocols
} //mainchain_potential
