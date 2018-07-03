// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/coupled_moves/CoupledMovesProtocol.cc
/// @brief executable for the CoupledMovesProtocol
/// @author Noah <nollikai@gmail.com>, refactored by Steven Lewis smlewi@gmail.com

#include <protocols/coupled_moves/CoupledMovesProtocol.hh>

//KIC
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.hh>
#include <protocols/kinematic_closure/perturbers/WalkingPerturber.hh>
#include <protocols/kinematic_closure/KicMover.hh>

#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/coupled_moves.OptionKeys.gen.hh>

static basic::Tracer TR( "apps.coupled_moves" );

void *
my_main( void* );

int
main( int argc, char * argv [] )
{
	try {

		OPT(in::path::database);
		OPT(in::ignore_unrecognized_res);
		OPT(out::nstruct);
		OPT(packing::resfile);
		OPT(in::file::native);
		OPT(constraints::cst_fa_weight);
		OPT(constraints::cst_fa_file);
		OPT(out::pdb_gz);
		OPT(coupled_moves::ntrials);
		OPT(coupled_moves::mc_kt);
		OPT(coupled_moves::boltzmann_kt);
		OPT(coupled_moves::mm_bend_weight);
		OPT(coupled_moves::trajectory);
		OPT(coupled_moves::trajectory_gz);
		OPT(coupled_moves::trajectory_stride);
		OPT(coupled_moves::trajectory_file);
		OPT(coupled_moves::output_fasta);
		OPT(coupled_moves::output_stats);
		OPT(coupled_moves::ligand_mode);
		OPT(coupled_moves::initial_repack);
		OPT(coupled_moves::min_pack);
		OPT(coupled_moves::save_sequences);
		OPT(coupled_moves::save_structures);
		OPT(coupled_moves::ligand_prob);
		OPT(coupled_moves::fix_backbone);
		OPT(coupled_moves::uniform_backrub);
		OPT(coupled_moves::bias_sampling);
		OPT(coupled_moves::bump_check);
		OPT(coupled_moves::ligand_weight);
		OPT(coupled_moves::output_prefix);
		//OPT(coupled_moves::walking_perturber_magnitude);

		// initialize Rosetta
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	TR << "********************d**o**n**e***********************" << std::endl;
	return 0;
}

void *
my_main( void* )
{

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// Make coupled moves
	protocols::coupled_moves::CoupledMovesProtocolOP coupled_moves( new protocols::coupled_moves::CoupledMovesProtocol );

	// Set command line options for coupled moves

	// Set backbone mover
	if ( option[ basic::options::OptionKeys::coupled_moves::backbone_mover ].user() ) {
		std::string backbone_mover = option[ basic::options::OptionKeys::coupled_moves::backbone_mover ];
		if ( backbone_mover == "kic" ) {
			coupled_moves->set_backbone_mover( backbone_mover );
			std::string kic_perturber = option[ basic::options::OptionKeys::coupled_moves::kic_perturber ];
			if ( kic_perturber == "fragment" ) {
				if ( ( option[ basic::options::OptionKeys::loops::frag_files ].user() ) && ( option[ basic::options::OptionKeys::loops::frag_sizes ].user() ) ) {
					utility::vector1<core::fragment::FragSetOP> frag_libs;
					protocols::loops::read_loop_fragments(frag_libs); // this function uses OptionKeys::loops::frag_sizes and OptionKeys::loops::frag_files to fill the frag_libs object, which is then used as an argument for the FragmentPerturber constructor.
					coupled_moves->set_perturber( protocols::kinematic_closure::perturbers::PerturberOP( new protocols::kinematic_closure::perturbers::FragmentPerturber(frag_libs) ) );
				} else {
					std::stringstream message;
					message << "[ ERROR - fragments ] Must specify the -loops:frag_sizes and -loops:frag_files " << std::endl;
					message << "[ ERROR - fragments ] options in order to use FragmentPerturber." << std::endl;
					throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
				}
			} else if ( kic_perturber == "walking" ) {
				coupled_moves->set_perturber( protocols::kinematic_closure::perturbers::PerturberOP( new protocols::kinematic_closure::perturbers::WalkingPerturber( option[ basic::options::OptionKeys::coupled_moves::walking_perturber_magnitude ] ) ) );
			}
			// Set kic loop size
			if ( !option[ basic::options::OptionKeys::coupled_moves::kic_loop_size ].user() ) {
				TR << TR.Red << "[ WARNING - kic_loop_size ] You did not specify -coupled_moves::kic_loop_size option." << std::endl;
				TR << TR.Red << "[ WARNING - kic_loop_size ] Using default, which is probably fine." << std::endl;
				TR << TR.Red << "[ WARNING - kic_loop_size ] See CoupledMoves wiki or doxygen for more information on kic_loop_size." << TR.Reset << std::endl;
			}
			coupled_moves->set_loop_size( option[ basic::options::OptionKeys::coupled_moves::kic_loop_size ] );
		} else if ( ( backbone_mover == "backrub" ) || ( backbone_mover == "" ) ) {
			coupled_moves->set_backbone_mover( backbone_mover );
		} else {
			std::stringstream message;
			message << "[ ERROR - backbone_mover ] Specified -backbone_mover '" << backbone_mover << "' is not recognized by Coupled Moves. Try 'kic' or 'backrub' instead." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
		}
	} else if ( !option[ basic::options::OptionKeys::coupled_moves::backbone_mover ].user() ) {
		// If no backbone mover is specified in command line
		TR << TR.Red << "[ WARNING - backbone_mover ] You did not specify -coupled_moves::backbone_mover option." << std::endl;
		TR << TR.Red << "[ WARNING - backbone_mover ] Defaulting to legacy behavior '-coupled_moves::backbone_mover backrub'." << std::endl;
		TR << TR.Red << "[ WARNING - backbone_mover ] Example usages:" << std::endl;
		TR << TR.Red << "[ WARNING - backbone_mover ]     '-coupled_moves::backbone_mover kic'" << std::endl;
		TR << TR.Red << "[ WARNING - backbone_mover ]     '-coupled_moves::backbone_mover backrub'" << std::endl;
		TR << TR.Red << "[ WARNING - backbone_mover ] See CoupledMoves wiki or doxygen for more information on available backbone movers." << TR.Reset << std::endl;
		coupled_moves->set_backbone_mover( "backrub" );
	}

	// Set repack_neighborhood option
	if ( option[ basic::options::OptionKeys::coupled_moves::repack_neighborhood ].user() ) {
		bool repack_neighborhood = option[ basic::options::OptionKeys::coupled_moves::repack_neighborhood ];
		coupled_moves->set_repack_neighborhood( repack_neighborhood );
	}

	protocols::jd2::JobDistributor::get_instance()->go( coupled_moves );

	return nullptr;
}
