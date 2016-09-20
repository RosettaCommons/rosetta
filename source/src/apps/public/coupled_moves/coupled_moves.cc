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

#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/coupled_moves.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.coupled_moves" );

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

		// initialize Rosetta
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	TR << "********************d**o**n**e***********************" << std::endl;
	return 0;
}

void *
my_main( void* )
{

	protocols::coupled_moves::CoupledMovesProtocolOP coupled_moves( new protocols::coupled_moves::CoupledMovesProtocol );
	protocols::jd2::JobDistributor::get_instance()->go( coupled_moves );

	return 0;
}
