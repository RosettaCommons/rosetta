// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   app/public/ligand_docking/ligand_dock.cc
///
/// @brief Ligand docking application
/// @author Rocco Moretti (rmoretti@u.washington.edu) Ian Davis (ian.w.davis@gmail.com)

// Unit Headers

#include <devel/init.hh>
#include <protocols/ligand_docking/ligand_dock_impl.hh>

// Utility Headers

#include <basic/Tracer.hh>

#include <basic/options/option.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

static basic::Tracer TR( "ligand_dock.main" );

int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		option.add_relevant(in::path::database);
		option.add_relevant(in::file::extra_res_fa);
		option.add_relevant(packing::unboundrot);
		option.add_relevant(packing::ex1::ex1);
		option.add_relevant(packing::ex1aro::ex1aro);
		option.add_relevant(packing::ex2::ex2);
		option.add_relevant(packing::extrachi_cutoff);
		option.add_relevant(packing::no_optH);
		option.add_relevant(packing::flip_HNQ);
		option.add_relevant(docking::ligand::soft_rep);
		option.add_relevant(docking::ligand::old_estat);

		option.add_relevant(in::file::s);
		option.add_relevant(in::file::native);
		option.add_relevant(out::nstruct);
		option.add_relevant(out::suffix);
		option.add_relevant(out::path::pdb);
		option.add_relevant(docking::uniform_trans);
		option.add_relevant(docking::randomize2);
		option.add_relevant(docking::dock_pert);
		option.add_relevant(docking::ligand::random_conformer);
		option.add_relevant(docking::ligand::improve_orientation);
		option.add_relevant(docking::ligand::minimize_ligand);
		option.add_relevant(docking::ligand::harmonic_torsions);
		option.add_relevant(docking::ligand::minimize_backbone);
		option.add_relevant(docking::ligand::harmonic_Calphas);
		option.add_relevant(docking::ligand::protocol);
		option.add_relevant(docking::ligand::start_from);
		option.add_relevant(docking::ligand::mutate_same_name3);
		option.add_relevant(enzdes::cstfile);

		// Parses command line options and inits RNG.
		// Doesn't seem to hurt to do it again if already done once (?)
		// Except in unit testing mode, where it wipes out e.g. -database
		devel::init(argc, argv);

		// Giving an AtomTreeDiff input file to ligand_dock is uncommon enough that it isn't worth special casing
		if ( option[ in::file::silent ].active() ) {
			TR.Warning << "Current version of ligand_dock expects regular silent file with -in:file:silent" << std::endl;
			TR.Warning << "         Use -in:file:atom_tree_diff for inputting Atom Tree Diff format files" << std::endl;
		}

		// Backward compatability default
		// If we've haven set any other JD2 output options, use silent.out with Atom Tree Diff format
		// (Options tested are a hacky repeat of cases in protocols::jd2::JobDistributorFactory::create_job_outputter() )
		if ( ! option[ out::file::silent ].user() && ! option[ out::pdb ].user() && ! option[ out::file::score_only ].user() &&
				! option[ jd2::no_output ].user() && ! option[ out::nooutput ].user() &&
				! option[ jd2::enzdes_out ].user() && ! option[ out::use_database ].user() ) {
			option[ out::file::atom_tree_diff ].value( "silent.out" );
		}

		// Backward compatibility options munging.
		if ( option[ out::file::silent ].user() && ! option[ out::file::silent_struct_type ].user() ) {
			TR.Warning << "For backward compatibility, by default ligand_dock will output Atom Tree Diff format files to -out:file:silent" << std::endl;
			TR.Warning << "         To output regular format silent file format, explicitly set -out:file:silent_struct_type" << std::endl;
			option[ out::file::atom_tree_diff ].value( option[ out::file::silent ] );
			option[ out::file::silent ].deactivate();
		}

		// Note: Another difference from the jd1 version of ligand_dock is that the job order for -s/-l isn't randomized

		return ligand_dock_main();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
