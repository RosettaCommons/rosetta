// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///
/// @brief
/// @author Gordon Lemmon

#include <devel/init.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

//Auto Headers
#include <basic/options/option.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/ligand_docking/GrowLigand.hh>

#include <utility/excn/Exceptions.hh>

/// This wrapper exists so David Kim's BOINC executable can call my real main() method.
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

	option.add_relevant(in::file::s);
	option.add_relevant(in::file::native);
	option.add_relevant(in::path::fragments);
	option.add_relevant(out::nstruct);
	option.add_relevant(out::suffix);
	option.add_relevant(out::path::pdb);
	option.add_relevant(docking::ligand::option_file);
	option.add_relevant(enzdes::cstfile);


	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	// Except in unit testing mode, where it wipes out e.g. -database
	devel::init(argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go(new protocols::ligand_docking::GrowLigand("X"));
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
    }
    return 0;
}

