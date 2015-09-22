// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file hotspot_hash.cc
/// @brief app to find potential hotspot residue placements on the surface of a protein target
///            Method: make a disembodied residue, with option to make residue backbone invisible to the scorefxn
///                    centroid dock residue to the target, then full-atom dock and repack the residue only (must use -norepack1)
///                    write out all docked residues that contact the target with even a somewhat reasonable energy
/// @author Jacob Corn
/// @created May 2008
/// @usage hotspot_hash -randomize1 -randomize2 -norepack1 [-score_interface] -prefix <Rosetta Residue name> -run::n_cycles <# lowres docks> -ntrials <# hires docks> [-sc_only] -user_tag <output suffix>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>

#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>


using basic::T;
using basic::Error;
using basic::Warning;

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace protocols::hotspot_hashing;

static THREAD_LOCAL basic::Tracer TR( "pilot_apps.jcorn.hh_test" );
FileOptionKey const hashfile( "hashfile" );
StringOptionKey const residue( "residue");

int
main( int argc, char * argv [] )
{
	try {

	using namespace scoring;
	option.add( hashfile, "PDB containing existing hash of target");
	option.add( residue, "The 3-letter Rosetta name for the residue to be docked." );

	devel::init(argc, argv);

	protocols::hotspot_hashing::HotspotStubSetOP stubset( new protocols::hotspot_hashing::HotspotStubSet );
	protocols::hotspot_hashing::HotspotStubSetOP new_set( new protocols::hotspot_hashing::HotspotStubSet );
	std::string hashin_fname = "";
	if (option[hashfile].user() )
	{
		hashin_fname = option[hashfile]();
		stubset->read( hashin_fname );
	}

	new_set = stubset->cluster();
	new_set->write_all( "clustered.stubs" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
 }
