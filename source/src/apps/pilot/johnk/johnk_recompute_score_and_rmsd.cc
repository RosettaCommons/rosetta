// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/rms_util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, ref_decoy )

static basic::Tracer TR( "apps.pilot.johnk_recompute_score_and_rmsd.main" );


/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	NEW_OPT( ref_decoy, "the structure to compute RMSD and relative score to", "" );

	devel::init(argc, argv);

	TR << "Starting recomputing scores and rmsds" << std::endl;

	std::string const ref_decoy_fname = option[ ref_decoy ];

	// scoring function
	scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
	//	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );

	// create pose from pdb
	pose::Pose ref_pose;
	core::import_pose::pose_from_pdb( ref_pose, ref_decoy_fname );
	(*scorefxn)(ref_pose);
	core::Real ref_score = ref_pose.energies().total_energies()[ total_score ];

	// Open output file, generate the header line (save it for printing in the log later), print to file
	std::string outfname = "score_vs_rmsd.out";
  utility::io::ozstream outstream;
	outstream.open(outfname, std::ios::out);

	outstream << "fname ca_rms allatom_rms relative_score" << std::endl;

	for (core::Size f=1; f <= basic::options::start_files().size(); f++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(f);
		TR << "Processing decoy " << curr_decoy_fname << std::endl;

		pose::Pose curr_pose;
		core::import_pose::pose_from_pdb( curr_pose, curr_decoy_fname );
		(*scorefxn)(curr_pose);

		core::Real score_diff = curr_pose.energies().total_energies()[ total_score ] - ref_score;
		core::Real CA_rms = rmsd_with_super( ref_pose, curr_pose, is_protein_CA );
		core::Real heavyatom_rms = rmsd_with_super( ref_pose, curr_pose, is_polymer_heavyatom );

		outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << score_diff << std::endl;

	}

	TR << "Done recomputing scores and rmsds" << std::endl;

	outstream.close();
	outstream.clear();

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}



