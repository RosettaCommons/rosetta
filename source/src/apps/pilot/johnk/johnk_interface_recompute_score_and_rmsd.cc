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
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
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

static basic::Tracer TR( "apps.pilot.johnk_interface_recompute_score_and_rmsd.main" );

// NOTE: GLOBAL VARIABLE - this would better be done with a static (and could be used to ensure it's filled, too...)
utility::vector1 <bool> interface;

void define_interface( core::pose::Pose const & pose ) {
	interface.resize( pose.total_residue(), false );
	core::Real interface_dist = 8.0;
	core::Size rb_jump = 1;
	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
	pack::task::TaskFactory tf;
	tf.push_back( new protocols::toolbox::task_operations::RestrictToInterface( rb_jump, interface_dist ) );
	pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );
	for ( core::Size i=1; i <= pose.total_residue(); ++i ) {
		if ( task->pack_residue(i) ) interface.at(i)=true;
	}
}

bool
is_interface_CA(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	if ( ! interface.at(resno) ) return false;
	return rsd.is_protein() && rsd.has("CA") && rsd.atom_index("CA") == atomno;
}

bool
is_interface_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	if ( ! interface.at(resno) ) return false;
	return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);
}


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
	scoring::ScoreFunctionOP scorefxn( getScoreFunction() );

	// create pose from pdb
	pose::Pose ref_pose;
	core::import_pose::pose_from_pdb( ref_pose, ref_decoy_fname );
	(*scorefxn)(ref_pose);
	core::Real ref_score = ref_pose.energies().total_energies()[ total_score ];

	define_interface( ref_pose );

	FArray1D_bool superpos_partner ( ref_pose.total_residue(), false );
	for ( Size i=1; i<= ref_pose.total_residue(); ++i ) {
		if ( ref_pose.pdb_info()->chain(i) == 'R' ) superpos_partner(i)=true;
	}

	// Open output file, generate the header line (save it for printing in the log later), print to file
	std::string outfname = "lig_score_vs_rmsd.out";
  utility::io::ozstream outstream;
	outstream.open(outfname, std::ios::out);

	outstream << "fname full_CA_rms interface_CA_rms interface_heavyatom_rms relative_score" << std::endl;

	for (core::Size f=1; f <= basic::options::start_files().size(); f++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(f);
		TR << "Processing decoy " << curr_decoy_fname << std::endl;

		pose::Pose curr_pose;
		core::import_pose::pose_from_pdb( curr_pose, curr_decoy_fname );
		(*scorefxn)(curr_pose);

		//		core::Real score_diff = curr_pose.energies().total_energies()[ total_score ] - ref_score;
		//		core::Real full_CA_rms = rmsd_with_super( ref_pose, curr_pose, is_protein_CA );
		//		core::Real interface_CA_rms = rmsd_with_super( ref_pose, curr_pose, is_interface_CA );
		//		core::Real interface_heavyatom_rms = rmsd_with_super( ref_pose, curr_pose, is_interface_heavyatom );

		core::Real score_diff = curr_pose.energies().total_energies()[ total_score ] - ref_score;
		core::Real full_CA_rms = rmsd_no_super_subset( ref_pose, curr_pose, superpos_partner, is_protein_CA );
		core::Real interface_CA_rms = rmsd_no_super_subset( ref_pose, curr_pose, superpos_partner, is_interface_CA );
		core::Real interface_heavyatom_rms = rmsd_no_super_subset( ref_pose, curr_pose, superpos_partner, is_interface_heavyatom );

		outstream << curr_decoy_fname << ' ' << full_CA_rms << ' ' << interface_CA_rms << ' ' << interface_heavyatom_rms << ' ' << score_diff << std::endl;

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



