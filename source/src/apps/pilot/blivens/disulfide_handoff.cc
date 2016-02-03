// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file disulfide_handoff.cc
/// @brief Removes disulfides from a structure and tries to replace them.
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created September 2008
/// @details
/// This app was written to experiment with techniques for preserving
/// disulfide bonds during the handoff between centroid and full atom.
/// Techniques tried include repacking, minimization, constraints.
/// @section cli Command Line
/// @code disulfide_staple -s input.pdb -o output.pdb -database db @endcode


//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>

//Core Chemistry
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>


//Command line Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Packing
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>

//Scoring
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <apps/pilot/blivens/disulfides.hh>

#include <utility>

using namespace core;
using namespace std;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::chemical;
using namespace core::conformation;

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

static THREAD_LOCAL basic::Tracer TR( "pilot_apps.blivens.disulfide_handoff" );

int
usage(char* msg)
{
	TR	<< "usage: disulfide_staple -s input.pdb -o output.pdb -database db" << endl
		<< msg << endl;
	exit(1);
}


int main( int argc, char * argv [] )
{
  try {
	string infile;
	string outfile;

	//init options system
	option.add_relevant( in::file::s );
	option.add_relevant( out::file::o );

	devel::init(argc, argv);

	if( option[ in::file::s ].user() ) {
		infile = basic::options::start_file();
	}
	else {
		return usage("No in file given: Use -s to designate a pdb file");
	}

	if( option[ out::file::o ].user() ) {
		outfile = option[ out::file::o ]();
	}
	else {
		return usage("No out file given: Use -o to designate a pdb file");
	}
	//done with options

	pose::PoseOP pose(new pose::Pose);
	core::import_pose::pose_from_file( *pose, infile , core::import_pose::PDB_file);

	scoring::ScoreFunctionOP sfxn = scoring::get_score_function();


	//initialize vectors of all disulf bonds
	vector1< pair<Size,Size> > disulfides;
	for(Size i(1); i<= pose->total_residue()-1;++i) {
		for(Size j(i+1); j<= pose->total_residue();++j) {
			if( actual_disulfide(*pose,i,j) ) {
				disulfides.push_back(make_pair(i,j));
			}
		}
	}
	for(vector1< pair<Size,Size> >::const_iterator ds_it = disulfides.begin(); ds_it != disulfides.end(); ++ds_it) {
		TR <<"Found DS at "<<ds_it->first<<" to "<<ds_it->second<<endl;
	}

	//////////////////////
	//Convert to centroid and back
	//core::util::switch_to_residue_type_set( cen_pose, chemical::CENTROID);
	//core::util::switch_to_residue_type_set( cen_pose, chemical::FA_STANDARD);
	//TR << "Writing fa version after conversion to "<< outfile+"_fa.pdb" << endl;
	//cen_pose.dump_scored_pdb(outfile+"_fa.pdb",*sfxn,"fa");

	// Force disulfides
	for(vector1<pair<Size, Size> >::const_iterator
			disulf(disulfides.begin()), end_disulf(disulfides.end());
			disulf != end_disulf; ++disulf)
	{
	  core::conformation::form_disulfide(*pose.conformation(), disulf->first, disulf->second);
	}

	// Setup Packer & Minimizer
	pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( *pose );
	task->initialize_from_command_line().or_include_current( true );
	task->restrict_to_repacking();

	kinematics::MoveMapOP mm(new kinematics::MoveMap);
	mm->set_bb( false );

	// Set up each residue individually
	for( Size i(1); i <= pose->total_residue(); ++i )
	{
		Residue const& res(pose->residue(i));
		if( !res.is_protein() )
			continue;

		// Determine if i is part of disulfides
		bool is_disulf = false;
		for(vector1<pair<Size, Size> >::const_iterator
				disulf(disulfides.begin()), end_disulf(disulfides.end());
				disulf != end_disulf; ++disulf)
		{
			if( i == disulf->first || i == disulf->second ) {
				is_disulf = true;
				break;
			}
		}

		if( is_disulf ) {
			// repack & minimize disulfides
			mm->set_chi(i, true);
		} else {
			// Other residues are unchanged
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}

	// Rebuild disulfides
	core::pose::rebuild_disulfide(*pose,disulfides, task, NULL, mm, NULL);

	//cen_pose.conformation().detect_disulfides();
	pose->dump_scored_pdb(outfile, *sfxn, "");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // end main

