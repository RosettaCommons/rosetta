// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file design_disulfide.cc derived from blivens/convert.cc
/// @brief Converts PDBs to polyCys-> detect disulfide bonds -> combine each disulfide bond into original full atoms PDB -> relax and score -> output scored PDB
/// @author Spencer Bliven <blivens@u.washington.edu>, Kui Chan <kuichan@u.washington.edu>
/// @date Created 2009
/// @details Given a structure, this program will try to find all the possible disulfide bonds.
/// @return a list of scored relaxed PDB with one disulfide each and polycys file.
/// @code design_disulfide -s input.pdb -o output.pdb -database db -detect_disulf



#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.cc.gen.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <core/util/disulfide_util.hh>
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

#include <core/chemical/util.hh>
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/relax/util.hh>


using namespace core;
using namespace std;
using utility::vector1;
using basic::options::option;
using namespace basic::options::OptionKeys;

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>


basic::Tracer TR( "pilot_apps.kuichan.design_disulfide" );

int
usage(char* msg)
{
	TR	<< "usage: design_disulfide -s input.pdb -o output.pdb -database db -detect_disulf" << endl
		<< "Example: design_disulfide.linuxiccrelease -database minirosetta_database -s input.pdb -o output.pdb -ex1 -ex2 -correct -detect_disulf" << endl
		<< msg << endl;
	exit(1);
}


int main( int argc, char * argv [] )
{
    try {
	//init options system
	option.add_relevant( in::file::s );
	option.add_relevant( out::file::o );
	option.add_relevant( in::file::residue_type_set );
	option.add_relevant( out::file::residue_type_set );
	option.add_relevant( in::file::fullatom );
	option.add_relevant( out::file::fullatom );
	option.add_relevant( in::file::centroid_input );
	devel::init(argc, argv);


	if( ! option[ in::file::s ].user() )
		return usage("No in file given: Use -s or -l to designate pdb files to search for disulfides");
	if( ! option[ out::file::o ].user() )
		return usage("No out file given: Use -o to designate an output file");

	string pdb = basic::options::start_file();
	string out = option[ out::file::o ]();


	string in_rsd_set = chemical::FA_STANDARD;

	chemical::ResidueTypeSetCAP rsd_set =
		chemical::ChemicalManager::get_instance()->residue_type_set(in_rsd_set);
	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, *rsd_set, pdb, false);

	pose::Pose start_pose;
	core::import_pose::pose_from_pdb( start_pose, *rsd_set, pdb, false);

	utility::vector1< core::Size > positions;
	for (Size ii = 1; ii <= pose.total_residue(); ii++){
		positions.push_back (ii);
	}


	string out_rsd_set = chemical::CENTROID;
	TR.Info << "Converting from " << in_rsd_set <<" to " << out_rsd_set << endl;
	core::util::switch_to_residue_type_set( pose, out_rsd_set);

	//Write out temp PDB
	pose.dump_pdb(out);

	//Change the residue type in the ATOM line to CYS
	//Write out a out file with .cys ext
	ifstream indata(out.c_str());
	string cenPDB;
	string cysfile = out + ".cys";
	ofstream outdata(cysfile.c_str());
	if (indata.is_open()){
		while (!indata.eof()){
			getline(indata,cenPDB);
			//TR.Info << cenPDB << endl;
			size_t found = cenPDB.find("ATOM");
			if (found!=string::npos){
				cenPDB.replace(17,3,"CYS");
				outdata << cenPDB << "\n";
			}
			//TR.Info << cenPDB << endl;
		}
		indata.close();
		outdata.close();
	}

	//read back the polyCys file
	in_rsd_set = chemical::CENTROID;
	out_rsd_set = chemical::FA_STANDARD;
	rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set(in_rsd_set);
	core::import_pose::pose_from_pdb( pose, *rsd_set, cysfile, false);

	option[ basic::options::OptionKeys::in::detect_disulf ].value( true );
	option[ basic::options::OptionKeys::run::rebuild_disulf ].value( true );

	//convert back to full atoms and build disufide bonds
	TR.Info << "Converting from " << in_rsd_set << " to " << out_rsd_set << endl;
	core::util::switch_to_residue_type_set( pose, out_rsd_set);

	utility::vector1< std::pair<core::Size,core::Size> >  NumDisulfides;
	core::conformation::disulfide_bonds ( pose.conformation(), NumDisulfides );

	for(vector1<pair<Size, Size> >::const_iterator
		disulf(NumDisulfides.begin()), end_disulf(NumDisulfides.end());
		disulf != end_disulf; ++disulf){
		TR.Info << "disulfide:" << disulf->first << " and " << disulf->second << "." << std::endl;
		TR.Info << "origPDB:" << start_pose.aa( disulf->first ) << " and " << start_pose.aa( disulf->second )<< "." << std::endl;
	}

	//write out full atoms pdb with disulfide bonds.
	pose.dump_pdb(out);

	//read in the orig pdb in full atoms mode
	pose::Pose pdbPose;
	in_rsd_set = chemical::FA_STANDARD;
	rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set(in_rsd_set);
	core::import_pose::pose_from_pdb( pdbPose, *rsd_set, pdb, false);


	option[ basic::options::OptionKeys::docking::symmetry::minimize_sidechains ].value( true );


	//incorporate each disulfide bond into original PDB
	//copy start_pose, relax, score and write out the disulfide PDB
	for(vector1<pair<Size, Size> >::const_iterator
		disulf(NumDisulfides.begin()), end_disulf(NumDisulfides.end());
		disulf != end_disulf; ++disulf){


	  if (core::conformation::is_disulfide_bond( pose.conformation(), disulf->first, disulf->second) ){

		TR.Info << "disulfide:" << disulf->first << " and " << disulf->second << "." << std::endl;
		TR.Info << "origPDB:" << start_pose.aa( disulf->first ) << " and " << start_pose.aa( disulf->second )<< "." << std::endl;

		stringstream disulf_out;
		disulf_out << disulf->first << "_" << disulf->second << "_" << pdb;
		TR.Info << "Create PDB:" << disulf_out.str() << std::endl;
		vector1< pair<Size,Size> > disulfides;
		Size s1 (disulf->first);
		Size s2 (disulf->second);
		disulfides.push_back(std::make_pair(s1,s2));

		pose::PoseOP disulf_pose ( new pose::Pose (pdbPose));
		disulf_pose->replace_residue(disulf->first, pose.residue(disulf->first),true);
		disulf_pose->replace_residue(disulf->second, pose.residue(disulf->second),true);
		disulf_pose->conformation().fix_disulfides( disulfides );

		scoring::ScoreFunctionOP sfxn = scoring::getScoreFunction();
		core::util:: rebuild_disulfide(*disulf_pose,disulfides,NULL,sfxn,NULL,sfxn);
		TR.Info << "relax:" << disulf_out.str() << std::endl;
		protocols::relax::relax_pose(*disulf_pose,sfxn,disulf_out.str());

		disulf_pose->dump_scored_pdb(disulf_out.str(),*sfxn,disulf_out.str());
		}

	}

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
} // end main

