// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/io/mpistream.hh>
#include <core/kinematics/MoveMap.hh>

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::id;
using namespace protocols::rigid;
using namespace protocols::simple_moves;


OPT_KEY( String, comparison_relax_pdb_num )
OPT_KEY( String, template_pdb_name )
OPT_KEY( String, source_contact_list )
OPT_KEY( String, template_contact_list )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.david_pocket_align_and_save.main" );

//set to store pdb info keys
std::list <std::string> source_interface;
std::list <std::string> template_interface;

bool
is_interface_heavyatom(
	core::pose::Pose const & pose,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	// ws get residue "key" for list
	std::ostringstream residuestream;
	residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
	std::string res_id = residuestream.str();

	core::conformation::Residue const & rsd = pose.residue(resno);
	bool found = false;
	if ( rsd.is_protein() && !rsd.atom_is_hydrogen(atomno) ) {
		for ( std::list<std::string>::iterator it=source_interface.begin(); it!=source_interface.end(); ++it ) {
			if ( !res_id.compare(*it) ) {
				std::cout<<*it<<" "<<res_id<<std::endl;
				found=true;
				break;
			}
		}
	}
	return found;
}

bool
is_interface_heavyatom_pair(
	core::pose::Pose const & pose,
	core::pose::Pose const & pose2,
	core::Size resno,
	core::Size resno2,
	core::Size atomno
)
{
	// ws get residue "key" for list
	std::ostringstream residuestream;
	residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
	std::string res_id = residuestream.str();
	std::ostringstream residuestream2;
	residuestream2 << pose2.pdb_info()->chain(resno2) << pose2.pdb_info()->number(resno2);
	std::string res_id2 = residuestream2.str();

	core::conformation::Residue const & rsd = pose.residue(resno);
	core::conformation::Residue const & rsd2 = pose2.residue(resno2);
	bool found = false;
	if ( rsd.is_protein() && !rsd.atom_is_hydrogen(atomno) ) {
		if ( rsd2.is_protein() && !rsd2.atom_is_hydrogen(atomno) ) {
			for ( std::list<std::string>::iterator it=source_interface.begin(), it2=template_interface.begin(); it!=source_interface.end() && it2!=source_interface.end(); ++it, ++it2 ) {
				if ( !res_id.compare(*it) ) {
					if ( !res_id2.compare(*it2) ) {
						found = true;
						break;
					}
				}
			}
		}
	}

	return found;
}
Real
iface_pdb_superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
)
{
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::BOGUS_ATOM_ID );
	for ( Size ii = 1; ii <= mod_pose.total_residue(); ++ii ) {
		if ( ! mod_pose.residue(ii).has("CA") ) continue;
		if ( ! mod_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= ref_pose.total_residue(); ++jj ) {
			if ( ! ref_pose.residue(jj).has("CA") ) continue;
			if ( ! ref_pose.residue(jj).is_protein() ) continue;
			if ( mod_pose.pdb_info()->chain(ii) != ref_pose.pdb_info()->chain(jj) ) continue;
			if ( mod_pose.pdb_info()->number(ii) != ref_pose.pdb_info()->number(jj) ) continue;
			if ( is_interface_heavyatom ( ref_pose, mod_pose, jj, 2) ) {
				id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
				id::AtomID const id2( ref_pose.residue(jj).atom_index("CA"), jj );
				atom_map.set( id1, id2 );
			}

			break;
		}

	}
	return superimpose_pose( mod_pose, ref_pose, atom_map );
}

Real
iface_pdb_superimpose_diff_prot(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
)
{
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::BOGUS_ATOM_ID );
	for ( Size ii = 1; ii <= mod_pose.total_residue(); ++ii ) {
		if ( ! mod_pose.residue(ii).has("CA") ) continue;
		if ( ! mod_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= ref_pose.total_residue(); ++jj ) {
			if ( ! ref_pose.residue(jj).has("CA") ) continue;
			if ( ! ref_pose.residue(jj).is_protein() ) continue;
			if ( is_interface_heavyatom_pair ( mod_pose, ref_pose, ii, jj, 2) ) {
				id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
				id::AtomID const id2( ref_pose.residue(jj).atom_index("CA"), jj );
				atom_map.set( id1, id2 );
				break;
			}

		}

	}
	return superimpose_pose( mod_pose, ref_pose, atom_map );
}


/// General testing code
int main( int argc, char * argv [] ) {

	try{

		NEW_OPT( comparison_relax_pdb_num, "comparison residue", "-1");
		NEW_OPT( template_pdb_name, "template pdb", "" );
		NEW_OPT ( source_contact_list, "File name for optional list of contact residues to align the source pdb on","");
		NEW_OPT ( template_contact_list, "File name for optional list of contact residues to align to the template pdb","");

		TR << "Calling init" << std::endl;
		//initializes Rosetta functions
		devel::init(argc, argv);

		TR << "done" << std::endl;
		//allows output when running the program
		//sets input residue numbers to resid and resid_c
		std::string const resid_c = option[comparison_relax_pdb_num];
		TR << "Starting pocket compare" << std::endl;

		// create pose for comparison pose from pdb
		std::string const comparison_pdb_name ( basic::options::start_file() );
		pose::Pose comparison_pose;
		core::import_pose::pose_from_file( comparison_pose, comparison_pdb_name , core::import_pose::PDB_file);
		TR << "set comparison pdb"<< "    Number of residues: " << comparison_pose.total_residue() << std::endl;

		std::string tag = "";
		if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
			tag = "." + option[ OptionKeys::out::output_tag ]();
		}

		std::string const template_fname ( option[ template_pdb_name ] );
		if ( template_fname != "" ) {
			std::cout << "template\n";
			//sets template input pdb name

			TR << "set template pdb" << template_fname << std::endl;
			//sets pdb as a Rosetta pose
			pose::Pose template_pose;
			core::import_pose::pose_from_file( template_pose, template_fname , core::import_pose::PDB_file);
			TR << "set template pdb" << "    Number of residues: " << template_pose.total_residue() << std::endl;


			std::string const cfilename = option[ source_contact_list ];
			if ( cfilename != "" ) {
				std::cout<<"source contact list\n";
				std::ifstream ifs(cfilename.c_str(), std::ifstream::in);
				if ( !ifs.is_open() ) {
					std::cout<< "Error opening source contact list file "<<cfilename<<std::endl;
					return -100;
				}
				//ifb.open (cfilename,std::ios::in);
				//std::ostream ios(&ifb);
				std::string intres;
				while ( ifs.good() ) {
					ifs >> intres;
					source_interface.push_back(intres);
				}

				std::string const tcfilename = option[ template_contact_list ];
				if ( tcfilename != "" ) {
					std::cout<<"template contact list\n";
					std::ifstream tifs(tcfilename.c_str(), std::ifstream::in);
					if ( !tifs.is_open() ) {
						std::cout<< "Error opening template contact list file "<<tcfilename<<std::endl;
						return -101;
					}
					//ifb.open (cfilename,std::ios::in);
					//std::ostream ios(&ifb);
					std::string tintres;
					while ( tifs.good() ) {
						tifs >> tintres;
						template_interface.push_back(tintres);
					}
					std::cout << "aligning\n";
					iface_pdb_superimpose_diff_prot( comparison_pose, template_pose);
				} else {
					std::cout<<"no template contact list\n";
					iface_pdb_superimpose_pose( comparison_pose, template_pose);
				}
			} else {
				std::cout<<"no source contact list\n";

				// align comparison pose to template pose
				// protocols::simple_moves::SuperimposeMoverOP sp_mover = new protocols::simple_moves::SuperimposeMover( template_pose );
				protocols::simple_moves::SuperimposeMoverOP sp_mover( new protocols::simple_moves::SuperimposeMover() );
				sp_mover->set_reference_pose( template_pose, 1, template_pose.total_residue() );
				sp_mover->set_target_range( 1, template_pose.total_residue() );
				sp_mover->apply( comparison_pose );
			}
		}
		std::cout << "done aligning\n";
		std::vector< conformation::ResidueCOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(comparison_pose, resid_c);
		// call function to make a grid around a target residue (seqpos)
		protocols::pockets::PocketGrid comparison_pg( residues );
		//call function to define the pocket
		comparison_pg.zeroAngle();
		comparison_pg.autoexpanding_pocket_eval( residues, comparison_pose ) ;
		//output the pocket in a pdb
		//comparison_pg.dumpGridToFile();
		TR << "Comparison pocket defined" << std::endl;
		std::cout << "Pocket score (unweighted) is: " << comparison_pg.netTargetPocketVolume() << std::endl;
		// dump to file, with name based on input pdb
		std::stringstream out_fname;
		out_fname << comparison_pdb_name << tag << ".pocket.pdb";
		std::cout<<out_fname.str() <<std::endl;
		comparison_pg.dumpTargetPocketsToPDB( out_fname.str() );
		if ( option[ OptionKeys::pocket_grid::pocket_dump_exemplars ]() ) {
			std::stringstream out_exfname;
			out_exfname << comparison_pdb_name << tag << ".exemplar.pdb";
			std::cout<<out_exfname.str() <<std::endl;
			comparison_pg.dumpExemplarToFile( out_exfname.str() );
		}
		std::stringstream out_pfname;
		out_pfname << comparison_pdb_name << tag << ".pocket";
		std::cout<<out_pfname.str() <<std::endl;
		comparison_pg.dumpTargetPocketsToFile( out_pfname.str() );

		//utility::io::ozstream fout;
		//fout.open("distance.txt", std::ios::out);
		//fout << d1 << std::endl;
		//fout.close();
		//  fout.clear();

		TR << "Done!" << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

