// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk + dj

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <sstream>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/pockets/PocketConstraint.cc>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/scoring/rms_util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, ref_decoy )
OPT_KEY( String, input_ligand_file )


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );

//set to store pdb info keys
std::set <std::string> interface;
//stores resid of the ligand residue
core::Size lig_res_num;

void define_interface( core::pose::Pose & ref_pose ) {
	lig_res_num =0;
	for ( int j = 1, resnum = ref_pose.total_residue(); j <= resnum; ++j ) {
		if ( !ref_pose.residue(j).is_protein() ) {
			lig_res_num = j;
			break;
		}
	}
	if ( lig_res_num == 0 ) {
		TR << "No ligand given in reference PDB structure.  Cannot identify interface."<<std::endl;
		exit (1);
	}

	TR <<"sele ";
	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	(*scorefxn)(ref_pose);
	EnergyGraph & energy_graph(ref_pose.energies().energy_graph());
	for ( graph::Graph::EdgeListIter
			iru  = energy_graph.get_node( lig_res_num )->lower_edge_list_begin(),
			irue = energy_graph.get_node( lig_res_num )->lower_edge_list_end();
			iru != irue; ++iru ) {
		EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
		Size const j( edge->get_first_node_ind() );

		// the pair energies cached in the link
		EnergyMap const & emap( edge->fill_energy_map());
		Real const attr( emap[ fa_atr ] );
		//TR<<"\n"<<j<<": "<<attr<<"\n";
		if ( attr < -.2 ) {
			// create string id to store in set
			std::ostringstream residuestream;

			TR << "resi "<< ref_pose.pdb_info()->number(j)<<" or ";

			residuestream << ref_pose.pdb_info()->chain(j) << ref_pose.pdb_info()->number(j);
			std::string res_id = residuestream.str();
			interface.insert(res_id);


		}
	}

	TR << std::endl;
	TR << lig_res_num<< std::endl;


	// ref_pose.delete_polymer_residue(lig_res_num);
}

bool
is_interface_heavyatom(
	core::pose::Pose const & pose,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	// ws get residue "key" for set
	std::ostringstream residuestream;
	residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
	std::string res_id = residuestream.str();

	core::conformation::Residue const & rsd = pose.residue(resno);
	if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);

	return false;
}

bool
is_interface_bbatom(
	core::pose::Pose const & pose,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	// ws get residue "key" for set
	std::ostringstream residuestream;
	residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
	std::string res_id = residuestream.str();

	core::conformation::Residue const & rsd = pose.residue(resno);
	if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && !rsd.atom_is_backbone(atomno);

	return false;
}

Real
calpha_pdb_superimpose_pose(
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
			id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
			id::AtomID const id2( ref_pose.residue(jj).atom_index("CA"), jj );
			atom_map.set( id1, id2 );
			break;
		}

	}
	return superimpose_pose( mod_pose, ref_pose, atom_map );
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


template< class T >
Real
interface_rmsd(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	T* predicate
)
{
	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( Size ii = 1; ii <= ref_pose.total_residue(); ++ii ) {
		if ( ! ref_pose.residue(ii).has("CA") ) continue;
		if ( ! ref_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= mod_pose.total_residue(); ++jj ) {
			if ( ! ref_pose.residue(ii).has("CA") ) continue;
			if ( ! ref_pose.residue(ii).is_protein() ) continue;
			if ( mod_pose.pdb_info()->chain(jj) != ref_pose.pdb_info()->chain(ii) ) continue;
			if ( mod_pose.pdb_info()->number(jj) != ref_pose.pdb_info()->number(ii) ) continue;
			Size num_atoms ( ref_pose.residue(ii).natoms() );

			for ( core::Size i = 1; i <= num_atoms; ++i ) {
				if ( predicate ( ref_pose, mod_pose, ii, i) ) {
					Size num_atoms2 ( mod_pose.residue(jj).natoms() );
					for ( core::Size j = 1; j <= num_atoms2; ++j ) {
						if ( !ref_pose.residue(ii).atom_name(i).compare(mod_pose.residue(jj).atom_name(j)) ) {
							p1_coords.push_back(ref_pose.residue(ii).xyz(i));
							p2_coords.push_back(mod_pose.residue(jj).xyz(j));
						}
					}
				}
			}

		}
	}
	assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}

	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );

}


/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT( ref_decoy, "the structure to compute RMSD and relative score to", "" );
		NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );

		core::init::init(argc, argv);

		TR << "Starting recomputing scores and rmsds" << std::endl;

		std::string const ref_decoy_fname = option[ ref_decoy ];
		std::string const input_ligand = option[ input_ligand_file ];

		// create pose from pdb
		pose::Pose ref_pose;
		core::import_pose::pose_from_file( ref_pose, ref_decoy_fname , core::import_pose::PDB_file);

		define_interface( ref_pose );

		TR << "Defined interface" << std::endl;

		std::string outfname;
		if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
			outfname = "alignedrmsd." + option[ OptionKeys::out::output_tag ]() + ".out";
		} else {
			outfname = "alignedrmsd.out";
		}
		//std::cout<<outfname<<" output_tag: "<<option[ OptionKeys::out::output_tag ]()<<std::endl;
		utility::io::ozstream outstream;
		outstream.open(outfname, std::ios::out);

		//outstream << "fname allatom_rms iface_rms" << std::endl;

		for ( core::Size f=1; f <= basic::options::start_files().size(); f++ ) {

			std::string const curr_decoy_fname = basic::options::start_files().at(f);
			TR << "Processing decoy " << curr_decoy_fname << std::endl;

			pose::Pose curr_pose;
			core::import_pose::pose_from_file( curr_pose, curr_decoy_fname , core::import_pose::PDB_file);
			core::Real CA_rms;
			CA_rms = iface_pdb_superimpose_pose( curr_pose, ref_pose);
			CA_rms = core::scoring::CA_rmsd( curr_pose, ref_pose );
			std::cout << "after superimpose to native. Rms to native: " << CA_rms << std::endl;


			//core::Real CA_rms = rmsd_with_super( ref_pose, curr_pose, is_protein_CA );
			//std::cout << "CA_rms to native: " << CA_rms << std::endl;
			//core::Real heavyatom_rms = rmsd_with_super( ref_pose, curr_pose, is_interface_heavyatom );
			core::Real heavyatom_rms = interface_rmsd( ref_pose, curr_pose, is_interface_heavyatom);
			core::Real bbatom_rms = interface_rmsd( ref_pose, curr_pose, is_interface_bbatom );
			std::cout << "Interface rmsd, fa: " << heavyatom_rms << " bb: "<< bbatom_rms<<std::endl;

			//outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << score_diff << ' ' << p_score_diff <<std::endl;
			outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' '<< bbatom_rms<<std::endl;

			//create 'tag' for output filenames
			//int dot_index1 = input_protein.rfind(".", input_protein.size());
			//int slash_index1 = input_protein.rfind("/", input_protein.size());
			//if (slash_index1 == -1) slash_index1 = 0;
			//assert(dot_index1 != -1 && "No dot found in filename");
			//std::string protein_name = input_protein.substr(slash_index1 + 1,dot_index1);
			int dot_index2 = input_ligand.rfind(".", input_ligand.size());
			int slash_index2 = input_ligand.rfind("/", input_ligand.size());
			//if (slash_index2 == -1) slash_index2 = 0;
			assert(dot_index2 != -1 && "No dot found in filename");
			std::string ligand_name = input_ligand.substr(slash_index2 +1,dot_index2);

			std::stringstream tagstream;
			tagstream<<ligand_name.c_str()<<"."<<f;
			std::string tag = tagstream.str();
			if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
				tag += "."+option[ OptionKeys::out::output_tag ]();
			}
			std::string complex_filename = "ALIGNEDCOMPLEX_" + tag + ".pdb";

			curr_pose.dump_pdb(complex_filename);

			std::ofstream PLfile(complex_filename.c_str(), std::ios_base::app);
			std::ifstream Lfile(input_ligand.c_str());
			std::string Llineread;

			if ( !Lfile ) {
				std::cout<< "Can't open Ligand-pose file " << input_ligand << std::endl;
				exit(1);
			}
			while ( std::getline(Lfile, Llineread) ) {
				if ( Llineread[0] == 'E' && Llineread[1] == 'N' && Llineread[2] == 'D' ) continue;
				if ( Llineread[0] == 'T' && Llineread[1] == 'E' && Llineread[2] == 'R' ) continue;
				PLfile << Llineread<<"\n";
			}
			PLfile <<"END\n";
			Lfile.close();
			PLfile.close();

		}
		TR << "Done recomputing scores and rmsds" << std::endl;

		outstream.close();
		outstream.clear();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}


