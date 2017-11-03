// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk + dj

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/pockets/PocketConstraint.hh>
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
OPT_KEY( String, contact_list )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );

//set to store pdb info keys
std::set <std::string> interface;
//stores resid of the ligand residue
core::Size lig_res_num;

void define_interface( core::pose::Pose & ref_pose ) {
	lig_res_num =0;
	for ( int j = 1, resnum = ref_pose.size(); j <= resnum; ++j ) {
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
	for ( utility::graph::Graph::EdgeListIter
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

	if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && rsd.atom_is_backbone(atomno);

	return false;
}

Real
calpha_pdb_superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
)
{
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::AtomID::BOGUS_ATOM_ID() );
	for ( Size ii = 1; ii <= mod_pose.size(); ++ii ) {
		if ( ! mod_pose.residue(ii).has("CA") ) continue;
		if ( ! mod_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= ref_pose.size(); ++jj ) {
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

	for ( Size ii = 1; ii <= ref_pose.size(); ++ii ) {
		if ( ! ref_pose.residue(ii).has("CA") ) continue;
		if ( ! ref_pose.residue(ii).is_protein() ) continue;
		for ( Size jj = 1; jj <= mod_pose.size(); ++jj ) {
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
		NEW_OPT ( contact_list, "File name for optional list of contact residues to check","");
		devel::init(argc, argv);

		TR << "Starting recomputing scores and rmsds" << std::endl;

		std::string const ref_decoy_fname = option[ ref_decoy ];


		// create pose from pdb
		pose::Pose ref_pose;
		core::import_pose::pose_from_file( ref_pose, ref_decoy_fname , core::import_pose::PDB_file);

		std::string const cfilename = option[ contact_list ];
		if ( cfilename != "" ) {
			std::ifstream ifs(cfilename.c_str(), std::ifstream::in);
			if ( !ifs.is_open() ) {
				std::cout<< "Error opening contact list file "<<cfilename<<std::endl;
				return -100;
			}
			//ifb.open (cfilename,std::ios::in);
			//std::ostream ios(&ifb);
			std::string intres;
			while ( ifs.good() ) {
				ifs >> intres;
				interface.insert(intres);
			}


		} else {
			define_interface(ref_pose);
		}

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

			core::Real CA_rms = calpha_pdb_superimpose_pose( curr_pose, ref_pose);
			CA_rms = core::scoring::CA_rmsd( curr_pose, ref_pose );
			std::cout << "superimpose to native. Rms to native: " << CA_rms << std::endl;


			//core::Real CA_rms = rmsd_with_super( ref_pose, curr_pose, is_protein_CA );
			//std::cout << "CA_rms to native: " << CA_rms << std::endl;
			//core::Real heavyatom_rms = rmsd_with_super( ref_pose, curr_pose, is_interface_heavyatom );
			core::Real heavyatom_rms = interface_rmsd( ref_pose, curr_pose, is_interface_heavyatom );
			core::Real bbatom_rms = interface_rmsd( ref_pose, curr_pose, is_interface_bbatom );
			std::cout << "Interface rmsd: " << heavyatom_rms << " " << bbatom_rms << std::endl;

			//outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << score_diff << ' ' << p_score_diff <<std::endl;
			outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << bbatom_rms << std::endl;

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


