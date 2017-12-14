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

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );

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
		auto * edge( static_cast< EnergyEdge *> (*iru) );
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

	/* core::scoring::TwelveANeighborGraph const & graph = ref_pose.energies().twelveA_neighbor_graph();
	for ( utility::graph::Graph::EdgeListConstIter
	iter = graph.get_node( lig_res_num )->const_edge_list_begin(),
	iter_end = graph.get_node( lig_res_num )->const_edge_list_end();
	iter != iter_end; ++iter ) {
	Size const neighbor_id( (*iter)->get_other_ind( lig_res_num ) );

	// create string id to store in set
	std::ostringstream residuestream;

	TR << "resi "<< ref_pose.pdb_info()->number(neighbor_id)<<" or ";

	residuestream << ref_pose.pdb_info()->chain(neighbor_id) << ref_pose.pdb_info()->number(neighbor_id);
	std::string res_id = residuestream.str();
	interface.insert(res_id);
	//allow_moving.at(neighbor_id) = true;
	}

	*/

	TR << std::endl;
	TR << lig_res_num<< std::endl;


	ref_pose.delete_polymer_residue(lig_res_num);
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

		// This is the residue we'll use for the pocket constraint
		/*        std::string resid(option[ OptionKeys::pocket_grid::central_relax_pdb_num ]);
		int  central_relax_pdb_number;
		char chain= ' ';
		core::Size central_relax_res = 0;
		std::size_t fpos( resid.find(':') );
		if ( fpos != std::string::npos ) {
		central_relax_pdb_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
		if (fpos != resid.size()-1 ) {
		chain = resid[ fpos+1 ];
		}
		} else {
		central_relax_pdb_number = ObjexxFCL::int_of( resid );
		}
		central_relax_res = 0;
		for ( int j = 1, resnum = ref_pose.size(); j <= resnum; ++j ) {
		if ( ref_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
		if (chain != ' '){
		if ( ref_pose.pdb_info()->chain(j) == chain ) {
		central_relax_res = j;
		}
		}else{
		central_relax_res = j;
		}
		}
		}

		*/


		/*
		int const central_relax_pdb_number = option[ OptionKeys::pocket_grid::central_relax_pdb_num ];
		//int const central_relax_pdb_number = option[ central_relax_pdb_num ];
		for ( int j = 1, resnum = ref_pose.size(); j <= resnum; ++j ) {
		if ( ref_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
		central_relax_res = j;
		}
		}
		*/
		//if ( central_relax_res == 0 ) {
		//        std::cerr << "ERROR!! Could not find residue measure pocket constraint" << std::endl;
		//        exit(1);
		//}

		// scoring function
		//scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
		//scoring::ScoreFunctionOP scorefxn( get_score_function() );

		// PocketConstraint set to unweighted (1)
		//scorefxn->set_weight( core::scoring::pocket_constraint, 1 );

		// Calculate score without Pocket Constraint
		//(*scorefxn)(ref_pose);
		//core::Real ref_score = ref_pose.energies().total_energies()[ total_score ];
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
			define_interface( ref_pose );
		}

		TR << "Defined interface" << std::endl;

		// Add pocket constraint, rescore
		//TR << "Add pocket scoring ref_pose" << std::endl;
		//core::scoring::constraints::PocketConstraintOP tmp( new core::scoring::constraints::PocketConstraint( ref_pose, central_relax_res ));
		//ref_pose.add_constraint( tmp );
		//ref_pose.add_constraint( new core::scoring::constraints::PocketConstraint( ref_pose) );
		//TR << "rescore" << std::endl;
		//(*scorefxn)(ref_pose);
		//core::Real const starting_pocket_score = ref_pose.energies().total_energies()[ pocket_constraint ];

		//TR << "P score: " << starting_pocket_score<< std::endl;
		// ref_pose.remove_constraint( tmp );

		//scorefxn->set_weight( core::scoring::pocket_constraint, 0 );
		// Open output file, generate the header line (save it for printing in the log later), print to file
		std::stringstream filename;
		if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
			filename << "score_vs_rmsd." << option[ OptionKeys::out::output_tag ]()<<".out";
		} else {
			filename<<"score_vs_rmsd.out";
		}
		utility::io::ozstream outstream;
		outstream.open(filename.str().c_str(), std::ios::out);

		outstream << "fname ca_rms allatom_rms relative_score relative_p_constr" << std::endl;

		for ( core::Size f=1; f <= basic::options::start_files().size(); f++ ) {

			std::string const curr_decoy_fname = basic::options::start_files().at(f);
			TR << "Processing decoy " << curr_decoy_fname << std::endl;

			pose::Pose curr_pose;
			core::import_pose::pose_from_file( curr_pose, curr_decoy_fname , core::import_pose::PDB_file);

			// This is the residue we'll use for PocketConstraint
			//central_relax_res = 0;
			//for ( int j = 1, resnum = curr_pose.size(); j <= resnum; ++j ) {
			//if ( curr_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
			//        central_relax_res = j;
			//}
			//}
			//TR << "Set constraint to 0" << std::endl;
			//scorefxn->set_weight( core::scoring::pocket_constraint, 0 );
			//(*scorefxn)(curr_pose);

			//TR << "Set constraint to 1" << std::endl;
			//core::Real score_diff = curr_pose.energies().total_energies()[ total_score ] - ref_score;
			//scorefxn->set_weight( core::scoring::pocket_constraint, 1 );

			//curr_pose.add_constraint( new protocols::constraints_additional::PocketConstraint( curr_pose) );

			// rescore, report new score
			//(*scorefxn)(curr_pose);
			//TR << "Done rescoring" << std::endl;
			//core::Real p_score_diff = curr_pose.energies().total_energies()[ pocket_constraint ];
			core::Real CA_rms = rmsd_with_super( ref_pose, curr_pose, is_protein_CA );
			core::Real heavyatom_rms = rmsd_with_super( ref_pose, curr_pose, is_interface_heavyatom );
			core::Real bb_rms = rmsd_with_super( ref_pose, curr_pose, is_interface_bbatom );

			//outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << score_diff << ' ' << p_score_diff <<std::endl;
			outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << bb_rms << std::endl;

		}

		TR << "Done recomputing scores and rmsds" << std::endl;

		outstream.close();
		outstream.clear();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}


