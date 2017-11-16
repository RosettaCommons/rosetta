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
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, ref_decoy )

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );

//set to store pdb info keys
std::set <std::string> interface;
//stores resid of the ligand residue

//mjo commenting out 'input_pose' because it is unused and causes a warning
void define_interface( core::pose::Pose & /*input_pose*/ ) {
}

/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);
		pose::Pose input_pose;

		//read in pdb file from command line
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( input_pose, input_pdb_name , core::import_pose::PDB_file);


		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = input_pose.size(); j <= resnum; ++j ) {
			if ( !input_pose.residue(j).is_protein() ) {
				lig_res_num = j;
				break;
			}
			//TR<<input_pose.residue(j).type().name();
			//TR<<" ";
		}
		//TR << "\n";
		if ( lig_res_num == 0 ) {
			TR << "No ligand given in reference PDB structure.  Cannot identify interface."<<std::endl;
			exit (1);
		}

		//TR <<"sele ";

		std::filebuf fb;
		std::stringstream filename;
		filename<<option[ OptionKeys::out::output_tag ]()<<".contacts";
		fb.open (filename.str().c_str(),std::ios::out);
		std::ostream os(&fb);

		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		(*scorefxn)(input_pose);

		EnergyGraph & energy_graph(input_pose.energies().energy_graph());
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

				TR << "resi "<< input_pose.pdb_info()->number(j)<<" or ";

				residuestream << input_pose.pdb_info()->chain(j) << input_pose.pdb_info()->number(j);
				std::string res_id = residuestream.str();
				interface.insert(res_id);
				os<<res_id<<std::endl;

			}
		}


		TR << std::endl;
		TR << lig_res_num<< std::endl;


		input_pose.delete_polymer_residue(lig_res_num);

		/*
		std::string outfname = "score_vs_rmsd.out";
		utility::io::ozstream outstream;
		outstream.open(outfname, std::ios::out);

		outstream << "fname ca_rms allatom_rms relative_score relative_p_constr" << std::endl;

		for (core::Size f=1; f <= basic::options::start_files().size(); f++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(f);
		TR << "Processing decoy " << curr_decoy_fname << std::endl;

		pose::Pose curr_pose;
		core::import_pose::pose_from_file( curr_pose, curr_decoy_fname , core::import_pose::PDB_file);

		// This is the residue we'll use for PocketConstraint
		central_relax_res = 0;
		for ( int j = 1, resnum = curr_pose.size(); j <= resnum; ++j ) {
		if ( curr_pose.pdb_info()->number(j) == central_relax_pdb_number ) {
		central_relax_res = j;
		}
		}
		TR << "Set constraint to 0" << std::endl;
		scorefxn->set_weight( core::scoring::pocket_constraint, 0 );
		(*scorefxn)(curr_pose);

		TR << "Set constraint to 1" << std::endl;
		core::Real score_diff = curr_pose.energies().total_energies()[ total_score ] - ref_score;
		scorefxn->set_weight( core::scoring::pocket_constraint, 1 );

		curr_pose.add_constraint( new core::scoring::constraints::PocketConstraint( curr_pose) );

		// rescore, report new score
		(*scorefxn)(curr_pose);
		TR << "Done rescoring" << std::endl;
		core::Real p_score_diff = curr_pose.energies().total_energies()[ pocket_constraint ];
		core::Real CA_rms = rmsd_with_super( input_pose, curr_pose, is_protein_CA );
		core::Real heavyatom_rms = rmsd_with_super( input_pose, curr_pose, is_interface_heavyatom );

		outstream << curr_decoy_fname << ' ' << CA_rms << ' ' << heavyatom_rms << ' ' << score_diff << ' ' << p_score_diff <<std::endl;

		}

		TR << "Done recomputing scores and rmsds" << std::endl;

		outstream.close();
		outstream.clear();
		*/

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}


