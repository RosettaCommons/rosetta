// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>


#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <ios>
#include <utility/io/izstream.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>


// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <basic/Tracer.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>

//protocols
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

#include <core/scoring/methods/RG_Energy_Fast.hh>
#include <ObjexxFCL/format.hh>

//#include "james_util.hh" //for calculation of burial
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>


double
reference_sasa(core::chemical::AA aminoacid){
	using namespace core::chemical;
	if ( aminoacid == core::chemical::aa_ala ) return  115.594; // 1 A
	else if ( aminoacid == aa_cys ) return  136.441; // 2 C
	else if ( aminoacid == aa_asp ) return  166.109; // 3 D
	else if ( aminoacid == aa_glu ) return  187.047; // 4 E
	else if ( aminoacid == aa_phe ) return  227.126; // 5 F
	else if ( aminoacid == aa_gly ) return   95.379; // 6 G
	else if ( aminoacid == aa_his ) return  174.877; // 7 H
	else if ( aminoacid == aa_ile ) return  175.089; // 8 I
	else if ( aminoacid == aa_lys ) return  215.766; // 9 K
	else if ( aminoacid == aa_leu ) return  178.404; // 10 L
	else if ( aminoacid == aa_met ) return  187.310; // 11 M
	else if ( aminoacid == aa_asn ) return  142.086; // 12 N
	else if ( aminoacid == aa_pro ) return  146.215; // 13 P
	else if ( aminoacid == aa_gln ) return  180.624; // 14 Q
	else if ( aminoacid == aa_arg ) return  216.548; // 15 R
	else if ( aminoacid == aa_ser ) return  125.275; // 16 S
	else if ( aminoacid == aa_thr ) return  138.319; // 17 T
	else if ( aminoacid == aa_val ) return  153.885; // 18 V
	else if ( aminoacid == aa_trp ) return  229.518; // 19 W
	else if ( aminoacid == aa_tyr ) return  215.013; // 20 Y r++ reference values
	else { return -1; }
}

/**
utility::vector1<double>
sasa_calpha_beta(core::pose::Pose p, double probe_size){
//computes sasa of increased probe size over a protein structure of backbone atoms and C-betas
core::pose::Pose copy(p);
//trim down to c-betas

//calculate sasa

}
**/

/**
void set_lk_dgfree_to_1(core::chemical::ResidueTypeSetOP rsd_set){
core::chemical::AtomTypeSets atom_set = rsd_set->
}
**/

int
main( int argc, char* argv [] )
{
	try {
		using namespace core;
		using namespace scoring;
		using namespace utility;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::id;
		using namespace core::io::silent;
		using namespace ObjexxFCL::format;

		devel::init( argc, argv );

		core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

		vector1<file::FileName> files;
		Real probe_radius = basic::options::option[pose_metrics::sasa_calculator_probe_radius]();
		core::io::silent::SilentFileData sfd;
		if ( basic::options::option[in::file::s].user() ) {
			files=option[in::file::s]();
		} else if ( option[in::file::l].user() ) {

			utility::vector1<file::FileName> list = basic::options::option[ in::file::l ]();
			for ( unsigned int h=1; h<=list.size(); h++ ) {
				utility::io::izstream pdbs(list[h]);
				std::string fname;
				while ( pdbs >> fname ) {
					files.push_back(fname);
				}
			}
		} else if ( option[ in::file::silent ].user() ) {

			std::string silentfilename = option[ in::file::silent ]()[1];
			sfd.set_filename( silentfilename );
			sfd.read_file( silentfilename );
			files = sfd.tags();
		}


		ScoreFunctionOP scfxn = get_score_function();

		// set_lk_dgfree_to_1();

		int width(10); int precision(6);

		core::scoring::methods::RG_Energy_Fast rg;

		std::cout << "SASA resnum sasa fractional_sasa mean_atomic_NO_sasa\n";
		for ( unsigned int f=1; f<=files.size(); f++ ) {
			pose::Pose pose;
			if ( option[ in::file::s ].user() || option[ in::file::l ].user() ) {
				core::import_pose::pose_from_file(pose, files[f], core::import_pose::PDB_file);
			} else {
				SilentStructOP ss = sfd[files[f]];
				ss->fill_pose( pose, *rsd_set );
			}

			(*scfxn)(pose);
			//calculate number of neighbors for each residue as an alternative way of
			//determining burial
			vector1<int> nbrs;
			core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

			core::Real total_num_nbrs = 0;
			nbrs.resize( pose.size() );
			for ( Size i=1; i<= pose.size(); ++i ) {

				{
					nbrs[i] = 1;
					for ( utility::graph::Graph::EdgeListConstIter
							ir  = tenA_neighbor_graph.get_node( i )->const_edge_list_begin(),
							ire = tenA_neighbor_graph.get_node( i )->const_edge_list_end();
							ir != ire; ++ir ) {
						core::Size const neighbor_id( (*ir)->get_other_ind( i ) );
						chemical::ResidueType const & nbr_rsd( pose.residue_type( neighbor_id ) );
						if ( nbr_rsd.is_protein() ) {
							nbrs[i] += 1;
						}
					}
				}
				total_num_nbrs += nbrs[i];
			}
			//end determination of num neighbors
			core::Real total_fa_atr = 0.0;
			core::Real total_fa_sol = 0.0;

			utility::vector1< core::Real > fa_atr;
			utility::vector1< core::Real > fa_sol;
			//get fa-sol and fa-atr
			for ( Size i = 1; i <= pose.size(); i++ ) {
				core::Real res_fa_atr = (pose.energies().residue_total_energies(i))[core::scoring::score_type_from_name("fa_atr")];
				core::Real res_fa_sol = (pose.energies().residue_total_energies(i))[core::scoring::score_type_from_name("fa_sol")];
				fa_atr.push_back(res_fa_atr);
				fa_sol.push_back(res_fa_sol);
				total_fa_atr += res_fa_atr;
				total_fa_sol += res_fa_sol;
			}

			//end get fa-sol and fa-atr


			id::AtomID_Map< Real > atom_sasa;
			utility::vector1< Real > rsd_sasa;
			calc_per_atom_sasa(pose,atom_sasa,rsd_sasa,probe_radius);
			std::cout << "SASA totals " << F( width, precision, calc_total_sasa(pose,probe_radius) )
				<< " " << F( width, precision, total_num_nbrs ) << " "
				<< F( width, precision, rg.calculate_rg_score( pose ) )
				<< " " << F( width, precision, total_fa_atr ) << " " << F( width, precision, total_fa_sol ) << std::endl;
			for ( unsigned int i=1; i<= rsd_sasa.size(); i++ ) {
				std::cout << "SASA " << i << " " << pose.residue(i).name3() << " "
					<< F( width, precision, rsd_sasa[i] ) << " "
					<<  F( width, precision, (rsd_sasa[i]/reference_sasa(pose.aa(i))) )<< " "
					<< nbrs[i] << " " << F( width, precision, fa_atr[i] )<< " " << F( width, precision, fa_sol[i] ) << " " << files[f] << "\n";
			}
		}
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
