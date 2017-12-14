// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file score_protein_ligand_interactions.cc
/// @brief score interactions with a given ligand
/// @author Noah Ollikainen
/// @author Roland A. Pache, PhD

#include <devel/init.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

//C++ Headers
#include <cmath>


static basic::Tracer TR( "pilot_apps.score_protein_ligand_interactions" );

//predicate for sorting by energy
bool min_energy(const std::pair<core::Real, std::string>& x, const std::pair<core::Real, std::string>& y) {
	return x.first < y.first;
}

//main
int main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// initialize Rosetta
		devel::init(argc, argv);

		//load input PDB into pose
		utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
		if ( input_jobs.size() != 1 ) {
			utility_exit_with_message( "Expected exactly one pdb to be specified from the -s or -l flags" );
		}
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, input_jobs[ 1 ]->input_tag() , core::import_pose::PDB_file);
		core::Size const nres( pose.size() );

		//define score function
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		//score pose
		(*score_fxn)(pose);
		score_fxn->show(TR, pose);
		/*score_fxn->set_weight(core::scoring::fa_sol,0);
		score_fxn->set_weight(core::scoring::fa_rep,0);
		score_fxn->show(TR, p);
		*/

		//calculate ligand interaction energies per residue separately for each ligand in the pose
		for ( core::Size i = 1; i <= nres; i++ ) {
			if ( pose.residue(i).is_ligand() ) {
				TR << std::endl << "Ligand: " << pose.residue(i).name3() << std::endl;
				std::vector<std::pair<core::Real, std::string> > ligand_interactions;
				std::vector<std::pair<core::Real, std::string> > ligand_interactions_2b;
				for ( core::Size j = 1; j <= nres; j++ ) {
					if ( i != j && pose.residue(j).is_protein() ) {
						//create new energymap, containing the one-body energies of the given residue and all two-body energies between the residue and ligand
						core::scoring::EMapVector energymap;
						core::scoring::EMapVector energymap_2b;
						//collect context-independent energies
						score_fxn->eval_ci_1b(pose.residue(j), pose, energymap);
						score_fxn->eval_ci_2b(pose.residue(i), pose.residue(j), pose, energymap);
						score_fxn->eval_ci_2b(pose.residue(i), pose.residue(j), pose, energymap_2b);
						//collect context-dependent energies
						score_fxn->eval_cd_1b(pose.residue(j), pose, energymap);
						score_fxn->eval_cd_2b(pose.residue(i), pose.residue(j), pose, energymap);
						score_fxn->eval_cd_2b(pose.residue(i), pose.residue(j), pose, energymap_2b);
						//apply weights to compute interaction energy
						core::Real total_energy = score_fxn->weights().dot(energymap);//
						core::Real two_body_interaction_energy = score_fxn->weights().dot(energymap_2b);
						if ( two_body_interaction_energy != 0.0 ) {
							//collect information about individual energies
							std::string interaction_string = pose.residue(i).name3() + ' '
								+ utility::to_string(pose.pdb_info()->number(i)) + ' '
								+ pose.residue(j).name3() + ' '
								+ utility::to_string(pose.pdb_info()->number(j)) + ' ';
							ligand_interactions.emplace_back(total_energy,interaction_string + utility::to_string(total_energy) + ' ' + energymap.weighted_string_of(score_fxn->weights()));
							ligand_interactions_2b.emplace_back(two_body_interaction_energy,interaction_string + utility::to_string(two_body_interaction_energy) + ' ' + energymap.weighted_string_of(score_fxn->weights()));
						}
					}
				}
				//sort ligand interactions by lowest energy
				std::sort(ligand_interactions.begin(), ligand_interactions.end(), min_energy);
				std::sort(ligand_interactions_2b.begin(), ligand_interactions_2b.end(), min_energy);
				//output ligand interactions
				/*TR << std::endl << "Interacting residues sorted by total energy:" << std::endl;
				for(core::Size index = 0; index < ligand_interactions.size(); index++) {
				TR << index+1 << ": " << ligand_interactions[index].second << std::endl;
				}
				TR << std::endl*/;
				TR << std::endl << "Interacting residues sorted by two-body interaction energy:" << std::endl;
				for ( core::Size index = 0; index < ligand_interactions_2b.size(); index++ ) {
					TR << index+1 << ": " << ligand_interactions_2b[index].second << std::endl;
				}
			}
		}

		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
