// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/steric_fusion_scan.cc
/// @brief concatenates poses together end-to-end while transforming them to avoid chainbreaks
/// @author frankdt (frankdt@email.unc.edu)

#include <apps/pilot/frankdt/steric_fusion_scan.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <algorithm>
#include <regex>
#include <boost/algorithm/string.hpp>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <core/conformation/Atom.hh>
#include <numeric/random/random.hh>
#include <numeric/HomogeneousTransform.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.frankdt.steric_fusion_scan" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

FileOptionKey const n_terminal_component_list( "n_terminal_component_list" );
FileOptionKey const c_terminal_component_list( "c_terminal_component_list" );
FileOptionKey const exclude_component_file_name( "exclude_component_file_name" );
FileOptionKey const output_filename( "output_file_name" );
FileOptionKey const n_terminal_range( "n_terminal_range" );
FileOptionKey const c_terminal_range( "c_terminal_range" );
//FileOptionKey const minimum_good_clashes( "minimum_good_clashes" );

namespace apps {
namespace pilot {
namespace frankdt {

steric_fusion_scan::steric_fusion_scan():
	utility::pointer::ReferenceCount()
{

}

steric_fusion_scan::~steric_fusion_scan(){}

steric_fusion_scan::steric_fusion_scan( steric_fusion_scan const & ) {

}



steric_fusion_scanOP
steric_fusion_scan::clone() const {
	return steric_fusion_scanOP( new steric_fusion_scan( *this ) );
}

int
main( int argc, char * argv [] )
{
	option.add( n_terminal_component_list, "filename containing the n-terminal PDBs");
	option.add( c_terminal_component_list, "filename containing the c-terminal PDBs");
	option.add( exclude_component_file_name, "filename containing structure to be sterically blocked");
	option.add( n_terminal_range, "what residues on the N-terminal component's C terminus to scan through, counted backwards from the terminus");
	option.add( c_terminal_range, "what residues on the C-terminal component's N terminus to scan through");
	option.add( output_filename, "what to name the output file");
	//option.add( minimum_good_clashes, "how few clashes with the blocked structure a fusion needs to count");
	core::Real clash_radius = 4.0;
	core::Real interaction_radius = 6.0;
	core::Size max_clashes = 0;
	core::Size min_good_clashes = 1; //std::stoi(option[minimum_good_clashes].value());
	core::Size min_interactions = 5; //std::stoi(option[minimum_good_clashes].value());
	devel::init(argc,argv);
	core::chemical::ResidueTypeSetCOP res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	utility::vector1<core::Size> best_interactions;
	utility::vector1<core::Size> best_good_clashes;
	utility::vector1<core::Size> best_interaction_n_terminal_resnums;
	utility::vector1<core::Size> best_interaction_c_terminal_resnums;
	utility::vector1<std::string> best_n_lines;
	utility::vector1<std::string> best_c_lines;
	TR << "Initialized" << std::endl;
	//
	std::string n_line;
	std::string c_line;
	core::pose::PoseOP n_terminal;
	core::pose::PoseOP c_terminal;
	core::pose::PoseOP exclude = core::import_pose::pose_from_file(option[exclude_component_file_name].value());

	std::string line;
	utility::vector1<std::string> tokens;
	utility::io::izstream n_terminal_pdb_source_file( option[n_terminal_component_list].value() );
	utility::io::izstream c_terminal_pdb_source_file( option[c_terminal_component_list].value() );
	while ( getline( n_terminal_pdb_source_file, n_line) ) {
		n_terminal = core::import_pose::pose_from_file(n_line);
		// now the n ranges
		line = option[n_terminal_range].value();
		boost::split( tokens, line, boost::is_any_of(",") );
		core::Size n_terminal_min = 1 + n_terminal->size() - std::stoi(tokens[2]);
		core::Size n_terminal_max = 1 + n_terminal->size() - std::stoi(tokens[1]);
		while ( getline( c_terminal_pdb_source_file, c_line) ) {
			c_terminal = core::import_pose::pose_from_file(c_line);
			// and the c ranges
			line = option[c_terminal_range].value();
			boost::split( tokens, line, boost::is_any_of(",") );
			core::Size c_terminal_min = std::stoi(tokens[1]);
			core::Size c_terminal_max = std::stoi(tokens[2]);
			//TR << "imported N and C poses" << std::endl;

			for ( core::Size n_terminal_resnum = n_terminal_min; n_terminal_resnum <= n_terminal_max; n_terminal_resnum++ ) {
				for ( core::Size c_terminal_resnum = c_terminal_min; c_terminal_resnum <= c_terminal_max; c_terminal_resnum++ ) {

					//TR <<"Superimposing N-term residue " << n_terminal_resnum << " on C-term residue " << c_terminal_resnum << std::endl;

					core::conformation::Residue stationary_basis_residue = c_terminal->residue(c_terminal_resnum);
					utility::vector1< core::conformation::Atom > & stationary_basis_atoms = stationary_basis_residue.atoms();
					//We'll define a coordinate frame for the mobile basis residue
					core::conformation::Residue mobile_basis_residue = n_terminal->residue(n_terminal_resnum);
					utility::vector1< core::conformation::Atom > & mobile_basis_atoms = mobile_basis_residue.atoms();
					//We'll transform all of the residues in all the connected segment to the local coordinate frame of the mobile basis residue and then to the global coordinate frame of the stationary basis residue
					numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_basis_atoms[ 3 ].xyz(), stationary_basis_atoms[ 1 ].xyz(), stationary_basis_atoms[ 2 ].xyz() );
					numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );
					numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
					numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;
					n_terminal->apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());
					//transform is made
					core::Size clashes = 0;
					core::Size good_clashes = 0;
					core::Size interactions = 0;
					for ( core::Size n_terminal_component_residue = 1; n_terminal_component_residue < n_terminal_resnum; n_terminal_component_residue++ ) {
						core::conformation::Atom n_term_CA = n_terminal->residue(n_terminal_component_residue).atom(2);
						for ( core::Size c_terminal_component_residue = c_terminal_resnum + 2; c_terminal_component_residue <= c_terminal->size(); c_terminal_component_residue++ ) {

							core::conformation::Atom c_term_CA = c_terminal->residue(c_terminal_component_residue).atom(2);

							core::Real current_distance = n_term_CA.xyz().distance(c_term_CA.xyz());
							if ( current_distance < clash_radius ) {
								clashes++;
							} else if ( current_distance < interaction_radius ) {
								interactions++;
							}

							if ( clashes > max_clashes ) {
								TR << "Too many clashes, aborting this fusion." << std::endl;
								interactions = 0;
								n_terminal_component_residue = n_terminal_resnum;
								c_terminal_component_residue = c_terminal->size()+1;
							}

						}
						// and the exclude check
						for ( core::Size exclude_component_residue = 1; exclude_component_residue <= exclude->size(); exclude_component_residue++ ) {
							core::conformation::Atom exclude_CA = exclude->residue(exclude_component_residue).atom(2);
							core::Real current_distance = n_term_CA.xyz().distance(exclude_CA.xyz());
							if ( current_distance < clash_radius ) {
								good_clashes++;
							}
						}
					}
					//TR << "Interactions Detected: " << interactions << std::endl;
					//TR << "Good Clashes: " << good_clashes << std::endl;
					if ( good_clashes >= min_good_clashes && interactions >= min_interactions ) {
						best_interactions.push_back(interactions);
						best_good_clashes.push_back(good_clashes);
						best_interaction_n_terminal_resnums.push_back(n_terminal_resnum);
						best_interaction_c_terminal_resnums.push_back(c_terminal_resnum);
						best_n_lines.push_back(n_line);
						best_c_lines.push_back(c_line);
						TR << "Interactions: " << best_interactions.back() << std::endl;
						TR << "Good Clashes: " << best_good_clashes.back() << std::endl;
						TR << "C Term Residue: " << best_interaction_n_terminal_resnums.back() << std::endl;
						TR << "N Term Residue: " << best_interaction_c_terminal_resnums.back() << std::endl;
					}
				}
			}
		}
	}



	//now write out the output structure

	if ( best_interactions.size() == 0 ) {
		TR << "No good fusion found." << std::endl;
		return 0;
	}

	for ( core::Size current_entry = 1; current_entry <= best_interactions.size(); current_entry++ ) {
		n_terminal = core::import_pose::pose_from_file(best_n_lines[current_entry]);
		c_terminal = core::import_pose::pose_from_file(best_c_lines[current_entry]);
		core::conformation::Residue stationary_basis_residue = c_terminal->residue(best_interaction_c_terminal_resnums[current_entry]);
		utility::vector1< core::conformation::Atom > & stationary_basis_atoms = stationary_basis_residue.atoms();
		//We'll define a coordinate frame for the mobile basis residue
		core::conformation::Residue mobile_basis_residue = n_terminal->residue(best_interaction_n_terminal_resnums[current_entry]);
		utility::vector1< core::conformation::Atom > & mobile_basis_atoms = mobile_basis_residue.atoms();
		numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_basis_atoms[ 3 ].xyz(), stationary_basis_atoms[ 1 ].xyz(), stationary_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
		numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;
		n_terminal->apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());

		core::pose::PoseOP output_pose = core::pose::PoseOP(new core::pose::Pose(*n_terminal,1,best_interaction_n_terminal_resnums[current_entry]));
		core::pose::PoseOP append_pose = core::pose::PoseOP(new core::pose::Pose(*c_terminal,best_interaction_c_terminal_resnums[current_entry]+1,c_terminal->size()));
		output_pose->append_pose_by_jump(*append_pose,output_pose->size(),"CA","CA");
		output_pose->fold_tree(core::kinematics::FoldTree(output_pose->size()));

		core::pose::PoseOP final_pose = core::pose::PoseOP(new core::pose::Pose());

		core::pose::make_pose_from_sequence(*final_pose, output_pose->sequence(), *res_type_set,false);

		final_pose->copy_segment(final_pose->size(),*output_pose,1,1);

		std::string full_output_file_name = option[output_filename].value();
		full_output_file_name = full_output_file_name + "_" + std::to_string(current_entry) + ".pdb";

		final_pose->dump_pdb(full_output_file_name);
	}
	return 0;
}


} //apps
} //pilot
} //frankdt

int
main( int argc, char * argv [] )
{
	apps::pilot::frankdt::main( argc, argv );
	return 0;
}
