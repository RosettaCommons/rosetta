// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/concatenate_poses.cc
/// @brief concatenates poses together end-to-end while transforming them to avoid chainbreaks
/// @author frankdt (frankdt@email.unc.edu)

#include <apps/pilot/frankdt/concatenate_poses.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <algorithm>
#include <regex>

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

static basic::Tracer TR( "apps.pilot.frankdt.concatenate_poses" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

FileOptionKey const component_file_name( "component_file" );
FileOptionKey const immobile_domain( "immobile_domain" );
FileOptionKey const output_filename( "output_filename" );
FileOptionKey const binding_partner( "binding_partner" );

namespace apps {
namespace pilot {
namespace frankdt {

concatenate_poses::concatenate_poses():
	utility::pointer::ReferenceCount()
{

}

concatenate_poses::~concatenate_poses(){}

concatenate_poses::concatenate_poses( concatenate_poses const & ) {

}



concatenate_posesOP
concatenate_poses::clone() const {
	return concatenate_posesOP( new concatenate_poses( *this ) );
}

int
main( int argc, char * argv [] )
{
	option.add( component_file_name, "what poses and linkers to join");
	option.add( immobile_domain, "which domain should not move?");
	option.add( output_filename, "output filename");
	option.add( binding_partner, "binding partner");
	devel::init(argc,argv);
	TR << "parsing pose components" << std::endl;
	utility::vector1<core::pose::PoseOP > poses;
	utility::io::izstream pdb_source_file( option[component_file_name].value() );
	utility::vector1<std::pair<core::Size,core::Size>> fixed_spans;
	core::Size immobile_residue=0;
	std::string line;
	std::string new_sequence;
	numeric::HomogeneousTransform< core::Real > immobile_ht;
	core::chemical::ResidueTypeSetCOP res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	while ( getline( pdb_source_file, line) ) {
		if ( line.compare(option[immobile_domain].value())==0 ) {
			TR << "Immobile domain found" << std::endl;
			immobile_residue = new_sequence.length()+1;
		}
		if ( line.size() > 4 && line.substr(line.length()-4,4) == ".pdb" ) {
			std::pair<core::Size,core::Size> new_span;
			poses.push_back(core::import_pose::pose_from_file(line));
			new_span.first = new_sequence.length()+1;
			new_sequence = new_sequence + poses.back()->sequence();
			new_span.second = new_sequence.length();
			fixed_spans.push_back(new_span);
			if ( line.compare(option[immobile_domain].value())==0 ) {
				core::conformation::Residue immobile_basis_residue = poses.back()->residue(1);
				utility::vector1< core::conformation::Atom > & immobile_basis_atoms = immobile_basis_residue.atoms();
				numeric::HomogeneousTransform< core::Real > new_ht( immobile_basis_atoms[ 3 ].xyz(), immobile_basis_atoms[ 1 ].xyz(), immobile_basis_atoms[ 2 ].xyz() );
				immobile_ht = new_ht;
			}
		} else {
			poses.push_back(core::pose::PoseOP(new core::pose::Pose));
			core::pose::make_pose_from_sequence(*(poses.back()), line, *res_type_set,false);
			for ( core::Size i = 1; i <= poses.back()->size(); i++ ) {
				poses.back()->set_phi(i,-180);
				poses.back()->set_psi(i,180);
				poses.back()->set_omega(i,180);
			}
			new_sequence = new_sequence + poses.back()->sequence();
		}
	}



	for ( core::Size current_pose_number = 1; current_pose_number <= poses.size(); current_pose_number++ ) {
		for ( core::Size current_residue_number = 1; current_residue_number<=poses[current_pose_number]->size(); current_residue_number++ ) {
			if ( poses[current_pose_number]->residue(current_residue_number).is_upper_terminus() ) {
				TR <<"found internal terminus" << std::endl;
				core::pose::remove_variant_type_from_pose_residue(*poses[current_pose_number], core::chemical::UPPER_TERMINUS_VARIANT, current_residue_number);
			}
			if ( poses[current_pose_number]->residue(current_residue_number).is_lower_terminus() ) {
				TR <<"found internal terminus" << std::endl;
				core::pose::remove_variant_type_from_pose_residue(*poses[current_pose_number], core::chemical::LOWER_TERMINUS_VARIANT, current_residue_number);
			}
		}
	}

	core::pose::PoseOP working_pose = core::pose::PoseOP(new core::pose::Pose(*poses[1]));

	for ( core::Size current_pose_number = 2; current_pose_number <= poses.size(); current_pose_number++ ) {
		poses[current_pose_number-1]->append_residue_by_bond(poses[current_pose_number]->residue(1),true);
		TR << "Appended residue 1" << std::endl;
		//HomogenousTransform everything to last residue of working pose here

		core::conformation::Residue stationary_basis_residue = poses[current_pose_number-1]->residue(poses[current_pose_number-1]->size());
		utility::vector1< core::conformation::Atom > & stationary_basis_atoms = stationary_basis_residue.atoms();
		//We'll define a coordinate frame for the mobile basis residue
		core::conformation::Residue mobile_basis_residue = poses[current_pose_number]->residue(1);
		utility::vector1< core::conformation::Atom > & mobile_basis_atoms = mobile_basis_residue.atoms();

		//We'll transform all of the residues in all the connected segment to the local coordinate frame of the mobile basis residue and then to the global coordinate frame of the stationary basis residue

		numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_basis_atoms[ 3 ].xyz(), stationary_basis_atoms[ 1 ].xyz(), stationary_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
		numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;
		TR << "Transform_complete 1" << std::endl;
		poses[current_pose_number]->apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());
		working_pose->append_pose_by_jump(*(poses[current_pose_number]),working_pose->size(),"CA","CA");
		//
	}



	working_pose->fold_tree(core::kinematics::FoldTree(working_pose->size()));

	if ( immobile_residue > 0 ) {
		TR <<"Transforming to immobile residue " << immobile_residue << std::endl;
		core::conformation::Residue mobile_basis_residue = working_pose->residue(immobile_residue);
		utility::vector1< core::conformation::Atom > & mobile_basis_atoms = mobile_basis_residue.atoms();
		numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );
		numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
		numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = immobile_ht * inverse_mobile_ht;
		working_pose->apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());
	}

	core::pose::PoseOP output_pose = core::pose::PoseOP(new core::pose::Pose);

	core::pose::make_pose_from_sequence(*output_pose, working_pose->sequence(), *res_type_set,false);

	output_pose->copy_segment(output_pose->size(),*working_pose,1,1);

	std::string fold_tree;
	fold_tree = fold_tree + "EDGE " + std::to_string(immobile_residue) + " 1 -1 ";
	fold_tree = fold_tree + "EDGE " + std::to_string(immobile_residue) + " " + std::to_string(output_pose->size()) + " -1 ";

	if ( option[binding_partner].active() ) {
		TR << "found binding partner" << std::endl;
		core::Size first_residue_of_partner = output_pose->size()+1;
		fold_tree = fold_tree + "EDGE " + std::to_string(immobile_residue) + " " + std::to_string(first_residue_of_partner) + " 1 ";
		output_pose->append_pose_by_jump(*(core::import_pose::pose_from_file(option[binding_partner].value())),output_pose->size(),"CA","CA");
		fold_tree = fold_tree + "EDGE " + std::to_string(first_residue_of_partner) + " " + std::to_string(output_pose->size()) + " -1 ";
	}

	for ( auto span : fixed_spans ) {
		TR << span.first << " " << span.second << std::endl;
	}
	TR << "FOLD TREE" << std::endl;
	TR << fold_tree << std::endl;


	output_pose->dump_pdb(option[output_filename].value());
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
