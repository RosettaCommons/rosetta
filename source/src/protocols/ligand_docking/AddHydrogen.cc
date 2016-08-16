// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/AddHydrogen.cc
///
/// @brief
/// @Gordon Lemmon

#include <protocols/ligand_docking/AddHydrogen.hh>

#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

#include <core/chemical/PatchOperation.hh>

#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>


namespace protocols {
namespace ligand_docking {


AddHydrogen::AddHydrogen():
	//utility::pointer::ReferenceCount(),
	Mover("AddHydrogen")
{
	Mover::type( "AddHydrogen" );
}

AddHydrogen::AddHydrogen(core::Size const residue_index, core::Size const connection_id):
	//utility::pointer::ReferenceCount(),
	Mover("AddHydrogen"),
	residue_index_(residue_index),
	connection_id_(connection_id)
{
	Mover::type( "AddHydrogen" );
}

AddHydrogen::AddHydrogen(AddHydrogen const & that):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( that ),
	residue_index_(that.residue_index_),
	connection_id_(that.connection_id_)
{}

AddHydrogen::~AddHydrogen() {}

std::string AddHydrogen::get_name() const{
	return "AddHydrogen";
}


void
AddHydrogen::apply( core::pose::Pose & pose )
{
	core::conformation::Residue const & res_to_fix= pose.residue(residue_index_);
	core::chemical::ResidueConnection const & res_connection= res_to_fix.residue_connection(connection_id_);
	core::chemical::AtomICoor const & new_i_coor= res_connection.icoor();

	core::chemical::ResidueTypeOP type_to_fix= res_to_fix.type().clone();
	type_to_fix->name( generate_unique_name() );
	core::Size res_conn_atom_index= type_to_fix->residue_connect_atom_index(connection_id_);
#if defined(WIN32) && !defined(WIN_PYROSETTA)
	core::chemical::AddAtomWIN32 aa(" HH ", "Hapo", "X", 0.09);
#else
	core::chemical::AddAtom aa(" HH ", "Hapo", "X", 0.09);
#endif
	aa.apply(*type_to_fix);
	core::chemical::AddBond ab(" HH ", res_to_fix.atom_name(res_conn_atom_index));
	ab.apply(*type_to_fix);

	core::Size stub_atom1= new_i_coor.stub_atom1().atomno();
	core::Size stub_atom2= new_i_coor.stub_atom2().atomno();
	core::Size stub_atom3= new_i_coor.stub_atom3().atomno();

	std::string name1= res_to_fix.atom_name(stub_atom1);
	std::string name2= res_to_fix.atom_name(stub_atom2);
	std::string name3= res_to_fix.atom_name(stub_atom3);

	core::chemical::SetICoor set_i_coor(
		"HH",/// name this in the style of the other Hs (H1,H2,H3, etc)
		new_i_coor.phi(),
		new_i_coor.theta(),
		1.11, ///TODO Lookup from bond-length table
		name1,
		name2,
		name3
	);
	set_i_coor.apply(*type_to_fix);

	type_to_fix->finalize();
	{
		core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
		core::chemical::ResidueTypeSet & rsd_set= cm->nonconst_residue_type_set( core::chemical::FA_STANDARD );
		rsd_set.add_custom_residue_type(type_to_fix);
	}
	utility::vector1< std::pair< std::string, std::string > > atom_pairs;
	atom_pairs.push_back(std::pair<std::string, std::string>(name1,name1) );
	atom_pairs.push_back(std::pair<std::string, std::string>(name2,name2) );
	atom_pairs.push_back(std::pair<std::string, std::string>(name3,name3) );

	core::conformation::Residue new_res(*type_to_fix, true);
	//type_to_fix_is_good
	pose.replace_residue(residue_index_, new_res, atom_pairs);
}

std::string generate_unique_name(std::string /*input_name*/){
	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::ResidueTypeSet & rsd_set= cm->nonconst_residue_type_set( core::chemical::FA_STANDARD );

	std::string new_name;

	do{
		new_name.clear();
		char a= numeric::random::random_range(65, 90); // ascii range for upper case letters
		char b= numeric::random::random_range(65, 90); // ascii range for upper case letters
		char c= numeric::random::random_range(65, 90); // ascii range for upper case letters


		new_name.append(1,a);
		new_name.append(1,b);
		new_name.append(1,c);

	} while( rsd_set.has_name(new_name));

	return new_name;

}

} // namespace ligand_docking
} // namespace protocols
