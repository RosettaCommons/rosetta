// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metalloproteins/util.cc
/// @brief  Utilities for working with metalloproteins.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


// Unit header
#include <core/util/metalloproteins_util.hh>

// Package headers
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/util.tmpl.hh>

// Project headers

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <stdio.h>

// External headers
#include <ObjexxFCL/string.functions.hh>
//#include <boost/functional/hash.hpp>
//#include <boost/foreach.hpp>


#define foreach BOOST_FOREACH


namespace core {
namespace util {

	static basic::Tracer TR("core.util.metalloproteins_util");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief: Adds an arbitrary covalent linkage between two atoms (resA_At and resB_At) in two residues (at positions resA_pos and resB_pos).
/// @details:  This is useful for adding covalent linkages between metal-binding side-chains and metal atoms.  This code was shamelessly
/// stolen from Florian's EnzConstraintParameters.cc in protocols/toolbox/match_enzdes_utils, and was modified to permit deletion of
/// unnecessary protons.  NOTE: THIS CODE MODIFIES THE RESIDUE TYPE LIST, AND IS CURRENTLY NOT THREADSAFE.
/// @author:  Vikram K. Mulligan (vmullig@uw.edu), Florian Richter (flosopher@gmail.com)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkage(
	core::pose::Pose & pose,
	Size const resA_pos,
	Size const resB_pos,
	Size const resA_At,
	Size const resB_At,
	bool const remove_hydrogens //Should extraneous hydrogens on the bonding atoms be removed?
)
{
	using namespace core::chemical;

	TR << "Adding covalent linkage between residue " << resA_pos << "'s " << pose.residue(resA_pos).atom_name(resA_At).c_str() <<
			" atom and residue " << resB_pos << "'s " << pose.residue(resB_pos).atom_name(resB_At).c_str() << " atom." << std::endl;

	//if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
	//	pose.dump_pdb("bef_resmod.pdb");
	//}
	std::string resA_base = residue_type_base_name( pose.residue_type(resA_pos) );
	std::string resB_base = residue_type_base_name( pose.residue_type(resB_pos) );
	std::string resA_var, resB_var;

	//std::string resA_type_set = pose.residue(resA_pos).residue_type_set().name();
	//std::string resB_type_set = pose.residue(resB_pos).residue_type_set().name();

	add_covalent_linkage_helper( pose, /*resA_,*/ resA_pos, resA_At, pose.residue(resB_pos).xyz(resB_At), /*core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( resA_type_set ), ntorsionA_, nangleA_, disAB_,*/ resA_var, remove_hydrogens);
	add_covalent_linkage_helper( pose, /*resB_,*/ resB_pos, resB_At, pose.residue(resA_pos).xyz(resA_At), /*core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( resB_type_set ), ntorsionB_, nangleB_, disAB_,*/ resB_var, remove_hydrogens);

	//if(basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ){
	//	pose.dump_pdb("after_resmod.pdb");
	//}

	std::string resA_atomname = pose.residue( resA_pos ).atom_name( resA_At /*(resA_->get_template_atoms_at_pos(pose, resA_pos))->atom1_[resA_At].atomno() */);
	std::string resB_atomname = pose.residue( resB_pos ).atom_name( resB_At /*(resB_->get_template_atoms_at_pos(pose, resB_pos))->atom1_[resB_At].atomno() */);

	//TR.Debug << "Adding chemical bond between " << pose.residue( resA_pos ).name() << " " << (resA_->get_template_atoms_at_pos(pose, resA_pos))->atom1_[resA_At].atomno() <<  " "<< resA_pos << " " << resA_atomname << " and "
	//				 << pose.residue( resB_pos ).name() << " " << resB_pos << " " << (resB_->get_template_atoms_at_pos(pose, resB_pos))->atom1_[resB_At].atomno() << " " << resB_atomname << std::endl;

	pose.conformation().declare_chemical_bond(
		resA_pos, resA_atomname,
		resB_pos, resB_atomname
	);

	//EnzdesCstParamCacheOP param_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block_ ) );
	//param_cache->covalent_connections_.push_back( new CovalentConnectionReplaceInfo(resA_base, resB_base, resA_var, resB_var, resA_pos, resB_pos, restype_set_ ) ); //new

} //add_covalent_linkage


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief: This is a helper function for the add_covalent_linkage function.  You probably don't want to call it directly, unless you
/// really know what you're doing.
/// @details:  This is useful for adding covalent linkages between metal-binding side-chains and metal atoms.  This code was shamelessly
/// stolen from Florian's EnzConstraintParameters.cc (colourful comments and all) in protocols/toolbox/match_enzdes_utils, and was
/// modified to permit deletion of unnecessary protons and to remove EnzDes-specific stuff.  NOTE: THIS CODE MODIFIES THE RESIDUE TYPE
/// LIST, AND IS CURRENTLY NOT THREADSAFE.
/// @author:  Vikram K. Mulligan (vmullig@uw.edu), Florian Richter (flosopher@gmail.com)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkage_helper(
	core::pose::Pose & pose,
	//EnzCstTemplateResOP template_res,
	core::Size const res_pos,
	core::Size const Atpos,
	numeric::xyzVector < core::Real > const partner_xyz, //Coordinates of the atom to which this will be bonded
	//core::Real itorsion,
	//core::Real iangle,
	//core::Real idis,
 	std::string & res_varname,
	bool const remove_hydrogens
 )
{
	//std::cout << "APL DEBUG EnzConstraintParameters.cc::make_constraint_covalent_helper begin" << std::endl;

	using namespace core::chemical;
	using namespace core::pack::dunbrack;

	std::string res_atom = pose.residue(res_pos).atom_name( Atpos /*(template_res->get_template_atoms_at_pos(pose, res_pos) )->atom1_[Atpos].atomno()*/ );
	core::Size res_atom_parent_index = pose.residue(res_pos).atom_base(Atpos);
	core::Size res_atom_parent_parent_index = pose.residue(res_pos).atom_base(res_atom_parent_index);
	//Special cases: Atpos == 1 or Atpos == 2
	if(Atpos==1) {
		res_atom_parent_index=2; res_atom_parent_parent_index=3;
	} else if (Atpos==2) {
		res_atom_parent_index=1; res_atom_parent_parent_index=3;
	}
	std::string res_atom_parent1 = pose.residue(res_pos).atom_name( res_atom_parent_index );
	std::string res_atom_parent2 = pose.residue(res_pos).atom_name( res_atom_parent_parent_index );


	//Calculate icoor information for the bond.  The connection point should be placed where the other atom goes.
	numeric::xyzVector < core::Real > Atxyz = pose.residue(res_pos).xyz(Atpos);
	numeric::xyzVector < core::Real > Atpxyz = pose.residue(res_pos).xyz(res_atom_parent_index);
	numeric::xyzVector < core::Real > Atppxyz = pose.residue(res_pos).xyz(res_atom_parent_parent_index);

	core::Real idis (Atxyz.distance(partner_xyz)); //The distance to the bonding partner.
	core::Real iangle ( numeric::angle_radians ( partner_xyz, Atxyz, Atpxyz) ); //The bond angle between the bonding atom, this atom, and its parent atom.
	core::Real itorsion ( numeric::dihedral_radians ( partner_xyz, Atxyz, Atpxyz, Atppxyz) ); //The torsion angle between the bonding atom, this atom, this atom's parent, and this atom's grandparent.

	//need to remove whitespace, why the f is this so clumsy in c++?!?
	int whitespace_pos = res_atom.find(" ");
	while( whitespace_pos != -1 ) {
		res_atom.erase(whitespace_pos, 1 );
		whitespace_pos = res_atom.find(" ");
	}
	whitespace_pos = res_atom_parent1.find(" ");
	while( whitespace_pos != -1 ) {
		res_atom_parent1.erase(whitespace_pos, 1 );
		whitespace_pos = res_atom_parent1.find(" ");
	}
	whitespace_pos = res_atom_parent2.find(" ");
	while( whitespace_pos != -1 ) {
		res_atom_parent2.erase(whitespace_pos, 1 );
		whitespace_pos = res_atom_parent2.find(" ");
	}

	std::string current_pose_type_basename( residue_type_base_name( pose.residue_type(res_pos) ) );
	std::string current_pose_type_patches_name( residue_type_all_patches_name( pose.residue_type(res_pos) ) );

	res_varname = "_connect" + res_atom;
	{// scope
		// Find a name for the new residue type / variant name that will be added to the existing
		// residue so that, if the existing residue already has this variant type, then the
		// new residue type will get one with a new name.
		core::chemical::ResidueType const & currres( pose.residue_type( res_pos ));
		Size count=0;
		while ( true ) {
			if ( count > 1000 ) {
				utility_exit_with_message( "Encountered infinite loop trying to find a new variant name for residue type " + currres.name() + " in EnzConstraintParameters.  Talk to Andrew.");
			}
			++count;
			if ( count == 1 ) {
				if ( ! currres.has_variant_type( res_varname )) break;
			} else {
				res_varname = "_"+utility::to_string(count)+"connect"+res_atom;
				if ( ! currres.has_variant_type( res_varname )) break;
			}
		}
	}
	std::string res_type_mod_name( current_pose_type_basename + res_varname + current_pose_type_patches_name );

	//check whether the modified residues have already been created earlier
	if( !pose.residue(res_pos).residue_type_set().has_name(res_type_mod_name) ){

		//holy jesus, we have to change the residue type set.
		//we not only need to clone, modify and add  the residue type
		//currently in the pose, but all similar ones (same basename )
		//in the ResidueTypeSet. We also need to create a copy of the
		// SingleLigandRotamerLibrary in case it exists.
		//reminds me of open heart surgery...

		//the following line is necessary to ensure that a ligand rotamer library exists
		//if this function is called before any scoring happened
		RotamerLibrary::get_instance().get_rsd_library( pose.residue_type( res_pos ));

		SingleLigandRotamerLibraryOP new_lrots = NULL;
		if( pose.residue_type(res_pos).is_ligand() &&
				RotamerLibrary::get_instance().rsd_library_already_loaded( pose.residue_type(res_pos) ) ) {
			new_lrots = new SingleLigandRotamerLibrary();
		}

		utility::pointer::access_ptr< ResidueTypeSet > mod_restype_set = & ChemicalManager::get_instance()->nonconst_residue_type_set( pose.residue(res_pos).residue_type_set().name() );

		//first get all residue types that correspond to the type in question
		ResidueTypeCOPs res_to_modify = mod_restype_set->name3_map( pose.residue_type(res_pos).name3() );

		for ( utility::vector1< ResidueTypeCOP >::iterator res_it = res_to_modify.begin(); res_it != res_to_modify.end(); ++res_it) {
			std::string const base_name( residue_type_base_name( *(*res_it) ) );
			//TR << "contemplating modification of residuetype " << (*res_it)->name() << " with basename " << base_name << std::endl; TR.flush(); //DELETE ME

			if( current_pose_type_basename == base_name ){

				ResidueTypeOP mod_res;
				core::Size con_res(0);
				//TR << " MODIFYING" << std::endl; TR.flush(); //DELETE ME
				std::string patches_name( residue_type_all_patches_name( *(*res_it) ) );
				std::string new_name( base_name + res_varname + patches_name );

				mod_res = (*res_it)->clone();
				con_res = mod_res->add_residue_connection( res_atom );

				mod_res->name( new_name );
				assert( ! mod_res->has_variant_type( res_varname ) );
				mod_res->add_variant_type( res_varname ); //necessary to restrict the packer to only use this residue variant in packing

				mod_res->set_icoor( "CONN"+ObjexxFCL::string_of( con_res ), itorsion, iangle, idis, res_atom, res_atom_parent1, res_atom_parent2, true );

				if(mod_res->is_metal()) {  //If this is a metal and I have fewer virts attached to the metal than there are connections, add virts.  Tethering the metal requires one virt per connection.
					core::Size virtcount = mod_res->n_virtual_atoms();
					while(true) {
						if(virtcount >= con_res) break;
						virtcount++;

						//Add a virt atom.  First, find a unique name for it:
						std::string virtname = "";
						core::Size virtnum = 1;
						while(true) {
							std::stringstream virtcandidate;
							virtcandidate << "V" << virtnum;
							if( !mod_res->has(virtcandidate.str()) ) {
								virtname = virtcandidate.str();
								break;
							}
							++virtnum;
						} //Finding unique name

						//Create a new virt atom, using the unique name found above:
#ifdef WIN32
						core::chemical::AddAtomWIN32 addvirt(virtname, "VIRT", "VIRT", 0.0);
#else
						core::chemical::AddAtom addvirt(virtname, "VIRT", "VIRT", 0.0);
#endif
						addvirt.apply( (*mod_res) );
						mod_res->add_bond(virtname, res_atom);
						mod_res->set_atom_base(virtname, res_atom);
						//TR << "Setting metal itorsion=" << itorsion << " iangle=" << iangle << " idis=" << idis << std::endl; TR.flush(); //DELETE ME
						mod_res->set_icoor( virtname,  itorsion, iangle, idis, res_atom, res_atom_parent1, res_atom_parent2, true ); //Place the virt atop the bonding partner.
						mod_res->finalize();
					}
				}

				if(remove_hydrogens) { //If we're stripping off hydrogens, loop through all the hydrogens and remove ONLY the first bonded hydrogen that we encounter.  Adjust charges accordingly.
					utility::vector1 < core::Size > bonded_indices = mod_res->bonded_neighbor( Atpos );
					mod_res->finalize(); //Extra finalize needed prior to atom_is_hydrogen() call, below.
					for(core::Size ia=1, iamax=bonded_indices.size(); ia<=iamax; ++ia) {
						//printf("This atom is bonded to index %lu.\n", bonded_indices[ia]); fflush(stdout); //DELETE ME
						if(mod_res->atom_is_hydrogen(bonded_indices[ia])) {
							std::string hname = mod_res->atom_name( bonded_indices[ia] );
							core::chemical::SetAtomicCharge setatchg( res_atom, (mod_res->atom(res_atom)).charge() - 1.0 + (mod_res->atom(bonded_indices[ia])).charge() ); //Adjust the charge
							setatchg.apply( (*mod_res) );
							core::chemical::SetAtomicCharge setatchg2 ( hname, 0.0); //Nullify the charge on the proton
							setatchg2.apply( (*mod_res) );
							core::chemical::SetAtomType set_virt1 ( hname, "VIRT");  //Make the proton to delete into a VIRT (probably less likely to break things than outright deleting an atom).
							core::chemical::SetMMAtomType set_virt2 ( hname, "VIRT");
							set_virt1.apply( (*mod_res) );
							set_virt2.apply( (*mod_res) );
							//mod_res->delete_atom( hname ); //Delete the hydrogen
							mod_res->finalize();
							//printf("Deleted a proton.\n"); fflush(stdout); //DELETE ME
							break; //Just delete one hydrogen!
						}
					}
				}

				//new_lrots is empty at the moment, but will be filled a couple of lines down
				if( pose.residue_type( res_pos ).is_ligand() ) {
					RotamerLibrary::get_instance().add_residue_library( *mod_res, new_lrots );
				}

				//finalize again just to make sure
				mod_res->nondefault(true);
				mod_res->base_restype_name( base_name );
				mod_res->finalize();

				mod_restype_set->add_residue_type( mod_res );
			}
		}

		//and last but not least we have to regenerate the rotamer library for the ligand
		if( pose.residue_type( res_pos ).is_ligand() ) {

			SingleLigandRotamerLibraryCAP old_lrots(
				static_cast< SingleLigandRotamerLibrary const * >
				( RotamerLibrary::get_instance().get_rsd_library( pose.residue_type( res_pos ))() ));

			if( old_lrots != 0 ){

				using namespace core::conformation;
				utility::vector1< ResidueOP > new_rotamers;
				new_rotamers.clear();
				utility::vector1< ResidueOP > const old_rotamers = old_lrots->get_rotamers();

				//std::cerr << "old rotamer library has " << old_rotamers.size() << "members";
				for( utility::vector1< ResidueOP>::const_iterator oldrot_it = old_rotamers.begin(); oldrot_it != old_rotamers.end(); ++oldrot_it){
					ResidueOP new_rot_res = new Residue( pose.residue(res_pos).residue_type_set().name_map(res_type_mod_name), true);
					//set the coordinates
					//1. we go over the atoms of the NEW residue on purpose, to make sure that no atom gets skipped
					for( core::Size at_ct = 1; at_ct <= new_rot_res->natoms(); at_ct++){
						if( !(*oldrot_it)->has( new_rot_res->atom_name( at_ct ) ) ){
							std::cerr << "Unexpected ERROR: when regenerating ligand rotamer library (for covalent constraints), one atom wasn't found in a template rotamer." << std::endl;
							utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
						}
						else{
							new_rot_res->set_xyz( at_ct, (*oldrot_it)->xyz( new_rot_res->atom_name( at_ct ) ) );
						}
					}
					//2. we also set the chis, to make sure everything is properly in place
					new_rot_res->chi( (*oldrot_it)->chi() );

					new_rotamers.push_back( new_rot_res );
				}
				new_lrots->set_reference_energy( old_lrots->get_reference_energy() );
				new_lrots->set_rotamers( new_rotamers );
				//std::cerr << "new rotlibrary has " << new_lrots->get_rotamers().size() << " members." << std::endl;
			}
		}
	}

	core::conformation::Residue new_res( pose.residue(res_pos).residue_type_set().name_map(res_type_mod_name), true);

	//Temporarily make a copy of the old residue:
	core::conformation::Residue old_res = (*pose.residue(res_pos).clone());

	//replacing the residue
	if(pose.residue(res_pos).is_metal()) {
		utility::vector1< std::pair< std::string, std::string > > atom_pairs;
		for(core::Size ia=1; ia<=3; ia++) atom_pairs.push_back(std::pair<std::string, std::string>(pose.residue(res_pos).atom_name(ia),new_res.atom_name(ia)) );
		pose.replace_residue( res_pos, new_res, atom_pairs); //If this is a metal, replace using only the first three atoms for alignment.
	} else {
		pose.replace_residue( res_pos, new_res, true);
	}

	//and resetting the xyz positions
	for( core::Size at_ct = 1, at_ctmax=pose.residue(res_pos).natoms(); at_ct <= at_ctmax; at_ct++) {
		if(old_res.has( pose.residue(res_pos).atom_name(at_ct) ) ) {
			pose.set_xyz( id::AtomID(at_ct, res_pos), old_res.xyz( pose.residue(res_pos).atom_name(at_ct) ) );
		}
	}

} //add_covalent_linkage_helper


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function that generates a list of metal-binding atoms that coordinate a metal in a protein.
/// @details This function generates the list by looping through all residues and checking all metal-binding atoms of all
/// metal-binding residues, so it's not super speedy.
/// Inputs:
///		pose (The pose that we'll operate on, unchanged by operation)
///		metal_postion (The residue number of the metal)
/// 	dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 < core::id::AtomID >
find_metalbinding_atoms (
	core::pose::Pose const &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier
) {
	if(!pose.residue(metal_position).is_metal()) {
		std::string message = "Error!  Asked to find metal-binding atoms coordinating something that isn't a metal.";
		utility_exit_with_message(message);
	}

	numeric::xyzVector < core::Real > const &metalatomxyz = pose.residue(metal_position).xyz(1); //Atom 1 should be the metal in any metal residue; this gets its xyz position.

	core::Size const nres = pose.n_residue(); //Number of residues in the pose

	if((metal_position < 1) || (metal_position > nres)) {
		std::string message = "Error!  Asked to find metal-binding atoms coordinating a residue that's not in the pose (metal_position < 1 or > n_residue.";
		utility_exit_with_message(message);
	}

	utility::vector1 < core::id::AtomID > outputlist; //The list of AtomIDs that we will output.

	for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues
		if (ir==metal_position) continue; //The metal doesn't bind itself.
		if (pose.residue(ir).is_metalbinding()) { //If this is a metal-binding residue
			utility::vector1 < core::Size > bindingatomlist;
			pose.residue(ir).get_metal_binding_atoms( bindingatomlist );
			if(bindingatomlist.size()>0) {
				for(core::Size ia=1, iamax=bindingatomlist.size(); ia<=iamax; ia++) { //Loop through all the metal-binding atoms in the residue.
					numeric::xyzVector < core::Real > bindingatomxyz = pose.residue(ir).xyz(bindingatomlist[ia]);
					core::Real distsq = metalatomxyz.distance_squared(bindingatomxyz);
					core::Real distcutoffsq = pow(dist_cutoff_multiplier*(pose.residue(metal_position).atom_type(1).lj_radius()+pose.residue(ir).atom_type(bindingatomlist[ia]).lj_radius()), 2);
					if ( distsq <= distcutoffsq ) {
						TR.Debug << "Residue " << ir << " atom " << pose.residue(ir).atom_name( bindingatomlist[ia] ) << " binds the residue " << metal_position << " metal." << std::endl;
						core::id::AtomID curatom( bindingatomlist[ia], ir);
						outputlist.push_back(curatom);
					} //if below distance cutoff
				}//for loop through all metal-binding atoms
			}//if bindingatomlist.size()>0
		}//if is_metalbinding
	}//for loop through all residues

	return outputlist;

} //find_metalbinding_atoms
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to add covalent linkages between a metal atom and all the liganding atoms provided in a vector of AtomIDs.
/// @details The inputs are:
///    pose (The pose to be modified)
///    metal_position (The residue number of the metal in the pose)
///    liganding_atomids (A list of AtomIDs on other residues that will be covalently linked to the metal)
///    remove_hydrogens (Should hydrogens on the liganding atoms be removed automatically?  Default true.)
/// This function uses core::pose::add_covalent_linkage, which can strip off extraneous hydrogens.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkages_to_metal (
	core::pose::Pose &pose,
	core::Size const metal_position,
	utility::vector1 < core::id::AtomID > &liganding_atomids,
	bool const remove_hydrogens
) {
	if(!pose.residue(metal_position).is_metal()) {
		std::string message = "Error!  Asked to add covalent bonds between metal-binding atoms and something that isn't a metal.";
		utility_exit_with_message(message);
	}

	for(core::Size ir=1; ir<=liganding_atomids.size(); ++ir) {
			add_covalent_linkage( pose, liganding_atomids[ir].rsd(), metal_position, liganding_atomids[ir].atomno(), 1, remove_hydrogens);
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to auto-detect and add covalent connections to metal ions.
/// @details This function calls find_metalbinding_atoms and then passes the list of metal binding atoms to add_covalent_linkages_to_metal.
/// Inputs:
///		pose (The pose that we'll operate on, unchanged by operation)
///		metal_postion (The residue number of the metal)
/// 	dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
///   remove_hydrogesn (Should hydrogens on the liganding atoms be auto-removed?  Default true.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_metal_bonds (
	core::pose::Pose &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier,
	bool const remove_hydrogens
) {
	utility::vector1 < core::id::AtomID > metalbinding_atomids = find_metalbinding_atoms(pose, metal_position, dist_cutoff_multiplier);
	add_covalent_linkages_to_metal (pose, metal_position, metalbinding_atomids, remove_hydrogens);

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to auto-detect and add covalent connections to ALL metal ions in a pose.
/// @details This function iteratively calls auto_setup_metal_bonds.
/// Inputs:
///		pose (The pose that we'll operate on, unchanged by operation)
/// 	dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
///   remove_hydrogesn (Should hydrogens on the liganding atoms be auto-removed?  Default true.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_bonds (
	core::pose::Pose &pose,
	core::Real const dist_cutoff_multiplier,
	bool const remove_hydrogens
) {

	TR << "Automatically setting covalent bonds between metal ions and metal-binding residues." << std::endl ;

	core::Size nres = pose.n_residue(); //Residue count
	if(nres > 0) {
		for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues.
			if(pose.residue(ir).is_metal()) auto_setup_metal_bonds (pose, ir, dist_cutoff_multiplier, remove_hydrogens); //If this is a metal, detect bonding atoms and add covalent connections to them.
		}
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// does not set the weights for the constraints terms in the scorefunction.
/// Inputs:
///		pose (The pose that we'll operate on, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints (
	core::pose::Pose &pose,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
	//core::scoring::ScoreFunctionOP sfxn
) {
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	using namespace core::id;

	TR << "Automatically setting up constraints between metal ions and metal-binding residues." << std::endl ;

	core::Size const nres = pose.n_residue();

	for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all residues
		if(pose.residue(ir).is_metal()) { //When we find a metal, check what it's bound to and set up distance constraints accordingly.
			core::Size const n_resconn = pose.residue(ir).n_residue_connections(); //Number of residue connections

			//Make a vector of indices of virtual atoms in this residue:
			utility::vector1 < core::Size > virtindices;
			for(core::Size ia=1, ia_max=pose.residue(ir).natoms(); ia<=ia_max; ++ia) {
				if(pose.residue(ir).is_virtual(ia)) virtindices.push_back(ia);
			}

			for (core::Size jr=1; jr<=n_resconn; ++jr) { //Loop through all of the metal's residue connections
				if(jr > virtindices.size()) { //We'll use virtual atoms for the constraints.  This means that if we have more connections than virtual residues, we can't do this.
					std::string message = "Error!  The number of connections to a metal is greater than the number of virtual atoms in that metal's residuetype.  Unable to continue.";
					utility_exit_with_message(message);
				}

				core::chemical::ResConnID const & curconnection = pose.residue(ir).connect_map(jr); //The current residue connection
				core::Size const otherres = curconnection.resid(); //The other residue to which this one is connected
				core::Size const otherres_conid = curconnection.connid(); //The other residue's connection id for the bond to the metal.
				core::Size const otherres_atom = pose.residue(otherres).residue_connection(otherres_conid).atomno(); //The atom index of the atom to which this metal is bonded.
				core::Size const otherres_atom_parent = pose.residue(otherres).type().atom_base(otherres_atom); //The atom index of the parent of the atom to which this metal is bonded in the other residue.

				AtomID metalID(1, ir);
				AtomID virtID(virtindices[jr], ir);
				AtomID otherID(otherres_atom, otherres);
				AtomID otherparentID(otherres_atom_parent, otherres);

				pose.set_xyz(virtID, pose.xyz(otherID)); //Move this residue's virt to the other residue's metal-binding atom position.

				//Setting up distance constraints:
				if(distance_constraint_multiplier > 1.0e-10) {
					HarmonicFuncOP hfunc = new HarmonicFunc( 0.0, 0.1 / sqrt(distance_constraint_multiplier)); //Harmonic function for constraining the position.
					AtomPairConstraintOP pairconst = new AtomPairConstraint(virtID, otherID, hfunc); //Atom pair constraint holding the virt at the position of the metal-binding atom.
					pose.add_constraint(pairconst); //Add the constraint to the pose, and carry on.
				}

				//Setting up angle constraints:
				if(angle_constraint_multiplier > 1.0e-10) {
					core::Real const ang1 = numeric::angle_radians( pose.residue(ir).xyz(1), pose.residue(otherres).xyz(otherres_atom),  pose.residue(otherres).xyz(otherres_atom_parent) ); //Angle between metal-bonding atom-bonding atom's parent.
					CircularHarmonicFuncOP circfunc1 = new CircularHarmonicFunc( ang1, 0.05/sqrt(angle_constraint_multiplier) ); //Circular harmonic function for constraining angles (works in RADIANS).
					AngleConstraintOP angleconst1 = new AngleConstraint( metalID, otherID, otherparentID, circfunc1 ); //Angle constraint holding the metal
					pose.add_constraint(angleconst1);
				}

			} //Loop through all of the metal's connections
		} //If this is a metal
	} //Loop through all residues

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// sets the weights for the constraints terms in the scorefunction to 1.0 if they're off, or scales the
/// constraints themselves appropriately if they're already on.
/// Inputs:
///		pose (The pose that we'll operate on, changed by operation)
///   sfxn (An owning pointer to the scorefunction, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints (
	core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
) {
	using namespace core::scoring;

	core::Real distmult = 1.0; //Constraint weight multiplier, for cases in which the atom_pair_constraint has already been turned on but is not equal to 1.0.
	if(sfxn->get_weight(atom_pair_constraint) < 1.0e-10) { //If this score term is turned off, turn it on.
		sfxn->set_weight(atom_pair_constraint, 1.0);
	} else { //Otherwise, if this score term is already turned on, set the weight multiplier appropriately.
		distmult = 1.0 / sfxn->get_weight(atom_pair_constraint); //Note -- we will take the square root of this later on.
	}

	core::Real angmult=1.0; //Constraint weight multiplier, for cases in which the atom_pair_constraint has already been turned on but is not equal to 1.0.
	if(sfxn->get_weight(angle_constraint) < 1.0e-10) { //If this score term is turned off, turn it on.
		sfxn->set_weight(angle_constraint, 1.0);
	} else { //Otherwise, if this score term is already turned on, set the weight multiplier appropriately.
		angmult = 1.0 / sfxn->get_weight(angle_constraint); //Note -- we will take the square root of this later on.
	}

	auto_setup_all_metal_constraints ( pose, distmult*distance_constraint_multiplier, angmult*angle_constraint_multiplier );

	return;
}

} // util
} // core
