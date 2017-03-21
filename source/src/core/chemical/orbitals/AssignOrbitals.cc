// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/orbitals/AssignOrbitals.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/NumericTraits.hh>

#include <basic/Tracer.hh>

namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS


namespace core {
namespace chemical {
namespace orbitals {

static THREAD_LOCAL basic::Tracer TR( "core::chemical::orbitals::AssignOrbitals" );

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

inline
std::string
strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to do this?
	return trimmed_name;
}
AssignOrbitals::AssignOrbitals(core::chemical::ResidueTypeOP const restype) :
	restype_(restype),
	n_orbitals_(0)

{

}

void AssignOrbitals::assign_orbitals( )
{
	core::chemical::ChemicalManager* chemical_manager = core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCOP atom_type_set = chemical_manager->atom_type_set("fa_standard");
	restype_->clear_orbitals(); // just in case.
	if ( restype_->aa() == aa_tyr || restype_->aa() == aa_phe || restype_->aa() == aa_trp ) {
		if ( restype_->actcoord_atoms().size() != 0 ) {

			utility::vector1<Size> act_atoms =restype_->actcoord_atoms();
			Aindex_ = restype_->actcoord_atoms()[2];
			core::chemical::AtomType const & atmtype(restype_->atom_type(Aindex_));
			if ( !atmtype.is_virtual() ) {


				core::Size atm_index2 =  restype_->bonded_neighbor(Aindex_)[1];
				core::Size atm_index3= restype_->bonded_neighbor(Aindex_)[2];
				AOdist_ =0.7;
				AOhybridization_=2;
				numeric::xyzVector<core::Real> acct_coord = ((restype_->atom(Aindex_).ideal_xyz() - restype_->atom(atm_index2).ideal_xyz())/2);

				//std::cout << restype_->name3() << " " << restype_->atom_name(Aindex_) << " " << restype_->atom_name(atm_index2) << " " << restype_->atom_name(atm_index3) << std::endl;
				//std::cout << acct_coord.x() << " " << acct_coord.y() << " " << acct_coord.z() << std::endl;

				numeric::xyzVector<core::Real> new_action;
				new_action.zero();
				for ( Size ii = 1; ii <= restype_->actcoord_atoms().size(); ++ii ) {
					new_action += restype_->atom(restype_->actcoord_atoms()[ii]).ideal_xyz();
				}
				new_action.x() /= restype_->actcoord_atoms().size();
				new_action.y() /= restype_->actcoord_atoms().size();
				new_action.z() /= restype_->actcoord_atoms().size();


				numeric::xyzVector<core::Real> vector_d( new_action - restype_->atom(atm_index2).ideal_xyz());
				numeric::xyzVector<core::Real> vector_f( new_action - restype_->atom(atm_index3).ideal_xyz());

				//Create an object of Class utility::vector1 to hold the xyz coordinates of orbitals(e.g., cross products)
				//Get two cross products of the two vectors, one is above, the other is below the plane defined by the two vectors
				utility::vector1< numeric::xyzVector<core::Real> > pi_orbital_xyz_vector;
				numeric::xyzVector<core::Real> xyz_right = cross_product(vector_d, vector_f);
				numeric::xyzVector<core::Real> xyz_left = cross_product(-vector_d, vector_f);
				pi_orbital_xyz_vector.push_back((xyz_left.normalized() *  AOdist_) + restype_->atom(atm_index3).ideal_xyz());
				pi_orbital_xyz_vector.push_back((xyz_right.normalized() *  AOdist_) + restype_->atom(atm_index3).ideal_xyz());


				for ( core::Size vector_index = 1; vector_index <= pi_orbital_xyz_vector.size(); ++vector_index ) {
					std::string p_orbital_type_full_name(make_orbital_type_name(atmtype, "pi", AOhybridization_) );
					std::string p_orbital_element_name( make_orbital_element_name() );
					set_orbital_type_and_bond(Aindex_, p_orbital_element_name, p_orbital_type_full_name);

					Vector const stub1_xyz = restype_->atom(Aindex_).ideal_xyz();
					Vector const stub2_xyz = restype_->atom(atm_index2).ideal_xyz();
					Vector const stub3_xyz = restype_->atom(atm_index3).ideal_xyz();

					core::Real const distance(pi_orbital_xyz_vector[vector_index].distance(stub1_xyz) );

					core::Real theta(0.0);
					core::Real phi(0.0);

					if ( distance <1e-2 ) {
						TR.Warning << "extremely small distance=" << distance << " for " <<
							p_orbital_element_name << " ,using 0.0 for theta and phi."<<
							" If you were not expecting this warning, something is very wrong" <<std::endl;
					} else {
						theta =  numeric::angle_radians<core::Real>(pi_orbital_xyz_vector[vector_index],stub1_xyz,stub2_xyz);
						if ( (theta < 1e-2) || (theta > numeric::NumericTraits<Real>::pi()-1e-2) ) {
							phi = 0.0;
						} else {
							phi = numeric::dihedral_radians<core::Real>(pi_orbital_xyz_vector[vector_index],stub1_xyz,stub2_xyz,stub3_xyz);
						}

					}

					std::string const stub1(strip_whitespace(restype_->atom_name(Aindex_)));
					std::string const stub2(strip_whitespace(restype_->atom_name(atm_index2)));
					std::string const stub3(strip_whitespace(restype_->atom_name(atm_index3)));
					core::Real const const_theta(theta);
					core::Real const const_phi(phi);
					std::string const const_name(p_orbital_element_name);
					restype_->set_orbital_icoor_id( const_name,const_phi,const_theta,distance,stub1,stub2,stub3);

				}
			}
		}

	}
	// Get the chemical atom_type for each atom by it index number in this residue
	for ( core::Size atm_index = 1; atm_index <= restype_->natoms(); ++atm_index ) {
		/*OrbInfo orbital_info;*/
		core::chemical::AtomType const & atmtype(restype_->atom_type(atm_index));


		// determine if atom has an orbital, if yes,
		if ( atmtype.atom_has_orbital() ) {
			//get hybridization state of atom_type and orbital type associated with atom type
			core::Size atom_type_index = atom_type_set->atom_type_index(atmtype.atom_type_name());
			core::Size parameter_hybridization = atom_type_set->extra_parameter_index("ORBITAL_HYBRIDIZATION");
			core::Size parameter_orbitaltypes = atom_type_set->extra_parameter_index("ORBITAL_TYPES");
			core::Size parameter_bohrradius = atom_type_set->extra_parameter_index("BOHR_RADIUS");

			Aindex_ = atm_index;
			AOhybridization_ = (core::Size)( atom_type_set->operator[](atom_type_index ).extra_parameter(parameter_hybridization) );
			Orbtype_ = (core::Size)( atom_type_set->operator[](atom_type_index).extra_parameter(parameter_orbitaltypes) );
			AOdist_ = atom_type_set->operator[](atom_type_index).extra_parameter(parameter_bohrradius);
			AObondedatoms_ = restype_->bonded_neighbor(atm_index);

			//very crappy hack. This whole code sucks. WTF? Needs to be rewritten!
			// This will skip C#N among other things, it will do it silently, this should bother you
			if ( atmtype.name() == "Nhis" && AObondedatoms_.size() < 2 ) {
				TR.Warning << "residue " << restype_->name() << " has an Nhis typed atom with < 2 bonds.  It is probably a C#N or something similar for which Rosetta has no reasonable atomtype.  This atom is not being assigned orbitals." <<std::endl;
			} else if ( atmtype.name() == "Nhis" && AObondedatoms_.size() > 1 ) {
				core::Size atm_index2(AObondedatoms_[2]);//atom index of the only bonded neighbor.
				//get the atom indices of the bonded neighbors of atm_index2.
				utility::vector1<core::Size> neighbor_bonded_atms2(restype_->bonded_neighbor(atm_index2));
				core::Size atm_index3(AObondedatoms_[1]);
				/*
				for(core::Size x=1; x<= neighbor_bonded_atms2.size(); ++x){
				//if atm_index is not equal to the index in neighbor bonded atoms to C2,do nothing,continue the loop.
				//crappy hack and causes the neighbor_bonded_atms2 to be somewhat random according to index
				//makes for hard debugging
				if(Aindex_ != neighbor_bonded_atms2[x])
				{
				atm_index3 = neighbor_bonded_atms2[x];
				}
				}
				*/


				numeric::xyzVector<core::Real> vector_a(restype_->atom(Aindex_).ideal_xyz() - restype_->atom(atm_index2).ideal_xyz()   );
				numeric::xyzVector<core::Real> vector_b( restype_->atom(Aindex_).ideal_xyz() - restype_->atom(atm_index3).ideal_xyz() );
				numeric::xyzVector<core::Real> vector_ab_norm = vector_a.normalized()+vector_b.normalized();
				//If this if statment evaluates false, Aindex_ represents the central atom in an azide.
				//Currently, we don't assign orbitals, but we can fix this
				if ( vector_ab_norm != 0.0 ) {
					utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vectors;
					orbital_xyz_vectors.push_back((vector_ab_norm.normalized()*AOdist_)+restype_->atom(Aindex_).ideal_xyz());

					//utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vectors = cross_product_helper(Aindex_,atm_index2,atm_index3,AOdist_);
					//add_orbitals_to_restype(atm_index2, atm_index3, /*orbital_info,*/ atmtype, "p", orbital_xyz_vectors);

					core::Real const phi(numeric::conversions::radians(180.0));
					core::Real const theta(numeric::conversions::radians(54.365));

					//core::Real const phi_De(numeric::conversions::radians(180.0));
					//core::Real const theta_De(numeric::conversions::radians(70.0));

					std::string orbital_type_full_name(make_orbital_type_name(atmtype, "p",AOhybridization_) );
					std::string const orbital_element_name( make_orbital_element_name() );

					set_orbital_type_and_bond(Aindex_, orbital_element_name, orbital_type_full_name);
					//set_orbital_type_and_bond(atm_index, orbital_element_name3, orbital_type_full_name);


					std::string const stub1(strip_whitespace(restype_->atom_name(Aindex_)));
					std::string const stub2(strip_whitespace(restype_->atom_name(atm_index2)));
					std::string const stub3(strip_whitespace(restype_->atom_name(atm_index3)));


					restype_->set_orbital_icoor_id( orbital_element_name,phi,theta,AOdist_,stub1,stub2,stub3);
				} else {
					TR.Warning << "residue " << restype_->name() << " has an Nhis typed atom that is the central N in an azide. It is not being assigned orbitals." <<std::endl;
				}


			}


			//std::cout << Aindex_ << " " << AOhybridization_  <<" " << Orbtype_ << " " << AOdist_ << std::endl;
			if ( AObondedatoms_.size()== 1 && AOhybridization_ == 2 ) {
				if ( Orbtype_ == 1 ) {
					// assign_only_pi_orbitals_to_atom(orbital_info, atmtype);
				} else if ( Orbtype_ == 4 ) {
					assign_sp2_sp_orbitals_to_one_bonded_atom(/*orbital_info,*/ atmtype);
				} else { //if orbitaltypes == 2
					assign_sp2_orbitals_to_one_bonded_atom(/*orbital_info,*/ atmtype);
					//utility_exit_with_message("Both P and Pi orbitals exist if # bonded atoms=1 and hybridization = 2 in case of C=0");
				}
			}
			if ( AObondedatoms_.size()== 2 && AOhybridization_ == 2 ) {
				//Consider this situation C=N-H, the C-N,N-H sigma bonds and the nitrogen lone-pair hybrid orbital are coplanar.
				//the un-hybridazed p orbital at a right anlge to the plane forming a pi bond.
				//utility_exit_with_message("Both P and Pi orbitals exist if # bonded atoms=2 and hybridization = 2, as seen in >C=N-");
				if ( Orbtype_ == 1 ) {
					assign_only_pi_orbitals_to_atom(/*orbital_info,*/ atmtype);
				}
				if ( Orbtype_ == 4 ) {
					//Currently, the atm_index is N1 which is sp2 hybridized and bonded to 2 atoms.N1 has both P and Pi orbitals.
					//Now we need to determine neighbor atoms in order to define a plane to place the orbitals.
					core::Size atm_index2(AObondedatoms_[1]);
					core::Size atm_index3(AObondedatoms_[2]);

					//We first want to add pi orbitals to the Nitrogen,
					assign_only_pi_orbitals_to_atom(/*orbital_info,*/ atmtype);

					//Place one lone pair of P orbitals on sp2 hybridized Nhis atom, which are bonded to two atoms.
					std::string const stub1(strip_whitespace(restype_->atom_name(atm_index)));
					std::string const stub2(strip_whitespace(restype_->atom_name(atm_index2)));
					std::string const stub3(strip_whitespace(restype_->atom_name(atm_index3)));

					core::Real const phi(numeric::conversions::radians(180.0));
					core::Real const theta(numeric::conversions::radians(54.0));

					std::string orbital_type_full_name(make_orbital_type_name(atmtype, "p", AOhybridization_) );
					std::string const orbital_element_name( make_orbital_element_name() );
					set_orbital_type_and_bond(atm_index, orbital_element_name, orbital_type_full_name);
					restype_->set_orbital_icoor_id( orbital_element_name,phi,theta,AOdist_,stub1,stub2,stub3);
				}
			}
			if ( AObondedatoms_.size()== 3 && AOhybridization_== 2 ) { //as seen in COO,NH2O and aroC.
				if ( Orbtype_ == 1 || Orbtype_ == 3 ) { //Assign pi orbitals
					core::Size atm_index2(AObondedatoms_[1]);
					core::Size atm_index3(AObondedatoms_[2]);
					utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vectors = cross_product_helper(atm_index,atm_index2,atm_index3,AOdist_);
					add_orbitals_to_restype(atm_index2, atm_index3, /*orbital_info,*/ atmtype, "pi",orbital_xyz_vectors);
				}
				/*    if (Orbtype_ == 2){
				utility_exit_with_message("P orbital does not exist if # bonded atoms=3 and hybridization = 2, as seen in >C=C<");
				}*/
			}
			if ( AObondedatoms_.size()== 1 && AOhybridization_== 3 ) {
				//this instance exists in PO4-, which has one P=O double bond and then three P-O single bonds.
				// The oxygens that are singly bonded have 3 sets of lone pairs and each has a -1 charge.

			}
			//sp3 hybrdization
			if ( AObondedatoms_.size()== 2 && AOhybridization_== 3 ) {
				/*    if(Orbtype_ != 2){
				utility_exit_with_message("Pi orbital does not exist if # bonded atoms=3 and hybridization = 3");
				}*/

				if ( Orbtype_ == 3 ) { //assign P orbitals to -O-, or -S-, such as -OH
					core::Size atm_index2(AObondedatoms_[1]);
					core::Size atm_index3(AObondedatoms_[2]);

					std::string const stub1(strip_whitespace(restype_->atom_name(atm_index)));
					std::string const stub2(strip_whitespace(restype_->atom_name(atm_index2)));
					std::string const stub3(strip_whitespace(restype_->atom_name(atm_index3)));

					core::Real const phi(numeric::conversions::radians(120.0));
					core::Real const theta(numeric::conversions::radians(70.0));

					//core::Real const phi_De(numeric::conversions::radians(180.0));
					//core::Real const theta_De(numeric::conversions::radians(70.0));

					std::string orbital_type_full_name(make_orbital_type_name(atmtype, "p",AOhybridization_) );
					std::string const orbital_element_name( make_orbital_element_name() );
					std::string const orbital_element_name2( make_orbital_element_name() );
					// amw cppcheck:orbital_element_name3 is not used, but commenting out its unuse is wrong!
					// you see, make_orbital_element_name() has side effects...
					make_orbital_element_name();
					//std::string const orbital_element_name3( make_orbital_element_name() );

					set_orbital_type_and_bond(atm_index, orbital_element_name, orbital_type_full_name);
					set_orbital_type_and_bond(atm_index, orbital_element_name2, orbital_type_full_name);
					//set_orbital_type_and_bond(atm_index, orbital_element_name3, orbital_type_full_name);

					restype_->set_orbital_icoor_id( orbital_element_name,phi,theta,AOdist_,stub1,stub2,stub3);
					restype_->set_orbital_icoor_id( orbital_element_name2,-phi,theta,AOdist_,stub1,stub2,stub3);
					//restype_->set_orbital_icoor_id( orbital_element_name3,phi_De,theta_De,AOdist_,stub1,stub2,stub3);

				}
			}
			if ( AObondedatoms_.size()== 3 && AOhybridization_ == 3 ) {
				/*    if(Orbtype_ != 2){
				utility_exit_with_message("Pi orbital does not exist if # bonded atoms=3 and hybridization = 3");
				}*/

				if ( Orbtype_ == 2 ) { // assign one p orbital to a sp3 N bonded to 3 atoms, as seen in -NH-
					core::Size atm_index2(AObondedatoms_[1]);
					core::Size atm_index3(AObondedatoms_[2]);
					core::Size atm_index4(AObondedatoms_[3]);

					utility::vector1< numeric::xyzVector<core::Real> >  orbital_xyz_vector = Coordinates_Tetrahedral_bondedto3atoms_helper(atm_index, atm_index2, atm_index3, atm_index4, AOdist_);
					for ( core::Size vector_index = 1; vector_index <= orbital_xyz_vector.size(); ++vector_index ) {
						//std::cout << "orb_index=" << vector_index << " x=" << orbital_xyz_vector[vector_index].x() << " y="<< orbital_xyz_vector[vector_index].y() << " z="<< orbital_xyz_vector[vector_index].z() << std::endl;
						std::string orbital_type_full_name(make_orbital_type_name(atmtype, "p", AOhybridization_) );
						std::string orbital_element_name( make_orbital_element_name() );
						set_orbital_type_and_bond(atm_index, orbital_element_name, orbital_type_full_name);
						calculate_orbital_icoor(orbital_xyz_vector[vector_index], atm_index, atm_index2,atm_index3, orbital_element_name);
					}
				}
			}

		}


	}
	restype_->finalize();
}


void AssignOrbitals::assign_only_pi_orbitals_to_atom(/*OrbInfo const & orbital_info,*/ core::chemical::AtomType const & atmtype){
	core::Size atm_index2(AObondedatoms_[1]);//atom index of the only bonded neighbor.
	//get the atom indices of the bonded neighbors of atm_index2.
	utility::vector1<core::Size> neighbor_bonded_atms2(restype_->bonded_neighbor(atm_index2));

	core::Size atm_index3(500);
	bool set_atm_index3 = false;

	for ( core::Size x=1; x<= neighbor_bonded_atms2.size(); ++x ) {
		//if atm_index is not equal to the index in neighbor bonded atoms to C2,do nothing,continue the loop.
		//crappy hack and causes the neighbor_bonded_atms2 to be somewhat random according to index
		//makes for hard debugging
		if ( Aindex_ != neighbor_bonded_atms2[x] ) {
			atm_index3 = neighbor_bonded_atms2[x];
			set_atm_index3 = true;
		}
	}
	if ( !set_atm_index3 ) {
		TR.Warning << "Unable to assign orbitals properly for atom of type " <<atmtype.name() << " On residue " << restype_->name() << " Orbital assignment for this atom is skipped" <<std::endl;
		return;
	}

	utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vectors = cross_product_helper(Aindex_,atm_index2,atm_index3,AOdist_);
	add_orbitals_to_restype(atm_index2, atm_index3, /*orbital_info,*/ atmtype, "pi", orbital_xyz_vectors);
}

void AssignOrbitals::assign_sp2_sp_orbitals_to_one_bonded_atom(/*OrbInfo const & orbital_info,*/ core::chemical::AtomType const & atmtype){
	//Consider this situation. >C1-C2=O1 (C single bond C double bond O). O1 needs two Pi and three P orbitals,
	core::Size atm_index2(AObondedatoms_[1]);//atom index of the only bonded neighbor.
	//get the atom indices of the bonded neighbors of atm_index2.
	utility::vector1<core::Size> neighbor_bonded_atms2(restype_->bonded_neighbor(atm_index2));

	core::Size atm_index3(500);

	for ( core::Size x=1; x<= neighbor_bonded_atms2.size(); ++x ) {
		//if atm_index is not equal to the index in neighbor bonded atoms to C2,do nothing,continue the loop, .
		if ( Aindex_ != neighbor_bonded_atms2[x] ) {
			atm_index3 = neighbor_bonded_atms2[x];
		}
	}
	utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vectors = cross_product_helper(Aindex_,atm_index2,atm_index3,AOdist_);
	add_orbitals_to_restype(atm_index2, atm_index3, /*orbital_info,*/ atmtype, "pi", orbital_xyz_vectors);


	//Place two lone pair of P orbitals on sp2 hybridized O atom, which are bonded to one atoms.
	orbital_xyz_vectors = Coordinates_TriganolPlanar_bondedto1atom_helper(Aindex_,atm_index2,atm_index3,AOdist_);
	add_orbitals_to_restype(atm_index2, atm_index3, /*orbital_info,*/ atmtype, "p", orbital_xyz_vectors);
}

void AssignOrbitals::assign_sp2_orbitals_to_one_bonded_atom(/*OrbInfo const & orbital_info,*/ core::chemical::AtomType const & atmtype){
	//Consider this situation. >C1-C2=O1 (C single bond C double bond O). O1 needs two Pi and three P orbitals,
	core::Size atm_index2(AObondedatoms_[1]);//atom index of the only bonded neighbor.
	//get the atom indices of the bonded neighbors of atm_index2.
	utility::vector1<core::Size> neighbor_bonded_atms2(restype_->bonded_neighbor(atm_index2));

	core::Size atm_index3(500);

	for ( core::Size x=1; x<= neighbor_bonded_atms2.size(); ++x ) {
		//if atm_index is not equal to the index in neighbor bonded atoms to C2,do nothing,continue the loop, .
		if ( Aindex_ != neighbor_bonded_atms2[x] ) {
			atm_index3 = neighbor_bonded_atms2[x];
		}
	}
	//Place two lone pair of P orbitals on sp2 hybridized O atom, which are bonded to one atoms.
	utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vectors = Coordinates_TriganolPlanar_bondedto1atom_helper(Aindex_,atm_index2,atm_index3,AOdist_);
	add_orbitals_to_restype(atm_index2, atm_index3, /*orbital_info,*/ atmtype, "p", orbital_xyz_vectors);
}


// To get a pair of pi orbitals, we calculate the cross products of two vectors both pointing towards the atom with an index of atm_index1.
// The two cross products are perpendicular to the plane; one is above and the other is below the plane.
// core::Real dist is the Bohr radius of H plus the Bohr radius of the first atom with atm_index1
utility::vector1< numeric::xyzVector<core::Real> > AssignOrbitals::cross_product_helper
(
	core::Size const atm_index1,
	core::Size const atm_index2,
	core::Size const atm_index3,
	core::Real const dist
)
{

	//std::cout << "atm_index1: " << atm_index1 << " index 2: " << atm_index2 << " index3 " << atm_index3 << std::endl;

	//define two vectors, both pointing back to the central atom with atm_index2
	numeric::xyzVector<core::Real> vector_d( restype_->atom(atm_index1).ideal_xyz() - restype_->atom(atm_index2).ideal_xyz());
	numeric::xyzVector<core::Real> vector_f( restype_->atom(atm_index1).ideal_xyz() - restype_->atom(atm_index3).ideal_xyz());

	//Create an object of Class utility::vector1 to hold the xyz coordinates of orbitals(e.g., cross products)
	//Get two cross products of the two vectors, one is above, the other is below the plane defined by the two vectors
	utility::vector1< numeric::xyzVector<core::Real> > pi_orbital_xyz_vector;
	numeric::xyzVector<core::Real> xyz_right = cross_product(vector_d, vector_f);
	numeric::xyzVector<core::Real> xyz_left = cross_product(-vector_d, vector_f);

	//Normalize the two new vectors, xyz_right and xyz_left to get a unit vector.
	//pi_orbital_xyz_vector now stores the new xyz coordinates of the pi orbitals.
	pi_orbital_xyz_vector.push_back((xyz_right.normalized() * dist) + restype_->atom(atm_index1).ideal_xyz());
	pi_orbital_xyz_vector.push_back((xyz_left.normalized() *  dist) + restype_->atom(atm_index1).ideal_xyz());

	return pi_orbital_xyz_vector;

}


void AssignOrbitals::add_orbitals_to_restype(
	core::Size const atm_index2,
	core::Size const atm_index3,
	//OrbInfo const & orbital_info,
	core::chemical::AtomType const & atmtype,
	std::string const & atom_hybridization,
	utility::vector1< numeric::xyzVector<core::Real> > const & orbital_xyz_vectors
){
	for ( core::Size vector_index = 1; vector_index <= orbital_xyz_vectors.size(); ++vector_index ) {
		std::string p_orbital_type_full_name(make_orbital_type_name(atmtype, atom_hybridization, AOhybridization_) );
		std::string p_orbital_element_name( make_orbital_element_name() );
		set_orbital_type_and_bond(Aindex_, p_orbital_element_name, p_orbital_type_full_name);
		calculate_orbital_icoor(orbital_xyz_vectors[vector_index],Aindex_, atm_index2, atm_index3, p_orbital_element_name);
	}
}

std::string AssignOrbitals::make_orbital_type_name
(
	AtomType const & atmtype,
	std::string const & orbitaltype,
	core::Size const hybridization
)
{
	std::string atm_element( atmtype.element() );

	std::string hyb;
	if ( hybridization == 1 ) {
		hyb = "sp";
	} else if ( hybridization ==2 ) {
		hyb="sp2";
	} else if ( hybridization ==3 ) {
		hyb="sp3";
	}

	std::string orbital_type_full_name(atm_element+"."+orbitaltype+"."+hyb);
	return orbital_type_full_name;
}

std::string AssignOrbitals::make_orbital_element_name()
{
	++n_orbitals_;
	std::string orbital_name("LP");
	std::string orb_index_string = utility::to_string<core::Size>(n_orbitals_) ;
	std::string orbital_element_name(orbital_name+orb_index_string);
	return orbital_element_name;
}


//Assign orbital types and bond information which will be passed to the function restype_->set_orbital_icoor_id to get
//icoord for all orbitals.
void AssignOrbitals::set_orbital_type_and_bond(
	core::Size atom_index,
	std::string orbital_element_name,
	std::string orbital_type_full_name

){
	// Orbital names are given by concatenate two strings:'LP" and the indices of the orbitals on the residue(restype);
	std::string atm_name(strip_whitespace(restype_->atom_name(atom_index)));

	restype_->add_orbital(orbital_element_name, orbital_type_full_name);
	restype_->add_orbital_bond(atm_name, orbital_element_name);

}


void AssignOrbitals::calculate_orbital_icoor(
	numeric::xyzVector<core::Real> const & orbital_xyz,
	core::Size const atm_index1,
	core::Size const atm_index2,
	core::Size const atm_index3,
	std::string const & orbital_element_name
)
{
	Vector const stub1_xyz = restype_->atom(atm_index1).ideal_xyz();
	Vector const stub2_xyz = restype_->atom(atm_index2).ideal_xyz();
	Vector const stub3_xyz = restype_->atom(atm_index3).ideal_xyz();

	core::Real const distance(orbital_xyz.distance(stub1_xyz) );

	core::Real theta(0.0);
	core::Real phi(0.0);

	if ( distance <1e-2 ) {
		TR.Warning << "extremely small distance=" << distance << " for " <<
			orbital_element_name << " ,using 0.0 for theta and phi."<<
			" If you were not expecting this warning, something is very wrong" <<std::endl;
	} else {
		theta =  numeric::angle_radians<core::Real>(orbital_xyz,stub1_xyz,stub2_xyz);
		if ( (theta < 1e-2) || (theta > numeric::NumericTraits<Real>::pi()-1e-2) ) {
			phi = 0.0;
		} else {
			phi = numeric::dihedral_radians<core::Real>(orbital_xyz,stub1_xyz,stub2_xyz,stub3_xyz);
		}

	}

	std::string const stub1(strip_whitespace(restype_->atom_name(atm_index1)));
	std::string const stub2(strip_whitespace(restype_->atom_name(atm_index2)));
	std::string const stub3(strip_whitespace(restype_->atom_name(atm_index3)));
	core::Real const const_theta(theta);
	core::Real const const_phi(phi);
	std::string const const_name(orbital_element_name);


	//restype_->add_orbital( orbital_element_name, );
	//restype_->add_orbital_bond(stub1, orbital_element_name);


	//tr << orbital << " " << stub_atom1 << " "<< stub_atom2 << " " <<stub_atom3 << " " <<distance << " " << phi << " " << theta <<std::endl;
	restype_->set_orbital_icoor_id( const_name,const_phi,const_theta,distance,stub1,stub2,stub3);
}


utility::vector1< numeric::xyzVector<core::Real> >  AssignOrbitals::Coordinates_TriganolPlanar_bondedto1atom_helper(
	core::Size const atm_index1,
	core::Size const atm_index2,
	core::Size const atm_index3,
	core::Real const dist

){

	//create a new vector to hold coordinates of P orbitals.
	utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vector;

	numeric::xyzVector<core::Real> vector_a = restype_->atom(atm_index1).ideal_xyz();
	numeric::xyzVector<core::Real> vector_b = restype_->atom(atm_index2).ideal_xyz();
	numeric::xyzVector<core::Real> vector_c = restype_->atom(atm_index3).ideal_xyz();

	//core::Real distance_xa = 01.0;
	core::Real angle_xab = numeric::constants::r::pi_over_3; //60 degrees
	//for one point it should be 180 for another it should be 0
	//core::Real dihedral_xabc = numeric::constants::r::pi;

	numeric::xyzVector<core::Real> a( (vector_a - vector_b).normalized());
	numeric::xyzVector<core::Real> b((vector_b - vector_c).normalized());
	numeric::xyzVector<core::Real> c(cross_product(a,b).normalized());
	numeric::xyzVector<core::Real> d(cross_product(a,c).normalized());

	core::Real dihedral_xabc1 = numeric::constants::r::pi;
	core::Real dihedral_xabc2 = 0;

	numeric::xyzVector<core::Real> v1 = a * std::cos(numeric::constants::r::pi - angle_xab);
	numeric::xyzVector<core::Real> v2 = c * std::sin(numeric::constants::r::pi - angle_xab);
	numeric::xyzVector<core::Real> v3 = d * std::sin(numeric::constants::r::pi - angle_xab);

	numeric::xyzVector<core::Real> x1
		(
		(
		v1 - v2 * std::sin(dihedral_xabc1) + v3 * std::cos(dihedral_xabc1)
		) * dist
	);

	numeric::xyzVector<core::Real> x2
		(
		(
		v1 - v2 * std::sin(dihedral_xabc2) + v3 * std::cos(dihedral_xabc2)
		) * dist
	);

	/*
	//add the degenerative orbital. Currently commented out as it screws with correct geometry
	numeric::xyzVector<core::Real> x3
	(
	(
	-(vector_a - vector_b).normalized()
	) * dist
	);
	*/

	orbital_xyz_vector.push_back(vector_a+x1);
	orbital_xyz_vector.push_back(vector_a+x2);
	//orbital_xyz_vector.push_back(vector_a+x3);

	return orbital_xyz_vector;
}


// Get the coordinates of one P orbital (x1)attached to a sp3 hybridized atom which is bonded to three atoms (e.g. -NH-).
utility::vector1< numeric::xyzVector<core::Real> >  AssignOrbitals::Coordinates_Tetrahedral_bondedto3atoms_helper(
	core::Size const atm_index1,
	core::Size const atm_index2,
	core::Size const atm_index3,
	core::Size const atm_index4,
	core::Real const dist

){

	utility::vector1< numeric::xyzVector<core::Real> > orbital_xyz_vector;

	numeric::xyzVector<core::Real> vector_a = restype_->atom(atm_index1).ideal_xyz();
	numeric::xyzVector<core::Real> vector_b = restype_->atom(atm_index2).ideal_xyz();
	numeric::xyzVector<core::Real> vector_c = restype_->atom(atm_index3).ideal_xyz();
	numeric::xyzVector<core::Real> vector_d = restype_->atom(atm_index4).ideal_xyz();

	numeric::xyzVector<core::Real> x
		(
		( ( vector_a - vector_b).normalized()
		+ ( vector_a - vector_c).normalized() + ( vector_a - vector_d).normalized()).normalized() * dist
	);


	orbital_xyz_vector.push_back(vector_a+x);

	return orbital_xyz_vector;
}


}//namespace
}//namespace
}//namespace


