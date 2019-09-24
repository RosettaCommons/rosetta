// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_core_chemical_orbitals_AssignOrbitals_hh
#define INCLUDED_core_chemical_orbitals_AssignOrbitals_hh

#include <numeric/xyzVector.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <utility/vector1.hh>
#include <core/chemical/ResidueGraphTypes.hh>

//Auto Headers
#include <core/chemical/AtomType.fwd.hh>


namespace core {
namespace chemical {
namespace orbitals {

/*struct OrbInfo{
core::Size atom_index;
core::Size hybridization;
core::Size orbitaltypes;
core::Real dist;
utility::vector1<core::Size> bondedatoms;

};*/

class AssignOrbitals {
public:


	AssignOrbitals(core::chemical::MutableResidueTypeOP const restype);
	void assign_orbitals();

	void assign_only_pi_orbitals_to_atom(/*OrbInfo const & orbital_info,*/ core::chemical::AtomType const & atmtype);

	void assign_sp2_sp_orbitals_to_one_bonded_atom(/*OrbInfo const & orbital_info*/ core::chemical::AtomType const & atmtype);
	void add_orbitals_to_restype(
		VD atm_vd2,
		VD atm_vd3,
		//OrbInfo const & orbital_info,
		core::chemical::AtomType const & atmtype,
		std::string const & atom_hybridization,
		utility::vector1< numeric::xyzVector<core::Real> > const & orbital_xyz_vectors
	);

	void assign_sp2_orbitals_to_one_bonded_atom(/*OrbInfo const & orbital_info,*/ core::chemical::AtomType const & atmtype);


	utility::vector1< numeric::xyzVector<core::Real> > cross_product_helper(
		VD const atm_vd1,
		VD const atm_vd2,
		VD const atm_vd3,
		core::Real const dist
	);

	void calculate_orbital_icoor(
		numeric::xyzVector<core::Real> const & orbital_xyz,
		VD const atm_vd1,
		VD const atm_vd2,
		VD const atm_vd3,
		std::string const & orbital_element_name
	);

	utility::vector1< numeric::xyzVector<core::Real> >  Coordinates_TriganolPlanar_bondedto1atom_helper(
		VD const atm_vd1,
		VD const atm_vd2,
		VD const atm_vd3,
		core::Real const dist

	);

	utility::vector1< numeric::xyzVector<core::Real> >  Coordinates_Tetrahedral_bondedto3atoms_helper(
		VD const atm_vd1,
		VD const atm_vd2,
		VD const atm_vd3,
		VD const atm_vd4,
		core::Real const dist

	);

private:

	core::chemical::VD Avd_;
	core::Size AOhybridization_;
	core::Size Orbtype_;
	core::Real AOdist_;
	utility::vector1<VD> AObondedatoms_;

	core::chemical::MutableResidueTypeOP restype_;
	core::Size n_orbitals_;

	std::string make_orbital_type_name
	(
		AtomType const & atmtype,
		std::string const & orbitaltype,
		core::Size const hybridization
	);

	std::string make_orbital_element_name();


	void set_orbital_type_and_bond(
		VD atom_vd,
		std::string orbital_element_name,
		std::string orbital_type_full_name
	);

};


}
}
}


#endif /* INCLUDED_core_chemical_orbitals_AssignOrbitals_hh */
