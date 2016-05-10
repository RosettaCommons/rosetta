// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraint_generator/HydrogenBondConstraintGenerator.hh
/// @brief
/// @author Tom Linsky ( tlinsky at uw dot edu )


#ifndef INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/HydrogenBondConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// Numeric headers

// C++ headers
#include <map>
#include <set>
#include <list>

namespace protocols {
namespace constraint_generator {

class HydrogenBondingAtom : public utility::pointer::ReferenceCount {

public:
	typedef std::list< core::Real > Dihedrals;

	HydrogenBondingAtom(
		std::string const & atom1,
		std::string const & atom2,
		std::string const & atom3,
		core::Real const ideal_distance,
		core::Real const ideal_angle,
		Dihedrals const & ideal_dihedrals );

	std::string const & hb_atom() const { return atom_; }
	std::string const & atom2() const { return atom2_; }
	std::string const & atom3() const { return atom3_; }

	core::Real distance() const { return distance_; }
	core::Real angle() const { return angle_; }
	Dihedrals const & dihedrals() const { return dihedrals_; }

private:
	HydrogenBondingAtom();

private:
	std::string atom_;
	std::string atom2_;
	std::string atom3_;
	core::Real distance_;
	core::Real angle_;
	Dihedrals dihedrals_;
};

typedef std::list< HydrogenBondingAtom > HydrogenBondingAtoms;

/// @brief Database to lookup atoms for hydrogen bonding
class HydrogenBondInfo : public utility::SingletonBase< HydrogenBondInfo > {

public:
	typedef std::map< std::string, HydrogenBondingAtoms > AtomNameMap;

	HydrogenBondInfo();

	HydrogenBondingAtoms const & atoms( std::string const & rsd_name ) const;

	static HydrogenBondInfo * create_singleton_instance();

	HydrogenBondingAtoms &
	create_residue( std::string const & rsd_name );

private:
	AtomNameMap atoms_;
	static HydrogenBondingAtoms const empty_;
};

/// @brief This constraint generator generates constraints for favoring formation of hydrogen bonds
class HydrogenBondConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	HydrogenBondConstraintGenerator();
	virtual ~HydrogenBondConstraintGenerator();

	virtual protocols::constraint_generator::ConstraintGeneratorOP clone() const;

	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data );

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const;

public:
	void
	set_residue_selector1( core::select::residue_selector::ResidueSelectorCOP selector );

	void
	set_residue_selector2( core::select::residue_selector::ResidueSelectorCOP selector );

	void
	set_atoms1( std::string const & atoms_str );

	void
	set_atoms1( std::set< std::string > const & atoms );

	void
	set_atoms2( std::string const & atoms_str );

	void
	set_atoms2( std::set< std::string > const & atoms );

	void
	set_atom_pair_func( std::string const & func_def );

	void
	set_atom_pair_func( core::scoring::func::FuncOP atompairfunc );

	void
	set_angle1_func( std::string const & func_def );

	void
	set_angle1_func( core::scoring::func::FuncOP angle_func );

	void
	set_angle2_func( std::string const & func_def );

	void
	set_angle2_func( core::scoring::func::FuncOP angle_func );

	void
	set_bounded( bool const bounded );

	void
	set_atom_pair_sd( core::Real const sd );

	void
	set_angle_sd( core::Real const sd );

private:

	core::scoring::constraints::ConstraintOP
	create_residue_constraint(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2 ) const;

	core::scoring::constraints::ConstraintOP
	create_residue_constraint(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		HydrogenBondingAtoms const & atoms1,
		HydrogenBondingAtoms const & atoms2 ) const;

	core::scoring::constraints::ConstraintOP
	create_residue_constraint(
		core::id::AtomID const & atomid1,
		core::id::AtomID const & parent_atomid1,
		core::id::AtomID const & parent2_atomid1,
		core::id::AtomID const & atomid2,
		core::id::AtomID const & parent_atomid2,
		core::id::AtomID const & parent2_atomid2,
		HydrogenBondingAtom const & a1,
		HydrogenBondingAtom const & a2 ) const;

	HydrogenBondingAtoms
	choose_atoms(
		core::conformation::Residue const & rsd,
		std::set< std::string > const & allowed_atoms ) const;

	core::scoring::func::FuncOP
	atom_pair_func( HydrogenBondingAtom const & a1,	HydrogenBondingAtom const & a2 ) const;

	core::scoring::func::FuncOP
	angle1_func( HydrogenBondingAtom const & a1 ) const;

	core::scoring::func::FuncOP
	angle2_func( HydrogenBondingAtom const & a2 ) const;

	core::scoring::func::FuncOP
	dihedral_func( core::Real const dihedral_value ) const;

private:
	core::select::residue_selector::ResidueSelectorCOP selector1_;
	core::select::residue_selector::ResidueSelectorCOP selector2_;
	std::set< std::string > atoms1_;
	std::set< std::string > atoms2_;
	core::scoring::func::FuncOP atom_pair_func_;
	core::scoring::func::FuncOP angle1_func_;
	core::scoring::func::FuncOP angle2_func_;
	bool bounded_;
	core::Real atom_pair_sd_;
	core::Real angle_sd_;

}; // HydrogenBondConstraintGenerator

} // namespace constraint_generator
} // namespace protocols

#endif // INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGenerator_HH
