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

// Numeric headers

// C++ headers

namespace protocols {
namespace constraint_generator {

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
	set_atoms1( utility::vector1< std::string > const & atoms );

	void
	set_atoms2( std::string const & atoms_str );

	void
	set_atoms2( utility::vector1< std::string > const & atoms );

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

private:

	core::scoring::constraints::ConstraintOP
	create_residue_constraint(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2 ) const;

	core::scoring::constraints::ConstraintOP
	create_residue_constraint(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		utility::vector1< core::Size > const & atoms1,
		utility::vector1< core::Size > const & atoms2 ) const;

	core::scoring::constraints::ConstraintOP
	create_residue_constraint(
		core::id::AtomID const & atomid1,
		core::id::AtomID const & parent_atomid1,
		core::id::AtomID const & atomid2,
		core::id::AtomID const & parent_atomid2 ) const;

	utility::vector1< core::Size >
	choose_atoms(
		core::conformation::Residue const & rsd,
		utility::vector1< std::string > const & atoms ) const;

private:
	core::select::residue_selector::ResidueSelectorCOP selector1_;
	core::select::residue_selector::ResidueSelectorCOP selector2_;
	utility::vector1< std::string > atoms1_;
	utility::vector1< std::string > atoms2_;
	core::scoring::func::FuncOP atom_pair_func_;
	core::scoring::func::FuncOP angle1_func_;
	core::scoring::func::FuncOP angle2_func_;

}; // HydrogenBondConstraintGenerator

} // namespace constraint_generator
} // namespace protocols

#endif // INCLUDED_protocols_constraint_generator_HydrogenBondConstraintGenerator_HH
