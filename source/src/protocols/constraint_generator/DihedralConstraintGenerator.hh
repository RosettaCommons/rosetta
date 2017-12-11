// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/DihedralConstraintGenerator.hh
/// @brief A cst generator that creates Dihedral constraints for specified residues using a residue selector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_constraint_generator_DihedralConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_DihedralConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/DihedralConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace constraint_generator {

///@brief A cst generator that creates Dihedral constraints for specified residues using a residue selector.
///@details
/// Uses CircularHarmonic constraints, since CircularGaussian func does not exist.
/// By default, works on Protein and carbohydrate BB dihedrals.
/// You must set either the set of atom names or phi/psi/omega.  Will only work on ONE type of dihedral angle.
///
class DihedralConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	DihedralConstraintGenerator();

	virtual ~DihedralConstraintGenerator();

	static std::string
	class_name();

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


public:

	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	///@brief Set the torsion type we will be working with.
	void
	set_torsion_type( core::id::MainchainTorsionType torsion );

	///@brief MANUALLY set the dihedrals we will work with.
	/// IGNORES Residue selector, but allows setting dihedrals for ANY torsion.
	void
	set_custom_dihedral(utility::vector1< core::id::AtomID > const & ids );

	///@brief Set the standard deviation for the dihedral constraint.  Default is 16.0.
	/// THis which was found by taking the mean SD of all dihedral angles of either
	/// PHI or PSI for each North CDR Cluster.  THis is a fairly tight constraint and allows a bit of movement while not changing overall struture much.
	///
	///@details
	/// Set this in Degrees.
	void
	set_sd_degree( core::Real sd_degree );

protected:

	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

private:

	core::id::MainchainTorsionType torsion_;

	core::Real sd_ = 16.0; //From SD of dihedrals within CDR Clusters.  THis is very conservative.

	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	std::string parsed_resnums_ = "";
	std::string parsed_atoms_ = "";
	bool parsed_custom_torsion_ = false;

	utility::vector1< core::id::AtomID > custom_torsion_;



};

///@brief Turn two strings, comma-separated with atom names and resnums into a custom dihedral vector.
utility::vector1< core::id::AtomID >
parse_custom_torsion( std::string const & atom_names, std::string const & resnums, core::pose::Pose const & pose );

///@brief Get the MainchainTorsionType from a string (phi/psi/omega/omega2/omega3).
core::id::MainchainTorsionType
parse_torsion_type( std::string const & torsion_name);


} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_DihedralConstraintGenerator_hh
