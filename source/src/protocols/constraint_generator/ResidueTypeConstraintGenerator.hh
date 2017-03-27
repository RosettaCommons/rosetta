// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ResidueTypeConstraintGenerator.hh
/// @brief Generates residue type constraints for a set of residues from the current or reference pose
/// @author Sharon Guffy (guffy@email.unc.edu)


#ifndef INCLUDED_protocols_constraint_generator_ResidueTypeConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_ResidueTypeConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/ResidueTypeConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
namespace protocols {
namespace constraint_generator {

///@brief Generates atom pair constraints for a set of residues in the current pose
class ResidueTypeConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {

public:
	ResidueTypeConstraintGenerator();

	~ResidueTypeConstraintGenerator() override;

	static std::string
	class_name() { return "ResidueTypeConstraintGenerator"; }

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

public:
	//Setters
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	void
	set_favor_native_bonus( core::Real bonus );

	void
	set_rsd_type_name3( std::string name3 );

	void
	set_reference_pose( core::pose::PoseCOP );

	//Getters (for testing purposes)
	core::select::residue_selector::ResidueSelectorCOP
	get_residue_selector() const;

	core::Real
	get_favor_native_bonus() const;

	std::string
	get_rsd_type_name3() const;

	core::pose::PoseCOP
	get_reference_pose() const;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::scoring::constraints::ConstraintCOPs
	generate_residue_type_constraints( core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset ) const;

	core::scoring::constraints::ConstraintCOPs
	generate_residue_type_constraints(
		core::pose::Pose const & pose,
		core::pose::Pose const & ref_pose,
		core::select::residue_selector::ResidueSubset const & subset,
		core::id::SequenceMapping const & seqmap ) const;


private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	core::Real favor_native_bonus_ = 1.0;
	std::string rsd_type_name3_ = "";
	core::pose::PoseCOP ref_pose_;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_ResidueTypeConstraintGenerator_hh
