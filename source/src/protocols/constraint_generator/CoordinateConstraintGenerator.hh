// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/CoordinateConstraintGenerator.hh
/// @brief Generates coodinate constraints
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_constraint_generator_CoordinateConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_CoordinateConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/CoordinateConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace constraint_generator {

///@brief Generates coodinate constraints
class CoordinateConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {
public:
	CoordinateConstraintGenerator();
	~CoordinateConstraintGenerator() override;

	static std::string
	class_name() { return "CoordinateConstraintGenerator"; }

	ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

public:
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	void
	set_sidechain( bool const sidechain );

	void
	set_ca_only( bool const ca_only );

	void
	set_ambiguous_hnq( bool const ambiguous_hnq );

	void
	set_sd( core::Real const sd );

	void
	set_bounded( bool const bounded );

	void
	set_bounded_width( core::Real const bounded_width );

	void
	set_reference_pose( core::pose::PoseCOP refpose );

private:
	core::pose::Pose const &
	select_constraint_target_pose( core::pose::Pose const & pose ) const;

	core::id::SequenceMapping
	generate_seqmap( core::pose::Pose const & pose ) const;

	core::pose::PoseOP
	prepare_constraint_target_pose(
		core::pose::Pose const & input,
		core::id::SequenceMapping const & seq_map ) const;

	void
	create_residue_constraints(
		core::scoring::constraints::ConstraintCOPs & csts,
		core::id::AtomID const & root_atomid,
		core::conformation::Residue const & pose_rsd,
		core::conformation::Residue const & ref_rsd ) const;

	core::scoring::func::FuncOP
	create_function() const;

	core::scoring::constraints::ConstraintOP
	create_coordinate_constraint(
		core::id::AtomID const & pose_atomid,
		core::id::AtomID const & ref_atom,
		core::id::AtomID const & root_atom,
		core::conformation::Residue const & ref_rsd ) const;

	core::scoring::constraints::ConstraintOP
	create_ambiguous_constraint(
		core::id::AtomID const & pose_atomid,
		core::conformation::Residue const & ref_rsd,
		core::id::AtomID const & root_atomid,
		std::pair< std::string, std::string > const & hnq_atoms ) const;

	void
	align_reference_pose(
		core::pose::Pose & reference,
		core::pose::Pose const & input_pose,
		core::id::SequenceMapping const & seq_map ) const;

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	core::pose::PoseCOP refpose_;

	bool bounded_;
	core::Real bounded_width_;
	core::Real sd_;
	bool ca_only_;
	bool sidechain_;
	bool ambiguous_hnq_;
	bool align_reference_;
};

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_CoordinateConstraintGenerator_hh

