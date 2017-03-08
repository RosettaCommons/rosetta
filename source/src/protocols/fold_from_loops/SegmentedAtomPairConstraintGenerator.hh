// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/SegmentedAtomPairConstraintGenerator.hh
/// @brief Given a set of non-continuous selected segments, generates differently scored atom pair constraints
///       for the resides in each segment and between segments.
/// @author Jaume Bonet (jaume.bonet@gmail.com)


#ifndef INCLUDED_protocols_fold_from_loops_SegmentedAtomPairConstraintGenerator_hh
#define INCLUDED_protocols_fold_from_loops_SegmentedAtomPairConstraintGenerator_hh

// Unit headers
#include <protocols/fold_from_loops/SegmentedAtomPairConstraintGenerator.fwd.hh>
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
namespace fold_from_loops {

struct ConstraintConditions {
	core::Real sd;
	core::Real weight;
	bool ca_only;
	bool use_harmonic;
	bool unweighted;
	core::Real max_distance;
	core::Size min_seq_sep;
};

///@brief Generates atom pair constraints for a set of residues from the current or reference pose
class SegmentedAtomPairConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {

public:
	SegmentedAtomPairConstraintGenerator();

	~SegmentedAtomPairConstraintGenerator() override;

	static std::string
	class_name() { return "SegmentedAtomPairConstraintGenerator"; }

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

public:
	void
	set_residue_selector( core::select::residue_selector::ResidueSelector const & selector );

	void
	set_inner_sd( core::Real const sd );

	void
	set_inner_ca_only( bool const ca_only );

	void
	set_inner_use_harmonic_function( bool const use_harmonic );

	void
	set_inner_unweighted_function( bool const unweighted );

	void
	set_inner_min_seq_sep( core::Size const min_seq_sep );

	void
	set_inner_weight( core::Real const weight );

	void
	set_outer_sd( core::Real const sd );

	void
	set_outer_ca_only( bool const ca_only );

	void
	set_outer_use_harmonic_function( bool const use_harmonic );

	void
	set_outer_unweighted_function( bool const unweighted );

	void
	set_outer_max_distance( core::Real const max_dist );

	void
	set_outer_weight( core::Real const weight );

	void
	set_reference_pose( core::pose::PoseCOP ref_pose );

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	inline static std::string add_sub_ct_name( std::string tag_name ) { return "contraint_generator_internal_" + tag_name + "_complex_type"; }

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	core::pose::PoseCOP reference_pose_;
	ConstraintConditions inner_;
	ConstraintConditions outer_;

};

} //protocols
} //fold_from_loops

#endif //INCLUDED_protocols_fold_from_loops_SegmentedAtomPairConstraintGenerator_hh
