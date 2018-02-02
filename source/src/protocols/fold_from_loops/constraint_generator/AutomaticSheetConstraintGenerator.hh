// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/AutomaticSheetConstraintGenerator.hh
/// @brief Finds angle constraints between the paired residues in betas
/// @author Jaume Bonet (jaume.bonet@gmail.com)


#ifndef INCLUDED_protocols_fold_from_loops_constraint_generator_AutomaticSheetConstraintGenerator_hh
#define INCLUDED_protocols_fold_from_loops_constraint_generator_AutomaticSheetConstraintGenerator_hh

// Unit headers
#include <protocols/fold_from_loops/constraint_generator/AutomaticSheetConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>
#include <core/scoring/func/Func.fwd.hh>

// Utility headers
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace fold_from_loops {
namespace constraint_generator {

///@brief Generates atom pair constraints for a set of residues from the current or reference pose
class AutomaticSheetConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {

public:
	AutomaticSheetConstraintGenerator();

	~AutomaticSheetConstraintGenerator() override;

	static std::string
	class_name() { return "AutomaticSheetConstraintGenerator"; }

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

public:

	void sd( core::Real const sd ) { sd_ = sd;};
	core::Real sd() {return sd_;};
	void distance( core::Real const dist ){distance_ = dist;};
	core::Real distance(){ return distance_;};


	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	inline static std::string add_sub_ct_name( std::string tag_name ) { return "contraint_generator_internal_" + tag_name + "_complex_type"; }
	static core::Real default_sd() {return 2.0;};
	static core::Real default_distance() {return 6.1;};

	core::scoring::func::FuncOP create_bb_angle_func( core::Real const ideal_angle ) const;
	core::scoring::func::FuncOP create_bb_dihedral_func( core::Real const ideal_dihedral ) const;
	core::scoring::func::FuncOP weighted_func( core::scoring::func::FuncOP func ) const;

private:
	core::Real sd_;
	core::Real distance_;
	core::Real angle_tolerance_;
	core::Real bb_dihedral_tolerance_;
	core::Real weight_;

};

}
} //protocols
} //fold_from_loops

#endif //INCLUDED_protocols_fold_from_loops_AutomaticSheetConstraintGenerator_hh
