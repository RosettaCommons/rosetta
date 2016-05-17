// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/fldsgn/SheetRemodelConstraintGenerator.hh
/// @brief Remodel constraint generator for adding sheet constraints
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_fldsgn_SheetRemodelConstraintGenerator_hh
#define INCLUDED_protocols_fldsgn_SheetRemodelConstraintGenerator_hh

// Unit header
#include <protocols/fldsgn/SheetRemodelConstraintGenerator.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// Protocol header
#include <protocols/fldsgn/SheetConstraintGenerator.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace fldsgn {

///@brief Remodel constraint generator for adding sheet constraints
class SheetRemodelConstraintGenerator : public protocols::forge::remodel::RemodelConstraintGenerator {
public:
	SheetRemodelConstraintGenerator();
	virtual ~SheetRemodelConstraintGenerator();

	protocols::moves::MoverOP
	clone() const;

	std::string
	get_name() const;

	static std::string
	class_name() { return "SheetCstGenerator"; }

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

protected:
	virtual void
	generate_remodel_constraints( core::pose::Pose const & pose );

private:
	SheetConstraintGeneratorOP cg_;
};


} //protocols
} //fldsgn

#endif //INCLUDED_protocols_fldsgn_SheetRemodelConstraintGenerator_hh
