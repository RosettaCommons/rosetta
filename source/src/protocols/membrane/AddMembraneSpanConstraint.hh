// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneSpanConstraintCreator.hh
/// @brief      add membrant span constraint
/// @author     Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_membrane_AddMembraneSpanConstraint_hh
#define INCLUDED_protocols_membrane_AddMembraneSpanConstraint_hh

// Unit Headers
#include <protocols/membrane/AddMembraneSpanConstraint.fwd.hh>
#include <protocols/membrane/AddMembraneSpanConstraintCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

class AddMembraneSpanConstraint : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	AddMembraneSpanConstraint();

	/// @brief Copy Constructor
	AddMembraneSpanConstraint( AddMembraneSpanConstraint const & src );

	/// @brief Assignment Operator
	AddMembraneSpanConstraint & operator = ( AddMembraneSpanConstraint const & src );

	/// @brief Destructor
	~AddMembraneSpanConstraint() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (AddMembraneSpanConstraint)
	// XRW TEMP  std::string get_name() const override;

	/// @brief Flip the downstream partner in the membrane
	void apply( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneSpanConstraint_hh
