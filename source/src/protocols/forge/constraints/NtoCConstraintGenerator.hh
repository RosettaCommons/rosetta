// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/constraints/NtoCConstraintGenerator.hh
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012

#ifndef INCLUDED_protocols_forge_constraints_NtoCConstraintGenerator_hh
#define INCLUDED_protocols_forge_constraints_NtoCConstraintGenerator_hh

// Unit Header
#include <protocols/forge/constraints/NtoCConstraintGenerator.fwd.hh>

// Package Header
#include <protocols/constraint_generator/TerminiConstraintGenerator.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// Proeject Header
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace forge {
namespace constraints {

class NtoCConstraintGenerator : public protocols::forge::remodel::RemodelConstraintGenerator{
public:

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

public:

	NtoCConstraintGenerator();

	NtoCConstraintGenerator( Real const dist, Real const coef );

	virtual ~NtoCConstraintGenerator();

	void
	parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;


	protocols::moves::MoverOP
	fresh_instance() const override;

	protocols::moves::MoverOP
	clone() const override;

	void
	generate_remodel_constraints( Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/*
	void set_weight( Real const coef );

	void set_distance( Real const dist );
	*/

private:
	protocols::constraint_generator::TerminiConstraintGenerator cg_;

}; //class NtoCConstraintGenerator


} //namespace constraints
} //namespace forge
} //namespace protocols


#endif // INCLUDED_protocols_forge_constraints_NtoCConstraintGenerator_HH
