// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ConstraintGenerator.hh
///
/// @brief  abstract base class that generates constraints
/// @author Tom Linsky (tlinsky at uw dot edu)


#ifndef INCLUDED_protocols_constraint_generator_ConstraintGenerator_hh
#define INCLUDED_protocols_constraint_generator_ConstraintGenerator_hh

// Unit headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>

// Core headers
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh> // WIN32 INCLUDE
#endif
#include <core/scoring/func/Func.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace constraint_generator {

/// @brief pure virtual base class
class ConstraintGenerator : public utility::pointer::ReferenceCount {

public: // constructors / destructors
	ConstraintGenerator( std::string const & class_name );
	virtual ~ConstraintGenerator();

private:
	ConstraintGenerator();

public: //overridden virtuals

	virtual ConstraintGeneratorOP
	clone() const = 0;

	/// @brief generates constraints and adds them to the pose
	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const = 0;

protected:
	/// @brief called by parse_my_tag -- should not be used directly
	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) = 0;

public:
	/// @brief parses XML tag -- calls protected parse_tag() function
	void
	parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	std::string const &
	id() const;

	void
	set_id( std::string const & id );

	std::string const &
	class_name() const;

public:
	void
	clear_stored_constraints() const;

	void
	store_constraints( core::scoring::constraints::ConstraintCOPs const & csts ) const;

	core::scoring::constraints::ConstraintCOPs const &
	retrieve_constraints() const;

private:
	std::string class_name_;
	std::string id_;

}; //class ConstraintGenerator

} //namespace constraint_generator
} //namespace protocols


#endif
