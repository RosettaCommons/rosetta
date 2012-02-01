// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief
/// @author Javier Castellanos ( javiercv@uw.edu )

#ifndef INCLUDED_constrained_sequence_design_SequenceContraintFactory_hh
#define INCLUDED_constrained_sequence_design_SequenceContraintFactory_hh

// Unit Headers

// Package Headers
#include <devel/constrained_sequence_design/SequenceConstraint.hh> 
#include <devel/constrained_sequence_design/SequenceConstraintCreator.hh> 

// Project Headers
#include <core/pose/Pose.fwd.hh>
// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.fwd.hh>

namespace devel {
namespace constrained_sequence_design {

// singleton class that generates different SequenceConstraint objects.
class SequenceConstraintFactory : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose Pose;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagPtr TagPtr;
	typedef std::map< std::string, SequenceConstraintCreatorOP > SequenceConstraintMap;
	
public:
	static SequenceConstraintFactory * get_instance();
	void factory_register( SequenceConstraintCreatorOP );



	/// @brief return a new SequenceConstraint object of constriant_type class
	SequenceConstraintOP newSequenceConstraint(const std::string & constraint_type);

	/// @brief return a new SequenceConstraint object by tag parsing
	SequenceConstraintOP
		newSequenceConstraint(
		TagPtr const,
		protocols::moves::DataMap &,
		Pose const &
	);
	
private:
	// singleton private constructor
	SequenceConstraintFactory();
	~SequenceConstraintFactory();

	// Unimplemented -- uncopyable
	SequenceConstraintFactory( SequenceConstraintFactory const & );
	SequenceConstraintFactory const & operator = ( SequenceConstraintFactory const & );

private:
	static SequenceConstraintFactory * instance_;
	SequenceConstraintMap creator_map_;


}; // class SequenceConstraintFactory

/// @brief Templated class used to register the SequenceConstraintCreator into
//  the SequenceContraintFactory
template < class T > 
class SequenceConstraintRegistrator : public utility::factory::WidgetRegistrator < SequenceConstraintFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< SequenceConstraintFactory, T > parent;
public:
	SequenceConstraintRegistrator() : parent() {}
}; // SequenceConstriantRegistator


} // namespace devel
} // namespace constrained_sequence_design 

#endif
