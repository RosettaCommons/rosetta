// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief
/// @author Javier Castellanos	(javiercv@uw.edu)
///

#ifndef INCLUDED_devel_constrained_sequence_design_SequenceConstraint_HH
#define INCLUDED_devel_constrained_sequence_design_SequenceConstraint_HH

// Unit header
#include <devel/constrained_sequence_design/SequenceConstraint.fwd.hh>
// Package headers

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace devel {
namespace constrained_sequence_design {


class SequenceConstraint : public utility::pointer::ReferenceCount {

public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::pose::PoseAP PoseAP;
	typedef core::pose::PoseCOP PoseCOP;
	typedef core::Real Real;
	typedef core::Size Size;

public:
	/// @brief constructor
	SequenceConstraint();
	/// @brief constructor
	SequenceConstraint( PoseOP pose );

	/// @brief destructor
	virtual ~SequenceConstraint();

	/// @brief reuturns a name for the constraint
	virtual std::string name() = 0;

	/// @brief adjust any data that needs to be cached to the actual
	/// state of the pose. should be called only once and after a
	/// design run. 
	/// IMPORTANT: if you overwrite this function allways call
	/// the parent's class update function.
	virtual void update(const Pose & pose) = 0;


	/// @brief returns -1 if rejects aminoacid AA for position
	/// pos, 1 if it accepts it and 0 if neither. this function
	/// is const because any internal data that should be cached
	/// should be calculated only after design cycles, not every 
	/// time apply is called. Any constraint that has to cache 
	/// data has to calculate and load the data from update().
	virtual Size apply(core::chemical::AA aa, Size pos) const  = 0;

	/// @brief sets the state of the object according to the tag
	virtual 
	void
	parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	Pose const & pose );
	// Setters

	/// @brief set constrained pose
	void set_pose(Pose& p);

	// Accessors

	/// @brief return the  sum of applying the apply functiond
	/// over the whole protein using the sequence present on the 
	/// pose.
	virtual Real raw_score() ;

	/// @brief return the weighted score (raw_score * weight)
	Real score() ;

	/// @brief return the weight
	inline Real weight() const { return weight_; } 

private:
	PoseAP	pose_;
	Real		weight_;
};

} // constrained_sequence_design
} // devel

#endif
