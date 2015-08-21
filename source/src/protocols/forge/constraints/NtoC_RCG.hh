// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/constraints/NtoC_RCG.hh
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012

#ifndef INCLUDED_protocols_forge_constraints_NtoC_RCG_hh
#define INCLUDED_protocols_forge_constraints_NtoC_RCG_hh

// Unit Header
#include <protocols/forge/constraints/NtoC_RCG.fwd.hh>

// Package Header
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

// Proeject Header
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace constraints {

class NtoC_RCG : public protocols::forge::remodel::RemodelConstraintGenerator{
public:

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

public:

	NtoC_RCG();

	NtoC_RCG( NtoC_RCG const & rval );

	NtoC_RCG( Real const dist, Real const coef );

	virtual ~NtoC_RCG();

	virtual void
	parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual std::string
	get_name() const;

	virtual protocols::moves::MoverOP
	fresh_instance() const;

	virtual protocols::moves::MoverOP
	clone() const;

	virtual
	void generate_remodel_constraints( Pose const & pose );

	void set_weight( Real const coef );

	void set_distance( Real const dist );

private:

	Real dist_;
	Real coef_;

}; //class NtoC_RCG


} //namespace constraints
} //namespace forge
} //namespace protocols


#endif // INCLUDED_protocols_forge_constraints_NtoC_RCG_HH
