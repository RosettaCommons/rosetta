// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009
/// @modified Tom Linsky, tlinsky@uw.edu

#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.forge.remodel.remodelconstraintgenerator" );

namespace protocols {
namespace forge {
namespace remodel {

/// @details Auto-generated virtual destructor
RemodelConstraintGenerator::~RemodelConstraintGenerator()
{}


RemodelConstraintGenerator::RemodelConstraintGenerator()
: protocols::moves::ConstraintGenerator(),
	seqmap_(/* NULL */),
	vlb_(/* NULL */)
{}

/// @brief This is called if this mover is instantiated from XML
void
RemodelConstraintGenerator::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	ConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
}

RemodelConstraintGenerator::VarLengthBuildAP
RemodelConstraintGenerator::vlb() const {
	return vlb_;
}

void
RemodelConstraintGenerator::set_vlb(
	VarLengthBuildAP vlb )
{
	vlb_ = vlb;
}

void
RemodelConstraintGenerator::set_seqmap(
	core::id::SequenceMappingCOP seqmap )
{
	seqmap_ = seqmap;
}

core::id::SequenceMappingCOP
RemodelConstraintGenerator:: seqmap() const
{
	return seqmap_;
}

} //namespace remodel
} //namespace forge
} //namespace protocols
