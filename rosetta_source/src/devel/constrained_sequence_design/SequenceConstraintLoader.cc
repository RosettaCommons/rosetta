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

// Unit Headers
#include <devel/constrained_sequence_design/SequenceConstraintLoader.hh>
#include <devel/constrained_sequence_design/SequenceConstraintLoaderCreator.hh>

// Package Headers
#include <devel/constrained_sequence_design/SequenceConstraintFactory.hh> 

// Project Headers
#include <basic/Tracer.hh>
#include <protocols/moves/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Boost Headers
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace devel {
namespace constrained_sequence_design {

static basic::Tracer TR( "devel.constrained_sequence_design.SequenceConstraintLoader" );
				
SequenceConstraintLoader::SequenceConstraintLoader() {}
SequenceConstraintLoader::~SequenceConstraintLoader() {}

void SequenceConstraintLoader::load_data(
	core::pose::Pose const & pose,
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data
) const
{
	using namespace utility::tag;
	using namespace utility;
	typedef utility::vector0< TagPtr > TagPtrs;


	foreach(TagPtr set_tag, tag->getTags() ){
		std::string const cset_name( set_tag->getName() );
	  TR << "Adding new sequence constriant set named " << cset_name << std::endl;
		SequenceConstraintSetOP constraints = new SequenceConstraintSet;
		foreach( TagPtr constraint_tag, set_tag->getTags() ) {
			std::string const type( constraint_tag->getName() );
			TR << "Adding new constriant type " << type << " to constraint set " << cset_name << std::endl;
			SequenceConstraintOP cntr = SequenceConstraintFactory::get_instance()->newSequenceConstraint( constraint_tag, data, pose );
			constraints->insert(cntr);
		}
		data.add("sequence_constraints", cset_name, constraints);
	}
}

// Loader Creator
protocols::jd2::parser::DataLoaderOP
SequenceConstraintLoaderCreator::create_loader() const { return new SequenceConstraintLoader; }

std::string
SequenceConstraintLoaderCreator::keyname() const { return "SEQUENCECONSTRAINTS"; }


} // namespace devel
} // namespace constrained_sequence_design
