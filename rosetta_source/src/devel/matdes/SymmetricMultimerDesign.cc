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
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )

// Unit headers
#include <devel/matdes/SymmetricMultimerDesign.hh>
#include <devel/matdes/SymmetricMultimerDesignMoverCreator.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
std::string
ConstrainedDesignMoverCreator::keyname() const
{
	return ConstrainedDesignMoverCreator::mover_name();
}

MoverOP
ConstrainedDesignMoverCreator::create_mover() const {
	return new ConstrainedDesign;
}

std::string
ConstrainedDesignMoverCreator::mover_name()
{
	return "ConstrainedDesign";
}
// -------------  Mover Creator -------------


SymmetricMultimerDesign::SymmetricMultimerDesign() {
} //SymmetricMultimerDesign()


SymmetricMultimerDesign::SymmetricMultimerDesign(const SymmetricMultimerDesign& rval) {
} // SymmetricMultimerDesign(const SymmetricMultimerDesign& rval) 

void 
SymmetricMultimerDesign::apply(Pose& pose) {
} // apply

MoverOP SymmetricMultimerDesign::clone() const {
} // clone


MoverOP
SymmetricMultimerDesign::fresh_instance() const {
} // fresh_instance


void 
SymmetricMultimerDesign::parse_my_tag( TagPtr const tag,
										 DataMap & data,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & ) {
} // parse_my_tag


} // devel
} // matdes
