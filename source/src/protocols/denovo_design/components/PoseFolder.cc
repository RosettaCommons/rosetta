// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/components/PoseFolder.cc
/// @brief Given a pose with all residues, assign a 3D conformation to that pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit header or inline function header
#include <protocols/denovo_design/components/PoseFolder.hh>

namespace protocols {
namespace denovo_design {
namespace components {

PoseFolder::PoseFolder( std::string const & type ):
	utility::pointer::ReferenceCount(),
	type_( type )
{}

// Defined to prevent pure virtual destructor error at run time.
PoseFolder::~PoseFolder()
{}

std::string const &
PoseFolder::type() const
{
	return type_;
}

void
PoseFolder::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	parse_tag( tag, data );
}

} //protocols
} //denovo_design
} //components
