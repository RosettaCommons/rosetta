// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Record.cc
/// @brief  Helper function definitions for Record data structures.
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/Record.hh>


namespace core {
namespace io {
namespace pdb {

std::ostream &
operator<<( std::ostream & os, Record const & record )
{
	for ( Record::const_iterator field = record.begin(), end = record.end(); field != end; ++field ) {
		os << "<Record>{" << field->first << ":" << field->second << "}" << std::endl;
	}
	return os;
}

}  // namespace pdb
}  // namespace io
}  // namespace core
