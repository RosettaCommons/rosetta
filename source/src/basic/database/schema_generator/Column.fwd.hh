// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/database/schema_generator/Column.fwd.hh
/// @brief forward hearder for the Column class in the schema generator framework
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_database_schema_generator_Column_fwd_hh
#define INCLUDED_basic_database_schema_generator_Column_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace basic {
namespace database {
namespace schema_generator {

class Column;
typedef utility::vector1< Column > Columns;

typedef utility::pointer::shared_ptr< Column > ColumnOP;
typedef utility::pointer::shared_ptr< Column const > ColumnCOP;

} //namesapce
} //namespace
} //namespace

#endif //include guard
