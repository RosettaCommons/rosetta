// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /protocols/filters/filters.fwd.hh
/// @brief  forward declaration for Filter base class
/// @author Florian Richter, floric@u.washington.edu, Sarel Fleishman sarelf@u.washington.edu


#ifndef INCLUDED_protocols_filters_BasicFilters_fwd_hh
#define INCLUDED_protocols_filters_BasicFilters_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>
// AUTO-REMOVED #include <string>


namespace protocols {
namespace filters {

// Forward
class TrueFilter;
class FalseFilter;
class CompoundFilter;
class CombinedFilter;
class MoveBeforeFilter;
class IfThenFilter;

// Types
typedef utility::pointer::shared_ptr< CompoundFilter > CompoundFilterOP;
typedef utility::pointer::shared_ptr< CompoundFilter const >  CompoundFilterCOP;

typedef utility::pointer::shared_ptr< CombinedFilter > CombinedFilterOP;
typedef utility::pointer::shared_ptr< CombinedFilter const >  CombinedFilterCOP;

typedef utility::pointer::shared_ptr< MoveBeforeFilter > MoveBeforeFilterOP;
typedef utility::pointer::shared_ptr< MoveBeforeFilter const >  MoveBeforeFilterCOP;

typedef utility::pointer::shared_ptr< IfThenFilter > IfThenFilterOP;
typedef utility::pointer::shared_ptr< IfThenFilter const >  IfThenFilterCOP;

// used by CompoundFilter
enum boolean_operations {
	AND=1,
	OR,
	XOR,
	NOR,
	NAND,
  ORNOT,
	ANDNOT,
	NOT
};


} // namespace protocols
} // namespace filters

#endif
