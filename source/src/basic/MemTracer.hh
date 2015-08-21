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
/// @author

#ifndef INCLUDED_basic_MemTracer_hh
#define INCLUDED_basic_MemTracer_hh

// C/C++ headers
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <ostream>
#include <sstream>
#include <vector>
#include <basic/Tracer.fwd.hh>


namespace basic {

void get_usage_from_procfilesystem( std::ostream& mem_report );

class MemTracer : public Tracer {
	typedef Tracer Parent;

public:
	MemTracer( TracerPriority priority=t_info, bool muted_by_default = true )
	: Tracer( "memory_usage", priority, muted_by_default ) {}

protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const &);

private:
	static bool single_line_;
};

extern MemTracer mem_tr;

}  // namespace basic

#endif
