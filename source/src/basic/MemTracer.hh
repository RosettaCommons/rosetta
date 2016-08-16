// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_basic_MemTracer_hh
#define INCLUDED_basic_MemTracer_hh

// Utility headers
#include <basic/Tracer.hh>
#include <iosfwd>


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
