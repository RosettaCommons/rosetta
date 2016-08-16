// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/gpu/Timer.hh
/// @brief  High-resolution Timer (ns resolution, but on *nix only)
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

#ifndef INCLUDED_basic_gpu_Timer_hh
#define INCLUDED_basic_gpu_Timer_hh

#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <basic/Tracer.hh>

namespace basic {
namespace gpu {

class Timer {
#ifdef WIN32
	clock_t start, end;
#else
	struct timeval start, end;
#endif
	const char *tag_;
	basic::Tracer::TracerProxy *t_;
public:
	Timer(const char *tag);
	Timer(basic::Tracer::TracerProxy& t, const char *tag =NULL);
	Timer();
	~Timer();
	void Report(const char *tag =NULL);
	void Reset();
	double GetTime();
};

} // gpu
} // basic

#endif // INCLUDED_basic_gpu_Timer_hh
