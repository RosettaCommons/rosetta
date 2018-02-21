// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/gpu/Timer.cc
/// @brief  High-resolution Timer (ns resolution)
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

#ifndef INCLUDED_basic_gpu_TIMER_cc
#define INCLUDED_basic_gpu_TIMER_cc

#include <basic/gpu/Timer.hh>
#include <iostream>

namespace basic {
namespace gpu {

Timer::Timer() {
	tag_ = nullptr;
	t_ = nullptr;
	Reset();
}

Timer::~Timer() {
	Report(tag_);
}

Timer::Timer(const char *tag) {
	t_ = nullptr;
	tag_ = tag;
	Reset();
}

Timer::Timer(basic::Tracer::TracerProxy& t, const char *tag) {
	t_ = &t;
	tag_ = tag;
	Reset();
}

#ifdef WIN32
void Timer::Reset() {
	start = clock();
}

double Timer::GetTime() {
	end = clock();
	double time = 1000 * (end-start) / CLOCKS_PER_SEC;
	return time;
}

#else

void Timer::Reset() {
	gettimeofday(&start, nullptr);
}

double Timer::GetTime() {
	gettimeofday(&end, nullptr);
	double time = ((end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.);
	return time;
}
#endif

void Timer::Report(const char *tag) {

	double time = GetTime();
	if ( t_ ) {
		if ( tag ) {
			(*t_) << "Time [" << tag << "]: ";
		} else {
			(*t_) << "Time: ";
		}
		(*t_) << time << " ms" << std::endl;
	} else {
		if ( tag_ ) {
			std::cout << "Time [" << tag << "]: ";
		} else {
			std::cout << "Time: ";
		}
		std::cout << time << " ms" << std::endl;
	}
}

} // gpu
} // basic

#endif // INCLUDED_basic_gpu_TIMER_cc
