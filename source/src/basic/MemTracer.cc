// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit header
#include <basic/MemTracer.hh>

#include <stdio.h>              // for fclose, fopen, fscanf, snprintf, EOF
#ifndef WIN32
#include <sys/resource.h>       // for getrusage, rusage, RUSAGE_SELF
#include <unistd.h>             // for getpagesize, getpid
#endif
#include <ObjexxFCL/format.hh>  // for I
#include <iomanip>              // for operator<<, setprecision
#include <iosfwd>               // for ostream, stringstream, string
#include <ostream>              // for operator<<, basic_ostream
#include <string>               // for allocator, char_traits

namespace basic {

// ---- from POSIX reference ---
// struct rusage {
//     struct timeval Ru_utime; /* user time used */
//     struct timeval ru_stime; /* system time used */
//     long   ru_maxrss;        /* maximum resident set size */
//     long   ru_ixrss;         /* integral shared memory size */
//     long   ru_idrss;         /* integral unshared data size */
//     long   ru_isrss;         /* integral unshared stack size */
//     long   ru_minflt;        /* page reclaims */
//     long   ru_majflt;        /* page faults */
//     long   ru_nswap;         /* swaps */
//     long   ru_inblock;       /* block input operations */
//     long   ru_oublock;       /* block output operations */
//     long   ru_msgsnd;        /* messages sent */
//     long   ru_msgrcv;        /* messages received */
//     long   ru_nsignals;      /* signals received */
//     long   ru_nvcsw;         /* voluntary context switches */
//     long   ru_nivcsw;        /* involuntary context switches */
// };

// it turns out LINUX kernel do not support getrusage() ...
void get_usage_from_procfilesystem( std::ostream& mem_report ) {
#if defined(_WIN32) || defined( __native_client__ )
	return;  // disabled on windows
#else
	char buf[30];
	int page_sz = getpagesize()/1024; //this value doesn't seem to be correct on BG

#ifdef MPICH_IGNORE_CXX_SEEK
	//now using MPICH_IGNORE_CXX_SEEK is a bit hacky, but this happens to be a flag active on BG/P
	page_sz = 1;
#endif

	unsigned pid = (unsigned)getpid();
	snprintf(buf, 30, "/proc/%u/statm", pid );
	FILE* pf = fopen(buf, "r");
	int const width( 4);
	if (pf) {
		unsigned size; //       total program size
		unsigned resident;//   resident set size
		unsigned share;//      shared pages
		unsigned text;//       text (code)
		unsigned lib;//        library
		unsigned data;//       data/stack
		// Width limits added to appease cppcheck:
		// 25 charachters should be sufficient for any number able to fit in an unsigned
		if (fscanf(pf, "%25u %25u %25u %25u %25u %25u", &size, &resident, &share, &text, &lib, &data ) == EOF){
			mem_report << "WARNING! End of file reached without assignments from fscanf!";
		}
		using namespace ObjexxFCL::format;
		mem_report << std::setprecision(4) << "Virt/Res/Share/Exe/Data "
							 << I( width, (int) size*4/1024 ) << "  "
							 << I( width, (int) resident*page_sz/1024)  << " "
							 << I( width, (int) share*4/1024 ) << " "
							 << I( width, (int) (lib+text)*4/1024 ) << "  "
							 << I( width, (int) data * 4/1024 ) << " MB";
		fclose(pf);
	} else {
		rusage usage;
		using namespace ObjexxFCL::format;
		getrusage( RUSAGE_SELF, &usage );
		mem_report << std::setprecision(4) << "Res/ixrss/rdrss/isrss/swap "

							 << I( width, (int) usage.ru_maxrss*page_sz/1024 ) //<< "  "
// 							 << I( width, (int) usage.ru_ixrss*4/1024  ) << " "
// 							 << I( width, (int) usage.ru_idrss*4/1024  ) << " "
// 							 << I( width, (int) usage.ru_isrss*4/1024  ) << "  "
//							 << I( width, (int) usage.ru_nswap*4/1024  )
							 << " MB";
	}
#endif
}

bool MemTracer::single_line_ = false;
void MemTracer::t_flush( std::string const &str ) {
#ifdef _WIN32
	return;  // disabled on windows
#else
	if ( visible() ) {
		std::stringstream mem_report;
		get_usage_from_procfilesystem( mem_report );
		mem_report << "   @ " << str;
		Parent::t_flush( mem_report.str() );
	} else {
		Parent::t_flush( str );
	}
#endif
}

MemTracer mem_tr;

} // basic
