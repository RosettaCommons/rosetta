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
/// @detailed
/// @author Oliver Lange
/// @author Mike Tyka
///


// Unit Headers
// AUTO-REMOVED #include <protocols/loops/util.hh>

// Package Headers
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>

#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS


static thread_local basic::Tracer TR( "protocols.loops.util" );
static thread_local basic::Tracer tr( "protocols.loops" );

namespace protocols {
namespace loops {
using namespace core;
using namespace pose;
using namespace kinematics;


void
loopfinder( core::pose::Pose & pose , loops::Loops & loops ){
	// this function will not find terminal loops
	// get dssp info
	core::Size lastres = pose.total_residue();
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	utility::vector1< utility::vector1< core::Size > >all_loops;
	utility::vector1< core::Size > ind_loops;
	for( core::Size ii = 1; ii < pose.total_residue(); ++ii ){
		if( pose.secstruct(ii) == 'L' && pose.secstruct(ii+1) == 'L' ) {
			ind_loops.push_back(ii);
		}
		if( pose.residue(ii).chain() != pose.residue(ii+1).chain() )  { // separate loops spanning multiple chains
			all_loops.push_back(ind_loops);
			ind_loops.clear();
			continue; // make sure we skip all further conditionals so that we don't double-add the loop
		}
		if( pose.secstruct(ii) == 'L' && pose.secstruct(ii+1) != 'L' ) {
			ind_loops.push_back(ii);
			all_loops.push_back(ind_loops);
			ind_loops.clear();
		}

	}

	if( pose.secstruct( lastres ) == 'L' ){
		ind_loops.push_back( lastres );
		all_loops.push_back(ind_loops);

		ind_loops.clear();
	}


	for( core::Size ii = 1; ii <= all_loops.size(); ++ii ){
		core::Size const lastlooppos = all_loops[ii].size();
		core::Size const chain_firstpos = pose.residue( all_loops[ii][1] ).chain();
		core::Size const chain_lastpos = pose.residue( all_loops[ii][lastlooppos] ).chain();
		runtime_assert( chain_firstpos == chain_lastpos ); // This should have been caught in the previous loop indexing.
		core::Size const chain_begin = pose.conformation().chain_begin( chain_firstpos );
		core::Size const chain_end = pose.conformation().chain_end( chain_firstpos );

		// dont include terminal loops
		if ( all_loops[ii][1] == 1 || all_loops[ii][lastlooppos] == lastres ) { continue; }
		if ( all_loops[ii][1] == chain_begin || all_loops[ii][lastlooppos] == chain_end ) { continue; } // skip chain begin/end for multimers
		else{

			// make sure loop is at least 3 residues long
			if( all_loops[ii][lastlooppos] - all_loops[ii][1] < 3 ){
				TR << "increasing loop from" << all_loops[ii][1] << " " << all_loops[ii][lastlooppos] << std::endl;
				TR << "increasing loop to  " << all_loops[ii][1]-1 << " " << all_loops[ii][lastlooppos]+1 << std::endl;
				// this may seem sloopy but it catches the case where loop size is 1 and other cases as well
				core::Size cut_point = (all_loops[ii][lastlooppos]+1 - all_loops[ii][1]-1)/2 + all_loops[ii][1]-1;
				TR << "cut_point" <<cut_point << std::endl;
				loops.add_loop( all_loops[ii][1]-1, all_loops[ii][lastlooppos]+1, cut_point, 0, false );
			} else {
				core::Size cut_point = (all_loops[ii][lastlooppos] - all_loops[ii][1])/2 + all_loops[ii][1];
				loops.add_loop( all_loops[ii][1], all_loops[ii][lastlooppos], cut_point, 0, false );
				TR << "cut_point" <<cut_point << std::endl;
			}

		}

	}

}


} // loops
} // protocols
