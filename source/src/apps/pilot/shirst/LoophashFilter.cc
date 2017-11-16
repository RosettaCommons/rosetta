// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/shirst/LoophashFilter.cc
/// @brief  trying to figure how to use loophash as a filter in folding with the TMHTopologySamplerClaimer
/// @author Stephanie DeLuca (stephanie.h.deluca@vanderbilt.edu)

//Loophash
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <utility/excn/EXCN_Base.hh>
//Scoring
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoreFunction.hh>

//Pose
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <sstream>

#include <devel/init.hh>

static basic::Tracer tr( "apps.pilot.LoophashFilter" );

//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	try {
		devel::init(argc, argv);

		using namespace basic::options::OptionKeys;
		using basic::options::option;

		//pose is read in from PDB for now
		utility::vector0<std::string> pdbs;
		pdbs = option[in::file::s]();
		std::string pdb=pdbs[0];
		core::pose::Pose pose; // starts NULL, coords *never* modified!
		core::import_pose::pose_from_file(pose, pdb, core::import_pose::PDB_file);
		core::Size nres = pose.size();

		//membrane stuff to figure out loop start and stop
		//At this point, assert if don't have a spanfile
		std::string spanfile = option[in::file::spanfile];
		runtime_assert (option[in::file::spanfile].user());

		//set up membrane topology
		core::scoring::MembraneTopologyOP topology( new core::scoring::MembraneTopology );
		pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topology );
		topology->initialize(spanfile);

		utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
		protocols::loophash::LoopHashLibraryOP library( new protocols::loophash::LoopHashLibrary( loop_sizes ) );

		//loophash stuff
		core::Size start = 0;
		core::Size stop = 0;

		bool rt_exists(false);

		// [1] get transform from start to stop
		numeric::geometry::hashing::Real6 loop_transform; //this is the RT for the loop in question (how far apart are the residues)

		library->load_db();

		// initialize some variables for setting up loops array
		core::Size const nspan = topology->tmhelix();
		core::Size const num_cut_loops = nspan-1;
		ObjexxFCL::FArray1D_int previous_span_begin((nspan-1),0);
		ObjexxFCL::FArray1D_int previous_span_end((nspan-1),0);
		ObjexxFCL::FArray1D_int span_begin((nspan-1),0);
		ObjexxFCL::FArray1D_int span_end((nspan-1),0);
		ObjexxFCL::FArray1D_int loop_begin(num_cut_loops,0);
		ObjexxFCL::FArray1D_int loop_end(num_cut_loops,0);
		core::Size span_index = 1;

		//Loop through spans to figure out loop defs. For each loop, figure out loop hashing
		for ( span_index = 1; span_index <= nspan-1; ++span_index ) {
			//need to know the beginning and end of the TMs to determine where the loops are
			previous_span_begin(span_index) = topology->span_begin(span_index);
			previous_span_end(span_index) = topology->span_end(span_index);
			span_begin(span_index) = topology->span_begin(span_index+1);
			span_end(span_index) = topology->span_end(span_index+1);

			loop_begin(span_index) = previous_span_end(span_index);
			loop_end(span_index) = span_begin(span_index);

			//if the predicted loop (that is, not span) is shorter than 3 res
			if ( loop_end(span_index) - loop_begin(span_index) < 2 ) {
				loop_begin(span_index) -= 1;
				loop_end(span_index) += 1;
			}

			//Print loop begin and end
			tr.Info << "span_index:  " << span_index << " loop_begin:  " << loop_begin(span_index) << " loop_end: "
				<< loop_end(span_index) << std::endl;

			assert(loop_begin(span_index) != 0);
			assert(loop_end(span_index) != 0);

			//loophash stuff
			start = loop_begin(span_index);
			stop = loop_end(span_index);

			if ( start > nres || start < 1 || stop > nres || stop < 1 ) {
				tr.Info << "ERROR!" << "residue range " << start << "," << stop << "unknown!" << std::endl;
			}
		}
		rt_exists = protocols::loophash::get_rt_over_leap_fast( pose, start, stop, loop_transform );

		core::Size radius_size(0);
		if ( option[lh::radius_size].user() ) {
			radius_size = (core::Size)(option[lh::radius_size].value());
		} else {
			radius_size = 2;
		}

		//Loop through loop library (hashes) and figure out, for the RT, does it return any hashes? (i.e., is there a loop in the DB that has this RT)
		for ( std::vector< core::Size >::const_iterator jt = library->hash_sizes().begin(); jt != library->hash_sizes().end(); ++jt ) {
			core::Size loop_size = *jt;
			if ( rt_exists==true ) {
				// Get the fragment bucket
				// leap_index_bucket contains loophash hits (see pilot app to extract the backbone)
				protocols::loophash::LoopHashMap &hashmap = library->gethash( loop_size );
				std::vector < core::Size > leap_index_bucket;
				tr.Info << "radius_size:  " << radius_size << "\tloop_size:  " << loop_size << "\tnumber of hits:  " << hashmap.radial_count(radius_size,loop_transform) << std::endl;
			}
		}
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;

}
