// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TestLoophash.cc
///
/// @brief

/// @author Tim Jacobs

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

// C++ headers
#include <iostream>

#include <devel/init.hh>

int
main( int argc, char * argv [] )
{

	using namespace protocols::loophash;
	using namespace std;
	
	// initialize core
	devel::init(argc, argv);
	
	utility::vector1<core::Size> loop_sizes;
	loop_sizes.push_back(3);
	LoopHashLibraryOP lh_library = new LoopHashLibrary( loop_sizes );
	lh_library->load_mergeddb();
	
	LoopHashMap &hashmap = lh_library->gethash( 3 );
	std::cout << "Number of loops in db of size 3: " << hashmap.n_loops() << std::endl;
	
	core::pose::PoseOP pose = core::import_pose::pose_from_pdb("1tua_clean_mpm.pdb");
	
	numeric::geometry::hashing::Real6 loop_transform;
	//TR << "Getting transform from residues " << lh_fragment_begin << " and " << lh_fragment_end+1 << std::endl;
	get_rt_over_leap_without_foldtree_bs( *pose, 5, 8, loop_transform );

	std::vector<core::Size> leap_index_bucket;
	hashmap.radial_lookup( 4, loop_transform, leap_index_bucket);
	
	cout << "Found " << leap_index_bucket.size() << " fragment for transform" << std::endl;

	LeapIndex cp = hashmap.get_peptide(leap_index_bucket[0]);
	BackboneSegment bs_seg;
	lh_library->backbone_database().get_backbone_segment(cp.index, cp.offset,3,bs_seg);

	bs_seg.print();

	BBData bb;
	BBExtraData extra_data;
	lh_library->backbone_database().get_protein(cp.index, bb);
	lh_library->backbone_database().get_extra_data(bb.extra_key, extra_data);
	std::string sequence = extra_data.sequence;

	cout << "PDB: " << extra_data.pdb_id << std::endl;
	cout << "Sequence " << sequence << std::endl;
	cout << "Offset " << cp.offset << std::endl;

	core::Size seq_offset = (cp.offset/3)+(0*3);//should I be subtracting 1 here?
	char aa_char = sequence[seq_offset+0];
	cout << "aa 1 " << aa_char << std::endl;
	aa_char = sequence[seq_offset+1];
	cout << "aa 2 " << aa_char << std::endl;
	aa_char = sequence[seq_offset+2];
	cout << "aa 3 " << aa_char << std::endl;

	return 0;
}
