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
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <protocols/sewing/sampling/LoophashAssemblyMover.hh>
#include <protocols/sewing/util/io.hh>
#include <protocols/sewing/conformation/Model.hh>
#include <protocols/sewing/conformation/ContinuousAssembly.hh>
#include <protocols/sewing/sampling/SewGraph.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <utility/io/ozstream.hh>

#include <protocols/loops/Loop.hh>

// C++ headers
#include <iostream>

#include <devel/init.hh>

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::sewing;

	// initialize core and read options
	devel::init(argc, argv);
	std::string model_file_name = option[sewing::model_file_name];
	std::string score_file_name = option[sewing::score_file_name];

	std::map< int, Model > models = read_model_file(model_file_name);
	utility::vector1<BasisPair> alignment_pairs = read_hashing_scores_from_file(score_file_name);
	SewGraphCOP graph = new SewGraph(models, alignment_pairs);

	utility::io::ozstream file("sew.json", std::ios::out | std::ios::binary);
	file << serialize_graph_json(graph, 200);

//	std::map< int, Model >::const_iterator it = models.begin();
//	std::map< int, Model >::const_iterator it_end = models.end();
//	for(; it != it_end; ++it) {
//		ContinuousAssembly assem(it->second);
//		assem.to_pose(core::chemical::CENTROID).dump_pdb(utility::to_string(it->first) + ".pdb");
//	}

//	using namespace protocols::loophash;
//	using namespace std;
//	
//	// initialize core
//	devel::init(argc, argv);
//
//	core::chemical::ResidueTypeSetCAP res_type_set =
//		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
//
//	core::pose::PoseOP pose = core::import_pose::pose_from_pdb(*res_type_set, "/Users/tjacobs2/PROJECTS/sewing2/refactor/discontinuous/1yzm_test/inputs/1yzm_no_loop.pdb");
//
//	for(core::Size i=1; i <= pose->total_residue(); ++i) {
//		core::conformation::remove_upper_terminus_type_from_conformation_residue( pose->conformation(), i );
//		core::conformation::remove_lower_terminus_type_from_conformation_residue( pose->conformation(), i );
//	}
//
//	protocols::sewing::LoophashAssemblyMover lh_mover;
//	lh_mover.add_single_loop(*pose, 18, 1, pose->total_residue());

	
//	utility::vector1<core::Size> loop_sizes;
//	loop_sizes.push_back(3);
//	LoopHashLibraryOP lh_library = new LoopHashLibrary( loop_sizes );
//	lh_library->load_mergeddb();
//	
//	LoopHashMap &hashmap = lh_library->gethash( 3 );
//	std::cout << "Number of loops in db of size 3: " << hashmap.n_loops() << std::endl;
//	
//	core::pose::PoseOP pose = core::import_pose::pose_from_pdb("1tua_clean_mpm.pdb");
//	
//	numeric::geometry::hashing::Real6 loop_transform;
//	//TR << "Getting transform from residues " << lh_fragment_begin << " and " << lh_fragment_end+1 << std::endl;
//	get_rt_over_leap_without_foldtree_bs( *pose, 5, 8, loop_transform );
//
//	std::vector<core::Size> leap_index_bucket;
//	hashmap.radial_lookup( 4, loop_transform, leap_index_bucket);
//	
//	cout << "Found " << leap_index_bucket.size() << " fragment for transform" << std::endl;
//
//	LeapIndex cp = hashmap.get_peptide(leap_index_bucket[0]);
//	BackboneSegment bs_seg;
//	lh_library->backbone_database().get_backbone_segment(cp.index, cp.offset,3,bs_seg);
//
//	bs_seg.print();
//
//	BBData bb;
//	BBExtraData extra_data;
//	lh_library->backbone_database().get_protein(cp.index, bb);
//	lh_library->backbone_database().get_extra_data(bb.extra_key, extra_data);
//	std::string sequence = extra_data.sequence;
//
//	cout << "PDB: " << extra_data.pdb_id << std::endl;
//	cout << "Sequence " << sequence << std::endl;
//	cout << "Offset " << cp.offset << std::endl;
//
//	core::Size seq_offset = (cp.offset/3)+(0*3);//should I be subtracting 1 here?
//	char aa_char = sequence[seq_offset+0];
//	cout << "aa 1 " << aa_char << std::endl;
//	aa_char = sequence[seq_offset+1];
//	cout << "aa 2 " << aa_char << std::endl;
//	aa_char = sequence[seq_offset+2];
//	cout << "aa 3 " << aa_char << std::endl;
//
//	return 0;
}
