// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ken Jung
/// @brief  lhalign
/// Given a partial threaded model pdb and a fasta, outputs the bin radius required to find a fragment
/// to fit inside each gap. Does not check for clashes between gaps for now.

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>


#include <core/scoring/rms_util.hh>


#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// C++ headers
#include <core/init.hh>
#include <map>
#include <ctime>

using namespace core;
using namespace io;
using namespace pdb;
using namespace chemical;
using namespace conformation;
using namespace protocols::loophash;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility;

using core::pose::PoseOP;
using core::pose::PoseAP;
using core::pose::Pose;
using utility::vector1;
using core::Size;
using std::map;
using std::string;

static basic::Tracer TR("main");

// lifted straight from protocols/idealize/idealize.cc
void find_next_gap( Pose & pose, Size & idx, Real & gap_distance ) {
	// squared distance at which bond is considered discontinuous
	Real const chain_break_cutoff = { 4.0 };
	Size const nres ( pose.total_residue() );

	// find chain breaks to add to gaplist
	pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	kinematics::FoldTree f( pose.fold_tree() );
	for ( Size i = idx + 1; i < nres; ++i ) {
		if ( f.is_cutpoint(i) ) continue;
		bool chain_break = false;
		Size j = i+1;
		if ( pdbinfo->number(i)+1 != pdbinfo->number(j) ) {
			TR.Info << "non-sequential at res nums " << i << '-' << j << std::endl;
			TR.Info << "non-sequential pdb res nums " << pdbinfo->number(i) << pdbinfo->chain(i) <<
					'-' << pdbinfo->number(j) << pdbinfo->chain(j) << std::endl;
			chain_break = true;
		} else {
			conformation::Residue const & rsd = pose.residue(i);
			conformation::Residue const & next_rsd = pose.residue(j);
			if (rsd.is_polymer() && next_rsd.is_polymer()) {
				Real dist_squared = rsd.atom( rsd.upper_connect_atom() ).xyz().distance_squared(next_rsd.atom( next_rsd.lower_connect_atom() ).xyz());
				gap_distance = std::sqrt(dist_squared);
				if (dist_squared > chain_break_cutoff) {
						//TR.Info << "chain break at res nums: " << i << '-' << j << ' ' << std::sqrt(dist_squared) << std::endl;
						//TR.Info << "chain break pdb res nums: " << pdbinfo->number(i) << pdbinfo->chain(i) <<
						//		'-' << pdbinfo->number(j) << pdbinfo->chain(j) << std::endl;
						chain_break = true;
				} else if ( dist_squared < 0.1 ) {
						//TR.Info << "zero length bond at res nums: " << i << '-' << j << std::endl;
						//TR.Info << "zero length bond pdb res nums: " << pdbinfo->number(i) << pdbinfo->chain(i) <<
						//	'-' << pdbinfo->number(j) << pdbinfo->chain(j) << std::endl;
						chain_break = true;
				}
			}
		}
		if ( chain_break ) {
			idx = i;
			return;
		}
	}
	idx = 0;
}

std::map< std::string, core::pose::Pose > poses_from_cmd_line(utility::vector1< std::string > const & fn_list){

  using utility::file::file_exists;
  using core::import_pose::pose_from_pdb;

  ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
  std::map< std::string, core::pose::Pose > poses;
  typedef utility::vector1< std::string >::const_iterator iter;
  for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
    if ( file_exists(*it) ) {
                        core::pose::Pose pose;
      core::import_pose::pose_from_pdb( pose, *rsd_set, *it );
                        std::string name = utility::file_basename( *it );
      name = name.substr( 0, 5 );
      poses[name] = pose;
    }
  }
  return poses;
}

int
main( int argc, char * argv [] )
{

	const Size max_loop_size = 14;
	const Size min_loop_size = 3;

	// initialize core
	core::init(argc, argv);
	
  // read in poses
  map< string, Pose > input_poses = poses_from_cmd_line(option[ in::file::s ]());

	// initialize lhlibrary
	// max gap size we can handle is 14
	utility::vector1 < core::Size > loop_sizes;
	for (Size i = min_loop_size; i<= max_loop_size; i++ ) {
		loop_sizes.push_back(i);
	}
	LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
	LoopHashLibraryOP lhlibrary = new LoopHashLibrary( loop_sizes );
	lhlibrary->load_mergeddb();
	LoopHashSamplerOP lhsampler = new LoopHashSampler( lhlibrary, simple_inserter );
	// make this more stringent?
	lhsampler->set_max_rms( 0.3 );
	lhsampler->set_nonideal( true );
	lhsampler->set_min_rms( 0.0 );
	
	long starttime = time(NULL);
  typedef map< string, Pose >::const_iterator iter;
	for( iter it = input_poses.begin(), end = input_poses.end(); it != end; ++it ) {

		TR << "Working on pose " << it->first << std::endl;
		// copy pose
		PoseOP working_pose = new Pose( it->second );
		Size idx = 1;

		// convert pose to centroid pose:
		if( working_pose->is_fullatom() ){
			core::util::switch_to_residue_type_set( *working_pose, core::chemical::CENTROID);
		}
		// Now go through each gap and try increasingly larger lh until something is returned
		Size next_gap = 0;
		Real gap_dist;
		find_next_gap( *working_pose, next_gap, gap_dist );
		while( next_gap != 0 ) {
			TR << "Attempting to fix gap following residue " << next_gap << std::endl;
			std::vector< Pose > lib_structs;
			
			// gogo loophash
			// increase loophash size until we get anything returned
			Size loop_size = std::max((Size)(gap_dist/3.5), min_loop_size); // no point in trying anything that cant reach across the gap
			while( lib_structs.size() == 0 && loop_size < max_loop_size ) {
				TR << "Trying loopsize " << ++loop_size << std::endl;
				lhsampler->set_start_res( next_gap + 3 < loop_size ? 0 : next_gap + 3 - loop_size );
				lhsampler->set_stop_res ( next_gap );
				lhsampler->close_gaps( *working_pose, lib_structs, loop_size );
			}
			if( lib_structs.size() != 0 ) {
				working_pose = new Pose (lib_structs[0]);
				string n = it->first + ".gap" + to_string(next_gap) + ".pdb";
				working_pose->dump_pdb(n);
			}

			find_next_gap( *working_pose, next_gap, gap_dist );
		}

	// dump resulting pdb
	string m = it->first + ".closed.pdb";
	working_pose->dump_pdb(m);
	}

	long endtime = time(NULL);
	TR << "Gap closing took " << to_string(endtime - starttime) << "seconds" << std::endl;
}


