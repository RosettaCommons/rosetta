// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//Part of the flatland docking protocol, symmetrizies sequence onto other modular parts of the building blocks. Rejects a building block if it is not modular. 
//Also generates the full lattice 

#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <cmath>
#include <iostream>
//#include "boost/filesystem.hpp" 
//#include <boost/algorithm/string.hpp>

#include <basic/datacache/DataMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <boost/foreach.hpp>

#include <devel/init.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/constants.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>

#include <core/id/AtomID_Map.hh>

#include <core/conformation/Residue.hh>
//#include <core/conformation/symmetry/SymmetricConformation.hh>
//#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/import_pose/import_pose.hh>

//#include <core/io/silent/ProteinSilentStruct.hh>
//#include <core/io/silent/BinaryProteinSilentStruct.hh>
//#include <core/io/silent/SilentFileData.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>

//#include <core/optimization/MinimizerOptions.hh>
//#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

//#include <core/pack/task/ResfileReader.hh>
//#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/operation/TaskOperations.hh>
//#include <core/pack/task/operation/util/interface_vector_calculate.hh>

#include <core/pose/motif/reference_frames.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/Remarks.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/util.hh>

//#include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/hbonds/HBondOptions.hh>
//#include <core/scoring/hbonds/hbonds.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
//#include <core/scoring/motif/util.hh>
//#include <core/scoring/motif/xfrags.hh>
//#include <core/scoring/motif/motif_hash_stuff.hh>
//#include <core/scoring/sasa.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
//#include <core/scoring/symmetry/SymmetricEnergies.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreFunction.hh>

#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

//#include <protocols/cryst/refineable_lattice.hh>
//#include <protocols/cryst/TwoDLattice.hh>
//#include <protocols/cryst/wallpaper.hh>

#include <protocols/moves/Mover.hh>

//#include <protocols/flxbb/LayerDesignOperation.hh>

//#include <protocols/filters/Filter.hh>

//#include <protocols/relax/FastRelax.hh>

//#include <protocols/simple_moves/MakePolyXMover.hh>
//#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
//#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
//#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

//#include <protocols/toolbox/SelectResiduesByLayer.hh>
//#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
//#include <protocols/toolbox/task_operations/SelectBySASAOperation.hh>

//#include <core/pack/interaction_graph/InteractionGraphBase.hh>


using namespace basic;
using namespace core;
using namespace core::pose;
using namespace core::conformation::symmetry;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//OPT_KEY( Integer, cn)
//OPT_KEY( String, wallpaper)
// OPT_KEY( Real, debug3)
//OPT_KEY( Boolean, dump_silent)
//OPT_KEY( Boolean, read_changes_from_cmdline)
OPT_KEY( String, src_pose_1_start)
OPT_KEY( String, target_pose_start)
OPT_KEY( Integer, span)
OPT_KEY( String, prefix)
OPT_KEY( Boolean, pad_3_heptads)
OPT_KEY( String, target_pose_start_2)
OPT_KEY( Boolean, align_bb)
// OPT_KEY( String, reference)

static basic::Tracer TR("apps.pilot.zibochen.pad_helices");

////////////////////////////////////////////////
int
main( int argc, char ** argv ) {

	NEW_OPT( src_pose_1_start, "starting indecies of the source pose 1", "");
	NEW_OPT( target_pose_start, "starting indecies of the target pose for src pose 2 ", "");
	NEW_OPT( span, "how long are the replaced regions?", 7);
	NEW_OPT( prefix, "prefix for the output file", "");
	NEW_OPT( pad_3_heptads, "pad 2 or 3 heptads?", false);
	NEW_OPT( target_pose_start_2, "starting indecies of the target pose for src pose 2 ", "");
	NEW_OPT( align_bb, "should backbones be aligned during residue replacement?", false);

	devel::init(argc, argv);

	// option[ in::preserve_crystinfo ].value(true); //need this line to add remarks to pdb
	// option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(false); //force this option to be false so Rosetta would not try to connect bonds between adjacent jumps

	utility::vector1< std::string > design_filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
	core::pose::Pose src_pose;
	core::import_pose::pose_from_file(src_pose, design_filenames[1] );
	core::pose::Pose target_pose;
	core::import_pose::pose_from_file(target_pose, design_filenames[2] );

	//get basename, used for outputing pose
	std::string base_name=design_filenames[2].substr( design_filenames[2].find_last_of( '/' ) +1 );
	const std::string ext(".pdb");
	base_name = base_name.substr(0, base_name.size() - ext.size());

	//string to size
	utility::vector1< Size > src_pose_1_start_ind;
	utility::vector1< Size > target_pose_start_ind;
	utility::vector1<std::string> src_pose_1_start_in=utility::string_split_multi_delim(option[src_pose_1_start],",");
	for (Size i=1;i<=src_pose_1_start_in.size();i++){
		src_pose_1_start_ind.push_back(utility::string2int(src_pose_1_start_in[i]));
	}
	utility::vector1<std::string> target_pose_start_in=utility::string_split_multi_delim(option[target_pose_start],",");
	for (Size i=1;i<=target_pose_start_in.size();i++){
		target_pose_start_ind.push_back(utility::string2int(target_pose_start_in[i]));
	}

	//replace residues
	for (Size i=1;i<=src_pose_1_start_ind.size();i++){
		for (int j=0;j<=option[span]-1;j++){
			target_pose.replace_residue(target_pose_start_ind[i]+j, src_pose.residue(src_pose_1_start_ind[i]+j), option[align_bb]);
			// target_pose.delete_polymer_residue(target_pose_start_ind[i]+j);
			// //target_pose.append_residue_by_bond(src_pose.residue(src_pose_1_start_ind[i]+j));
			// target_pose.append_residue_by_jump(src_pose.residue(src_pose_1_start_ind[i]+j),target_pose_start_ind[i]+j-1);
		}
		//target_pose.copy_segment(option[span],src_pose,target_pose_start_ind[i],src_pose_1_start_ind[i]);
	}

	//check if there is an additional heptad to be padded
	if (option[pad_3_heptads]){
		//string to size
		utility::vector1< Size > target_pose_start_2_ind;
		utility::vector1<std::string> target_pose_start_2_in=utility::string_split_multi_delim(option[target_pose_start_2],",");
		for (Size i=1;i<=target_pose_start_2_in.size();i++){
			target_pose_start_2_ind.push_back(utility::string2int(target_pose_start_2_in[i]));
		}
		core::pose::Pose src_pose_2;
		core::import_pose::pose_from_file(src_pose_2, design_filenames[3] );
		//replace residues
		for (Size i=1;i<=target_pose_start_2_ind.size();i++){
			for (int j=0;j<=option[span]-1;j++){
				target_pose.replace_residue(target_pose_start_2_ind[i]+j, src_pose_2.residue(src_pose_1_start_ind[i]+j), option[align_bb]);
				// target_pose.delete_polymer_residue(target_pose_start_ind[i]+j);
				// //target_pose.append_residue_by_bond(src_pose.residue(src_pose_1_start_ind[i]+j));
				// target_pose.append_residue_by_jump(src_pose.residue(src_pose_1_start_ind[i]+j),target_pose_start_ind[i]+j-1);
			}
			//target_pose.copy_segment(option[span],src_pose,target_pose_start_ind[i],src_pose_1_start_ind[i]);
		}
	}
	//dump pdb
	target_pose.dump_pdb(option[prefix]+base_name+"_padded.pdb");
	return 0;
}//main
