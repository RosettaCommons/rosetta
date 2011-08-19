// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <utility/file/file_sys_util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <protocols/ligand_docking/ligand_dock_impl.hh>
#include <protocols/jobdist/JobDistributors.hh>

#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

#include <core/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //for addding constraints if demanded by user
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/moves/Mover.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>

#include <ctime>
#include <fstream>

//Auto Headers
#include <core/svn_version.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/grid/CartGrid.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/raw_data/DecoyFileData.hh>
#include <core/io/raw_data/ScoreFileData.hh>
#include <core/io/silent/EnergyNames.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <utility/stream_util.hh>
#include <basic/datacache/CacheableData.hh>
#include <protocols/checkpoint/Checkpoint.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/io/izstream.hh>
#include <utility/keys/Key2Tuple.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <time.h>



class LigandDockTest : public CxxTest::TestSuite {

public:

	void setUp() {
		// Stupid C++!  Without at least one string object, can't do concatenation.
		std::string empty = "";
		core_init_with_additional_options(empty
			+" -run:constant_seed"
			+" -run:rng mt19937"
			+" -in:file:s apps/public/ligand_docking/7cpa_7cpa_input.pdb"
			+" -in:file:native apps/public/ligand_docking/7cpa_7cpa_native.pdb"
			+" -out:path:pdb apps/public/ligand_docking"
			+" -packing:no_optH"
			// omitted for speed of execution, though it doesn't seem to matter much:
			//+" -packing:ex1"
			//+" -packing:ex1aro"
			//+" -packing:ex2"
			+" -docking:randomize2"
			+" -docking:uniform_trans 5"
			+" -docking:ligand:start_from  -1.731  32.589  -5.039"
			+" -docking:ligand:minimize_ligand"
			+" -docking:ligand:harmonic_torsions 10"
		);
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCAP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("ZN1")) params_files.push_back("apps/public/ligand_docking/ZN1.params");
		if(!residue_set.has_name("CP1")) params_files.push_back("apps/public/ligand_docking/7cpa.params");
		if(!residue_set.has_name("AQ1")) params_files.push_back("apps/public/ligand_docking/1aq1.params");
		residue_set.read_files(params_files,
			ChemicalManager::get_instance()->atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->element_set( FA_STANDARD ),
			ChemicalManager::get_instance()->mm_atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->orbital_type_set(FA_STANDARD));//,
			//ChemicalManager::get_instance()->csd_atom_type_set( FA_STANDARD ));
	}

	void tearDown() {}

	// This test is unreasonably sensitive to differences between hardware, leading to wildly different results.
	// As a result, this test is always broken. Comparing floats within some delta is not enough to fix the problem.
	// Therefore, I'm removing this test.  -IWD
	//
	//void test_ligand_docking() {
	//	using namespace utility::file;
	//	std::string silent_out = "apps/public/ligand_docking/silent.out";
	//	if( file_exists(silent_out) ) file_delete(silent_out);
	//	TS_ASSERT_EQUALS( 0, ligand_dock_main() );
	//	TS_ASSERT_FILE_EQ( "apps/public/ligand_docking/silent.out.ref", silent_out.c_str() )
	//}

	void test_select_best_poses() {
		core::import_pose::atom_tree_diffs::AtomTreeDiff atdiff("apps/public/ligand_docking/1dm2_1aq1_10poses.out");
		std::set< std::string > best_tags;

		protocols::ligand_docking::select_best_poses(atdiff, best_tags);
		TS_ASSERT_EQUALS(best_tags.size(), 1);
		TS_ASSERT_EQUALS(best_tags.count("flexbb0_5_1dm2_1aq1_0009"), 1);
	}
};

