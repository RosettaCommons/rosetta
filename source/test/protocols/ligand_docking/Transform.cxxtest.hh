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

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <core/pose/util.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/ligand_docking/Transform.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
//#include <test/core/init_util.hh>


static basic::Tracer TR("protocols.ligand_docking.StartFrom.cxxtest");


class TransformTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init_with_additional_options("-extra_res_fa protocols/ligand_docking/ZNx.params protocols/ligand_docking/7cpa.params");
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		//  using namespace core::chemical;
		//  utility::vector1< std::string > params_files;
		//  ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		//  ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		//  if ( !residue_set.has_name("ZNx") ) params_files.push_back("protocols/ligand_docking/ZNx.params");
		//  if ( !residue_set.has_name("CP1") ) params_files.push_back("protocols/ligand_docking/7cpa.params");
		//  residue_set.read_files_for_custom_residue_types(params_files);
	}

	void tearDown() {}

	//inline void core_init_with_additional_options( std::string const & commandline_in );

	void test_initial_perturb() {
		using namespace protocols::ligand_docking;

		numeric::random::rg().set_seed( "mt19937", time(0) );

		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "protocols/ligand_docking/7cpa_7cpa_native.pdb" , core::import_pose::PDB_file);

		Transform mover;

		core::Size chain_id=core::pose::get_chain_id_from_chain('X', pose);
		core::Size begin=pose.conformation().chain_begin(chain_id);

		core::conformation::UltraLightResidue test_ligand(pose.residue(begin).get_self_ptr());
		core::conformation::UltraLightResidue start_ligand(test_ligand);
		core::Vector start_center(start_ligand.center());

		//Test initial perturb of residue
		core::Size rejected = 0;
		core::Size accepted = 0;

		for ( core::Size i=0; i <= 500; i++ ) {
			mover.randomize_ligand(test_ligand, 5, 360);
			core::Real distance = test_ligand.center().distance(start_center);

			if ( distance > 5.0 ) {
				rejected++;
			} else {
				accepted++;
			}

			test_ligand = start_ligand;
		}

		TS_ASSERT_EQUALS(rejected, 0.0);

		//Test conformer change of residue
		rejected = 0;
		accepted = 0;
		core::Real deviation = 0;

		mover.setup_conformers(pose, begin);

		for ( core::Size i=0; i <= 500; i++ ) {
			mover.change_conformer(test_ligand);
			core::Real distance = test_ligand.center().distance(start_center);

			if ( distance > 5.0 ) {
				rejected++;
			} else {
				accepted++;
			}

			std::cout << "conformer distance: " << distance << std::endl;
			deviation = 0;

			utility::vector1<core::PointPosition > target_coords = start_ligand.coords_vector();
			utility::vector1<core::PointPosition > copy_coords = test_ligand.coords_vector();

			for ( core::Size i=1; i <= copy_coords.size(); ++i ) {
				core::Real deviation_x = ((copy_coords[i][0]-target_coords[i][0]) * (copy_coords[i][0]-target_coords[i][0]));
				core::Real deviation_y = ((copy_coords[i][1]-target_coords[i][1]) * (copy_coords[i][1]-target_coords[i][1]));
				core::Real deviation_z = ((copy_coords[i][2]-target_coords[i][2]) * (copy_coords[i][2]-target_coords[i][2]));

				core::Real total_dev = deviation_x + deviation_y + deviation_z;
				deviation += total_dev;
			}

			deviation /= (core::Real)copy_coords.size();
			deviation = sqrt(deviation);

			std::cout << "RMSD: " << deviation << std::endl;

			test_ligand = start_ligand;
		}

		TS_ASSERT_EQUALS(rejected, 0.0);



	}
};

