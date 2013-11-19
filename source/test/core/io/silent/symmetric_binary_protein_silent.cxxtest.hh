// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/io/silent/symmtric_protein_silent.cxxtest.hh
/// @brief  test suite for symmtric protein silent-file format
/// @author Lin Jiang

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <test/UTracer.hh>

// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// AUTO-REMOVED #include <numeric/random/random.hh>
#include <utility/file/file_sys_util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>
#include <iostream>

//Auto Headers
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("test.core.io.silent.symm_protein_silent");

using namespace core;

class Symmetric_Binary_Protein_Silent_File_Tests : public CxxTest::TestSuite
{
public:
	Symmetric_Binary_Protein_Silent_File_Tests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-symmetry::symmetry_definition core/scoring/symmetry/fibril_symm.dat -in::file::silent_struct_type binary -in::file::fullatom" );

		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCAP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if(!residue_set.has_name("GTP")) params_files.push_back("core/io/GTP.params");
		residue_set.read_files(params_files,
			ChemicalManager::get_instance()->atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->element_set( FA_STANDARD ),
			ChemicalManager::get_instance()->mm_atom_type_set( FA_STANDARD ),
			ChemicalManager::get_instance()->orbital_type_set(FA_STANDARD));
	}

	// Shared finalization goes here.
	void tearDown() {
	}


void test_save_and_restore()
{
  using namespace core::conformation::symmetry;

	double rms_threshold = 1e-2;
	double score_threshold = 1e-1;

	pose::Pose ref_pose, restored_pose;
	core::chemical::ResidueTypeSetCAP rsd =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::import_pose::pose_from_pdb( ref_pose, *rsd, std::string("core/scoring/symmetry/fibril_in.pdb"));
  pose::set_ss_from_phipsi( ref_pose );
	core::pose::symmetry::make_symmetric_pose( ref_pose );

  // Write out the silent-file
  core::io::silent::SilentFileData sfd;
  std::string silent_outfile = "core/io/sym_bin_silentfile_test.out";
  utility::file::file_delete( silent_outfile );
  core::io::silent::BinaryProteinSilentStruct pss,new_pss;
  pss.fill_struct( ref_pose, "ref_pose" );
  sfd.write_silent_struct( pss, silent_outfile );

	// Read the ProteinSilentStruct from the silent-file
  core::io::silent::SilentFileData sfe;
  utility::vector1 < std::string > tags;
  tags.push_back( "ref_pose" );
//	sfd.read_file( silent_outfile ); // read file w/ non-ideal geometry
	sfe.read_file( silent_outfile ); // read file w/ non-ideal geometry
	TS_ASSERT( sfe.size() > 0 );
	core::io::silent::SilentFileData::iterator iter = sfe.begin();
	iter->fill_pose( restored_pose, *rsd );

	// test symmetry info difference
	std::ostringstream si_ref, si_restored;
  si_ref << *core::pose::symmetry::symmetry_info( ref_pose );
  si_restored << *core::pose::symmetry::symmetry_info( restored_pose );
	if( si_ref.str() == si_restored.str() ) {
		TR << "SYMMETRY_INFO no difffernce" <<std::endl;
  }
	TS_ASSERT( si_ref.str() == si_restored.str() );


	// test rms difference
	Real rms_to_restored = scoring::CA_rmsd( ref_pose, restored_pose );
	TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
	TS_ASSERT( rms_to_restored < rms_threshold );

 // test score3 difference
/*  core::scoring::ScoreFunctionOP scorefxn3_sym =
        core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
  scorefxn3_sym = new core::scoring::symmetry::SymmetricScoreFunction(*scorefxn3_sym);
  Real score3_ref = (*scorefxn3_sym)(ref_pose);
  Real score3_restored = (*scorefxn3_sym)(restored_pose);
  Real score3_del = std::fabs( score3_restored - score3_ref );
  scorefxn3_sym->show( std::cout, ref_pose );
  scorefxn3_sym->show( std::cout, restored_pose );
  TR << "Score difference: " << score3_del << std::endl;
  TS_ASSERT( score3_del < score_threshold );
*/
	// test score13 difference
	core::scoring::ScoreFunctionOP scorefxn13_sym =
				core::scoring::ScoreFunctionFactory::create_score_function( "score13_env_hb" );
  scorefxn13_sym = core::scoring::symmetry::symmetrize_scorefunction(*scorefxn13_sym);
	Real score13_ref = (*scorefxn13_sym)(ref_pose);
	Real score13_restored = (*scorefxn13_sym)(restored_pose);
	Real score13_del = std::fabs( score13_restored - score13_ref );
	TR << "Score difference: " << score13_del << std::endl;
	TS_ASSERT( score13_del < score_threshold );

}

};
