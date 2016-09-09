// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// I think almost all of this is deprecated now. -- rhiju, 2014.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/Energies.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/types.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>

#include <basic/options/option_macros.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>

//RNA stuff.
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/movers/RNA_Minimizer.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/farna/setup/RNA_DeNovoPoseInitializer.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <protocols/farna/util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/RNA_LJ_BaseEnergy.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap

#include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <time.h>

using namespace core;
using namespace core::pose::rna;
using namespace basic;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options::OptionKeys;

using utility::vector1;

using ObjexxFCL::format::A;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using ObjexxFCL::string_of;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, create_vall_torsions )
OPT_KEY( Boolean, create_jump_database )
OPT_KEY( Boolean, icoord_test )
OPT_KEY( Boolean, print_internal_coord )
OPT_KEY( Boolean, extract )
OPT_KEY( Boolean, fullatom_score_test )
OPT_KEY( Boolean, fullatom_minimize )
OPT_KEY( Boolean, fullatom_multiscore )
OPT_KEY( Boolean, fullatom_minimize_silent)
OPT_KEY( Boolean, o2prime_test )
OPT_KEY( Boolean, skip_o2prime_pack )
OPT_KEY( Boolean, lores_score )
OPT_KEY( Boolean, lores_score_silent )
OPT_KEY( Boolean, env_test )
OPT_KEY( Boolean, rna_design )
OPT_KEY( Boolean, rna_design_gap )
OPT_KEY( Boolean, rna_idealize )
OPT_KEY( Boolean, idl_close_chainbreaks )
OPT_KEY( Boolean, rna_assemble )
OPT_KEY( Boolean, rna_stats )
OPT_KEY( Boolean, rna_torsion_check )
OPT_KEY( Boolean, sum_lores_plus_hires )
OPT_KEY( Boolean, fa_standard )
OPT_KEY( Boolean, close_chainbreaks_test)
OPT_KEY( Boolean, filter_base_pairs)
OPT_KEY( Boolean, filter_lores_base_pairs)
OPT_KEY( Boolean, rna_jumping )
OPT_KEY( Boolean, create_benchmark )
OPT_KEY( Boolean, minimize_rna )
OPT_KEY( Boolean, relax_rna )
OPT_KEY( Boolean, simple_relax )
OPT_KEY( Boolean, ignore_secstruct )
OPT_KEY( Boolean, chain_closure_test )
OPT_KEY( Boolean, backbone_rebuild_test )
OPT_KEY( Boolean, crazy_minimize )
OPT_KEY( Boolean, crazy_fold_tree )
OPT_KEY( Boolean, sasatest )
OPT_KEY( Boolean, close_loops )
OPT_KEY( Boolean, output_lores_silent_file )
OPT_KEY( Boolean, heat )
OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, convert_to_native )
OPT_KEY( Boolean, disable_o2prime_rotamers )
OPT_KEY( Boolean, disable_include_current )
OPT_KEY( Boolean, sample_chi )
OPT_KEY( Boolean, pymol_struct_type )
OPT_KEY( Boolean, print_hbonds )
OPT_KEY( Boolean, calc_rmsd )
OPT_KEY( Boolean, dinucleotide )
OPT_KEY( Boolean, vary_geometry )
OPT_KEY( Boolean, more_rotamers )
OPT_KEY( Boolean, quick_test )
OPT_KEY( Boolean, rotamerize_test )
OPT_KEY( Boolean, build_next_nucleotide )
OPT_KEY( Boolean, prepend_residue )
OPT_KEY( Boolean, sugar_geometry_test )
OPT_KEY( Boolean, sugar_frag_test )
OPT_KEY( Boolean, color_by_geom_sol )
OPT_KEY( Boolean, color_by_rna_lj_base )
OPT_KEY( Boolean, files_for_openMM )
OPT_KEY( Boolean, print_torsions )
OPT_KEY( Boolean, print_secstruct )
OPT_KEY( Boolean, skip_coord_constraints )
OPT_KEY( Real, fa_stack_weight )
OPT_KEY( Real, temperature )
OPT_KEY( Real, jump_change_frequency )
OPT_KEY( Real, rmsd_cutoff )
OPT_KEY( Integer, cycles )
OPT_KEY( String,  vall_torsions )
OPT_KEY( String,  assemble_file )
OPT_KEY( String,  basepair_file )
OPT_KEY( String,  jump_library_file )
OPT_KEY( String,  params_file )
OPT_KEY( String,  data_file )
OPT_KEY( String,  cst_file )
OPT_KEY( String,  rsd_type_set )
OPT_KEY( Real, atom_pair_constraint_weight )
OPT_KEY( Real, coordinate_constraint_weight )


// ///////////////////////////////////////////////////////////////////////////////
// void
// rna_basepair_jump_atoms( chemical::AA const res,
//              Size & atom1, Size & atom2, Size & atom3){
//  // using namespace param_aa;
//  // assert( is_RNA( res ));
//  using namespace chemical;

//  if ( res == na_rad ){
//   atom1 = 18; // C6
//   atom2 = 13; // N1
//   atom3 = 14; // C2
//  }

//  if ( res == na_rcy ){
//   atom1 = 17; // C4
//   atom2 = 16; // N3
//   atom3 = 14; // C2
//  }

//  if ( res == na_rgu ){
//   atom1 = 19; // C6
//   atom2 = 13; // N1
//   atom3 = 14; // C2
//  }

//  if ( res == na_ura ){
//   atom1 = 17; // C4
//   atom2 = 16; // N3
//   atom3 = 14; // C2
//  }
// }

// //////////////////////////////////////////
// // COPY FROM PHIL'S DEVEL CODE
// // DOESN'T WORK FOR RNA?
// void
// set_na_jump_atoms( pose::Pose & pose )
// {
//  using conformation::Residue;
//  using namespace id;

//  kinematics::FoldTree f( pose.fold_tree() );

//  // anchor intra-dna jumps at fourth chi1 atom (out in the base)
//  //
//  for ( Size i=1; i<= f.num_jump(); ++i ) {
//   Residue const & rsd1( pose.residue( f.  upstream_jump_residue( i ) ) );
//   Residue const & rsd2( pose.residue( f.downstream_jump_residue( i ) ) );
//   if ( rsd1.is_NA() && rsd2.is_NA() ) {
//    f.set_jump_atoms( i, rsd1.atom_name( rsd1.chi_atoms(1)[4] ), rsd2.atom_name( rsd2.chi_atoms(1)[4] ) );
//   } else if ( rsd1.is_NA() && rsd2.is_protein() ) {
//    f.set_jump_atoms( i, rsd1.atom_name( rsd1.chi_atoms(1)[4] ),  "CA" );
//   } else if ( rsd2.is_NA() && rsd1.is_protein() ) {
//    f.set_jump_atoms( i, "CA", rsd2.atom_name( rsd2.chi_atoms(1)[4] ) );
//   }
//  }

//  pose.fold_tree( f );

//  // tinker with atomorder in atomtree so that we'll get the jump stubs we want
//  // useful for graphics -- keeps root base fixed in space
//  //
//  // I'm not actually sure that this is so important anymore
//  // I mean, it seems likely that this will in fact be the stub... but can't remember the details...
//  //
//  for ( Size i=1; i<= pose.size(); ++i ) {
//   Residue const & rsd( pose.residue(i) );
//   if ( rsd.is_NA() ) {
//    pose.conformation().set_jump_atom_stub_id( StubID ( AtomID( rsd.chi_atoms(1)[4], i ),
//                              AtomID( rsd.chi_atoms(1)[3], i ),
//                              AtomID( rsd.chi_atoms(1)[2], i ) ) );
//   }
//  }


// }


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
figure_out_icoord_test( ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Create an extended pose from scratch.
	pose::Pose extended_pose;
	std::string sequence = "ggguuu";

	std::cout << "ABOUT TO MAKE EXTENDED POSE" << std::endl;

	//  ResidueTypeCOPs const & rsd_type_list( rsd_set->aa_map( aa_from_oneletter_code( 'a' )  ) );
	//  std::cout << "SIZE " << rsd_type_list.size() << std::endl;
	//  for ( Size j=1; j<= rsd_type_list.size(); ++j ) {
	//   ResidueType const & rsd_type( *(rsd_type_list[j]) );
	//   std::cout << "ISPOLYMER" << j << " " << rsd_type.is_polymer() << std::endl;
	//  }

	make_pose_from_sequence(
		extended_pose,
		sequence,
		*rsd_set );

	dump_pdb( extended_pose, "extended.pdb" );

	std::string infile  = option[ in ::file::s ][1];

	pose::Pose pose,start_pose;
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);
	/////////////////////////////////////////
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
	/////////////////////////////////////////

	std::cout << "READ POSE FROM PDB" << std::endl;

	dump_pdb( pose, "test.pdb" );

	std::cout << "DUMPED PDB" << std::endl;

	int const res_num = 1;

	start_pose = pose;

	//Hmm, do I understand how to address the atom tree degrees of freedom?
	for ( Size j=1; j <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
		pose = start_pose;
		id::TorsionID my_ID( res_num /* res num*/, id::BB, j );
		id::AtomID id1,id2,id3,id4;
		bool fail = pose.conformation().get_torsion_angle_atom_ids( my_ID, id1, id2, id3, id4 );
		if ( !fail ) {
			std::cout << " BB: " << j << "  " <<
				id1.rsd()  << ' ' << pose.residue_type( id1.rsd() ).atom_name( id1.atomno() ) << "   " <<
				id2.rsd()  << ' ' << pose.residue_type( id2.rsd() ).atom_name( id2.atomno() ) << "   " <<
				id3.rsd()  << ' ' << pose.residue_type( id3.rsd() ).atom_name( id3.atomno() ) << "   " <<
				id4.rsd()  << ' ' << pose.residue_type( id4.rsd() ).atom_name( id4.atomno() ) << " : " <<
				pose.torsion( my_ID ) <<
				std::endl;

			pose.set_torsion( my_ID, 180.0 );

			dump_pdb( pose, "testBB"+string_of(j)+".pdb" );
			std::cout << "DUMPED PDB" << j << std::endl;
		}
	}


	for ( Size j=1; j <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {

		pose = start_pose;

		id::TorsionID my_ID( res_num /* res num*/, id::CHI, j );

		pose.set_torsion( my_ID, 180.0 );

		id::AtomID id1,id2,id3,id4;
		pose.conformation().get_torsion_angle_atom_ids( my_ID, id1, id2, id3, id4 );
		std::cout << "CHI: " << j << "  " <<
			id1.rsd()  << ' ' << pose.residue_type( id1.rsd() ).atom_name( id1.atomno() ) << "   " <<
			id2.rsd()  << ' ' << pose.residue_type( id2.rsd() ).atom_name( id2.atomno() ) << "   " <<
			id3.rsd()  << ' ' << pose.residue_type( id3.rsd() ).atom_name( id3.atomno() ) << "   " <<
			id4.rsd()  << ' ' << pose.residue_type( id4.rsd() ).atom_name( id4.atomno() ) <<  " : " <<
			pose.torsion( my_ID ) <<
			std::endl;

		dump_pdb( pose, "testCHI"+string_of(j)+".pdb" );
		std::cout << "DUMPED PDB" << j << std::endl;
	}


}


///////////////////////////////////////////////////////////////////////////////
void
rna_fullatom_score_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::io::silent;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	if ( option[fa_standard] ) rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	for ( Size i = 1; i <= pdb_files.size(); i++ ) {
		std::string const pdb_file = pdb_files[i];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		//import_pose::pose_from_file( pose, pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		std::cout << "Check it! SEQUENCE " << pose.sequence() << std::endl;

		if ( option[ data_file].user() ) {
			core::io::rna::RNA_DataReader rna_data_reader( option[ data_file ] );
			rna_data_reader.fill_rna_data_info( pose );
		}


		//Score these suckers.
		ScoreFunctionOP scorefxn = get_score_function();

		(*scorefxn)(pose);
		scorefxn->show( std::cout, pose );

		RNA_SilentStruct s( pose, pdb_file );
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );

		//Testing HB env dep...
		if ( false ) {
			//    std::cout << "*****************************" << std::endl;
			//    std::cout << "TURNING OFF HB ENV DEP" << std::endl;
			//    std::cout << "*****************************" << std::endl;
			//    scorefxn->energy_method_options().hbond_options()->use_hb_env_dep( false );
			//    (*scorefxn)( pose );
			//    scorefxn->show( std::cout, pose );

			//    std::cout << "*****************************" << std::endl;
			//    std::cout << "TURNING ON HB ENV DEP" << std::endl;
			//    std::cout << "*****************************" << std::endl;
			//    scorefxn->energy_method_options().hbond_options()->use_hb_env_dep( true );
			//    (*scorefxn)( pose );
			//    scorefxn->show( std::cout, pose );

			//    std::cout << "*****************************" << std::endl;
			//    std::cout << "TURNING OFF HB ENV DEP" << std::endl;
			//    std::cout << "*****************************" << std::endl;
			//    scorefxn->energy_method_options().hbond_options()->use_hb_env_dep( false );
			//    (*scorefxn)( pose );
			//    scorefxn->show( std::cout, pose );

			using namespace core::scoring::methods;
			std::cout << "*****************************" << std::endl;
			std::cout << "TURNING ON HB ENV DEP2" << std::endl;
			std::cout << "*****************************" << std::endl;
			EnergyMethodOptions options( scorefxn->energy_method_options() );
			options.hbond_options().use_hb_env_dep( true );
			scorefxn->set_energy_method_options( options );
			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );

			std::cout << "*****************************" << std::endl;
			std::cout << "TURNING OFF HB ENV DEP2" << std::endl;
			std::cout << "*****************************" << std::endl;
			options.hbond_options().use_hb_env_dep( false );
			scorefxn->set_energy_method_options( options );
			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );

			//Rescore
			std::cout << "*****************************" << std::endl;
			std::cout << " RESCORE " << std::endl;
			std::cout << "*****************************" << std::endl;
			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );

			using namespace core::scoring::methods;
			std::cout << "*****************************" << std::endl;
			std::cout << "TURNING ON HB ENV DEP2" << std::endl;
			std::cout << "*****************************" << std::endl;
			options.hbond_options().use_hb_env_dep( true );
			scorefxn->set_energy_method_options( options );
			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );

			std::cout << "*****************************" << std::endl;
			std::cout << "TURNING OFF HB ENV DEP2" << std::endl;
			std::cout << "*****************************" << std::endl;
			options.hbond_options().use_hb_env_dep( false );
			scorefxn->set_energy_method_options( options );
			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );

		}

		//  pose.dump_pdb( "score.pdb" );

		//   for (Size i = 1; i <= pose.size(); i++ ){
		//    id::TorsionID torsion_id( i, id::CHI, 4 );
		//    id::DOF_ID const dof_id(  pose.conformation().dof_id_from_torsion_id( torsion_id ) );
		//    std::cout << i << " 2'-OH-TORSION: " <<  pose.torsion( torsion_id ) << " " << numeric::conversions::degrees( pose.conformation().dof( dof_id ) ) << std::endl;
		//   }

		//////////////////////////////////
		//////////////////////////////////
		//test
		//    for (Real x = -180.0; x <= 180.0; x++ ){
		//      pose::Pose my_pose;
		//      my_pose = pose;
		//      pose.set_torsion( id::TorsionID( 2, id::CHI, 4), x );
		//      Real my_score = (*scorefxn)( my_pose );
		//      RNA_SilentStruct s( my_pose, string_of( x ) );
		//      silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
		//     }
		//////////////////////////////////
		//////////////////////////////////


	}

}

//////////////////////////////////////////////////////////////////////////////
void
add_coordinate_constraints( pose::Pose & pose ) {
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	ConstraintSetOP cst_set( new ConstraintSet() );

	Real const coord_sdev( 2.0 );
	Size const my_anchor( 1 ); //anchor atom on first residue?

	Size const nres( pose.size() );
	for ( Size i=1; i<= nres;  ++i ) {

		Residue const & i_rsd( pose.residue(i) );

		for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
			core::scoring::func::FuncOP fx( new HarmonicFunc( 0.0, coord_sdev ) );
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( AtomID(ii,i), AtomID(1,my_anchor), i_rsd.xyz(ii), fx ) ) ) );
		}
	}

	pose.constraint_set( cst_set );
}

///////////////////////////////////////////////////////////////////////////////
void
pack_o2prime( core::pose::Pose & pose, core::scoring::ScoreFunction const & scorefxn ) {

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !pose.residue_type(i).is_RNA() ) continue;
		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		task->nonconst_residue_task(i).or_ex4( true );
		task->nonconst_residue_task(i).or_include_current( true );
	}

	std::cout << "Packing 2' hydroxyls..." << std::endl;
	pack::pack_rotamers( pose, scorefxn, task);
}

//////////////////////////////////////////////////////////////////////////////
//Ahem, copied from phil.
void
setup_rna_chainbreak_constraints(
	pose::Pose & pose
) {
	using namespace scoring::constraints;
	using namespace scoring::func;
	using namespace chemical;
	using namespace conformation;
	using namespace id;
	using numeric::conversions::radians;

	Real const O3_P_distance( 1.608 );
	Real const O3_angle( 119.8 );
	Real const  P_angle( 103.4 );

	Real const distance_stddev( 0.3 ); // amber is 0.0659
	Real const angle_stddev_degrees( 35 ); // amber is 8.54 (P angle), 5.73 (O3 angle)

	ConstraintSetOP cst_set( pose.constraint_set()->clone() );
	assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();

	FuncOP const distance_func( new HarmonicFunc( O3_P_distance, distance_stddev ) );
	FuncOP const O3_angle_func( new HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees ) ) );
	FuncOP const  P_angle_func( new HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees ) ) );

	for ( Size i=1; i< pose.size(); ++i ) {
		ResidueType const & rsd1( pose.residue_type( i   ) );
		ResidueType const & rsd2( pose.residue_type( i+1 ) );
		if ( rsd1.is_NA() && !rsd1.is_upper_terminus() && rsd2.is_NA() && !rsd2.is_lower_terminus() ) {
			//   tt << "adding dna chainbreak constraint between residues " << i << " and " << i+1 << std::endl;

			AtomID const C3_id( rsd1.atom_index( "C3'" ), i   );
			AtomID const O3_id( rsd1.atom_index( "O3'" ), i   );
			AtomID const  P_id( rsd2.atom_index( "P"   ), i+1 );
			AtomID const O5_id( rsd2.atom_index( "O5'" ), i+1 );

			// distance from O3' to P
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( O3_id, P_id, distance_func ) ) ) );

			// angle at O3'
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) ) ) );

			// angle at P
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) ) ) );
		}
	}

	pose.constraint_set( cst_set );

}

//////////////////////////////////////////////////////////////////////////////
// Probably belongs in util.hh
void
get_basepair_atoms( pose::Pose & pose,
	Size const & i, Size const & j,
	std::string & atom1, std::string & atom2 ) {

	using namespace core::chemical;
	char aa1 = pose.residue(i).name1();
	char aa2 = pose.residue(j).name1();

	if ( aa1=='a' && aa2=='u' ) {
		atom1 = " N1 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1=='u' && aa2=='a' ) {
		atom1 = " N3 ";
		atom2 = " N1 ";
		return;
	}
	if ( aa1=='g' && aa2=='c' ) {
		atom1 = " N1 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1=='c' && aa2=='g' ) {
		atom1 = " N3 ";
		atom2 = " N1 ";
		return;
	}
	if ( aa1=='g' && aa2=='u' ) {
		atom1 = " O6 ";
		atom2 = " N3 ";
		return;
	}
	if ( aa1=='u' && aa2=='g' ) {
		atom1 = " N3 ";
		atom2 = " O6 ";
		return;
	}

	//SPECIAL HACK -- force a stack.
	if ( aa1=='a' && aa2=='g' ) {
		std::cout << "BIG PROBLEM!" << std::endl;
		std::cout << "CRAZY STACK HACK!" << std::endl;
		atom1 = " C2'";
		atom2 = " O4'";
	}

	// assert( 1+1 == 3);
	return;

}


//////////////////////////////////////////////////////////////////////////////
void
setup_rna_base_pair_constraints( pose::Pose & pose ){

	using namespace scoring::constraints;
	using namespace scoring::func;
	using namespace conformation;
	using namespace options;
	using namespace id;

	std::string const basepair_filename = option[ basepair_file ];
	std::cout << "READING IN BASE PAIR FILE: " << basepair_filename << std::endl;
	utility::io::izstream data_stream( basepair_filename );

	ConstraintSetOP cst_set( pose.constraint_set()->clone() );
	assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();

	std::string line;

	utility::vector1 < std::pair< Size,Size> > basepair_list;

	Size pos1, pos2;
	while (  getline(data_stream, line) ) {
		std::istringstream line_stream( line );
		line_stream >> pos1;

		while ( !line_stream.fail() ) {
			line_stream >> pos2;
			if ( !line_stream.fail() )  {
				basepair_list.push_back( std::make_pair( pos1, pos2 ) );
			}
		}
	}

	Real const WC_distance( 2.9 );
	Real const distance_stddev( 0.25 );
	FuncOP const distance_func( new HarmonicFunc( WC_distance, distance_stddev ) );

	std::string atom1, atom2;
	for ( Size n = 1; n <= basepair_list.size(); n++ ) {
		Size const i = basepair_list[n].first;
		Size const j = basepair_list[n].second;
		get_basepair_atoms( pose, i, j,  atom1, atom2 );

		Distance const dist = (pose.residue(i).xyz( atom1 ) - pose.residue(j).xyz( atom2 ) ).length();
		std::cout << "Adding BASEPAIR constraint: " << pose.residue(i).name1() << I(3,i) << " <-->  " <<
			pose.residue(j).name1() << I(3,j) << "   " << atom1 << " <--> " << atom2 << ".  Current distance: " << dist << std::endl;

		Size const atom1_index = pose.residue(i).atom_index( atom1 );
		Size const atom2_index = pose.residue(j).atom_index( atom2 );
		cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( AtomID( atom1_index, i ), AtomID( atom2_index, j), distance_func ) ) ) );
	}

	pose.constraint_set( cst_set );

}

///////////////////////////////////////////////////////////////////////////////
void
rna_fullatom_multiscore_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::id;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	if ( option[fa_standard] ) rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	pose::Pose pose;
	std::string const pdb_file = pdb_files[1];

	import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
	/////////////////////////////////////////
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
	/////////////////////////////////////////

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Pose const pose_start = pose;

	ScoreFunctionOP scorefxn = get_score_function();

	clock_t const time_start( clock() );

	Size numscores = 10000;
	if ( option[ out::nstruct ].user() ) numscores = option[ out::nstruct ]();


	kinematics::MoveMap mm;
	mm.set_bb(  true );
	mm.set_chi( true );
	mm.set_jump( true );
	MinimizerMap min_map;
	min_map.setup( pose, mm );
	scorefxn->setup_for_minimizing( pose, min_map );
	// setup the function that we will pass to the low-level minimizer
	AtomTreeMultifunc f( pose, min_map, *scorefxn, false, false );
	Multivec dofs( min_map.nangles() );

	for ( Size count = 1; count <= numscores; count++ ) {

		if ( count % 1000 == 0 ) std::cout << "Done with score number: " << count << std::endl;

		for ( Size i = 1; i <= pose.size(); i++ ) {
			for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {
				pose.set_xyz( AtomID(j,i), pose_start.xyz( AtomID(j,i) ) + numeric::xyzVector<Real>( 0.0, 0.0, 0.001) * count  );
			}
		}

		//( *scorefxn )( pose );

		min_map.copy_dofs_from_pose( pose, dofs );
		Multivec G( dofs.size() );
		f.dfunc(  dofs,  G);
	}


	std::cout << "Total time for " << numscores << " score calls: " <<
		static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;


}


///////////////////////////////////////////////////////////////////////////////
void
convert_to_native_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	if ( option[fa_standard] ) rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	for ( Size i = 1; i <= pdb_files.size(); i++ ) {
		std::string const pdb_file = pdb_files[i];

		pose::Pose pose, start_pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		//Yes this is ridiculous
		std::string const output_file = "native_"+pdb_file;
		std::cout << "Outputting ==> " << output_file << std::endl;
		dump_pdb( pose, output_file );

	}

}


///////////////////////////////////////////////////////////////////////////////
void
rna_fullatom_minimize_silent_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	std::string const infile  = option[ in::file::silent  ][1];
	std::string const outfile = option[ out::file::silent  ]();

	// Silent file setup?
	SilentFileData silent_file_data_in, silent_file_data_out;
	silent_file_data_in.read_file( infile );

	pose::Pose native_pose;
	std::string native_pdb_file  = option[ in::file::native ];
	import_pose::pose_from_file( native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);
	/////////////////////////////////////////
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( native_pose );
	/////////////////////////////////////////

	utility::vector1< std::string > tags_done;
	std::map< std::string, bool > tag_is_done;
	if ( utility::file::file_exists( outfile ) ) {
		tags_done = silent_file_data_in.read_tags_fast( outfile );
		for ( utility::vector1< std::string >::const_iterator iter = tags_done.begin(); iter != tags_done.end(); iter++ ) {
			std::cout << "Already done? " << *iter << std::endl;
			tag_is_done[ *iter ] = true;
		}
	}

	pose::Pose ideal_pose;
	bool const use_input_pose = option[ in::file::s ].active();
	if ( use_input_pose ) {
		std::string ideal_pdb_file  = option[ in::file::s ][1];
		import_pose::pose_from_file( ideal_pose, *rsd_set, ideal_pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( ideal_pose );
		/////////////////////////////////////////
	}

	//Would be nice to set up some operators to "add" and "multiply" weight sets.
	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );

	for ( core::io::silent::SilentFileData::iterator iter = silent_file_data_in.begin(), end = silent_file_data_in.end(); iter != end; ++iter ) {

		std::string const tag = iter->decoy_tag();
		std::string const out_tag =  "minimize_" + tag;
		std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
		std::cout << "Let us do: " << tag << ". Already done? " << tag_is_done[ out_tag  ] << std::endl;
		if ( tag_is_done[ out_tag ] ) continue;

		pose::Pose pose = ideal_pose;
		iter->fill_pose( pose, *rsd_set /*, use_input_pose*/ );

		protocols::farna::movers::RNA_Minimizer rna_minimizer;
		rna_minimizer.apply( pose );

		RNA_SilentStruct s( pose, out_tag );

		Real const rmsd = all_atom_rmsd( native_pose, pose );
		std::cout << "All atom rmsd: " << rmsd  << std::endl;
		s.add_energy( "rms", rmsd );

		//finish up -- calculate lo res terms, and save values in silent struct.
		pose::Pose pose_for_lores_scoring  = pose;
		(*lores_scorefxn)( pose_for_lores_scoring );
		for ( Size i = 1; i <= n_score_types; i++ ) {
			ScoreType const t = ScoreType( i );
			if ( lores_scorefxn->get_weight( t ) != 0.0 ) {
				s.add_energy( name_from_score_type( t ), pose_for_lores_scoring.energies().total_energies()[ t ] );
			}
		}

		std::cout << "Outputting " << out_tag << " to silent file: " << outfile << std::endl;
		silent_file_data_out.write_silent_struct( s, outfile, false /*write score only*/ );


		//  dump_pdb( pose, out_tag+".pdb" );

	}

}
///////////////////////////////////////////////////////////////////////////////
void
rna_o2prime_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	std::string const silent_file = option[ out::file::silent  ]();
	// Silent file setup?
	SilentFileData silent_file_data;

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	for ( Size i = 1; i <= pdb_files.size(); i++ ) {
		std::string const pdb_file = pdb_files[i];

		pose::Pose pose, start_pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////
		start_pose = pose;

		std::cout << "Check it! SEQUENCE " << pose.sequence() << std::endl;

		ScoreFunctionOP scorefxn = get_score_function();

		scorefxn->show( std::cout, pose );

		{ //output initial structure.
			RNA_SilentStruct s( pose, pdb_file );
			std::cout << "Outputting " << pdb_file << " to silent file: " << silent_file << std::endl;
			silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
		}

		//2'-OH rotamer trials?!
		pack_o2prime( pose, *scorefxn );

		scorefxn->show( std::cout, pose );

		std::string const out_file =  "zzz_"+pdb_file;
		RNA_SilentStruct s( pose, out_file );
		std::cout << "Outputting " << out_file << " to silent file: " << silent_file << std::endl;
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

		dump_pdb( pose, out_file );

	}

}

///////////////////////////////////////////////////////////////////////////////
void
rna_lores_score_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );


	//Score these suckers.
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	for ( Size i = 1; i <= pdb_files.size(); i++ ) {
		std::string const pdb_file = pdb_files[i];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		std::cout << "Check it! SEQUENCE " << pose.sequence() << std::endl;

		scorefxn->show( std::cout, pose );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
rna_lores_score_silent_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::io::silent;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	//Score these suckers.
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );

	//Need to refactor silent file -- why don't we give it a filename here?
	SilentFileData silent_file_data_in, silent_file_data_out;
	std::string const infile  = option[ in::file::silent  ][1];
	std::string const outfile = option[ out::file::silent  ]();
	silent_file_data_in.read_file( infile );

	//Native pose setup
	pose::Pose native_pose;
	bool const use_native = option[in::file::native].active();
	if ( use_native )  {
		import_pose::pose_from_file( native_pose, *rsd_set, option( in::file::native ) , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( native_pose );
		/////////////////////////////////////////
	}
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	pose::Pose ideal_pose;
	bool const use_input_pose = option[ in::file::s ].active();
	if ( use_input_pose ) {
		std::string ideal_pdb_file  = option[ in::file::s ][1];
		import_pose::pose_from_file( ideal_pose, *rsd_set, ideal_pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( ideal_pose );
		/////////////////////////////////////////
	}


	for ( core::io::silent::SilentFileData::iterator iter = silent_file_data_in.begin(), end = silent_file_data_in.end(); iter != end; ++iter ) {

		std::string const tag = iter->decoy_tag();
		std::string const out_tag =  tag;

		pose::Pose pose = ideal_pose;
		iter->fill_pose( pose, *rsd_set /*, use_input_pose*/ );

		std::cout << "Check it! SEQUENCE " << pose.sequence() << std::endl;

		scorefxn->show( std::cout, pose );

		RNA_SilentStruct s( pose, out_tag );
		if ( use_native ) {
			Real const rmsd = all_atom_rmsd( native_pose, pose );
			std::cout << "All atom rmsd: " << rmsd  << std::endl;
			s.add_energy( "rms", rmsd );
		}
		std::cout << "Outputting " << out_tag << " to silent file: " << outfile << std::endl;
		silent_file_data_out.write_silent_struct( s, outfile, true /*write score only*/ );

	}

}


///////////////////////////////////////////////////////////////////////////////
void
output_struct_type(
	utility::io::ozstream & pymol_out,
	pose::Pose const & pose,
	ObjexxFCL::FArray1D_int & struct_type,
	int const & match_type)
{

	bool start_list( false );

	Size const nres( pose.size() );

	for ( Size i = 1; i <= nres; i++ )  {
		if ( struct_type(i) == match_type ) {
			if ( start_list ) pymol_out << '+';
			pymol_out << pose.pdb_info()->number(i);
			start_list = true;
		}
	}

	pymol_out << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
pymol_struct_type_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	pose::Pose pose;
	std::string infile  = option[ in::file::s ][1];
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);

	Size const nres = pose.size();
	FArray1D_int struct_type( nres, -1 );
	protocols::farna::check_base_pair( pose, struct_type );

	utility::io::ozstream pymol_out( infile+".pml" );

	pymol_out << "reinitialize "  << std::endl;
	pymol_out << "load " << infile << ",native  "  << std::endl;
	pymol_out << "set cartoon_ring_mode, 1" << std::endl;
	pymol_out << "cartoon oval" << std::endl;
	pymol_out << "set cartoon_oval_width, 0.4" << std::endl;
	pymol_out << "set cartoon_oval_length, 0.8" << std::endl;
	pymol_out << "show cartoon, native " << std::endl;
	pymol_out << "color white, native " << std::endl;
	pymol_out << "bg_color white" << std::endl;


	// bool start_list;
	pymol_out << "select SS, resi ";
	output_struct_type( pymol_out, pose, struct_type, 0);
	pymol_out << "color red, SS " << std::endl;

	pymol_out << "select DS, resi ";
	output_struct_type( pymol_out, pose, struct_type, 1);
	pymol_out << "color gray50, DS " << std::endl;

	pymol_out << "select TS, resi ";
	output_struct_type( pymol_out, pose, struct_type, 2);
	pymol_out << "color red, SS " << std::endl;

	pymol_out.close();

}

///////////////////////////////////////////////////////////////////////////////
void
rna_design_gap_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::methods;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	Real const design_temp( option[ temperature ]() );

	ScoreFunctionOP scorefxn = get_score_function();

	// scorefxn->energy_method_options().exclude_DNA_DNA( exclude_DNA_DNA );
	EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );

	std::string const in_path = option[ in::path::path ]()[1];

	for ( Size i = 1; i <= pdb_files.size(); i++ ) {
		std::string const pdb_file = pdb_files[i];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, in_path + pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			//Hmmm, extras.
			task->nonconst_residue_task( ii ).and_extrachi_cutoff( 0 );
			task->nonconst_residue_task( ii ).or_ex4( true );
			task->nonconst_residue_task( ii ).or_include_current( true );
		}

		// Special -- set low temp by hand
		task->low_temp( design_temp );
		// Special -- disallow quench.
		task->disallow_quench( false );

		pack::task::PackerTaskOP task_design = task->clone();
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			task_design->nonconst_residue_task( ii ).allow_aa( na_rad );
			task_design->nonconst_residue_task( ii ).allow_aa( na_ura );
			task_design->nonconst_residue_task( ii ).allow_aa( na_rgu );
			task_design->nonconst_residue_task( ii ).allow_aa( na_rcy );
			assert( task->design_residue(ii) );
		}

		Size pos( pdb_file.find( ".pdb" ) );

		// Special -- do ten repacks
		Size const nstruct = option[ out::nstruct ];
		utility::vector1< std::pair< Real, std::string > > results;
		utility::vector1< pose::PoseOP > pose_list;
		pack::pack_rotamers_loop( pose, *scorefxn, task, nstruct, results, pose_list);
		{
			std::string outfile( pdb_file );
			outfile.replace( pos, 4, ".pack.txt" );
			protocols::farna::export_packer_results( results, pose_list, scorefxn, outfile );
		}

		// Do ten designs.
		results.clear();
		pose_list.clear();
		pack::pack_rotamers_loop( pose, *scorefxn, task_design, nstruct, results, pose_list);
		{
			std::string outfile( pdb_file );
			outfile.replace( pos, 4, ".design.txt" );
			protocols::farna::export_packer_results( results, pose_list, scorefxn, outfile );
		}
	}

}


///////////////////////////////////////////////////////////////////////////////
void
print_internal_coord_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::id;
	using numeric::conversions::degrees;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( option[ rsd_type_set]()  );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
	/////////////////////////////////////////
	if ( option[ rsd_type_set]() == core::chemical::FA_STANDARD ) {
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );
	}

	kinematics::FoldTree f( pose.size() );
	Size start( 2 ), end( 6 ), cutpos( 3 );
	// Size start( 8 ), end( 9 );
	//  Size start( 2 ), end( 5 );
	f.new_jump( start, end, cutpos );


	//  f.set_jump_atoms( 1,
	//           core::chemical::rna::chi1_torsion_atom( pose.residue( start ) ),
	//           core::chemical::rna::chi1_torsion_atom( pose.residue( end ) ) );

	//  f.set_jump_atoms( 1,
	//           core::chemical::rna::chi1_torsion_atom( pose.residue( end ) ),
	//           " O2'" );

	f.set_jump_atoms( 1, " Y  ", " Y  " );
	pose.fold_tree( f );
	//  set_na_jump_atoms( pose );


	pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
	pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

	pose.set_xyz( id::NamedAtomID( "OVL1", cutpos  ),  pose.residue( cutpos+1 ).xyz( " P  " ) );
	pose.set_xyz( id::NamedAtomID( "OVL2", cutpos  ),  pose.residue( cutpos+1 ).xyz( " S  " ) );
	pose.set_xyz( id::NamedAtomID( "OVU1", cutpos+1),  pose.residue( cutpos   ).xyz( " S  " ) );


	//Special trick to figure out where to put VO4' virtual atom for sugar ring closure.
	//  for (Size i= 1; i <= pose.size(); i++ ) {
	//   conformation::Residue const & rsd( pose.residue( i ) ) ;
	//   for (Size j = 1; j <= rsd.natoms(); j++ ) {
	//    if ( rsd.atom_name( j ) == "VO4'" ) {
	//     pose.set_xyz( id::AtomID(j,i ), rsd.xyz( "O4'") );
	//    }
	//   }
	//  }

	protocols::farna::print_internal_coords( pose );


	if ( option[ rsd_type_set ]() == "coarse_rna" ) {
		{
			conformation::Residue const & rsd = pose.residue(cutpos);
			if ( rsd.has( "OVL1" ) ) {
				std::cout << "OVL1 torsion: " << numeric::conversions::degrees(  numeric::dihedral_radians( rsd.xyz( "OVL1" ), rsd.xyz( "S" ), rsd.xyz( "P" ), rsd.xyz( "CEN" ) ) ) << std::endl;
				std::cout << "OVL2 torsion: " << numeric::conversions::degrees(  numeric::dihedral_radians( rsd.xyz( "OVL2" ), rsd.xyz( "OVL1" ), rsd.xyz( "S" ), rsd.xyz( "P" ) ) ) << std::endl;
			}
		}

		{
			conformation::Residue const & rsd = pose.residue(cutpos+1);
			if ( rsd.has( "OVU1" ) ) {
				std::cout << "OVU1 torsion: " << numeric::conversions::degrees(  numeric::dihedral_radians( rsd.xyz( "OVU1" ), rsd.xyz( "P" ), rsd.xyz( "S" ), rsd.xyz( "CEN" ) ) ) << std::endl;
			}
		}
	}


}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
set_ideal_geometry( pose::Pose & pose, pose::Pose const & extended_pose, chemical::ResidueTypeSetCOP & rsd_set ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::s ].active() ) {

		pose::Pose ideal_pose;
		std::string ideal_pdb_file  = option[ in::file::s ][1];
		import_pose::pose_from_file( ideal_pose, *rsd_set, ideal_pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( ideal_pose );
		/////////////////////////////////////////

		pose = ideal_pose;
		for ( Size i = 1; i <= pose.size(); i++ ) {
			for ( Size j=1; j <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
				id::TorsionID my_ID( i, id::BB, j );
				pose.set_torsion( my_ID, extended_pose.torsion( my_ID ) );
			}

			//Is it really necessary to hard code 6 backbone torsion angles and 4 chi angles?
			for ( Size j=1; j <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {
				id::TorsionID my_ID( i, id::CHI, j );
				pose.set_torsion( my_ID, extended_pose.torsion( my_ID ) );
			}
		}

	}

}


//////////////////////////////////////////////////////////////////////////////
void
copy_rna_torsions( Size const new_pos, Size const src_pos, pose::Pose & new_pose, pose::Pose & src_pose)
{
	for ( Size j=1; j <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
		id::TorsionID    my_ID( new_pos, id::BB, j );
		id::TorsionID input_ID( src_pos, id::BB, j );
		new_pose.set_torsion( my_ID, src_pose.torsion( input_ID ) );
	}

	for ( Size j=1; j <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {
		id::TorsionID    my_ID(  new_pos, id::CHI, j );
		id::TorsionID input_ID(  src_pos, id::CHI, j );
		new_pose.set_torsion( my_ID, src_pose.torsion( input_ID ) );
	}
}




//////////////////////////////////////////////////////////////////////////////
void
rna_idealize_test() {

	using namespace options;

	pose::Pose pose;
	utility::vector1 <std::string> pdb_files ( option[ in::file::s ]() );

	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	bool const close_chainbreaks = option[ idl_close_chainbreaks ];

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {

		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		if ( !close_chainbreaks ) core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );

		pose::Pose const start_pose( pose );

		pose = start_pose;

		protocols::idealize::IdealizeMover idealizer;

		// set some options
		if ( option[ coordinate_constraint_weight ].user() ) {
			idealizer.coordinate_constraint_weight( option[ coordinate_constraint_weight ] ) ;
		}
		if ( option[ atom_pair_constraint_weight ].user() ) {
			idealizer.atom_pair_constraint_weight( option[ atom_pair_constraint_weight ] );
		}
		idealizer.fast( false /* option[ fast ] */ );

		idealizer.apply( pose );

		// confirm that nothing changes:
		// actually something *does* change!
		idealizer.apply( pose );

		//test this.
		//  pose::Pose refold_pose;
		//  std::string refold_sequence = pose.sequence();
		//  refold_sequence.erase( refold_sequence.size()-1 );
		//  std::cout << "ABOUT TO MAKE POSE FROM SEQUENCE " << refold_sequence << std::endl;
		//  core::chemical::make_pose_from_sequence( refold_pose, refold_sequence, *rsd_set );
		//  std::cout << "HEY! " << pose.size() << " " << refold_pose.size() << std::endl;
		//  refold_pose.dump_pdb( "extended.pdb" );
		//  for (Size i = 1; i <= refold_pose.size(); i++ ) copy_rna_torsions( i, i, refold_pose, pose );
		//  refold_pose.dump_pdb( "refold.pdb" );

		pose.dump_pdb( "idealize_"+pdb_file );

	}
}


///////////////////////////////////////////////////////////////////////
void
print_torsions_check( pose::Pose & pose )
{
	for ( Size i = 1 ; i < pose.size(); i ++ ) {
		std::cout << "POSITION " << I(3,i);
		for ( Size j=1; j <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
			id::TorsionID    my_ID( i, id::BB, j );
			std::cout << " BB" << j << " " << F(8,3,pose.torsion(my_ID) ) << ";";
		}

		for ( Size j=1; j <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {
			id::TorsionID    my_ID(  i, id::CHI, j );
			std::cout << " CHI" << j << " " << F(8,3,pose.torsion(my_ID) ) << ";";
		}
		std::cout << std::endl;
	}

}

///////////////////////////////////////////////////////////////////////
void
rna_torsion_check_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	pose::Pose pose;

	core::sequence::Sequence fasta_sequence = *(core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]);
	make_pose_from_sequence(
		pose,
		fasta_sequence.sequence(),
		*rsd_set );

	std::cout << " --------- ORIGINAL ---------" << std::endl;
	print_torsions_check( pose );
	pose.dump_pdb( "S1.pdb" );

	Size const nres( pose.size() );

	pose::Pose pose2 = pose;
	kinematics::FoldTree f( nres );
	f.new_jump(nres/4,(3*nres/4),nres/2);
	pose2.fold_tree( f );

	//Force a refold.
	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	(*lores_scorefxn)( pose2 );

	std::cout << " --------- NEW FOLDTREE------" << std::endl;
	print_torsions_check( pose2 );
	pose.dump_pdb( "S2.pdb" );

	for ( Size i =1; i<=nres; i++ ) copy_rna_torsions(i,i,pose2,pose);
	(*lores_scorefxn)( pose2 );

	std::cout << " --------- COPY TORSIONS-----" << std::endl;
	print_torsions_check( pose2 );
	pose.dump_pdb( "S3.pdb" );


}


///////////////////////////////////////////////////////////////////////
void
rna_close_chainbreaks_test(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::optimization;

	ResidueTypeSetCOP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Just a check.
	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	lores_scorefxn->set_weight( linear_chainbreak, 1.0 );

	std::string infile  = option[ in::file::s ][1];
	Pose pose;
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);
	/////////////////////////////////////////
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
	/////////////////////////////////////////

	pose.dump_pdb( "start.pdb" );

	lores_scorefxn->show( std::cout, pose );

	// Make copy of pose to save, and screw with a new fold tree
	Pose start_pose( pose );
	FoldTree f( pose.size() );
	f.new_jump( 2, 7, 6 );
	pose.fold_tree( f );

	// Add chainbreak overlap atoms?
	add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, 6 );
	add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, 7 );

	pose.dump_pdb( "overlap.pdb" );

	// score?
	lores_scorefxn->show( std::cout, pose );

	// Copy bb torsions over.
	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j = 1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
			id::TorsionID torsion_id( i, id::BB, j );
			pose.set_torsion( torsion_id, start_pose.torsion( torsion_id ) ) ;
		}
	}


	pose.dump_pdb( "set_torsion.pdb" );
	lores_scorefxn->show( std::cout, pose );

	//minimize?
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;

	mm.set_bb(  true );
	mm.set_chi( true );
	mm.set_jump( true );

	minimizer.run( pose, mm, *(lores_scorefxn), options );
	lores_scorefxn->show( std::cout, pose );

	pose.dump_pdb( "minimize.pdb" );

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
Size
is_regular_helix( pose::Pose const & pose,
	Size const & seqpos ){

	using namespace core::chemical::rna;
	using namespace core::scoring::rna;
	using namespace core::chemical;
	using namespace core::conformation;

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	EnergyBasePairList const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	bool forms_canonical_base_pair( false );
	// bool forms_non_canonical_base_pair( false );

	Size k( 0 ), m( 0 );
	Size partner( 0 ), canonical_partner( 0 );
	for ( EnergyBasePairList::const_iterator it = scored_base_pair_list.begin();
			it != scored_base_pair_list.end(); ++it ) {

		BasePair const base_pair = it->second;

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		if ( i == seqpos ) {
			k = base_pair.edge1();
			m = base_pair.edge2();
			partner = j;
		} else if ( j == seqpos ) {
			k = base_pair.edge2();
			m = base_pair.edge1();
			partner = i;
		} else {
			continue;
		}

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
				&& base_pair.orientation() == 1 )  &&
				possibly_canonical( rsd_i.aa(), rsd_j.aa() ) ) { //&&
			//     pose.torsion( id::TorsionID( i, id::CHI, 1 ) ) > 0  && //Need to check syn/anti
			//     pose.torsion( id::TorsionID( j, id::CHI, 1 ) ) > 0     //Need to check syn/anti

			std::string atom1, atom2;
			get_watson_crick_base_pair_atoms( rsd_i.type(), rsd_j.type(), atom1, atom2 );
			if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() < 3.5 ) {
				forms_canonical_base_pair = true;
				canonical_partner = partner;
				break;
			}  //else {
			//   forms_non_canonical_base_pair = true;
			//   break;
			//  }
		}
	}

	if ( forms_canonical_base_pair ) { //&& !forms_non_canonical_base_pair) {
		return canonical_partner;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
namespace benchmark_ns {
ObjexxFCL::FArray1D < Size > vanilla_helix;
ObjexxFCL::FArray1D < Size > domain_map;
Size num_domain( 0 );
utility::vector1< utility::vector1< Size > > domain_residue_list;
Size nres;
FArray2D <bool> hbonded;
utility::vector1< std::string > colors;
}

///////////////////////////////////////////////////////////////////////////////
void
figure_out_domain_neighbors( Size const & i, core::pose::Pose const & pose  ){

	using namespace benchmark_ns;
	if ( vanilla_helix( i ) ) return;
	if ( domain_map( i ) > 0 ) {
		//Check!!
		if ( domain_map( i ) != num_domain ) utility_exit_with_message( "Problem at seqpos "+I(3,i)+" ==> "+I(3,num_domain)+" "+I(3,domain_map( i )) );
		return;
	}

	// std::cout << "Assigning domain: " <<  i << " " << num_domain << std::endl;
	domain_map( i ) = num_domain;
	// std::cout << "Figure out i-1 " << std::endl;
	if ( i > 1 ) figure_out_domain_neighbors( i - 1, pose );
	// std::cout << "Figure out i+1 " << std::endl;
	if ( i < nres ) figure_out_domain_neighbors( i + 1, pose );

	// RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	// RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	// EnergyBasePairList const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	// Reach across to anything base paired.
	// for ( EnergyBasePairList::const_iterator it = local_scored_base_pair_list.begin();
	//    it != local_scored_base_pair_list.end(); ++it ){
	//  BasePair const base_pair = it->second;
	//  Size const pos1 = base_pair.res1;
	//  Size const pos2 = base_pair.res2;
	//  if ( pos1 == i ) figure_out_domain_neighbors( pos2 );
	//  if ( pos2 == i ) figure_out_domain_neighbors( pos1 );
	// }

	// Reach across to anything hydrogen bonded.
	// core::scoring::hbonds::HBondSet const & hbond_set
	//  ( static_cast< core::scoring::hbonds::HBondSet const & >
	//   ( pose.energies().data().get( util::HBOND_SET )));
	// for ( int n = 1; n <= hbond_set.nhbonds(); ++n ) {
	//  core::scoring::hbonds::HBond const & hbond( hbond_set.hbond(n) );
	//  if (i == Size( hbond.don_res() ) ) figure_out_domain_neighbors( hbond.acc_res(), pose );
	//  if (i == Size( hbond.acc_res() ) ) figure_out_domain_neighbors( hbond.don_res(), pose );
	// }

	// the neighbor/energy links
	using namespace core::scoring;
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( utility::graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node(i)->const_edge_list_begin(),
			irue = energy_graph.get_node(i)->const_edge_list_end();
			iru != irue; ++iru ) {

		EnergyEdge const * edge( static_cast< EnergyEdge const *> ( *iru ) );

		Size const j( edge->get_other_ind(i) );

		if ( (*edge)[ hbond_sc ] < 0.0 ) {
			figure_out_domain_neighbors( j, pose );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
figure_out_domain_map( pose::Pose const & pose ) {
	using namespace benchmark_ns;
	domain_map.dimension( nres );
	domain_map = 0;
	num_domain = 0;
	for ( Size i = 1; i <= nres; i++ ) {
		if ( vanilla_helix( i ) )  continue;
		if ( domain_map(i) > 0 ) continue;
		num_domain++;
		figure_out_domain_neighbors( i, pose );
	}

}


///////////////////////////////////////////////////////////////////////////////
void
absorb_surrounded_canonicals(){
	//Absorb a few a canonicals in, if they're surrounded.
	using namespace benchmark_ns;
	for ( Size i = 2; i < nres; i++ ) {
		if ( domain_map(i) == 0 &&  //canonical helix
				domain_map(i-1) > 0 && // surrounded by noncanonicals?
				domain_map(i-1) == domain_map(i+1)  ) {
			domain_map( i ) = domain_map(i-1);
			domain_map( vanilla_helix( i ) ) = domain_map( i ); //This could lead to disaster.
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void setup_domain_list(){
	using namespace benchmark_ns;
	// Check it out.
	domain_residue_list.clear();
	for ( Size n = 1; n <= num_domain; n++ ) {
		utility::vector1< Size > temp_list;
		for ( Size i = 1; i <= nres; i++ )  {
			if ( domain_map(i)==n ) {
				temp_list.push_back( i );
			}
		}
		domain_residue_list.push_back( temp_list );
	}
}

///////////////////////////////////////////////////////////////////////////
void
fill_hbond_neighbors( pose::Pose const & pose  ){

	using namespace benchmark_ns;

	hbonded.dimension( nres, nres );
	hbonded = false;

	//Look through list of hydrogen bonded neighbors -- add to stem list.
	for ( Size i = 1; i <= nres; i++ ) {
		using namespace core::scoring;
		EnergyGraph const & energy_graph( pose.energies().energy_graph() );

		for ( utility::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {

			EnergyEdge const * edge( static_cast< EnergyEdge const *> ( *iru ) );

			Size const j( edge->get_other_ind(i) );

			if ( (*edge)[hbond_sc] < 0.0 ) {
				hbonded(i,j) = true;
			}
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
output_benchmark_stuff( pose::Pose const & pose,
	std::string const & tag,
	utility::vector1<Size> const & domain_res_list,
	std::list<Size> const & stem_res_list,
	Size const & n,
	utility::io::ozstream & pymol_out )
{

	using namespace benchmark_ns;

	// which residues to output?
	std::list <Size > all_res_list = stem_res_list;
	for ( Size m = 1; m <= domain_res_list.size(); m++ ) all_res_list.push_back( domain_res_list[m] );

	all_res_list.sort();
	all_res_list.unique();

	// Slice out these residues.
	pose::Pose mini_pose;
	Size count( 0 );
	std::map< Size, Size > res_map;
	for ( std::list<Size>::iterator it = all_res_list.begin(),
			it_end = all_res_list.end(); it != it_end; ++ it ) {
		res_map[ *it ] = ++count;
		mini_pose.append_residue_by_bond( pose.residue( *it ) );
	}

	mini_pose.dump_pdb( tag );

	///////////////////////////////////////////////////
	Size pos( tag.find( "RNA.pdb" ) );
	std::string param_file( tag );
	param_file.replace( pos, 7, ".prm" );
	utility::io::ozstream params_out( param_file );

	//Demarcate chainbreaks.
	core::pose::rna::figure_out_reasonable_rna_fold_tree( mini_pose );
	if ( mini_pose.num_jump() > 0 ) {
		params_out << "CUTPOINT_OPEN " ;
		for ( Size j = 1; j < mini_pose.size(); j++ ) {
			if ( mini_pose.fold_tree().is_cutpoint( j ) ) {
				params_out << " " << I(3,j);
			}
		}
		params_out << std::endl;
	}

	///////////////////////////////////////////////////////
	// Some pymol crap.
	pymol_out << "select  chunk"<< lead_zero_string_of(n,3) << ", resi ";
	for ( std::list<Size>::iterator it = all_res_list.begin(),
			it_end = all_res_list.end(); it != it_end; ++ it ) {
		pymol_out << *it << "+";
	}
	pymol_out << "0" << std::endl;
	pymol_out << "color " << colors[ ( ( n-1 ) % colors.size() ) + 1 ] <<",  chunk"<< lead_zero_string_of(n,3) << std::endl;

	pymol_out << "label native and resi " << domain_res_list[1] << " and name P, " << "\"chunk"<<lead_zero_string_of(n,3)<<"\"" << std::endl;

	pymol_out << "select  stems"<< lead_zero_string_of(n,3) << ", resi ";
	for ( std::list<Size>::const_iterator it = stem_res_list.begin(),
			it_end = stem_res_list.end(); it != it_end; ++ it ) {
		pymol_out << *it << "+";
	}
	pymol_out << "0" << std::endl;
	pymol_out << "color purple,  stems"<< lead_zero_string_of(n,3) << std::endl;

	/////////////////////////////////////////////////////////////
	// Figure out stems to force.
	FArray1D <bool> stem_residue_written_out( nres, false );
	// Note that this doesn't (yet) group stem residues together --
	//  may be necessary for pseudoknots!
	for ( std::list<Size>::const_iterator it = stem_res_list.begin(),
			it_end = stem_res_list.end(); it != it_end; ++ it ) {
		Size const i( *it );
		if ( !stem_residue_written_out( i ) ) {
			params_out << "STEM PAIR " << res_map[i] << " " << res_map[vanilla_helix(i)] << " W W A " << std::endl;
			stem_residue_written_out( i ) = true;
			stem_residue_written_out( vanilla_helix( i ) ) = true;
		}
	}


	/////////////////////////////////////////////////////////////////
	// Need to do a quick check -- what if two segments separated by
	// a chainbreak have no potential canonical pairing? Could
	// be a problem -- create a "stem" with appropriate tertiary pairings?


	params_out.close();


	///////////////////////////////////////////////////
	std::string fasta_file( tag );
	pos = tag.find( "RNA.pdb" );
	fasta_file.replace( pos, 7, ".fasta" );
	utility::io::ozstream fasta_out( fasta_file );
	fasta_out << ">" << tag << std::endl;
	fasta_out << mini_pose.sequence() << std::endl;
	fasta_out.close();


}

///////////////////////////////////////////////////////////////////////////////
void
figure_out_stems( pose::Pose & pose, std::string const & infile, utility::io::ozstream & pymol_out )
{
	using namespace benchmark_ns;

	for ( Size n = 1; n <= num_domain; n++ ) {

		std::list < Size > stem_residues;
		std::cout << "DOMAIN " << I(3,n) << "==> ";

		if ( domain_residue_list[n].size() < 3 ) continue;

		for ( Size m = 1; m <= domain_residue_list[n].size(); m++ ) {

			Size const i = domain_residue_list[n][m];
			std::cout << " " << i;

			for ( Size j = 1; j <= nres; j++ ) {
				if ( hbonded(i,j) && vanilla_helix(j) ) {
					if ( !domain_map( j) || !domain_map( vanilla_helix(j) ) ) {
						stem_residues.push_back( j );
						stem_residues.push_back( vanilla_helix(j) );
					}
				}
			}

			if ( i > 1  && vanilla_helix(i-1) ) {
				Size const j = i-1;
				if ( !domain_map( j) || !domain_map( vanilla_helix(j) ) ) {
					stem_residues.push_back( j );
					stem_residues.push_back( vanilla_helix(j) );
				}
			}

			if ( i < nres  && vanilla_helix(i+1) ) {
				Size const j = i+1;
				if ( !domain_map( j) || !domain_map( vanilla_helix(j) ) ) {
					stem_residues.push_back( j );
					stem_residues.push_back( vanilla_helix(j) );
				}
			}

		}
		std::cout << std::endl;


		stem_residues.sort();
		stem_residues.unique();

		std::cout << "  ASSOCIATED STEM RESIDUES ==> ";
		for ( std::list<Size>::iterator it = stem_residues.begin(),
				it_end = stem_residues.end(); it != it_end; ++ it ) {
			std::cout << " " << ( *it );
		}
		std::cout << std::endl;

		std::string const tag = "chunk" + lead_zero_string_of( n, 3 ) + "_" + infile;
		output_benchmark_stuff( pose, tag, domain_residue_list[n], stem_residues, n, pymol_out );

	}

}

/////////////////////////////////////////////////////////////////////////
void
initialize_pymol_colors(){
	using namespace benchmark_ns;
	colors.clear();
	colors.push_back( "blue" );
	colors.push_back( "green" );
	colors.push_back( "red" );
	colors.push_back( "cyan" );
	colors.push_back( "yellow" );
	colors.push_back( "orange" );
	colors.push_back( "lightblue" );
	colors.push_back( "forest" );
	colors.push_back( "pink" );
}

///////////////////////////////////////////////////////////////////////////////
void
create_rna_benchmark_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;

	// hmmm.
	using namespace benchmark_ns;

	initialize_pymol_colors();

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	utility::vector1< std::string > const infiles( option[ in::file::s ]() );

	for ( Size i = 1; i <= infiles.size(); i++ ) {

		std::string infile  = infiles[i];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);
		/////////////////////////////////////////
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		nres = pose.size();

		core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );

		// Need some stuff to figure out which residues are base paired. First score.
		ScoreFunctionOP scorefxn( new ScoreFunction );
		scorefxn->set_weight( rna_base_pair, 1.0 );
		scorefxn->set_weight( hbond_sc, 1.0 );
		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		//Divide up the structure based on "modules" separated by canonical RNA base pair helices.
		// Set up list of residues that make canonical base pairs and nothing else.
		vanilla_helix.dimension( nres );
		vanilla_helix = 0;
		std::cout << "Vanilla helix setup... " << std::endl;
		for ( Size i = 1; i <= nres; i++ ) vanilla_helix( i ) = is_regular_helix( pose, i );

		std::cout << "Domain map setup... " << std::endl;
		// March through each residue.
		fill_hbond_neighbors( pose );
		figure_out_domain_map( pose );

		absorb_surrounded_canonicals();

		setup_domain_list();

		utility::io::ozstream pymol_out( infile+".pml" );
		pymol_out << "reinitialize "  << std::endl;
		pymol_out << "load " << infile << ",native  "  << std::endl;
		pymol_out << "set cartoon_ring_mode, 1" << std::endl;
		pymol_out << "cartoon oval" << std::endl;
		pymol_out << "set cartoon_oval_width, 0.4" << std::endl;
		pymol_out << "set cartoon_oval_length, 0.8" << std::endl;
		pymol_out << "show cartoon, native " << std::endl;
		pymol_out << "color white, native " << std::endl;
		pymol_out << "bg_color white" << std::endl;


		figure_out_stems( pose, infile, pymol_out );

		std::string vall_file( infile );
		Size const pos = infile.find( "_RNA.pdb" );
		vall_file.replace( pos, 8, ".torsions" );
		utility::vector1< Size > blank_list;
		protocols::farna::create_rna_vall_torsions( pose, vall_file, blank_list );

		pymol_out << "color black, resi 1 " << std::endl;
		pymol_out << "label native and resi 1 and name P, \"" << infile << "\"" << std::endl;
		pymol_out << "set label_size, 18 " << std::endl;
		pymol_out.close();

	}
}


///////////////////////////////////////////////////////////////
void
rna_chain_closure_test()
{
	using namespace chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	std::string infile  = option[ in::file::s ][1];

	pose::Pose pose;
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);
	/////////////////////////////////////////
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
	/////////////////////////////////////////

	//Save a copy.
	pose::Pose start_pose;
	start_pose = pose;

	//Put in chain break
	Size const cutpoint( 8 );
	kinematics::FoldTree f( pose.size() );
	f.new_jump( cutpoint, cutpoint+1, cutpoint );
	f.set_jump_atoms( 1,
		core::chemical::rna::chi1_torsion_atom( pose.residue_type( cutpoint) ),
		core::chemical::rna::chi1_torsion_atom( pose.residue_type( cutpoint+1) )   );
	pose.fold_tree( f );

	add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint );
	add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint+1 );

	for ( Size i = cutpoint; i <= cutpoint+1; i++ ) {
		for ( Size j=1; j <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
			id::TorsionID my_ID( i, id::BB, j );
			pose.set_torsion( my_ID, start_pose.torsion( my_ID ) );
		}
	}
	pose.dump_pdb( "start.pdb" );

	//Perturb torsions
	utility::vector1< id::TorsionID > tor_ids;
	for ( Size j = 5; j <=6; ++ j ) {
		tor_ids.push_back( id::TorsionID( cutpoint, id::BB, j ) );
	}
	for ( Size j = 1; j <=4; ++ j ) {
		tor_ids.push_back( id::TorsionID( cutpoint+1, id::BB, j ) );
	}

	for ( Size n = 1; n <= tor_ids.size(); n++ ) {
		//Size const n(1);
		Real const current_val = pose.torsion( tor_ids[n] );
		pose.set_torsion( tor_ids[n], current_val + 20.0);
	}

	//Save "perturbed" configuration
	pose.dump_pdb( "perturb.pdb" );

	protocols::farna::movers::RNA_LoopCloser rna_ccd_closer;
	Real mean_dist_err = rna_ccd_closer.apply( pose, cutpoint );
	std::cout << "CCD closure at residue " << cutpoint << " ==> final mean distance error: " << mean_dist_err << std::endl;

	pose.dump_pdb( "final.pdb" );

}

////////////////////////////////////////////////////////////////////
Real
get_jump_distance( kinematics::Jump const & jump1, kinematics::Jump const & jump2 ){
	Real dist, theta;
	kinematics::jump_distance( jump1, jump2, dist, theta );
	return ( dist + 4.0 * theta );
}

//////////////////////////////////////////////////////////////////////////////////////
void
setup_crazy_fold_tree( pose::Pose & pose, core::chemical::ResidueTypeSetCOP & rsd_set )
{
	using namespace chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::conformation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::kinematics;
	using namespace core::id;

	pose::Pose original_pose = pose;
	core::pose::rna::figure_out_reasonable_rna_fold_tree( original_pose );

	std::string rna_sequence = pose.sequence();
	Size const nres_real( pose.size() );
	ResidueOP rsd( ResidueFactory::create_residue( rsd_set->name_map( "VRT" ) ) );
	pose.append_residue_by_jump( *rsd, 1 );

	//Crazy fold tree.
	Size const num_jumps( nres_real );
	ObjexxFCL::FArray2D <int> jump_points( 2, num_jumps );
	ObjexxFCL::FArray1D <int> cuts( num_jumps );
	for ( Size n = 1; n <= num_jumps; n++ ) {
		jump_points( 1, n ) = n;
		jump_points( 2, n ) = nres_real+1;
		cuts( n ) = n;
	}

	FoldTree f( nres_real+1 );
	f.tree_from_jumps_and_cuts( nres_real + 1, num_jumps,
		jump_points, cuts, 1, false /*verbose*/ );
	f.reorder( nres_real + 1 );
	for ( Size n = 1; n <= num_jumps; n++ ) {
		f.set_jump_atoms( n,
			"ORIG",
			core::chemical::rna::chi1_torsion_atom( pose.residue_type( n ) ) );

	}

	pose.fold_tree( f );

	for ( Size n = 1; n < nres_real; n++ ) {
		if ( !original_pose.fold_tree().is_cutpoint( n ) ) {
			add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, n );
			add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, n+1 );
		}
	}

}

///////////////////////////////////////////////////////////////
void
rna_backbone_rebuild_test()
{
	using namespace chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::conformation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::kinematics;
	using namespace core::id;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	std::string infile  = option[ in::file::s ][1];

	//Input pose, tack on virtual atom.
	pose::Pose pose;
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);
	/////////////////////////////////////////
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );
	/////////////////////////////////////////

	std::string rna_sequence = pose.sequence();
	setup_crazy_fold_tree( pose, rsd_set );

	pose.dump_pdb( "start.pdb" );

	//Save a copy.
	pose::Pose start_pose;
	start_pose = pose;

	//Now create idealized pose.
	pose.clear();

	ResidueOP rsd( ResidueFactory::create_residue( rsd_set->name_map( "VRT" ) ) );
	make_pose_from_sequence( pose, rna_sequence, *rsd_set );
	pose.append_residue_by_jump( *rsd, 1 );
	pose.fold_tree( start_pose.fold_tree() );
	Size const nres_real( rna_sequence.size()  );
	Size const num_jumps( nres_real );
	for ( Size n = 1; n <= num_jumps; n++ ) {
		pose.set_jump( n, start_pose.jump( n ) );
	}

	for ( Size n = 1; n < nres_real; n++ ) {
		add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, n );
		add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, n+1 );
	}

	pose.dump_pdb( "backbone_extend.pdb" );

	///////////////////////////////////////////////////////////////////////
	// Now read in vall_torsions
	utility::io::izstream vall_in( option[ vall_torsions ] );
	std::string line, tag;

	utility::vector1< utility::vector1< Real > >  torsion_sets;

	while (  getline( vall_in, line) ) {
		std::istringstream line_stream( line );
		char dummychar( ' ' );
		Real val( 0.0 );
		line_stream >> dummychar;
		utility::vector1 < Real > torsion_set;
		for ( Size i = 1; i <= core::chemical::rna::NUM_RNA_TORSIONS; i++ ) {
			line_stream  >> val;
			torsion_set.push_back( val );
		}
		torsion_sets.push_back( torsion_set );
	}

	/////////////////////////////////////////////////////////////////////////////
	// Figure out a mapping between torsion sets and jumps
	pose::Pose mini_pose;
	make_pose_from_sequence( mini_pose, "uu", *rsd_set );
	FoldTree mini_f( 2 );
	mini_f.new_jump( 1, 2, 1);
	mini_f.set_jump_atoms( 1, " C2 ", " C2 " );
	// mini_pose.fold_tree( mini_f );
	// add_variant_type_to_pose_residue( mini_pose, chemical::CUTPOINT_LOWER, 1 );
	// add_variant_type_to_pose_residue( mini_pose, chemical::CUTPOINT_UPPER, 2 );

	utility::vector1< Jump > base_base_jumps, forwardconnect_base_jumps;
	std::cout << "NUM TORSION SETS " << torsion_sets.size() << std::endl;

	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( hbond_sc, 1.0 );

	for ( Size n = 1; n < torsion_sets.size(); n++ ) {
		mini_pose.fold_tree( FoldTree( 2 ) ); //simple.
		for ( Size j = 1; j <= NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
			mini_pose.set_torsion( id::TorsionID( 1, id::BB, j), torsion_sets[n][j] );
			mini_pose.set_torsion( id::TorsionID( 2, id::BB, j), torsion_sets[n+1][j] );
		}
		for ( Size j = 1; j <= NUM_RNA_CHI_TORSIONS; ++ j ) {
			mini_pose.set_torsion( id::TorsionID( 1, id::CHI, j), torsion_sets[n][j+ NUM_RNA_MAINCHAIN_TORSIONS] );
			mini_pose.set_torsion( id::TorsionID( 2, id::CHI, j), torsion_sets[n+1][j+NUM_RNA_MAINCHAIN_TORSIONS] );
		}

		//Need to force a refold?
		(*scorefxn)( mini_pose );
		//  mini_pose.dump_pdb( "S_"+string_of( n )+".pdb" );

		mini_pose.fold_tree( mini_f );

		id::AtomID base_atom1_id(       mini_pose.residue(1).atom_index( " C1'" ), 1 );
		tree::AtomCOP base_atom1 ( mini_pose.atom_tree().atom( base_atom1_id ).get_self_ptr() );

		id::AtomID forward_connect1_id( mini_pose.residue(1).atom_index( " O3'" ), 1 );
		tree::AtomCOP forward_connect1 ( mini_pose.atom_tree().atom( forward_connect1_id ).get_self_ptr() );

		id::AtomID base_atom2_id(       mini_pose.residue(2).atom_index( " C1'" ), 2 );
		tree::AtomCOP base_atom2 ( mini_pose.atom_tree().atom( base_atom2_id ).get_self_ptr() );

		base_base_jumps.push_back( Jump( base_atom1->get_stub(), base_atom2->get_stub() ) );
		forwardconnect_base_jumps.push_back( Jump( forward_connect1->get_stub(), base_atom2->get_stub() ) );

	}

	/////////////////////////////////////////////////////////////////////////////
	// Find best jump to superimpose onto current residue -- do residue 1 first.
	{
		Size const i( 1 );

		id::AtomID base_atom1_id(       pose.residue(i).atom_index( " C1'" ), i );
		tree::AtomCOP base_atom1 ( pose.atom_tree().atom( base_atom1_id ).get_self_ptr() );
		id::AtomID base_atom2_id(       pose.residue(i+1).atom_index( " C1'" ), i+1 );
		tree::AtomCOP base_atom2 ( pose.atom_tree().atom( base_atom2_id ).get_self_ptr() );
		Jump current_jump( base_atom1->get_stub(), base_atom2->get_stub() );

		Real best_jump_distance( 999.9 );
		Size best_n( 1 );
		for ( Size n = 1; n <= base_base_jumps.size(); n++ ) {
			Real const jump_distance = get_jump_distance( current_jump, base_base_jumps[ n ] );
			//  std::cout <<  "JUMP DIST ==> " << n << " " << jump_distance << " " << base_base_jumps[n] << std::endl;
			if ( jump_distance < best_jump_distance ) {
				best_jump_distance = jump_distance;
				best_n = n;
				//    std::cout << "Found better jump: " << best_n  << " " << best_jump_distance << std::endl;
			}
		}

		for ( Size j = 4; j <= NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
			pose.set_torsion( id::TorsionID( i, id::BB, j), torsion_sets[best_n][j] );
		}
		for ( Size j = 1; j <= 4; ++j ) {
			pose.set_torsion( id::TorsionID( i+1, id::BB, j), torsion_sets[best_n+1][j] );
		}
		for ( Size j = 1; j <= NUM_RNA_CHI_TORSIONS; ++ j ) {
			pose.set_torsion( id::TorsionID( i, id::CHI, j), torsion_sets[best_n][j+NUM_RNA_MAINCHAIN_TORSIONS] );
			pose.set_torsion( id::TorsionID( i+1, id::CHI, j), torsion_sets[best_n+1][j+NUM_RNA_MAINCHAIN_TORSIONS] );
		}

		pose.dump_pdb( "close1.pdb" );
	}


	/////////////////////////////////////////////////////////////////////////////
	// Now close the rest of the residues...
	for ( Size i = 2; i < nres_real; i++ ) {

		id::AtomID base_atom1_id(       pose.residue(i).atom_index( " O3'" ), i );
		tree::AtomCOP base_atom1 ( pose.atom_tree().atom( base_atom1_id ).get_self_ptr() );
		id::AtomID base_atom2_id(       pose.residue(i+1).atom_index( " C1'" ), i+1 );
		tree::AtomCOP base_atom2 ( pose.atom_tree().atom( base_atom2_id ).get_self_ptr() );
		Jump current_jump( base_atom1->get_stub(), base_atom2->get_stub() );

		Real best_jump_distance( 999.9 );
		Size best_n( 1 );
		for ( Size n = 1; n <= base_base_jumps.size(); n++ ) {
			Real const jump_distance = get_jump_distance( current_jump, forwardconnect_base_jumps[ n ] );
			//  std::cout <<  "JUMP DIST ==> " << n << " " << jump_distance << " " << base_base_jumps[n] << std::endl;
			if ( jump_distance < best_jump_distance ) {
				best_jump_distance = jump_distance;
				best_n = n;
			}
		}
		std::cout << "Found best jump: " << best_n  << " " << best_jump_distance << std::endl;

		for ( Size j = 5; j <= NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
			pose.set_torsion( id::TorsionID( i, id::BB, j), torsion_sets[best_n][j] );
		}
		for ( Size j = 1; j <= 4; ++j ) {
			pose.set_torsion( id::TorsionID( i+1, id::BB, j), torsion_sets[best_n+1][j] );
		}
		for ( Size j = 1; j <= NUM_RNA_CHI_TORSIONS; ++ j ) {
			pose.set_torsion( id::TorsionID( i+1, id::CHI, j), torsion_sets[best_n+1][j+NUM_RNA_MAINCHAIN_TORSIONS] );
		}
	}
	pose.dump_pdb( "close.pdb" );

	// CCD close?
	protocols::farna::movers::RNA_LoopCloser rna_loop_closer;
	for ( Size i = 1; i < nres_real; i++ ) {

		//Real mean_dist_err = rna_ccd_close( pose, i );

		Real mean_dist_err = rna_loop_closer.apply( pose, i );
		std::cout << "CCD CLOSURE ==> " << i << " : " << mean_dist_err << std::endl;

	}
	pose.dump_pdb( "ccd_close.pdb" );

}

///////////////////////////////////////////////////////////////////////////////
void
output_sasa( std::ostream & out,
	pose::Pose & pose,
	Size const rsdno,
	id::AtomID_Map< Real > & atom_sasa,
	utility::vector1< std::string> const & atom_names,
	bool const & output )
{
	if ( !output ) {
		out << F(8,3, 0.0 );
		return;
	}

	bool atom_found( true );
	Real atom_sasa_tot( 0.0 );

	for ( Size m = 1; m <= atom_names.size(); m++ ) {
		Size const atom_i( pose.residue( rsdno ).atom_index( atom_names[m] ) );
		if ( atom_i > 0 ) {
			id::AtomID atom_id( atom_i, rsdno );
			atom_sasa_tot += atom_sasa[ atom_id ];
			//out << F(8, 3, atom_sasa[ atom_id ] );
		} else {
			atom_found = false;
			//out << F(8, 3, 0.0 );
		}
	}

	if ( atom_found ) {
		out << F(8, 3, atom_sasa_tot );
	} else {
		out << F(8,3, 0.0 );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
sasa_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );


	utility::vector1< std::string > atom_names_sugar, atom_names_sugar_hydrogen, atom_names_phosphate,
		atom_names_dms, atom_names_kethoxal;
	atom_names_sugar.push_back( " C1'");
	atom_names_sugar.push_back( " C2'");
	atom_names_sugar.push_back( " C3'");
	atom_names_sugar.push_back( " C4'");
	atom_names_sugar.push_back( " C5'");

	atom_names_sugar_hydrogen.push_back( " H1'");
	atom_names_sugar_hydrogen.push_back( " H2'");
	atom_names_sugar_hydrogen.push_back( " H3'");
	atom_names_sugar_hydrogen.push_back( " H4'");
	atom_names_sugar_hydrogen.push_back( " H5'");
	atom_names_sugar_hydrogen.push_back( "H5''");

	atom_names_phosphate.push_back( " OP2" );
	atom_names_phosphate.push_back( " OP1" );

	atom_names_dms.push_back( " N1 " );
	atom_names_dms.push_back( " N3 " );

	atom_names_kethoxal.push_back( " N1 " );

	Real const probe_radius( 1.4 );

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		//  protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );

		id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, true );

		for ( Size n = 1; n <= pose.size(); n++ ) {

			std::cout << I(3,n) << " " << F(8,3,rsd_sasa[n]);

			std::string resname = pose.residue(n).name3();

			//std::cout << resname;

			output_sasa( std::cout, pose, n, atom_sasa, atom_names_sugar, true );
			output_sasa( std::cout, pose, n, atom_sasa, atom_names_sugar_hydrogen, true );
			output_sasa( std::cout, pose, n, atom_sasa, atom_names_phosphate, true );
			output_sasa( std::cout, pose, n, atom_sasa, atom_names_dms, (resname == "  A" || resname == "  C")  );
			output_sasa( std::cout, pose, n, atom_sasa, atom_names_kethoxal, (resname == "  G") );

			// OP2/OP1
			// A or C watson/crick N (DMS)
			// G watson/crick N


			std::cout << std::endl;

		}

	}


}

///////////////////////////////////////////////////////////////////////////////
void
env_sugar_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( option[ rsd_type_set ]()  );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );


	utility::vector1< std::string > atom_names_sugar;
	if ( option[ rsd_type_set ]() == "coarse_rna" ) {
		atom_names_sugar.push_back( " S  ");
	} else {
		atom_names_sugar.push_back( " C1'");
		atom_names_sugar.push_back( " C2'");
		atom_names_sugar.push_back( " C3'");
		atom_names_sugar.push_back( " C4'");
		atom_names_sugar.push_back( " C5'");
	}

	Real const probe_radius_bin_size( 2.0 );
	Real const max_probe_radius( 16.0 );
	Size const numbins( static_cast<int>(max_probe_radius/probe_radius_bin_size) );
	Size const num_sugar_atoms( atom_names_sugar.size() );

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		//  protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );

		for ( Size i = 1; i <= pose.size(); i++ ) {

			utility::vector1< numeric::xyzVector<Real> > sugar_atoms_xyz;
			numeric::xyzVector< Real > sugar_atom_centroid( 0.0 );
			for ( Size j = 1; j <= num_sugar_atoms; j++ ) {
				numeric::xyzVector< Real > const & sugar_atom( pose.residue(i).xyz( atom_names_sugar[j] ) );
				sugar_atoms_xyz.push_back( sugar_atom );
				sugar_atom_centroid += sugar_atom;
			}
			sugar_atom_centroid /= num_sugar_atoms;
			sugar_atoms_xyz.push_back( sugar_atom_centroid );

			FArray2D< Size > nbrs(  numbins, num_sugar_atoms+2 );
			nbrs = 0;

			for ( Size k = 1; k <= pose.size(); k++ ) {
				conformation::Residue const & rsd( pose.residue(k) );

				for ( Size m = 1; m <= rsd.nheavyatoms(); m++ ) {

					numeric::xyzVector< Real > const & other_atom( rsd.xyz( m ) );
					FArray1D< bool > found_neighbor( numbins, false );

					// Is this atom near any of the sugar atoms in the residue of interest?
					for ( Size j = 1; j <= sugar_atoms_xyz.size(); j++ ) {

						Real const dist = ( other_atom - sugar_atoms_xyz[j] ).length();

						if ( dist < max_probe_radius ) {
							Size const bin = static_cast<Size>( dist / probe_radius_bin_size ) + 1;
							//for ( Size b = 1; b <= bin; b++ ) {
							nbrs( bin, j )++;
							found_neighbor( bin ) = true;
							//       }
						}
					}

					// the "any sugar atom" bin.
					for ( Size bin = 1; bin <= numbins; bin++ ) {
						if ( found_neighbor( bin ) ) {
							nbrs(  bin, num_sugar_atoms+2 )++;
						}
					}

				}

			}

			std::cout << "ENV_SUGAR ";

			for ( Size j = 1; j <= num_sugar_atoms+2; j++ ) {
				std::cout << "     ";
				for ( Size bin = 1; bin <= numbins; bin++ ) {
					std::cout <<  ' '<< nbrs( bin, j );
				}
			}

			std::cout << std::endl;
		}

	}


}


///////////////////////////////////////////////////////////////////////////////
void
print_hbonds_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( option[ rsd_type_set]() );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	ScoreFunctionOP scorefxn = get_score_function();

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);

		(*scorefxn)(pose);
		hbonds::HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
		hbond_options->use_hb_env_dep( false );
		hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( *hbond_options ) );

		hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

		core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
		for ( Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
			hbonds::HBond const & hbond( hbond_set->hbond( i ) );

			Size const don_res_num = hbond.don_res();
			Size const don_hatm = hbond.don_hatm();

			Size const acc_res_num = hbond.acc_res();
			Size const acc_atm = hbond.acc_atm();


			std::cout << "HBOND: " <<
				pdb_info->chain( don_res_num ) << ":" <<
				pose.residue( don_res_num ).name1() <<pdb_info->number( don_res_num ) << " " <<
				pose.residue( don_res_num ).atom_name( don_hatm ) << " --- " <<
				pdb_info->chain( acc_res_num ) << ":" <<
				pose.residue( acc_res_num).name1() << pdb_info->number( acc_res_num ) << " " <<
				pose.residue( acc_res_num ).atom_name( acc_atm ) <<
				" ==> " << hbond.energy()
				<< std::endl;

		}

	}
}


///////////////////////////////////////////////////////////////////////////////
void
create_dihedral_constraint( core::scoring::constraints::ConstraintSetOP & cst_set,
	core::pose::Pose & pose,
	core::id::TorsionID torsion_id )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;

	id::AtomID id1,id2,id3,id4;
	pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );

	Real torsion_value( pose.torsion( torsion_id ) );

	HarmonicFuncOP harm_func( new HarmonicFunc( numeric::conversions::radians( torsion_value ),
		numeric::conversions::radians( 20.0 ) ) );

	ConstraintOP dihedral1( new DihedralConstraint( id1, id2, id3, id4,
		harm_func,
		rna_torsion ) );
	cst_set->add_constraint( dihedral1 );
}


///////////////////////////////////////////////////////////////////////////////
void
get_backbone_rotamers( utility::vector1< utility::vector1 <Real > > & backbone_rotamers,
	PuckerState const & pucker1,
	PuckerState const & pucker2 ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

	utility::vector1< Real > torsion_samples;
	if ( option[ more_rotamers ] ) {
		torsion_samples.push_back( -0.5 );
		torsion_samples.push_back( 0.5 );
	} else  {
		torsion_samples.push_back( 0 );
	}

	backbone_rotamers.clear();

	Real delta1, nu2_1, nu1_1, epsilon1, zeta1;
	Real alpha2, beta2, gamma2, delta2, nu2_2, nu1_2;

	chemical::rna::GaussianParameter gp( 0.0, 0.0, 0.0);

	for ( int d1 = 1; d1 <= 2; d1++ ) {

		if ( pucker1 > 0 &&  pucker1 != d1 ) continue;

		if ( d1 == 1 ) {
			delta1 = rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center;
			nu2_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center;
			nu1_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center;
		} else {
			delta1 = rna_fitted_torsion_info.gaussian_parameter_set_delta_south()[1].center;
			nu2_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_south()[1].center;
			nu1_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_south()[1].center;
		}

		for ( Size e1 = 1; e1 <= 2; e1++ ) {

			for ( Size e1_std = 1; e1_std <= torsion_samples.size(); e1_std++ ) {

				if ( d1 == 1 ) gp = rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[e1];
				else         gp = rna_fitted_torsion_info.gaussian_parameter_set_epsilon_south()[e1];

				epsilon1 = gp.center + torsion_samples[ e1_std ] * gp.width;

				//zeta depends on alpha of next residue...
				for ( Size a2 = 1; a2 <= 3; a2++ ) {

					for ( Size a2_std = 1; a2_std <= torsion_samples.size(); a2_std++ ) {

						gp = rna_fitted_torsion_info.gaussian_parameter_set_alpha()[a2];
						alpha2 = gp.center + torsion_samples[ a2_std ] * gp.width;

						for ( Size z1 = 1; z1 <= 2; z1++ ) {

							for ( Size z1_std = 1; z1_std <= torsion_samples.size(); z1_std++ ) {

								if ( a2 == 1 )    gp = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[z1];
								else if ( a2==2 ) gp = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_plus()[z1];
								else            gp = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_ap()[z1];

								zeta1 = gp.center + torsion_samples[ z1_std ] * gp.width;


								for ( Size b2 = 1; b2 <= 1; b2++ ) {

									for ( Size b2_std = 1; b2_std <= torsion_samples.size(); b2_std++ ) {

										gp = rna_fitted_torsion_info.gaussian_parameter_set_beta()[b2];
										beta2 = gp.center + torsion_samples[ b2_std ] * gp.width;

										for ( Size g2 = 1; g2 <= 3; g2++ ) {

											for ( Size g2_std = 1; g2_std <= torsion_samples.size(); g2_std++ ) {

												gp = rna_fitted_torsion_info.gaussian_parameter_set_gamma()[g2];
												gamma2 = gp.center + torsion_samples[ g2_std ] * gp.width;

												for ( int d2 = 1; d2 <= 2; d2++ ) {

													if ( pucker2 > 0 &&  pucker2 != d2 ) continue;

													if ( d2 == 1 ) {
														delta2 = rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center;
														nu2_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center;
														nu1_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center;
													} else {
														delta2 = rna_fitted_torsion_info.gaussian_parameter_set_delta_south()[1].center;
														nu2_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_south()[1].center;
														nu1_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_south()[1].center;
													}

													utility::vector1 < Real >  backbone_rotamer; //Would be better to make this a class, RNA_SuiteToSuiteRotamer
													backbone_rotamer.push_back( delta1 );
													backbone_rotamer.push_back( nu2_1 );
													backbone_rotamer.push_back( nu1_1 );
													backbone_rotamer.push_back( epsilon1 );
													backbone_rotamer.push_back( zeta1 );
													backbone_rotamer.push_back( alpha2 );
													backbone_rotamer.push_back( beta2 );
													backbone_rotamer.push_back( gamma2 );
													backbone_rotamer.push_back( delta2 );
													backbone_rotamer.push_back( nu2_2 );
													backbone_rotamer.push_back( nu1_2 );

													backbone_rotamers.push_back( backbone_rotamer );

												} // delta2
											} // gamma2_samples
										} // gamma2
									} // beta2_samples
								} // beta2
							} // zeta1_samples
						} // zeta1
					} // alpha2_samples
				} // alpha2
			} // epsilon1_samples
		} // epsilon1
	} // delta1

}

///////////////////////////////////////////////////////////////////////////////
void
apply_backbone_rotamer( pose::Pose & pose, Size const i, utility::vector1< Real > backbone_rotamer,
	bool const pucker1 = true, bool const pucker2 = true)
{
	using namespace core::id;

	if ( pucker1 ) {
		pose.set_torsion( TorsionID( i  , id::BB,  4 ), backbone_rotamer[1] );
		pose.set_torsion( TorsionID( i  , id::CHI, 2 ), backbone_rotamer[2] );
		pose.set_torsion( TorsionID( i  , id::CHI, 3 ), backbone_rotamer[3] );
	}

	pose.set_torsion( TorsionID( i  , id::BB,  5 ), backbone_rotamer[4] );
	pose.set_torsion( TorsionID( i  , id::BB,  6 ), backbone_rotamer[5] );
	pose.set_torsion( TorsionID( i+1, id::BB,  1 ), backbone_rotamer[6] );
	pose.set_torsion( TorsionID( i+1, id::BB,  2 ), backbone_rotamer[7] );
	pose.set_torsion( TorsionID( i+1, id::BB,  3 ), backbone_rotamer[8] );

	if ( pucker2 ) {
		pose.set_torsion( TorsionID( i+1, id::BB,  4 ), backbone_rotamer[9] );
		pose.set_torsion( TorsionID( i+1, id::CHI, 2 ), backbone_rotamer[10] );
		pose.set_torsion( TorsionID( i+1, id::CHI, 3 ), backbone_rotamer[11] );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
dinucleotide_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	pose::Pose pose;
	std::string sequence = "cc";

	core::pose::make_pose_from_sequence(
		pose,
		sequence,
		*rsd_set );
	dump_pdb( pose, "start.pdb" );

	chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;
	for ( Size i = 1; i <=2; i++ ) { // starting values for torsions.

		pose.set_torsion( TorsionID( i, id::BB, 1 ), rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 2 ), rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 3 ), rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 4 ), rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 5 ), rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
		pose.set_torsion( TorsionID( i, id::BB, 6 ), rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );

		pose.set_torsion( TorsionID( i, id::CHI, 1 ), rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
		pose.set_torsion( TorsionID( i, id::CHI, 2 ), rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center );
		pose.set_torsion( TorsionID( i, id::CHI, 3 ), rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center );
		pose.set_torsion( TorsionID( i, id::CHI, 4 ), 0.0 );


	}

	pose::Pose current_pose( pose );

	protocols::viewer::add_conformation_viewer( current_pose.conformation(), "current", 400, 400 );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	scorefxn->set_weight( dihedral_constraint, 1.0 );

	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	// Now cycle through one of the torsions. Score and output.
	Size count( 1 );

	utility::vector1< utility::vector1 <Real > > backbone_rotamers;
	get_backbone_rotamers( backbone_rotamers, ANY_PUCKER, ANY_PUCKER );

	for ( Size n = 1;  n <= backbone_rotamers.size(); n++ )  {
		apply_backbone_rotamer( pose, 1, backbone_rotamers[n] );

		current_pose = pose;

		ConstraintSetOP cst_set( new ConstraintSet() );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::BB, 4 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::CHI, 1 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::CHI, 2 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::CHI, 3 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::BB, 5 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 1, id::BB, 6 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 1 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 2 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 3 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::BB, 4 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::CHI, 1 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::CHI, 2 ) );
		create_dihedral_constraint( cst_set, current_pose, TorsionID( 2, id::CHI, 3 ) );
		current_pose.constraint_set( cst_set );

		AtomTreeMinimizer minimizer;
		float const dummy_tol( 0.0000025);
		bool const use_nblist( true );
		MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
		options.nblist_auto_update( true );

		kinematics::MoveMap mm;
		mm.set_bb(  true );
		mm.set_chi( true );
		mm.set_jump( true );

		minimizer.run( current_pose, mm, *(scorefxn), options );
		//   (*scorefxn)(current_pose);

		std::string const tag( "S_"+lead_zero_string_of(count++, 3) );
		RNA_SilentStruct s( current_pose, tag );
		silent_file_data.write_silent_struct(s, silent_file, false /*write score only*/);
		dump_pdb( current_pose, tag+".pdb" );
		dump_pdb( pose, "nomin_"+tag+".pdb" );
	}


}

////////////////////////////////////////////////////
void
vary_bond_length( pose::Pose & pose,
	core::id::TorsionID & tor_id,
	core::kinematics::MoveMap & mm,
	core::scoring::constraints::ConstraintSetOP & cst_set )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	AtomID id1,id2,id3,id4, my_ID;
	bool const failure = pose.conformation().get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 );
	if ( failure ) return;

	core::kinematics::tree::AtomCOP atom2 ( pose.atom_tree().atom( id2 ).get_self_ptr() );
	core::kinematics::tree::AtomCOP atom3 ( pose.atom_tree().atom( id3 ).get_self_ptr() );

	DOF_ID dof_id;
	if ( atom2->parent() == atom3 ) {
		my_ID = id2;
	} else if ( atom3->parent() == atom2 ) {
		my_ID = id3;
	} else  {
		utility_exit_with_message( "Problem with atoms " );
	}

	dof_id = DOF_ID( my_ID, D );
	// std::cout << "Attempt to vary bond length for resno " << my_ID.rsd() << " atom: " << pose.residue( my_ID.rsd() ).atom_name( my_ID.atomno() ) << std::endl;

	mm.set( dof_id, true );

	mm.set( DOF_ID( my_ID, PHI) , true );
	mm.set( DOF_ID( my_ID, THETA ), true );

	cst_set->add_dof_constraint( dof_id, core::scoring::func::FuncOP( new HarmonicFunc( (atom2->xyz() - atom3->xyz() ).length() , 0.01 ) ) );

}


/////////////////////////////////////////////////////////////////////////////////
void vary_geometry_RNA( pose::Pose & pose, kinematics::MoveMap & mm )
{

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::scoring::constraints;

	ConstraintSetOP cst_set(  pose.constraint_set()->clone() ) ;

	//Change this to also include D DOF's for sidechains.
	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size n = 1; n <= NUM_RNA_MAINCHAIN_TORSIONS; n++ ) {
			TorsionID tor_id( i, id::BB, n );
			vary_bond_length( pose, tor_id, mm, cst_set );
		}
		for ( Size n = 1; n <= NUM_RNA_CHI_TORSIONS; n++ ) {
			TorsionID tor_id( i, id::CHI, n );
			vary_bond_length( pose, tor_id, mm, cst_set );
		}
	}

	pose.constraint_set( cst_set );

}


///////////////////////////////////////////////////////////////////////////////
void
build_next_nucleotide_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::chemical::rna;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Read in reference
	pose::Pose pose_start, pose_reference;

	bool const prepend_res = option[ prepend_residue] ;

	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_file( pose_start, *rsd_set, infile , core::import_pose::PDB_file);

	//This is a bit complicated...
	if ( option[in::file::native].user() ) {
		import_pose::pose_from_file( pose_reference, *rsd_set, option[in::file::native] , core::import_pose::PDB_file);
		if ( prepend_res ) {
			ResidueType rsd_type =  pose_reference.residue( 1 ).type();
			ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type ) ;
			pose_start.prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
		} else {
			remove_upper_terminus_type_from_pose_residue(pose_start,
				pose_start.size() );

			Size nres( pose_reference.size() );
			chemical::AA res_aa =  pose_reference.residue( nres ).aa();
			ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->get_representative_type_aa( res_aa )) ) ;
			pose_start.append_residue_by_bond( *new_rsd, true );

			add_upper_terminus_type_to_pose_residue(pose_start,
				pose_start.size() );

		}
	} else {
		pose_reference = pose_start;
	}

	pose_start.dump_pdb( "start.pdb" );

	pose::Pose pose = pose_start;
	Size const nres( pose.size() );

	// Now erase last residue, and rebuild last suite-to-suite, from scratch.
	Size const which_res = prepend_res ? 1 : nres;
	ResidueType rsd_type =  pose.residue( which_res   ).type();
	//pose.delete_polymer_residue( which_res );
	ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type ) ;

	utility::vector1< std::pair< std::string, std::string> > atom_pairs;
	if ( prepend_res ) {
		atom_pairs.push_back( std::make_pair( " O3'", " O3'" ) );
		atom_pairs.push_back( std::make_pair( " C3'", " C3'" ) );
		atom_pairs.push_back( std::make_pair( " C4'", " C4'" ) );
	} else {
		atom_pairs.push_back( std::make_pair( " P  ", " P  " ) );
		atom_pairs.push_back( std::make_pair( " O5'", " O5'" ) );
		atom_pairs.push_back( std::make_pair( " C5'", " C5'" ) );
	}

	pose.replace_residue( which_res, *new_rsd, atom_pairs );

	chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

	pose.set_torsion( TorsionID( which_res, id::BB, 1 ), rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 2 ), rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 3 ), rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 4 ), rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 5 ), rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
	pose.set_torsion( TorsionID( which_res, id::BB, 6 ), rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );

	pose.set_torsion( TorsionID( which_res, id::CHI, 1 ), rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
	pose.set_torsion( TorsionID( which_res, id::CHI, 2 ), rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center );
	pose.set_torsion( TorsionID( which_res, id::CHI, 3 ), rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center );
	pose.set_torsion( TorsionID( which_res, id::CHI, 4 ), 0.0 );

	pose.dump_pdb( "start2.pdb" );

	pose::Pose current_pose( pose );

	////////////////////////////////////////////////////
	protocols::viewer::add_conformation_viewer( current_pose.conformation(), "current", 400, 400 );

	///////////////////////////////////////////////////
	// Minimizer setup
	///////////////////////////////////////////////////

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	}
	scorefxn->set_weight( dihedral_constraint, 1.0 );
	scorefxn->set_weight( atom_pair_constraint, 1.0 );

	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	utility::vector1< utility::vector1 <Real > > backbone_rotamers;
	if ( prepend_res ) {
		get_backbone_rotamers( backbone_rotamers, ANY_PUCKER, NORTH );
	} else {
		get_backbone_rotamers( backbone_rotamers, NORTH, ANY_PUCKER );
	}

	Size const suite_res = prepend_res ? 1: (nres-1);

	bool quick = option[ quick_test ];

	Real rmsd( 0.0 );
	Real const RMSD_CUTOFF = option[ rmsd_cutoff ];

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;

	//mm.set_bb(  which_res, true );
	//mm.set_chi( which_res, true );
	//if (prepend_res) {
	// mm.set_bb(  which_res+1, true );
	//} else {
	// mm.set_bb(  which_res-1, true );
	//}

	//Just testing...
	mm.set_bb( true );
	mm.set_chi( true );

	///////////////////////////////////////////////////
	// Initial pose...
	Size count( 0 );
	current_pose = pose_reference;
	if ( !quick ) minimizer.run( current_pose, mm, *(scorefxn), options );
	(*scorefxn)(current_pose);
	rmsd = all_atom_rmsd( pose_reference, current_pose );
	std::string const tag( "S_"+lead_zero_string_of(count, 3) );
	RNA_SilentStruct s( current_pose, tag );
	s.add_energy( "rms", rmsd );
	silent_file_data.write_silent_struct(s, silent_file, false /*write score only*/);

	//////////////////////////////////////////////////////////
	for ( Size n = 1;  n <= backbone_rotamers.size(); n++ )  {

		count++;

		if ( prepend_res ) {
			apply_backbone_rotamer( pose, suite_res, backbone_rotamers[n],
				true /* pucker1 */, false /* pucker2 */);
		} else {
			apply_backbone_rotamer( pose, suite_res, backbone_rotamers[n],
				false /* pucker1 */, true /* pucker2 */);
		}

		current_pose = pose;

		if ( option[ vary_geometry] )   vary_geometry_RNA( pose, mm );

		//  mm.set_jump( true );

		rmsd = all_atom_rmsd( pose_reference, current_pose );
		std::cout << count << " All atom rmsd: " << rmsd  << std::endl;

		if ( rmsd > RMSD_CUTOFF ) continue;

		if ( !quick ) minimizer.run( current_pose, mm, *(scorefxn), options );

		(*scorefxn)(current_pose);

		rmsd = all_atom_rmsd( pose_reference, current_pose );

		std::string const tag( "S_"+lead_zero_string_of(count, 5) );
		RNA_SilentStruct s( current_pose, tag );
		s.add_energy( "rms", rmsd );

		silent_file_data.write_silent_struct(s, silent_file, false /*write score only*/);
		if ( !quick ) dump_pdb( current_pose, tag+".pdb" );

		//dump_pdb( pose, "nomin_"+tag+".pdb" );
	}
}

///////////////////////////////////////////////////////////////////////////////
typedef utility::vector1< GaussianParameter > GaussianParameterSet;

void
copy_rotamerized_torsions( pose::Pose & pose,
	pose::Pose & source_pose,
	Size const i,
	Size const rna_torsion_number,
	GaussianParameterSet const & gaussian_parameter_set )
{

	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	conformation::Residue rsd( source_pose.residue( i ) );

	//All this converting between torsion number <--> mainchain/chi is perhaps unnecessary?
	id::TorsionID  torsion_id( i, id::BB, rna_torsion_number );
	if ( rna_torsion_number > NUM_RNA_MAINCHAIN_TORSIONS ) {
		Size const chino( rna_torsion_number - NUM_RNA_MAINCHAIN_TORSIONS );
		torsion_id = id::TorsionID( i, id::CHI, chino  );
	}

	id::AtomID id1,id2,id3,id4;
	bool fail = source_pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
	if ( fail ) return;

	Real torsion_value( source_pose.torsion( torsion_id ) );

	// If there are multiple harmonic tethers available choose the "closest" one.
	assert( gaussian_parameter_set.size() > 0 );

	Real best_center( 0.0 ), best_sigma2( 10000000.0 );

	utility::vector1< Real > torsion_samples;
	torsion_samples.push_back( 0 );
	if ( option[ more_rotamers ] ) {
		torsion_samples.push_back( -1 );
		torsion_samples.push_back( 1 );
	}

	for ( Size n = 1; n <= gaussian_parameter_set.size(); n++ ) {

		Real const width = gaussian_parameter_set[n].width;
		Real const weight = 1.0/( width * width );

		for ( Size k = 1; k <= torsion_samples.size(); k++ ) {
			Real center = gaussian_parameter_set[n].center + torsion_samples[k] * width;

			Real deviation = numeric::principal_angle_degrees( torsion_value - center );
			Real const sigma2 = deviation * deviation * weight;

			if ( sigma2 < best_sigma2 ) {
				best_center = center;
				best_sigma2 = sigma2;
			}
		}
	}

	std::cout << "Residue " << i << "  torsion_num " << rna_torsion_number << "  " << F(8,3,torsion_value) << " " <<
		F(8,3,best_center) << "  ==> " << F(8,3,torsion_value-best_center) << std::endl;
	pose.set_torsion( torsion_id, best_center );
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
rotamerize_structure( pose::Pose & pose, pose::Pose & source_pose )
{
	using namespace core::scoring;
	using namespace core::chemical::rna;

	chemical::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

	Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );
	Size const nres( pose.size() );

	for ( Size i = 1; i <=nres; i++ ) {

		conformation::Residue const & rsd( source_pose.residue( i ) );
		copy_rotamerized_torsions( pose, source_pose, i, ALPHA, rna_fitted_torsion_info.gaussian_parameter_set_alpha());
		copy_rotamerized_torsions( pose, source_pose, i, BETA, rna_fitted_torsion_info.gaussian_parameter_set_beta());
		copy_rotamerized_torsions( pose, source_pose, i, GAMMA, rna_fitted_torsion_info.gaussian_parameter_set_gamma());

		Real const & delta( rsd.mainchain_torsion( DELTA ) );

		if ( delta <= DELTA_CUTOFF ) { // North, or 3'-endo sugar pucker, favored by RNA.
			copy_rotamerized_torsions( pose, source_pose, i, DELTA, rna_fitted_torsion_info.gaussian_parameter_set_delta_north());
		} else { // South, or 2'-endo sugar pucker.
			copy_rotamerized_torsions( pose, source_pose, i, DELTA, rna_fitted_torsion_info.gaussian_parameter_set_delta_south());
		}

		if ( delta <= DELTA_CUTOFF ) { // North, or 3'-endo sugar pucker, favored by RNA.
			copy_rotamerized_torsions( pose, source_pose, i, EPSILON, rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north());
		} else { // South, or 2'-endo sugar pucker.
			copy_rotamerized_torsions( pose, source_pose, i, EPSILON, rna_fitted_torsion_info.gaussian_parameter_set_epsilon_south());
		}

		Real next_alpha( -60.0 );
		if ( i < nres && source_pose.residue(i+1).is_RNA() ) next_alpha = source_pose.residue( i+1 ).mainchain_torsion( ALPHA );

		if ( next_alpha > -120.0 && next_alpha <= 0.0 ) { //default A-form, alpha sc-
			copy_rotamerized_torsions( pose, source_pose, i, ZETA, rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus());
		} else if ( next_alpha > 0.0 && next_alpha < 100.0 ) { // alpha sc+
			copy_rotamerized_torsions( pose, source_pose, i, ZETA, rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_plus());
		} else { // alpha ap
			copy_rotamerized_torsions( pose, source_pose, i, ZETA, rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_ap());
		}


		if ( delta <= DELTA_CUTOFF ) { // North, or 3'-endo sugar pucker, favored by RNA.
			copy_rotamerized_torsions( pose, source_pose, i, core::chemical::rna::CHI, rna_fitted_torsion_info.gaussian_parameter_set_chi_north());
		} else { // South, or 2'-endo sugar pucker.
			copy_rotamerized_torsions( pose, source_pose, i, core::chemical::rna::CHI, rna_fitted_torsion_info.gaussian_parameter_set_chi_south());
		}

		if ( delta <= DELTA_CUTOFF ) { // North, or 3'-endo sugar pucker, favored by RNA.
			copy_rotamerized_torsions( pose, source_pose, i, NU2, rna_fitted_torsion_info.gaussian_parameter_set_nu2_north());
		} else { // South, or 2'-endo sugar pucker.
			copy_rotamerized_torsions( pose, source_pose, i, NU2, rna_fitted_torsion_info.gaussian_parameter_set_nu2_south());
		}

		if ( delta <= DELTA_CUTOFF ) { // North, or 3'-endo sugar pucker, favored by RNA.
			copy_rotamerized_torsions( pose, source_pose, i, NU1, rna_fitted_torsion_info.gaussian_parameter_set_nu1_north());
		} else { // South, or 2'-endo sugar pucker.
			copy_rotamerized_torsions( pose, source_pose, i, NU1, rna_fitted_torsion_info.gaussian_parameter_set_nu1_south());
		}

	}

}


///////////////////////////////////////////////////////////////////////////////
void
rotamerize_rna_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace protocols::farna;
	using namespace core::chemical::rna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	pose::Pose source_pose,pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_file( source_pose, *rsd_set, infile , core::import_pose::PDB_file);

	std::string sequence = source_pose.sequence();

	core::pose::make_pose_from_sequence(
		pose,
		sequence,
		*rsd_set );
	dump_pdb( pose, "start.pdb" );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	///////////////////////////////////////////////////////
	rotamerize_structure( pose, source_pose );

	pose.dump_pdb( "rotamerized.pdb");
	std::cout << "All-atom RMSD " <<   all_atom_rmsd( source_pose, pose ) << std::endl;


	//////////////////////////////////////////////////
	// Stay near this rotamer
	ConstraintSetOP cst_set( new ConstraintSet() );
	for ( Size i = 1; i < pose.size(); i++ ) {
		create_dihedral_constraint( cst_set, pose, TorsionID( 1, id::BB, 4 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 1, id::CHI, 1 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 1, id::CHI, 2 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 1, id::CHI, 3 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 1, id::BB, 5 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 1, id::BB, 6 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 2, id::BB, 1 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 2, id::BB, 2 ) );
		create_dihedral_constraint( cst_set, pose, TorsionID( 2, id::BB, 3 ) );
	}


	//////////////////////////////////////////////////
	// base pairing constraints.
	Real const WC_distance( 1.9 );
	Real const distance_stddev( 0.25 ); //Hmm. Maybe try linear instead?
	FuncOP const distance_func( new HarmonicFunc( WC_distance, distance_stddev ) );

	utility::vector1< std::string > atom_ids1, atom_ids2;
	Size i( 1 );
	Size j( pose.size() );
	get_watson_crick_base_pair_atoms( pose.residue_type(i), pose.residue_type(j), atom_ids1, atom_ids2 );

	for ( Size p = 1; p <= atom_ids1.size(); p++ ) {

		Size const atom1 = pose.residue_type(i).atom_index( atom_ids1[p] ) ;
		Size const atom2 = pose.residue_type(j).atom_index( atom_ids2[p] ) ;

		std::cout << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue_type(i).name1() << I(3,i) << " <-->  " <<
			pose.residue_type(j).name1() << I(3,j) << "   " <<
			atom_ids1[p] << " <--> " <<
			atom_ids2[p] << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;

		cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint(
			id::AtomID(atom1,i),
			id::AtomID(atom2,j),
			distance_func ) ) ) );
	}

	//////////////////////////////////////////////////
	pose.constraint_set( cst_set );

	///////////////////////////////////////////////////
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	}
	scorefxn->set_weight( dihedral_constraint, 1.0 );
	scorefxn->set_weight( atom_pair_constraint, 1.0 );


	///2'-OH
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !pose.residue_type(i).is_RNA() ) continue;
		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
	}
	//  pack::pack_rotamers( pose, *packer_scorefxn, task);
	pack::rotamer_trials( pose, *scorefxn, task);


	//////////////////////////////////////////////////
	// minimize...
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;
	mm.set_bb(  true );
	mm.set_chi( true );
	mm.set_jump( true );

	if ( option[ vary_geometry] ) {
		vary_geometry_RNA( pose, mm );
	}
	scorefxn->set_weight( dof_constraint, 1.0 );

	minimizer.run( pose, mm, *(scorefxn), options );

	scorefxn->show( std::cout, pose );
	std::cout << "All-atom RMSD " <<   all_atom_rmsd( source_pose, pose ) << std::endl;

	pose.dump_pdb( "minimized.pdb" );

}

//////////////////////////////////////////////////////////////////////
void
calc_rmsd_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	pose::Pose native_pose;
	std::string native_pdb_file  = option[ in::file::native ];
	import_pose::pose_from_file( native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( native_pose );

	Real rmsd( 0.0 );
	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );

		rmsd = all_atom_rmsd( native_pose, pose );
		std::cout << " RMSD " << pdb_file << " " << rmsd << std::endl;
	}

}

/////////////////////////////////////////////////////////////////
void
output_sugar_geometry_parameters(
	pose::Pose & pose,
	utility::vector1< utility::vector1< std::string> > & bond_length_atoms,
	utility::vector1< utility::vector1< std::string> > & bond_angle_atoms,
	// utility::vector1< utility::vector1< std::string> > & bond_torsions
	utility::io::ozstream & out
) {

	using namespace core::conformation;
	using namespace core::chemical::rna;
	Size const nres ( pose.size() );

	for ( Size i = 1; i <= nres; i++ ) {
		Residue const & rsd( pose.residue( i ) );

		out << rsd.mainchain_torsion( 4 );

		for ( Size n = 1; n <= bond_length_atoms.size(); n++ ) {
			std::string a1 = bond_length_atoms[n][1];
			std::string a2 = bond_length_atoms[n][2];
			if ( a1 == "BASE" ) a1 = first_base_atom( rsd.type() );
			if ( a2 == "BASE" ) a2 = first_base_atom( rsd.type() );
			//std::cout << "LENGTH " << a1 << " " << a2 << " " << ( rsd.xyz( a1 ) - rsd.xyz( a2 ) ).length() << std::endl;
			out << ' ' << ( rsd.xyz( a1 ) - rsd.xyz( a2 ) ).length();
		}

		for ( Size n = 1; n <= bond_angle_atoms.size(); n++ ) {
			std::string a1 = bond_angle_atoms[n][1];
			std::string a2 = bond_angle_atoms[n][2];
			std::string a3 = bond_angle_atoms[n][3];
			if ( a1 == "BASE" ) a1 = first_base_atom( rsd.type() );
			if ( a2 == "BASE" ) a2 = first_base_atom( rsd.type() );
			if ( a3 == "BASE" ) a3 = first_base_atom( rsd.type() );
			//std::cout << "LENGTH " << a1 << " " << a2 << " " << ( rsd.xyz( a1 ) - rsd.xyz( a2 ) ).length() << std::endl;
			out << ' ' << numeric::conversions::degrees( angle_radians( rsd.xyz( a1 ), rsd.xyz( a2 ), rsd.xyz( a3 ) ) );
		}

		//////////////////////////////////////////////
		Real const a = (rsd.xyz( " O4'" ) - rsd.xyz( " C4'" )).length();
		Real const b = (rsd.xyz( " C3'" ) - rsd.xyz( " C4'" )).length();
		Real const c = (rsd.xyz( "VO4'" ) - rsd.xyz( " C3'" )).length();
		Real const theta = numeric::conversions::degrees( std::acos( ( a*a + b*b - c*c )/(2 * a *b ) ) );
		out << ' ' << theta;

		out << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////
void
put_it_in_list(
	utility::vector1< utility::vector1< std::string> > & my_list,
	std::string const  a1, std::string const  a2 ){
	utility::vector1< std::string > vec;
	vec.push_back( a1 );
	vec.push_back( a2 );
	my_list.push_back( vec );
}
void
put_it_in_list(
	utility::vector1< utility::vector1< std::string> > & my_list,
	std::string const  a1, std::string const  a2, std::string const a3 ){
	utility::vector1< std::string > vec;
	vec.push_back( a1 );
	vec.push_back( a2 );
	vec.push_back( a3 );
	my_list.push_back( vec );
}
void
put_it_in_list(
	utility::vector1< utility::vector1< std::string> > & my_list,
	std::string const  a1, std::string const  a2,
	std::string const  a3, std::string const  a4 ){
	utility::vector1< std::string > vec;
	vec.push_back( a1 );
	vec.push_back( a2 );
	vec.push_back( a3 );
	vec.push_back( a4 );
	my_list.push_back( vec );
}

/////////////////////////////////////////////////////////////////////
void
fill_sugar_atom_list( utility::vector1< std::string> & sugar_atom_list )
{
	sugar_atom_list.push_back( " C2'" );
	sugar_atom_list.push_back( " O2'" );
	sugar_atom_list.push_back( " C1'" );
	sugar_atom_list.push_back( " O4'" );
	sugar_atom_list.push_back( "BASE" );
}

/////////////////////////////////////////////////////////////////////
void
fill_bond_atoms(
	utility::vector1< utility::vector1< std::string> > & bond_length_atoms ,
	utility::vector1< utility::vector1< std::string> > & bond_angle_atoms
)
{
	put_it_in_list( bond_length_atoms, " C3'", " C4'" );
	put_it_in_list( bond_length_atoms, " C2'", " C3'" );
	put_it_in_list( bond_length_atoms, " C1'", " C2'" );
	put_it_in_list( bond_length_atoms, " C4'", " O4'" );
	put_it_in_list( bond_length_atoms, " O4'", " C1'" );
	put_it_in_list( bond_length_atoms, " C1'", "BASE" );

	put_it_in_list( bond_angle_atoms, " C3'", " C4'", " O4'" );
	put_it_in_list( bond_angle_atoms, " C2'", " C3'", " C4'" );
	put_it_in_list( bond_angle_atoms, " C1'", " C2'", " C3'" );
	put_it_in_list( bond_angle_atoms, " C2'", " C1'", " O4'" );
	put_it_in_list( bond_angle_atoms, " C4'", " O4'", " C1'" );
	put_it_in_list( bond_angle_atoms, " C2'", " C1'", "BASE" );
	put_it_in_list( bond_angle_atoms, " O4'", " C1'", "BASE" );

}

////////////////////////////////////////////////////////////////
void
output_sugar_internal_dof( pose::Pose & pose, utility::vector1< std::string > const & sugar_atom_list, utility::io::ozstream & out )
{
	using namespace core::id;
	using namespace numeric::conversions;

	for ( Size i = 1; i <= pose.size(); i++ ) {
		core::conformation::Residue const & rsd( pose.residue( i ) );

		std::cout << " DELTA " << rsd.mainchain_torsion( 4 ) << std::endl;
		out << rsd.mainchain_torsion( 4 );

		for ( Size n = 1; n <= sugar_atom_list.size(); n++ ) {

			std::string atom_name(  sugar_atom_list[ n ] );
			if ( atom_name == "BASE" ) atom_name = core::chemical::rna::first_base_atom( pose.residue_type(i) );

			Size const j = rsd.atom_index( atom_name );

			core::kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID(j,i) ).get_self_ptr() );

			std::cout << A( 5, rsd.atom_name( j )) << " " <<
				F(11,6, degrees( pose.atom_tree().dof( DOF_ID( current_atom->id(), id::PHI ) ) ) )  << " " <<
				F(11,6, degrees(  pose.atom_tree().dof( DOF_ID( current_atom->id(), id::THETA ) ) ) ) << " " <<
				F(11,6, pose.atom_tree().dof( DOF_ID( current_atom->id(), id::D ) ) )    << "  " << std::endl;

			out << " " <<
				F(11,6, degrees( pose.atom_tree().dof( DOF_ID( current_atom->id(), id::PHI ) ) ) )  << " " <<
				F(11,6, degrees(  pose.atom_tree().dof( DOF_ID( current_atom->id(), id::THETA ) ) ) ) << " " <<
				F(11,6, pose.atom_tree().dof( DOF_ID( current_atom->id(), id::D ) ) );

		}
		out << std::endl;
	}

}


//////////////////////////////////////////////////////////////////////
void
replace_torsion_angles( pose::Pose & extended_pose, pose::Pose & fixed_pose, pose::Pose & pose )
{

	using namespace core::conformation;
	using namespace core::id;

	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j=1; j <= core::chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
			id::TorsionID my_ID( i, id::BB, j );
			extended_pose.set_torsion( my_ID, pose.torsion( my_ID ) );
		}
		for ( Size j=1; j <= core::chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {
			id::TorsionID my_ID( i, id::CHI, j );
			extended_pose.set_torsion( my_ID, pose.torsion( my_ID ) );
		}
	}

	fixed_pose = extended_pose;

	for ( Size i = 1; i <= pose.size(); i++ ) {
		Residue rsd = pose.residue( i );
		kinematics::Stub const input_stub( rsd.xyz( " C3'" ), rsd.xyz( " C3'" ), rsd.xyz( " C4'" ), rsd.xyz( " C5'" ) );

		Residue rsd_fixed = fixed_pose.residue( i );
		kinematics::Stub const input_stub_fixed( rsd_fixed.xyz( " C3'" ), rsd_fixed.xyz( " C3'" ), rsd_fixed.xyz( " C4'" ), rsd_fixed.xyz( " C5'" ) );

		utility::vector1< std::string > my_atoms;
		my_atoms.push_back( " C2'" );
		my_atoms.push_back( " C1'" );
		my_atoms.push_back( " O4'" );
		//my_atoms.push_back( " C4'" );

		utility::vector1< Vector > start_vectors;
		utility::vector1< utility::vector1< Real > > new_dof_sets;

		//What DOFS do I need to get the ring atoms where I want them?
		for ( Size n = 1; n <= my_atoms.size(); n++  ) {
			Size const j = rsd.atom_index( my_atoms[ n ] );
			Vector v1 = input_stub.global2local( rsd.xyz( my_atoms[ n ] ) );
			Vector v2 = input_stub_fixed.local2global( v1 );
			Vector v = rsd_fixed.xyz( my_atoms[ n ]  );
			start_vectors.push_back( v );
			fixed_pose.set_xyz( id::AtomID( j,i ), v2 );

			utility::vector1< Real > dof_set;
			dof_set.push_back( pose.atom_tree().dof( DOF_ID( AtomID( j,i) , id::D ) ) );
			dof_set.push_back( pose.atom_tree().dof( DOF_ID( AtomID( j,i) , id::PHI ) ) );
			dof_set.push_back( pose.atom_tree().dof( DOF_ID( AtomID( j,i) , id::THETA ) ) );
			new_dof_sets.push_back( dof_set );

		}


		// Now put the ring atoms in the desired spots, but by changing internal DOFS --
		// rest of the atoms (e.g., in base and 2'-OH) will scoot around as well, preserving
		// ideal bond lengths and angles.
		for ( Size n = 1; n <= my_atoms.size(); n++  ) {
			Size const j = rsd.atom_index( my_atoms[ n ] );
			fixed_pose.set_xyz( id::AtomID( j,i ), start_vectors[n] );
		}
		for ( Size n = 1; n <= my_atoms.size(); n++  ) {
			Size const j = rsd.atom_index( my_atoms[ n ] );
			fixed_pose.set_dof(  DOF_ID( AtomID( j,i) , id::D ), new_dof_sets[n][1] );
			fixed_pose.set_dof(  DOF_ID( AtomID( j,i) , id::PHI ), new_dof_sets[n][2] );
			fixed_pose.set_dof(  DOF_ID( AtomID( j,i) , id::THETA ), new_dof_sets[n][3] );
		}

	}

}

////////////////////////////////////////////////////////
void
fix_sugar_bond_angles_EMPIRICAL( pose::Pose & pose )
{
	using namespace core::id;

	Real theta( 0.0 ), phi( 0.0 );

	for ( Size i = 1; i <= pose.size(); i++ ) {

		core::conformation::Residue const & rsd( pose.residue( i ) );
		Real const delta = rsd.mainchain_torsion( 4 );

		{
			std::string const atom_name = " C2'";
			Size const j = rsd.atom_index( atom_name );
			core::kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID(j,i) ).get_self_ptr() );
			if ( delta < 100.0 ) {
				theta  = -0.138 * delta + 89.4;
			} else {
				theta  =  0.0423 * delta + 71.4;
			}
			pose.set_dof( DOF_ID( current_atom->id(), id::THETA )  , numeric::conversions::radians( theta ) );
		}


		{
			std::string const atom_name = " O4'";
			Size const j = rsd.atom_index( atom_name );
			core::kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID(j,i) ).get_self_ptr() );
			if ( delta < 100.0 ) {
				theta  =  0.132 * delta + 59.5;
				phi = 0.0118 * delta - 118.0;
			} else {
				theta  =  0.0168 * delta + 68.2;
				phi = 0.0757 * delta - 131.0;
			}
			pose.set_dof( DOF_ID( current_atom->id(), id::THETA )  , numeric::conversions::radians( theta ) );
			pose.set_dof( DOF_ID( current_atom->id(), id::PHI )  , numeric::conversions::radians( phi ) );
		}


	}

}

////////////////////////////////////////////////////////
void
fix_sugar_bond_angles_CLOSE_BOND( pose::Pose & pose )
{
	using namespace core::id;

	for ( Size i = 1; i <= pose.size(); i++ ) {

		core::conformation::Residue const & rsd( pose.residue( i ) );

		Real const a = (rsd.xyz( " O4'" ) - rsd.xyz( " C4'" )).length();
		Real const b = (rsd.xyz( " C3'" ) - rsd.xyz( " C4'" )).length();
		Real const c = (rsd.xyz( "VO4'" ) - rsd.xyz( " C3'" )).length();
		Real const theta = std::acos( ( c*c - a*a - b*b )/(2 * a *b ) );

		{
			kinematics::Stub const input_stub( rsd.xyz( " C4'"), rsd.xyz( " C4'"),
				rsd.xyz( " C3'"), rsd.xyz( " O4'") );

			Vector v1 = input_stub.global2local( rsd.xyz( " O4'" ) );
			//   std::cout << "VEC " << v1.x() << " " << v1.y() << " " << v1.z() << std::endl;

			std::string atom_name = " O4'";
			Size const j = rsd.atom_index( atom_name );
			pose.set_xyz( id::AtomID(j,i), input_stub.spherical( 0.0, theta, a ) );

			//v1 = input_stub.global2local( rsd.xyz( " O4'" ) );
			//   std::cout << "VEC_NEW " << v1.x() << " " << v1.y() << " " << v1.z() << std::endl;

		}

		//  std::cout << "LENGTH_O4'_VO4': " << ( pose.residue(i).xyz( "VO4'" ) - pose.residue(i).xyz(" C3'") ).length() << " " << ( pose.residue(i).xyz( " O4'" ) - pose.residue(i).xyz(" C3'") ).length() << " " << std::endl;


		std::string atom_name = " C2'";
		Size j = rsd.atom_index( atom_name );
		{
			kinematics::Stub const input_stub( rsd.xyz( " C3'"), rsd.xyz( " C3'"),
				rsd.xyz( " C4'"), rsd.xyz( " C2'") );


			Vector const v1 = input_stub.global2local( pose.residue(i).xyz( " O4'" ) );
			Vector const v2 = input_stub.global2local( pose.residue(i).xyz( "VO4'" ) );
			//Real const theta_difference = std::atan2( v2.y(), v2.x() )  - std::atan2( v1.y(), v1.x() ) ;

			//   std::cout << "BLAH  " << input_stub.global2local( pose.residue(i).xyz( " O4'" ) ).x() << " " <<
			//    input_stub.global2local( pose.residue(i).xyz( "VO4'" ) ).x() << std::endl;

			Real const d = sqrt( v2.x() * v2.x() + v2.y()*v2.y() );
			Real const alpha = std::atan2( v2.y(), v2.x() );
			Real const theta_difference = std::acos( v1.x()/d ) - alpha;

			//   std::cout << "THETA DIFFERENCE " << alpha << " " << theta_difference << std::endl;
			id::DOF_ID dof_id( DOF_ID( id::AtomID( j, i), THETA ) );
			pose.set_dof(  dof_id, pose.dof( dof_id ) + theta_difference );

			//   std::cout << "BLAH2 " << input_stub.global2local( pose.residue(i).xyz( " O4'" ) ).x() << " " <<
			//    input_stub.global2local( pose.residue(i).xyz( "VO4'" ) ).x() << std::endl;

		}

		{
			core::kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID(j,i) ).get_self_ptr() );
			core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
			core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );
			kinematics::Stub const input_stub( input_stub_atom1->xyz(), input_stub_atom1->xyz(), input_stub_atom2->xyz(), input_stub_atom3->xyz());

			Vector const v1 = input_stub.global2local( pose.residue(i).xyz( " O4'" ) );
			Vector const v2 = input_stub.global2local( pose.residue(i).xyz( "VO4'" ) );
			Real const phi_difference = std::atan2( v2.z(), v2.y() )  - std::atan2( v1.z(), v1.y() ) ;
			//std::cout << "PHI_DIFFERENCE " << phi_difference << std::endl;
			id::DOF_ID dof_id( DOF_ID( id::AtomID( j, i), PHI ) );
			Real const phi = pose.dof( dof_id );
			//std::cout << "PHI " << phi << std::endl;
			pose.set_dof(  dof_id, phi - phi_difference );
		}

		//{
		//Vector const v1 = input_stub.global2local( pose.residue(i).xyz( " O4'" ) );
		//Vector const v2 = input_stub.global2local( pose.residue(i).xyz( "VO4'" ) );
		//Real const phi_difference = std::atan2( v2.z(), v2.y() )  - std::atan2( v1.z(), v1.y() ) ;
		//std::cout << "PHI_DIFFERENCE_NEW " << phi_difference << std::endl;
		//}

		//  std::cout << "LENGTH: " << ( pose.residue(i).xyz( "VO4'" ) - pose.residue(i).xyz( " O4'" ) ).length() << std::endl;


	}

}


//////////////////////////////////////////////////////////////////////
void
sugar_geometry_RNA_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );


	utility::vector1< utility::vector1< std::string> > bond_length_atoms;
	utility::vector1< utility::vector1< std::string> > bond_angle_atoms;
	fill_bond_atoms( bond_length_atoms, bond_angle_atoms );

	utility::vector1< std::string > sugar_atom_list;
	fill_sugar_atom_list( sugar_atom_list );

	utility::io::ozstream out( "sugar_internal_dof.txt" );

	utility::io::ozstream orig_out( "sugar_geometry_orig.txt" );
	utility::io::ozstream ideal_out( "sugar_geometry_ideal.txt" );
	utility::io::ozstream fix_out( "sugar_geometry_fix.txt" );

	pose::Pose extended_pose, pose, fixed_pose;

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];

		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );


		output_sugar_internal_dof( pose, sugar_atom_list, out ) ;
		output_sugar_geometry_parameters( pose, bond_length_atoms, bond_angle_atoms, orig_out );

		make_pose_from_sequence( extended_pose, pose.sequence(), *rsd_set );
		replace_torsion_angles( extended_pose, fixed_pose, pose );
		output_sugar_geometry_parameters( extended_pose, bond_length_atoms, bond_angle_atoms, ideal_out );

		//fix_sugar_bond_angles( extended_pose );
		output_sugar_geometry_parameters( fixed_pose, bond_length_atoms, bond_angle_atoms, fix_out );

	}

	std::cout << "Numerical values of dofs in sugar_internal_dof.txt " << std::endl;
}

//////////////////////////////////////////////////////////////////////
void
sugar_frag_RNA_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	pose::Pose pose;
	make_pose_from_sequence( pose, "au", *rsd_set );
	pose.dump_pdb( "start.pdb" );
	kinematics::FoldTree f( pose.size() );
	f.new_jump(1,2,1);
	f.set_jump_atoms( 1,
		core::chemical::rna::chi1_torsion_atom( pose.residue_type( 1 ) ),
		core::chemical::rna::chi1_torsion_atom( pose.residue_type( 2 ) )   );
	pose.fold_tree( f );
	std::cout << f << std::endl;

	{
		pose::Pose pose_temp = pose;
		for ( Size  n = 1; n <= pose.size(); n++ )  {
			pose::rna::apply_ideal_c2endo_sugar_coords( pose_temp, n /*, true  */ );
		}
		pose_temp.dump_pdb( "final_WORKS.pdb" );
	}
	{
		pose::Pose pose_temp = pose;
		for ( Size  n = 1; n <= pose.size(); n++ )  {
			pose::rna::apply_ideal_c2endo_sugar_coords( pose_temp, n /*, false */ );
		}
		pose_temp.dump_pdb( "final_NOWORKS.pdb" );
	}
}

////////////////////////////////////////////////////////
void
color_by_geom_sol_RNA_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	utility::vector1< std::string > pdb_files( option[ in::file::s ]() );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	geometric_solvation::GeometricSolEnergyEvaluator geometric_sol_energy_method( scorefxn->energy_method_options() );

	Pose pose;

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );

		PDBInfoOP pdb_info( new PDBInfo( pose, true ) );

		(*scorefxn)( pose );

		for ( Size i = 1; i <= pose.size(); i++ ) {
			for ( Size j = 1; j <= pose.residue( i ) .natoms(); j++ ) {
				pdb_info->temperature( i, j, geometric_sol_energy_method.eval_atom_energy( AtomID( j, i ), pose ) );
			}
		}

		pose.pdb_info( pdb_info );

		Size pos( pdb_file.find( ".pdb" ) );
		std::string outfile( pdb_file );
		outfile.replace( pos, 4, "_color_geomsol.pdb" );
		pose.dump_pdb( outfile );
	}

}
////////////////////////////////////////////////////////
void
color_by_lj_base_RNA_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::farna;
	using namespace core::scoring::etable;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( option[ rsd_type_set ] );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );

	EtableOP etable_ptr( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
		EtableOptions() ) );
	// core::chemical::rna::RNA_LJ_BaseEnergy rna_lj_base_energy( *etable_ptr );

	Pose pose;

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );

		PDBInfoOP pdb_info( new PDBInfo( pose, true ) );

		(*scorefxn)( pose );

		for ( Size i = 1; i <= pose.size(); i++ ) {
			for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ) {
				Real const score = 0.0; //rna_lj_base_energy.eval_atom_energy( AtomID( j, i ), pose );
				std::cout << i <<  " " << j << " " <<  score << std::endl;
				pdb_info->temperature( i, j, score );
			}
		}

		pose.pdb_info( pdb_info );

		Size pos( pdb_file.find( ".pdb" ) );
		std::string outfile( pdb_file );
		outfile.replace( pos, 4, "_color_lj_base.pdb" );
		pose.dump_pdb( outfile );
	}

}

////////////////////////////////////////////////////////
void
rna_stats_test()
{
	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );

	std::string const silent_file = option[ out::file::silent  ]();
	protocols::farna::RNA_DeNovoProtocol rna_de_novo_protocol;

	pose::PoseOP native_pose_OP( new pose::Pose );
	pose::Pose & native_pose = *native_pose_OP;

	if ( option[ in::file::native ].active() ) {
		std::string native_pdb_file  = option[ in::file::native ];
		import_pose::pose_from_file( native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( native_pose );
		rna_de_novo_protocol.set_native_pose( native_pose_OP );
	}

	Pose pose;

	for ( Size n = 1; n <= pdb_files.size(); n++ ) {
		std::string const pdb_file = pdb_files[n];
		import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
		protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose );

		utility::vector1< core::pose::rna::BasePair > base_pair_list;
		utility::vector1< bool > is_bulged;
		core::pose::rna::classify_base_pairs( pose, base_pair_list, is_bulged );

		rna_de_novo_protocol.output_to_silent_file( pose, silent_file, pdb_files[n], true /*score_only*/ );

	}

}


///////////////////////////////////////////////////////
void
files_for_openMM_test(){

	// Read in pdb.
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::farna;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );
	std::string const pdb_file = pdb_files[1];
	Pose pose_start;
	import_pose::pose_from_file( pose_start, *rsd_set, pdb_file , core::import_pose::PDB_file);
	protocols::farna::ensure_phosphate_nomenclature_matches_mini( pose_start );

	Pose pose;
	make_pose_from_sequence( pose, pose_start.sequence(), *rsd_set );
	//apply_ideal_coordinates_for_alternative_pucker( pose_start, pose );


	pose_start.dump_pdb( "rna_start.pdb" );
	pose.dump_pdb( "rna_reference.pdb" );

	std::map< AtomID, Size > global_index_map;

	// xyz coordinates.
	utility::io::ozstream out1( "xyz.txt" );
	Size count( 0 );
	for ( Size i = 1; i <= pose_start.size(); i++ ) {
		for ( Size j = 1; j <= pose_start.residue_type(i).natoms(); j++ ) {
			numeric::xyzVector< Real > const & v = pose_start.residue(i).xyz(j);
			out1 << 0.1*v(1) << ' ' << 0.1*v(2) << ' ' << 0.1*v(3) << std::endl;

			global_index_map[ AtomID(j,i) ] = count;
			count++;
		}
	}
	out1.close();

	// bonds
	utility::io::ozstream out2( "bonds.txt" );
	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j = 1; j <= pose.residue_type(i).natoms(); j++ ) {

			AtomID atom_id1( j,i);
			utility::vector1< AtomID > atom_ids2 = pose.conformation().bonded_neighbor_all_res( atom_id1 );
			for ( Size n = 1; n <= atom_ids2.size(); n++ ) {
				AtomID const & atom_id2 = atom_ids2[ n ] ;
				if ( global_index_map[ atom_id1 ] > global_index_map[ atom_id2 ] ) {
					Real const atom_atom_distance = ( pose.xyz( atom_id1 ) - pose.xyz( atom_id2 ) ).length();
					out2 << global_index_map[ atom_id1 ] << ' ' <<  global_index_map[ atom_id2 ] << ' ' << 0.1 * atom_atom_distance << std::endl;
				}
			}

		}
	}
	out2.close();

	// angles
	utility::io::ozstream out3( "angles.txt" );
	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j = 1; j <= pose.residue_type(i).natoms(); j++ ) {

			AtomID atom_id1( j,i);
			utility::vector1< AtomID > atom_ids2 = pose.conformation().bonded_neighbor_all_res( atom_id1 );

			for ( Size n = 1; n <= atom_ids2.size(); n++ ) {
				AtomID const & atom_id2 = atom_ids2[ n ] ;

				for ( Size q = (n+1); q <= atom_ids2.size(); q++ ) {
					AtomID const & atom_id3 = atom_ids2[ q ] ;
					Real const angle = angle_radians( pose.xyz( atom_id2 ), pose.xyz( atom_id1 ), pose.xyz( atom_id3 ) );
					out3 << global_index_map[ atom_id2] << ' ' << global_index_map[ atom_id1 ] << ' ' << global_index_map[ atom_id3 ] << ' ' << angle << std::endl;
				}
			}

		}
	}
	out3.close();

	utility::vector1< Size > resnum_for_donor, resnum_for_acceptor;

	//donors
	utility::io::ozstream out4( "donors.txt" );
	for ( Size i = 1; i <= pose.size(); i++ ) {
		AtomIndices const & hpos_polar = pose.residue_type(i).Hpos_polar();
		for ( Size n = 1; n <= hpos_polar.size(); n++ ) {
			AtomID atom_id( hpos_polar[n], i );
			AtomID atom_id_base( pose.residue_type(i).atom_base( hpos_polar[n] ), i );
			out4 << global_index_map[ atom_id ] << ' ' << global_index_map[ atom_id_base ] << std::endl;
			resnum_for_donor.push_back( i );
		}
	}
	out4.close();


	//donors
	utility::io::ozstream out5( "acceptors.txt" );
	for ( Size i = 1; i <= pose.size(); i++ ) {
		AtomIndices const & acceptors = pose.residue_type(i).accpt_pos();
		for ( Size n = 1; n <= acceptors.size(); n++ ) {
			AtomID atom_id( acceptors[n], i );
			AtomID atom_id_base( pose.residue_type(i).atom_base( acceptors[n] ), i );
			out5 << global_index_map[ atom_id ] << ' ' << global_index_map[ atom_id_base ] << std::endl;
			resnum_for_acceptor.push_back( i );
		}
	}
	out5.close();

	//donors
	utility::io::ozstream out6( "donor_acceptor_exclude.txt" );
	for ( Size i = 1; i <= resnum_for_donor.size(); i++ ) {
		for ( Size j = 1; j <= resnum_for_acceptor.size(); j++ ) {
			if ( resnum_for_donor[i] == resnum_for_acceptor[j] ) {
				out6 << i-1 << ' ' << j-1 << std::endl;
			}
		}
	}
	out6.close();


}

/////////////////////////////////////////////////
void
print_secstruct_test(){
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	ResidueTypeSetCOP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);

	utility::vector1< std::pair<Size, Size> > base_pairing_list;
	protocols::farna::get_base_pairing_list( pose, base_pairing_list );

	std::cout << "WATSON-CRICK BASE PAIRS: " << std::endl;
	for ( Size n = 1; n <= base_pairing_list.size(); n++ ) {
		std::cout << base_pairing_list[n].first << ' ' << base_pairing_list[n].second << std::endl;
	}

}

/////////////////////////////////////////////////
void
print_all_torsions_test(){
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace core::id;
	using namespace chemical::rna;

	ResidueTypeSetCOP rsd_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);

	std::cout << "  N    ALPHA     BETA    GAMMA    DELTA  EPSILON     ZETA         CHI" << std::endl;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		std::cout << I( 3, n )
			<< ' ' << F( 8, 3, pose.torsion( TorsionID( n, id::BB, ALPHA ) ) )
			<< ' ' << F( 8, 3, pose.torsion( TorsionID( n, id::BB, BETA ) ) )
			<< ' ' << F( 8, 3, pose.torsion( TorsionID( n, id::BB, GAMMA ) ) )
			<< ' ' << F( 8, 3, pose.torsion( TorsionID( n, id::BB, DELTA ) ) )
			<< ' ' << F( 8, 3, pose.torsion( TorsionID( n, id::BB, EPSILON ) ) )
			<< ' ' << F( 8, 3, pose.torsion( TorsionID( n, id::BB, ZETA ) ) )
			<< "   " << F( 8, 3, pose.torsion( TorsionID( n, id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ) ) )
			<< std::endl;
	}

}


///////////////////////////////////////////////////////////////////////////////
void
rna_denovo_test()
{
	std::cout << std::endl;
	std::cout << "RNA de novo structure prediction is now in       " << std::endl;
	std::cout << " the rna_denovo application. " << std::endl;
	std::cout << "  [code is in src/apps/public/stepwise/rna/rna_denovo.cc] " << std::endl;
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
rna_design_test()
{
	std::cout << std::endl;
	std::cout << "RNA design is now in       " << std::endl;
	std::cout << " the rna_design application. " << std::endl;
	std::cout << "  [code is in src/apps/public/stepwise/rna/rna_design.cc] " << std::endl;
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions_test()
{
	std::cout << std::endl;
	std::cout << "RNA vall torsions code is now in       " << std::endl;
	std::cout << " the rna_database application. " << std::endl;
	std::cout << "  [code is in src/apps/public/stepwise/rna/rna_database.cc] " << std::endl;
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
create_bp_jump_database_test()
{
	std::cout << std::endl;
	std::cout << "Jump Database code is now in       " << std::endl;
	std::cout << " the rna_database application. " << std::endl;
	std::cout << "  [code is in src/apps/public/stepwise/rna/rna_database.cc] " << std::endl;
	std::cout << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
void
extract_pdbs_test()
{
	std::cout << std::endl;
	std::cout << "RNA silent extractions code is now in       " << std::endl;
	std::cout << " the rna_database application. " << std::endl;
	std::cout << "  [code is in src/apps/public/stepwise/rna/rna_extract.cc] " << std::endl;
	std::cout << std::endl;
}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	if ( option[ create_vall_torsions ] ) {
		create_rna_vall_torsions_test();
	} else if ( option[ create_jump_database ] ) {
		create_bp_jump_database_test();
	} else if ( option[ icoord_test ] ) {
		figure_out_icoord_test();
	} else if ( option[ extract ] ) {
		extract_pdbs_test();
	} else if ( option[ print_internal_coord ] ) {
		print_internal_coord_test();
	} else if ( option[ fullatom_score_test ] ) {
		rna_fullatom_score_test();
	} else if ( option[ o2prime_test ] ) {
		rna_o2prime_test();
	} else if ( option[ fullatom_multiscore ] ) {
		rna_fullatom_multiscore_test();
	} else if ( option[ fullatom_minimize_silent ] ) {
		rna_fullatom_minimize_silent_test();
	} else if ( option[ lores_score ] ) {
		rna_lores_score_test();
	} else if ( option[ lores_score_silent ] ) {
		rna_lores_score_silent_test();
	} else if ( option[ rna_design ] ) {
		rna_design_test();
	} else if ( option[ rna_design_gap ] ) {
		rna_design_gap_test();
	} else if ( option[ rna_idealize ] ) {
		rna_idealize_test();
	} else if ( option[ rna_torsion_check ] ) {
		rna_torsion_check_test();
	} else if ( option[ rna_stats ] ) {
		rna_stats_test();
	} else if ( option[ close_chainbreaks_test ] ) {
		rna_close_chainbreaks_test();
	} else if ( option[ create_benchmark ] ) {
		create_rna_benchmark_test();
	} else if ( option[ chain_closure_test ] ) {
		rna_chain_closure_test();
	} else if ( option[ backbone_rebuild_test ] ) {
		rna_backbone_rebuild_test();
	} else if ( option[ convert_to_native ] ) {
		convert_to_native_test();
	} else if ( option[ pymol_struct_type ] ) {
		pymol_struct_type_test();
	} else if ( option[ env_test ] ) {
		env_sugar_test();
	} else if ( option[ sasatest ] ) {
		sasa_test();
	} else if ( option[ print_hbonds ] ) {
		print_hbonds_test();
	} else if ( option[ calc_rmsd ] ) {
		calc_rmsd_test();
	} else if ( option[ dinucleotide ] ) {
		dinucleotide_test();
	} else if ( option[ rotamerize_test ] ) {
		rotamerize_rna_test();
	} else if ( option[ build_next_nucleotide ] ) {
		build_next_nucleotide_test();
	} else if ( option[ sugar_geometry_test ] ) {
		sugar_geometry_RNA_test();
	} else if ( option[ sugar_frag_test ] ) {
		sugar_frag_RNA_test();
	} else if ( option[ color_by_geom_sol ] ) {
		color_by_geom_sol_RNA_test();
	} else if ( option[ color_by_rna_lj_base ] ) {
		color_by_lj_base_RNA_test();
	} else if ( option[ files_for_openMM ] ) {
		files_for_openMM_test();
	} else if ( option[ print_secstruct ] ) {
		print_secstruct_test();
	} else if ( option[ print_torsions ] ) {
		print_all_torsions_test();
	} else {
		//  virtual_test();
		rna_denovo_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;

		//Uh, options?
		NEW_OPT( create_vall_torsions, "Generate a torsions file from a big RNA file", false );
		NEW_OPT( create_jump_database, "Generate a database of jumps extracted from base pairings from a big RNA file", false );
		NEW_OPT( icoord_test, "Rhiju testing generation of ideal RNA coordinates", false );
		NEW_OPT( print_internal_coord, "Rhiju: generate ideal RNA coordinates", false );
		NEW_OPT( extract, "Extract RNA pdbs", false );
		NEW_OPT( fullatom_score_test, "Test RNA full atom score", false );
		NEW_OPT( fullatom_minimize, "Test RNA full minimize", false );
		NEW_OPT( fullatom_multiscore, "Test RNA fullatom scoring speed", false );
		NEW_OPT( fullatom_minimize_silent, "Test RNA full minimize on silent file", false );
		NEW_OPT( o2prime_test, "Test RNA 2'-OH rotamer trials", false );
		NEW_OPT( skip_o2prime_pack, "Turn off RNA 2'-OH packing during minimize", false );
		NEW_OPT( lores_score, "Test RNA low resolution score", false );
		NEW_OPT( lores_score_silent, "Test RNA low resolution score on silent file", false );
		NEW_OPT( env_test, "environment", false );
		NEW_OPT( rna_idealize, "RNA idealize", false );
		NEW_OPT( idl_close_chainbreaks, "RNA idealize, close chain breaks", false );
		NEW_OPT( rna_assemble, "RNA assemble", false );
		NEW_OPT( rna_design, "Test RNA design", false );
		NEW_OPT( rna_design_gap, "Test RNA design gap", false );
		NEW_OPT( rna_stats, "Check rms, no. shared base pairs, etc. for decoys", false );
		NEW_OPT( rna_torsion_check, "Test RNA torsions, reverse fold tree", false );
		NEW_OPT( sum_lores_plus_hires, "Sum in low res score after hi res minimize", false );
		NEW_OPT( fa_standard, "standard residue types",false );
		NEW_OPT( close_chainbreaks_test, "Close chainbreaks after fragment insertion in a pdb --test",false );
		NEW_OPT( rna_jumping, "Jumping test",false );
		NEW_OPT( create_benchmark, "Create RNA benchmark",false );
		NEW_OPT( minimize_rna, "Minimize RNA after fragment assembly",false );
		NEW_OPT( relax_rna, "Relax RNA after fragment assembly",false );
		NEW_OPT( simple_relax, "Relax by minimizing after any fragment insertion",false );
		NEW_OPT( color_by_geom_sol, "Color by geometric solvation",false );
		NEW_OPT( color_by_rna_lj_base, "Color by Lennard-Jones energy (just between bases)",false );
		NEW_OPT( ignore_secstruct, "Ignore sec struct in input file",false );
		NEW_OPT( chain_closure_test, "RNA CCD test",false );
		NEW_OPT( filter_base_pairs, "Filter silent file for models that satisfy structure parameters",false );
		NEW_OPT( filter_lores_base_pairs, "Filter for models that satisfy structure parameters",false );
		NEW_OPT( backbone_rebuild_test, "RNA backbone rebuild test",false );
		NEW_OPT( crazy_minimize, "RNA crazy minimize test",false );
		NEW_OPT( crazy_fold_tree, "RNA star fold tree",false );
		NEW_OPT( sasatest, "SASA calculator",false );
		NEW_OPT( fa_stack_weight, "RNA full atom stacking potential weight", 1.0 );
		NEW_OPT( cycles, "Default number of Monte Carlo cycles", 10000 );
		NEW_OPT( vall_torsions, "Torsions file?", "sampling/rna/1jj2.torsions" );
		NEW_OPT( temperature, "temperature", 0.3 );
		NEW_OPT( jump_change_frequency, "jump change frequency", 0.1 );
		NEW_OPT( close_loops, "close loops during frag insertion and jump mover", false );
		NEW_OPT( output_lores_silent_file, "output lores stuff", false );
		NEW_OPT( heat, "Heat (random frag insertions)", false );
		NEW_OPT( dump, "Dump pdb", false );
		NEW_OPT( pymol_struct_type, "struct type", false );
		NEW_OPT( convert_to_native, "Convert input pdb to rosetta-ordered native",false);
		NEW_OPT( disable_o2prime_rotamers, "In designing, don't sample 2'-OH",false);
		NEW_OPT( disable_include_current, "In designing, don't include current",false);
		NEW_OPT( sample_chi, "In designing RNA, chi torsion sample",false);
		NEW_OPT( print_hbonds, "Read in PDB output H-BONDS",false);
		NEW_OPT( calc_rmsd, "rmsd",false);
		NEW_OPT( dinucleotide, "dinucleotide",false);
		NEW_OPT( vary_geometry, "vary geometry",false);
		NEW_OPT( more_rotamers, "more rotamers",false);
		NEW_OPT( quick_test, "more rotamers",false);
		NEW_OPT( sugar_geometry_test, "sugar geometry",false);
		NEW_OPT( sugar_frag_test, "sugar frag",false);
		NEW_OPT( rmsd_cutoff, "rmsd_cutoff", 999.99 );
		NEW_OPT( rotamerize_test, "rotamerize test",false);
		NEW_OPT( build_next_nucleotide, "rotamerize test",false);
		NEW_OPT( prepend_residue, "rotamerize test",false);
		NEW_OPT( skip_coord_constraints, "skip coord constraints",false);
		NEW_OPT( atom_pair_constraint_weight, "atompair constraint weight", 0.0 );
		NEW_OPT( coordinate_constraint_weight, "coordinate constraint weight", 0.0 );
		NEW_OPT( assemble_file, "Input file for RNA assembly", "assemble.txt" );
		NEW_OPT( basepair_file, "Input file for RNA base pair definition", "default.basepairs" );
		NEW_OPT( jump_library_file, "Input file for jumps", "sampling/rna/1jj2_RNA_jump_library.dat" );
		NEW_OPT( params_file, "Input file for pairings", "default.prm" );
		NEW_OPT( data_file, "Input file for pairings", "default.prm" );
		NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );
		NEW_OPT( rsd_type_set, "Input file for RNA assembly", core::chemical::FA_STANDARD );
		NEW_OPT( files_for_openMM, "get files ready for openMM", false );
		NEW_OPT( print_secstruct, "print secondary structure", false );
		NEW_OPT( print_torsions, "print RNA torsions", false );


		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);


		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
