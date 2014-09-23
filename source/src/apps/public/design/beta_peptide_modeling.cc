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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options::OptionKeys;
using namespace basic::options;

OPT_KEY( String, algorithm )
OPT_KEY( String, force_field )
OPT_KEY( Boolean, apply_dihedral_cst )
OPT_KEY( Boolean, no_symmetry )
OPT_KEY( IntegerVector, repack_res )
OPT_KEY( Integer, n_repeat )
OPT_KEY( Integer, repeat_size )

///////////////////////////////////////////////////////////////////////////////

void
figure_out_fold_tree( pose::Pose & pose )
{
	using namespace core::conformation;

	//Look for chainbreaks in PDB.
	Size const nres = pose.total_residue();
	kinematics::FoldTree f( nres );
	Real const dist2_cutoff = 1.7 * 1.7;

	Size m( 0 );

	for (Size i=1; i < nres; ++i) {

		Residue const & current_rsd( pose.residue( i   ) ) ;
		Residue const &    next_rsd( pose.residue( i+1 ) ) ;
		Size atom_C = current_rsd.atom_index( " C  " );
		Size atom_N =    next_rsd.atom_index( " N  " );
		Real const dist2 =
			( current_rsd.atom( atom_C ).xyz() - next_rsd.atom( atom_N ).xyz() ).length_squared();

		if ( dist2 > dist2_cutoff ){
			std::cout << "Jump from " << i << " to " << i+1 << std::endl;
			f.new_jump( i, i+1, i );
			m++;
		}

	}

	pose.fold_tree( f );
}
///////////////////////////////////////////////////////////////////

void
minimize_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace pose;
	using namespace optimization;
	using namespace kinematics;
	using namespace id;

	//Setup scoring function
	std::string const & force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );

	ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	Pose pose;
	std::string pdbname = option[in::file::native]();
	pdbname.replace( pdbname.rfind(".pdb", pdbname.length() ), 4, "" );
	import_pose::pose_from_pdb ( pose, *rsd_set, option[in::file::native]() );
	Pose start_pose = pose;
	pose.dump_pdb( pdbname + "_start.pdb" );
	std::cout << "Loaded in PDB file..." << std::endl;

	////////////////////////////////
	// Personal fold tree.
/*
	Size const nres = pose.total_residue();

	int my_anchor = 20;

	int const num_jump_in( 7 );
	ObjexxFCL::FArray2D_int jump_point( 2, num_jump_in, 0);
	ObjexxFCL::FArray1D_int cuts( num_jump_in, 0);
	cuts(1) = 12;
	cuts(2) = 24;
	cuts(3) = 36;
	cuts(4) = 48;
	cuts(5) = 60;
	cuts(6) = 72;
	cuts(7) = 84;

	jump_point( 1, 1) =  5;
	jump_point( 2, 1) = 16;

	jump_point( 1, 2) = 28;
	jump_point( 2, 2) = 41;

	jump_point( 1, 3) = 53;
	jump_point( 2, 3) = 64;

	jump_point( 1, 4) = 76;
	jump_point( 2, 4) = 89;

	jump_point( 1, 5) = 20;
	jump_point( 2, 5) = 32;

	jump_point( 1, 6) = 68;
	jump_point( 2, 6) = 80;

	jump_point( 1, 7) = 32;
	jump_point( 2, 7) = 68;

	FoldTree f;
	f.tree_from_jumps_and_cuts( nres, num_jump_in, jump_point, cuts );
	f.reorder( my_anchor );
	pose.fold_tree( f );
*/


	Size const nres = pose.total_residue();
	Size const my_anchor = 1;
	figure_out_fold_tree( pose );

	////////////////////////////////////////////
	//Need to set up coordinate constraints
	////////////////////////////////////////////

	ConstraintSetOP cst_set( new ConstraintSet() );

	Real const coord_sdev( 2.0 );

	for ( Size i=1; i<= nres; ++i ) {
		Residue const & i_rsd( pose.residue(i) );
		for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
			core::scoring::func::FuncOP fx = new core::scoring::func::HarmonicFunc( 0.0, coord_sdev );
			cst_set->add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,my_anchor), i_rsd.xyz(ii), fx ) );
		}
	}


	Pose pose_start = pose;

	Real dihedral_sdev( numeric::conversions::radians( 200.0 ) );

	for ( Size i=1; i<= nres; ++i ) {
		//How about going through each atom in the tree, and finding its stub atoms?
		Residue const & i_rsd( pose.residue(i) );
		for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {

			kinematics::tree::AtomCOP current_atom ( pose.atom_tree().atom( AtomID(ii,i) ).get_self_ptr() );
			if ( current_atom->is_jump() ) continue;

			kinematics::tree::AtomCOP stub_atom1( current_atom->input_stub_atom1() );

			if ( stub_atom1->is_jump() ) continue;

			kinematics::tree::AtomCOP stub_atom2( current_atom->input_stub_atom2() );
			kinematics::tree::AtomCOP stub_atom3( current_atom->input_stub_atom3() );

			Real angle = numeric::dihedral_radians
				( current_atom->xyz(), stub_atom1->xyz(),
					stub_atom2->xyz(), stub_atom3->xyz() );

			core::scoring::func::FuncOP fx = new core::scoring::func::HarmonicFunc( angle, dihedral_sdev );
			cst_set->add_constraint( new DihedralConstraint( current_atom->id(), stub_atom1->id(),
															 stub_atom2->id(), stub_atom3->id(), fx ) );

		}

	}

	pose.constraint_set( cst_set );

	if ( option[ apply_dihedral_cst ] () ) {
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
		scorefxn->set_weight(  dihedral_constraint, 1.0 );
		scorefxn->set_weight(     angle_constraint, 1.0 );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	Energy score_orig = (*scorefxn)( pose );
	scorefxn->set_weight( coordinate_constraint, 1.0 );

	/////////////////////////////////////////////////////
	// minimzer
	AtomTreeMinimizer minimizer;
	bool const use_mm_potential = ! ( scorefxn -> get_weight( mm_twist ) == 0 );
	float const dummy_tol( 0.0000001 );
	//Default: use_nb_list= false, but MM potential fails in minimization when use_nb_list = false so use true instead
	bool const use_nblist( use_mm_potential ? true : false );
	bool const deriv_check( false );
	MinimizerOptions options( "dfpmin_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check);

	kinematics::MoveMap mm;


	mm.set_bb(  false );
	mm.set_chi( false );
	mm.set_jump( true );

	std::cout << "Minimizing... jumps" << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );

	mm.set_bb(  true );
	mm.set_chi( false );
	mm.set_jump( true );

	std::cout << "Minimizing... backbone+chi" << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );


	mm.set_bb(  true );
	mm.set_chi( true );
	mm.set_jump( true );

	std::cout << "Minimizing... jumps+backbone+chi" << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );

	scorefxn->set_weight( coordinate_constraint, 0.0 );
 	minimizer.run( pose, mm, *scorefxn, options );

	Energy minimize_score = (*scorefxn)( pose );

	std::cout << "Completed minimize_test() with new score: " << minimize_score << " vs orig: " << score_orig << std::endl;

	scorefxn->show( start_pose );
	scorefxn->show( pose );
	pose.dump_pdb( pdbname + "_minimize.pdb" );

}

///////////////////////////////////////////////////////////////////////////////
void
repack_test () {
	using namespace conformation;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace pose;
	using namespace pack;

	//Setup scoring function
	std::string const & algorithm_name = option[algorithm];
	std::string const & force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );

	ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	Pose pose;
	std::string pdbname = option[in::file::native]();
	pdbname.replace( pdbname.rfind(".pdb", pdbname.length() ), 4, "" );
	import_pose::pose_from_pdb ( pose, *rsd_set, option[in::file::native]() );
	Pose start_pose = pose;
	pose.dump_pdb( pdbname + "_start.pdb" );
	std::cout << "Loaded in PDB file..." << std::endl;

	pack::task::PackerTaskOP designtask( pack::task::TaskFactory::create_packer_task( pose ));
	designtask->initialize_from_command_line();
	if ( algorithm_name == "repack" ) {
		designtask->restrict_to_repacking();
	} else {
		utility::vector1< bool > const allowed_aas (chemical::num_canonical_aas, false);
		std::string const beta_peptide_names [21] = {"B3A", "B3C", "B3D", "B3E", "B3F", "B3G", "B3H", "B3I", "B3K", "B3L", "B3M", "B3N",
		"B3O","B3P", "B3Q", "B3R", "B3S", "B3T", "B3V", "B3W", "B3Y"};

		for(Size i = 1; i <= pose.total_residue(); ++i) {
			for (Size j = 0; j != 21; ++j) {
				designtask->nonconst_residue_task(i).allow_noncanonical_aa( beta_peptide_names[j] );
			}
			designtask->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
		}
	}

	Size const nres = pose.total_residue();
	std::cout << "NRES  = " << nres << std::endl;

	//Read in pack_residue file and setup symmetry links
	bool const is_no_symmetry = option[no_symmetry]();
	utility::vector1< core::Size > const repack_res_list = option[repack_res]();
	Size const repeats = option[n_repeat]();
	Size const repeat_unit = option[repeat_size]();

	utility::vector1< bool > residues_to_repack( nres, false );
	rotamer_set::RotamerLinksOP links( new rotamer_set::RotamerLinks() );
	links->resize( nres );
	std::cout << "residue rebuilding:" << std::endl;

	for (Size i = 1; i <= repack_res_list.size(); ++i) {
		Size const first_res = repack_res_list[i];
		utility::vector1< Size > linked_res;
		for (Size j = 0; j < repeats; ++j) {
			Size const curr_res = first_res + repeat_unit * j;
			if (curr_res > nres) {
				std::cerr << "ERROR!!! Residue " << curr_res << " outside the pose total residues!!!" << std::endl;
				exit(1);
			}
			std::cout << curr_res << ' ';
			residues_to_repack[ curr_res ] = true;
			for (Size k = 1; k <= linked_res.size(); ++k) {
				links->set_equiv(linked_res[k], curr_res);
				links->set_equiv(curr_res, linked_res[k]);
			}
			linked_res.push_back( curr_res );
		}
		std::cout << std::endl;
	}
	designtask->restrict_to_residues( residues_to_repack );
	if (repeats > 1 && (! is_no_symmetry)) designtask->rotamer_links( links );

/*

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Only repack core residues -- in Schepartz designs these are all beta3_leucine.
	// This is totally hacky, of course, but waiting for Andrew Leaver-Fay ResFile reader...
	utility::vector1< bool > residues_to_repack( nres, false );

	for ( Size ii = 1; ii <= nres ; ++ii ) {
		if ( ii > 96 ) continue;
		if ( ii % 3 == 2 ) residues_to_repack[ ii ] = true;
	}

	designtask->restrict_to_residues( residues_to_repack );

	// setup residue couplings
	rotamer_set::RotamerLinksOP links( new rotamer_set::RotamerLinks() );
	links->resize( nres );
	Size const res_list1 [4] = {2, 5, 8, 11};
	for (Size i = 0; i != 4; ++i) {
		Size const curr_res = res_list1[i];
		Size const link_res_list [4] = {curr_res, curr_res+36, curr_res+48, curr_res+84};
		for (Size j = 0; j != 4; ++j) {
			Size const res1 = link_res_list[j];
			for (Size k = 0; k != 4; ++k) {
				Size const res2 = link_res_list[k];
				if (res1 != res2) links->set_equiv(res1, res2);
			}
		}
	}

	Size const res_list2 [4] = {14, 17, 20, 23};
	for (Size i = 0; i != 4; ++i) {
		Size const curr_res = res_list2[i];
		Size const link_res_list [4] = {curr_res, curr_res+12, curr_res+48, curr_res+60};
		for (Size j = 0; j != 4; ++j) {
			Size const res1 = link_res_list[j];
			for (Size k = 0; k != 4; ++k) {
				Size const res2 = link_res_list[k];
				if (res1 != res2) links->set_equiv(res1, res2);
			}
		}
	}

	designtask->rotamer_links( links );
*/
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	Energy score_orig = (*scorefxn)( pose );

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// design
	pack::pack_rotamers( pose, *scorefxn, designtask);
	Energy design_score = (*scorefxn)( pose );
	scorefxn->weights();

	std::cout << "Completed pack_rotamers_test() with new score: " << design_score << " vs orig: " << score_orig << std::endl;
	scorefxn->show(start_pose);
	scorefxn->show(pose);

	pose.dump_pdb( pdbname + "_" + algorithm_name + ".pdb" );
}
///////////////////////////////////////////////////////////////////////////////
void
my_main ( )
{
	std::string const & algorithm_name = option[algorithm];

	if ( algorithm_name == "repack" || algorithm_name == "redesign" ) {
		repack_test();
	} else if (algorithm_name == "minimize") {
		minimize_test();
	}
}
///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	using namespace basic::options;
	utility::vector1< Size > blank_size_vector;

	NEW_OPT ( algorithm, "", "" );
	NEW_OPT ( force_field, "score_file", "" );
	NEW_OPT ( apply_dihedral_cst, "", true );
	NEW_OPT ( no_symmetry, "", false );
	NEW_OPT ( repack_res, "", blank_size_vector );
	NEW_OPT ( n_repeat, "", 1 );
	NEW_OPT ( repeat_size, "", 0 );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	my_main();

	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
