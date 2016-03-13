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

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/Methods.hh>

//#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/OutputMovers.hh>
//#include <protocols/moves/PackMover.hh>
//#include <protocols/rigid/RigidBodyMover.hh>
//#include <protocols/moves/rigid_body_moves.hh>
//#include <protocols/moves/TrialMover.hh>


#include <protocols/viewer/viewers.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>

#include <core/util/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/keys/OptionKeys.hh>

#include <core/options/after_opts.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/options/util.hh>//option.hh>
//#include <core/options/after_opts.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//REMOVE LATER!
//#include <utility/io/izstream.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;


using namespace core;
using namespace protocols;

using utility::vector1;



///////////////////////////////////////////////////////////////////////////////
Real
rhiju_pack( pose::Pose & pose, scoring::ScoreFunctionOP & scorefxn )
{

	using namespace scoring;
	using namespace  pose;

	int const nres = pose.total_residue();

	utility::vector1< bool > residues_to_repack( nres, true );

	pack::task::PackerTaskOP packertask( pack::task::TaskFactory::create_packer_task( pose ));
	packertask->initialize_from_command_line();//.or_include_current( true );
	packertask->restrict_to_residues( residues_to_repack );
	packertask->restrict_to_repacking();

	Energy score_orig = (*scorefxn)( pose );
	EnergyMap emap_orig = pose.energies().total_energies();
	emap_orig *= scorefxn->weights();

	clock_t starttime = clock();
	pack::pack_rotamers( pose, *scorefxn, packertask);
	clock_t stoptime = clock();
	std::cout << "pack_rotamers with design took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
	Energy design_score = (*scorefxn)( pose );
	EnergyMap emap_final = pose.energies().total_energies();
	emap_final *= scorefxn->weights();


	std::cout << "Completed pack_rotamers_test() with new score: " << design_score << " vs orig: " << score_orig << std::endl;
	//	std::cout << "emap_orig" << std::endl;
	//	emap_orig.print();
	std::cout << "emap_final:" << std::endl;
	emap_final.print();
	//std::cout << "diff" << std::endl;
	EnergyMap ediff( emap_final );
	ediff -= emap_orig;
	//	ediff.print();

	return design_score;
}

///////////////////////////////////////////////////////////////////////////////
Real
soft_pack( pose::Pose & pose )
{
	using namespace scoring;
	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( SOFT_REP_WTS, SCORE12_PATCH ) );
	return rhiju_pack( pose, scorefxn );
}

///////////////////////////////////////////////////////////////////////////////
Real
hard_pack( pose::Pose & pose )
{
	using namespace scoring;
	ScoreFunctionOP scorefxn( get_score_function() );
	return rhiju_pack( pose, scorefxn );
}

///////////////////////////////////////////////////////////////////////////////
void
hard_minimize( pose::Pose & pose ){

	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace pose;

	using namespace optimization;
	using namespace kinematics;
	using namespace id;

	Size const nres = pose.total_residue();

	// Use constraints to disallow big movements.
	ConstraintSetOP cst_set( new ConstraintSet() );

	Real const coord_sdev( 2.0 );

	int const my_anchor = nres/2;

	for ( Size i=1; i<= nres; ++i ) {
		Residue const & i_rsd( pose.residue(i) );
		for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
			cst_set->add_constraint( new CoordinateConstraint( AtomID(ii,i),
																												 AtomID(1,my_anchor),
																												 i_rsd.xyz(ii),
																												 new HarmonicFunc( 0.0, coord_sdev ) ) );
		}
	}

// Turn on hard fa_rep
// Minimize chi angles
// Minimize backbone angles
// Minimize rigid body

	ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );

	Energy score_orig = (*scorefxn)( pose );
	EnergyMap emap_orig = pose.energies().total_energies();
	emap_orig *= scorefxn->weights();

	scorefxn->set_weight( coordinate_constraint, 1.0 );

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.000001 );
	bool const use_nblist( true ), deriv_check( false );
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check);

	kinematics::MoveMap mm;

	mm.set_bb(  false );
	mm.set_chi( true );
	mm.set_jump( false );

	std::cout << "Minimizing... chi" << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );

	mm.set_bb(  true );
	mm.set_chi( true );
	mm.set_jump( false );

	std::cout << "Minimizing... backbone+chi" << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );


	mm.set_bb(  true );
	mm.set_chi( true );
	//	mm.set_jump( true );

	std::cout << "Minimizing... jumps+backbone+chi" << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );

	scorefxn->set_weight( coordinate_constraint, 0.0 );
 	minimizer.run( pose, mm, *scorefxn, options );


}

///////////////////////////////////////////////////////////////////////////////
void
get_subpose( pose::Pose & pose, pose::Pose & subpose,
						 Size const start1, Size const end1,
						 Size const start2, Size const end2)
{

	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	using namespace kinematics;

	utility::vector1< Size > positions;

	subpose.clear();

	std::cout << "About to create a subpose " << std::endl;

	for (Size i=start1; i <= end1; ++i ) {
		conformation::Residue const & rsd( pose.residue( i ) );
		subpose.append_residue_by_jump( rsd, 1 );
	}

	for (Size i=start2; i <= end2; ++i ) {
		conformation::Residue const & rsd( pose.residue( i ) );
		subpose.append_residue_by_jump( rsd, 1 );
	}

	int const nres = subpose.total_residue();
	subpose.conformation().insert_chain_ending( nres-1 );

	std::cout << "SUBPOSE TOTAL RES " << nres << std::endl;
	FoldTree f( nres );
	std::cout << "SUBPOSE TOTAL RES " << nres << std::endl;
	subpose.fold_tree( f );
	//	create_subpose( pose, positions, f, subpose );

	std::cout << "Finished creating a subpose " << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
void
find_closest_points( pose::Pose & pose,
										 Size const start1, Size const end1,
										 Size const start2, Size const end2,
										 Size & closest_i, Size & closest_j )
{
	using namespace chemical;
	using namespace conformation;
	using namespace pose;

	Real closest_dist2 = 9999.9;
	for ( Size i=start1; i<= end1; ++i ) {
		Residue const & i_rsd( pose.residue(i) );

		for ( Size j=start2; j<= end2; ++j ) {
			Residue const & j_rsd( pose.residue(j) );

			Real dist2 = ( i_rsd.xyz( i_rsd.atom_index("CA") )
										 - j_rsd.xyz( j_rsd.atom_index("CA") )  ).length_squared();

			if (dist2 < closest_dist2 ){
				closest_i = i;
				closest_j = j;
				closest_dist2 = dist2;
			}
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
get_useful_fold_trees_for_gcc185( pose::Pose & pose,
											kinematics::FoldTree & f_minimize,
											kinematics::FoldTree & f_repack)
// Set up fold tree
{
	int const num_jump_in( 3 );
	FArray2D_int jump_point( 2, num_jump_in, 0);
	FArray1D_int cuts( num_jump_in, 0);

	// MAGIC NUMBERS. In principle this could be much more flexible,
	// looking at chainbreaks within input PDB to figure out what's going on.

	Size gtpase1_start,	gtpase1_end  ,	coiledcoil1_start,	coiledcoil1_end  ,	gtpase2_start,	gtpase2_end  ,	coiledcoil2_start,	coiledcoil2_end  ;

	if ( false )
	{ // appropriate for mini_3bbp
	gtpase1_start = 1;
	gtpase1_end   = 6;
	coiledcoil1_start = 7;
	coiledcoil1_end   = 13;

	gtpase2_start = 14;
	gtpase2_end   = 19;
	coiledcoil2_start = 20;
	coiledcoil2_end   = 26;
	} else
	{ // appropriate for 3bbp_threecoiledcoilturns
	gtpase1_start = 1;
	gtpase1_end   = 161;
	coiledcoil1_start = 162;
	coiledcoil1_end   = 189;

	gtpase2_start = 190;
	gtpase2_end   = 350;
	coiledcoil2_start = 351;
	coiledcoil2_end   = 378;
	}

	Size const nres = pose.total_residue();

	//////////////////////////////////////////////////////////////
	// This connects GTPase1 <--> coil1 <--> coil2 <--> GTPase2
	//////////////////////////////////////////////////////////////
	Size i,j;

	find_closest_points( pose,
											 gtpase1_start, gtpase1_end,
											 coiledcoil1_start, coiledcoil1_end,
											 i, j );
	jump_point( 1, 1) = i;
	jump_point( 2, 1) = j;
	cuts( 1 ) = gtpase1_end;

	find_closest_points( pose,
											 coiledcoil1_start, coiledcoil1_end,
											 coiledcoil2_start, coiledcoil2_end,
											 i, j );
	jump_point( 1, 2) = i;
	jump_point( 2, 2) = j;
	cuts( 2 ) = coiledcoil1_end;

	find_closest_points( pose,
											 gtpase2_start, gtpase2_end,
											 coiledcoil2_start, coiledcoil2_end,
											 i, j );
	jump_point( 1, 3) = i;
	jump_point( 2, 3) = j;
	cuts( 3 ) = gtpase2_end;

	f_minimize.tree_from_jumps_and_cuts( nres, num_jump_in, jump_point, cuts );


	//////////////////////////////////////////////////////////////
	// This connects GTPase2 <--> GTPase1 <--> coil1 <--> coil2
	//////////////////////////////////////////////////////////////
	find_closest_points( pose,
											 gtpase1_start, gtpase1_end,
											 coiledcoil1_start, coiledcoil1_end,
											 i, j );
	cuts( 1 ) = gtpase1_end;

	find_closest_points( pose,
											 coiledcoil1_start, coiledcoil1_end,
											 coiledcoil2_start, coiledcoil2_end,
											 i, j );
	jump_point( 1, 2) = i;
	jump_point( 2, 2) = j;
	cuts( 2 ) = gtpase2_end;

	find_closest_points( pose,
											 gtpase1_start, gtpase1_end,
											 gtpase2_start, gtpase2_end,
											 i, j );
	jump_point( 1, 3) = i;
	jump_point( 2, 3) = j;
	cuts( 3 ) = coiledcoil1_end;

	f_repack.tree_from_jumps_and_cuts( nres, num_jump_in, jump_point, cuts );

	return;
}

///////////////////////////////////////////////////////////////////////////////
void
get_ddG( pose::Pose & pose, std::string const pdb ) {

	using namespace pose;

	Pose pose1, pose2;

	partition_pose_by_jump( pose, 1, pose1, pose2 );

	Real const Etot     = hard_pack( pose );
	Real const E1       = hard_pack( pose1 );
	Real const E2       = hard_pack( pose2 );

	dump_pdb( pose1, pdb+"_partner1.pdb" );
	dump_pdb( pose2, pdb+"_partner2.pdb" );

	std::cout << "DELDELG " << pdb << " " << Etot << " " <<
		E1 << " " << E2 << " " << (Etot - E1 - E2) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
void
repack_minimize_test( ){

	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	using namespace optimization;
	using namespace kinematics;
	using namespace id;

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;

	//Read in PDB
	Pose pose;
	//std::string const pdb = "mini_3bbp";

	std::string pdb  = option[ in ::file::s ][1];

	//	std::string const pdb = "3bbp_threecoiledcoilturns";

	io::pdb::pose_from_file( pose, pdb+".pdb" , core::import_pose::PDB_file);

	int const nres = pose.total_residue();
	std::cout << "NRES: " << nres << std::endl;
	dump_pdb( pose, pdb+"_start.pdb" );

	FoldTree f_minimize, f_repack;
	//THIS HAS MAGIC NUMBERS!!!
	get_useful_fold_trees_for_gcc185( pose, f_minimize, f_repack );

	////////////////////////////////////////////
	// Soft repack
	////////////////////////////////////////////
	pose.fold_tree( f_repack );
	soft_pack( pose );
	dump_pdb( pose, pdb+"_repack.pdb" );

	////////////////////////////////////////////
	// Hard minimize
	////////////////////////////////////////////
	pose.fold_tree( f_minimize );
	hard_minimize( pose );
	dump_pdb( pose, pdb+"_minimize.pdb" );


	///////////////////////////////////////////////////////////////////////////////////////////
	// Calculate ddG, by parsing out coiled-coil, GTPase, combined system ... do hard repacks.
	///////////////////////////////////////////////////////////////////////////////////////////
	pose.fold_tree( f_repack );
	get_ddG( pose, pdb );

	// Save to outfile?


}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;
	using namespace core::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	devel::init(argc, argv);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	repack_minimize_test();
	exit(0);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
