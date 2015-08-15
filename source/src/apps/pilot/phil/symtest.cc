// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/pilot/phil/test1.cc
/// @brief  Some simple examples of how to use basic functionality + some DNA functionality
/// @author Phil Bradley (pbradley@fhcrc.org)

// libRosetta headers
#include <apps/pilot/phil/phil.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

// #include <protocols/frags/VallData.hh>
// #include <protocols/frags/TorsionFragment.hh>

// #include <core/scoring/dna/setup.hh>
// #include <core/scoring/dna/base_geometry.hh>
// #include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/Energies.hh>
// #include <core/scoring/EnergyGraph.hh>
// #include <core/scoring/GenBornPotential.hh>
// #include <core/scoring/LREnergyContainer.hh>
// #include <core/scoring/methods/Methods.hh>
// #include <utility/excn/Exceptions.hh>

// #include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
// #include <protocols/moves/MonteCarlo.hh>
// #include <protocols/moves/Mover.hh>
// #include <protocols/simple_moves/PackRotamersMover.hh>
// #include <protocols/moves/MoverContainer.hh>
// #include <protocols/moves/OutputMovers.hh>
// #include <protocols/simple_moves/RotamerTrialsMover.hh>
// #include <protocols/rigid/RigidBodyMover.hh>
// //#include <protocols/moves/rigid_body_moves.hh>
// #include <protocols/moves/TrialMover.hh>
// #include <protocols/moves/RepeatMover.hh>

// #include <protocols/loops/loop_closure/ccd/ccd_closure.hh>
// #include <protocols/loops/loops_main.hh>

// #include <protocols/viewer/viewers.hh>

// #include <core/types.hh>

// #include <core/scoring/sasa.hh>

// #include <basic/prof.hh> // profiling
// #include <core/pose/datacache/CacheableDataType.hh> // profiling
// #include <basic/datacache/BasicDataCache.hh>

// #include <core/sequence/DerivedSequenceMapping.hh>
// #include <core/sequence/util.hh>


// #include <core/chemical/AA.hh>
// #include <core/conformation/Residue.hh>
// #include <core/conformation/ResidueMatcher.hh>
// #include <core/pack/rotamer_set/RotamerCouplings.hh>
// #include <core/pack/rtmin.hh>
// #include <core/pack/min_pack.hh>
// #include <core/chemical/ResidueTypeSet.hh>
// #include <core/chemical/ResidueTypeSelector.hh>
// #include <core/conformation/ResidueFactory.hh>
// #include <core/chemical/VariantType.hh>

// #include <core/chemical/ChemicalManager.hh>

// #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// #include <core/pack/dunbrack/RotamerLibrary.hh>
// #include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
// #include <core/scoring/hbonds/HBondSet.hh>
// #include <core/scoring/hbonds/hbonds.hh>

// #include <core/pack/rotamer_trials.hh>
// #include <core/pack/pack_rotamers.hh>
// #include <core/pack/task/PackerTask.hh>
// #include <core/pack/task/TaskFactory.hh>
// #include <core/pack/task/ResfileReader.hh>

// #include <core/kinematics/FoldTree.hh>
// #include <protocols/viewer/visualize.hh>
// #include <core/kinematics/MoveMap.hh>
// #include <core/id/AtomID_Map.hh>


// #include <core/optimization/AtomTreeMinimizer.hh>
// #include <core/optimization/MinimizerOptions.hh>

// #include <core/pose/Pose.hh>
// #include <core/pose/util.hh>

// #include <basic/options/util.hh>//option.hh>
// //#include <basic/options/after_opts.hh>

// #include <basic/basic.hh>


// #include <core/init/init.hh>
#include <devel/init.hh>

// #include <core/io/pdb/pose_io.hh>

// #include <utility/vector1.hh>

// #include <numeric/xyzVector.hh>
// #include <numeric/constants.hh>

// #include <ObjexxFCL/format.hh>
// #include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// //#include <utility/io/izstream.hh>


// // C++ headers
// //#include <cstdlib>
// #include <fstream>
// #include <iostream>
// #include <string>

#include <basic/Tracer.hh>

// // option key includes

// #include <basic/options/keys/out.OptionKeys.gen.hh>
// #include <basic/options/keys/phil.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>

// #include <core/import_pose/import_pose.hh>
// #include <utility/vector0.hh>
// #include <numeric/xyz.functions.hh>

// //Auto Headers
// #include <core/pose/util.tmpl.hh>


////////////////////////////////////////////////
// danger USING ////////////////////////////////
// using namespace core;
// using namespace protocols;
// using namespace id;
// using namespace pose;
// using namespace conformation;
// using namespace chemical;
// using namespace scoring;
// using namespace basic::options;
// using namespace optimization;
// using utility::vector1;
// using std::string;
// using namespace ObjexxFCL;
// using namespace ObjexxFCL::format;
// using core::import_pose::pose_from_pdb;
// using io::pdb::dump_pdb; // deprecated though
static thread_local basic::Tracer TR( "apps.pilot.phil.symtest" );


///////////////////////////////////////////////////////////////////////////////
// rhino01 source$ python src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type cn -nsub 2
// symmetry_name c2
// subunits 2
// recenter
// number_of_interfaces  1
// E = 2*VRT0001 + 1*(VRT0001:VRT0002)
// anchor_residue COM
// virtual_transforms_start
// start -1,0,0 0,1,0 0,0,0
// rot Rz 2
// virtual_transforms_stop
// connect_virtual JUMP1 VRT0001 VRT0002
// set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)


void
setup_symmetric_dimer(
											Vector const & axis,
											Vector const & center,
											Size const monomer_root_residue,
											Pose & pose // a dimer, structurally symmetric, but not a "symmetric pose" yet
											)
{
	Size const nres_protein( pose.total_residue() ), nres_monomer( nres_protein/2 );
	runtime_assert( nres_protein == nres_monomer * 2 );
	runtime_assert( !pose::symmetry::is_symmetric( pose ) ); // not yet
	runtime_assert( fabs( 1.0 - axis.length() ) < 1e-3 ); // should be a unit vector

	// add two virtual residues along the symm axis
	ResidueOP vrtrsd
		( conformation::ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ) );
	vrtrsd->set_xyz("ORIG",center);

	// orient first vrt with x-axis toward monomer 1, or maybe away from
	Vector centroid(0,0,0);
	for ( Size i=1; i<= nres_monomer; ++i ) centroid += pose.residue(i).nbr_atom_xyz();
	centroid /= nres_monomer;

	Vector const z( axis ), y( ( z.cross( centroid - center ) ).normalized() ), x( y.cross(z) );
	vrtrsd->set_xyz( "X", center - x ); // minus here to match the symmdef file ??
	vrtrsd->set_xyz( "Y", center + y );

	pose.append_residue_by_jump( *vrtrsd, 1 );
	pose.conformation().insert_chain_ending( nres_protein );

	// now rotate vrtrsd 180 degrees about the "z" axis
	vrtrsd->set_xyz( "X", center + x );
	vrtrsd->set_xyz( "Y", center - y );
	pose.append_residue_by_jump( *vrtrsd, 1 );

	Size const vrt1( nres_protein+1 ), vrt2( nres_protein+2 );

	// set up a better foldtree
	{
		FoldTree f( pose.total_residue() );
		// jump from vrt1 to monomer1
		f.new_jump( monomer_root_residue, vrt1, nres_monomer );
		// jump from vrt2 to monomer2
		f.new_jump( monomer_root_residue+nres_monomer, vrt2, nres_protein );
		// jump from vrt1 to vrt2
		f.new_jump( vrt1, vrt2, vrt1 );
		f.reorder( vrt1 );
		pose.fold_tree( f );
	}


	// now symmetry setup
	conformation::symmetry::SymmetryInfo symminfo;

	for ( Size i=1; i<= nres_monomer; ++i ) {
		symminfo.add_bb_clone( i, i+nres_monomer );
		symminfo.add_chi_clone( i, i+nres_monomer );
	}

	symminfo.add_jump_clone( 1, 2, 0.0 );

	/// jumps that build the monomers have all dofs flexible except z-trans and z-rot
	///
	{
		using core::conformation::symmetry::SymDof;
		map< Size, SymDof > symdofs;
		SymDof symdof;
		/// this part I'm not so sure about...
		/// it may be useful to have more dofs here, for minimization, even though theoretically there are only 4 dofs
		symdof.read( "x angle_x angle_y angle_z" );
		//symdof.read( "x y angle_x angle_y angle_z" ); // better flexibility?
		//symdof.read( "x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)" ); // set ranges of motion?
		symdofs[ 1 ] = symdof; // setting for jump#1, which goes from vrt1 to monomer1
		symminfo.set_dofs( symdofs );
	}


	symminfo.num_virtuals( 2 );
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( pose.total_residue(), 1 );

	for ( Size i=1; i<= nres_protein; ++i ) {
		if ( i < nres_monomer ) symminfo.set_score_multiply( i, 2 );
		else symminfo.set_score_multiply( i, 1 );
	}

	symminfo.update_score_multiply_factor();

	pose::symmetry::make_symmetric_pose( pose, symminfo );

	runtime_assert( pose::symmetry::is_symmetric( pose ) ); // should be now
}

void
bk_test()
{
	// load monomer
	Pose monomer;
	pose_from_pdb( monomer, start_file() );

	monomer.dump_pdb("monomer.pdb");

	// define axis, center
	Vector axis( 0,0,1 ), center(0,0,0);

	Pose dimer( monomer );

	// transform monomer
	Matrix R( rotation_matrix( axis, numeric::constants::d::pi ) );
	Vector const v(0,0,0); // Rx+v fixes center in this special case so we are OK

	monomer.apply_transform_Rx_plus_v( R, v );

	for ( Size i=1; i<= monomer.total_residue(); ++i ) {
		if ( i==1 ) dimer.append_residue_by_jump( monomer.residue(i), 1 );
		else dimer.append_residue_by_bond( monomer.residue(i) );
	}

	setup_symmetric_dimer( axis, center, 1, dimer ); // would normally choose a smarter position than 1

	dimer.dump_pdb("dimer.pdb");

	{ // minimize
		ScoreFunctionOP scorefxn( get_score_function() );
		runtime_assert( pose::symmetry::is_symmetric( *scorefxn ) ); // requires stoopid command line flag
		scorefxn->set_weight( fa_rep, 0.01 );

		Real const start_score( (*scorefxn)( dimer ) );

		TR.Trace << "start_score: " << F(9,3,start_score) << ' ' <<
			dimer.energies().total_energies().weighted_string_of( scorefxn->weights() ) << endl;

		kinematics::MoveMapOP movemap = new kinematics::MoveMap;
		movemap->set_jump(1,true);
		protocols::simple_moves::symmetry::SymMinMoverOP min_mover
			( new protocols::simple_moves::symmetry::SymMinMover(movemap, scorefxn, "dfpmin", 0.00001, true ) );

		min_mover->apply( dimer );
		dimer.dump_pdb("dimer_aftermin.pdb");
	}

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
		using namespace basic::options;

		devel::init(argc, argv);

		bk_test();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
