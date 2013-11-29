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
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/scoring/dunbrack/DunbrackRotamer.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>

#include <protocols/viewer/viewers.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//StepWise!
#include <protocols/swa/StepWiseClusterer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/ccd_closure.hh>
#include <protocols/loops/CcdLoopClosureMover.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>

//clustering
#include <protocols/cluster/cluster.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/optimizeH.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/pose_stream/ExtendedPoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>
#include <core/util/datacache/BasicDataCache.hh>
#include <core/util/datacache/CacheableString.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>


#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
//RNA stuff.
//#include <protocols/rna/RNA_FragmentsClasses.hh>
//#include <protocols/rna/RNA_DeNovoProtocol.hh>
//#include <protocols/rna/RNA_StructureParameters.hh>

//Job dsitributor
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <vector>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/cluster.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, cluster_test )
OPT_KEY( StringVector, s1 )
OPT_KEY( StringVector, s2 )
OPT_KEY( StringVector, silent1 )
OPT_KEY( StringVector, silent2 )
OPT_KEY( StringVector, tags1 )
OPT_KEY( StringVector, tags2 )
OPT_KEY( IntegerVector, slice_res1 )
OPT_KEY( IntegerVector, slice_res2 )
OPT_KEY( IntegerVector, input_res1 )
OPT_KEY( IntegerVector, input_res2 )
OPT_KEY( IntegerVector, loop_res )
OPT_KEY( String, cst_file )
OPT_KEY( String, pack_weights )
OPT_KEY( Boolean, cluster_by_all_atom_rmsd )
OPT_KEY( Boolean, output_start )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Boolean, auto_tune )
OPT_KEY( Boolean, parse_pathway )


///////////////////////////////////////////////////////////////////////////////
void
loop_closure_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::io::pose_stream;
	using namespace core::pose;
	using namespace core::pack;
	using namespace protocols::swa;
	using namespace protocols::loops;
	using namespace numeric::kinematic_closure;

	//////////////////////
	// POSE SETUP
	//////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	Pose pose, pose_input;
	io::pdb::pose_from_pdb( pose_input, *rsd_set, option[ in::file::s]()[1] );

	pose = pose_input;
	//	make_pose_from_sequence( pose, pose_input.sequence(), *rsd_set );

	/////////////////////
	// Read in native
	PoseOP native_pose;
	bool native_exists( false );
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		std::string native_pdb_file  = option[ in::file::native ];
		io::pdb::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );
	}

	utility::vector1< Size > const loop_residues = option[ loop_res ]();

	Size min_loop_res = pose.total_residue();
	Size max_loop_res = 1;
	for ( Size i = 1; i <= loop_residues.size(); i++ ) {
		if ( loop_residues[i] < min_loop_res ) min_loop_res = loop_residues[ i ];
		if ( loop_residues[i] > max_loop_res ) max_loop_res = loop_residues[ i ];
	}

	FoldTree f = pose.fold_tree();
	Size const cutpoint = max_loop_res-1;
	f.new_jump( min_loop_res-1, max_loop_res+1, cutpoint );
	pose.fold_tree( f );
	pose.dump_pdb( "init.pdb" );

	if ( true ) {
		// make extended chain
		for ( Size i = 1; i <= loop_residues.size(); i++ ) {
			Size const pos = loop_residues[ i ];
			if ( ! pose.residue(pos).is_protein() ) continue;
			pose.set_phi( pos, -150 );
			pose.set_psi( pos, 150);
			pose.set_omega( pos, 180 );
		}
	}


	pose.dump_pdb( "extended.pdb" );

// 	std::map< Size, Size > res_map;
// 	FArray1D<bool> is_loop( pose.total_residue(), false );
// 	for ( Size i = 1; i <= loop_residues.size(); i++ ) is_loop( loop_residues[ i ] ) = true;
// 	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
// 		if (!is_loop( n ) ) {
// 			res_map[ n ] = n;
// 			std::cout << "will copy dofs " << n << std::endl;
// 		}
// 	}
// 	copy_dofs( pose, pose_input, res_map );
// 	pose.dump_pdb( "copy_dofs.pdb" );

	Pose start_pose = pose;

	add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint     );
	add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint + 1 );


	Loop loop( min_loop_res, max_loop_res, cutpoint );


	///////////////////////////////////////
	// CCD closure
	///////////////////////////////////////
	if ( true  ) {
		MoveMapOP movemap( new MoveMap);
		movemap->set_bb( false );
		movemap->set_chi( false );
		movemap->set_jump( false );
		for ( Size i = 1; i <= loop_residues.size(); i++ ) {
			movemap->set_bb( loop_residues[ i ], true );
		}
		movemap->set_bb( true );

		///// CCD close
		//CcdLoopClosureMover ccd_closer( loop, movemap );
		//	ccd_closer.set_ccd_cycles( 1000 );
		//	ccd_closer.apply( pose );
		Real forward_deviation, backward_deviation, torsion_delta, rama_delta;
		fast_ccd_loop_closure( pose, *movemap, loop.start(), loop.stop(), loop.cut(), 10000, 0.0000001, false, 0.5, 10, 50, 75,
													 forward_deviation, backward_deviation, torsion_delta, rama_delta );
		pose.dump_pdb( "ccd_closed.pdb" );
	}

	std::cout << "About to kick off KIC!" << std::endl;
	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	utility::vector1<utility::vector1<Real> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	// outputs from loop closure
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, save_t_ang, save_b_len, save_b_ang, R0 (3);
	utility::vector1<Real> dummy_t_ang, dummy_b_ang, dummy_b_len;
	utility::vector1<conformation::ResidueOP> save_residues;
 	Size ind;

	Size const start_res_ = loop.start();
	Size const middle_res_ = loop.start() + 1;
	Size const end_res_ = loop.stop();
	Size const middle_offset = middle_res_ - start_res_; // is used to set central pivot atom
	Size const seg_len = end_res_ - start_res_ + 1;
	atoms.resize((seg_len + 2) * 3); // one extra residue on each side to establish the geometric frame

// 	std::cout << "About to fill atoms" << std::endl;
// 	ind = 1;
// 	for (Size i =  start_res_ - 1;  i <= end_res_ + 1;		 i++) {
// 		std::cout << "Filling residue " << i << std::endl;
// 		conformation::Residue res = start_pose.residue(i);
// 		for (Size j=1; j<=3; j++) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
// 			atoms[ind].resize(3);
// 			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());
// 			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());
// 			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());
// 			ind++;
// 		}
// 	}

// 	std::cout << "About to run chainTORS" << std::endl;

// 	chainTORS(atoms.size(), atoms, save_t_ang , save_b_ang, save_b_len, R0, Q0);

	std::cout << "About to run chainTORS" << std::endl;
	ind = 1;
	for (Size i =  start_res_ - 1;  i <= end_res_ + 1;		 i++) {
		std::cout << "Filling residue " << i << std::endl;
		conformation::Residue res = pose.residue(i);
		for (Size j=1; j<=3; j++) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());
			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());
			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());
			ind++;
		}
	}

	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);

	order[1]=1;
	order[2]=2;
	order[3]=3;
	// Set the pivot atoms
	Size pvatom1=5; // second C-alpha
	Size pvatom2=5 + (3 * middle_offset); // middle res C-alpha
	Size pvatom3=(3 * (seg_len+1)) - 1; // second-to-last C-alpha
	pivots[1]=pvatom1;
	pivots[2]=pvatom2;
	pivots[3]=pvatom3;

// 	//idealize?
// 	if ( false ){
// 		// set all bond lengths, angles, and omegas for closure to ideal values
// 		Real idl_C_N_CA_(121.7), 	idl_N_CA_C_(111.2), 	idl_CA_C_N_(116.2), 	idl_C_N_(1.32869), 	idl_N_CA_(1.458), 	idl_CA_C_(1.52326);
// 		for (Size i=1; i<=atoms.size(); i+=3) {
// 			db_ang[i]=idl_C_N_CA_;
// 			db_ang[i+1]=idl_N_CA_C_;
// 			db_ang[i+2]=idl_CA_C_N_;
// 			db_len[i]=idl_N_CA_;
// 			db_len[i+1]=idl_CA_C_;
// 			db_len[i+2]=idl_C_N_;
// 			//dt_ang[i+2]=OMEGA_MEAN_;
// 		}
// 	}


	std::cout << "About to run bridgeObjects" << std::endl;

	///////////////////////////////////
	// Perform loop closure
	///////////////////////////////////
	bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);

	std::cout << "Finished bridgeObjects" << std::endl;

	for (Size i = 1; i <= t_ang.size(); i++) {
		for( core::Size res = 0; res < seg_len; res++ ){
			pose.set_phi( start_res_ + res, t_ang[ i ][ (3*(res+1)) + 1 ] );
			pose.set_psi( start_res_ + res, t_ang[ i ][ (3*(res+1)) + 2 ] );
		}
		pose.dump_pdb( "KIC_"+string_of( i )+".pdb" );
	}


}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	loop_closure_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	//Uh, options?
	NEW_OPT( cluster_test, "cluster", false );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );

	NEW_OPT( s1, "input file(s)", blank_string_vector );
	NEW_OPT( s2, "input file(s)", blank_string_vector );
	NEW_OPT( silent1, "input file", blank_string_vector );
	NEW_OPT( silent2, "input file", blank_string_vector );
	NEW_OPT( tags1, "input tag(s)", blank_string_vector );
	NEW_OPT( tags2, "input tag(s)", blank_string_vector );
	NEW_OPT( slice_res1, "Residues to slice out of starting file", blank_size_vector );
	NEW_OPT( slice_res2, "Residues to slice out of starting file", blank_size_vector );
	NEW_OPT( input_res1, "Residues already present in starting file", blank_size_vector );
	NEW_OPT( input_res2, "Residues already present in starting file2", blank_size_vector );
	NEW_OPT( pack_weights, "weights for green packing", "standard.wts" );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT( cluster_by_all_atom_rmsd, "cluster by all atom rmsd", false );
	NEW_OPT( output_start, "output starting pdb", false );
	NEW_OPT( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT( parse_pathway, "parse the pathway", false );

	NEW_OPT( loop_res, "Loop residues to remodel", blank_size_vector );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

	exit( 0 );

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
