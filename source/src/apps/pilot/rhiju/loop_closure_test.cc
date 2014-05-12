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
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
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
#include <protocols/farna/util.hh>

#include <protocols/viewer/viewers.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>

//StepWise!
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/ccd_closure.hh>
#include <protocols/loops/CcdLoopClosureMover.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>


#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>
#include <core/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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
//#include <protocols/farna/RNA_FragmentsClasses.hh>
//#include <protocols/farna/RNA_DeNovoProtocol.hh>
//#include <protocols/farna/RNA_StructureParameters.hh>

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


using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;
using core::util::Error;
using core::util::Warning;
using ObjexxFCL::format::F;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
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
OPT_KEY( Integer, cut )
OPT_KEY( Integer, sample_res )
OPT_KEY( String, cst_file )
OPT_KEY( String, pack_weights )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Real, bin_width )
OPT_KEY( Boolean, auto_tune )
OPT_KEY( Boolean, sampling )

///////////////////////////////////////////////////////////////////////////////
void
extend_loop( pose::Pose & pose, utility::vector1< Size > const & loop_residues ){

	// make extended chain
	for ( Size i = 1; i <= loop_residues.size(); i++ ) {
		Size const pos = loop_residues[ i ];
		if ( ! pose.residue(pos).is_protein() ) continue;
		pose.set_phi( pos, -150 );
		pose.set_psi( pos, 150);
		pose.set_omega( pos, 180 );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
read_native_pose(  pose::PoseOP & native_pose,
									 core::chemical::ResidueTypeSetCAP & rsd_set ){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::pose;

	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		std::string native_pdb_file  = option[ in::file::native ];
		io::pdb::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );
	}
}

////////////////////////////////////////////////////////////////////
void
setup_pose_with_loop( pose::Pose & pose,
											utility::vector1< Size > const & loop_residues,
											protocols::loops::Loop & loop ){

	using namespace core::kinematics;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace protocols::loops;

	Size min_loop_res = pose.total_residue();
	Size max_loop_res = 1;
	for ( Size i = 1; i <= loop_residues.size(); i++ ) {
		if ( loop_residues[i] < min_loop_res ) min_loop_res = loop_residues[ i ];
		if ( loop_residues[i] > max_loop_res ) max_loop_res = loop_residues[ i ];
	}

	FoldTree f = pose.fold_tree();
	Size cutpoint = max_loop_res-1;
	if ( option[ cut ].user() ) cutpoint = option[ cut ]();

	f.new_jump( min_loop_res-1, max_loop_res+1, cutpoint );
	pose.fold_tree( f );
	pose.dump_pdb( "init.pdb" );

	add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint     );
	add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint + 1 );

	loop = Loop( min_loop_res, max_loop_res, cutpoint );

}


///////////////////////////////////////////////////////////////////////////////
	//void
	//copy_dofs_outside_loop(){
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



///////////////////////////////////////////////////////////////////////////////
void
ccd_loop_close( pose::Pose & pose,
								utility::vector1< Size > const & loop_residues,
								protocols::loops::Loop const& loop )
{

	using namespace core::kinematics;
	using namespace protocols::loops;

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
	fast_ccd_loop_closure( pose, *movemap, loop.start(), loop.stop(), loop.cut(),
												 10000, 0.0000001, false, 0.5, 10, 50, 75,
												 forward_deviation, backward_deviation, torsion_delta, rama_delta );
	pose.dump_pdb( "ccd_closed.pdb" );
}



///////////////////////////////////////////////////////////////////////////
void
output_chainTORS( utility::vector1< core::Real > const & dt_ang,
									utility::vector1< core::Real > const & db_ang,
									utility::vector1< core::Real > const & db_len ) {

	std::cout << "------  chainTORS output ---- " << std::endl;
	for (Size i = 1; i <= ( dt_ang.size()/3) ; i++) {

		std::cout << "TORSIONS: ";
		for (Size j = 1; j <= 3; j++) std::cout << F(8,3,dt_ang[ 3*(i-1)+ j ]) << " ";

		std::cout << "   BOND_ANGLES: ";
		for (Size j = 1; j <= 3; j++) std::cout << F(8,3,db_ang[ 3*(i-1)+ j ]) << " ";

		std::cout << "   BOND_LENGTHS: ";
		for (Size j = 1; j <= 3; j++) std::cout << F(8,3,db_len[ 3*(i-1)+ j ]) << " ";

		std::cout << std::endl;

	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
fill_chainTORS_info( pose::Pose const & pose,
					 utility::vector1<utility::vector1<Real> > & atoms,
					 utility::vector1<Real> & dt_ang,
					 utility::vector1<Real> & db_ang,
					 utility::vector1<Real> & db_len,
					 Size const & start_res_ ,
					 Size const & end_res_,
					 bool const verbose = true ) {

	using namespace numeric::kinematic_closure;

	if ( verbose ) std::cout << "About to run chainTORS" << std::endl;

	Size ind = 1;
	for (Size i =  start_res_ - 1;  i <= end_res_ + 1;		 i++) {
		if ( verbose ) std::cout << "Filling residue " << i << std::endl;
		conformation::Residue res = pose.residue(i);
		for (Size j=1; j<=3; j++) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());
			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());
			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());
			ind++;
		}
	}

	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> R0 (3);

	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);

	if ( verbose )	output_chainTORS( dt_ang, db_ang, db_len );

}

///////////////////////////////////////////////////////////////////////////////////////
void
KIC_loop_close( pose::Pose & pose,
								protocols::loops::Loop const & loop,
								utility::vector1< pose::PoseOP > & pose_list,
								bool const verbose = false ){

	using namespace core::kinematics;
	using namespace protocols::loops;
	using namespace numeric::kinematic_closure;

	if ( verbose ) std::cout << "About to kick off KIC!" << std::endl;
	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	utility::vector1<utility::vector1<Real> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<Real> dt_ang, db_len, db_ang, save_t_ang, save_b_len, save_b_ang;
	utility::vector1<Real> dummy_t_ang, dummy_b_ang, dummy_b_len;
	utility::vector1<conformation::ResidueOP> save_residues;

	Size const start_res_ = loop.start();
	Size const middle_res_ = loop.start() + 1;
	Size const end_res_ = loop.stop();
	Size const middle_offset = middle_res_ - start_res_; // is used to set central pivot atom
	Size const seg_len = end_res_ - start_res_ + 1;
	atoms.resize( (seg_len + 2) * 3); // one extra residue on each side to establish the geometric frame

	//	fill_chainTORS_info( reference_pose, atoms, save_t_ang, save_b_ang, save_b_len, start_res_, end_res_, verbose );
	fill_chainTORS_info( pose, atoms, dt_ang, db_ang, db_len, start_res_, end_res_, verbose );

	// Need to fix bond lengths and angles at cutpoint
	static const Real idl_CA_C_N_(116.2);
	static const Real idl_C_N_CA_(121.7);
	static const Real idl_C_N_(1.32869);
	static const Real OMEGA_MEAN_(179.8);
	Size const cut_offset_ = loop.cut() - start_res_;
	dt_ang[ 3 + 3*cut_offset_ + 3 ] = OMEGA_MEAN_;
	db_len[ 3 + 3*cut_offset_ + 3 ] = idl_C_N_;
	db_ang[ 3 + 3*cut_offset_ + 3 ] = idl_CA_C_N_;
	db_ang[ 3 + 3*cut_offset_ + 4 ] = idl_C_N_CA_;

	//cheat.
	//fill_chainTORS_info( reference_pose, atoms, dt_ang, db_ang, db_len, start_res_, end_res_ );

	///////////////////////////////////////////////////////////////////////////////
	// omega from CCD chain closure looks totally weirdo. Manually reset to 180.0
	///////////////////////////////////////////////////////////////////////////////
	for (Size i = 1; i <= seg_len; i++)		dt_ang[ 3*i + 3 ] = 180.0;

	if ( verbose ) 	output_chainTORS( dt_ang, db_ang, db_len );

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

	///////////////////////////////////
	// Perform loop closure
	///////////////////////////////////
	// outputs from loop closure
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	if ( verbose ) std::cout << "About to run bridgeObjects" << std::endl;
	bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
	if ( verbose ) std::cout << "Finished bridgeObjects" << std::endl;

	//	std::cout << "Found number of solutions " << t_ang.size() << std::endl;
	for (Size i = 1; i <= t_ang.size(); i++) {
		for ( core::Size res = 0; res < seg_len; res++ ){

			Real const phi   = t_ang[ i ][ (3*(res+1)) + 1 ];
			Real const psi   = t_ang[ i ][ (3*(res+1)) + 2 ];
			Real const omega = t_ang[ i ][ (3*(res+1)) + 3 ];

			Size const resnum = start_res_ + res;
			//			std::cout << resnum << " " << phi << " " << psi << " " << omega << std::endl;
			pose.set_phi  ( resnum,  phi);
			pose.set_psi  ( resnum,  psi);
			pose.set_omega( resnum,  omega);
		}
		if ( verbose ) {
			pose.dump_pdb( "KIC_"+ ObjexxFCL::string_of( i )+".pdb" );
		}

		// Save it.
		pose::PoseOP pose_op = new pose::Pose;
		*pose_op = pose;
		pose_list.push_back( pose_op );

		//		protocols::viewer::clear_conformation_viewers();		exit( 0 );

	}

	// just for output.
	if ( verbose ) fill_chainTORS_info( pose, atoms, dt_ang, db_ang, db_len, start_res_, end_res_ );

}

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

	//////////////////////
	// POSE SETUP
	//////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	Pose pose, pose_input;
	io::pdb::pose_from_pdb( pose_input, *rsd_set, option[ in::file::s]()[1] );
	pose = pose_input;

	// Read in native
	PoseOP native_pose;
	read_native_pose( native_pose, rsd_set );

	Loop loop;
	utility::vector1< Size > const loop_residues = option[ loop_res ]();
	setup_pose_with_loop( pose, loop_residues, loop );
	extend_loop( pose, loop_residues );
	pose.dump_pdb( "extended.pdb" );

	//copy_dofs_outside_loop();
	Pose start_pose = pose;

	// CCD closure
	ccd_loop_close( pose, loop_residues, loop );

	// KIC closure
	utility::vector1< pose::PoseOP > pose_list;
	KIC_loop_close( pose, loop, pose_list, true /*verbose*/ );

}



void
output_to_silent( Size const count,
									pose::Pose const & pose,
									Size const nsol,
									Size const & sample_residue,
									utility::vector1< Size > const & loop_residues,
									pose::PoseOP & native_pose,
									core::io::silent::SilentFileData & silent_file_data,
									std::string const & silent_file_out ){

	using namespace core::io::silent;
	using namespace core::scoring;

	std::string const tag = "S_" + lead_zero_string_of( count, 5);
	BinarySilentStruct s( pose, tag );
	s.add_energy( "nsol", nsol );

	Size k( 1 );
	s.add_energy( "tau"+string_of( k++ ),  pose.phi( sample_residue ) );
	s.add_energy( "tau"+string_of( k++ ),  pose.psi( sample_residue ) );

	for ( Size m = 1; m <= loop_residues.size(); m++ ) {
		s.add_energy( "tau"+string_of( k++ ),  pose.phi( loop_residues[m] ) );
		s.add_energy( "tau"+string_of( k++ ),  pose.psi( loop_residues[m] ) );
	}

	if ( native_pose ){
		s.add_energy( "backbone_rms", rmsd_no_super( pose, *native_pose, is_protein_backbone_including_O ) );
	}

	//silent_file_data.write_silent_struct( s, silent_file_out, true /*just scores*/ );
	silent_file_data.write_silent_struct( s, silent_file_out, true /*just scores*/ );
	silent_file_data.add_structure( s );

}


///////////////////////////////////////////////////////////////////////////////
void
sampling_closure_test(){

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

	//////////////////////
	// POSE SETUP
	//////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	Pose pose, pose_input;
	io::pdb::pose_from_pdb( pose_input, *rsd_set, option[ in::file::s]()[1] );
	pose = pose_input;

	// Read in native
	PoseOP native_pose;
	read_native_pose( native_pose, rsd_set );

	Loop loop;
	utility::vector1< Size > const loop_residues = option[ loop_res ]();
	setup_pose_with_loop( pose, loop_residues, loop );
	extend_loop( pose, loop_residues );
	pose.dump_pdb( "extended.pdb" );

	Pose start_pose = pose;

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	SilentFileDataOP silent_file_data = new SilentFileData;
	std::string const silent_file_out = option[ out::file::silent]();
	std::string const score_all_file_out = "SCORE_ALL_" + silent_file_out;

	Real const bin_size( option[ bin_width]()  );
	Size const num_bins = static_cast<Size>(360.0/bin_size);
	Size const sample_residue = option[ sample_res ]();
	ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );

	Size count( 0 );
	(*scorefxn)( pose );
	output_to_silent( count, pose, 0, sample_residue, loop_residues, native_pose, *silent_file_data, score_all_file_out  );

	for ( Size i = 1; i <= num_bins; i++ ){

		std::cout << "Scanning " << i << " out of " << num_bins << std::endl;
		Real const phi  = i * ( 360.0 / num_bins );
		pose.set_phi( sample_residue, phi );

		for ( Size j = 1; j <= num_bins; j++ ){

			Real const psi  = j * ( 360.0 / num_bins );
			pose.set_psi( sample_residue, psi );

			/////////////////////////////////////////
			// KIC closure
			/////////////////////////////////////////
			utility::vector1< pose::PoseOP > pose_list;
			KIC_loop_close( pose, loop, pose_list, false /*verbose*/ );
			/////////////////////////////////////////

			Size const nsol = pose_list.size();
			for ( Size n = 1; n <= nsol; n++ ){

					Real const score = (*scorefxn)( *pose_list[n] );
					count++;
					output_to_silent( count, *pose_list[n], nsol, sample_residue, loop_residues, native_pose, *silent_file_data, score_all_file_out  );
				}

		}
	}

	StepWiseLegacyClusterer stepwise_clusterer( silent_file_data );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );
	utility::vector1< Size > calc_rms_res;
	calc_rms_res.push_back( sample_residue );
	for ( Size n = 1; n <= loop_residues.size(); n++ ) calc_rms_res.push_back( loop_residues[n] );
	stepwise_clusterer.set_calc_rms_res( calc_rms_res );
	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius	);
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	// Do it!
	stepwise_clusterer.cluster();
	stepwise_clusterer.output_silent_file( silent_file_out );


}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ sampling ] ){
		sampling_closure_test();
	} else {
		loop_closure_test();
	}


	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

	return 0;
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
	NEW_OPT( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT( cut, "cutpoint", 0 );
	NEW_OPT( sample_res, "cutpoint", 0 );
	NEW_OPT( sampling, "sampling over 1 residue, close the others...", false );
	NEW_OPT( loop_res, "Loop residues to remodel", blank_size_vector );
	NEW_OPT( bin_width, "width of bins", 10.0 );

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
		return -1;
	}

}
