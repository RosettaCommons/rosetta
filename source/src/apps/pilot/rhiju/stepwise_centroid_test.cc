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
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/util.hh>
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
#include <protocols/swa/StepWiseFilterer.hh>
#include <protocols/swa/StepWiseClusterer.hh>
#include <protocols/swa/StepWisePoseMinimizer.hh>
#include <protocols/swa/StepWisePoseSetup.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/swa/StepWiseResidueSampler.hh>

//clustering
#include <protocols/cluster/cluster.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>

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
#include <core/io/pose_stream/PoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
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

#include <ObjexxFCL/format/formatted.o.hh>
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
#include <core/options/keys/abinitio.OptionKeys.gen.hh>
#include <core/options/keys/frags.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;
typedef std::map< std::string, core::pose::PoseOP > PoseList;

OPT_KEY( Boolean, start_from_scratch )
OPT_KEY( Boolean, n_terminus )
OPT_KEY( Boolean, c_terminus )
OPT_KEY( Boolean, cluster_by_all_atom_rmsd )
OPT_KEY( Boolean, cluster_test )
OPT_KEY( Boolean, centroid_screen )
OPT_KEY( Integer, max_input )
OPT_KEY( Integer, min_res )
OPT_KEY( Integer, max_res )
OPT_KEY( Integer, insert_res )
OPT_KEY( String, cst_file )
OPT_KEY( String, centroid_weights )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( FileVector, s2 )


///////////////////////////////////////////////////////////////////////////////
void
output_centroid_silent_struct(
															core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
															core::io::silent::SilentFileDataOP & sfd, std::string const & tag,
															std::string const & parent_tag,
															std::string const silent_file = "" ){

	using namespace core::io::silent;
	using namespace core::scoring;

	// Will this give the intuitive phi,psi,omega --basic non-fullatom silent struct?
	ProteinSilentStruct s( pose, tag );

	if ( native_pose_op != 0 ){

		core::pose::Pose pose_for_superpose = pose;
		core::pose::Pose const & native_pose_for_superpose( *native_pose_op );

		Real const rmsd = CA_rmsd( pose_for_superpose, native_pose_for_superpose );
		Real const backbone_rmsd = rmsd_with_super(  pose_for_superpose, native_pose_for_superpose, is_protein_backbone_including_O);
		Real const all_rmsd = rms_at_corresponding_heavy_atoms( pose_for_superpose, native_pose_for_superpose );

		s.add_energy( "rms", rmsd );
		s.add_energy( "all_rms", all_rmsd );
		s.add_energy( "backbone_rms", backbone_rmsd );
	}


	s.add_comment( "PARENT_TAG", parent_tag );

	sfd->add_structure( s );
	if (silent_file.size() > 0 )	sfd->write_silent_struct( s, silent_file, false /*write score only*/ );

}


////////////////////////////////////////////////////////////////////////////////////////////////////
Size
figure_out_nested_positions(
														std::string const & inside_sequence,
														std::string const & desired_sequence,
														Size const min_start_res = 0)
{

	Size const max_start_res = desired_sequence.size() - inside_sequence.size()+1;

	for (Size potential_start_res = 1; potential_start_res <= max_start_res; potential_start_res++ ) {

		if ( potential_start_res < min_start_res ) continue;

		bool found_match = true;
		//Look for exact sequence match.
		for ( Size i = 0; i < inside_sequence.size(); i++ ) {
			if ( inside_sequence[i] != desired_sequence[i - 1 + potential_start_res] ) {
				found_match = false;
				break;
			}
		}
		if ( found_match ) {
			return potential_start_res;
		}
	}

	return 0;

}

///////////////////////////////////////////////////////////////////////////////
void
slice_pdb( core::pose::Pose & pose,
					 Size const & start_res,
					 Size const & end_res  )
{
	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::conformation;

	Pose new_pose;
	new_pose.clear();

	for ( Size i = start_res; i <= end_res; i++ )  {
		ResidueOP residue_to_add = pose.residue( i ).clone() ;
		new_pose.append_residue_by_bond(  *residue_to_add  ) ;
	}

	pose = new_pose;

}

///////////////////////////////////////////////////////////////////////////////
void
slice_sequence( std::string & sequence,
					 Size const & start_res,
					 Size const & end_res  )
{
	std::string const sequence_new = sequence.substr( start_res-1,   end_res - start_res + 1 );
	sequence = sequence_new;
}

///////////////////////////////////////////////////////////////////////////////
// might be useful in a util.cc somewhere
core::scoring::constraints::ConstraintSetOP
constraint_set_slice( core::scoring::constraints::ConstraintSetOP & cst_set, Size const & min_res, Size const & max_res )
{

	using namespace core::scoring::constraints;
	using namespace core::scoring;
	using namespace core::id;

	ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );

	ConstraintCOPs csts( cst_set->get_all_constraints() );

	for ( Size n = 1; n <= csts.size(); n++ ) {

		ConstraintCOP const & cst( csts[n] );

		if ( cst->score_type() == atom_pair_constraint)  { // currently only defined for pairwise distance constraints.
			Size const i = cst->atom( 1 ).rsd();
			Size const j = cst->atom( 2 ).rsd();
			//			Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
			//			if ( dist  > separation_cutoff ) continue;
			if (i < min_res || i > max_res ) continue;
			if (j < min_res || j > max_res ) continue;

			AtomID atom1_new( cst->atom(1).atomno(),  cst->atom(1).rsd() - min_res + 1);
			AtomID atom2_new( cst->atom(2).atomno(),  cst->atom(2).rsd() - min_res + 1);
			ConstraintOP cst_new = new AtomPairConstraint( atom1_new, atom2_new,
																										 cst->get_func().clone() /*is this defined?*/, cst->score_type() );

			cst_set_new->add_constraint( cst_new );

		}

	}


	std::cout << "NUM CONSTRAINTS " << cst_set_new->get_all_constraints().size() << " out of " <<
		csts.size() << std::endl;

	return cst_set_new;
}

///////////////////////////////////////////////////////////////////////////////
void
apply_all_frags(
	 pose::Pose & pose,
	 pose::PoseOP & native_pose,
	 Size const & insert_pos,
	 std::string const & parent_tag,
	 Real const & score_cut,
	 Real & score_min,
	 Size & count,
	 core::io::silent::SilentFileDataOP & sfd,
	 core::fragment::ConstantLengthFragSetOP & fragset,
	 core::scoring::ScoreFunctionOP & scorefxn )
{

	using namespace core::fragment;

	std::cout << "STARTING FROM POSE " << parent_tag << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	// Find all fragments for first residue. Maybe this should be done outside the loop.
	FrameList frames;
	fragset->frames( insert_pos, frames );

	// apply fragments one at a time, score, and save in silent file.
	std::cout << "NUMBER OF FRAMES: " << frames.size() << std::endl;

	FrameOP & frame( frames[1] );

	std::cout << "NUMBER OF FRAGSS " << frame->nr_frags() << std::endl;

	for ( Size q = 1; q <= frame->nr_frags(); q++ ) {

		frame->apply( q, pose );

		Real const score = (*scorefxn)( pose );

		count += 1;
		std::string const tag(  "S_" + string_of( count - 1) );

		if ( count == 1 ) score_min = score;
		if ( score < score_min ) score_min = score;

		if ( score < score_min + score_cut ) {
			output_centroid_silent_struct( pose, native_pose, sfd, tag, parent_tag  /*, "junk.out"*/ );
		}
	}
	//			pose.dump_pdb( tag+".pdb");

}

////////////////////////////////////////////////////////////////////////
Size
setup_pose( pose::Pose & pose,
						pose::Pose const & start_pose,
						std::string const & desired_sequence,
						Size const min_start_res = 0 )
{
	Size const start_res = figure_out_nested_positions( start_pose.sequence(), desired_sequence, min_start_res );
	Size const end_res = start_res +  start_pose.sequence().size() - 1;

	//Now actually copy into the pose.
	Size count( 0 );
	for ( Size n = start_res; n <= end_res; n++ ) {
		count++;
		pose.set_phi( n, start_pose.phi( count ) );
		pose.set_psi( n, start_pose.psi( count ) );
		pose.set_omega( n, start_pose.omega( count ) );
		pose.set_secstruct( n, start_pose.secstruct( count ) );
	}

	return end_res;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
rebuild_centroid_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::io::pose_stream;
	using namespace core::fragment;

	////////////////////////////////////////////////////
	Size const min_res_( option[ min_res ]() );
	Size const max_res_( option[ max_res ]() );
	Size insert_pos( option[ insert_res ]() );

	////////////////////////////////////////////////////
	//Read in sequence information and native
	////////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( CENTROID );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( option[ centroid_weights ]()  );

	//Read in desired fasta. [for now copy from native]
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string desired_sequence = fasta_sequence->sequence();
	std::string const desired_sequence_full( desired_sequence );
	if ( min_res_ > 0 ) 	slice_sequence( desired_sequence, min_res_, max_res_ );

	PoseOP native_pose;
	bool native_exists( false );
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		std::string native_pdb_file  = option[ in::file::native ];
		io::pdb::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );
		native_exists = true;
		if ( min_res_ > 0 ) slice_pdb( *native_pose, min_res_, max_res_ );
		//		native_pose->dump_pdb( "NATIVE.pdb" );
	}

	////////////////////////////////////////////////////
	// Actual read in of any starting poses.
	////////////////////////////////////////////////////

	// Can have up to two input streams.
	PoseInputStreamOP input1, input2;
	bool start_from_scratch( false );

	if ( option[ in::file::s].user() ) {
		// pdb input(s).

		input1 = new PDBPoseInputStream( option[ in::file::s ]() );

		if (option[ s2 ].user() ) {
			input2 = new PDBPoseInputStream( option[ s2 ]() );
		}

	} else if ( option[ in::file::silent].user() ) {
		// silent input(s).

		utility::vector1< std::string > const & silent_files_in = option[ in::file::silent ]();

		if ( silent_files_in.size() > 1 ) {
			assert( silent_files_in.size() == 2);
			utility::vector1< std::string > silent_files_in1, silent_files_in2;

			silent_files_in1.push_back( silent_files_in[ 1 ] );
			silent_files_in2.push_back( silent_files_in[ 2 ] );

			input1 = new SilentFilePoseInputStream( silent_files_in1 );
			input2 = new SilentFilePoseInputStream( silent_files_in2 );

		} else {
			if ( option[ in::file::tags].user() ){
				input1 = new SilentFilePoseInputStream( silent_files_in ,
																							 option[ in::file::tags ]() );
			} else {
				input1 = new SilentFilePoseInputStream( option[ in::file::silent ]()[ 1 ] );
			}
		}
	} else {
		start_from_scratch = true;
	}

	// Current allowing:
	// Either 0 inputs ( start from scratch ), 1 input (typically fragment insertions at termini), or
	//   2 inputs (internal fragment insertions between patches poses)
	Size const num_inputs ( static_cast<Size>(input1 != 0) +
													static_cast<Size>(input2 != 0)  );

	///////////////////////////////////
	// pose initialize.
	Pose pose;
	make_pose_from_sequence( pose, desired_sequence, *rsd_set, false /*auto_termini*/);

	// make extended chain
	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
		if ( ! pose.residue(pos).is_protein() ) continue;
		pose.set_phi( pos, -150 );
		pose.set_psi( pos, 150);
		pose.set_omega( pos, 180 );
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	/////////////////////////////
	// Constraints...
	ConstraintSetOP cst_set;
	if ( option[ cst_file ].user() ) {
		Pose full_pose; /* Need to make this full_pose for constraint machinery atom ID checking */
		make_pose_from_sequence( full_pose, desired_sequence_full, *rsd_set, false /*auto_termini*/);
		cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, full_pose );
		if ( min_res_ > 0 ) cst_set = constraint_set_slice( cst_set, min_res_, max_res_ );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
	}
	pose.constraint_set( cst_set );

	///////////////////////////////////
	// Read in fragment library.
	std::string const frag_file  = option[ in::file::frag_files ]()[ 1 ];
	core::fragment::ConstantLengthFragSetOP fragset;
	fragset = new ConstantLengthFragSet( 0 /*frag_length ... is reset by reader*/, frag_file );
	if ( min_res_ > 0 ) {
		fragment_set_slice( fragset, min_res_, max_res_ );
		insert_pos = insert_pos - min_res_ + 1;
	}

	Size count( 0 );
	Real const score_cut = option[ score_diff_cut ]();
	Real score_min( 0.0 );
	SilentFileDataOP sfd = new SilentFileData;
	std::string parent_tag("");
	Size const max_input_poses( option[ max_input ] );

	///////////////////////////////////
	// main loop.
	///////////////////////////////////
	if ( num_inputs == 0 ){
		///////////////////////
		// start from scratch
		assert(  pose.total_residue() == fragset->max_frag_length() );

		insert_pos = 1;
		parent_tag = "START_FROM_SCRATCH";
		apply_all_frags( pose, native_pose,  insert_pos, parent_tag, score_cut, score_min, count, sfd, fragset, scorefxn );

	} else if  (num_inputs == 1 ){

		///////////////////////
		// input pose.

		Pose start_pose;

		while( input1->has_another_pose() ) {
			input1->fill_pose( start_pose, *rsd_set );
			setup_pose( pose, start_pose, desired_sequence );
			parent_tag = tag_from_pose( start_pose );

			apply_all_frags( pose, native_pose,  insert_pos, parent_tag, score_cut, score_min, count, sfd, fragset, scorefxn );
		}

	} else {

		assert( num_inputs == 2);
		pose::Pose start_pose1, start_pose2;

		Size num_pose1( 0 );
		while( input1->has_another_pose() ) {

			num_pose1++;
			if ( max_input_poses > 0 && num_pose1 > max_input_poses ) break;

			input1->fill_pose( start_pose1, *rsd_set );
			Size const end_res1 = setup_pose( pose, start_pose1, desired_sequence );
			std::string const parent_tag1 = tag_from_pose( start_pose1 );

			Size num_pose2( 0 );
			while( input2->has_another_pose() ) {

				num_pose2++;
				if ( max_input_poses > 0 && num_pose2 > max_input_poses ) break;

				input2->fill_pose( start_pose2, *rsd_set );
				setup_pose( pose, start_pose2, desired_sequence, end_res1 /*start_pose2 better be after the first one.*/ );
				std::string const parent_tag2 = tag_from_pose( start_pose2 );

				parent_tag = parent_tag1+"_"+parent_tag2;

				apply_all_frags( pose, native_pose,  insert_pos, parent_tag, score_cut, score_min, count, sfd, fragset, scorefxn );

			}

			input2->reset();

		}
	}

	std::cout << "NUM DECOYS ======> " << sfd->size() << std::endl;


	protocols::swa::StepWiseClusterer stepwise_clusterer( sfd );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );
	stepwise_clusterer.set_cluster_radius(	option[ OptionKeys::cluster::radius ]()	);
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] ); // false by default
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
	stepwise_clusterer.cluster();

	stepwise_clusterer.output_silent_file( option[ out::file::silent] );

}

///////////////////////////////////////////////////////////////
void
cluster_outfile_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	protocols::swa::StepWiseClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 2500 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );

	stepwise_clusterer.set_cluster_radius(	option[ OptionKeys::cluster::radius ]()	);
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] );
	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	stepwise_clusterer.cluster();

	std::string const silent_file_out( option[ out::file::silent  ]() );
	stepwise_clusterer.output_silent_file( silent_file_out );

}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ cluster_test ] ){
		cluster_outfile_test();
	} else {
		rebuild_centroid_test();
	}

	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	//Uh, options?
	NEW_OPT( start_from_scratch, "start from scratch", false );
	NEW_OPT( n_terminus, "build N terminus", false );
	NEW_OPT( c_terminus, "build C terminus", false );
	NEW_OPT( cluster_test, "cluster", false );
	NEW_OPT( cluster_by_all_atom_rmsd, "cluster by all atom rmsd", false );
	NEW_OPT( centroid_screen, "centroid screen", false );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );
	NEW_OPT( min_res, "min. residue for fragment library slicing", 0 );
	NEW_OPT( max_res, "max. residue for fragment library slicing", 0 );
	NEW_OPT( max_input, "max. number of input poses from, e.g., silent files", 0 );
	NEW_OPT( insert_res, "where to insert fragment", 0 );
	NEW_OPT( centroid_weights, "centroid score function", "score3.wts" );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 10.0 );

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
