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
#include <protocols/farna/RNA_ProtocolUtil.hh>

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

//StepWiseProtein!
#include <protocols/stepwise/StepWiseClusterer.hh>
#include <protocols/stepwise/protein/StepWiseProteinFilterer.hh>
#include <protocols/stepwise/protein/StepWiseProteinPoseMinimizer.hh>
#include <protocols/stepwise/protein/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/protein/StepWiseProteinScreener.hh>
#include <protocols/stepwise/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/protein/StepWiseProteinResidueSampler.hh>
#include <protocols/stepwise/protein/MainChainTorsionSet.hh>

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
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/pose_stream/PoseInputStream.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>
#include <core/util/datacache/BasicDataCache.hh>
#include <core/util/datacache/CacheableString.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>
////
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>


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
OPT_KEY( Boolean, repack )
OPT_KEY( Boolean, fullatom )
OPT_KEY( Boolean, minimize )
OPT_KEY( Boolean, rebuild )
OPT_KEY( Boolean, start_from_scratch )
OPT_KEY( Boolean, make_ideal_helix )
OPT_KEY( Boolean, cluster_test )
OPT_KEY( Boolean, cluster_by_all_atom_rmsd )
OPT_KEY( Boolean, sample_trp )
OPT_KEY( Boolean, sample_trp_tyr )
OPT_KEY( Boolean, score12_plot )
OPT_KEY( Boolean, filter_native_big_bins )
OPT_KEY( Boolean, rename_tags )
OPT_KEY( Boolean, peptide_plane )
OPT_KEY( Boolean, n_terminus )
OPT_KEY( Boolean, c_terminus )
OPT_KEY( Boolean, add_peptide_plane )
OPT_KEY( Boolean, minimize_test )
OPT_KEY( Boolean, deriv_check )
OPT_KEY( Boolean, centroid_screen )
OPT_KEY( Boolean, no_sample_junction )
OPT_KEY( Boolean, entropy_test )
OPT_KEY( Boolean, centroid_output )
OPT_KEY( Boolean, ghost_loops )
OPT_KEY( Boolean, trans_omega )
OPT_KEY( Boolean, color_by_lj )
OPT_KEY( Real, filter_rmsd )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Real, centroid_score_diff_cut )
OPT_KEY( Integer, sample_residue )
OPT_KEY( Integer, start_res )
OPT_KEY( Integer, end_res )
OPT_KEY( Integer, n_sample )
OPT_KEY( Integer, nstruct_centroid )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( String, pack_weights )
OPT_KEY( String, cst_file )
OPT_KEY( String, centroid_weights )

///////////////////////////////////////////////////////////////////////////////
void
pdbslice( core::pose::Pose & pose,
					Size const & start_res,
					Size const & end_res  )
{
	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::conformation;

	Pose new_pose;
	new_pose.clear();

	for ( Size i = start_res; i <= end_res; i++ )  {
		//std::cout << "About to append " << i << std::endl;
		ResidueOP residue_to_add = pose.residue( i ).clone() ;
		//if ( i == start_res ) residue_to_add.
		new_pose.append_residue_by_bond(  *residue_to_add  ) ;
	}

	pose = new_pose;

}

///////////////////////////////////////////////////////////////////////////////
// now deprecated!!!??
void
sample_rama_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::io::silent;
	using namespace ObjexxFCL::format;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace protocols::stepwise::protein;

	// Read in protein
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	//rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( CENTROID );

	pose::PoseOP pose_op( new pose::Pose );
	pose::Pose & pose( *pose_op );

	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, *rsd_set, pdb_file );

	Size n1 = option[ sample_residue ]();
	if ( option[ start_res ].user() ){
		assert( option[ end_res ].user() );
		pdbslice( pose, option[ start_res ](), option[ end_res ]()  );
		pose.dump_pdb( "slice.pdb" ) ;
		n1 = n1 - option[ start_res ]() + 1;
	}

	utility::vector1< Size > moving_residues;
	moving_residues.push_back( n1 );
	moving_residues.push_back( n1+1 );
	//	moving_residues.push_back( n1+2 );

	pose::PoseOP native_pose = pose_op;
	pose::Pose fa_pose = pose;
	switch_to_residue_type_set( pose, CENTROID );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "cen_std.wts"  );
	scorefxn->set_weight( rama, 1.0 );
	scorefxn->set_weight( hbond_lr_bb, 1.0 );
	scorefxn->set_weight( hbond_sr_bb, 1.0 );

	// Quick checkaroo on phi,psi,rama.
	(*scorefxn)( pose );
	std::cout << "RAMA CHECK!" << std::endl;
	Ramachandran ramachandran;
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		std::cout << I(3,i) << F(8,2,pose.phi(i)) << " " << F(8,2,pose.psi(i)) << " " <<
			F(8,2,pose.omega(i)) << " " <<
			F(8,3,pose.energies().onebody_energies( i )[ rama ]) << " " <<
			ramachandran.eval_rama_score_residue( pose.aa( i ), pose.phi( i ), pose.psi( i ) ) <<
			std::endl;
	}


	//	Size count( 0 );
	//	Real rmsd( 0.0 );
	PoseList pose_list;
	std::string const silent_file = option[ out::file::silent  ]();
	std::string const tag( "NATIVE" );
	//std::string const dummy_file( "dummy.out" );
	//	quick_output( pose, fa_pose, tag, moving_residues, dummy_file, native_pose, pose_list );

	// Can we recreate the pose from scratch?
	if ( option[ make_ideal_helix ]() ) {
		make_pose_from_sequence( pose, pose.sequence(), *rsd_set );
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			pose.set_phi( i, -70.0 );
			pose.set_psi( i, -40.0 );
			pose.set_omega( i, 180.0 );
		}
	}

	pose.dump_pdb( "start.pdb" );

	///////////////
	///////////////
	// DEPRECATED!!
	///////////////
	///////////////
// 	StepWiseProteinScreener stepwise_screener( moving_residues );
// 	stepwise_screener.set_rmsd_cutoff( option[ filter_rmsd ]() );
// 	stepwise_screener.set_n_sample( option[ n_sample ]() );
// 	stepwise_screener.set_native_pose( native_pose );

// 	stepwise_screener.apply( pose );
// 	utility::vector1< MainChainTorsionSetList >const & main_chain_torsion_set_lists
// 		= stepwise_screener.main_chain_torsion_set_lists();

// 	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	StepWiseProteinResidueSampler stepwise_residue_sampler( moving_residues, main_chain_torsion_set_lists );
// 	ScoreFunctionOP pack_scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );
// 	stepwise_residue_sampler.set_native_pose( native_pose );
// 	stepwise_residue_sampler.set_scorefxn( pack_scorefxn );
// 	stepwise_residue_sampler.set_silent_file( silent_file );
// 	stepwise_residue_sampler.apply( pose );

	std::cout << "--- DONE! --" << std::endl;

}



///////////////////////////////////////////////////////////////////////////////
void
minimizer_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::io::pose_stream;
	using namespace core::pose;
	using namespace protocols::stepwise::protein;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );


	PoseOP native_pose_op;
	if ( option[ in::file::native ].user() ) {
		std::string native_pdb_file  = option[ in::file::native ];
		native_pose_op = new Pose;
		io::pdb::pose_from_pdb( *native_pose_op, *rsd_set, native_pdb_file );
	}


	bool const deriv_check_( option[ deriv_check] );
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );
	kinematics::MoveMap mm;
	mm.set_bb( true );
	mm.set_chi( true );

	pose::Pose pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	ScoreFunctionOP scorefxn = getScoreFunction();

	if ( option[ in::file::s].user() ) {

		std::string pdb_file  = option[ in::file::s ][1];
		io::pdb::pose_from_pdb( pose, *rsd_set, pdb_file );

		if ( option[ trans_omega ]() ) for ( Size n = 1; n <= pose.total_residue(); n++ ) pose.set_omega( n, 180.0 );

		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		exit( 0 );

		minimizer.run( pose, mm, *scorefxn, options );

		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		std::string const out_file =  "minimize.pdb";
		dump_pdb( pose, out_file );


	} else {

		utility::vector1< std::string > const silent_files( option[ in::file::silent ]() );
		PoseInputStreamOP input  = new SilentFilePoseInputStream( silent_files );

		std::string const silent_file_out = "minimize_" + silent_files[1];
		SilentFileData silent_file_data;

		while ( input->has_another_pose() ) {

			input->fill_pose( pose, *rsd_set );

			(*scorefxn)( pose );
			scorefxn->show( std::cout, pose );

			exit( 0 );

			minimizer.run( pose, mm, *scorefxn, options );

			//			SilentStructOP s = new BinaryProteinSilentStruct( pose, tag_from_pose( pose ) );
			//			silent_file_data.write_silent_struct( *s, silent_file_out, false /*write score only*/ );

			output_silent_struct( pose, native_pose_op, silent_file_out, tag_from_pose( pose ) );

		}

	}


}


///////////////////////////////////////////////////////////////////////
// repack -- this is silliness.
void
repack_test(){

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::optimization;
	using namespace core::pose;
	using namespace core::options;
	using namespace core::options::OptionKeys;


	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	//rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( CENTROID );
	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, *rsd_set, pdb_file );

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->restrict_to_repacking();

	static ScoreFunctionOP fa_scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );

	fa_scorefxn->show( std::cout, pose );
	pack::pack_rotamers( pose, *fa_scorefxn, task);
	fa_scorefxn->show( std::cout, pose );

	pose.dump_pdb( "repack.pdb" );

}


///////////////////////////////////////////////////////////////////////
Real
get_sidechain_rmsd( pose::Pose const & pose, pose::Pose const & start_pose, Size const & n )
{

	chemical::ResidueType rsd_type( pose.residue_type( n ) );

	Real d2( 0.0 );
	Size count( 0 );
	for ( Size j = rsd_type.first_sidechain_atom(); j <= rsd_type.nheavyatoms(); j++ ) {
		d2 += ( pose.residue( n ).xyz( j ) - start_pose.residue( n ) .xyz( j ) ).length_squared();
		count += 1;
	}

	return std::sqrt( d2 / count );

}

///////////////////////////////////////////////////////////////////////
void
sample_trp_test()
{

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::optimization;
	using namespace core::pose;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace protocols::moves;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;


	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose, *rsd_set, pdb_file );

	Size n = option[ sample_residue ](); //tryptophan.

	protocols::simple_moves::UserDefinedGroupDiscriminatorOP user_defined_group_discriminator( new UserDefinedGroupDiscriminator);
	utility::vector1< Size > group_ids;
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( i == n ) {
			group_ids.push_back( 0 );
		} else {
			group_ids.push_back( 1 );
		}
	}
	user_defined_group_discriminator->set_group_ids( group_ids );
	protocols::simple_moves::GreenPackerOP green_packer( new GreenPacker );
	green_packer->set_group_discriminator( user_defined_group_discriminator );

	static ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );
	green_packer->set_scorefunction( *scorefxn );

	TaskFactoryOP initial_task_factory( new TaskFactory );
	initial_task_factory->push_back( new InitializeFromCommandline );
	initial_task_factory->push_back( new RestrictToRepacking );
	PreventRepackingOP prevent_repacking( new PreventRepacking );
	prevent_repacking->include_residue( n );
	initial_task_factory->push_back( prevent_repacking );
	green_packer->set_reference_round_task_factory( initial_task_factory );

	TaskFactoryOP general_task_factory( new TaskFactory );
	general_task_factory->push_back( new InitializeFromCommandline );
	general_task_factory->push_back( new RestrictToRepacking );
	general_task_factory->push_back( prevent_repacking );
	green_packer->set_task_factory( general_task_factory );

	std::string const silent_file = option[ out::file::silent  ]();

	PoseOP native_pose_op = new Pose;
	std::string native_pdb_file  = option[ in::file::native ];
	io::pdb::pose_from_pdb( *native_pose_op, *rsd_set, native_pdb_file );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Size const N_SAMPLE( option[ n_sample]() );
	Size count( 0 );
	Pose start_pose = pose;

	for (Size i1 = 1; i1 <= N_SAMPLE; i1++ ) {
		for (Size j1 = 1; j1 <= N_SAMPLE; j1++ ) {

			pose.set_chi( 1, n, protocols::stepwise::protein::get_rotamer_angle( i1, N_SAMPLE ) );
			pose.set_chi( 2, n, protocols::stepwise::protein::get_rotamer_angle( j1, N_SAMPLE ) );
			green_packer->apply( pose );

			( *scorefxn )( pose );
			setPoseExtraScores( pose, "trp_rms", get_sidechain_rmsd( pose, start_pose, n )  );
			std::string const tag = "S_"+ lead_zero_string_of( ++count, 5 );
			protocols::stepwise::protein::output_silent_struct( pose, native_pose_op, silent_file, tag );

		}
	}

}


///////////////////////////////////////////////////////////////////////
void
add_chi_tags( core::pose::Pose & pose, core::pose::Pose const & native_pose, Size const n_trp, Size const n_tyr ){

	if ( pose.residue_type( n_trp ).nchi() > 1 ) {
		setPoseExtraScores( pose, "trp_chi1", pose.chi( 1, n_trp ) );
		setPoseExtraScores( pose, "trp_chi2", pose.chi( 2, n_trp ) );
	} else {
		setPoseExtraScores( pose, "trp_chi1", native_pose.chi( 1, n_trp) );
		setPoseExtraScores( pose, "trp_chi2", native_pose.chi( 2, n_trp) );
	}

	if ( pose.residue_type( n_tyr ).nchi() > 1 ) {
		setPoseExtraScores( pose, "tyr_chi1", pose.chi( 1, n_tyr ) );
		setPoseExtraScores( pose, "tyr_chi2", pose.chi( 2, n_tyr ) );
	} else {
		setPoseExtraScores( pose, "tyr_chi1", native_pose.chi( 1, n_tyr) );
		setPoseExtraScores( pose, "tyr_chi2", native_pose.chi( 2, n_tyr) );
	}

}

///////////////////////////////////////////////////////////////////////
void
sample_trp_tyr_test()
{

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::optimization;
	using namespace core::pose;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace protocols::moves;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose, pose_trp_tyr, pose_trp, pose_tyr, pose_null, pose_input;
	std::string pdb_file  = option[ in::file::s ][1];
	io::pdb::pose_from_pdb( pose_input, *rsd_set, pdb_file );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//std::string const silent_file = option[ out::file::silent ]();
	std::string const silent_file_trp_tyr = "sample_trp_tyr.out";
	std::string const silent_file_trp     = "sample_trp.out";
	std::string const silent_file_tyr     = "sample_tyr.out";
	std::string const silent_file_null     = "sample_null.out";

	// substitute everything except n1, n2 to ala/gly.
	std::string const input_sequence = pose_input.sequence();
	std::string new_sequence_trp_tyr, new_sequence_trp, new_sequence_tyr, new_sequence_null;
	Size n_trp( 0 ), n_tyr(  0 );
	std::map< Size, Size > res_map;

	for (Size n = 1; n <= pose_input.total_residue(); n++ ){
		res_map[ n ] = n;
		if ( input_sequence[n-1] == 'W' ) {
			n_trp = n;
		} else if ( input_sequence[n-1] == 'Y' ) {
			n_tyr = n;
		} else if ( input_sequence[n-1] != 'G' ) {
			new_sequence_trp_tyr += 'A'; continue;
		}
		new_sequence_trp_tyr += input_sequence[ n-1 ];
	}

	for (Size n = 1; n <= pose_input.total_residue(); n++ ){
		if ( n == n_tyr )
			new_sequence_trp += 'A';
		else
			new_sequence_trp += new_sequence_trp_tyr[ n-1 ];
	}

	for (Size n = 1; n <= pose_input.total_residue(); n++ ){
		if ( n == n_trp )
			new_sequence_tyr += 'A';
		else
			new_sequence_tyr += new_sequence_trp_tyr[ n-1 ];
	}

	for (Size n = 1; n <= pose_input.total_residue(); n++ ){
		if ( n == n_trp || n == n_tyr )
			new_sequence_null += 'A';
		else
			new_sequence_null += new_sequence_trp_tyr[ n-1 ];
	}


	std::cout << "Found trp and tyr at " << n_trp << " and " << n_tyr << std::endl;

	make_pose_from_sequence( pose_trp_tyr, new_sequence_trp_tyr, *rsd_set );
	make_pose_from_sequence( pose_trp, new_sequence_trp, *rsd_set );
	make_pose_from_sequence( pose_tyr, new_sequence_tyr, *rsd_set );
	make_pose_from_sequence( pose_null, new_sequence_null, *rsd_set );

	std::cout << new_sequence_trp_tyr << std::endl;
	std::cout << new_sequence_trp << std::endl;
	std::cout << new_sequence_tyr << std::endl;
	std::cout << new_sequence_null << std::endl;

	copy_dofs_match_atom_names( pose_trp_tyr, pose_input, res_map );
	copy_dofs_match_atom_names( pose_trp, pose_input, res_map );
	copy_dofs_match_atom_names( pose_tyr, pose_input, res_map );
	copy_dofs_match_atom_names( pose_null, pose_input, res_map );



	ScoreFunctionOP scorefxn = getScoreFunction();

	PoseOP start_pose_op;
	//Pose & start_pose = *start_pose_op;
	//start_pose.dump_pdb( "START.pdb" );

	Pose native_pose = pose_input;

	// no sampling -- just gly/ala
	std::string silent_file = silent_file_null;
	{
		std::string tag = "S_0";
		protocols::stepwise::protein::output_silent_struct( pose_null, start_pose_op, silent_file, tag );
	}
	pose = pose_null;
	(*scorefxn)( pose );
	scorefxn->show( pose );

	Size const N_SAMPLE( option[ n_sample]() );
	Size count( 0 );

	//trp only
	pose = pose_trp;
	pose.dump_pdb( "TRP.pdb" );
	(*scorefxn)( pose );
	scorefxn->show( pose );
	silent_file = silent_file_trp;
	for (Size i1 = 1; i1 <= N_SAMPLE; i1++ ) {
		for (Size j1 = 1; j1 <= N_SAMPLE; j1++ ) {

			pose.set_chi( 1, n_trp, protocols::stepwise::protein::get_rotamer_angle( i1, N_SAMPLE ) );
			pose.set_chi( 2, n_trp, protocols::stepwise::protein::get_rotamer_angle( j1, N_SAMPLE ) );

			( *scorefxn )( pose );
			add_chi_tags( pose, native_pose, n_trp, n_tyr );
			std::string const tag = "S_"+ lead_zero_string_of( ++count, 5 );
			protocols::stepwise::protein::output_silent_struct( pose, start_pose_op, silent_file, tag );

		}
	}

	//tyr only
	count = 0;
	pose = pose_tyr;
	(*scorefxn)( pose );
	scorefxn->show( pose );
	silent_file = silent_file_tyr;
	for (Size i1 = 1; i1 <= N_SAMPLE; i1++ ) {
		for (Size j1 = 1; j1 <= N_SAMPLE; j1++ ) {

			pose.set_chi( 1, n_tyr, protocols::stepwise::protein::get_rotamer_angle( i1, N_SAMPLE ) );
			pose.set_chi( 2, n_tyr, protocols::stepwise::protein::get_rotamer_angle( j1, N_SAMPLE ) );

			( *scorefxn )( pose );
			add_chi_tags( pose, native_pose, n_trp, n_tyr );
			std::string const tag = "S_"+ lead_zero_string_of( ++count, 5 );
			protocols::stepwise::protein::output_silent_struct( pose, start_pose_op, silent_file, tag );

		}
	}


	//trp and tyr
	count = 0;
	pose = pose_trp_tyr;
	(*scorefxn)( pose );
	scorefxn->show( pose );

	return; // early exit.
	silent_file = silent_file_trp_tyr;
	for (Size i1 = 1; i1 <= N_SAMPLE; i1++ ) {
		for (Size j1 = 1; j1 <= N_SAMPLE; j1++ ) {

			pose.set_chi( 1, n_trp, protocols::stepwise::protein::get_rotamer_angle( i1, N_SAMPLE ) );
			pose.set_chi( 2, n_trp, protocols::stepwise::protein::get_rotamer_angle( j1, N_SAMPLE ) );

			for (Size p1 = 1; p1 <= N_SAMPLE; p1++ ) {
				for (Size q1 = 1; q1 <= N_SAMPLE; q1++ ) {

					pose.set_chi( 1, n_tyr, protocols::stepwise::protein::get_rotamer_angle( p1, N_SAMPLE ) );
					pose.set_chi( 2, n_tyr, protocols::stepwise::protein::get_rotamer_angle( q1, N_SAMPLE ) );

					( *scorefxn )( pose );
					add_chi_tags( pose, native_pose, n_trp, n_tyr );
					std::string const tag = "S_"+ lead_zero_string_of( ++count, 5 );
					protocols::stepwise::protein::output_silent_struct( pose, start_pose_op, silent_file, tag );

				}
			}

		}
	}

	return;

}


///////////////////////////////////////////////////////////////////////
std::string
get_file_name( std::string const & silent_file, std::string const & tag )
{
	Size pos( silent_file.find( ".out" ) );
	std::string silent_file_sample( silent_file );
	silent_file_sample.replace( pos, 4, tag+".out" );
	return silent_file_sample;

}


///////////////////////////////////////////////////////////////////////
// [this should eventually be a MOVER]
///////////////////////////////////////////////////////////////////////
void
rebuild_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::pack;
	using namespace protocols::stepwise::protein;

	// A lot of the following might be better handled by a JobDistributor!?

	////////////////////////////////////////////////////
	//Read in sequence information and native
	////////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	//Read in desired fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const desired_sequence = fasta_sequence->sequence();

	PoseOP native_pose;
	bool native_exists( false );
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		std::string native_pdb_file  = option[ in::file::native ];
		io::pdb::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );
		native_exists = true;
	}

	////////////////////////////////////////////////////
	// Actual read in of any starting poses.
	////////////////////////////////////////////////////
	Pose pose;

	utility::vector1< std::string > start_tags, silent_files_in;
	if ( !option[start_from_scratch] ){
		if ( option[ in::file::s].user() ) {
			start_tags = option[ in::file::s ]();
		} else {
			// There better be some kind of silent file specified...
			if ( !option[ in::file::silent ].user() ) {
				utility_exit_with_message( "Specify -start_from_scratch, -s <pdb>, or -in:file:silent <file> -tags <tags>" );
			}
			silent_files_in = option[ in::file::silent ]();
			start_tags = option[ in::file::tags ]();
			if ( silent_files_in.size() != start_tags.size() ) utility_exit_with_message( "one tag per silent file" );
		}
	}

	////////////////////////////////////////////////////
	// Setup
	////////////////////////////////////////////////////
	StepWiseProteinPoseSetup stepwise_pose_setup( desired_sequence, start_tags /*could be empty vector*/, silent_files_in);
	if ( option[ n_terminus ] ) stepwise_pose_setup.set_n_terminus( true );
	if ( option[ c_terminus ] ) stepwise_pose_setup.set_c_terminus( true );
	if ( option[ add_peptide_plane ] ) stepwise_pose_setup.set_add_peptide_plane( true );
	if ( option[ no_sample_junction ] ) stepwise_pose_setup.set_sample_junction( false );

	stepwise_pose_setup.apply( pose );
	utility::vector1 < Size > moving_residues = stepwise_pose_setup.moving_residues();

	pose.dump_pdb( "after_setup.pdb" );

	// Constraints...
	ConstraintSetOP cst_set;
	if ( option[ cst_file ].user() ) {
		cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, pose );
		pose.constraint_set( cst_set );
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	std::string const silent_file = option[ out::file::silent  ]();
	std::string const silent_file_sample = get_file_name( silent_file, "_pack" );
	std::string const silent_file_centroid = get_file_name( silent_file, "_centroid" );

	/////////////////////////////////////////////
	// Screen that predefines (phi, psi, omega)  for moving residues
	// --> input for later mover carries out all the green packer moves.
	/////////////////////////////////////////////
// 	ScoreFunctionOP centroid_scorefxn = ScoreFunctionFactory::create_score_function( option[centroid_weights] );
// 	StepWiseProteinScreener stepwise_screener( moving_residues );
// 	stepwise_screener.set_rmsd_cutoff( option[ filter_rmsd ]() );
// 	stepwise_screener.set_n_sample( option[ n_sample ]() );
// 	stepwise_screener.set_filter_native_big_bins( option[ filter_native_big_bins ]  );
// 	if (native_exists) stepwise_screener.set_native_pose( native_pose );
// 	if ( option[ centroid_output ] ) stepwise_screener.set_silent_file( silent_file_centroid );
// 	if ( option[ centroid_screen] ) {
// 		stepwise_screener.set_centroid_screen( true );
// 		stepwise_screener.set_centroid_score_diff_cut( option[ centroid_score_diff_cut ] );

// 		if ( option[ ghost_loops ] ){
// 			// Trying a mode where loops are not included in scoring.
// 			// the idea was to cut out poses where secondary structure elements
// 			// are in contact -- energy less than a reference pose in which
// 			// the secondary structure elements are really far apart.
// 			stepwise_screener.set_ghost_loops( true );
// 			stepwise_screener.set_centroid_score_diff_cut( 0.0 );
// 			// if reference is "expanded", rg doesn't make sense.
// 			centroid_scorefxn->set_weight( rg, 0.0 );
// 			// disallow steric clashes beyond what it in expanded pose.
// 			stepwise_screener.set_apply_vdw_cut( true );
// 		}

// 		stepwise_screener.set_centroid_scorefxn( centroid_scorefxn );
// 		stepwise_screener.set_nstruct_centroid( option[ nstruct_centroid ]() );

// 	}

// 	stepwise_screener.apply( pose );

// 	pose.dump_pdb( "after_screen.pdb" );

// 	utility::vector1< MainChainTorsionSetList >const & main_chain_torsion_set_lists
// 		= stepwise_screener.main_chain_torsion_set_lists();

// 	//////////////////////////////////////////////////////////////////////////
// 	//////////////////////////////////////////////////////////////////////////
// 	// StepWiseProteinResidueSampler -- iterative and enumerative sampling
// 	//   of backbone degrees of freedom.
// 	//////////////////////////////////////////////////////////////////////////
// 	//////////////////////////////////////////////////////////////////////////
// 	ScoreFunctionOP pack_scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );
// 	if (pose.constraint_set()->has_constraints() )	pack_scorefxn->set_weight( atom_pair_constraint, 1.0 );
// 	StepWiseProteinResidueSampler stepwise_residue_sampler( moving_residues, main_chain_torsion_set_lists );
// 	if (native_exists) stepwise_residue_sampler.set_native_pose( native_pose );
// 	stepwise_residue_sampler.set_scorefxn( pack_scorefxn );
// 	stepwise_residue_sampler.set_silent_file( silent_file_sample /*useful for checkpointing*/ );
// 	stepwise_residue_sampler.apply( pose );


// 	/////////////////////////////
// 	// Cluster...
// 	/////////////////////////////
// 	// Have an option to read structure back in from disk?  may be useful for checkpointing.
// 	// For now just take in silent structs prepared by the stepwise_residue_sampler.
// 	protocols::stepwise::StepWiseClusterer stepwise_clusterer(  stepwise_residue_sampler.silent_file_data() );
// 	Size max_decoys( 400 );
// 	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
// 	stepwise_clusterer.set_max_decoys( max_decoys );
// 	stepwise_clusterer.set_cluster_radius(	option[ OptionKeys::cluster::radius ]()	);
// 	stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] ); // false by default
// 	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
// 	stepwise_clusterer.cluster();

// 		// Perhaps we should output decoys into a silent file at this point -- for checkpointing.

// 	if ( option[ minimize ] ){
// 		// Minimize...
// 		PoseList minimize_pose_list = stepwise_clusterer.clustered_pose_list();

// 		//StepWiseProteinFilterer stepwise_filterer;
// 		//		stepwise_filterer.set_final_number( option[ n_minimize ]()  );
// 		//		stepwise_filterer.filter( pose_list, minimize_pose_list );
// 		//		std::cout << "FILTER " << pose_list.size() << " " << minimize_pose_list.size() << std::endl;

// 		StepWiseProteinPoseMinimizer stepwise_pose_minimizer( minimize_pose_list, moving_residues );
//     ScoreFunctionOP minimize_scorefxn( core::scoring::getScoreFunction() );
// 		if (pose.constraint_set()->has_constraints() )	minimize_scorefxn->set_weight( atom_pair_constraint, 1.0 );
// 		stepwise_pose_minimizer.set_scorefxn( minimize_scorefxn );
// 		stepwise_pose_minimizer.set_silent_file( silent_file );
// 		//stepwise_pose_minimizer.set_constraint_set( cst_set );
// 		if (native_exists) stepwise_pose_minimizer.set_native_pose( native_pose );

// 		// also outputs to silent file.
// 		stepwise_pose_minimizer.apply( pose );

// 	} else {
// 		stepwise_clusterer.output_silent_file( silent_file );
// 	}

	//	protocols::viewer::clear_conformation_viewers();
}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// Don't really need this one anymore.
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void
cluster_outfile_test_OLD(){
	//eventually will allow for readin of several files, cluster across all of them.

	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace protocols::cluster;
	using namespace core::io::pose_stream;
	using namespace core::chemical;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	/////////////////////////////
	//static ScoreFunctionOP scorefxn = getScoreFunction();
	static ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	ClusterPhilStyleOP clustering( new ClusterPhilStyle );
	PoseInputStreamOP input  = new SilentFilePoseInputStream( option[ in::file::silent ]() );
	clustering->set_score_function( scorefxn );

	//	clustering->set_cluster_by_protein_backbone( true );
	//	clustering->set_cluster_by_all_atom( true );

	while ( input->has_another_pose() ) {
		core::pose::Pose pose;
		input->fill_pose( pose, *rsd_set );
		clustering->apply( pose );
	}

	clustering->set_cluster_radius(		option[ OptionKeys::cluster::radius ]()	);

	clustering->do_clustering();
	clustering->do_redistribution();
	clustering->sort_each_group_by_energy();
	//	clustering->sort_groups_by_energy();

	if ( option[ OptionKeys::cluster::remove_singletons ].user() ) clustering->remove_singletons();

	/////////////////////////////////
	// output to screen
	/////////////////////////////////
	clustering->print_summary();

	/////////////////////////////////
	// create silent file with the
	// lowest-energy poses from each cluster
	/////////////////////////////////
	std::vector< std::deque< int > > const & cluster_list( clustering->get_cluster_list() );
	std::vector< core::pose::Pose > & pose_list( clustering->get_pose_list() );
	std::vector< std::string > const & tag_list( clustering->get_tag_list() );

	BinaryProteinSilentStructOP s;

	static const SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	for ( Size n = 0 ; n < cluster_list.size(); n++ ) {
		Size const cluster_size = cluster_list[ n ].size();
		Size pose_list_index( 0 );
		if ( cluster_size == 1 ){
			pose_list_index = cluster_list[ n ][ 0 ];
		} else {
			pose_list_index = cluster_list[ n ][ 1 ];
		}

		pose::Pose & pose( pose_list[ pose_list_index ] );
		(*scorefxn)( pose );
		std::string const & tag( tag_list[ pose_list_index ] );

		s = new BinaryProteinSilentStruct( pose, tag );

		silent_file_data.write_silent_struct( *s, silent_file, false /*write score only*/ );

	}



}


///////////////////////////////////////////////////////////////
void
peptide_plane_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::io::pose_stream;
	using namespace core::pose;
	using namespace core::pack;
	using namespace core::id;


	// setup residue types
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	PoseInputStreamOP input  = new SilentFilePoseInputStream( option[ in::file::silent ]() );

	while ( input->has_another_pose() ) {
		core::pose::Pose pose;
		input->fill_pose( pose, *rsd_set );
		Residue const & rsd( pose.residue(1) );
		Real angle = angle_radians( rsd.xyz( "H" ), rsd.xyz( "N" ), rsd.xyz( "CP" ) );
		Real const angle_check = numeric::conversions::to_degrees( angle );
		std::cout << tag_from_pose( pose ) << " " << angle_check << std::endl;
	}

}


///////////////////////////////////////////////////////////////
void
peptide_plane_test_OLD(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::pack;
	using namespace core::id;


	// setup residue types
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	std::string const desired_sequence = "W";

	ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );

	Pose pose;
	make_pose_from_sequence( pose, desired_sequence, *rsd_set, false /*auto_termini*/);
	pose.set_chi( 1, 1, 60.0 );
	pose.set_chi( 2, 1, 60.0 );

	std::cout << " NO TERMINI" << std::endl;
	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );
	pose.dump_pdb( "no_termini.pdb" );

	{
		Pose termini_pose = pose;

		chemical::add_lower_terminus_type_to_pose_residue( termini_pose, 1   );
		chemical::add_upper_terminus_type_to_pose_residue( termini_pose, termini_pose.total_residue()   );

		std::cout << " ADDED TERMINI" << std::endl;
		(*scorefxn)( termini_pose );
		scorefxn->show( std::cout, termini_pose );
		termini_pose.dump_pdb( "add_termini.pdb" );

		// NEW! Now can only add N-acetylation, C-methylamidation to terminus variants.
		//		pose = termini_pose;

	}


	Pose peptide_plane_pose = pose;

	{

		if ( false ){
			// Might be useful in a util.cc
			Residue residue( pose.residue(1) );
			Size const atom_no( residue.atom_index( " H  " ));
			std::cout << "GOT N ACETYLATION? " << pose.residue(1).has_variant_type( "N_ACETYLATION" ) << std::endl;
			std::cout<< " " << 	pose.residue(1).type().icoor( atom_no ).stub_atom1().atomno() <<
				" " << pose.residue(1).type().icoor( atom_no ).stub_atom2().atomno() <<
				" " << pose.residue(1).type().icoor( atom_no ).stub_atom3().atomno() <<
				" " << pose.residue(1).type().icoor( atom_no ).is_internal() <<
				" " << pose.residue(1).type().icoor( atom_no ).depends_on_polymer_lower()
							 << std::endl;
		}

		chemical::add_variant_type_to_pose_residue( peptide_plane_pose, "N_ACETYLATION", 1   );

		if ( false ){
			//ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( 'W' ).set_property( "N_ACETYLATION" ).select( *rsd_set )[1] );
			ResidueTypeCOP new_rsd_type( rsd_set->get_residue_type_with_variant_added( pose.residue(1).type(), "N_ACETYLATION" ) );

			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );

			Size const atom_no( new_rsd->atom_index( " H  " ));
			std::cout << " X " << 	new_rsd->type().icoor( atom_no ).stub_atom1().atomno() <<
				" " << new_rsd->type().icoor( atom_no ).stub_atom2().atomno() <<
				" " << new_rsd->type().icoor( atom_no ).stub_atom3().atomno() <<
				" " << new_rsd->type().icoor( atom_no ).is_internal() <<
				" " << new_rsd->type().icoor( atom_no ).depends_on_polymer_lower()
								<< std::endl;

		}

		if ( false ){
			// Might be useful in a util.cc
			Residue residue( peptide_plane_pose.residue(1) );
			Size const atom_no( residue.atom_index( " H  " ));

			std::cout << "GOT N ACETYLATION? " << peptide_plane_pose.residue(1).has_variant_type( "N_ACETYLATION" ) << std::endl;
			std::cout<< " " << 	peptide_plane_pose.residue(1).type().icoor( atom_no ).stub_atom1().atomno() <<
				" " << peptide_plane_pose.residue(1).type().icoor( atom_no ).stub_atom2().atomno() <<
				" " << peptide_plane_pose.residue(1).type().icoor( atom_no ).stub_atom3().atomno() <<
				" " << peptide_plane_pose.residue(1).type().icoor( atom_no ).is_internal() <<
				" " << peptide_plane_pose.residue(1).type().icoor( atom_no ).depends_on_polymer_lower()
							 << std::endl;

			peptide_plane_pose.set_xyz( AtomID( atom_no, 1 ),
																	residue.build_atom_ideal( atom_no, peptide_plane_pose.conformation() ) );

		}


		//peptide_plane_pose.append_residue_by_bond( *new_rsd );

		std::cout << " ADDED PEPTIDE_PLANE" << std::endl;
		(*scorefxn)( peptide_plane_pose );
		scorefxn->show( std::cout, peptide_plane_pose );

		std::cout << "PHI " << peptide_plane_pose.phi( 1 ) << std::endl;
		peptide_plane_pose.dump_pdb( "add_peptide_plane.pdb" );

		peptide_plane_pose.set_phi( 1, 0.0 );
		(*scorefxn)( peptide_plane_pose );
		scorefxn->show( std::cout, peptide_plane_pose );
		std::cout << "PHI " << peptide_plane_pose.phi( 1 ) << std::endl;
		peptide_plane_pose.dump_pdb( "add_peptide_plane_setphi.pdb" );
	}

	{

		chemical::add_variant_type_to_pose_residue( peptide_plane_pose, "C_METHYLAMIDATION", 1   );

		std::cout << " ADDED PEPTIDE_PLANE2" << std::endl;
		(*scorefxn)( peptide_plane_pose );
		scorefxn->show( std::cout, peptide_plane_pose );

		std::cout << "PSI " << peptide_plane_pose.psi( 1 ) << std::endl;
		std::cout << "OMEGA " << peptide_plane_pose.omega( 1 ) << std::endl;
		peptide_plane_pose.dump_pdb( "add_peptide_plane2.pdb" );

		peptide_plane_pose.set_psi( 1, -50.0 );
		peptide_plane_pose.set_omega( 1, 180.0 );
		(*scorefxn)( peptide_plane_pose );
		scorefxn->show( std::cout, peptide_plane_pose );

		std::cout << "PSI " << peptide_plane_pose.psi( 1 ) << std::endl;
		std::cout << "OMEGA " << peptide_plane_pose.omega( 1 ) << std::endl;
		peptide_plane_pose.dump_pdb( "add_peptide_plane2_setpsi_omega.pdb" );
		}


	{

		chemical::remove_variant_type_from_pose_residue( peptide_plane_pose, "N_ACETYLATION", 1   );

		std::cout << " REMOVED PEPTIDE_PLANE3" << std::endl;
		(*scorefxn)( peptide_plane_pose );
		scorefxn->show( std::cout, peptide_plane_pose );

		peptide_plane_pose.dump_pdb( "remove_peptide_plane3.pdb" );

		chemical::add_variant_type_to_pose_residue( peptide_plane_pose, "N_ACETYLATION", 1   );
		(*scorefxn)( peptide_plane_pose );
		scorefxn->show( std::cout, peptide_plane_pose );

		peptide_plane_pose.dump_pdb( "add_peptide_plane4.pdb" );

	}

}


///////////////////////////////////////////////////////////////
void
cluster_outfile_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	protocols::stepwise::StepWiseClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
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


////////////////////////////////////////////////////////////////////
// This loop creates scorefiles that I'm going to plot in MATLAB
// so that I can visualize  all these screwy non-bonded terms
// in Rosetta's score12.
////////////////////////////////////////////////////////////////////
void
score12_plot_test()
{
	//eventually will allow for readin of several files, cluster across all of them.

	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace protocols::cluster;
	using namespace core::io::pose_stream;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	/////////////////////////////
	//ScoreFunctionOP scorefxn = getScoreFunction();
	static ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();

	///////////////////////////////////////
  // Go through each possible sequence
	///////////////////////////////////////

	std::string const protein_sequence = "ACDEFGHIKLMNPQRSTVWY";

	for ( Size i = 1; i <= 20 ; i++ ){
		//		AA aa( i );

		std::stringstream seqstream;
		seqstream << "G" << protein_sequence[i-1] << "G";
		std::string sequence;
		seqstream >> sequence;

		char const seqchar( protein_sequence[i-1] );
		//		ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( seqchar ).exclude_variants().select( *rsd_set )[1] );
		//		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );

		std::cout << "DOING: " << seqchar << std::endl;

		Pose pose;

		//		pose.append_residue_by_jump( *rsd, 1 );

		//		kinematics::FoldTree f( pose.total_residue() );
		//		f.new_jump( 1, 2, 1);
		//		pose.fold_tree( f );

		make_pose_from_sequence( pose, sequence, *rsd_set );

		// main loop
		Size const N_SAMPLE( option[ n_sample]() );
		Size count( 0 );
		for (Size i1 = 1; i1 <= N_SAMPLE; i1++ ) {
			for (Size j1 = 1; j1 <= N_SAMPLE; j1++ ) {

				Real const phi = protocols::stepwise::protein::get_rotamer_angle( i1, N_SAMPLE );
				Real const psi = protocols::stepwise::protein::get_rotamer_angle( j1, N_SAMPLE );
				pose.set_phi( 2, phi );
				pose.set_psi( 2, psi );

				//pose.set_phi( 1, phi );
				//pose.set_psi( 1, psi );

				//std::cout << "CHECK " << pose.phi( 2 ) << " " << pose.psi( 2 ) << std::endl;

				//Need to force a refold somehow
				pose.dump_pdb( "blah.pdb" );

				( *scorefxn )( pose );
				//scorefxn->show( std::cout, pose );

				std::string const tag = "S_"+ lead_zero_string_of( ++count, 5 );
				BinaryProteinSilentStruct s( pose, tag );
				s.add_energy( "phi", phi );
				s.add_energy( "psi", psi );

				static const SilentFileData silent_file_data;

				std::stringstream ss;
				std::string filename;
				ss  << seqchar << "_phi_psi.sc";
				ss  >> filename;
				silent_file_data.write_silent_struct( s, filename, true /*write score only*/ );

			}
		}

	}


	// phi/psi terms.

	////////////////////////////////////////
	// Then freeze at particular phi/psi,
	// and sample through chi's.

}



////////////////////////////////////////////////////////
// Move ths out of stepwise_protein_test.cc!!
////////////////////////////////////////////////////////
////////////////////////////////////////////////
void
eval_lj(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real const & d2,
	Real & fa_atr_score,
	Real & fa_rep_score,
	scoring::etable::Etable const & etable
)
{

	ObjexxFCL::FArray3D< Real > const & ljatr_ = etable.ljatr();
	ObjexxFCL::FArray3D< Real > const & ljrep_ = etable.ljrep();

	//	Real temp_score( 0.0 );
	fa_atr_score = 0.0;
	fa_rep_score = 0.0;

	Real const & safe_max_dis2_ =  etable.get_safe_max_dis2();
	Real const & get_bins_per_A2_ =  etable.get_bins_per_A2();

	if ( ( d2 < safe_max_dis2_) && ( d2 != Real(0.0) ) ) {

		Real const d2_bin = d2 * get_bins_per_A2_;
		int	disbin = static_cast< int >( d2_bin ) + 1;
		Real	frac = d2_bin - ( disbin - 1 );

		// l1 and l2 are FArray LINEAR INDICES for fast lookup:
		// [ l1 ] == (disbin  ,attype2,attype1)
		// [ l2 ] == (disbin+1,attype2,attype1)

		{
			int const l1 = ljatr_.index( disbin, atom2.type(), atom1.type() );
			int const l2 = l1 + 1;
			fa_atr_score = ( (1.-frac)* ljatr_[ l1 ] + frac * ljatr_[ l2 ]);
		}


		{
			int const l1 = ljrep_.index( disbin, atom2.type(), atom1.type() );
			int const l2 = l1 + 1;
			fa_rep_score = ( (1.-frac)* ljrep_[ l1 ] + frac * ljrep_[ l2 ]);
		}

	}

}

//////////////////////////////////////////////////////////////////////////////////////////
Real
get_lj_atom_score( core::id::AtomID const atom_id, core::pose::Pose const & pose,
									 scoring::etable::Etable const & etable ){

	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	Vector const heavy_atom_i( rsd1.xyz( m ) );

	Real score = 0.0;

	for ( Size j = 1; j < pose.total_residue(); j ++ ) {

		if ( i == j ) continue;

		conformation::Residue const & rsd2( pose.residue( j ) );

		using namespace scoring::etable::count_pair;
		CountPairFunctionOP cpfxn =
			CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

		for ( Size n = 1; n <= rsd2.natoms(); ++n ) {

			Real cp_weight = 1.0;

			if ( cpfxn->count(m, n, cp_weight ) ) {

				Vector const heavy_atom_j( rsd2.xyz( n ) );
				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Real const d2 = d_ij.length_squared();
				//Real const d = std::sqrt( d2 );
				Vector const d_ij_norm = d_ij.normalized();

				Real fa_atr( 0.0 ), fa_rep( 0.0 );
				eval_lj( rsd1.atom(m), rsd2.atom(n), d2, fa_atr, fa_rep, etable );

				// In principle could pass in an emap and weight the components
				// by the Emap.
				score += cp_weight * ( fa_atr + fa_rep );
      }
		}
	}

	return score;
}

////////////////////////////////////////////////////////
// Move ths out of stepwise_protein_test.cc!!
////////////////////////////////////////////////////////
void
color_by_lj_test()
{
	// Read in pdb.
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::io::pose_stream;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::farna;
	using namespace core::scoring::etable;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	ScoreFunctionOP scorefxn = getScoreFunction();

	EtableOP etable_ptr
		( new Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ),
									EtableOptions() ) );

	Pose pose;

	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
			);
		} else {
			input = new SilentFilePoseInputStream( option[ in::file::silent ]() );
		}
	} else {
		input = new PDBPoseInputStream( option[ in::file::s ]() );
	}

	utility::vector1< Size > sample_residues = option[ sample_res ]();

	Size count( 0 );
	while ( input->has_another_pose() ){

		input->fill_pose( pose, *rsd_set );

		PDBInfoOP pdb_info(  new PDBInfo( pose, true ) );

		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		if ( sample_residues.size() == 0 ) {
			for ( Size n = 1; n <= pose.total_residue(); n++ ) sample_residues.push_back( n );
		}

		for ( Size n = 1; n <= sample_residues.size(); n++ ) {

			Size const i = sample_residues[ n ];

			for (Size j = 1; j <= pose.residue( i ) .natoms(); j++ ) {
				Real const score =	get_lj_atom_score(  id::AtomID( j, i ), pose, *etable_ptr );
				//std::cout << i <<  " " << j << " " <<  score << std::endl;
				pdb_info->temperature( i, j, score );

			}
		}

		pose.pdb_info( pdb_info );

		std::string pdb_name =  "COLOR_" + string_of(count++) + ".pdb";
		std::cout << "About to output: " << pdb_name  << std::endl;
		pose.dump_pdb( pdb_name );

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ rebuild ] ){
		rebuild_test();
	} else if ( option[ cluster_test ] ){
		cluster_outfile_test();
	} else if ( option[ sample_trp ] ){
		sample_trp_test();
	} else if ( option[ sample_trp_tyr ] ){
		sample_trp_tyr_test();
	} else if ( option[ repack ] ){
		repack_test();
	} else if ( option[ score12_plot ] ){
		score12_plot_test();
	} else if ( option[ peptide_plane ] ){
		peptide_plane_test_OLD();
	} else if ( option[ minimize_test ] ){
		minimizer_test();
	} else if ( option[ color_by_lj ] ){
		color_by_lj_test();
	} else {
		sample_rama_test();
	}

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
	NEW_OPT( fullatom, "fullatom", false );
	NEW_OPT( minimize, "minimize", false );
	NEW_OPT( rebuild, "rebuild", false );
	NEW_OPT( start_from_scratch, "start from scratch", false );
	NEW_OPT( make_ideal_helix, "make_ideal_helix", false );
	NEW_OPT( cluster_test, "cluster", false );
	NEW_OPT( cluster_by_all_atom_rmsd, "cluster by all atom rmsd", false );
	NEW_OPT( sample_trp, "sample trp in trp cage", false );
	NEW_OPT( sample_trp_tyr, "sample trp an tyr in trp cage", false );
	NEW_OPT( repack, "repack", false );
	NEW_OPT( centroid_output, "output centroid structure during screening", false );
	NEW_OPT( sample_residue, "residue to sample", 5 );
	NEW_OPT( start_res, "starting residue", 0 );
	NEW_OPT( end_res, "ending residue", 0 );
	NEW_OPT( n_sample, "number of samples per torsion angle", 18 );
	NEW_OPT( filter_rmsd, "for fast sampling", -1.0 );
	NEW_OPT( trans_omega, "omega->180.0", false );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT( centroid_score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT( filter_native_big_bins, "Figure out various terms for score12", false );
	NEW_OPT( rename_tags, "After clustering, rename tags S_0, S_1, etc.", false );
	NEW_OPT( peptide_plane, "testing peptide plane extensions", false );
	NEW_OPT( n_terminus, "build N terminus", false );
	NEW_OPT( c_terminus, "build C terminus", false );
	NEW_OPT( minimize_test, "Minimize test", false );
	NEW_OPT( deriv_check, "derivative check for minimize test", false );
	NEW_OPT( pack_weights, "weights for green packing", "standard.wts" );
	NEW_OPT( centroid_weights, "weights for centroid filter", "score3.wts" );
	NEW_OPT( score12_plot, "Figure out various terms for score12", false );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );
	NEW_OPT( centroid_screen, "Centroid Screen", false );
	NEW_OPT( no_sample_junction, "disable sampling of residue at junction inherited from start pose", false );
	NEW_OPT( add_peptide_plane, "Include N-acetylation and C-methylamidation caps at termini", false );
	NEW_OPT( entropy_test, "Entropy test", false );
	NEW_OPT( ghost_loops, "Virtualize loops in centroid screening", false );
	NEW_OPT( nstruct_centroid, "Number of decoys to output from centroid screening", 0 );
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector ); //I am here.

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
