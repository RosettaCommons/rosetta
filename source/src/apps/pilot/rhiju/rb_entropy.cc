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

// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/rna/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/recces/util.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/basic_sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/stepwise/modeler/util.hh> //has euler angle stuff.
#include <protocols/toolbox/sample_around/util.hh> //has euler angle stuff.
#include <protocols/toolbox/rigid_body/util.hh>

#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/sample_around.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/constants.hh>

#include <utility/stream_util.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "rb_entropy" );

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;

using utility::vector1;
typedef  numeric::xyzMatrix< Real > Matrix;

// all helper functions moved to protocols/toolbox/sample_around/util.cc
using namespace protocols::toolbox::sample_around;

OPT_KEY( String, nucleobase )
OPT_KEY( Boolean, twodimensional )
OPT_KEY( Real, xyz_size )

OPT_KEY( Boolean, recces_mode )
OPT_KEY( Integer, n_cycle )
OPT_KEY( Real, a_form_range )
OPT_KEY( RealVector, temps )
OPT_KEY( RealVector, st_weights )
OPT_KEY( Real, rmsd_cutoff )
OPT_KEY( Real, translation_mag )
OPT_KEY( Real, rotation_mag )
OPT_KEY( String, out_prefix )
OPT_KEY( Boolean, save_score_terms )
OPT_KEY( Boolean, dump_pdb )
OPT_KEY( Integer, n_intermediate_dump )
OPT_KEY( Boolean, block_stack )

//////////////////////////////////////////////////////////////////
//
// Can we calculate free energy of forming a base pair
//  to start a helix ("init" in the nearest neighbor rules)?
//
// Basic tests of:
//
//  1. Fast RMSD calculation
//  2. Analytical calculation of RMSD distributions at "1 M" reference state.
//  3. Movers to install in thermal_sampler.
//
//        -- rhiju, 2016
//
//////////////////////////////////////////////////////////////////
core::Real
calc_base_centroid_rmsd( core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2 )
{
	Real rmsd( 0.0 ); Size numatoms( 0 );
	for ( Size i = rsd1.first_sidechain_atom() + 1; i <= rsd1.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.
		if ( rsd1.is_virtual( i ) ) continue;
		if ( rsd1.is_repulsive( i ) ) continue;
		Vector dist = ( rsd1.xyz(i) - rsd2.xyz(i) );
		rmsd += dist.length_squared();
		numatoms++;
	}
	rmsd = sqrt( rmsd/numatoms );
	return rmsd;
}

void
print_base_centroid_atoms( core::conformation::Residue const & rsd, std::string filename_xyz  )
{
	utility::io::ozstream out_xyz;
	out_xyz.open( filename_xyz );
	for ( Size i = rsd.first_sidechain_atom() + 1; i <= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.
		//  TR << rsd.atom_name( i ) << std::endl;
		if ( rsd.is_virtual( i ) ) continue;
		out_xyz << rsd.xyz(i).x() << " " << rsd.xyz(i).y() << " " << rsd.xyz(i).z() << std::endl;
	}
	out_xyz.close();
	TR << "Outputted xyz values into " << filename_xyz << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////
void
rb_entropy_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;
	using namespace core::kinematics;
	using numeric::random::rg;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose;
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	if ( option[ in::file::s ].user() ) {
		std::string infile  = option[ in::file::s ][1];
		import_pose::pose_from_file( pose, *rsd_set_op, infile , core::import_pose::PDB_file);
	} else {
		std::string const sequence = option[ nucleobase ]();
		runtime_assert( sequence.size() == 1 );
		make_pose_from_sequence( pose, sequence, *rsd_set_op, false /*auto_termini*/ );
		core::pose::rna::apply_Aform_torsions( pose, 1 );
	}

	pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, 1 );
	pose::add_variant_type_to_pose_residue( pose, VIRTUAL_RIBOSE, 1 );

	rotate_into_nucleobase_frame( pose );
	std::string const out_prefix = option[ out::prefix]();
	pose.dump_pdb( std::string(option[ out::path::path]()) + "/" + out_prefix +  "a_rotated.pdb" );

	add_virtual_res( pose );

	kinematics::FoldTree f ( pose.fold_tree() );
	std::cout << pose.annotated_sequence() << std::endl;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );


	Size cycles( 1000000 );
	core::pose::Pose start_pose = pose;
	Vector start_base_centroid = core::chemical::rna::get_rna_base_centroid( start_pose.residue( 1 ) );
	Size const probe_jump_num( 1 );
	Stub const upstream_stub = pose.conformation().upstream_jump_stub( probe_jump_num );
	kinematics::Jump start_jump( start_pose.jump( probe_jump_num ) );
	TR << "original centroid: " << start_base_centroid.x() << " " << start_base_centroid.y() << " " << start_base_centroid.z() << std::endl;

	// to compute moments of inertia in MATLAB
	print_base_centroid_atoms( pose.residue(1), std::string( option[ out::path::path]() ) + "/xyz.txt" );
	Matrix M;
	Real const & box_size = option[ xyz_size ]();

	utility::io::ozstream out;
	std::string filename( std::string( option[ out::path::path]() ) + "/rmsd.txt" );
	out.open( filename );
	for ( Size n = 1; n <= cycles; n++ ) {

		kinematics::Jump jump( start_jump );

		// rotation
		if ( option[ twodimensional]() ) {
			// 2D rotations, just for testing...
			M = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), 2 * numeric::constants::d::pi * rg().uniform() );
		} else {
			// 3D rotation -- thank you will sheffler & quaternions
			M = numeric::random::random_rotation();
		}
		jump.rotation_by_matrix( upstream_stub, start_base_centroid, M );

		// translation
		Vector random_translation( 2.0 * rg().uniform() - 1.0,
			2.0 * rg().uniform() - 1.0,
			2.0 * rg().uniform() - 1.0 );
		random_translation *= box_size;
		jump.set_translation( jump.get_translation() +  random_translation );

		pose.set_jump( probe_jump_num, jump );

		Real rmsd = calc_base_centroid_rmsd( pose.residue(1), start_pose.residue(1) );
		out << rmsd << std::endl;

		// tests that base_centroid is in same place:
		// base_centroid = get_rna_base_centroid( pose.residue( 1 ) );
		// TR << "new centroid:     " << base_centroid.x() << " " << base_centroid.y() << " " << base_centroid.z() << std::endl;

		if ( n <= 10 ) pose.dump_pdb( std::string( option[ out::path::path]() ) + "/test"+ObjexxFCL::lead_zero_string_of( n, 4 )+".pdb" );

	}
	out.close();
	TR << "Outputted " << cycles << " RMSD values into " << filename << std::endl;

}


//////////////////////////////////////////////////////////////////////////////
// horrific -- this is a direct copy of recces_turner code from Fang, from
//  the `solve_challenges` branch from early in 2016 -- essentially the
//  exact code that Fang checked in. There has been some work by AMW to refactor
//  the code since then, but it appears fairly broken. So I am going to
//  develop separately and then ask AMW to integrate.
//////////////////////////////////////////////////////////////////////////////
void
MC_run () {
	using namespace protocols::moves;
	using namespace scoring;
	using namespace core::chemical;
	using namespace core::chemical::rna;
	using namespace core::kinematics;
	using namespace core::pose;
	using namespace numeric;
	using namespace numeric::random;
	using namespace protocols::recces;

	clock_t const time_start( clock() );

	utility::vector1<Real> const & temps_( option[ temps ]() );
	runtime_assert( temps_.size() != 0 );

	utility::vector1<Real> weights_;
	utility::vector1<Real> const & orig_weights( option[ st_weights ]() );
	if ( temps_.size() != orig_weights.size() ) {
		weights_.push_back( 0 );
	}
	weights_.insert( weights_.end(), orig_weights.begin(),
		orig_weights.end() );
	runtime_assert( temps_.size() == weights_.size() );

	Size const n_cycle_( option[n_cycle]() );

	// Score function setup
	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/turner_no_zeros.wts" );
	}

	// Pose setup
	Pose pose;
	// should be starter base pair (e.g., "cg.pdb") -- only 2 residues!
	std::string infile  = option[ in::file::s ][1];
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	import_pose::pose_from_file( pose, *rsd_set_op, infile , core::import_pose::PDB_file);
	TR << "Annotated sequence of pose from " << infile << ": " << pose.annotated_sequence() << std::endl;
	kinematics::FoldTree f( 2 );
	f.new_jump( 1, 2, 1 );
	f.set_jump_atoms( 1, default_jump_atom(pose.residue_type(1)), default_jump_atom(pose.residue_type(2) )) ;
	pose.fold_tree( f );
	Size const probe_jump_num( 1 );
	Size const moving_rsd( 2 );
	print_base_centroid_atoms( pose.residue(moving_rsd), std::string( option[ out::path::path]() ) + "/xyz.txt" );
	for ( Size n = 1; n <= 2; n++ ) {
		pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, n );
		pose::add_variant_type_to_pose_residue( pose, VIRTUAL_RIBOSE, n );
		if ( option[ block_stack ]() ) {
			pose::add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, n );
			pose::add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, n );
		}
	}

	// Simulated Tempering setup
	SimulatedTempering tempering( pose, scorefxn, temps_, weights_ );
	// tempering.set_rep_cutoff( 100 );
	// set_gaussian_stdev( sampler, tempering, bp_rsd, dangling_rsd );

	// Setup for data saving for output
	Size curr_counts( 1 );
	utility::vector1<float> scores;
	update_scores( scores, pose, scorefxn );
	utility::vector1<float> const null_arr_;
	utility::vector1<utility::vector1<float> > data(
		temps_.size(), null_arr_ );

	Real const min( -100.05 ), max( 800.05 ), spacing( 0.1 );
	Histogram null_hist( min, max, spacing);
	utility::vector1<Histogram> hist_list( temps_.size(), null_hist );

	// Useful coounters and variables during the loop
	Size n_accept_total( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	Size const n_t_jumps( n_cycle_ / t_jump_interval );
	Size temp_id( tempering.temp_id() );
	bool const is_save_scores( option[save_score_terms]() );

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );
	Pose start_pose = pose;

	// SimulatedTempering doesn't handle restoration of pose when it rejects!!!
	Pose accepted_pose = pose;

	std::cout << "Start the main sampling loop." << std::endl;
	if ( option[dump_pdb]() ) pose.dump_pdb( "init.pdb" );

	Size const n_dump = option[n_intermediate_dump]();
	Size curr_dump = 1;

	// Main sampling cycle
	Real const base_centroid_rmsd_cutoff = option[ rmsd_cutoff ]() /* 2.0 */;
	Real const translation_mag_( option[ translation_mag]() /* 0.01 */);
	Real const rotation_mag_( option[ rotation_mag ]() /* 1.0 degrees */ );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	for ( Size n = 1; n <= n_cycle_; ++n ) {

		// let's do the move -- later stick into a StepWiseSampler for inclusion in RECCES.
		// (Perhaps StepWiseRigidBodySampler? I think that uses base centroids.)
		{
			kinematics::Jump jump( pose.jump( probe_jump_num ) );
			Vector base_centroid = get_rna_base_centroid( pose.residue( moving_rsd ) );
			// following few lines are from protocols/rigid/UniformRigidBodyMover.cc -- I'm not sure if
			//  its "kosher", however. Would be better to write myself in terms of, say, Euler angles or quaternions.
			core::Real theta = random_rotation_angle<core::Real>( rotation_mag_, numeric::random::rg() );
			xyzVector<core::Real> axis = random_point_on_unit_sphere<core::Real>( numeric::random::rg() );
			xyzMatrix<core::Real> delta_rot = rotation_matrix_radians( axis, theta );
			Stub const upstream_stub = pose.conformation().upstream_jump_stub( probe_jump_num );
			jump.rotation_by_matrix( upstream_stub, base_centroid, delta_rot );
			// translation
			xyzVector<core::Real> delta_trans = random_translation( translation_mag_, numeric::random::rg() );
			jump.set_translation( jump.get_translation() +  delta_trans );
			pose.set_jump( probe_jump_num, jump );
		}

		// later, may need to create 'ghost' rsd at ideal location, or use local2global stuff.
		// (for speed, could also use moment-of-inertia-based short cut)
		Real rmsd( calc_base_centroid_rmsd( pose.residue(moving_rsd), start_pose.residue(moving_rsd) ) );
		//  TR << "RMSD: " << rmsd << " SCORE: " << (*scorefxn)( pose ) << std::endl;
		if ( rmsd > base_centroid_rmsd_cutoff ) tempering.force_next_move_reject(); // artificial "wall"

		if ( ( tempering.check_boltzmann( pose ) ) || n == n_cycle_ ) {
			if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
			++n_accept_total;
			hist_list[temp_id].add( scores[1], curr_counts );
			update_scores( scores, pose, scorefxn );
			accepted_pose = pose;
			if ( n == n_cycle_ ) break;
			// sampler.update(); // What was this for?
			curr_counts = 1;
			if ( option[dump_pdb]() && scores[1] < min_score ) {
				min_score = scores[1];
				min_pose = pose;
			}
			if ( n_dump != 0 && n * (n_dump + 1) / double(n_cycle_) >= curr_dump ) {
				std::ostringstream oss;
				oss << "intermediate" << '_' << curr_dump << ".pdb";
				pose.dump_pdb(oss.str());
				++curr_dump;
			}
		} else {
			++curr_counts;
			pose = accepted_pose;
		}

		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;
			if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
			hist_list[temp_id].add( scores[1], curr_counts );
			curr_counts = 1;
			//   set_gaussian_stdev( sampler, tempering, bp_rsd, dangling_rsd );
			temp_id = tempering.temp_id();
		}
	}
	if ( option[dump_pdb]() ) {
		pose.dump_pdb( "end.pdb" );
		min_pose.dump_pdb( "min.pdb" );
		scorefxn->show( min_pose );
	}

	///////////////////////////////////////////////////////////////////////////
	// aargh -- hacky helper function to figure out scoretypes & weights to use in python --
	// AMW & cooper will fix this --
	utility::vector1<core::scoring::ScoreType> const & score_types( scorefxn->get_nonzero_weighted_scoretypes() );
	for ( core::Size i = 1; i<= score_types.size(); ++i ) TR << score_types[i] << std::endl;
	TR << " weights [" ;
	for ( core::Size i = 1; i<= score_types.size(); ++i ) {
		TR << " " << scorefxn->get_weight( score_types[ i ] );
		if ( i < score_types.size() ) TR << ",";
	}
	TR << " ] " << std::endl;
	///////////////////////////////////////////////////////////////////////////

	// Output simple statistics and the data
	//pose.dump_pdb("final.pdb");

	std::cout << "n_cycles: " << n_cycle_ << std::endl;
	std::cout << "Accept rate: " << double( n_accept_total ) / n_cycle_
		<< std::endl;
	std::cout << "T_jump accept rate: " << double( n_t_jumps_accept ) / n_t_jumps
		<< std::endl;
	Real const time_in_test = static_cast<Real>( clock() - time_start )
		/ CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;

	for ( Size i = 1; i <= temps_.size(); ++i ) {
		if ( is_save_scores ) {
			std::ostringstream oss;
			oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
				<< temps_[i] << ".bin.gz";
			Size const data_dim2( scorefxn->get_nonzero_weighted_scoretypes().size() + 2 );
			Size const data_dim1( data[i].size() / data_dim2 );
			vector2disk_in2d( oss.str(), data_dim1, data_dim2, data[i] );
		}
		std::ostringstream oss;
		oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
			<< temps_[i] << ".hist.gz";
		utility::vector1<Size> const & hist( hist_list[i].get_hist() );
		utility::vector1<Real> const & scores( hist_list[i].get_scores() );
		vector2disk_in1d( oss.str(), hist );

		std::ostringstream oss1;
		oss1 << option[out_prefix]() << "_hist_scores.gz";
		vector2disk_in1d( oss1.str(), scores );
	}
}






///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	if ( option[ recces_mode ]() ) {
		MC_run();
	} else {
		rb_entropy_test();
	}
	protocols::viewer::clear_conformation_viewers();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		utility::vector1< Real > null_real_vector;

		NEW_OPT( nucleobase, "nucleobase to sample around", "a" );
		NEW_OPT( twodimensional, "choose random rotation about Z-axis (2D)", false );
		NEW_OPT( xyz_size, "box half-diameter (max in x,y,z)", 1.0 );

		// for monte carlo
		NEW_OPT( recces_mode, "try RECCES-style sampling of a base pair", false );
		NEW_OPT( n_cycle, "cycle number for random sampling", 10000 );
		NEW_OPT( rmsd_cutoff, "base-centroid RMSD cutoff", 2.0 );
		NEW_OPT( translation_mag, "magnitude of random translation moves", 0.01 );
		NEW_OPT( rotation_mag, "magnitude of random rotation moves (in degrees)", 1.0 );
		NEW_OPT( temps, "Simulated tempering temperatures", null_real_vector );
		NEW_OPT( st_weights, "Simulated tempering weights", null_real_vector );
		NEW_OPT( out_prefix, "prefix for the out file", "turner" );
		NEW_OPT( save_score_terms,
			"Save scores and individual score terms"
			" of all sampled conformers", false );
		NEW_OPT( dump_pdb, "Dump pdb files", false );
		NEW_OPT( n_intermediate_dump,
			"Number of intermediate conformations to be dumped", 0 );
		NEW_OPT( block_stack, "add block stack pseudoatoms to block stacking", false );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
