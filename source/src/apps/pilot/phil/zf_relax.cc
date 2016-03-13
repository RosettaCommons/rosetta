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
#include <devel/dna/ProteinDNA_Relax.hh>
#include <devel/dna/relax_util.hh>
#include <devel/dna/base_movers.hh>
#include <devel/cartesian_frags/DNA_FragLib.hh>

#include <protocols/loops/ccd_closure.hh>
#include <utility/excn/Exceptions.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>


#include <core/chemical/VariantType.hh>


#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>


#include <core/pose/Pose.hh>

#include <basic/options/util.hh>

#include <basic/prof.hh> // profiling

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>


#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>


// // C++ headers

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <ObjexxFCL/format.hh>

//#include <basic/options/keys/specificity.OptionKeys.gen.hh>


using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
//using namespace protocols;

using utility::vector1;
using std::string;


static THREAD_LOCAL basic::Tracer tt( "demo.phil.zf_relax", basic::t_trace );
static THREAD_LOCAL basic::Tracer td( "demo.phil.zf_relax", basic::t_debug );
static THREAD_LOCAL basic::Tracer ti( "demo.phil.zf_relax", basic::t_info );
static THREAD_LOCAL basic::Tracer tw( "demo.phil.zf_relax", basic::t_warning );


///////////////////////////////////////////////////////////////////////////////
namespace zf_relax_ns {

	Size const protein_root_mpos( 3 ); // dna->protein jump has downstream pos at motif pos=3
	Size const   dna_anchor_mpos( 2 ); // anchored at position 2 in each triplet

	Size const   cutpoint_mpos( 13 ); // gly
	Size const loop_begin_mpos( 12 ); // ser
	Size const   loop_end_mpos( 17 ); // phe
}


///////////////////////////////////////////////////////////////////////////////
int
protein_motif_offset( Size const finger, pose::Pose const & pose )
{
	return pose.fold_tree().downstream_jump_residue( finger ) - zf_relax_ns::protein_root_mpos;
}

///////////////////////////////////////////////////////////////////////////////
Real
zf_chainbreak_distance( Size const finger, pose::Pose const & pose )
{
	using namespace zf_relax_ns;

	Real const target_CN_distance( 1.33 ); // from ALA.params, icoor for upper,lower
	Size const pmotif_offset( protein_motif_offset( finger, pose ) );
	Size const cutpoint  ( pmotif_offset + cutpoint_mpos );
	assert( pose.fold_tree().is_cutpoint( cutpoint ) );
	Real const current_CN_distance( pose.residue(cutpoint).xyz("C").distance( pose.residue( cutpoint+1 ).xyz("N") ) );
	return std::abs( current_CN_distance - target_CN_distance );
}

///////////////////////////////////////////////////////////////////////////////
// returns true for success, false for failure


bool
ccd_close_zf_chainbreak( Size const finger, pose::Pose & pose )
{
	using namespace zf_relax_ns;

	/// ccd closure params
	Size const ccd_cycles( 100 );
	Real const tolerance( 0.1 );
	bool const rama_check( true );
	Real const max_rama_score_increase( 2.0 );
	Real const max_tdH( 30.0 ), max_tdS( 45.0 ), max_tdL( 60.0 );

	Size const pmotif_offset( protein_motif_offset( finger, pose ) );

	Size const cutpoint  ( pmotif_offset + cutpoint_mpos );
	Size const loop_begin( pmotif_offset + loop_begin_mpos );
	Size const loop_end  ( pmotif_offset + loop_end_mpos );

	assert( pose.fold_tree().is_cutpoint( cutpoint ) );

	kinematics::MoveMap mm;
	for ( Size i=loop_begin; i<= loop_end; ++i ) {
		mm.set( id::TorsionID( i, id::BB, 1 ), true ); // phi
		mm.set( id::TorsionID( i, id::BB, 2 ), true ); // psi
	}

	// try to close using ccd closure
	Real fwd_dev, bwd_dev, tor_delta, rama_delta;
	protocols::loops::fast_ccd_loop_closure( pose, mm, loop_begin, loop_end, cutpoint, ccd_cycles, tolerance, rama_check,
																			 max_rama_score_increase, max_tdH, max_tdS, max_tdL, fwd_dev, bwd_dev,
																			 tor_delta, rama_delta );

	tt << "after ccd: " << fwd_dev << ' ' << bwd_dev << ' ' << tor_delta << ' ' << rama_delta << std::endl;

	if ( fwd_dev < tolerance*2 && bwd_dev < tolerance*2 ) {
		// successfully closed!
		return true;
	} else {
		return false;
	}
}


///////////////////////////////////////////////////////////////////////////////

class ZF_PatchupMover : public protocols::moves::Mover {
public:
	typedef protocols::moves::MoverOP MoverOP;


public:
	//// ctor
	ZF_PatchupMover( MoverOP my_mover ):
		my_mover_( my_mover )
	{
		protocols::moves::Mover::type( my_mover_->type()+"PatchupMover" );
	}

	void
	apply( pose::Pose & pose )
	{
		Real const chainbreak_distance_delta_threshold( 10.0 ); // too big
		pose::Pose const start_pose( pose );
		Real const start_chainbreak1( zf_chainbreak_distance( 1, pose ) );
		Real const start_chainbreak2( zf_chainbreak_distance( 2, pose ) );

		Size ntries(0), max_tries( 20 );
		while ( true ) {
			++ntries;
			td << "ZF_PatchupMover::apply ntries: " << ntries << std::endl;

			if ( ntries > max_tries ) {
				tw << "too many tries in ZF_PatchupMover::apply! " << ntries << std::endl;
				break;
			}

			my_mover_->apply( pose );

			Real const delta1( zf_chainbreak_distance( 1, pose ) - start_chainbreak1 );
			Real const delta2( zf_chainbreak_distance( 2, pose ) - start_chainbreak2 );

			if ( delta1 > chainbreak_distance_delta_threshold || delta2 > chainbreak_distance_delta_threshold ) {
				// undo the move, try again
 				pose = start_pose;
				continue;
			}

			if ( delta1 > 0.1 ) {
				// the mover perturbed this chainbreak:
				if ( !ccd_close_zf_chainbreak( 1, pose ) ) {
					pose = start_pose;
					continue;
				}
			}

			if ( delta2 > 0.1 ) {
				// the mover perturbed this chainbreak:
				if ( !ccd_close_zf_chainbreak( 2, pose ) ) {
					pose = start_pose;
					continue;
				}
			}

			// if we got here we successfully closed both chainbreaks!
			break;
		}
	}

private:
	MoverOP my_mover_;


}; // ZF_PatchupMover


///////////////////////////////////////////////////////////////////////////////
void
setup_zf_pose(
	Size const nfinger,
	Size const dna_motif_offset,
	vector1< Size > const & prot_motif_offset,
	pose::Pose & pose
)
{
	using namespace zf_relax_ns;
	using namespace devel::dna;

	scoring::dna::set_base_partner( pose );

	// first let's delete unpaired dna residues
	delete_unpaired_bases( pose );

	//pose.dump_pdb( "no_unpaired.pdb" );

	Size const nres( pose.total_residue() );

	kinematics::FoldTree f( nres );

	Size protein_dna_cutpoint(0);
	for ( Size i=1; i< nres; ++i ) {
		if ( ( pose.residue(i  ).is_DNA() && pose.residue(i+1).is_protein() ) ||
				 ( pose.residue(i+1).is_DNA() && pose.residue(i  ).is_protein() ) ) {
			protein_dna_cutpoint = i;
			break;
		}
	}


	// first the nfinger jumps to the fingers:
	for ( Size i=1; i<= nfinger; ++i ) {
		Size const dna_anchor( dna_motif_offset + 3 * ( nfinger-i ) + dna_anchor_mpos );
		Size const protein_root( prot_motif_offset[i] + protein_root_mpos );
		if ( i<nfinger ) {
			// chainbreak after glycine
			f.new_jump( protein_root, dna_anchor, prot_motif_offset[i] + cutpoint_mpos );

		} else {
			// chainbreak at end of protein
			f.new_jump( protein_root, dna_anchor, protein_dna_cutpoint );
		}

		if ( i == (nfinger+1)/2 ) {
			f.reorder( dna_anchor );
		}
	}

	// now the intra-dna jumps
	add_dna_base_jumps_to_fold_tree( pose, f );

	pose.fold_tree( f );

	set_dna_jump_atoms( pose );

	{ // add cutpoint variants
		using namespace chemical;
		for ( Size i=1; i< nfinger; ++i ) {
			Size const cutpoint( protein_motif_offset(i,pose) + cutpoint_mpos );
			core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, cutpoint   );
			core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
		}
	}

	tt << "setup zf pose: final foldtree: " << pose.fold_tree() << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
/**

	 Test zinc finger relaxation protocol.

	 for starters, hard code 1aay, later extend to other guys, read Zinc atoms, calculate motif positions...


	 setup a foldtree for DNA fragment insertions with jumps from triplet position 2 to each finger,
	 and chainbreaks between the three fingers, after the glycine residue.

	 moves:
	  1. rigid-body perturbation of one of the fingers, followed by closure of the loop
		2. internal dna flex-move followed by closure of one of the fingers, if necessary

	 experiment with ccd closure ? to close the finger chainbreaks, or fragment-based closure? or
	 torsion/angle minimization as above?


**/

void
zf_relax_test()
{
	using namespace zf_relax_ns;
	using namespace pose;
	using namespace kinematics;
	using namespace scoring;
	using namespace id;

	using namespace devel::dna; // base movers
	using namespace devel::cartesian_frags; // DNA-fraglib

	// parameters
	bool const repack_dna( true ); // change ME !!!!!!!!!!!!
	bool const exclude_DNA_DNA( false );
	Real const energycut( 0.1 );
	Real const min_tol( 0.001 );
	Real rb_mover_trans_mag( 0.5 );
	Real rb_mover_rot_mag( 3.0 );
	Real base_mover_frag_dev_threshold( 2.0 );
	Size base_mover_max_tries( 20);
	Real base_mover_max_score_increase( 2.0 );

	// dont forget to write the code for flipping the dna foldtree around !!!!!!!!!!!!!!

	//// retrieve commandline options /////////////////
	std::string output_tag, score_function_file;
	Size nstruct( 10 );
	Size ninner( 20 );
	Size nouter( 10 );
	vector1< string > fraglib_files;//, input_files;

	{ // scope
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace basic::options::OptionKeys::dna;
		output_tag = option[ out::output_tag ];
		score_function_file = option[ specificity::score_function ];
		fraglib_files = option[ specificity::frag_files ]();
		//input_files = start_files();
		if ( option[ out::nstruct              ].user() ) nstruct = option[ out::nstruct ];
		if ( option[ specificity::n_inner ].user() ) ninner  = option[ specificity::n_inner ];
		if ( option[ specificity::n_outer ].user() ) nouter  = option[ specificity::n_outer ];
		if ( option[ specificity::params ].user() ) {
			vector1< Real > const params( option[ specificity::params ]() );
			base_mover_frag_dev_threshold = params[1];
			base_mover_max_tries          = Size( params[2] );
			base_mover_max_score_increase = params[3];
			rb_mover_trans_mag            = params[4];
			rb_mover_rot_mag              = params[5];
		}
	}

	// randomize the order of the input files:
	//numeric::random::random_permutation( input_files, numeric::random::rg() );

	// the simulation pose
	Pose pose;
	core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);

	protocols::viewer::add_conformation_viewer( pose.conformation(), "zf_relax_pose" );

	//assert( pose.residue( 1 ).is_protein() && pose.residue( nres ).is_DNA() ); // temporary

	// hard-coded positions -- all wrt the pdb after unpaired bases have been removed
	Size const nfinger(3);
	Size const dna_motif_offset( 85 ); // the start of the nfinger*3 basepair dna motif minus 1

	vector1< Size > pmotif_offset;
	pmotif_offset.push_back( 16 ); // motif position "0" in finger 1
	pmotif_offset.push_back( 44 ); //                              2
	pmotif_offset.push_back( 72 ); //                              3

	setup_zf_pose( nfinger, dna_motif_offset, pmotif_offset, pose );

	assert( pose.residue( dna_motif_offset+4 ).name1() == 't' );
	assert( pose.residue( dna_motif_offset+5 ).name1() == 'g' );
	assert( pose.residue( dna_motif_offset+6 ).name1() == 'g' );
	assert( pose.residue( pmotif_offset[1] + 7 ).name1() == 'H' );
	assert( pose.residue( pmotif_offset[2] + 7 ).name1() == 'H' );
	assert( pose.residue( pmotif_offset[3] + 7 ).name1() == 'H' );

	////// setup ///////////////////////////////////

	//// the constraints in the pose
	setup_dna_chainbreak_constraints( pose );
	//setup_zf_chainbreak_constraints( pose );


	//// fragments
	DNA_FragLibOP fraglib( new DNA_FragLib() );
	build_frag_libraries( fraglib_files, *fraglib );

	//// scorefunctions
	scoring::ScoreFunctionOP cst_scorefxn( new ScoreFunction() );
	cst_scorefxn->set_weight( atom_pair_constraint, 1.0 );
	cst_scorefxn->set_weight(     angle_constraint, 1.0 );
	cst_scorefxn->set_weight(           chainbreak, 0.5 ); // tune this!

	scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );

	//	scorefxn->energy_method_options().exclude_DNA_DNA( exclude_DNA_DNA );
	methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( exclude_DNA_DNA );
	scorefxn->set_energy_method_options( options );

	scorefxn->add_weights_from_file( score_function_file );
	scorefxn->set_weight( atom_pair_constraint, 1.0 );
	scorefxn->set_weight(     angle_constraint, 1.0 );
	scorefxn->set_weight(           chainbreak, 1.0 );

	Real const start_score( (*scorefxn)( pose ) );
	ti << "start_scores:" << std::endl;
	scorefxn->show( ti, pose );
	ti << std::endl; // flush??

	Pose const start_pose( pose );
	Size const nres( pose.total_residue() );

	for ( Size n=1; n<= nstruct; ++n ) {
		using namespace protocols;
		using namespace protocols::moves;

		pose.clear();
		pose = start_pose;

		////// more setup ///////////////////////////////////

		//// montecarlo
		MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, 0.8 ) );

		//// movemaps

		// rbmover movemap
		MoveMapOP protein_jumps_mm( new MoveMap() );
		for ( Size i=1; i<= nfinger; ++i ) {
			protein_jumps_mm->set_jump( i, true );
		}

		// minimizer movemap
		// freeze protein backbone except in loops during minimization
		MoveMapOP min_mm( new MoveMap() );
		min_mm->set_jump( true );
		min_mm->set_chi( true );
		min_mm->set_bb( false );
		for ( Size i=1; i<= nres; ++i ) {
			if ( pose.residue(i).is_DNA() ) {
				min_mm->set_bb( i, true );
			} else {
				min_mm->set_bb( i, false );
			}
		}
		for ( Size ii=1; ii< nfinger; ++ii ) {
			Size const pmo( protein_motif_offset( ii, pose ) );
			for ( Size i= pmo + loop_begin_mpos; i<= pmo + loop_end_mpos; ++i ) {
				min_mm->set( TorsionID( i, BB, 1 ), true ); // phi
				min_mm->set( TorsionID( i, BB, 2 ), true ); // psi
			}
		}

		//// tasks
		pack::task::PackerTaskOP pack_task( pack::task::TaskFactory::create_packer_task( pose ));

		pack_task->initialize_from_command_line();
		for ( Size i = 1; i <= nres; ++i ) {
			if ( pose.residue(i).is_protein() ) {
				pack_task->nonconst_residue_task( i ).restrict_to_repacking();
			} else {
				if ( repack_dna ) {
					pack_task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					pack_task->nonconst_residue_task( i ).prevent_repacking();
				}
			}
		}
		pack_task->or_include_current( true );

		// use more rotamers for rotamer trials, since its faster
		pack::task::PackerTaskOP rottrial_task( pack_task->clone() );

		for ( Size i = 1; i <= nres; ++i ) {
			if ( pose.residue(i).is_protein() ) {
				rottrial_task->nonconst_residue_task(i).or_ex1( true );
				rottrial_task->nonconst_residue_task(i).or_ex2( true );
				rottrial_task->nonconst_residue_task(i).or_ex1aro( true );
				rottrial_task->nonconst_residue_task(i).or_ex1aro_sample_level( pack::task::EX_THREE_THIRD_STEP_STDDEVS );
				rottrial_task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
			}
		}

		// movers
		MoverOP bp_mover( new ZF_PatchupMover( new BasePairMover( fraglib,
																															base_mover_frag_dev_threshold,
																															base_mover_max_tries,
																															base_mover_max_score_increase,
																															cst_scorefxn ) ) );
		MoverOP bs_mover( new ZF_PatchupMover( new BaseStepMover( fraglib,
																															base_mover_frag_dev_threshold,
																															base_mover_max_tries,
																															base_mover_max_score_increase,
																															cst_scorefxn ) ) );
		MoverOP rb_mover( new ZF_PatchupMover( new devel::dna::RB_Mover( protein_jumps_mm, rb_mover_trans_mag,
																																		 rb_mover_rot_mag ) ) );
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( min_mm, scorefxn, "lbfgs_armijo_nonmonotone", min_tol, true ) );
		protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( scorefxn, pack_task, 25 ) );

		// rotamer trials w/ energycut
		protocols::simple_moves::EnergyCutRotamerTrialsMoverOP rottrial_mover
			( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn, *rottrial_task, mc, energycut ) );

		protocols::simple_moves::RotamerTrialsMoverOP full_rottrial_mover
			( new protocols::simple_moves::RotamerTrialsMover( scorefxn, *rottrial_task ) );

		// trials:
		TrialMoverOP   rb_min_trial = devel::dna::setup_MCM_trial(   rb_mover,      rottrial_mover, min_mover, mc );
		TrialMoverOP   bp_min_trial = devel::dna::setup_MCM_trial(   bp_mover,      rottrial_mover, min_mover, mc );
		TrialMoverOP   bs_min_trial = devel::dna::setup_MCM_trial(   bs_mover,      rottrial_mover, min_mover, mc );
		TrialMoverOP pack_min_trial = devel::dna::setup_MCM_trial( pack_mover, full_rottrial_mover, min_mover, mc );

		TrialMoverOP min_trial = new TrialMover( min_mover, mc );

		// initial min-trial
		min_trial->apply( pose );
		core::Real const aftermin_rms  ( scoring::all_atom_rmsd( pose, start_pose ) );
		core::Real const aftermin_score( (*scorefxn)( pose ) );
		ti << "aftermin_scores:" << std::endl;
		scorefxn->show( ti, pose );
		ti << std::endl;

		// initial pack
		pack_min_trial->apply( pose );
		ti << "afterpack_scores:" << std::endl;
		scorefxn->show( ti, pose );
		ti << std::endl;

		basic::prof_reset();
		for ( Size outer=1; outer<= nouter; ++outer ) {
			if ( outer>1 ) pack_min_trial->apply( pose );

			for ( Size inner=1; inner<= ninner; ++inner ) {
				bp_min_trial->apply( pose );
				bs_min_trial->apply( pose );
				rb_min_trial->apply( pose );
			}

			mc->recover_low( pose );

			// flip the foldtree
			//core::Real const score_before( (*scorefxn)( pose ) );
			//setup_dna_only_fold_tree( pose, outer%2 == 1 );
			//core::Real const score_after( (*scorefxn)( pose ) );
			//assert( std::abs( score_before - score_after ) < 1e-1 );
			//mc->reset( pose );

			basic::prof_show();
		} // outer

		// write the pdb
		string const outfilename( output_tag + "_final" + lead_zero_string_of( n, 4 )+".pdb" );
		pose.dump_scored_pdb( outfilename, *scorefxn );
		basic::prof_show();

		// status output
		mc->show_counters();
		ti << "final_scores:" << std::endl;
		scorefxn->show( ti, pose );
		ti << std::endl;
		core::Real const final_rms( scoring::all_atom_rmsd( pose, start_pose ) );
		core::Real const final_score( (*scorefxn)(pose) );
		{
			using namespace ObjexxFCL::format;
			std::cout << "final_rmsd_and_score: " << outfilename << F(9,3,aftermin_rms) << F(9,3,final_rms) <<
				F(9,3,start_score) << F(9,3,aftermin_score) << F(9,3,final_score) << std::endl;
		}

	} // nstruct loop

}

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	// silly
	zf_relax_test();
	exit(0);
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
