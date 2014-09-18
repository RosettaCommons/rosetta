// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Frank DiMaio

#include <protocols/hybridization/MRMover.hh>
#include <protocols/hybridization/HybridizeProtocol.hh>
#include <protocols/electron_density/util.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <protocols/relax/util.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/electron_density/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCD.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.hh>

#include <protocols/relax/RelaxProtocolBase.hh>

#include <protocols/rbsegment_relax/AutoRBRelaxMover.hh>
#include <protocols/rbsegment_relax/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/util.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/pose/symmetry/util.hh>


#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/Residue.hh>
#include <utility/string_util.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;


namespace protocols {
namespace hybridization {

static thread_local basic::Tracer TR( "protocols.electron_density.util" );

using namespace protocols;
using namespace core;

// helper func
core::Size
parse_res( core::pose::Pose const &pose, std::string resnum ) {
	core::Size num;
	char chain;
	std::string::const_iterator input_end = resnum.end(), number_start = resnum.begin(), number_end = resnum.begin();
	while( number_end != input_end && *number_end >= '0' && *number_end <= '9' )
		++number_end;
	std::string num_str(number_start,number_end);
	num = std::atoi( num_str.c_str() );
	if (number_end == input_end) {
		chain = pose.pdb_info()->chain(1);
	} else {
		chain = *number_end;
	}
	return pose.pdb_info()->pdb2pose( chain, num );
}

//
//

MRMover::MRMover() :
		fragments_big_(NULL),
		fragments_small_(NULL) {
	init();
}


void
MRMover::init(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	max_gaplength_to_model_ = 4;
	cen1_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("score3");
	cen2_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
	fa_scorefxn_ = core::scoring::get_score_function();

	// default is a relatively short relax
	relax_max_iter_ = 200;
	relax_cycles_ = 2;
	censcale_ = 1.0;

	if (option[ OptionKeys::relax::default_repeats ].user())
		relax_cycles_ = option[ OptionKeys::relax::default_repeats ]();
	if (option[ OptionKeys::optimization::default_max_cycles ].user())
		relax_max_iter_ = option[ OptionKeys::optimization::default_max_cycles ]();

	// initialize scorefunctions
	cen1_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.25 );
	cen2_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.25 );
	fa_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.0 );

	// use reasonable defaults
	if (option[ OptionKeys::edensity::mapfile ].user()) {
		cen1_scorefxn_->set_weight( core::scoring::elec_dens_fast, 4.0 );
		cen2_scorefxn_->set_weight( core::scoring::elec_dens_fast, 4.0 );
		fa_scorefxn_->set_weight( core::scoring::elec_dens_window, 1.0 );
		core::scoring::electron_density::getDensityMap().setWindow( 5 );
	}
}

//
// apply()
void MRMover::apply( Pose &pose ) {
	using namespace protocols::loops;
	using namespace protocols::jd2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");

	bool threaded = true;
	protocols::comparative_modeling::ThreadingJobCOP job = dynamic_cast< protocols::comparative_modeling::ThreadingJob const*  >(
		JobDistributor::get_instance()->current_job()->inner_job().get() );
	if ( !job ) {
		if (option[ OptionKeys::in::file::fasta ].user()) {
			utility_exit_with_message(
				"CORE ERROR: You must use the ThreadingJobInputter with the LoopRelaxThreadingMover "
				"- did you forget the -in:file:template_pdb option?" );
		}

		threaded = false;

		// make sure gaps were read correctly
		protocols::simple_moves::MissingDensityToJumpMover read_miss_dens;
		read_miss_dens.apply( pose );
	}

	// we do a nonideal relax so make sure fa scorefunction is setup for that
	if (   fa_scorefxn_->get_weight( core::scoring::cart_bonded ) == 0
	    && fa_scorefxn_->get_weight( core::scoring::cart_bonded_angle ) == 0
	    && fa_scorefxn_->get_weight( core::scoring::cart_bonded_length ) == 0
	    && fa_scorefxn_->get_weight( core::scoring::cart_bonded_torsion ) == 0 ) {
		fa_scorefxn_->set_weight( core::scoring::cart_bonded, 0.5 );
		fa_scorefxn_->set_weight( core::scoring::pro_close, 0.0 );
	}

	// set up initial loop build
	LoopsOP my_loops;
	if (threaded) {
		core::Size nres = pose.total_residue();
		while (!pose.residue(nres).is_polymer()) nres--;
		my_loops = new Loops( job->loops( nres ) );

		if ( max_gaplength_to_model_ < 999 ) {
			trim_target_pose( pose, *my_loops, max_gaplength_to_model_ );
		}
	}

	if (max_gaplength_to_model_ > 0) {
		to_centroid.apply(pose);

		// setup & call single-template hybridize
		// 1 - remove loops from the input pose, add as template
		utility::vector1< int > pdb_numbering;
		utility::vector1< char > pdb_chains;
		core::pose::PoseOP template_pose = new core::pose::Pose;
		bool add_by_jump = true;
		for (Size i=1; i<=pose.total_residue(); ++i) {
			if (!threaded || !my_loops->is_loop_residue(i)) {
				if (add_by_jump) {
					if (template_pose->total_residue() > 0
								 && !template_pose->residue(template_pose->total_residue()).is_upper_terminus()
								 && template_pose->residue(template_pose->total_residue()).is_polymer())
						core::pose::add_upper_terminus_type_to_pose_residue( *template_pose, template_pose->total_residue() );

					template_pose->append_residue_by_jump( pose.residue(i), template_pose->total_residue(), "", "", true );
					add_by_jump = (!pose.residue(i).is_polymer() || pose.residue(i).is_upper_terminus());
					if (template_pose->residue(template_pose->total_residue()).is_polymer()
								 && !template_pose->residue(template_pose->total_residue()).is_lower_terminus() )
						core::pose::add_lower_terminus_type_to_pose_residue( *template_pose, template_pose->total_residue() );
				} else if ( !pose.residue(i).is_polymer() ) {
					if (template_pose->total_residue() > 0
								 && !template_pose->residue(template_pose->total_residue()).is_upper_terminus()
								 && template_pose->residue(template_pose->total_residue()).is_polymer())
						core::pose::add_upper_terminus_type_to_pose_residue( *template_pose, template_pose->total_residue() );
					template_pose->append_residue_by_jump( pose.residue(i), template_pose->total_residue(), "", "", true );
					add_by_jump = true;
				} else {
					template_pose->append_residue_by_bond( pose.residue(i), false );
					add_by_jump = (pose.residue(i).is_upper_terminus());
				}

				pdb_numbering.push_back( i );
				pdb_chains.push_back( 'A' );
			} else {
				add_by_jump = true;
			}
		}
		if (template_pose->total_residue() > 0
				 && !template_pose->residue(template_pose->total_residue()).is_upper_terminus()
				 && template_pose->residue(template_pose->total_residue()).is_polymer())
		core::pose::add_upper_terminus_type_to_pose_residue( *template_pose, template_pose->total_residue() );
		core::pose::PDBInfoOP new_pdb_info = new core::pose::PDBInfo( *template_pose );

		// pdbinfo
		new_pdb_info->set_numbering( pdb_numbering );
		new_pdb_info->set_chains( pdb_chains );
		template_pose->pdb_info( new_pdb_info );
		template_pose->pdb_info()->obsolete( false );

		// pose must be ideal going into hybrid
		//  only the foldtree+sequence is used from input pose
		for (Size i=1; i<=pose.total_residue(); ++i) {
			core::conformation::idealize_position(i, pose.conformation());
		}


		protocols::hybridization::HybridizeProtocol rebuild;
		rebuild.add_template( template_pose, "AUTO", symm_def_file_);
		if (fragments_big_trim_) rebuild.add_big_fragments( fragments_big_trim_ );
		if (fragments_small_trim_) rebuild.add_small_fragments( fragments_small_trim_ );
		rebuild.set_stage1_scorefxn( cen1_scorefxn_ );
		rebuild.set_stage2_scorefxn( cen2_scorefxn_ );
		rebuild.set_stage1_increase_cycles( threaded ? 1.0*censcale_ : 0.0 );
		rebuild.set_stage2_increase_cycles( 0.5*censcale_ );
		rebuild.set_batch_relax( 0 ); // centroid only
		rebuild.apply( pose );
	}

	if (relax_cycles_ > 0) {
		to_fullatom.apply(pose);

		// setup disulfides
		if (disulfs_.size() > 0) {
			utility::vector1< std::pair<Size,Size> > disulfides;
			core::Size ndisulf = disulfs_.size();
			for (core::Size i=1; i<=ndisulf; ++i) {
				utility::vector1<std::string> pair_i = utility::string_split( disulfs_[i], ':');
				runtime_assert( pair_i.size() == 2 );
				core::Size lres=parse_res( pose, pair_i[1] );
				core::Size ures=parse_res( pose, pair_i[2] );
				disulfides.push_back(std::make_pair(lres,ures));
			}
			pose.conformation().fix_disulfides( disulfides );
		} else {
			pose.conformation().detect_disulfides();
		}

		// apply fullatom constraints
		//setup_fullatom_constraints( pose, templates_, template_weights_, "AUTO", "NONE" );

		// relax with flexible angles & jumps
		protocols::relax::RelaxProtocolBaseOP relax_prot = new protocols::relax::FastRelax( fa_scorefxn_, relax_cycles_ );
		relax_prot->set_current_tag( get_current_tag() );
		core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
		mm->set_bb( true ); mm->set_chi( true ); mm->set_jump( true );
		mm->set( core::id::THETA, true );
		relax_prot->set_movemap( mm );
		relax_prot->min_type("lbfgs_armijo_nonmonotone");
		relax_prot->max_iter( relax_max_iter_ );
		relax_prot->apply( pose );
	}
}


//
// repack any missing sidechains after threading (threading does not do this)
void MRMover::pack_missing_sidechains( Pose & pose ) {
	utility::vector1< bool > needToRepack( pose.total_residue() , false );
	bool needToRepackAny = false;  // unused ~Labonte
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		// if there is missing density in the sidechain, then we need to repack
		// check:
		//    (a) CA-CB distance > 3A
		//    (b) CB-any distance > 10A
		if ( pose.residue_type(i).is_protein() && pose.residue_type(i).has("CA") ) {
			numeric::xyzVector< core::Real> ca_pos = pose.residue(i).atom("CA").xyz();
			numeric::xyzVector< core::Real> cb_pos = ca_pos;

			if ( !pose.residue_type(i).has("CB") ) continue;

			cb_pos = pose.residue(i).atom("CB").xyz();
			if ((ca_pos - cb_pos).length() > 3) {
				needToRepack[i] = true;
				needToRepackAny = true;  // unused ~Labonte
			} else {
				for (int j=(int)pose.residue(i).first_sidechain_atom()+1; j<=(int)pose.residue(i).natoms(); ++j) {
					if ( (cb_pos - pose.residue(i).atom(j).xyz()).length() > 10 ) {
						needToRepack[i] = true;
						needToRepackAny = true;
						break;
					}
				}
			}
		}
	}

	if (!needToRepackAny) return; //?? fpd: it's used right here

	core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose );
	taskstd->restrict_to_repacking();
	taskstd->or_include_current(false);
	core::pose::symmetry::make_residue_mask_symmetric( pose, needToRepack );
	taskstd->restrict_to_residues(needToRepack);

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		protocols::simple_moves::symmetry::SymPackRotamersMover pack1( fa_scorefxn_, taskstd );
		pack1.apply( pose );
	} else {
		protocols::simple_moves::PackRotamersMover pack1( fa_scorefxn_, taskstd );
		pack1.apply( pose );
	}
}


//
// trim long gaps from input pose
void MRMover::trim_target_pose( Pose & query_pose, protocols::loops::Loops &loops , core::Size max_gaplength ) {
	using namespace protocols::comparative_modeling;
	using namespace core::sequence;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< bool > to_trim(query_pose.total_residue(), false) ;
	protocols::loops::Loops new_loops;

	//
	for (int i=1; i<=(int)loops.size(); ++i) {
		if (loops[i].size() > max_gaplength ) {
			for (int j=(int)loops[i].start()+1; j<= (int)loops[i].stop()-1; ++j) {
				to_trim[j] = true;
			}
			//fpd check bb connectivity
			bool trim_start, trim_stop;
			trim_stop  = loops[i].stop() == query_pose.total_residue() || query_pose.fold_tree().is_cutpoint(loops[i].stop());
			trim_start = loops[i].start() == 1 || query_pose.fold_tree().is_cutpoint(loops[i].start()-1);
			if (!trim_start) {
				numeric::xyzVector< core::Real > x0 = query_pose.residue( loops[i].start() ).atom( "N" ).xyz();
				numeric::xyzVector< core::Real > x1 = query_pose.residue( loops[i].start()-1 ).atom( "C" ).xyz();
				if ( (x0-x1).length() > 4 )
					trim_start = true;
			}
			if (!trim_stop) {
				numeric::xyzVector< core::Real > x0 = query_pose.residue( loops[i].stop() ).atom( "C" ).xyz();
				numeric::xyzVector< core::Real > x1 = query_pose.residue( loops[i].stop()+1 ).atom( "N" ).xyz();
				if ( (x0-x1).length() > 4 )
					trim_stop = true;
			}
			to_trim[loops[i].start()] = trim_start;
			to_trim[loops[i].stop()] = trim_stop;
		} else {
			new_loops.push_back( loops[i] );
		}
	}

	// # residues in new pose
	Size new_nres = 0;
	for ( Size i = 1; i <= query_pose.total_residue(); ++i ) {
		if (!to_trim[i]) new_nres++;
	}

	if (new_nres == query_pose.total_residue()) return;  // nothing to do
	if (new_nres < 2) {
		std::cerr << "Error: not enough aligned residues! trying to continue" << std::endl;
		return;
	}

	// build new pose, generate seq mapping
	core::Size old_root = query_pose.fold_tree().root();
	core::pose::Pose new_query_pose;

	Size out_ctr=1;
	bool add_by_jump = true;
	core::id::SequenceMapping new_mapping(new_nres, query_pose.total_residue());
	core::id::SequenceMapping new_invmapping(query_pose.total_residue() , new_nres);

	for ( Size i = 1; i <= query_pose.total_residue(); ++i ) {
		if (!to_trim[i]) {
			if (add_by_jump) {
				if (new_query_pose.total_residue() > 0
							 && new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
					core::pose::add_upper_terminus_type_to_pose_residue( new_query_pose, new_query_pose.total_residue() );

				new_query_pose.append_residue_by_jump( query_pose.residue(i), new_query_pose.total_residue(), "", "", true );
				add_by_jump = !query_pose.residue(i).is_polymer();

				if (new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
					core::pose::add_lower_terminus_type_to_pose_residue( new_query_pose, new_query_pose.total_residue() );

			} else if ( !query_pose.residue(i).is_polymer() ) {
				if (new_query_pose.total_residue() > 0
							 && new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
					core::pose::add_upper_terminus_type_to_pose_residue( new_query_pose, new_query_pose.total_residue() );
				new_query_pose.append_residue_by_jump( query_pose.residue(i), new_query_pose.total_residue(), "", "", true );
				add_by_jump = true;
			} else {
				new_query_pose.append_residue_by_bond( query_pose.residue(i), false );
			}

			new_mapping[out_ctr] = i;
			new_invmapping[i] = out_ctr;
			out_ctr++;
		} else {
			add_by_jump = true;
		}
	}
	if (new_query_pose.total_residue() > 0
		 && new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
		core::pose::add_upper_terminus_type_to_pose_residue( new_query_pose, new_query_pose.total_residue() );

	core::kinematics::FoldTree f = new_query_pose.fold_tree();
	if (new_invmapping[old_root] != 0)
		f.reorder( new_invmapping[old_root] );

	//PDBInfo stuff
	core::pose::PDBInfoOP pdb_info( query_pose.pdb_info() );
	core::pose::PDBInfoOP new_pdb_info( new core::pose::PDBInfo(new_nres) );
	utility::vector1< int > pdb_numbering;
	utility::vector1< char > pdb_chains;

	for ( Size i(1); i <= query_pose.total_residue(); ++i ) {
		if (new_invmapping[i] != 0) {
			pdb_numbering.push_back( i );
			pdb_chains.push_back( 'A' );
		}
	}

	// set pdb-wide information
	new_pdb_info->set_numbering( pdb_numbering );
	new_pdb_info->set_chains( pdb_chains );
	new_query_pose.pdb_info( new_pdb_info );
	new_query_pose.pdb_info()->obsolete( false );

	// if pose is rooted on VRT update jump point
	if ( new_query_pose.residue( f.root() ).aa() == core::chemical::aa_vrt ) {
		// find residue closest to center-of-mass
		numeric::xyzVector< core::Real > massSum(0.0,0.0,0.0);
		int nAtms=0, nres=new_query_pose.total_residue();
		for ( int i=1; i<= nres; ++i ) {
			core::conformation::Residue const & rsd( new_query_pose.residue(i) );
			if (rsd.aa() == core::chemical::aa_vrt) continue;
			for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
				core::conformation::Atom const & atom( rsd.atom(j) );
				massSum += atom.xyz();
				nAtms++;
			}
		}
		massSum /= nAtms;

		// try to avoid putting the vrt too close to termini
		int i_min = 1;
		int r_start = (int)std::floor(   nres/3. );
		int r_end   = (int)std::ceil ( 2.*nres/3. );
		core::Real d_min = 99999, this_d;
		for ( int i=r_start; i<=r_end; ++i ) {
			core::conformation::Residue const & rsd( new_query_pose.residue(i) );
			if (!rsd.is_protein() ) continue;

			core::conformation::Atom const & atom( rsd.atom("CA") );
			this_d = (atom.xyz() - massSum).length();
			if (this_d < d_min) {
				d_min = this_d;
				i_min = i;
			}
		}

		f.renumber_jumps();
		f.slide_jump( 1, new_invmapping[old_root], i_min );
		f.reorder( new_invmapping[old_root] );
	}
	new_query_pose.fold_tree( f );

	// remap loops
	for (int i=1; i<=(int)new_loops.size(); ++i) {
		new_loops[i].set_start( new_invmapping[new_loops[i].start()] );
		new_loops[i].set_stop( new_invmapping[new_loops[i].stop()] );

		//fpd
		if (new_invmapping[new_loops[i].start()] == new_invmapping[new_loops[i].stop()]
				&& new_invmapping[new_loops[i].start()] != 0
				&& new_invmapping[new_loops[i].stop()] != 0 ) {
			if (new_invmapping[new_loops[i].start()] != 1)
				new_loops[i].set_start( new_invmapping[new_loops[i].start()]-1 );
			if (new_invmapping[new_loops[i].stop()] != new_nres)
				new_loops[i].set_stop( new_invmapping[new_loops[i].stop()]+1 );
		}
	}
	loops = new_loops;

	// copy pose
	query_pose = new_query_pose;

	if (fragments_big_) {
		// remap fragments
		core::fragment::FragSetOP new_big_frags( fragments_big_->empty_clone() );

		// iterate over frames, clone if mapped
		for ( core::fragment::ConstFrameIterator f=fragments_big_->begin(); f != fragments_big_->end(); ++f ) {
			core::Size start_res = f->start();
			core::Size end_res = f->end();

			// if any residue is unmapped, remove the frame
			bool keepthis = true;
			for (Size j=start_res; j<=end_res; ++j) keepthis &= ( new_invmapping[ j ] != 0 );
			if ( keepthis ) {
				core::fragment::FrameOP new_f = f->clone_with_frags();
				new_f->align(new_invmapping);
				new_big_frags->add( new_f );
			}
		}
		fragments_big_trim_ = new_big_frags;
	}

	if (fragments_small_) {
		core::fragment::FragSetOP new_small_frags( fragments_small_->empty_clone() );
		for ( core::fragment::ConstFrameIterator f=fragments_small_->begin(); f != fragments_small_->end(); ++f ) {
			core::Size start_res = f->start();
			core::Size end_res = f->end();

			// if any residue is unmapped, remove the frame
			bool keepthis = true;
			for (Size j=start_res; j<=end_res; ++j) keepthis &= ( new_invmapping[ j ] != 0 );
			if ( keepthis ) {
				core::fragment::FrameOP new_f = f->clone_with_frags();
				new_f->align(new_invmapping);
				new_small_frags->add( new_f );
			}
		}
		// update
		fragments_small_trim_ = new_small_frags;
	}

}

}
}

