/// @file
/// @brief


#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <protocols/relax/util.hh>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCD.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>

#include <protocols/relax/RelaxProtocolBase.hh>

#include <protocols/rbsegment_relax/AutoRBRelaxMover.hh>
#include <protocols/rbsegment_relax/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/util.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/constraints/util.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

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
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//
#include <iostream>
#include <string>
#include <sstream>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/comparative_modeling/AlignmentSet.fwd.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>


OPT_1GRP_KEY(Integer, MR, max_gaplength_to_model)
OPT_1GRP_KEY(Integer, MR, nrebuildcycles)
OPT_1GRP_KEY(Real, MR, cen_dens_wt)
OPT_1GRP_KEY(StringVector, MR, disulf)
OPT_1GRP_KEY(Boolean, MR, debug)
OPT_1GRP_KEY(Boolean, MR, smart_foldtree)
OPT_1GRP_KEY(String, MR, mode)


static basic::Tracer TR("rosetta_MR");
///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////

class MRMover : public protocols::moves::Mover {
private:
	core::scoring::ScoreFunctionOP densonly_scorefxn_;
	core::scoring::ScoreFunctionOP fadens_scorefxn_;
	core::scoring::ScoreFunctionOP cen_scorefxn_;
	core::scoring::ScoreFunctionOP cendens_scorefxn_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;

	utility::vector1< core::fragment::FragSetOP > frag_libs_;
	utility::vector1<std::string> disulfs_;
	core::Size nrebuildcycles;

public:
	MRMover(){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		densonly_scorefxn_ = new core::scoring::ScoreFunction();
		fadens_scorefxn_ = core::scoring::getScoreFunction();
		cen_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("cen_std","score4L");
		cendens_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("cen_std","score4L");
		fa_scorefxn_ = core::scoring::getScoreFunction();

		// if not specified should set patterson wt to some value
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *fadens_scorefxn_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *densonly_scorefxn_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *cendens_scorefxn_ );

		if (option[ OptionKeys::MR::disulf ].user()) {
			disulfs_ = option[ OptionKeys::MR::disulf ]();
		}

		// use fastdens
		if (option[ OptionKeys::edensity::mapfile ].user())
			cen_scorefxn_->set_weight( core::scoring::elec_dens_fast, option[ OptionKeys::MR::cen_dens_wt ]() );

		// set cst weight (fa only)
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *cen_scorefxn_  );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fadens_scorefxn_  );

		//nrebuildcycles =  option[ OptionKeys::MR::nrebuildcycles ]();
		if ( option[ OptionKeys::loops::frag_files ].user() )
			protocols::loops::read_loop_fragments( frag_libs_ );
	}

	virtual std::string get_name() const {return "MRMover";}

	void pack_missing_sidechains( Pose & pose ) {
		utility::vector1< bool > needToRepack( pose.total_residue() , false );
		//bool needToRepackAny = false;  // unused ~Labonte
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
					//needToRepackAny = true;  // unused ~Labonte
				} else {
					for (int j=(int)pose.residue(i).first_sidechain_atom()+1; j<=(int)pose.residue(i).natoms(); ++j) {
						if ( (cb_pos - pose.residue(i).atom(j).xyz()).length() > 10 ) {
							needToRepack[i] = true;
							//needToRepackAny = true;  // unused ~Labonte
							break;
						}
					}
				}
			}
		}

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

	void fast_loopclose( Pose &pose, protocols::loops::LoopsOP const loops ) {
		using namespace protocols::loops;

		// for all loops
		core::kinematics::FoldTree f_orig = pose.fold_tree();
		for ( Loops::iterator it=loops->v_begin(), it_end=loops->v_end(); it != it_end; ++it ) {
			Loop buildloop( *it );

			set_single_loop_fold_tree( pose, buildloop );
			//  --> idealize+extend
			set_extended_torsions( pose, buildloop );
			//  --> ccd close
			bool chainbreak_present = !pose.residue(  buildloop.start() ).is_lower_terminus() &&
	                                  !pose.residue(  buildloop.stop() ).is_upper_terminus();
			if ( chainbreak_present ) {
				core::kinematics::MoveMapOP mm_one_loop = new core::kinematics::MoveMap();
				set_move_map_for_centroid_loop( buildloop, *mm_one_loop );
				loop_mover::perturb::fast_ccd_close_loops( pose, buildloop,  *mm_one_loop );
			}

			//  --> restore foldtree
			pose.fold_tree( f_orig );
		}

		// now call the kinematic loop sampler
		protocols::comparative_modeling::LoopRelaxMoverOP lr_mover( new protocols::comparative_modeling::LoopRelaxMover );
		lr_mover->scorefxns( cen_scorefxn_, fadens_scorefxn_ );   // centroid threading doesn't use density score
		lr_mover->loops( loops );
		lr_mover->remodel( "perturb_kic" );
		lr_mover->cmd_line_csts( true );
		lr_mover->rebuild_filter( 999 );
		lr_mover->n_rebuild_tries( 10 );
		lr_mover->copy_sidechains( true );
		lr_mover->set_current_tag( get_current_tag() );
		//lr_mover->apply( pose );
		//remove_cutpoint_variants( pose );

		// to fullatom
		protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");
		to_fullatom.apply(pose);

		// finally, cartesian relax
		core::kinematics::MoveMapOP mm_relax = new core::kinematics::MoveMap();
		mm_relax->set_chi(true);
		mm_relax->set_bb(true);
		for ( Loops::iterator it=loops->v_begin(), it_end=loops->v_end(); it != it_end; ++it )
			for (core::Size i=it->start(); i<=it->stop(); ++i)
				mm_relax->set_bb(i, true);

		protocols::relax::RelaxProtocolBaseOP relax_prot = protocols::relax::generate_relax_from_cmd();
		relax_prot->set_current_tag( get_current_tag() );
		relax_prot->cartesian( true );
		relax_prot->min_type("lbfgs_armijo_nonmonotone");
		relax_prot->set_movemap(mm_relax);

		core::scoring::ScoreFunctionOP cart_fa = new core::scoring::ScoreFunction(*fadens_scorefxn_);
		cart_fa->set_weight( core::scoring::cart_bonded_angle, 0.5 );
		cart_fa->set_weight( core::scoring::cart_bonded_length, 0.5 );
		cart_fa->set_weight( core::scoring::cart_bonded_torsion, 0.5 );
		relax_prot->set_scorefxn( cart_fa );

		relax_prot->apply( pose );
	}


	void trim_target_pose( Pose & query_pose, protocols::loops::Loops &loops , core::Size max_gaplength ) {
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
					       && !new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
						core::pose::add_upper_terminus_type_to_pose_residue( new_query_pose, new_query_pose.total_residue() );

					new_query_pose.append_residue_by_jump( query_pose.residue(i), new_query_pose.total_residue(), "", "", true );
					add_by_jump = !query_pose.residue(i).is_polymer();

					if (!new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
						core::pose::add_lower_terminus_type_to_pose_residue( new_query_pose, new_query_pose.total_residue() );

				} else if ( !query_pose.residue(i).is_polymer() ) {
					if (new_query_pose.total_residue() > 0
					       && !new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
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
			 && !new_query_pose.residue(new_query_pose.total_residue()).is_polymer())
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
			//int resid_num;  // unused ~Labonte
			if (new_invmapping[i] != 0) {
				pdb_numbering.push_back( i );
				pdb_chains.push_back( 'A' );
			}
		}

		// set pdb-wide information
		new_pdb_info->set_numbering( pdb_numbering );
		new_pdb_info->set_chains( pdb_chains );
		new_query_pose.pdb_info( new_pdb_info );

		/////
		//new_query_pose.conformation().fill_missing_atoms( missing );

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
			int r_start = (int)std::floor(   nres/3 );
			int r_end   = (int)std::ceil ( 2*nres/3 );
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

		// remap fragments
		utility::vector1< core::fragment::FragSetOP > new_frag_libs;
		for (int i=1; i<=(int)frag_libs_.size(); ++i) {
			// empty_clone
			core::fragment::FragSetOP new_frag_set( frag_libs_[i]->empty_clone() );

			// iterate over frames, clone if mapped
			for ( core::fragment::ConstFrameIterator f=frag_libs_[i]->begin(); f != frag_libs_[i]->end(); ++f ) {
				core::Size start_res = f->start();
				if ( new_invmapping[ start_res ] != 0 ) {
					core::fragment::FrameOP new_f = f->clone_with_frags();
					new_f->align(new_invmapping);
					new_frag_set->add( new_f );
				}
			}
			new_frag_libs.push_back( new_frag_set );
		}
		frag_libs_ = new_frag_libs;
		query_pose = new_query_pose;
	}

	void apply_autorb( Pose &pose ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::moves;
		using namespace core::pack::task;

		MoverOP autorb( new protocols::rbsegment_relax::AutoRBMover );
		autorb->apply( pose );
	}

	void apply_relax( Pose &pose ) {
		using namespace protocols::loops;
		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMover pre_mover;
			pre_mover.apply( pose );
		} else {
			core::pose::addVirtualResAsRoot( pose );
		}

		// optionally set "smart" foldtree
		// TO DO symmetric version
		if ( option[ MR::smart_foldtree ].user() ) {
			using namespace protocols::rbsegment_relax;
			utility::vector1< RBSegment > rigid_segs, rb_chunks;
			utility::vector1< core::Size > jumps;
			protocols::loops::Loops loops;
			core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();

			guess_rbsegs_from_pose( pose, rigid_segs, rb_chunks, loops );
			jumps = setup_pose_rbsegs_keep_loops( pose,  rigid_segs , loops,  movemap );
		}

		if ( option[ OptionKeys::edensity::mapfile ].user() ) {
			protocols::moves::MoverOP dens_dock( new protocols::electron_density::SetupForDensityScoringMover );
			dens_dock->apply( pose );
		}

		protocols::relax::RelaxProtocolBaseOP relax_prot = protocols::relax::generate_relax_from_cmd();
		relax_prot->set_current_tag( get_current_tag() );
		relax_prot->apply( pose );
	}

	void apply_loopmodel( Pose &pose ) {
		using namespace protocols::loops;
		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		bool debug = option[ OptionKeys::MR::debug ]();

		// set up initial loop build
		core::Size nres = pose.total_residue();
		while (!pose.residue(nres).is_polymer()) nres--;

		// pack all missing SCs
		pack_missing_sidechains( pose );

		// setup for symmetry
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMover pre_mover;
			pre_mover.apply( pose );
		} else {
			core::pose::addVirtualResAsRoot( pose );
		}
		if ( option[ OptionKeys::edensity::mapfile ].user() ) {
			protocols::moves::MoverOP dens_dock( new protocols::electron_density::SetupForDensityScoringMover );
			dens_dock->apply( pose );
		}

		if (debug) pose.dump_pdb( "pre_loopmodel.pdb" );

		// if a loopfile is given, use it
		// otherwise, use model
		protocols::loops::LoopsOP to_rebuild;
		if ( option[ OptionKeys::loops::loop_file ].user() ) {
			to_rebuild = new protocols::loops::Loops( true );
		} else {
			// autoselect
			if (fadens_scorefxn_->get_weight(core::scoring::patterson_cc) > 0) {
				to_rebuild = new Loops( protocols::electron_density::findLoopFromPatterson( pose, 10, (core::Size)std::ceil(pose.total_residue()/40), false ) );
			} else {
				to_rebuild = new Loops( protocols::electron_density::findLoopFromDensity( pose, 0.3, -1, 1 ) );
			}
		}
		to_rebuild->choose_cutpoints( pose );

		// now do loopbuilding
		protocols::comparative_modeling::LoopRelaxMoverOP lr_mover( new protocols::comparative_modeling::LoopRelaxMover );
		lr_mover->loops( to_rebuild );
		lr_mover->scorefxns( cendens_scorefxn_, fadens_scorefxn_ );
		lr_mover->frag_libs( frag_libs_ );
		lr_mover->relax( option[ OptionKeys::loops::relax ]() );
		lr_mover->remodel( option[ OptionKeys::loops::remodel ]() );
		lr_mover->cmd_line_csts( true );
		lr_mover->rebuild_filter( option[ OptionKeys::cm::loop_rebuild_filter ]() );
		lr_mover->n_rebuild_tries( option[ OptionKeys::cm::max_loop_rebuild ]() );
		lr_mover->copy_sidechains( true );
		lr_mover->set_current_tag( get_current_tag() );
		lr_mover->apply(pose);
		remove_cutpoint_variants( pose );

		// score pose
		if (pose.is_fullatom())
			(*fadens_scorefxn_)(pose);
		else
			(*cen_scorefxn_)(pose);
	}

	void apply_threading( Pose &pose, bool superfast=false ) {
		using namespace protocols::loops;
		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		bool debug = option[ OptionKeys::MR::debug ]();

		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
		protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");

		protocols::comparative_modeling::ThreadingJobCOP job = dynamic_cast< protocols::comparative_modeling::ThreadingJob const*  >(
			JobDistributor::get_instance()->current_job()->inner_job().get() );
		if ( !job ) {
			utility_exit_with_message(
				"CORE ERROR: You must use the ThreadingJobInputter with the LoopRelaxThreadingMover "
				"- did you forget the -in:file:template_pdb option?" );
		}

		// set up initial loop build
		core::Size nres = pose.total_residue();
		while (!pose.residue(nres).is_polymer()) nres--;
		LoopsOP my_loops = new Loops( job->loops( nres ) );

		if ( option[ MR::max_gaplength_to_model ].user() ) {
			trim_target_pose( pose, *my_loops, option[ MR::max_gaplength_to_model ]() );
		}

		// find disulfides conserved from template
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

		// pack all missing SCs
		pack_missing_sidechains( pose );

		// setup for symmetry
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMover pre_mover;
			pre_mover.apply( pose );
		} else {
			core::pose::addVirtualResAsRoot( pose );
		}

		if ( option[ OptionKeys::edensity::mapfile ].user() ) {
			protocols::moves::MoverOP dens_dock( new protocols::electron_density::SetupForDensityScoringMover );
			dens_dock->apply( pose );
		}

		if (debug) pose.dump_pdb( "pre_loopmodel.pdb" );

		my_loops->choose_cutpoints( pose );

		if ( superfast ) {  // no loop modeling
			// fix loops initially
			fast_loopclose( pose, my_loops );
		} else if (option[ MR::max_gaplength_to_model ]() == 0) {  // no loop modeling
			// fastrelax only
			apply_relax( pose );
		} else {
		  protocols::comparative_modeling::LoopRelaxMoverOP lr_mover( new protocols::comparative_modeling::LoopRelaxMover );
			lr_mover->scorefxns( cen_scorefxn_, fadens_scorefxn_ );   // centroid threading doesn't use density score
			lr_mover->frag_libs( frag_libs_ );
			lr_mover->loops( my_loops );
			lr_mover->relax( option[ OptionKeys::loops::relax ]() );
			lr_mover->remodel( option[ OptionKeys::loops::remodel ]() );
			lr_mover->cmd_line_csts( true );
			lr_mover->rebuild_filter( option[ OptionKeys::cm::loop_rebuild_filter ]() );
			lr_mover->n_rebuild_tries( option[ OptionKeys::cm::max_loop_rebuild ]() );
			lr_mover->copy_sidechains( true );
			lr_mover->set_current_tag( get_current_tag() );

			lr_mover->apply( pose );
			remove_cutpoint_variants( pose );
		}

		// rescore pose
		if (pose.is_fullatom())
			(*fadens_scorefxn_)(pose);
		else
			(*cen_scorefxn_)(pose);
	}

	void apply_makefrags( Pose &pose ) {
		using namespace protocols::loops;
		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::fragment;
		using namespace utility::file;
		protocols::comparative_modeling::ThreadingJobCOP job = dynamic_cast< protocols::comparative_modeling::ThreadingJob const*  >(
		  JobDistributor::get_instance()->current_job()->inner_job().get() );
		if ( !job ) {
			utility_exit_with_message(
				"CORE ERROR: You must use the ThreadingJobInputter with the LoopRelaxThreadingMover "
				"- did you forget the -in:file:template_pdb option?" );
		}

		// set up initial loop build
		core::Size nres = pose.total_residue();
		while (!pose.residue(nres).is_polymer()) nres--;
		Loops my_loops( job->loops( nres ) );

		if ( !option[ MR::max_gaplength_to_model ].user() ) return;   // nothing to do

		trim_target_pose( pose, my_loops, option[ MR::max_gaplength_to_model ]() );

		// dump remapped fragfiles
		utility::vector1<int> frag_sizes( option[ OptionKeys::loops::frag_sizes ] );
		FileVectorOption frag_files( option[ OptionKeys::loops::frag_files ] );
		runtime_assert( frag_sizes.size() == frag_libs_.size());
		for (int i=1; i<=(int)frag_libs_.size(); ++i) {
			// don't dump 1mers
			if (frag_sizes[i] == 1 || frag_files[i] == std::string("none") ) continue;
			FragmentIO().write_data( option[ OptionKeys::out::file::frag_prefix ]+ ObjexxFCL::string_of(frag_sizes[i],0) , *frag_libs_[i] );
		}
	}

	void apply( Pose & pose ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string mode = option[ OptionKeys::MR::mode ]();
		if ( mode == "cm" ) {
			apply_threading( pose );
		} else if ( mode == "fastcm" ) {
			apply_threading( pose, true );
		} else if ( mode == "cmfrags" ) {
			apply_makefrags( pose ); //
		} else if ( mode == "loopmodel" ) {
			apply_loopmodel( pose ); //
		} else if ( mode == "relax" ){
			apply_relax( pose );
		} else if ( mode == "autorb" ) {
			apply_autorb( pose ); //
		} else {
			TR.Error << "UNKNOWN MODE " << mode << std::endl;
		}
	}

	protocols::moves::MoverOP clone() const { return new MRMover( *this ); }
	protocols::moves::MoverOP fresh_instance() const { return new MRMover; }
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( new MRMover() );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	try {
	// initialize option and random number system
	NEW_OPT(MR::max_gaplength_to_model, "max gaplength to rebuild", true);
    NEW_OPT(MR::mode, "mode", "cm");
    NEW_OPT(MR::cen_dens_wt, "centroid density weight", 4.0);
    NEW_OPT(MR::disulf, "disulf patterning", utility::vector1<std::string>());
    NEW_OPT(MR::smart_foldtree, "mode", false);
	NEW_OPT(MR::debug, "output debug pdbs?", false);
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	}
	catch ( utility::excn::EXCN_Base const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}
