// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file cenrot_jd2.cc
/// @brief test centroid rot model
/// @author Yuan Liu


// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

//
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/Minimizer.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/Mover.fwd.hh>

//for fragment insertion
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <protocols/simple_moves/GunnCost.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>

//sampling
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/ReplicaExchangeMC.hh>

#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>

//docking
#include <utility/tools/make_vector1.hh>
#include <protocols/docking/util.hh>
#include <protocols/docking/types.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
// #include <protocols/simple_filters/RmsdEvaluator.hh>
// #include <protocols/evaluation/EvaluatorFactory.hh>
// #include <protocols/evaluation/PoseEvaluator.hh>

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

/////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::import_pose;
using namespace core::io::pdb;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core::pose;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::kinematics;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;

using namespace core::pack::task;
using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;

//////////////////////////////////////////////////////////////////
static numeric::random::RandomGenerator RG(62331911);
basic::Tracer TR("pilot.wendao.cenrotjd2");

utility::vector1<core::Size> nrecovery(20);
utility::vector1<core::Size> n_total(20);

void* my_main( void* );

///////////////////////////////////////////////////////////////////
OPT_KEY(Boolean, fit_best_rotamer)
OPT_KEY(Boolean, switch_to_centroid)
OPT_KEY(Boolean, input_cenrot_pdb)
OPT_KEY(Boolean, output_cenrot_intcoord)
OPT_KEY(Boolean, output_bestrot_err)
OPT_KEY(Boolean, output_cenrot_pdb)
OPT_KEY(Boolean, repack_cenrot)
OPT_KEY(Boolean, relax_cenrot)
OPT_KEY(Boolean, frag_repack)
OPT_KEY(Boolean, opt_after_relax)
OPT_KEY(Boolean, docking_cenrot)
OPT_KEY(Boolean, rbrelax_cenrot)
OPT_KEY(Boolean, docking_skip_repack)
OPT_KEY(Boolean, repack_min)
OPT_KEY(Boolean, relax_bb)
OPT_KEY(Boolean, keep_silent_header)
OPT_KEY(Real, repack_vdw_scale)

OPT_KEY(Integer, repack_buried_cutoff)
OPT_KEY(Real, repack_bfactor_cutoff)
OPT_KEY(Integer, repack_ncycle)
OPT_KEY(String, output_cenrot_score)
OPT_KEY(String, output_cenrot_dir)
OPT_KEY(String, output_cenrot_prefix)

OPT_KEY(Boolean, cenrot_canonical)
OPT_KEY(Boolean, canonical_recover)
OPT_KEY(Real, canonical_sc_prob)

OPT_KEY(String, cenrot_score)

OPT_KEY(Real, relax_temp)
OPT_KEY(Integer, relax_step_per_cycle)
OPT_KEY(Integer, relax_cycle_number)

OPT_1GRP_KEY(Boolean, min, cenrot) //add min mover alone
OPT_1GRP_KEY(Boolean, min, debug)
OPT_1GRP_KEY(Boolean, min, cartesian)

int main( int argc, char * argv [] ) {
	NEW_OPT(fit_best_rotamer, "fit the exact centroid to the closest in lib", false);
	NEW_OPT(switch_to_centroid, "switch the fa pdb to the old centroid", false);
	NEW_OPT(output_cenrot_intcoord, "output the internal coordinate, for building lib", false);
	NEW_OPT(output_bestrot_err, "output the distance between the native rot and best fitted one", false);
	NEW_OPT(input_cenrot_pdb, "input centroid pdbs for scoring", false);
	NEW_OPT(output_cenrot_pdb, "output centroid pdbs for building database", false);
	NEW_OPT(output_cenrot_score, "score the centroid pdbs", "cenrot_score.out");
	NEW_OPT(output_cenrot_dir, "dir for output centroid pdbs", ".");
	NEW_OPT(output_cenrot_prefix, "prefix for pdbs", "idealized_");

	NEW_OPT(repack_cenrot, "repack the centroid rotamer model", false);
	NEW_OPT(repack_bfactor_cutoff, "count repack side-chain with Bfactor lower than default 100", 100.0);
	NEW_OPT(repack_buried_cutoff, "count repack side-chain with buried cutoff", 0);
	NEW_OPT(repack_ncycle, "how many times to repack", 1);
	NEW_OPT(cenrot_score, "cenrot score weight file", "test.wts");

	NEW_OPT(cenrot_canonical, "relax the centroid rotamer model, bbg+sc", false);
	NEW_OPT(canonical_recover, "recover_low at the end of canonical mover", false);
	NEW_OPT(canonical_sc_prob, "sampling probability of sidechain", 0.8);

	NEW_OPT(docking_cenrot, "docking", false);
	NEW_OPT(docking_skip_repack, "docking without repack", false);
	NEW_OPT(rbrelax_cenrot, "rbrelax_cenrot", false);
	NEW_OPT(keep_silent_header, "keep old scoring lines in silent file", false);

	NEW_OPT(relax_cenrot, "relax the centroid rotamer model", false);
	NEW_OPT(frag_repack, "frag insertion and repack", false);
	NEW_OPT(opt_after_relax, "opt after relax", false);
	NEW_OPT(repack_min, "repack+min", false);
	NEW_OPT(relax_bb, "repack backbone", false);
	NEW_OPT(repack_vdw_scale, "repack using scaled vdw", 1.0);
	NEW_OPT(relax_temp, "temp", 1.0);
	NEW_OPT(relax_step_per_cycle, "step", 100);
	NEW_OPT(relax_cycle_number, "cycle", 1);

	NEW_OPT(min::cenrot, "add min mover alone", false);
	NEW_OPT(min::debug, "debug derivs?", false);
	NEW_OPT(min::cartesian, "cartesian minimization?", false);

	try {
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );
	}
	catch ( utility::excn::EXCN_Base const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////
class OutputCenrotIntCoord : public Mover
{
public:
	OutputCenrotIntCoord() {}
	virtual std::string get_name() const {return "OutputCenrotIntCoordMover";}

	void apply( Pose & p )
	{
		for ( core::Size i=1; i<= p.total_residue(); ++i )
		{
			Residue const & rsd( p.residue(i) );
			/// for each residue, find out the dof-id of CEN
			id::DOF_ID id_dis(id::AtomID(p.residue(i).atom_index("CEN"), i), id::D);
			id::DOF_ID id_ang(id::AtomID(p.residue(i).atom_index("CEN"), i), id::THETA);
			id::DOF_ID id_dih(id::AtomID(p.residue(i).atom_index("CEN"), i), id::PHI);

			//output internal coordinates of centroids
			//as well as phi/psi angle
			if ( !rsd.is_terminus() ) {
				TR << "CEN-INT: " << rsd.name3() << " " << i << " "
				<< p.dof(id_dis) << " "
				<< p.dof(id_ang) << " "
				<< p.dof(id_dih) << " "
				<< p.psi(i) << " "
				<< p.phi(i) << std::endl;
			}
		}
	}
};

class ClearPoseHeader : public Mover
{
public:
	ClearPoseHeader() {}
	virtual std::string get_name() const {return "ClearPoseHeaderMover";}
	void apply( Pose &pose )
	{
		clearPoseExtraScores(pose);
	}
};

class RepackCenrotMover : public Mover {
private:
	core::scoring::ScoreFunctionOP scfxn_;

public:
	virtual std::string get_name() const {return "RepackCenrotMover";}

	void set_scfxn(core::scoring::ScoreFunctionOP const &in_score)
	{
		scfxn_ = in_score;
	}

	void apply( Pose & p ) {
		using namespace core::pack::task;
		TaskFactoryOP main_task_factory = new TaskFactory;
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );

		protocols::simple_moves::PackRotamersMover packrotamersmover;
		packrotamersmover.task_factory(main_task_factory);
		packrotamersmover.score_function(scfxn_);

		packrotamersmover.apply(p);
	}
};

//////////////////////////////////////////////////////////
//bbg+sc
class CenRotSidechainMover : public Mover
{
private:
	core::pack::dunbrack::RotamerLibrary const & rotamer_library_;

	CentroidRotamerSampleData const get_random_rotamer( Residue const &mobile_res )
	{
		ResidueType const &residue_type( mobile_res.type() );

		SingleResidueRotamerLibraryCAP residue_rotamer_library(
			rotamer_library_.get_rsd_library(residue_type)
		);

		SingleResidueCenrotLibraryCAP residue_cenrot_library(
			dynamic_cast< SingleResidueCenrotLibrary const * >(residue_rotamer_library.get())
		);

		assert( residue_cenrot_library );

		utility::vector1< CentroidRotamerSampleData > samples(
			residue_cenrot_library->get_rotamer_samples(mobile_res)
		);

		//randomly select a rotamer, so last_proposal_density would be 1
		Size nrot = samples.size();
		assert( nrot>0 );
		if (nrot==1) {
			return samples[1];
		}
		else {
			return samples[RG.random_range(1,nrot)];
		}
	}

public:
	CenRotSidechainMover(): rotamer_library_( RotamerLibrary::get_instance() ) {}

	virtual std::string get_name() const {return "CenRotSidechainMover";}

	virtual core::Real last_proposal_density_ratio() {return 1.0;}

	void apply( Pose &pose )
	{
		//radomly select a residue
		Size const nres( pose.n_residue() );
		Size mobile_index;

		//pick random res other than GLY/ALA
		do {
			mobile_index = RG.random_range(1, nres);
		} while(pose.residue(mobile_index).aa()==aa_gly || pose.residue(mobile_index).aa()==aa_ala);

		//get random sample
		CentroidRotamerSampleData const sampledata(get_random_rotamer(pose.residue(mobile_index)));

		//assign
		id::DOF_ID id_dis(id::AtomID(pose.residue(mobile_index).atom_index("CEN"), mobile_index), id::D);
		id::DOF_ID id_ang(id::AtomID(pose.residue(mobile_index).atom_index("CEN"), mobile_index), id::THETA);
		id::DOF_ID id_dih(id::AtomID(pose.residue(mobile_index).atom_index("CEN"), mobile_index), id::PHI);
		pose.set_dof(id_dis, sampledata.distance());
		pose.set_dof(id_ang, sampledata.angle());
		pose.set_dof(id_dih, sampledata.dihedral());
	}
};

//////////////////////////////////////////////////////////
//bbg+sc
class CenRotCanonicalMover : public Mover
{
private:
	Size mc_steps_;
	Real mc_temp_;
	Real sc_prob_;

	bool first_run_;

	core::scoring::ScoreFunctionOP scorefxn_;
	MonteCarloOP mc_;
	BBG8T3AMoverOP bbgmover_;

	typedef utility::pointer::owning_ptr< CenRotSidechainMover > CenRotSidechainMoverOP;
	CenRotSidechainMoverOP sidechainmover_;

public:
	CenRotCanonicalMover():first_run_(true)
	{
		using namespace core::pack::task;

		scorefxn_ = core::scoring::getScoreFunction();

		//bbgmover
		bbgmover_ = new BBG8T3AMover();
		//scmove
		sidechainmover_ = new CenRotSidechainMover();

		mc_steps_ = option[relax_step_per_cycle];

		sc_prob_ = option[canonical_sc_prob];
	}

	void apply( Pose &pose )
	{
		mc_temp_ = option[relax_temp];

		if (first_run_) {
			mc_ = new MonteCarlo(pose, *scorefxn_, mc_temp_);
			first_run_ = false;
		}
		else {
			mc_->reset(pose);
		}

		for (Size i=1; i<=mc_steps_; i++) {
			std::string move_type("fake");
			Real proposal_density_ratio=1.0;
			core::Real prob = RG.uniform();

			if ( prob > sc_prob_ ) { //bbg
				bbgmover_->apply(pose);

				//fill ss info
				//without fragment_insertion
				//no ss is assigned if we want to use that for scoring
				core::scoring::dssp::Dssp dssp( pose );
				dssp.insert_ss_into_pose( pose );

				move_type = bbgmover_->type();
				proposal_density_ratio = bbgmover_->last_proposal_density_ratio();
				mc_->boltzmann(pose, move_type, proposal_density_ratio);
			}
			else {
				sidechainmover_->apply(pose);
				move_type = sidechainmover_->type();
				proposal_density_ratio = sidechainmover_->last_proposal_density_ratio();
				mc_->boltzmann(pose, move_type, proposal_density_ratio);
			}
		}

		if (option[canonical_recover]) mc_->recover_low(pose);
		mc_->show_counters();
	}

	virtual std::string get_name() const {return "CenRotCanonicalMover";}
};

//////////////////////////////////////////////////////////
//Rigid body relax (repack and min)
class CenRotRBRelaxMover : public Mover
{
private:
	core::scoring::ScoreFunctionOP scorefxn_dock_;
	core::scoring::ScoreFunctionOP scorefxn_repack_;

	Real repack_score_scale_;

public:
	CenRotRBRelaxMover()
	{
		//setup score functions and movers
		//scorefxn_dock_ = core::scoring::getScoreFunction();
		scorefxn_dock_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_cartmin");
		scorefxn_repack_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_repack");

		repack_score_scale_ = option[repack_vdw_scale];
		Real pack_vdw_weight = scorefxn_repack_->get_weight(core::scoring::vdw);
		if (pack_vdw_weight>0) {
			scorefxn_repack_->set_weight(core::scoring::vdw,pack_vdw_weight*repack_score_scale_);
		}
	}

	void apply( Pose &pose )
	{
		init(pose);

		//packer
		RepackCenrotMover repack;
		repack.set_scfxn(scorefxn_repack_);
		repack.apply(pose);

		//rigid minimizer
		kinematics::MoveMap cstmm;
		cstmm.set_bb(false);
		cstmm.set_chi(false);
		cstmm.set_jump(true);
		core::optimization::AtomTreeMinimizer cstmin;
		core::optimization::MinimizerOptions cstoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, false, false );
		cstoptions.nblist_auto_update(true);
		cstoptions.max_iter(option[run::max_min_iter]);

		core::Real max_vdw = scorefxn_dock_->get_weight( core::scoring::vdw );
		int ncyc = option[ relax::default_repeats ]();

		for (int i=1; i<=ncyc; ++i) {
			scorefxn_dock_->set_weight( core::scoring::vdw, max_vdw/(pow(2.0,Real(ncyc-i))) );
			scorefxn_repack_->set_weight( core::scoring::vdw, repack_score_scale_*max_vdw/(pow(2.0,Real(ncyc-i))) );

			(*scorefxn_repack_)(pose);
			repack.apply(pose);
			cstmin.run( pose, cstmm, *scorefxn_dock_, cstoptions );
			scorefxn_dock_->show( TR.Debug, pose  );
		}
	}

	void init( Pose &pose ) {
		using namespace protocols::docking;

		TR << "Initializing the pose:" << std::endl;

		//foldtree
		//default only handle the first jump
		DockJumps movable_jumps(utility::tools::make_vector1<core::SSize>(1));
		std::string partners="_";
		setup_foldtree(pose, partners, movable_jumps);
		TR << "new fold tree: " << pose.fold_tree();

		//fill dssp info
		//ss info will change during docking
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
	}

	virtual std::string get_name() const {return "CenRotRBRelaxMover";}
};

//////////////////////////////////////////////////////////
//Rigid body perturb and repack
class CenRotDockingMover : public Mover
{
private:
	core::scoring::ScoreFunctionOP scorefxn_dock_;
	core::scoring::ScoreFunctionOP scorefxn_repack_;

	moves::MonteCarloOP mc_;

	simple_moves::PackRotamersMoverOP pack_rotamers_;
	protocols::rigid::RigidBodyPerturbNoCenterMoverOP rb_mover_;
	moves::TrialMoverOP trial_;
	moves::SequenceMoverOP combo_;


	Size inner_cycles_;
	Size outer_cycles_;

	Real temperature_;
	bool first_run_;

	Real rot_magnitude_;
	Real trans_magnitude_;

	bool do_repack_;

public:
	CenRotDockingMover():first_run_(true)
	{
		//setup score functions and movers
		scorefxn_dock_ = core::scoring::getScoreFunction();
		scorefxn_repack_ = core::scoring::getScoreFunction();

		//repack using high vdw can get better recovery power
		Real pack_vdw_weight = scorefxn_repack_->get_weight(core::scoring::vdw);
		if (pack_vdw_weight>0) {
			scorefxn_repack_->set_weight(core::scoring::vdw,pack_vdw_weight*option[repack_vdw_scale]);
		}

		//repack
		pack_rotamers_ = new protocols::simple_moves::PackRotamersMover();
		TaskFactoryOP main_task_factory = new TaskFactory;
		main_task_factory->push_back( new operation::RestrictToRepacking );
		pack_rotamers_->task_factory(main_task_factory);
		pack_rotamers_->score_function(scorefxn_repack_);

		inner_cycles_ = option[relax_step_per_cycle];
		outer_cycles_ = option[relax_cycle_number];

		rot_magnitude_ = 8.0;
		trans_magnitude_ = 3.0;

		do_repack_ = true;
	}

	void skip_repack() {do_repack_=false;}

	void slide_into_contact( Pose &pose ) {
		using namespace moves;

		rigid::RigidBodyTransMover mover( pose, 1 );
		( *scorefxn_dock_ )( pose );

		TR << "Moving away" << std::endl;
		core::Size const counter_breakpoint( 500 );
		core::Size counter( 0 );
		// first try moving away from each other
		//std::cout << "Init vdw: " << pose.energies().total_energies()[ scoring::vdw ] << std::endl;
		// cenrot has bigger intra-chain vdw value
		Real old_vdw = 9999.0;
		Real curr_vdw = pose.energies().total_energies()[ scoring::vdw ];
		while ( curr_vdw < old_vdw-0.0001 && counter <= counter_breakpoint ) {
			mover.apply( pose );
			( *scorefxn_dock_ )( pose );
			++counter;

			old_vdw = curr_vdw;
			curr_vdw = pose.energies().total_energies()[ scoring::vdw ];
		}
		if( counter > counter_breakpoint ){
			TR<<"failed moving away with original vector. Aborting DockingSlideIntoContact."<<std::endl;
			set_current_tag( "fail" );
			return;
		}
		counter = 0;
		// then try moving towards each other
		TR << "Moving together" << std::endl;
		mover.trans_axis().negate();
		while ( counter <= counter_breakpoint && pose.energies().total_energies()[ scoring::vdw ] < 0.1 ) {
			mover.apply( pose );
			( *scorefxn_dock_ )( pose );
			++counter;
		}
		if( counter > counter_breakpoint ){
			TR<<"moving together failed. Aborting DockingSlideIntoContact."<<std::endl;
			set_current_tag( "fail" );
			return;
		}
		// move away again until just touching
		mover.trans_axis().negate();
		mover.apply( pose );
	}

	void init( Pose &pose ) {
		using namespace protocols::docking;

		TR << "Initializing the pose:" << std::endl;

		//foldtree
		TR << "##############################" << std::endl;
		TR << "1. Setting up docking foldtree" << std::endl;
		TR << "old fold tree: " << pose.fold_tree();
		//default only handle the first jump
		DockJumps movable_jumps(utility::tools::make_vector1<core::SSize>(1));
		std::string partners="_";
		setup_foldtree(pose, partners, movable_jumps);
		TR << "new fold tree: " << pose.fold_tree();

		//perturb
		TR << "##############################" << std::endl;
		TR << "2. Sliding into contact" << std::endl;
		slide_into_contact(pose);

		//fill dssp info
		//ss info will change during docking
		TR << "##############################" << std::endl;
		TR << "3. Setup SS info" << std::endl;
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
	}

	void setup( Pose &pose )
	{
		// mover
		core::kinematics::MoveMap mm; mm.set_jump(true);
		rb_mover_ = new rigid::RigidBodyPerturbNoCenterMover( pose, mm, rot_magnitude_, trans_magnitude_, protocols::rigid::n2c );
		combo_ = new moves::SequenceMover();
		combo_->add_mover(rb_mover_);

		if (do_repack_) combo_->add_mover(pack_rotamers_);

		// Monte Carlo
		mc_ = new moves::MonteCarlo( pose, (*scorefxn_dock_), temperature_ );

		trial_ = new moves::TrialMover(combo_, mc_);

		first_run_ = false;
	}

	virtual std::string get_name() const {return "CenRotDockingMover";}

	void apply( Pose &pose ) {
		// only setup once for jd2
		temperature_ = option[relax_temp]; // init temp
		if (first_run_) setup(pose);

		// setup foldtree, gen init conforamtion (perturb, contact)
		init(pose);

		if (!first_run_) mc_->reset(pose); // for multiple runs

		for (Size i=1; i<=outer_cycles_; i++) {
			moves::RepeatMover( trial_, inner_cycles_ ).apply(pose);

			Real accept_rate = trial_->acceptance_rate();
			if ( accept_rate<0.3 ) {
				temperature_ *= 1.1;
				rot_magnitude_ *= 0.9;
				trans_magnitude_ *= 0.9;
				rb_mover_->rot_magnitude( rot_magnitude_ );
				rb_mover_->trans_magnitude( trans_magnitude_ );
			}
			else {
				temperature_  *= 0.9;
				rot_magnitude_ *= 1.1;
				trans_magnitude_ *= 1.1;
				rb_mover_->rot_magnitude( rot_magnitude_ );
				rb_mover_->trans_magnitude( trans_magnitude_ );
			}

			mc_->recover_low(pose); //restore
			mc_->set_temperature(temperature_); //reset temp
		}
	}
};

//////////////////////////////////////////////////////////
//intermediate level sampling, frag-insertion+repack -> MC
class SmoothFragRepackMover : public Mover {
private:
	simple_moves::ClassicFragmentMoverOP sms_;
	core::fragment::FragSetCOP fragset_small_;
	simple_moves::PackRotamersMoverOP pack_rotamers_;
	moves::TrialMoverOP smooth_trial_small_pack_;
	moves::SequenceMoverOP combo_smooth_;

	core::scoring::ScoreFunctionOP scorefxn_;
	moves::MonteCarloOP mc_;

	Size inner_cycles_;
	Size outer_cycles_;

public:
	SmoothFragRepackMover(){
		using namespace core::fragment;
		using simple_moves::SmoothFragmentMover;
		using simple_moves::GunnCost;

		scorefxn_ = core::scoring::getScoreFunction();

		kinematics::MoveMapOP movemap = new kinematics::MoveMap;
		movemap->set_bb( true );

		fragset_small_ = FragmentIO(option[ abinitio::number_3mer_frags ] ).read_data( option[in::file::frag3] );
		sms_ = new SmoothFragmentMover ( fragset_small_, movemap, new GunnCost );

		//repack
		pack_rotamers_ = new protocols::simple_moves::PackRotamersMover();
		TaskFactoryOP main_task_factory = new TaskFactory;
		main_task_factory->push_back( new operation::RestrictToRepacking );
		pack_rotamers_->task_factory(main_task_factory);
		pack_rotamers_->score_function(scorefxn_);

		combo_smooth_ = new moves::SequenceMover();
		combo_smooth_->add_mover(sms_);
		combo_smooth_->add_mover(pack_rotamers_);

		inner_cycles_ = option[relax_step_per_cycle];
		outer_cycles_ = option[relax_cycle_number];
	}

	virtual std::string get_name() const {return "SmoothFragRepackMover";}

	void apply( Pose &p )
	{
		Real temperature = option[relax_temp]; // init temp
		mc_ = new moves::MonteCarlo( p, (*scorefxn_), temperature );
		smooth_trial_small_pack_ = new moves::TrialMover(combo_smooth_, mc_);

		for (Size i=1; i<=outer_cycles_; i++) {
			moves::RepeatMover( smooth_trial_small_pack_, inner_cycles_ ).apply(p);

			Real accept_rate = smooth_trial_small_pack_->acceptance_rate();
			if ( accept_rate<0.3 ) {
				temperature *= 1.2;
			}
			else {
				temperature /= 1.5;
			}

			mc_->recover_low(p); //restore
			mc_->set_temperature(temperature); //reset temp
		}
	}
};


class RescoreCenrot : public Mover {
private:
	core::scoring::ScoreFunctionOP scfxn_;
	pose::PoseOP native_pose_;

public:
	RescoreCenrot(){}
	virtual std::string get_name() const {return "RescoreCenrotMover";}

	void set_scfxn(core::scoring::ScoreFunctionOP const &in_score) {
		scfxn_ = in_score;
	}
	void set_native(pose::PoseOP const &in_pose) {
		native_pose_ = in_pose;
	}

	void apply( Pose & p ) {
		//setup the ss correctly
		core::scoring::dssp::Dssp dssp( p );
		dssp.insert_ss_into_pose( p );

		(*scfxn_)(p);

		if ( option[OptionKeys::evaluation::gdttm] && option[in::file::native].user() ) {
			// protocols::evaluation::MetaPoseEvaluatorOP evaluator = new protocols::evaluation::MetaPoseEvaluator;
			// protocols::evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);
			// evaluator->add_evaluation(
			// 	new protocols::simple_filters::SelectRmsdEvaluator( native_pose_, "_native" )
			// );

			Real gdttm_sc, gdtha_sc;
			core::scoring::CA_gdttm( *native_pose_, p, gdttm_sc, gdtha_sc );
			protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
			job->add_string_real_pair("gdttm", gdttm_sc);
			job->add_string_real_pair("gdtha", gdtha_sc);
		}
	}
};

class RepackMinCenrotMover : public Mover {
public:
	virtual std::string get_name() const {return "RepackMinCenrotMover";}

	void apply( Pose & p ) {
		using namespace core::id;
		using namespace core::pack::task;

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();

		TaskFactoryOP main_task_factory = new TaskFactory;
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );

		protocols::simple_moves::PackRotamersMover packrotamersmover;
		packrotamersmover.task_factory(main_task_factory);
		packrotamersmover.score_function(scorefxn);

		packrotamersmover.apply(p);

		// min option
		bool debug_derivs = option[ OptionKeys::min::debug ]();
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, debug_derivs, debug_derivs );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( option[run::max_min_iter] );

		core::kinematics::MoveMap mm;
		mm.set_bb  ( false );
		mm.set_chi ( false );
		mm.set_jump( false );

		if (option[min::cartesian]) {
			//setup movemap for cartmin: only alows CEN move
			Size const n_res( p.n_residue() );
			for (Size i=1; i<=n_res; i++) {
				core::conformation::Residue const &res_i = p.residue(i);
				//it's a hack, for telling minimizer move only CEN, rather than the whole sidechain
				mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), i ), core::id::RB1 ), true );
			}

			//cart_min sc
			core::optimization::CartesianMinimizer minimizer;
			minimizer.run( p, mm, *scorefxn, minoptions );
		}
		else {
			//setup movemap for atomtree: D,THETA,CHI of CEN
			mm.set_chi(true);
			for ( Size ii = 1; ii <= p.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = p.residue(ii);
				if ( res_i.aa()!=chemical::aa_gly
				     && res_i.aa()!=chemical::aa_ala
				     && res_i.type().has( "CEN")) {
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::D ), true );
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::THETA ), true );
				}
			}

			core::optimization::AtomTreeMinimizer minimizer;
			minimizer.run( p, mm, *scorefxn, minoptions );
		}
	}
};

class MinCenrotMover : public Mover {
public:
	MinCenrotMover(){}
	virtual std::string get_name() const { return "MinCenrotMover"; }

	void apply( core::pose::Pose & pose) {
		using namespace core;
		using namespace core::id;

		bool debug_derivs = option[ OptionKeys::min::debug ]();

		//movemap
		core::kinematics::MoveMap mm;
		if ( option[relax_bb] ) mm.set_bb( true );
		mm.set_chi( false ); mm.set_jump( true );

		// min option
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, debug_derivs, debug_derivs );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( option[run::max_min_iter] );

		// scorefunction1
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_cartmin");

		if (option[min::cartesian]) {
			core::optimization::CartesianMinimizer minimizer;

			//cart boned term, fix CB
			if ( scorefxn->get_weight(core::scoring::cart_bonded)==0 ) {
				scorefxn->set_weight( core::scoring::cart_bonded, 0.1 );
			}

			minimizer.run( pose, mm, *scorefxn, minoptions );
		}
		else {
			mm.set_chi( true );

			core::optimization::AtomTreeMinimizer minimizer;

			//setup movemap for sidechain
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.aa()!=chemical::aa_gly
				     && res_i.aa()!=chemical::aa_ala
				     && res_i.type().has( "CEN")) {
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::D ), true );
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::THETA ), true );
				}
			}

			minimizer.run( pose, mm, *scorefxn, minoptions );
		}

		scorefxn->show( TR.Debug, pose  );
	}
};

class CenRotRelaxMover : public Mover {
private:
	core::scoring::ScoreFunctionOP scorefxn_min_;
	core::scoring::ScoreFunctionOP scorefxn_repack_;

	Real repack_score_scale_;

public:
	CenRotRelaxMover()
	{
		scorefxn_min_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_cartmin");
		scorefxn_repack_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_cenrot_repack");

		repack_score_scale_ = option[repack_vdw_scale];
		Real pack_vdw_weight = scorefxn_repack_->get_weight(core::scoring::vdw);
		if (pack_vdw_weight>0) {
			scorefxn_repack_->set_weight(core::scoring::vdw,pack_vdw_weight*repack_score_scale_);
		}
	}

	virtual std::string get_name() const { return "CenRotRelaxMover"; }

	void apply( core::pose::Pose & pose) {
		using namespace core;
		using namespace core::id;

		bool debug_derivs = option[ OptionKeys::min::debug ]();

		////////////////////////
		//align setup
		////////////////////////

		//movemap
		core::kinematics::MoveMap mm;
		if (option[relax_bb]) mm.set_bb( true );
		mm.set_chi ( true );
		mm.set_jump( true );

		// min option
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, debug_derivs, debug_derivs );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( option[run::max_min_iter] );

		// ramp vdw term
		core::Real max_vdw_min = scorefxn_min_->get_weight( core::scoring::vdw );
		core::Real max_vdw_pack = scorefxn_repack_->get_weight( core::scoring::vdw );
		int ncyc = option[ relax::default_repeats ]();

		//repack in case read from centroid, with wrong sidechain
		RepackCenrotMover repack;
		repack.set_scfxn(scorefxn_repack_);
		repack.apply(pose);

		//min type
		if (option[OptionKeys::min::cartesian]()) {
			core::optimization::CartesianMinimizer minimizer;

			// only for cart_min, set cart_bonded
			// scorefunction0 -- fix bad chainbreaks
			// core::scoring::ScoreFunctionOP scorefxn0 = core::scoring::getScoreFunction();
			// scorefxn0->reset();
			// scorefxn0->set_weight( core::scoring::vdw, 0.1 );
			// scorefxn0->set_weight( core::scoring::cart_bonded, 0.1 );
			// core::scoring::methods::EnergyMethodOptions options0(scorefxn0->energy_method_options());
			// options0.set_cartesian_bonded_linear(true);
			// options0.set_cartesian_bonded_parameters(10,1,0,0,0);
			// scorefxn0->set_energy_method_options(options0);
			// minimizer.run( pose, mm, *scorefxn0, minoptions );

			//make sure cart_bonded is turned on
			if ( scorefxn_min_->get_weight(core::scoring::cart_bonded)==0 ) {
				TR << "Warning: cart_bonded is turned on automatically to 0.1" << std::endl;
				scorefxn_min_->set_weight( core::scoring::cart_bonded, 0.1 );
			}

			TR << "cart_min " << ncyc << std::endl;

			for (int i=1; i<=ncyc; ++i) {
				scorefxn_min_->set_weight( core::scoring::vdw, max_vdw_min/(pow(2.0,Real(ncyc-i))) );
				scorefxn_repack_->set_weight( core::scoring::vdw, max_vdw_pack/(pow(2.0,Real(ncyc-i))) );
				//(*scorefxn_repack_)(pose);
				repack.apply(pose);
				minimizer.run( pose, mm, *scorefxn_min_, minoptions );
				scorefxn_min_->show( TR.Debug, pose  );
			}
		}
		else {
			core::optimization::AtomTreeMinimizer minimizer;

			//setup movemap for sidechain
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.aa()!=chemical::aa_gly
				     && res_i.aa()!=chemical::aa_ala
				     && res_i.type().has( "CEN")) {
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::D ), true );
					mm.set( DOF_ID( AtomID( res_i.atom_index("CEN"), ii ), core::id::THETA ), true );
				}
			}

			TR << "atomtree_min " << ncyc << std::endl;

			for (int i=1; i<=ncyc; ++i) {
				scorefxn_min_->set_weight( core::scoring::vdw, max_vdw_min/(pow(2.0,Real(ncyc-i))) );
				scorefxn_repack_->set_weight( core::scoring::vdw, max_vdw_pack/(pow(2.0,Real(ncyc-i))) );
				//(*scorefxn_repack_)(pose);
				repack.apply(pose);
				minimizer.run( pose, mm, *scorefxn_min_, minoptions );
				scorefxn_min_->show( TR.Debug, pose  );
			}
		}
	}
};


///////////////////////////////////////////////////////////////////
void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves;

	SequenceMoverOP do_cenrot( new SequenceMover() );

	//switch
	//input pose is either centroid or fullatom
	if (option[switch_to_centroid]()) {
		do_cenrot->add_mover(new SwitchResidueTypeSetMover("centroid"));
	}
	else {
		TR << "Switch to CenRot model" << std::endl;
		do_cenrot->add_mover(new SwitchResidueTypeSetMover("centroid_rot"));
	}

	if (!option[keep_silent_header]) do_cenrot->add_mover(new ClearPoseHeader());

	//output intcoord
	if (option[output_cenrot_intcoord]()) {
		do_cenrot->add_mover(new OutputCenrotIntCoord());
	}
	else {
		//setup the scorefxn
		// TODO: if score:weights is not specified, load default cenrot rather than talaris2013
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

		//repack
		if (option[repack_cenrot]) {
			utility::pointer::owning_ptr< RepackCenrotMover > repack = new RepackCenrotMover();
			repack->set_scfxn(score_fxn);
			do_cenrot->add_mover(repack);
		}

		if (option[min::cenrot]) {
			utility::pointer::owning_ptr< MinCenrotMover > mincenrot = new MinCenrotMover();
			do_cenrot->add_mover(mincenrot);
		}

		if (option[repack_min]) {
			utility::pointer::owning_ptr< RepackMinCenrotMover > repackmincenrot = new RepackMinCenrotMover();
			do_cenrot->add_mover(repackmincenrot);
		}

		//relax
		if (option[relax_cenrot]) {
			utility::pointer::owning_ptr< CenRotRelaxMover > relax = new CenRotRelaxMover();
			do_cenrot->add_mover(relax);
		}

		if (option[frag_repack]) {
			utility::pointer::owning_ptr< SmoothFragRepackMover > fragrepack = new SmoothFragRepackMover();
			do_cenrot->add_mover(fragrepack);
		}

		if (option[docking_cenrot]) {
			utility::pointer::owning_ptr< CenRotDockingMover > cenrotdock = new CenRotDockingMover();
			if (option[docking_skip_repack]) cenrotdock->skip_repack();
			do_cenrot->add_mover(cenrotdock);
		}

		if (option[rbrelax_cenrot]) {
			utility::pointer::owning_ptr< CenRotRBRelaxMover > rbrelax = new CenRotRBRelaxMover();
			do_cenrot->add_mover(rbrelax);
		}

		if (option[cenrot_canonical]) {
			utility::pointer::owning_ptr< CenRotCanonicalMover > cenrotcan = new CenRotCanonicalMover();
			do_cenrot->add_mover(cenrotcan);
		}

		//rescore anyway
		utility::pointer::owning_ptr< RescoreCenrot > rescore = new RescoreCenrot();
		rescore->set_scfxn(score_fxn);
		if (option[in::file::native].user()) {
			ResidueTypeSetCAP rsd_set;
			rsd_set=ChemicalManager::get_instance()->residue_type_set( "centroid" );
			core::pose::PoseOP native_pose;
			native_pose = new Pose();
			core::import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
			rescore->set_native(native_pose);
		}
		do_cenrot->add_mover(rescore);
	}

	protocols::jd2::JobDistributor::get_instance()->go( do_cenrot );

	return 0;
}
