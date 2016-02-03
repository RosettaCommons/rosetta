// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/krishna/pH_docking.cc
/// @Use variable protonation states during/after docking

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/pHEnergy.hh>
#include <core/scoring/Energies.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/relax/ClassicRelax.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>
#include <string>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

//neighbors & hbonds
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>

//temp includes krishna
#include <iostream>
#include <fstream>
#include <utility/string_util.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// option key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.krishna.PhDocking" );

class PhDocking : public protocols::moves::Mover {

public:
	PhDocking()
	{
		init_options();
	}

	virtual ~PhDocking(){};

	void
	set_default(){
//		pka_value_ = 1.0;
//		ipka_ = 1.0;
		pH_mode_ = false;
		resfile_ = false;
		unbound_avail_ = false;
	}

	void
	init_options(){
		using namespace basic::options;

		set_default();

		if (option[ OptionKeys::in::file::s ].user()){
			set_file_name(option[ OptionKeys::in::file::s ]()[1]);
		}
		if (option[ OptionKeys::docking::partners ].user()){
			set_partners(option[ OptionKeys::docking::partners ]());
		}
		if (option[ OptionKeys::pH::pH_mode ].user()){
			set_pH_mode(true);
		}
		if (option[ OptionKeys::packing::resfile ].user()){
			set_resfile(true);
		}
		if (option[ OptionKeys::pH::pH_unbound ].user()){
			set_unbound_pdbs(option[ OptionKeys::pH::pH_unbound ]());
		}

	}

	void
	set_file_name(std::string pdb_file_name){
		pdb_file_name_ = pdb_file_name;
	}

	void
	set_partners(std::string partner_info){
		partner_info_ = partner_info;
	}

	void
	set_pH_mode(bool pH_mode){
		pH_mode_ = pH_mode;
	}

	void
	set_resfile(bool resfile){
		resfile_ = resfile;
	}

	void
	set_unbound_pdbs(utility::vector1< utility::file::FileName > unbound_pdbs){
		unbound_avail_ = true;
		unbound_pdbs_ = unbound_pdbs;
	}

	void
	print_hbonds(core::pose::Pose & pose){
		using namespace core;
		using namespace scoring;

		hbonds::HBondDatabaseCOP hb_database(hbonds::HBondDatabase::get_database());
		hbonds::HBondSet hb_set;
		pose.update_residue_neighbors();
		hbonds::fill_hbond_set( pose, false, hb_set, false );
		hb_set.show(pose, std::cout);
	}

	void
	setup_foldtree( core::pose::Pose & pose, std::string partner_chainID)
	{

		using namespace core;
		using namespace protocols;

		core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
		char second_chain = '_';
		core::Size cutpoint = 0;

		using namespace kinematics;
		FoldTree f( pose.fold_tree() );

		for (Size i=1; i<=partner_chainID.length()-1; i++){ //identify second chain from input partner_chainID
			if (partner_chainID[i-1] == '_') second_chain = partner_chainID[i];
		}
		for ( Size i=2; i<= pose.total_residue(); ++i ) {
			if(pdb_info->chain( i ) == second_chain){ //identify cutpoint corresponding to second chain in partner_chainID
				cutpoint = i-1;
				break;
			}
		}

		Size jump_pos1 ( geometry::residue_center_of_mass( pose, 1, cutpoint ) );
		Size jump_pos2 ( geometry::residue_center_of_mass( pose, cutpoint+1, pose.total_residue() ) );

		//setup fold tree based on cutpoints and jump points
		f.clear();
		f.simple_tree( pose.total_residue() );
		f.new_jump( jump_pos1, jump_pos2, cutpoint);
		movable_jumps_.clear();
		movable_jumps_.push_back( 1 );

		Size chain_begin(0), chain_end(0);

		//rebuild jumps between chains N-terminal to the docking cutpoint
		chain_end = cutpoint;
		chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
		while (chain_begin != 1){
			chain_end = chain_begin-1;
			f.new_jump( chain_end, chain_begin, chain_end);
			chain_begin = pose.conformation().chain_begin( pose.chain(chain_end) );
		}

		//rebuild jumps between chains C-terminal to the docking cutpoint
		chain_begin = cutpoint+1;
		chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
		while (chain_end != pose.total_residue()){
			chain_begin = chain_end+1;
			f.new_jump( chain_end, chain_begin, chain_end);
			chain_end = pose.conformation().chain_end( pose.chain(chain_begin) );
		}

		f.reorder( 1 );
		f.check_fold_tree();
		pose.fold_tree( f );
//		f.show( std::cout );

	}

	void
	setup_task_and_repack( core::pose::Pose & pose ){

		using namespace core;
		using namespace scoring;

		core::scoring::ScoreFunctionOP std_scorefxn( ScoreFunctionFactory::create_score_function( "score12" ) );
		pack::task::PackerTaskOP my_task( pack::task::TaskFactory::create_packer_task( pose ));
		my_task->initialize_from_command_line();
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			my_task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		protocols::simple_moves::RotamerTrialsMinMoverOP prepack_mover( new protocols::simple_moves::RotamerTrialsMinMover( std_scorefxn, *my_task ) );
		prepack_mover->apply( pose );

	}


	utility::vector1< core::Real >
	calc_energies( const core::pose::Pose & pose, core::scoring::ScoreFunctionOP scorefxn, utility::vector1_int movable_jumps){

		using namespace core;
		using namespace scoring;

		utility::vector1< core::Real > my_energies;
		Real interaction_energy(0);
		Real binding_energy(0);
		Real unbound_energy(0);

		core::scoring::ScoreFunctionOP docking_scorefxn;
		core::pose::Pose complex_pose = pose;

		docking_scorefxn = new core::scoring::ScoreFunction( *scorefxn ) ;
	    docking_scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.0 );

		// calculate energy of complexed pose
		Real const complex_energy = (*docking_scorefxn)( complex_pose );

		for( utility::vector1_int::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {

			Size const rb_jump = *it;
			core::pose::Pose separated_pose = complex_pose;
			Real trans_magnitude = 1000;
			protocols::rigid::RigidBodyTransMoverOP translate_away ( new protocols::rigid::RigidBodyTransMover( separated_pose, rb_jump ) );
			translate_away->step_size( trans_magnitude );
			translate_away->apply( separated_pose );

			//Before repacking -> Interaction Energy
			Real const separated_energy = (*docking_scorefxn)( separated_pose );
			interaction_energy += (complex_energy - separated_energy);

			//Repack the residues -> Binding Energy
			setup_task_and_repack (separated_pose);
			Real const separated_repack_energy = (*docking_scorefxn)( separated_pose );
			binding_energy += (complex_energy - separated_repack_energy);

		}

		//calculate energies of unbound structures if available
		if (unbound_avail_){
			unbound_energy = ( complex_energy - calc_unbound_score(docking_scorefxn));
		}

		my_energies.push_back( interaction_energy );
		my_energies.push_back( binding_energy );
		my_energies.push_back( unbound_energy );
		return my_energies;
	}


	//Setup if pH_mode is turned on
	void
	setup_pH_mode(core::pose::Pose & pose, core::scoring::ScoreFunctionOP score_fxn ){

		using namespace core;
		using namespace chemical;
		using namespace conformation;
		using namespace scoring;

		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line();

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue( ii ).type().aa() == aa_asp ||
					pose.residue( ii ).type().aa() == aa_glu ||
					pose.residue( ii ).type().aa() == aa_his ||
					pose.residue( ii ).type().aa() == aa_tyr ||
					pose.residue( ii ).type().aa() == aa_lys)
				task->nonconst_residue_task( ii ).restrict_to_repacking();
			else
				task->nonconst_residue_task( ii ).prevent_repacking();
		}

		//std::cout << *task << std::endl;
		core::scoring::ScoreFunctionOP pack_score_fxn( core::scoring::get_score_function() );
		TR << "Tot_E of pose1:" << pdb_file_name_ << "\t" << "\tweighted_scores:\t" << pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << "\ttotal_score:\t" << pose.energies().total_energies().dot( score_fxn->weights() ) << std::endl;

		protocols::simple_moves::RotamerTrialsMinMoverOP pack_mover( new protocols::simple_moves::RotamerTrialsMinMover( pack_score_fxn, *task ) );
		pack_mover->apply(pose);
		(*score_fxn)(pose);

		TR << "Tot_E of pose2:" << pdb_file_name_ << "\t" << "\tweighted_scores:\t" << pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << "\ttotal_score:\t" << pose.energies().total_energies().dot( score_fxn->weights() ) << std::endl;
	}


	//calculate energy of binding defined as U_sc = E_complex - ( E_ligand_unbound + E_receptor_unbound )
	core::Real
	calc_unbound_score(core::scoring::ScoreFunctionOP score_fxn){

		core::Real total_unbound_sc(0);
		core::pose::Pose r_ubd_pose, l_ubd_pose;

		core::import_pose::pose_from_file(r_ubd_pose, unbound_pdbs_[1], core::import_pose::PDB_file);
		core::import_pose::pose_from_file(l_ubd_pose, unbound_pdbs_[2], core::import_pose::PDB_file);

		setup_task_and_repack (r_ubd_pose);
		setup_task_and_repack (l_ubd_pose);

		total_unbound_sc = (*score_fxn)(r_ubd_pose) + (*score_fxn)(l_ubd_pose);
		return total_unbound_sc;
	}


	virtual
	void
	apply( core::pose::Pose & pose ){

		using namespace core;
		using namespace chemical;
		using namespace conformation;
		using namespace scoring;

		setup_foldtree(pose, partner_info_);

		core::Real score12_native(0);
		core::Real score12_mut(0);
		core::Real score12_ddG(0);
		core::scoring::ScoreFunctionOP mut_score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( "score12" ) );
		core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( "docking" ) );
		(*score_fxn)(pose);

		if (pH_mode_){
			setup_pH_mode(pose, score_fxn);
		}

		//setup for ddG calculation for mutants given a resfile
		if (resfile_){
			score12_native = (*mut_score_fxn)(pose);
			pack::task::TaskFactoryOP mut_task_factory = new pack::task::TaskFactory;
			mut_task_factory->push_back( new pack::task::operation::ReadResfile );
			pack::task::PackerTaskOP mut_task = mut_task_factory->create_task_and_apply_taskoperations( pose );
			protocols::simple_moves::RotamerTrialsMinMoverOP mut_pack_mover( new protocols::simple_moves::RotamerTrialsMinMover( mut_score_fxn, *mut_task ) );
			mut_pack_mover->apply(pose);
		}

		//Do a rigid body linmin before feeding the complex to score calculations
		core::kinematics::MoveMapOP movemap ( new core::kinematics::MoveMap );
		movemap->set_bb( false );
		movemap->set_chi( true );
		for( utility::vector1_int::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
			movemap->set_jump(*it, true);
		}
		protocols::simple_moves::MinMoverOP minmover ( new protocols::simple_moves::MinMover(movemap, score_fxn, "linmin", 0.1, false /*use_nblist*/) );
		minmover->apply( pose );

		//calculate ddG for mutants given a resfile
		if (resfile_){
			(*mut_score_fxn)(pose);
			score12_mut = (*mut_score_fxn)(pose);
			score12_ddG = score12_mut - score12_native;
		}

		//Calculate and print out interaction / binding scores
		utility::vector1< core::Real > my_scores = calc_energies( pose, score_fxn, movable_jumps_);
		protocols::jd2::JobOP job = protocols::jd2::JobDistributor::get_instance()->current_job();
		job->add_string_real_pair("I_sc", my_scores[1]);
		job->add_string_real_pair("B_sc", my_scores[2]);
		job->add_string_real_pair("U_sc", my_scores[3]);
		job->add_string_real_pair("ddG", score12_ddG);

	}

	std::string get_name() const { return "PhDocking"; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new PhDocking;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

private:
	std::string pdb_file_name_;
	std::string partner_info_;
	utility::vector1_int movable_jumps_;
	bool pH_mode_;
	bool resfile_;
	bool unbound_avail_;
	utility::vector1< utility::file::FileName > unbound_pdbs_;
};

typedef utility::pointer::owning_ptr< PhDocking > PhDockingOP;

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	devel::init(argc, argv);
	PhDockingOP pH_test(new PhDocking);
	protocols::jd2::JobDistributor::get_instance()->go(pH_test);
}
