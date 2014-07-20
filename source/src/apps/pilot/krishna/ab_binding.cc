// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file src/apps/pilot/krishna/ab_binding.cc
/// @ab-ant binding prediction

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

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

//neighbors, hbonds & sasa calculator
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/sasa.hh>

//temp includes krishna
#include <iostream>
#include <fstream>
#include <core/io/pdb/pose_io.hh>
#include <utility/string_util.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingPrepackProtocol.hh>
#include <protocols/antibody2/LHRepulsiveRamp.hh>

// option key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//#include <apps/pilot/krishna/metrics.cc>

static basic::Tracer TR("apps.pilot.krishna.abbinding");

using namespace core;

class abbinding : public protocols::moves::Mover {

public:
	abbinding()
	{
		init_options();
	}

	virtual ~abbinding(){};

	void
	set_default(){
		refine_ = false;
		cognate_ = false;
		output_raw_scores_ = false;
		frmk_rms_ = 0.0;
		pre_process_ = false;
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
		if (option[ OptionKeys::docking::docking_local_refine ].user()){
			set_refine(option[ OptionKeys::docking::docking_local_refine ]());
		}
		if (option[ OptionKeys::pH::pre_process ].user()){
			set_pre_process(option[ OptionKeys::pH::pre_process ]());
		}
		if (option[ OptionKeys::pH::output_raw_scores ].user()){
			set_output_raw_scores(option[ OptionKeys::pH::output_raw_scores ]());
		}
		if (option[ OptionKeys::pH::cognate_partners ].user()){
			set_cognate_partners(option[ OptionKeys::pH::cognate_partners ]());
		}
		if (option[ OptionKeys::pH::cognate_pdb ].user()){
			set_cognate_pdb(option[ OptionKeys::pH::cognate_pdb ]());
		}

	}

	void
	set_file_name(utility::file::FileName pdb_file_name){
		pdb_file_name_ = pdb_file_name;
	}

	void
	set_partners(std::string partner_info){
		partner_info_ = partner_info;
	}

	void
		set_pre_process(bool pre_process){
		pre_process_ = pre_process;
	}

	void
	set_refine(bool refine){
		refine_ = refine;
	}

	void
	set_output_raw_scores(bool output_raw_scores){
		output_raw_scores_ = output_raw_scores;
	}

	void
	set_cognate_partners(std::string cognate_partner_info){
		cognate_partner_info_ = cognate_partner_info;
	}

	void
	set_cognate_pdb(utility::file::FileName cognate_pdb){
		cognate_pdb_ = cognate_pdb;
		cognate_ = true;
	}

	void
	setup_foldtree(
		core::pose::Pose & pose, std::string partner_chainID
	){

		using namespace core;
		using namespace protocols;

		core::pose::PDB_InfoCOP pdb_info = pose.pdb_info();
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
		f.show( std::cout );

	}

	utility::vector1< core::Real >
	return_raw_scores(
		const core::pose::Pose & pose
	){
		using namespace core;
		using namespace scoring;
		utility::vector1< core::Real > raw_scores;
		raw_scores.push_back(pose.energies().total_energies()[dslf_ca_dih]);
		raw_scores.push_back(pose.energies().total_energies()[dslf_cs_ang]);
		raw_scores.push_back(pose.energies().total_energies()[dslf_ss_dih]);
		raw_scores.push_back(pose.energies().total_energies()[dslf_ss_dst]);
		raw_scores.push_back(pose.energies().total_energies()[fa_atr]);
		raw_scores.push_back(pose.energies().total_energies()[fa_dun]);
		raw_scores.push_back(pose.energies().total_energies()[fa_pair]);
		raw_scores.push_back(pose.energies().total_energies()[fa_rep]);
		raw_scores.push_back(pose.energies().total_energies()[fa_sol]);
		raw_scores.push_back(pose.energies().total_energies()[hack_elec]);
		raw_scores.push_back(pose.energies().total_energies()[hbond_bb_sc]);
		raw_scores.push_back(pose.energies().total_energies()[hbond_lr_bb]);
		raw_scores.push_back(pose.energies().total_energies()[hbond_sc]);
		raw_scores.push_back(pose.energies().total_energies()[hbond_sr_bb]);
		return raw_scores;
	}

	void
	superimpose_antigen(
		core::pose::Pose & pose1,
		core::pose::Pose const & pose2,
		std::string const & ant1_chains,
		std::string const & ant2_chains
	){
		using namespace core;
		using core::id::AtomID;
		using core::id::AtomID_Map;

		std::string atom_name("CA");
		AtomID_Map< AtomID > atom_map;
	    pose::initialize_atomid_map(atom_map, pose1, id::BOGUS_ATOM_ID);

		Size ant1_start_chain = pose::get_chain_id_from_chain(ant1_chains[0], pose1);
		Size ant1_stop_chain = pose::get_chain_id_from_chain(ant1_chains[ant1_chains.length()-1], pose1);
		Size ant2_start_chain = pose::get_chain_id_from_chain(ant2_chains[0], pose2);
		Size ant2_stop_chain = pose::get_chain_id_from_chain(ant2_chains[ant1_chains.length()-1], pose2);

		Size ant1_start_res = pose1.conformation().chain_begin(ant1_start_chain);
		Size ant1_stop_res = pose1.conformation().chain_end(ant1_stop_chain);
		Size ant2_start_res = pose2.conformation().chain_begin(ant2_start_chain);
		Size ant2_stop_res = pose2.conformation().chain_end(ant2_stop_chain);

		Size jj = ant2_start_res;
		for (Size ii = ant1_start_res; ii <= ant1_stop_res; ++ii){
			if (! pose1.residue(ii).has(atom_name)) break;
			AtomID const id1(pose1.residue(ii).atom_index(atom_name), ii);
			AtomID const id2(pose2.residue(jj).atom_index(atom_name), jj);
			atom_map.set(id1, id2);
			jj++;
		}

		core::scoring::superimpose_pose(pose1, pose2, atom_map);
	}

	core::Real
	rmsd_frmk_chothia_num_ab(
		pose::Pose const & pose1,
		pose::Pose const & pose2,
		std::string const & ab1_chains,
		std::string const & ab2_chains
	){
		using namespace core;
		using core::id::AtomID;

		std::string atom_name("CA");
		std::map< core::id::AtomID, core::id::AtomID > atom_id_map;

		Size light_frmk_segments(4);
		Size light_start [] = {1, 35, 57, 98};
		Size light_end [] = {23, 49, 88, 103};

		for(Size i = 0; i < light_frmk_segments; i++){
			for(Size j = light_start[i]; j <= light_end[i]; j++){
				Size ab1_res = pose1.pdb_info()-> pdb2pose(ab1_chains[0],j);
				Size ab2_res = pose2.pdb_info()-> pdb2pose(ab2_chains[0],j);
				if (ab1_res == 0 || ab2_res == 0) continue;
				if (! pose1.residue(ab1_res).has(atom_name)) continue;
				if (! pose2.residue(ab2_res).has(atom_name)) continue;
				/*TR << "Light Chain Alignment: " << j << ab1_chains[0] << " (" << ab1_res << ") -" <<
						j << ab2_chains[0] << " (" << ab2_res << ")"<< std::endl;*/
				AtomID const id1(pose1.residue(ab1_res).atom_index(atom_name), ab1_res);
				AtomID const id2(pose2.residue(ab2_res).atom_index(atom_name), ab2_res);
				atom_id_map[id1] = id2;
			}
		}


		Size heavy_frmk_segments(4);
		Size heavy_start [] = {1, 36, 66, 103};
		Size heavy_end [] = {25, 49, 94, 108};

		for(Size m = 0; m < heavy_frmk_segments; m++){
			for(Size n = heavy_start[m]; n <= heavy_end[m]; n++){
				Size ab1_res = pose1.pdb_info()-> pdb2pose(ab1_chains[1],n);
				Size ab2_res = pose2.pdb_info()-> pdb2pose(ab2_chains[1],n);
				if (ab1_res == 0 || ab2_res == 0) continue;
				if (! pose1.residue(ab1_res).has(atom_name)) continue;
				if (! pose2.residue(ab2_res).has(atom_name)) continue;
				/*TR << "Heavy Chain Alignment: " << n << ab1_chains[1] << " (" << ab1_res << ") -" <<
						n << ab2_chains[1] << " (" << ab2_res << ")"<< std::endl;*/
				AtomID const id1(pose1.residue(ab1_res).atom_index(atom_name), ab1_res);
				AtomID const id2(pose2.residue(ab2_res).atom_index(atom_name), ab2_res);
				atom_id_map[id1] = id2;
			}
		}

		Size natoms(0);
		Real sum(0.0);
		for (std::map< core::id::AtomID, core::id::AtomID >::const_iterator iter = atom_id_map.begin(); iter != atom_id_map.end(); iter++){
			assert (pose1.residue((iter->first).rsd()).atom_name((iter->first).atomno()) ==
					pose2.residue((iter->second).rsd()).atom_name((iter->second).atomno()));

			Vector const & p1(pose1.xyz(iter->first));
			Vector const & p2(pose2.xyz(iter->second));

			sum += (p1-p2).length_squared();
			natoms++;
		}
		return std::sqrt(sum/natoms);
	}

	utility::vector1< std::string >
	partners_from_info(
		std::string const partner_info
	){
		utility::vector1< std::string > partners;   //Vec(partner1,partner2)
		char split = '_';

		for (Size i = 0; i <= partner_info.length()-1; i++){
			if (partner_info[i] == '_'){
				partners.push_back(partner_info.substr(0,i));
				partners.push_back(partner_info.substr(i+1,partner_info.length()-i-1));
				break;
			}
		}
		return partners;
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
		protocols::simple_moves::PackRotamersMoverOP prepack_mover( new protocols::simple_moves::PackRotamersMover( std_scorefxn, my_task ) );
		prepack_mover->apply( pose );
	}

	utility::vector1< core::Real >
	calc_energies_and_sasa(
		const core::pose::Pose & pose,
		core::scoring::ScoreFunctionOP scorefxn,
		utility::vector1_int movable_jumps,
		bool output_raw_scores
	){
		using namespace core;
		using namespace scoring;

		utility::vector1< core::Real > energies_and_sasa;
		utility::vector1< core::Real > raw_interface_energies;
		Real interaction_energy(0);
		Real binding_energy(0);
		Real del_sasa(0);

		core::scoring::ScoreFunctionOP docking_scorefxn;
		core::pose::Pose complex_pose = pose;

		docking_scorefxn = new core::scoring::ScoreFunction( *scorefxn ) ;
		docking_scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.0 );

		// calculate energy of complexed pose
		Real const complex_energy = (*docking_scorefxn)(complex_pose);
		utility::vector1< core::Real > const raw_complex_energies = return_raw_scores(complex_pose);

		Real bound_sasa = core::scoring::calc_total_sasa(complex_pose, 1.4);

		for( utility::vector1_int::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {

			Size const rb_jump = *it;
			core::pose::Pose separated_pose = complex_pose;
			Real trans_magnitude = 1000;
			protocols::rigid::RigidBodyTransMoverOP translate_away ( new protocols::rigid::RigidBodyTransMover( separated_pose, rb_jump ) );
			translate_away->step_size( trans_magnitude );
			translate_away->apply( separated_pose );

			//Calculate del_sasa before repacking
			Real unbound_sasa = core::scoring::calc_total_sasa(separated_pose, 1.4);
			del_sasa = unbound_sasa - bound_sasa;

			//Before repacking -> Interaction Energy
			Real const separated_energy = (*docking_scorefxn)(separated_pose);
			interaction_energy += (complex_energy - separated_energy);

			if (output_raw_scores){
				utility::vector1< core::Real > const raw_separated_energies = return_raw_scores(separated_pose);
				for ( Size n=1; n<= raw_complex_energies.size(); ++n ){
					raw_interface_energies.push_back(raw_complex_energies[n]-raw_separated_energies[n]);
				}
			}

			//Repack the residues -> Binding Energy
			setup_task_and_repack (separated_pose);
			Real const separated_repack_energy = (*docking_scorefxn)( separated_pose );
			binding_energy += (complex_energy - separated_repack_energy);
			//core::io::pdb::dump_pdb( separated_pose, "separated_pose.pdb" );
		}

		energies_and_sasa.push_back(interaction_energy);
		energies_and_sasa.push_back(binding_energy);
		energies_and_sasa.push_back(del_sasa);
		//output raw scores if needed
		if (output_raw_scores){
			for ( Size n=1; n<= raw_complex_energies.size(); ++n ){
				energies_and_sasa.push_back(raw_interface_energies[n]);
			}
		}
		return energies_and_sasa;
	}

	protocols::antibody2::LHRepulsiveRampOP
	setup_RepulsiveRampMover(core::pose::Pose & pose){

		using namespace core;
		using namespace pack::task;
		using namespace pack::task::operation;

		core::scoring::ScoreFunctionCOP refine_score_fxn1(core::scoring::ScoreFunctionFactory::create_score_function("docking", "docking_min"));
		core::scoring::ScoreFunctionCOP refine_score_fxn2(core::scoring::ScoreFunctionFactory::create_score_function("standard"));

		//movemap
		core::kinematics::MoveMapOP movemap = new kinematics::MoveMap();
		movemap->set_chi(false);
		movemap->set_bb(false);
		for( utility::vector1_int::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it){
			movemap->set_jump(*it, true);
		}

		//taskfactory
		core::pack::task::TaskFactoryOP tf = new TaskFactory;
		tf->clear();
		tf->push_back(new OperateOnCertainResidues(new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN")));
		tf->push_back(new InitializeFromCommandline);
		tf->push_back(new IncludeCurrent);
		tf->push_back(new RestrictToRepacking);
		tf->push_back(new NoRepackDisulfides);

	    using namespace protocols::toolbox::task_operations;
		tf->push_back(new RestrictToInterface());

		// unbound_rot
	    pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new pack::rotamer_set::UnboundRotamersOperation();
		unboundrot->initialize_from_command_line();
		operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet(unboundrot);
		tf->push_back(unboundrot_operation);
		core::pack::dunbrack::load_unboundrot(pose);

		protocols::antibody2::LHRepulsiveRampOP RR_mover = new protocols::antibody2::LHRepulsiveRamp(movable_jumps_, refine_score_fxn1, refine_score_fxn2 );
		RR_mover->set_rot_mag(5.0);
		RR_mover->set_task_factory(tf);
		RR_mover->set_move_map(movemap);

		return RR_mover;
	}

	virtual
	void
	apply(
		core::pose::Pose & pose
	){

		using namespace core;
		using namespace chemical;
		using namespace conformation;
		using namespace scoring;
		using namespace protocols;

		setup_foldtree(pose, partner_info_);

		//Dock local refine before feeding the complex to score calculations
		if (pre_process_){
			if (refine_) {
				antibody2::LHRepulsiveRampOP dock_refine_mover = setup_RepulsiveRampMover(pose);
				dock_refine_mover->apply(pose);
			}
			else {
				docking::DockingProtocolOP dock_refine_mover = new docking::DockingProtocol(movable_jumps_, false, false, false);
				dock_refine_mover->apply(pose);
			}
		}

		core::scoring::ScoreFunctionOP score_fxn(core::scoring::ScoreFunctionFactory::create_score_function("docking"));
		//(*score_fxn)(pose);

		//RMSD calculations using cognate Ab framework
		if (cognate_){
			utility::vector1< std::string > cognate_chains = partners_from_info(cognate_partner_info_);
			utility::vector1< std::string > decoy_chains = partners_from_info(partner_info_);
			core::pose::Pose cognate_pose;
			core::import_pose::pose_from_pdb(cognate_pose, cognate_pdb_);
			superimpose_antigen(cognate_pose, pose, cognate_chains[2], decoy_chains[2]);
			frmk_rms_ = rmsd_frmk_chothia_num_ab(cognate_pose, pose, cognate_chains[1], decoy_chains[1]);
		}

/*		//Do a rigid body linmin before feeding the complex to score calculations
		core::kinematics::MoveMapOP movemap (new core::kinematics::MoveMap);
		movemap->set_bb(false);
		movemap->set_chi(true);
		for( utility::vector1_int::const_iterator it = movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
			movemap->set_jump(*it, true);
		}
		protocols::simple_moves::MinMoverOP minmover ( new protocols::simple_moves::MinMover(movemap, score_fxn, "linmin", 0.1, false use_nblist) );
		minmover->apply(pose);*/

		//Calculate and print out interaction / binding scores
		utility::vector1< core::Real > scores_and_sasa = calc_energies_and_sasa( pose, score_fxn, movable_jumps_, output_raw_scores_);
		protocols::jd2::JobOP job = protocols::jd2::JobDistributor::get_instance()->current_job();
		job->add_string_real_pair("I_sc", scores_and_sasa[1]);
		job->add_string_real_pair("B_sc", scores_and_sasa[2]);
		job->add_string_real_pair("dSASA", scores_and_sasa[3]);
		if (cognate_){
			job->add_string_real_pair("Frms", frmk_rms_);
		}
		if (output_raw_scores_){
			job->add_string_real_pair("rint.dslf_ca_dih", scores_and_sasa[4]);
			job->add_string_real_pair("rint.dslf_cs_ang", scores_and_sasa[5]);
			job->add_string_real_pair("rint.dslf_ss_dih", scores_and_sasa[6]);
			job->add_string_real_pair("rint.dslf_ss_dst", scores_and_sasa[7]);
			job->add_string_real_pair("rint.fa_atr", scores_and_sasa[8]);
			job->add_string_real_pair("rint.fa_dun", scores_and_sasa[9]);
			job->add_string_real_pair("rint.fa_pair", scores_and_sasa[10]);
			job->add_string_real_pair("rint.fa_rep", scores_and_sasa[11]);
			job->add_string_real_pair("rint.fa_sol", scores_and_sasa[12]);
			job->add_string_real_pair("rint.hack_elec", scores_and_sasa[13]);
			job->add_string_real_pair("rint.hbond_bb_sc", scores_and_sasa[14]);
			job->add_string_real_pair("rint.hbond_lr_bb", scores_and_sasa[15]);
			job->add_string_real_pair("rint.hbond_sc", scores_and_sasa[16]);
			job->add_string_real_pair("rint.hbond_sr_bb", scores_and_sasa[17]);
		}
	}

	std::string get_name() const { return "abbinding"; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new abbinding;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

private:
	utility::file::FileName pdb_file_name_;
	std::string partner_info_;
	std::string cognate_partner_info_;
	utility::file::FileName cognate_pdb_;
	utility::vector1_int movable_jumps_;
	bool refine_;
	bool pre_process_;
	bool output_raw_scores_;
	bool cognate_;
	core::Real frmk_rms_;
	utility::vector1< utility::file::FileName > unbound_pdbs_;
};

typedef utility::pointer::owning_ptr< abbinding > abbindingOP;

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	devel::init(argc, argv);
	abbindingOP ab_test(new abbinding);
	protocols::jd2::JobDistributor::get_instance()->go(ab_test);
}
