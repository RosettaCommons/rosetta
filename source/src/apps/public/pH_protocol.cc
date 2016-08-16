// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/public/pH_protocol.cc
/// @Calculates the pKa value for a residue

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/pHEnergy.hh>
#include <core/scoring/Energies.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>
#include <string>
#include <utility/vector1.hh>

//neighbors & hbonds
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>

//temp includes krishna
#include <iostream>
#include <fstream>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>


// option key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.krishna.PhProtocol" );

class PhProtocol : public protocols::moves::Mover {

public:
	PhProtocol()
	{
		init_options();
	}

	virtual ~PhProtocol(){};

	void
	set_default(){
		refine_stage_ = false;
		shift_pH_ = 1.0 ;
		pka_value_ = 1.0;
		ipka_ = 1.0;
		neighbor_count_ = 0;
		pka_all_ = false;
		pH_prepack_ = false;
		pH_neighbor_pack_ = false;
		titration_successful_ = false;
		pH_relax_ = false;
		rotamer_prot_stats_ = false;
		method_code_ = "S"; //Default: Site_repack
	}

	void
	init_options(){
		using namespace basic::options;

		set_default();

		if ( option[ OptionKeys::in::file::s ].user() ) {
			set_file_name(option[ OptionKeys::in::file::s ]()[1]);
		}
		if ( option[ OptionKeys::pH::calc_pka::pka_for_resnos].user() ) {
			set_res_no(option[ OptionKeys::pH::calc_pka::pka_for_resnos]());
		}
		if ( option[ OptionKeys::pH::calc_pka::pka_for_chainno].user() ) {
			set_chain_no(option[ OptionKeys::pH::calc_pka::pka_for_chainno]());
		}
		if ( option[ OptionKeys::pH::calc_pka::pka_all].user() ) {
			set_pka_all(true);
		}
		if ( option[ OptionKeys::pH::calc_pka::pH_prepack].user() ) {
			set_pH_prepack(true);
		}
		if ( option[ OptionKeys::pH::calc_pka::pH_relax].user() ) {
			set_pH_relax(true);
		}
		if ( option[ OptionKeys::pH::calc_pka::pH_neighbor_pack].user() ) {
			set_pH_neighbor_pack(true);
		}
		if ( option[ OptionKeys::pH::calc_pka::pka_rad].user() ) {
			set_pka_rad(option[ OptionKeys::pH::calc_pka::pka_rad]());
		}
		if ( option[ OptionKeys::pH::calc_pka::rotamer_prot_stats].user() ) {
			set_rotamer_prot_stats(true);
		}
	}

	void
	set_file_name(std::string pdb_file_name){
		pdb_file_name_ = pdb_file_name;
	}

	void
	set_res_no(utility::vector1< core::Real > pdb_res_nos){
		pdb_res_nos_ = pdb_res_nos;
	}

	void
	set_chain_no(std::string pdb_chain_no){
		pdb_chain_no_ = pdb_chain_no;
	}

	void
	set_pka_all(bool pka_all){
		pka_all_ = pka_all;
	}

	void
	set_pH_prepack(bool pH_prepack){
		pH_prepack_ = pH_prepack;
		method_code_ = "P";
	}

	void
	set_pH_relax(bool pH_relax){
		pH_relax_ = pH_relax;
	}

	void
	set_pH_neighbor_pack(bool pH_neighbor_pack){
		pH_neighbor_pack_ = pH_neighbor_pack;
		method_code_ = "N";
	}

	void
	set_pka_rad(core::Real pka_rad){
		pka_rad_ = pka_rad;
	}

	void
	set_rotamer_prot_stats(bool rotamer_prot_stats){
		rotamer_prot_stats_ = rotamer_prot_stats;
	}

	void
	print_hbonds(core::pose::Pose & pose){
		using namespace core;
		using namespace scoring;

		hbonds::HBondDatabaseCOP hb_database(hbonds::HBondDatabase::get_database());
		hbonds::HBondSet hb_set;
		pose.update_residue_neighbors();
		hbonds::fill_hbond_set( pose, false, hb_set, false );
		hb_set.show(pose, true, std::cout);
	}

	void
	prepack_pose(core::pose::Pose & pose, std::string pdbname ){
		using namespace core;
		using namespace scoring;

		scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
		//  scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) );
		pack::task::PackerTaskOP my_task( pack::task::TaskFactory::create_packer_task( pose ));
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			my_task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		protocols::simple_moves::PackRotamersMoverOP prepack_mover( new protocols::simple_moves::PackRotamersMover( score_fxn, my_task ) );

		(*score_fxn)(pose);
		TR << " SCORE FOR " << pdbname << " BEFORE PREPACK " << (*score_fxn)(pose) << std::endl;

		prepack_mover->apply(pose);

		(*score_fxn)(pose);
		TR << " SCORE FOR " << pdbname << " AFTER PREPACK " << (*score_fxn)(pose) << std::endl;
	}


	void
	relax_pose(core::pose::Pose & pose, std::string pdbname){
		using namespace core;
		using namespace scoring;

		scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
		//  scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) );
		protocols::relax::ClassicRelax my_relax( score_fxn );

		(*score_fxn)(pose);
		TR << " SCORE FOR " << pdbname << " BEFORE RELAX " << (*score_fxn)(pose) << std::endl;

		my_relax.apply(pose);

		(*score_fxn)(pose);
		TR << " SCORE FOR " << pdbname << " AFTER RELAX " << (*score_fxn)(pose) << std::endl;
	}


	//loop over a range of pH values to find the value of pH at which the deprotonated form of a residue
	//will be preferred over the protonated form
	void
	titrate_pH(core::pose::Pose & pose, core::Size res_no, core::Real start_pH, core::Real stop_pH, core::Real step_size){

		for ( core::Real curr_pH = start_pH; curr_pH <= stop_pH; curr_pH += step_size ) {

			using namespace core;
			using namespace chemical;
			using namespace conformation;

			curr_pose_ = pose;
			TR << " CURRENTLY SIMULATING AT pH " << curr_pH << std::endl;
			TR << " CURRENTLY SELECTED RESIDUE " << pose.residue( res_no ).name() << std::endl;
			core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( curr_pose_ ));
			task->initialize_from_command_line();

			if ( !pH_neighbor_pack_ ) {
				for ( Size ii = 1; ii <= curr_pose_.total_residue(); ++ii ) {
					if ( ii == res_no ) {
						task->nonconst_residue_task( ii ).restrict_to_repacking();
					} else {
						task->nonconst_residue_task( ii ).prevent_repacking();
					}
				}
			} else {
				core::conformation::Residue const & rsd1 = curr_pose_.residue( res_no );
				for ( Size ii = 1; ii <= curr_pose_.total_residue(); ++ii ) {
					core::conformation::Residue const & rsd2 = curr_pose_.residue( ii );
					if ( rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom()) ) < pka_rad_ ) {
						task->nonconst_residue_task( ii ).restrict_to_repacking();
					} else {
						task->nonconst_residue_task( ii ).prevent_repacking();
					}
				}
			}

			//   std::cout << *task << std::endl;

			core::scoring::ScoreFunctionOP score_fxn( core::scoring::get_score_function() );
			//   scoring::ScoreFunctionOP score_fxn2( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) );
			core::scoring::methods::pHEnergy::set_pH ( curr_pH );

			//   protocols::simple_moves::RotamerTrialsMinMoverOP pack_mover( new protocols::simple_moves::RotamerTrialsMinMover( score_fxn, *task ) );
			protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( score_fxn, task ) );
			pack_mover->apply(curr_pose_);

			(*score_fxn)(curr_pose_);

			//   TR << "Res_E of pose2:" << pdb_file_name_ << "\t" << curr_pose_.residue( res_no ).name() << "\t" << curr_pose_.pdb_info()->pose2pdb( res_no ) << "\tweighted_scores:\t" << curr_pose_.energies().residue_total_energies( res_no ).weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< curr_pose_.energies().residue_total_energies( res_no ).dot( score_fxn->weights() ) << std::endl;

			if ( (!refine_stage_) && (!rotamer_prot_stats_) ) {

				if ( ( old_pose_.residue( res_no ).name() != curr_pose_.residue( res_no ).name() ) &&
						(!curr_pose_.residue( res_no ).has_variant_type(core::chemical::PROTONATED)) && ( curr_pH != 1.0 ) ) {

					core::conformation::Residue const & res1 = curr_pose_.residue( res_no );
					for ( Size ii = 1; ii <= curr_pose_.total_residue(); ++ii ) {
						core::conformation::Residue const & res2 = curr_pose_.residue( ii );
						/*      if ( res1.xyz( res1.nbr_atom() ).distance( res2.xyz( res2.nbr_atom()) ) <  10 )
						neighbor_count_++;*/
						if ( ii != res_no ) {
							for ( Size jj=1; jj<= res2.nheavyatoms(); ++jj ) {
								//count all heavy atoms within 6A of the nbr_atom as neighbors
								if ( res1.xyz( res1.nbr_atom() ).distance( res2.xyz(jj) ) < 6.0 ) ++neighbor_count_;
							}
						}
					}
					shift_pH_ = curr_pH-1;
					titration_successful_ = true;
					break;
				}
				old_pose_ = curr_pose_;
			}

			if ( refine_stage_ ) {

				if ( (old_pose_.residue( res_no ).name() != curr_pose_.residue( res_no ).name()) &&
						(!curr_pose_.residue( res_no ).has_variant_type(core::chemical::PROTONATED)) ) {

					/*
					TR << "Res_E of pose1:" << pdb_file_name_ << "\t" << old_pose_.residue( res_no ).name() << "\t" << old_pose_.pdb_info()->pose2pdb( res_no ) << "\tweighted_scores:\t" << old_pose_.energies().residue_total_energies( res_no ).weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< old_pose_.energies().residue_total_energies( res_no ).dot( score_fxn->weights() ) << std::endl;
					TR << "Res_E of pose2:" << pdb_file_name_ << "\t" << curr_pose_.residue( res_no ).name() << "\t" << curr_pose_.pdb_info()->pose2pdb( res_no ) << "\tweighted_scores:\t" << curr_pose_.energies().residue_total_energies( res_no ).weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< curr_pose_.energies().residue_total_energies( res_no ).dot( score_fxn->weights() ) << std::endl;
					TR << "Tot_E of pose1:" << pdb_file_name_ << "\t" << old_pose_.residue( res_no ).name() << "\t" << old_pose_.pdb_info()->pose2pdb( res_no ) << "\tweighted_scores:\t" << old_pose_.energies().total_energies().weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< old_pose_.energies().total_energies().dot( score_fxn->weights() ) << std::endl;
					TR << "Tot_E of pose2:" << pdb_file_name_ << "\t" << curr_pose_.residue( res_no ).name() << "\t" << curr_pose_.pdb_info()->pose2pdb( res_no ) << "\tweighted_scores:\t" << curr_pose_.energies().total_energies().weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< curr_pose_.energies().total_energies().dot( score_fxn->weights() ) << std::endl;
					*/

					/*
					TR << "Printing all energies" << std::endl;
					TR << "Old Pose" << std::endl;
					for ( Size ii = 1; ii <= old_pose_.total_residue(); ++ii ) {
					TR << old_pose_.residue(ii).name() << "\t" << old_pose_.pdb_info()->pose2pdb(ii) << "\tweighted_scores:\t" << old_pose_.energies().residue_total_energies(ii).weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< old_pose_.energies().residue_total_energies(ii).dot( score_fxn->weights() ) << std::endl;
					}
					TR << "Current Pose" << std::endl;
					for ( Size ii = 1; ii <= curr_pose_.total_residue(); ++ii ) {
					TR << curr_pose_.residue(ii).name() << "\t" << curr_pose_.pdb_info()->pose2pdb(ii) << "\tweighted_scores:\t" << curr_pose_.energies().residue_total_energies(ii).weighted_string_of( score_fxn->weights() ) <<"\ttotal_score:\t"<< curr_pose_.energies().residue_total_energies(ii).dot( score_fxn->weights() ) << std::endl;
					}
					*/
					/*

					TR << "# Hbond info for\t" << "protonated\t" <<  old_pose_.residue( res_no ).name() << "\t" << old_pose_.pdb_info()->pose2pdb( res_no ) << "\t" << pdb_file_name_ << std::endl;
					print_hbonds(old_pose_);
					TR << "# Hbond info for\t" << "deprotonated\t" <<  curr_pose_.residue( res_no ).name() << "\t" << curr_pose_.pdb_info()->pose2pdb( res_no ) << "\t" << pdb_file_name_ << std::endl;
					print_hbonds(curr_pose_);
					*/

					pka_value_ = curr_pH;
					break;
				}
				old_pose_ = curr_pose_;
			}
		}
	}

	void
	finalize_res_list(core::pose::Pose & pose){
		using namespace core;
		using namespace scoring;
		if ( pka_all_ ) {
			for ( Size i=1; i<=pose.total_residue(); ++i ) {
				std::string res_name3 = pose.residue(i).name3();
				if ( (pose.residue(i).is_protein()) && (res_name3 == "ASP" || res_name3 == "GLU" || res_name3 == "HIS" || res_name3 == "TYR" || res_name3 == "LYS") ) {
					final_res_list_.push_back(i);
				}
			}
		} else {
			for ( Size n=1; n<=pdb_res_nos_.size(); ++n ) {
				Size res_no = pose.pdb_info()->pdb2pose( pdb_chain_no_[0], static_cast< core::SSize >(pdb_res_nos_[n]));
				final_res_list_.push_back(res_no);
			}
		}
	}


	virtual
	void
	apply( core::pose::Pose & pose ){

		using namespace core;
		using namespace chemical;
		using namespace conformation;
		using namespace scoring;

		core::pose::Pose my_pose;
		std::ofstream PKA_LIST;

		finalize_res_list(pose);
		if ( final_res_list_[1] == 0 ) {
			TR << "PLEASE ENTER RESIDUE NOS FOR WHICH PKA IS TO BE CALCULATED AND TRY AGAIN" << std::endl;
			return;
		}

		if ( pH_relax_ ) relax_pose(pose, pdb_file_name_);
		if ( pH_prepack_ ) prepack_pose(pose, pdb_file_name_);

		utility::io::ozstream pka_output("pkas.txt");
		pka_output << "#Methods: S: Site Repack; N: Neighbor Repack; P: Prepack Complete Structure; *: Did Not Titrate;" << std::endl;
		pka_output << "#Fields: Residue N Chain" << "\t" << "IpKa" << "\t" << "pKa" << "\t" << "Method" << "\t" << "File" << std::endl;

		for ( Size n=1; n<= final_res_list_.size(); ++n ) {

			my_pose = pose;

			init_options();

			std::string pdb_res = my_pose.pdb_info()->pose2pdb( final_res_list_[n] );
			std::string res_name = my_pose.residue( final_res_list_[n] ).name();
			std::string res_name3 = my_pose.residue( final_res_list_[n] ).name3();

			//ipka_ details
			switch ( my_pose.residue( final_res_list_[n] ).type().aa() ){
			case aa_asp : ipka_ = 4.0; break;
			case aa_glu : ipka_ = 4.4; break;
			case aa_his : ipka_ = 6.3; break;
			case aa_lys : ipka_ = 10.4; break;
			case aa_tyr : ipka_ = 10.0; break;
			default : ipka_ = 0.0; break;
			}

			TR << " RES BEING PACKED IS " << res_name3 << " " << pdb_res << pdb_file_name_ << std::endl;
			old_pose_ = pose;  //initialize the pose

			//no of neighbors of the residue - determine burial - damn seg faults - no clue why
			//   scoring::TenANeighborGraph const & neighbor_graph ( old_pose_.energies().tenA_neighbor_graph() );
			//   neighbor_count_2 = neighbor_graph.get_node( res_no )->num_neighbors_counting_self();

			titrate_pH(my_pose, final_res_list_[n], 1.0, 14.0, 1.0);
			if ( rotamer_prot_stats_ ) continue;
			while ( !titration_successful_ ) {
				if ( !pH_neighbor_pack_ ) {
					pH_neighbor_pack_ = true;
					method_code_ = "N";
					titrate_pH(my_pose, final_res_list_[n], 1.0, 14.0, 1.0);
				} else {
					pH_neighbor_pack_ = false;
					method_code_ = "P";
					prepack_pose(my_pose, pdb_file_name_);
					titrate_pH(my_pose, final_res_list_[n], 1.0, 14.0, 1.0);
					continue;
				}
			}
			refine_stage_ = true;
			titrate_pH(my_pose, final_res_list_[n], shift_pH_ + 0.1, shift_pH_ + 1.1, 0.1);
			pka_output << res_name3 << " " << pdb_res << "\t" << ipka_ << "\t" << pka_value_ << "\t" << method_code_ << "\t" << pdb_file_name_ << std::endl;
			if ( !titration_successful_ ) {
				pka_output << res_name3 << " " << pdb_res << "\t" << ipka_ << "\t" << ipka_ << "\t" << "*" << "\t" << pdb_file_name_ << std::endl;
			}

			//   PKA_LIST.open("pka_list.txt", std::ios::app);
			//   TR << " THE PREDICTED PKA VALUE FOR " << res_name << pdb_res_nos_ << " FROM CHAIN " << pdb_chain_no_ << " FROM " << pdb_file_name_ << " IS " << pka_value_ << std::endl;
			//   TR << "NEIGHBOR COUNT FOR\t" << res_name3 << "\t" << pdb_chain_no_ << "\t" << pdb_res_nos_ << "\t" << pdb_file_name_ << "\t" << neighbor_count_ << std::endl;
			//   TR << "NEIGHBOR COUNT2 FOR\t" << res_name3 << "\t" << pdb_chain_no_ << "\t" << pdb_res_nos_ << "\t" << pdb_file_name_ << "\t" << neighbor_count_2 << std::endl;
			//   PKA_LIST.close();
		}
	}

	std::string get_name() const { return "PhProtocol"; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new PhProtocol );
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

private:
	core::pose::Pose old_pose_;
	core::pose::Pose curr_pose_;
	bool pka_all_;
	bool refine_stage_;
	bool pH_prepack_;
	bool pH_relax_;
	bool pH_neighbor_pack_;
	bool rotamer_prot_stats_;
	bool titration_successful_;
	core::Real shift_pH_;
	core::Real pka_value_;
	core::Real ipka_;
	core::Real pka_rad_;
	std::string pdb_file_name_;
	utility::vector1< core::Real > pdb_res_nos_;
	utility::vector1< core::Size > final_res_list_;
	std::string pdb_chain_no_;
	std::string method_code_;
	core::Size neighbor_count_;
};

typedef utility::pointer::shared_ptr< PhProtocol > PhProtocolOP;

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{

		devel::init(argc, argv);
		PhProtocolOP pH_test( new PhProtocol );
		protocols::jd2::JobDistributor::get_instance()->go(pH_test);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "Caught Exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
