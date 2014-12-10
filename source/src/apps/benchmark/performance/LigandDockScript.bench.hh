// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/LigandDock.bench.hh
///
/// @brief Dock the ligand in the 7cpa complex.
/// Use all options (flexible ligand, flexible backbone)
/// @author Gordon Lemmon


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//Auto Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/ligand_docking/FinalMinimizer.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/Rotate.hh>
#include <protocols/ligand_docking/SlideTogether.hh>
#include <protocols/ligand_docking/StartFrom.hh>
#include <protocols/ligand_docking/Translate.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>


core::scoring::ScoreFunctionOP
make_scorefxn(std::string weights_tag){
	using namespace core::scoring;

	ScoreFunctionOP sfxn( new ScoreFunction() );
	sfxn->reset();

	// manipulate EnergyMethodOptions here
	methods::EnergyMethodOptions options( sfxn->energy_method_options() );

	options.exclude_protein_protein_fa_elec( true ); // Are we doing this in the parser version?
	sfxn->set_energy_method_options( options );
	sfxn->add_weights_from_file( basic::database::full_name( "scoring/weights/"+weights_tag+".wts" ) );

	if( sfxn->has_zero_weight( fa_intra_rep ) ) sfxn->set_weight( fa_intra_rep, 0.004 ); // from standard.wts

	// For some reason, electrostatics is not in the .wts files...
	// fa_elec has a different dielectric constant than Rosetta++ (10r vs. 6r in ++)
	// It also includes all atom pairs instead of only ligand-protein interactions.
	if( sfxn->has_zero_weight( fa_elec ) ) sfxn->set_weight( fa_elec, 0.42 ); // from Meiler & Baker 2006

	sfxn->set_weight( hbond_sc, 1.30 ); // from Lin Jiang
	sfxn->set_weight( hbond_bb_sc, 1.30 ); // from Lin Jiang
	sfxn->set_weight( rama, 0.2);

	if( sfxn->has_zero_weight( rama ) ) sfxn->set_weight( rama, 0.2 ); // from score12.wts_patch

	return sfxn;

}

using namespace protocols::ligand_docking;

class LigandDockBench : public protocols::moves::Mover {
private:
	StartFrom start_from_;
	Translate_info translate_info_;
	Rotate_info rotate_info_;
	SlideTogether slide_together_;
	HighResDocker high_res_docker_;
	FinalMinimizer final_minimizer_;

public:

	void setup(){
		start_from_.chain("X");
		start_from_.coords(core::Vector(-1.731,32.589,-5.039), "default");
		translate_info_.distribution= get_distribution("uniform");
		translate_info_.angstroms = 5.0;
		translate_info_.cycles = 50;
		rotate_info_.chain = "X";
		rotate_info_.cycles = 200;
		rotate_info_.degrees = 360;
		rotate_info_.distribution = get_distribution("uniform");

		//// Creation of high_res_docker
		LigandAreaOP docking_sidechain( new LigandArea() );
		docking_sidechain->chain_ = 'X';
		docking_sidechain->cutoff_= 6.0;
		docking_sidechain->add_nbr_radius_ = true;
		docking_sidechain->all_atom_mode_ = true;
		docking_sidechain->minimize_ligand_ = 10;

		InterfaceBuilderOP sc_4_docking( new InterfaceBuilder(utility::vector1<LigandAreaOP>(1, docking_sidechain)) );
		MoveMapBuilderOP docking_movemap( new MoveMapBuilder(sc_4_docking, 0, true) );
		high_res_docker_ = HighResDocker(6, 3, std::vector<std::string>(1,"X"), make_scorefxn("ligand_soft_rep") , docking_movemap);

		/////// Creation of final_minimizer
		LigandAreaOP final_sidechain( new LigandArea() );
		final_sidechain->chain_ = 'X';
		final_sidechain->cutoff_= 6.0;
		final_sidechain->add_nbr_radius_ = true;
		final_sidechain->all_atom_mode_ = true;
		//////////////////////////////
		LigandAreaOP final_backbone( new LigandArea() );
		final_backbone->chain_ = 'X';
		final_backbone->cutoff_= 7.0;
		final_backbone->add_nbr_radius_ = false;
		final_backbone->all_atom_mode_ = true;
		final_backbone->Calpha_restraints_ = 0.3;
		/////////////////////////////
		InterfaceBuilderOP sc_4_final( new InterfaceBuilder(utility::vector1<LigandAreaOP>(1, final_sidechain)) );
		InterfaceBuilderOP bb_4_final( new InterfaceBuilder(utility::vector1<LigandAreaOP>(1, final_backbone), 3) );
		MoveMapBuilderOP final_movemap( new MoveMapBuilder(sc_4_final, bb_4_final, true) );
		final_minimizer_ = FinalMinimizer( make_scorefxn("ligand"), final_movemap);
	}

	LigandDockBench():
			protocols::moves::Mover("LigandDockBench"), slide_together_("X")
	{
		//setup(); // have to read in flags first
	}
	virtual ~LigandDockBench(){}
	LigandDockBench(LigandDockBench const & that):
		protocols::moves::Mover( that ),
		translate_info_(that.translate_info_),
		rotate_info_(that.rotate_info_),
		slide_together_(that.slide_together_),
		high_res_docker_(that.high_res_docker_),
		final_minimizer_(that.final_minimizer_)
	{}

	virtual protocols::moves::MoverOP clone() const{
		return protocols::moves::MoverOP( new LigandDockBench( *this ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const{
		return protocols::moves::MoverOP( new LigandDockBench );
	}
	virtual std::string get_name() const{ return "LigandDockBench";}

	virtual void apply( core::pose::Pose & pose ){
		start_from_.apply(pose);
		translate_info_.chain_id = core::pose::get_chain_id_from_chain("X", pose);
		translate_info_.jump_id = core::pose::get_jump_id_from_chain_id(translate_info_.chain_id, pose);
		Translate translate(translate_info_);
		translate.apply(pose);
		rotate_info_.chain_id = core::pose::get_chain_id_from_chain(rotate_info_.chain, pose);
		rotate_info_.jump_id = core::pose::get_jump_id_from_chain_id(rotate_info_.chain_id, pose);
		Rotate rotate(rotate_info_);
		rotate.apply(pose);
		slide_together_.apply(pose);
		high_res_docker_.apply(pose);
		final_minimizer_.apply(pose);
	}

};
typedef utility::pointer::shared_ptr<LigandDockBench> LigandDockBenchOP;


class LigandDockScriptBenchmark : public PerformanceBenchmark
{
public:
	LigandDockScriptBenchmark(std::string name) : PerformanceBenchmark(name) {};

	core::pose::PoseOP ligand_dock_pose_;
	LigandDockBenchOP ligand_dock_protocol_;

	virtual void setUp() {
		basic::options::option.load_options_from_file("ligand_dock/ligand_dock_script_flags.txt");

		ligand_dock_protocol_ = LigandDockBenchOP( new LigandDockBench);
		ligand_dock_protocol_->setup(); // Can't call this from the constructor because flags are needed by score function

		std::string pdb_file_name= basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
		ligand_dock_pose_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*ligand_dock_pose_, pdb_file_name);
	};

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(1 * scaleFactor) );
		if( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scale factor
		for(core::Size i=0; i< reps; i++) {
			ligand_dock_protocol_->apply(*ligand_dock_pose_);
		}

	};

	virtual void tearDown() {};
};
