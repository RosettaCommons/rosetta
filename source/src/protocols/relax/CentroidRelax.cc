// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Centroid-based relax, using Frank Dimaio's updated centroid stats
/// @brief CentroidRelax Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

//Project Headers
#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

//Core Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>

//Protocol Headers
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MonteCarlo.hh>

//Option Headers
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Utility Headers
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>

static basic::Tracer TR( "protocols.relax.CentroidRelax" );
namespace protocols {
namespace relax {
using namespace basic::options;
using namespace core::scoring;
using namespace protocols::simple_moves;
using namespace core::pack::task;
using namespace protocols::moves;
using namespace core::kinematics;
using core::pose::Pose;
using core::Real;
using std::string;


CentroidRelax::CentroidRelax():
	RelaxProtocolBase( std::string("CentroidRelax") )
{
	string def_score = option [OptionKeys::relax::centroid::weights]();
	ScoreFunctionOP score(ScoreFunctionFactory::create_score_function(def_score));
	set_score_function(score);
	set_defaults();

}

CentroidRelax::CentroidRelax(MoveMapOP mm):
	RelaxProtocolBase( std::string("CentroidRelax") )
{
	string def_score = option [OptionKeys::relax::centroid::weights]();
	ScoreFunctionOP score(ScoreFunctionFactory::create_score_function(def_score));
	set_score_function(score);
	set_defaults();
	set_movemap(mm);
}

CentroidRelax::CentroidRelax(MoveMapOP mm, ScoreFunctionOP cen_scorefxn_in):
	RelaxProtocolBase( std::string("CentroidRelax") )
{
	set_score_function(cen_scorefxn_in);
	set_defaults();
	set_movemap(mm);

}

//Mover Methods
CentroidRelax::~CentroidRelax()= default;

protocols::moves::MoverOP
CentroidRelax::clone() const {
	return protocols::moves::MoverOP( new CentroidRelax(*this) );
}

string
CentroidRelax::get_name() const {
	return "CentroidRelax";
}

void
CentroidRelax::set_score_function(ScoreFunctionOP cen_score){
	ScoreFunctionOP score = cen_score->clone();
	//Both Base class and CentroidRelax class have owning pointers to the centroid scorefunction.
	set_scorefxn(score);
	cen_scorefxn_=get_scorefxn();
	cen_scorefxn_->set_weight(chainbreak, 100.00); //I hate chainbreaks.
}

void
CentroidRelax::set_fa_score_function(ScoreFunctionOP fa_score){
	ScoreFunctionOP score = fa_score->clone();
	fa_scorefxn_ = score;
}
void
CentroidRelax::set_defaults(){

	set_rounds(option [OptionKeys::relax::default_repeats]());
	set_ramp_vdw(option [OptionKeys::relax::centroid::ramp_vdw]());
	set_ramp_rama(option [OptionKeys::relax::centroid::ramp_rama]());
	set_cartesian(option[OptionKeys::relax::cartesian]());
	set_fa_score_function(core::scoring::get_score_function());
	set_use_rama2b(option [OptionKeys::score::ramaneighbors]());
	set_use_increased_vdw_radii(option [OptionKeys::relax::centroid::increase_vdw_radii]());
	do_final_repack(option [OptionKeys::relax::centroid::do_final_repack]());
	read_default_parameters();
}


void
CentroidRelax::set_rounds(core::Size rounds){
	rounds_ = rounds;
}

void
CentroidRelax::set_use_increased_vdw_radii(bool use){
	use_increased_vdw_radii_ = use;
}

void
CentroidRelax::set_use_rama2b(bool use){
	//Should be already in score function - But it's not.
	Real rama_weight = cen_scorefxn_->get_weight(rama);
	Real rama2b_weight = cen_scorefxn_->get_weight(rama2b);

	if ( use ) {
		rama_type_=rama2b;
		if ( rama2b_weight==0 ) {
			cen_scorefxn_->set_weight(rama2b, rama_weight);
			cen_scorefxn_->set_weight(rama, 0.0);
		}
	} else {
		rama_type_=rama;
		if ( rama_weight==0 ) {
			cen_scorefxn_->set_weight(rama, rama2b_weight);
			cen_scorefxn_->set_weight(rama2b, 0.0);
		}
	}
}

void
CentroidRelax::set_ramp_rama(bool use){
	ramp_rama_ = use;
}
void
CentroidRelax::setup_class_movemap_and_constraints(Pose & pose){

	movemap_= get_movemap()->clone();
	if ( core::pose::symmetry::is_symmetric( pose )  )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );
	}

	initialize_movemap(pose, *movemap_);
	make_dna_rigid(pose, *movemap_);
	set_movemap(movemap_);
	set_up_constraints(pose, *movemap_);
	cen_scorefxn_->show(TR, pose);

}

void
CentroidRelax::setup_increased_vdw_radii(){
	///Courtesy of Frank Dimaio:
	TR<<"Increasing BB VDW Radii"<<std::endl;
	core::scoring::methods::EnergyMethodOptions lowres_options(cen_scorefxn_->energy_method_options());
	lowres_options.atom_vdw_atom_type_set_name("centroid_min");
	cen_scorefxn_->set_energy_method_options(lowres_options);
}

void
CentroidRelax::set_ramp_vdw(bool use){
	ramp_vdw_ = use;
}

void
CentroidRelax::set_cartesian(bool cart){
	cartesian( cart );

	if ( cart ) {
		set_min_type("lbfgs_armijo_nonmonotone");
		ScoreFunctionOP score(ScoreFunctionFactory::create_score_function("score4_smooth_cart"));
		set_score_function(score);
	} else {
		set_min_type(option [OptionKeys::relax::min_type]());
		string def_score = option [OptionKeys::relax::centroid::weights]();
		ScoreFunctionOP score(ScoreFunctionFactory::create_score_function(def_score));
		set_score_function(score);
	}
}
void
CentroidRelax::set_min_type(string min){
	min_type( min );
}

void
CentroidRelax::read_default_parameters(){
	string const fname = option [ OptionKeys::relax::centroid::parameters]();
	utility::io::izstream param_stream;
	basic::database::open(param_stream, fname);
	Real vdw_ramp; Real rama_ramp; Real min; Real cst;
	core::Size line_count = 0;
	while ( ! param_stream.eof() ) {
		++line_count;
		param_stream >> vdw_ramp >> rama_ramp >> min >> cst;
		def_parameters.vdw_params.push_back(vdw_ramp);
		def_parameters.rama_params.push_back(rama_ramp);
		def_parameters.min_params.push_back(min);
		def_parameters.cst_params.push_back(cst);

	}
	param_stream.close();

}

void
CentroidRelax::do_final_repack(bool repack_sc){
	repack_sc_=repack_sc;
}

void
CentroidRelax::apply(Pose& pose){


	if ( use_increased_vdw_radii_ ) {
		setup_increased_vdw_radii();
	}


	//Get the rama_weight.  If it's rama2b, get that.
	Real rama_weight = cen_scorefxn_->get_weight(rama);
	if ( rama_weight==0 ) {
		rama_weight = cen_scorefxn_->get_weight(rama2b);
		if ( rama_weight==0 ) {
			set_ramp_rama(false);
		}
	}

	//Get vdw weight.
	Real vdw_weight = cen_scorefxn_->get_weight(vdw);
	if ( vdw_weight==0 ) {
		set_ramp_vdw(false);
	}


	MonteCarloOP cen_mc;


	ReturnSidechainMoverOP recover_sc;


	bool passed_centroid;
	if ( ! pose.is_centroid() ) {
		TR << "FullAtom Score::"<<std::endl;

		fa_scorefxn_->show(TR, pose);

		recover_sc = ReturnSidechainMoverOP( new ReturnSidechainMover(pose) );

		SwitchResidueTypeSetMoverOP to_cen( new SwitchResidueTypeSetMover("centroid") );

		to_cen->apply(pose);

		passed_centroid = false;

		TR << "Centroid Score::"<<std::endl;
		//cen_scorefxn_->show(TR, pose);
		cen_mc = MonteCarloOP( new MonteCarlo(pose, *cen_scorefxn_, 1.0) );

	} else {

		passed_centroid = true;
		cen_scorefxn_->show(TR, pose);
		cen_mc = MonteCarloOP( new MonteCarlo(pose, *cen_scorefxn_, 1.0) );

	}


	setup_class_movemap_and_constraints(pose);

	//Initialize MinMover
	MinMoverOP minmover;
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		minmover = MinMoverOP( new symmetry::SymMinMover( movemap_, cen_scorefxn_, min_type(), def_parameters.min_params[1], true ) );
	} else {
		minmover = MinMoverOP( new MinMover( movemap_, cen_scorefxn_, min_type(), def_parameters.min_params[1], true ) );
	}
	if ( cartesian() ) {
		minmover->cartesian( true );
		minmover->min_type("lbfgs_armijo_nonmonotone");
	}

	//Test to make sure something is ramped - If not, run basic.
	if ( !ramp_vdw_ && !ramp_rama_ ) {
		TR << "Running BASIC Relax" <<std::endl;
		for ( core::Size i=1; i<=rounds_; i++ ) {
			minmover->apply(pose);
			TR <<"Cen Score: "<<(*cen_scorefxn_)(pose)<<std::endl;
			cen_mc->boltzmann(pose);
		}

		cen_mc->recover_low(pose);
	} else {
		TR <<"Starting RAMPED Protocol"<<std::endl;
		EnergyMap cen_emap = cen_scorefxn_->weights();
		TR <<"Total Rounds : " <<rounds_ <<std::endl;
		TR <<"Total Ramps  : " <<def_parameters.vdw_params.size() <<std::endl;
		//Outer loop
		for ( core::Size i=1; i<= rounds_; i++ ) {
			TR <<"Starting Round "<<i<<std::endl;

			//Inner loop
			for ( core::Size i2 = 1; i2 <=def_parameters.vdw_params.size(); i2++ ) {

				//Set Weights
				TR<< "Ramp "<<i2<<std::endl;
				if ( ramp_vdw_ ) {
					cen_scorefxn_->set_weight(vdw, cen_emap[vdw]*def_parameters.vdw_params[i2]);
				}
				if ( ramp_rama_ ) {
					cen_scorefxn_->set_weight(rama_type_, cen_emap[rama_type_]*def_parameters.rama_params[i2]);
				}
				//Removing ramp, as without any coordinate constraints structures core will compress
				//if (constrain_coords_ && ramp_down_constraints_){
				//cen_scorefxn_->set_weight(coordinate_constraint, cen_emap[coordinate_constraint]*def_parameters.cst_params[i2]);
				//}
				minmover->tolerance(def_parameters.min_params[i2]);
				minmover->apply(pose);
				TR << (*cen_scorefxn_)(pose) << std::endl;
				//Use MC on last run->Check how many times we ramp from the vdw_params ->Doesn't matter which param we check.
				if ( i2==def_parameters.vdw_params.size() ) {
					TR <<"Cen Energy with full weights :"<<(*cen_scorefxn_)(pose)<<std::endl;
					cen_scorefxn_->show(TR, pose);
					cen_mc->boltzmann(pose);
				}
			}
			//Recover lowest scoring decoy if last round.
			if ( i==rounds_ ) {
				cen_mc->recover_low(pose);
			}
		}
	}

	if ( !passed_centroid ) {
		//TR <<"Scoring Final Pose :"<<pose.is_centroid()<<std::endl;
		cen_scorefxn_->show(TR, pose);
		recover_sc->apply(pose);
		fa_scorefxn_->show(TR, pose);
		if ( repack_sc_ ) {
			TR <<"Repacking Sidechains"<<std::endl;
			PackerTaskOP task = TaskFactory::create_packer_task( pose );
			task->restrict_to_repacking();
			task->temporarily_fix_everything();
			for ( core::Size i=1; i<=pose.size(); ++i ) {
				if ( movemap_->get_chi(i) ) {
					task->temporarily_set_pack_residue(i, true);
				}
			}
			protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover(fa_scorefxn_, task) );
			packer->apply(pose);
			fa_scorefxn_->show(TR, pose);
		}

		TR << "Centroid relax complete.  Returning all-atom Structure."<<std::endl;
		return;
	} else {

		TR << "Centroid relax complete.  Returning centroid Structure." <<std::endl;
		cen_scorefxn_->show(TR, pose);
		return;
	}


}//END apply method
}
}
