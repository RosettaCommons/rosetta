// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file Centroid-based relax, using Frank Dimaio's updated centroid stats
/// @brief CentroidRelax Protocol.  
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>

//Core
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/conformation/Residue.hh>

//Movers
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

//Options
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Basic
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>

static basic::Tracer TR("protocols.relax.CentroidRelax");
namespace protocols {
namespace relax {
    using namespace basic::options;
    using namespace core::scoring;
    using namespace protocols::simple_moves;
    using namespace core::pack::task;
    using namespace protocols::moves;
    using core::pose::Pose;
    using core::Real;
    using std::string;
    

    CentroidRelax::CentroidRelax():
        RelaxProtocolBase("CentroidRelax")
    {
        string def_score = option [OptionKeys::relax::centroid::weights]();
        ScoreFunctionOP score(ScoreFunctionFactory::create_score_function(def_score));
        set_score_function(score);
        initialize_objects();
        set_defaults();
        
    }
    
    CentroidRelax::CentroidRelax(MoveMapOP mm):
        RelaxProtocolBase("CentroidRelax")
    {
        string def_score = option [OptionKeys::relax::centroid::weights]();
        ScoreFunctionOP score(ScoreFunctionFactory::create_score_function(def_score));
        set_score_function(score);
        initialize_objects();
        set_defaults();
        set_movemap(mm);
    }
    
    CentroidRelax::CentroidRelax(MoveMapOP mm, ScoreFunctionOP cen_scorefxn_in):
        RelaxProtocolBase("CentroidRelax")
    {
        set_score_function(cen_scorefxn_in);
        initialize_objects();
        set_defaults();
        set_movemap(mm);
        
    }
    
    //Mover Methods
    CentroidRelax::~CentroidRelax(){}
    
    protocols::moves::MoverOP
    CentroidRelax::clone() const {
        return new CentroidRelax(*this);
    }
    
    string
    CentroidRelax::get_name() const {
        return "CentroidRelax";
    }
    void
    CentroidRelax::initialize_objects(){
        to_fa_ = new SwitchResidueTypeSetMover("fa_standard");
        to_cen_ = new SwitchResidueTypeSetMover("centroid");
        
    }
    
    void
    CentroidRelax::set_score_function(ScoreFunctionOP cen_score){
        ScoreFunctionOP score = cen_score->clone();
        scorefxn_ = score;
        scorefxn_->set_weight(chainbreak, 100.00); //I hate chainbreaks.
    }
    
    void
    CentroidRelax::set_fa_score_function(ScoreFunctionOP fa_score){
        ScoreFunctionOP score = fa_score->clone();
        fa_scorefxn_ = score;
    }
    void
    CentroidRelax::set_defaults(){
        set_centroid_only_mode(option [OptionKeys::relax::centroid::centroid_only_mode]());
        
        set_rounds(option [OptionKeys::relax::default_repeats]());
        set_ramp_vdw(option [OptionKeys::relax::centroid::ramp_vdw]());
        set_ramp_rama(option [OptionKeys::relax::centroid::ramp_rama]());
        set_cartesian(option[OptionKeys::relax::cartesian]());
        set_fa_score_function(core::scoring::getScoreFunction());
        set_use_rama2b(option [OptionKeys::score::ramaneighbors]());
        set_default_parameters();
        sc_recoverable_=false;
    }
    
    void
    CentroidRelax::set_centroid_only_mode(bool cen_only){
        centroid_only_mode_ = cen_only;
        if (cen_only) {
            set_repack_intermittent(false);
            set_final_min_sidechains(false);
            set_final_repack_sidechains(false);
            set_final_min_fullatom(false);
        }
        else {
            set_repack_intermittent(option[OptionKeys::relax::centroid::do_inter_repacks]() );
            set_final_min_sidechains(option [OptionKeys::relax::centroid::do_final_min_sc]());
            set_final_repack_sidechains(option [OptionKeys::relax::centroid::do_final_repack]());
            set_final_min_fullatom(option [OptionKeys::relax::centroid::do_final_min_fa]());
        }
    }
    
    void
    CentroidRelax::set_rounds(core::Size rounds){
        rounds_ = rounds;
    }
    void
    CentroidRelax::set_use_rama2b(bool use){
        //Should be already in score function - But it's not.
        Real rama_weight = scorefxn_->get_weight(rama);
        Real rama2b_weight = scorefxn_->get_weight(rama2b);
        
        if (use){
            rama_type_=rama2b;
            if (rama2b_weight==0){
                scorefxn_->set_weight(rama2b, rama_weight);
                scorefxn_->set_weight(rama, 0.0);
            }
        }
        else{
            rama_type_=rama;
            if (rama_weight==0){
                scorefxn_->set_weight(rama, rama2b_weight);
                scorefxn_->set_weight(rama2b, 0.0);
            }
        }
    }
        
    void
    CentroidRelax::set_ramp_rama(bool use){
        ramp_rama_ = use;
    }
    void
    CentroidRelax::setup_class_movemap(Pose & pose){
        
        movemap_= get_movemap()->clone();
        if ( core::pose::symmetry::is_symmetric( pose )  )  {
            core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );
        }
        
        initialize_movemap(pose, *movemap_); //From RelaxProtocolBase
        makeDnaRigid(pose);
        set_up_constraints(pose, *movemap_); //FromRelaxProtocolBase
        set_movemap(movemap_);
        ///From RelaxProtocolBase
        if (constrain_coords_){   
            if ( scorefxn_->get_weight( coordinate_constraint ) == 0 ) {
		scorefxn_->set_weight( coordinate_constraint, 0.5 );
            }
            if ( fa_scorefxn_->get_weight( coordinate_constraint ) == 0 ) {
                fa_scorefxn_->set_weight( coordinate_constraint, 0.5);
            }
        }

    }
    
    void
    CentroidRelax::set_ramp_vdw(bool use){
        ramp_vdw_ = use;
    }
    void
    CentroidRelax::set_repack_intermittent(bool use){
        repack_intermittent_ = use;
    }
    void
    CentroidRelax::set_recover_sidechains_mover(ReturnSidechainMoverOP recover_sc){
        recover_sc_ = recover_sc;
        sc_recoverable_=true;
    }
    void
    CentroidRelax::set_cartesian(bool cart){
        cartesian_ = cart;
        
        if (cart){
            set_min_type("lbfgs_armijo_nonmonotone");
            ScoreFunctionOP score(ScoreFunctionFactory::create_score_function("score4_smooth_cart"));
            set_score_function(score);
        }
        else {
            set_min_type(option [OptionKeys::relax::min_type]());
            string def_score = option [OptionKeys::relax::centroid::weights]();
            ScoreFunctionOP score(ScoreFunctionFactory::create_score_function(def_score));
            set_score_function(score);
        }
    }
    void
    CentroidRelax::set_min_type(string min){
        min_type_ = min; 
    }
    void
    CentroidRelax::set_return_full_atom_from_centroid(bool fa_ret){
        fa_return_ = fa_ret;
    }
    void
    CentroidRelax::set_final_min_fullatom(bool min){
        final_min_fa_ = min;
    }
    void
    CentroidRelax::set_final_min_sidechains(bool min){
        final_min_sc_ = min;
    }
    void
    CentroidRelax::set_final_repack_sidechains(bool repack){
        final_repack_sc_ = repack;
    }
    
    void
    CentroidRelax::set_default_parameters(){
        string const fname = option [ OptionKeys::relax::centroid::parameters]();
        utility::io::izstream param_stream;
        basic::database::open(param_stream, fname);
        Real vdw_ramp; Real rama_ramp; Real min; Real cst;
        core::Size line_count = 0;
        while (! param_stream.eof() ){
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
    CentroidRelax::makeDnaRigid(Pose & pose){
        //This is exactly from FastRelax.
	using namespace core::conformation;
	//if DNA present set so it doesn't move
	for ( core::Size i=1; i<=pose.total_residue() ; ++i )      {
		if( pose.residue(i).is_DNA()){
			TR << "turning off DNA bb and chi move" << std::endl;
			movemap_->set_bb( i, false );
			movemap_->set_chi( i, false );
		}
	}
	//end DNA rigidity
    }
    
    
///Main CentroidRelax.

    void
    CentroidRelax::movemap_to_taskfactory(Pose & pose, PackerTaskOP task){
        //From RelaxProtocolBase.cc
        utility::vector1<bool> allow_repack(pose.total_residue(), false);
        for ( core::Size i = 1; i<= pose.total_residue() ; ++i ) {
            allow_repack[i] = movemap_->get_chi(i);
        }
        task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
    }
    
    
    void
    CentroidRelax::apply(Pose& pose){
        setup_class_movemap(pose);

        
        
        //Get the rama_weight.  If it's rama2b, get that.
        Real rama_weight = scorefxn_->get_weight(rama);
        if (rama_weight==0){
            rama_weight = scorefxn_->get_weight(rama2b);
            if (rama_weight==0){
                set_ramp_rama(false);
            }
        }
        
        //Get vdw weight.
        Real vdw_weight = scorefxn_->get_weight(vdw);
        if (vdw_weight==0){
            set_ramp_vdw(false);
        }
        
        //Initialize MinMover
        MinMoverOP minmover;
        if ( core::pose::symmetry::is_symmetric( pose ) )  {
            minmover = new symmetry::SymMinMover( movemap_, scorefxn_, min_type_, def_parameters.min_params[1], true );}
        else {
            minmover = new MinMover( movemap_, scorefxn_, min_type_, def_parameters.min_params[1], true );
        }
	minmover->cartesian( cartesian_ );
        
        TR <<"MinMover Initialized"<<std::endl;
        

        
        
        
        //Initialize Packer
          //Declaring the objects just in case.
                // For memory, this maybe should be changed.
                // For speed, they only need to be declared once. 

        

        //Initialize Starting residue set + MonteCarlo for both.
        MonteCarloOP fa_mc;
        MonteCarloOP cen_mc;
        PackRotamersMoverOP packer;  
        PackerTaskOP task;
        
        bool passed_centroid;
        if (! pose.is_centroid()){
            TR << "FullAtom Score::"<<std::endl;
            task = TaskFactory::create_packer_task( pose );
            movemap_to_taskfactory(pose, task);
            
            fa_scorefxn_->show(TR, pose);
            fa_mc = new MonteCarlo(pose, *fa_scorefxn_, 1.0);
            
            recover_sc_ = new ReturnSidechainMover(pose);
            to_cen_->apply(pose);
            
            passed_centroid = false;
            sc_recoverable_=true;
            
            TR << "Centroid Score::"<<std::endl;
            scorefxn_->show(TR, pose);
            cen_mc = new MonteCarlo(pose, *scorefxn_, 1.0);
        
        } else if (!centroid_only_mode_ && !recover_sc_){
            TR << "Please pass a Full Atom representation of the pose or give a ReturnSidechainMoverOP"<<std::endl;
            return;
        } else {
            
            passed_centroid = true;
            scorefxn_->show(TR, pose);
            cen_mc = new MonteCarlo(pose, *scorefxn_, 1.0);
            
            if (sc_recoverable_){
                recover_sc_->apply(pose);
                task = TaskFactory::create_packer_task( pose );
                movemap_to_taskfactory(pose, task);
                fa_mc = new MonteCarlo(pose, *fa_scorefxn_, 1.0);
                recover_sc_ = new ReturnSidechainMover(pose);
                to_cen_->apply(pose);
            }
        }
        
        TR << "MonteCarlo Initialized"<<std::endl;
        
        //Create packer task, set up the PackRotamers Mover
        

        if (sc_recoverable_) {
            if ( core::pose::symmetry::is_symmetric( pose ) )  {
                packer = new simple_moves::symmetry::SymPackRotamersMover( fa_scorefxn_, task );
            }
            else{
                packer = new PackRotamersMover(fa_scorefxn_, task);
            }
        }
        //Initialize MonteCarlo (Not screwing with kT yet)
        
        TR <<"Packer Initialized"<<std::endl;

        //Test to make sure something is ramped - If not, run basic.
            //This could be rounds_*parameters.rama_params.Size() to make it very much like other protocol.  
            //Though, using the basic, we can run it at least 4x quicker probably (check Energy benefits).
        if (!ramp_vdw_ && !ramp_rama_ && !constrain_coords_){
            TR << "Running BASIC Relax" <<std::endl;
            for (core::Size i=1; i<=rounds_; i++){
                minmover->apply(pose);
                TR <<"Cen Score: "<<(*scorefxn_)(pose)<<std::endl;
                if (repack_intermittent_ && sc_recoverable_){
                    recover_sc_->apply(pose);
                    packer->apply(pose);
                    TR <<"FA  Score: "<<(*fa_scorefxn_)(pose)<<std::endl;
                    fa_mc->boltzmann(pose);
                    recover_sc_ = new ReturnSidechainMover(pose);
                    to_cen_->apply(pose);
                }
                else{ cen_mc->boltzmann(pose); }
            }
            if (repack_intermittent_ && sc_recoverable_) {
                fa_mc->recover_low(pose);
                to_cen_->apply(pose);
            } 
            else { cen_mc->recover_low(pose);}
            rounds_=0;
        }
        else{
            TR <<"Starting RAMPED Protocol"<<std::endl;
            EnergyMap cen_emap = scorefxn_()->weights();
            EnergyMap fa_emap  = fa_scorefxn_()->weights();
            TR <<"Total Rounds : " <<rounds_ <<std::endl;
            TR <<"Total Ramps  : " <<def_parameters.vdw_params.size() <<std::endl;
            //Outer loop
            for (core::Size i=1; i<= rounds_; i++){
                TR <<"Starting Round "<<i<<std::endl;
                
                //Inner loop
                for (core::Size i2 = 1; i2 <=def_parameters.vdw_params.size(); i2++){
                    
                    //Set Weights
                    TR<< "Ramp "<<i2<<std::endl;
                    if (ramp_vdw_){
                        scorefxn_->set_weight(vdw, cen_emap[vdw]*def_parameters.vdw_params[i2]);
                        fa_scorefxn_->set_weight(fa_rep, fa_emap[fa_rep]*def_parameters.vdw_params[i2]);
                    }
                    if (ramp_rama_){
                        scorefxn_->set_weight(rama_type_, cen_emap[rama_type_]*def_parameters.rama_params[i2]);
                    }
                    if (constrain_coords_){
                        scorefxn_->set_weight(coordinate_constraint, cen_emap[coordinate_constraint]*def_parameters.cst_params[i2]);
                        fa_scorefxn_->set_weight(coordinate_constraint, fa_emap[coordinate_constraint]*def_parameters.cst_params[i2]);
                    }
                    minmover->tolerance(def_parameters.min_params[i2]);
                    minmover->apply(pose);
                    TR << (*scorefxn_)(pose) << std::endl;
                    //Use MC on last run.
                    if (i2==def_parameters.vdw_params.size()){ 
                        if (repack_intermittent_  && sc_recoverable_){
                            TR <<"Cen Energy with full weights :"<<(*scorefxn_)(pose)<<std::endl;
                            recover_sc_->apply(pose);
                            
                            if (repack_intermittent_){
                                packer->score_function(fa_scorefxn_);
                                packer->apply(pose);
                            }
                            
                            TR <<"FA Energy with full weights :"<<(*fa_scorefxn_)(pose)<<std::endl;
                            fa_mc->boltzmann(pose);
                            recover_sc_ = new ReturnSidechainMover(pose);
                            to_cen_->apply(pose);
                        }
                        else {
                            TR <<"Cen Energy with full weights :"<<(*scorefxn_)(pose)<<std::endl;
                            cen_mc->boltzmann(pose);
                        }
                    }
                }
                //Recover lowest scoring decoy if last round.
                if (i==rounds_){
                    if (repack_intermittent_   && sc_recoverable_){
                        fa_mc->recover_low(pose);
                        recover_sc_ = new ReturnSidechainMover(pose);
                        to_cen_->apply(pose);
                    }
                    else{
                        cen_mc->recover_low(pose);
                    }
                }
            }

        
            TR <<"Final Modeling"<<std::endl;
            //movemap_->show(TR, pose.total_residue());
            //End of protocol - Convert to FA, repack/Minimize SC or return centroid pose.
            if ((!passed_centroid) || (fa_return_ && sc_recoverable_)){
                //TR <<"Scoring Final Pose :"<<pose.is_centroid()<<std::endl;
                scorefxn_->show(TR, pose);
                recover_sc_->apply(pose);
                
                if (final_min_fa_){
                    MoveMapOP mm = movemap_->clone();
                    minmover->score_function(fa_scorefxn_);
                    minmover->movemap(mm);
                    minmover->apply(pose);
                }
                if (final_repack_sc_){
                    packer->score_function(fa_scorefxn_);
                    packer->apply(pose);
                }
            
                if (final_min_sc_){
                
                    MoveMapOP mm = movemap_->clone();
                    mm->set_bb(false);
                    minmover->score_function(fa_scorefxn_);
                    minmover->movemap(mm);
                    minmover->apply(pose);
                }
                fa_scorefxn_->show(TR, pose);
                TR << "Centroid relax complete.  Returning all-atom Structure."<<std::endl;
                return;
            
            } else {
            
                TR << "Centroid relax complete.  Returning centroid Structure." <<std::endl;
                scorefxn_->show(TR, pose);
                return;}
        

            }
        }
    }
}
