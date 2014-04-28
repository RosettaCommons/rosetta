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
/// @author Yifan Song

#include <protocols/hybridization/BackboneTorsionSampler.hh>
#include <protocols/hybridization/BackboneTorsionSamplerCreator.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/import_pose/import_pose.hh>
#include <core/sequence/util.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/hybridization/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <ObjexxFCL/format.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.hybridization.BackboneTorsionSampler" );
static numeric::random::RandomGenerator RG(5694684);

namespace protocols {
namespace hybridization {

BackboneTorsionSampler::BackboneTorsionSampler() {
    init();
}

BackboneTorsionSampler::~BackboneTorsionSampler() {}

void BackboneTorsionSampler::init() {
    set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" ) );
    increase_cycles_ = 1.0;
    temperature_ = 1.0;
    recover_low_ = true;
    local_ = 0;
    dump_snapshots_ = false;
    snapshot_interval_ = 10;
}

void BackboneTorsionSampler::local_perturb(core::pose::Pose pose, core::Real max_delta_torsion) {
    Size n_torsion(0);
    for (Size ires=1; ires <= pose.total_residue(); ++ires) {
        n_torsion += 2; // for now
    }
    Size perturbed_torsion = RG.random_range(1, n_torsion);
    Size perturbed_residue = (perturbed_torsion + 1)/2;
    core::Real delta = (2.* RG.uniform() - 1.) * max_delta_torsion;
    if ( (perturbed_torsion - (perturbed_residue-1)*2) == 1 ) {
        pose.set_phi(perturbed_residue, pose.phi(perturbed_residue) + delta);
    }
    else {
        //core::Real psi = (2.* RG.uniform() - 1.) * max_delta_torsion + pose.psi(perturbed_residue);
        pose.set_psi(perturbed_residue, pose.psi(perturbed_residue) + delta);
    }
    
    core::Real variance(3.);
    int next_torsion = perturbed_torsion;
    while (next_torsion == (int) perturbed_torsion) {
        next_torsion = floor( perturbed_res_ + variance*RG.gaussian() + 0.5 );
    }
    
    if ( next_torsion >= 1 && next_torsion <= (int) pose.total_residue() * 2 ) {
        Size next_perturbed_residue = (next_torsion + 1)/2;
        if ( (next_torsion - (next_perturbed_residue-1)*2) == 1 ) {
            pose.set_phi(next_perturbed_residue, pose.phi(next_perturbed_residue) - delta);
        }
        else {
            //core::Real psi = (2.* RG.uniform() - 1.) * max_delta_torsion + pose.psi(next_perturbed_residue);
            pose.set_psi(next_perturbed_residue, pose.psi(next_perturbed_residue) - delta);
        }
    }
}
    
void BackboneTorsionSampler::perturb(core::pose::Pose & pose,
                                     core::Size level, // level 1 is the base
                                     core::Real max_delta_torsion,
                                     core::Size local,
                                     bool rama_biased,
                                     bool repack,
                                     bool minimize) {
    if (level == 1){
        perturbed_res_ = RG.random_range(1, pose.total_residue());
    }
    for (core::Size ires=1; ires <= pose.total_residue() ; ++ires) {
        if (local != 0) {
            if ( abs((int) ires - (int) perturbed_res_) > (int) ((local - 1) * level/2)) continue;
        }
        
        if (pose.residue_type(ires).is_protein()) {
            core::Real phi(0), psi(0);
            if (level == 1 && rama_biased) {
                core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();
                
                rama.random_phipsi_from_rama( pose.residue_type(ires).aa(), phi, psi);
            }
            else {
                phi = (2.* RG.uniform() - 1.) * max_delta_torsion + pose.phi(ires);
                psi = (2.* RG.uniform() - 1.) * max_delta_torsion + pose.psi(ires);
            }
            pose.set_phi(ires, phi);
            pose.set_psi(ires, psi);
        }
    }
    
    pose.conformation().detect_disulfides();
    if (repack) {
        pack_full_repack_->apply(pose);
    }
    if (minimize) {
        mm_.set_bb(true);
        minimizer_->run( pose, mm_, *scorefxn_, *options_ );
    }
}
    
void BackboneTorsionSampler::apply( core::pose::Pose & pose ) {
    //core::pack::task::TaskFactoryOP local_tf = new core::pack::task::TaskFactory();
    //local_tf->push_back(new core::pack::task::operation::RestrictToRepacking());
    // if present, task_factory_ always overrides/regenerates task_
	core::pack::task::PackerTaskOP task;
    if ( task_factory_ != 0 ) {
		task = task_factory_->create_task_and_apply_taskoperations( pose );
        task->restrict_to_repacking();
	} else {
        task = core::pack::task::TaskFactory::create_packer_task( pose );
        task->initialize_from_command_line().restrict_to_repacking();
	}
    pack_full_repack_ = new protocols::simple_moves::PackRotamersMover( scorefxn_, task );
    //task->show_all_residue_tasks();
    
    pose.conformation().detect_disulfides();

    // initialize monte carlo
    utility::vector1<protocols::moves::MonteCarloOP> mc(n_nested_+1);
    utility::vector1<core::Size> ncycles(n_nested_+1);
    utility::vector1<core::Size> counters(n_nested_+1, 0);
    
    perturbed_res_ = 0;
    for (core::Size i_nest=1; i_nest<=n_nested_+1; ++i_nest) {
        mc[i_nest]= new protocols::moves::MonteCarlo( pose, *scorefxn_, temperature_ );
        mc[i_nest]->set_autotemp( false, temperature_ );
        if (i_nest == 1) {
            ncycles[i_nest] = 10 * pose.total_residue() * increase_cycles_;
        }
        else if (i_nest == n_nested_+1) {
            ncycles[i_nest] = 10 * increase_cycles_;
        }
        else {
            ncycles[i_nest] = 5 * increase_cycles_;
        }
    }

    minimizer_ = new core::optimization::AtomTreeMinimizer();
    //minimizer_ = new core::optimization::CartesianMinimizer();
	options_ =  new core::optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	
	options_->max_iter(10);
	mm_.set_bb  ( true );
	mm_.set_chi ( true );
	mm_.set_jump( false );

    if (dump_snapshots_) {
        std::string pdb_fn = snapshot_prefix_;
        for (core::Size i_nest=1; i_nest<=n_nested_+1; ++i_nest) {
            pdb_fn += "_" + ObjexxFCL::format::I(4,4,counters[i_nest]);
        }
        pdb_fn += ".pdb";
        pose.dump_pdb(pdb_fn);
    }

    Size snapshot_counter = 0;
    for (core::Size i_nest=1; i_nest<n_nested_; ++i_nest) {
        if (i_nest == 1) {
            perturb(pose, i_nest, 0., local_, true, false, false);
        }
        else {
            perturb(pose, i_nest, 10., local_, false, false, false);
        }
    }
    while (true) {
        Size i_nest = n_nested_+1;
        ++counters[i_nest];

        // lowest level move
        if (n_nested_ == 0) {
            perturb(pose, i_nest, 0., local_, true, false, false);
        }
        else {
            if (counters[i_nest] < ncycles[i_nest]/2) {
                perturb(pose, i_nest, 5., 0, false, true, false);
            }
            else {
                perturb(pose, i_nest, 5., 0, false, true, true);
            }
        }
        
        core::Real score=(*scorefxn_)(pose);
        if (native_ && native_->total_residue()) {
            core::Real gdtmm(0.);
            core::sequence::SequenceAlignmentOP native_aln;
			gdtmm = get_gdtmm(*native_, pose, native_aln);
            using namespace ObjexxFCL::format;
			TR << "Trace: " << F(8,3,gdtmm) << " " << F(11,3,score) << std::endl;
		}

        mc[i_nest]->boltzmann(pose, "BackboneTorsionSampler");
        
        if (n_nested_ != 0) {
            while (true) {
                if (counters[i_nest] >= ncycles[i_nest]) {
                    counters[i_nest] = 0;
                    mc[i_nest]->show_scores();
                    mc[i_nest]->show_counters();
                    
                    mc[i_nest]->recover_low(pose);

                    --i_nest;
                    ++counters[i_nest];
                    
                    (*scorefxn_)(pose);
                    mc[i_nest]->boltzmann(pose, "BackboneTorsionSampler");
                    
                    if (i_nest == 1) {
                        perturb(pose, i_nest, 0., local_, true, false, false);
                    }
                    else {
                        perturb(pose, i_nest, 10., local_, false, false, false);
                    }

                    if (i_nest == 1) break; // do not evaluate the bottom level counter
                }
                else {
                    break;
                }
            }

            for (core::Size j_nest=n_nested_+1; j_nest>=1; --j_nest) {
                if (counters[j_nest] == 0) {
                    mc[j_nest]->reset(pose);
                }
                else {
                    break;
                }
            }
            
            if (counters[1] >= ncycles[1]) {
                break;
            }
        }
        //TR << counters[1] << " " << counters[2] << std::endl;
        
        if (dump_snapshots_) {
            ++snapshot_counter;
            if ( snapshot_counter % snapshot_interval_ == 0 ) {
                std::string pdb_fn = snapshot_prefix_;
                for (core::Size i_nest=1; i_nest<=n_nested_+1; ++i_nest) {
                    pdb_fn += "_" + ObjexxFCL::format::I(4,4,counters[i_nest]);
                }
                pdb_fn += ".pdb";
                pose.dump_pdb(pdb_fn);
            }
        }
    }
    
    for (core::Size i_nest=1; i_nest<=n_nested_+1; ++i_nest) {
        mc[i_nest]->show_scores();
        mc[i_nest]->show_counters();
    }
    if (recover_low_) mc[1]->recover_low(pose);
    if (dump_snapshots_) {
        std::string pdb_fn = snapshot_prefix_;
        for (core::Size i_nest=1; i_nest<=n_nested_+1; ++i_nest) {
            pdb_fn += "_" + ObjexxFCL::format::I(4,4,counters[i_nest]);
        }
        pdb_fn += ".pdb";
        pose.dump_pdb(pdb_fn);
    }
}

void BackboneTorsionSampler::task_factory( core::pack::task::TaskFactoryCOP tf )
{
    runtime_assert( tf );
    task_factory_ = tf;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
BackboneTorsionSampler::parse_my_tag(
    TagCOP tag,
    basic::datacache::DataMap & datamap,
    Filters_map const & ,
    moves::Movers_map const & /*movers*/,
    Pose const & pose
) {
    core::Size start_res = 1;
    core::Size stop_res = pose.total_residue();
    
    if( tag->hasOption( "start_res" ) ) start_res = tag->getOption< core::Size >( "start_res" );
    if( tag->hasOption( "stop_res" ) ) stop_res = tag->getOption< core::Size >( "stop_res" );
    for (core::Size ires=start_res; ires<=stop_res; ++ires) {
            residue_list_.push_back(ires);
    }
    if (tag->hasOption( "native") ) {
        native_ = new core::pose::Pose;
        core::import_pose::pose_from_pdb( *native_, tag->getOption< std::string >( "native" ) );
    }
    
    //String const  user_defined_mover_name_( tag->getOption< String >( "mover_name" ,""));
	//Movers_map::const_iterator  find_mover ( movers.find( user_defined_mover_name_ ));
	//if( find_mover == movers.end() && user_defined_mover_name_ != "" ) {
	//	TR.Error << "ERROR !! mover not found in map: \n" << tag << std::endl;
	//	runtime_assert( find_mover != movers.end() );
	//}

    if( tag->hasOption( "increase_cycles" ) ) increase_cycles_ = tag->getOption< core::Real >( "increase_cycles" );
    if( tag->hasOption( "recover_low" ) ) recover_low_ = tag->getOption< bool >( "recover_low" );
    if( tag->hasOption( "temp" ) ) temperature_ = tag->getOption< core::Real >( "temp" );
    local_ = tag->getOption< core::Size >( "local" , 0);
    n_nested_ = tag->getOption< core::Size >( "nested" , 2);
    
    if( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		set_scorefunction ( (datamap.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name ))->clone() );
	}

    core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0) return;
	task_factory( new_task_factory );

    if( tag->hasOption( "dump_snapshots" ) ) {
        dump_snapshots_ = tag->getOption< core::Size >( "dump_snapshots" );
        snapshot_prefix_ = tag->getOption< std::string >( "snapshot_prefix", "snapshot");
        snapshot_interval_ = tag->getOption< core::Size >( "snapshot_interval" , 100);
    }
}

moves::MoverOP BackboneTorsionSampler::clone() const {
	return new BackboneTorsionSampler( *this );
}
moves::MoverOP BackboneTorsionSampler::fresh_instance() const {
	return new BackboneTorsionSampler;
}

std::string
BackboneTorsionSampler::get_name() const {
    return "BackboneTorsionSampler";
}

protocols::moves::MoverOP
BackboneTorsionSamplerCreator::create_mover() const {
	return new BackboneTorsionSampler;
}

std::string
BackboneTorsionSamplerCreator::keyname() const {
	return BackboneTorsionSamplerCreator::mover_name();
}

std::string
BackboneTorsionSamplerCreator::mover_name() {
	return "BackboneTorsionSampler";
}

} // moves
} // protocols
