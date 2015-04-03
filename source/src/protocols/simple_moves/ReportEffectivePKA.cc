// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Yifan Song

#include <protocols/simple_moves/ReportEffectivePKA.hh>
#include <protocols/simple_moves/ReportEffectivePKACreator.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/cryst/util.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <utility/vector1.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

static thread_local basic::Tracer TR( "protocols.simple_moves.ReportEffectivePKA" );

namespace protocols {
namespace simple_moves {

using namespace core;

ReportEffectivePKA::ReportEffectivePKA() : moves::Mover() {
	init();
}

void ReportEffectivePKA::init() {
    IonizableResidue his = IonizableResidue( "HIS", 6.1, 1.0 );
    his.add_neutral_restype("HIS");
    his.add_neutral_restype("HIS_D");
    his.add_ionized_restype("HIS_P");
    ionizables_.push_back(his);

    IonizableResidue asp = IonizableResidue( "ASP", 3.9, -1.0 );
    asp.add_neutral_restype("ASP_P1");
    asp.add_neutral_restype("ASP_P2");
    asp.add_ionized_restype("ASP");
    ionizables_.push_back(asp);

    IonizableResidue glu = IonizableResidue( "GLU", 4.2, -1.0 );
    glu.add_neutral_restype("GLU_P1");
    glu.add_neutral_restype("GLU_P2");
    glu.add_ionized_restype("GLU");
    ionizables_.push_back(glu);

    IonizableResidue lys = IonizableResidue( "LYS", 10.1, 1.0 );
    lys.add_neutral_restype("LYS_D");
    lys.add_ionized_restype("LYS");
    ionizables_.push_back(lys);

}

void ReportEffectivePKA::apply(core::pose::Pose & pose) {
    
    if (task_factory() != 0) {
        task_ = task_factory()->create_task_and_apply_taskoperations( pose );
    }
    
	core::conformation::symmetry::SymmetryInfoOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}
	for (core::Size i = 1; i <= pose.total_residue(); ++i) {
		if (symm_info && !symm_info->bb_is_independent( i ) ) continue;
		if (task_factory() != 0 && !task_->residue_task(i).being_designed()) {
            // if task operation is used, only calculated pka for the residue that is being designed
            continue;
        }
        core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
        
        for (Size i_restype =1 ; i_restype<=ionizables_.size(); ++i_restype) {
            if (ionizables_[i_restype].name3() == rsd_i.type().name3()) {
                core::pose::Pose pose_copy(pose);
                
                chemical::ResidueTypeSet const& restype_set( pose_copy.residue(i).residue_type_set() );
                
                core::pose::Pose ref_pose;
                ref_pose.append_residue_by_bond(rsd_i);
                
                Real score_neutral(0.), score_neutral_ref(0.);
                for (Size i_neutral_type = 1; i_neutral_type<=ionizables_[i_restype].neutral_restypes().size(); ++i_neutral_type) {

                    // Create the new residue and replace it
                    conformation::ResidueOP new_res = conformation::ResidueFactory::create_residue(
                                                                                                   restype_set.name_map(ionizables_[i_restype].neutral_restypes()[i_neutral_type] ), pose_copy.residue(i),
                                                                                                   pose_copy.conformation());
                    // Make sure we retain as much info from the previous res as possible
                    conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose_copy.residue(i),
                                                                                     *new_res, pose_copy.conformation() );
                    
                    pose_copy.replace_residue(i, *new_res, false );
                    ref_pose.replace_residue(1, *new_res, false );

                    if (i_neutral_type == 1) {
                        score_neutral = (*scorefxn_)(pose_copy);
                        score_neutral_ref = (*scorefxn_)(ref_pose);

                    }
                    else {
                        Real score = (*scorefxn_)(pose_copy);
                        if (score < score_neutral) {
                            score_neutral = score;
                        }

                        Real score_ref = (*scorefxn_)(ref_pose);
                        if (score_ref < score_neutral_ref) {
                            score_neutral_ref = score_ref;
                        }
                    }
                    
                }
                
                Real score_ionized(0.), score_ionized_ref(0.);
                for (Size i_ionized_type = 1; i_ionized_type<=ionizables_[i_restype].ionized_restypes().size(); ++i_ionized_type) {
                    
                    // Create the new residue and replace it
                    conformation::ResidueOP new_res = conformation::ResidueFactory::create_residue(
                                                                                                   restype_set.name_map(ionizables_[i_restype].ionized_restypes()[i_ionized_type] ), pose_copy.residue(i),
                                                                                                   pose_copy.conformation());
                    // Make sure we retain as much info from the previous res as possible
                    conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose_copy.residue(i),
                                                                                     *new_res, pose_copy.conformation() );
                    
                    pose_copy.replace_residue(i, *new_res, false );
                    ref_pose.replace_residue(1, *new_res, false );
                    
                    if (i_ionized_type == 1) {
                        score_ionized = (*scorefxn_)(pose_copy);
                        score_ionized_ref = (*scorefxn_)(ref_pose);
                        
                    }
                    else {
                        Real score = (*scorefxn_)(pose_copy);
                        if (score < score_ionized) {
                            score_ionized = score;
                        }
                        
                        Real score_ref = (*scorefxn_)(ref_pose);
                        if (score_ref < score_ionized_ref) {
                            score_ionized_ref = score_ref;
                        }
                    }
                }
                
                using namespace ObjexxFCL::format;
                //TR.Debug << "Effective pKa of " << rsd_i.type().name3() << I(4,i) << " is: " << F(8,3, score_ionized) << F(8,3, score_neutral) << F(8,3, score_ionized_ref) << F(8,3, score_neutral_ref) << std::endl;
                
                core::Size output_resi = i;
                if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
                    output_resi = pose.pdb_info()->number( i );
                }

                
                TR << "Effective pKa of " << rsd_i.type().name3() << I(4,output_resi) << " is: " << F(8,3,ionizables_[i_restype].ref_pKa()-ionizables_[i_restype].acid_base_coefficient()*(score_ionized-score_neutral-score_ionized_ref+score_neutral_ref)/1.36) << ", pKa shift is " << F(8,3, -ionizables_[i_restype].acid_base_coefficient()*(score_ionized-score_neutral-score_ionized_ref+score_neutral_ref)/1.36) << std::endl;


                // tag
                core::pose::RemarkInfo remark;
                std::stringstream oss;
                
                protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());

                oss << "Effective pKa of " << rsd_i.type().name3() << I(4,output_resi) << " is: " << F(8,3,ionizables_[i_restype].ref_pKa()-ionizables_[i_restype].acid_base_coefficient()*(score_ionized-score_neutral-score_ionized_ref+score_neutral_ref)/1.36) << ", pKa shift is " << F(8,3, -ionizables_[i_restype].acid_base_coefficient()*(score_ionized-score_neutral-score_ionized_ref+score_neutral_ref)/1.36);
                job->add_string(oss.str());

            }
        }
	}
}

    
/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
ReportEffectivePKA::parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap & datamap,
		Filters_map const &,
		moves::Movers_map const &,
		Pose const &)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

    scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", datamap, "talaris2013" )->clone();
	if( tag->hasOption( "task_operations" ) ) {
		TR << "WARNING: task_operations only active for proteins" << std::endl;
        task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	}
}

protocols::moves::MoverOP
ReportEffectivePKACreator::create_mover() const {
	return protocols::moves::MoverOP( new ReportEffectivePKA );
}

std::string
ReportEffectivePKACreator::keyname() const
{
	return ReportEffectivePKACreator::mover_name();
}

std::string
ReportEffectivePKACreator::mover_name()
{
	return "ReportEffectivePKA";
}

}
}
