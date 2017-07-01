// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelEnzdesCstModule.hh
///
/// @brief this file handles merging constraint defined by enzdes type cstfile
/// @brief and blueprint definition of positions and add them to the pose
/// @author Possu Huang, possu@u.washington.edu, Jan 2010

#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelWorkingSet.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

#include <core/chemical/ResidueType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

using namespace basic::options;

namespace protocols {
namespace forge {
namespace remodel {

static THREAD_LOCAL basic::Tracer TR( "protocols.forge.remodel.RemodelEnzdesCstModule" );

RemodelEnzdesCstModule::RemodelEnzdesCstModule(RemodelData external_data)
: protocols::toolbox::match_enzdes_util::EnzConstraintIO(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ))
{
	scorefxn_ = core::scoring::get_score_function();
	remodel_data_ = external_data;

}

RemodelEnzdesCstModule::~RemodelEnzdesCstModule(){}

void
RemodelEnzdesCstModule::apply(core::pose::Pose & pose)
{
	using namespace protocols::toolbox::match_enzdes_util;

	//set up constraints (read cstfile, do mapping, etc, then add to pose)
	if ( option[OptionKeys::enzdes::cstfile].user() ) {
		enable_constraint_scoreterms(scorefxn_);

		//tmp hack -- from florian for cstcashe observer initialization
		toolbox::match_enzdes_util::EnzdesCstCacheOP cst_cache = toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache();
		if ( !cst_cache ) {
			TR << "cst_cache nonexistant; make new instance." << std::endl;
			toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_cst_cache( toolbox::match_enzdes_util::EnzdesCstCacheOP( new EnzdesCstCache( get_self_ptr(), cst_pairs_.size() ) ) );
			cst_cache = toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache();
		}
		//tmp hack over

		blueprint_cst_definition(pose);

		bool not_packed = true;
		for ( core::Size block = 1 ; block <= cstblocksize_; ++block ) {
			//initialize block internals
			cst_pairs_[ block ]->set_mcfi( this->mcfi_list( block )->mcfi( 1 ) );


			if ( !backbone_only_ && (!cst_pairs_[block]->resA()->is_backbone() || !cst_pairs_[block]->resB()->is_backbone()) ) { //s-b, s-s don't apply if not after design

				TR << "replace sidechain" << std::endl;

				RemodelWorkingSet working_model;
				working_model.workingSetGen(pose, remodel_data_);

				if ( not_packed ) {
					//unfortunately needed this to make the sidechains for fullatom cst
					RemodelDesignMover designMover(remodel_data_, working_model, scorefxn_);
					designMover.set_state("stage");
					designMover.apply(pose);
					not_packed= false;
				}

				TR  << "applying sidechain constraints" << std::endl;
				add_constraints_to_pose_for_block_without_clearing_and_header_processing(pose, scorefxn_, block);
			} else if ( (cst_pairs_[block]->resA()->is_backbone() && cst_pairs_[block]->resB()->is_backbone()) ) { //bb to bb always apply
				TR  << "applying backbone constraints" << std::endl;
				add_constraints_to_pose_for_block_without_clearing_and_header_processing(pose, scorefxn_, block);
			} else {
				TR << "no constraint applied (only sidechain cst was defined in bb stage?) EnzdesCstModule Apply" << std::endl;
			}
		}
		// (*scorefxn_)( pose );
	}

}

void
RemodelEnzdesCstModule::enable_constraint_scoreterms(core::scoring::ScoreFunctionOP scorefxn){
	TR << "turning on constraint weights" << std::endl;
	scorefxn->set_weight(core::scoring::coordinate_constraint, 1.0 );
	scorefxn->set_weight(core::scoring::atom_pair_constraint, 1.0 );
	scorefxn->set_weight(core::scoring::angle_constraint, 1.0 );
	scorefxn->set_weight(core::scoring::dihedral_constraint, 1.0 );
	scorefxn->set_weight(core::scoring::res_type_constraint, 1.0);

}


void
RemodelEnzdesCstModule::blueprint_cst_definition(core::pose::Pose & pose ){

	utility::vector1<core::Size> cstblock;
	utility::vector1<std::string> role;
	utility::vector1<core::Size> position;


	for ( core::Size i = 0, ie = remodel_data_.blueprint.size(); i < ie ; i++ ) {
		if ( remodel_data_.blueprint[i].has_constraints ) {
			for ( std::vector<std::string>::iterator it = remodel_data_.blueprint[i].constraint_definition.begin(), end = remodel_data_.blueprint[i].constraint_definition.end(); it != end; ++it ) {
				//casting, sort of....
				core::Size found_start = (*it).find_first_of("0123456789");
				core::Size found_end = (*it).find_last_of("0123456789");
				std::string buffer((*it).substr(found_start,found_end-found_start+1));
				//TR << "buffer number: " << buffer << std::endl;
				std::istringstream bufferstream(buffer);
				core::Size temp;
				bufferstream >> temp;
				cstblock.push_back(temp);

				core::Size bufferLength = (*it).length();

				role.push_back((*it).substr(bufferLength-1,1));
				position.push_back(remodel_data_.blueprint[i].index);
			}
		}
	}
	TR<< "cstblock.size = " << cstblock.size() << std::endl;
	TR<< "cst_pair_.size = " << cst_pairs_.size() << std::endl;

	runtime_assert ( cstblock.size()/2 == cst_pairs_.size());
	cstblocksize_ = cstblock.size()/2;


	for ( core::Size i = 1; i <= cstblock.size(); ++i ) {
		if ( role[i] == "A" ) {
			TR<< "adding constraint to " << position[i]  << " in " << cstblock[i] << " residue: " << pose.residue_type(position[i]).name3() << std::endl;
			cst_pairs_[cstblock[i]]->nonconst_resA()->set_external_position(position[i]);
			cst_pairs_[cstblock[i]]->nonconst_resA()->find_in_pose_if_missing_from_header( pose );
			// cst_cache->param_cache( cstblock[i] )->template_res_cache( 1 )->add_position_in_pose( position[i] );//new
			// if (!backbone_only_ && !cst_pairs_[cstblock[i]]->resA()->is_backbone()){
			//  TR << "replace Sidechain to " <<  remodel_data_.blueprint[position[i]-1].aminoAcidList[0] << std::endl;//-1 on position array for 0 based switch
			//  core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			//  pose.replace_residue( position[i], *core::conformation::ResidueFactory::create_residue(residue_set->name_map(name_from_aa(remodel_data_.blueprint[position[i]-1].aminoAcidList[0]))),true) ;
			// }
		} else if ( role[i] == "B" ) {
			TR<< "adding constraint to " << position[i] << " in " << cstblock[i] << " residue: " << pose.residue_type(position[i]).name3() << std::endl;
			cst_pairs_[cstblock[i]]->nonconst_resB()->set_external_position(position[i]);
			cst_pairs_[cstblock[i]]->nonconst_resB()->find_in_pose_if_missing_from_header( pose );
			// cst_cache->param_cache( cstblock[i] )->template_res_cache( 2 )->add_position_in_pose( position[i] );//new
			// if (!backbone_only_ && !cst_pairs_[cstblock[i]]->resB()->is_backbone()){
			//  TR << "replace Sidechain to " <<  remodel_data_.blueprint[position[i]-1].aminoAcidList[0] << std::endl;
			//  core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			//  pose.replace_residue( position[i], *core::conformation::ResidueFactory::create_residue(residue_set->name_map(name_from_aa(remodel_data_.blueprint[position[i]-1].aminoAcidList[0]))),true) ;
			// }
		} else {
			TR << "mistake in merging blueprint cst and cst_pairs data" << std::endl;
		}
	}
}

void
RemodelEnzdesCstModule::use_backbone_only_blocks(){
	backbone_only_ = true;
}

void
RemodelEnzdesCstModule::use_all_blocks(){
	backbone_only_ = false;
}


} //namespace protocols
} //namespace forge
} //namespace remodel


