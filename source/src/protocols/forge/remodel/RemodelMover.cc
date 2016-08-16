// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelMover.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Gabriel Rocklin (grocklin@gmail.com) (disulfides)

// unit headers
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>
#include <protocols/forge/remodel/RemodelLigandHandler.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.fwd.hh>
#include <protocols/forge/remodel/RemodelMoverCreator.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>

// project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/io/Remarks.hh>
#include <core/pose/util.hh> // for pdbinfo
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>

#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelMoverCreator.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.fwd.hh>
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelAccumulator.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>

#include <protocols/relax/util.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/DisulfideEntropyFilter.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh> // dihedral constraint
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelAccumulator.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>


// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

// C++ headers
#include <utility>


//#define FILE_DEBUG 1

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;
using namespace core;
using namespace core::scoring;

namespace protocols {
namespace forge {
namespace remodel {


static THREAD_LOCAL basic::Tracer TR( "protocols.forge.remodel.RemodelMover" );

// parser
std::string
RemodelMoverCreator::keyname() const {
	return RemodelMoverCreator::mover_name();
}

protocols::moves::MoverOP
RemodelMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RemodelMover );
}

std::string
RemodelMoverCreator::mover_name() {
	return "RemodelMover";
}

///
/// @brief
/// Default constructor. Checks values of options packing::soft_rep_design, and remodel::dr_cycles
/// Creates and sets the centroid and fullatom score functions.
///
RemodelMover::RemodelMover() :
	Super( "RemodelMover" )
{
	//use_fullmer_ = false;
	//use_sequence_bias_ = false;
	//max_linear_chainbreak_ = 0.15;
	max_linear_chainbreak_ =option[OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak];
	//centroid_loop_mover_str_ = "quick_ccd";
	centroid_loop_mover_str_ = "RemodelLoopMover";
	redesign_loop_neighborhood_ = false;
	dr_cycles_ =option[OptionKeys::remodel::dr_cycles];
	centroid_sfx_ = core::scoring::ScoreFunctionFactory::create_score_function(option[OptionKeys::remodel::cen_sfxn] );
	fullatom_sfx_ = core::scoring::get_score_function();

	register_user_options();

	if ( option[OptionKeys::packing::soft_rep_design] ) {
		TR << "SWITCHING FULLATOM FUNCITON TO SOFT_REP_DESIGN" << std::endl;
		fullatom_sfx_ = ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS );
	}

	blueprint_ = "";
	rosetta_scripts_quick_and_dirty_ = false;
	rosetta_scripts_build_disulfide_ = false;
	rosetta_scripts_fast_disulfide_ = false;
	rosetta_scripts_bypass_fragments_ = false;
	rosetta_scripts_include_current_ds_ = false;
	rosetta_scripts_keep_current_ds_ = false;
	rosetta_scripts_match_rt_limit_ = 1;
	rosetta_scripts_min_disulfides_ = 1;
	rosetta_scripts_max_disulfides_ = 1;
	rosetta_scripts_min_loop_ = 1;
	last_input_pose_ = NULL;
	rosetta_scripts_ = false;
	relax_bb_for_disulf_ = false;
	use_match_rt_ = true;
	use_disulf_fa_score_ = false;
	disulf_fa_max_ = -0.25;

}

///
/// @brief
/// Copy constructor.
///
RemodelMover::RemodelMover( RemodelMover const & rval )
: //utility::pointer::ReferenceCount(),
	Super( rval ),

	//manager_( rval.manager_ ),
	design_info_( rval.design_info_ ),
	//use_fullmer_( rval.use_fullmer_ ),
	//use_sequence_bias_( rval.use_sequence_bias_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	centroid_loop_mover_str_( rval.centroid_loop_mover_str_ ),
	redesign_loop_neighborhood_( rval.redesign_loop_neighborhood_ ),
	//resfile_( rval.resfile_ ),
	dr_cycles_( rval.dr_cycles_ ),
	centroid_sfx_( rval.centroid_sfx_->clone() ),
	fullatom_sfx_( rval.fullatom_sfx_->clone() ),
	blueprint_( rval.blueprint_ ),
	rosetta_scripts_quick_and_dirty_( rval.rosetta_scripts_quick_and_dirty_),
	rosetta_scripts_build_disulfide_( rval.rosetta_scripts_build_disulfide_),
	rosetta_scripts_fast_disulfide_( rval.rosetta_scripts_fast_disulfide_),
	rosetta_scripts_bypass_fragments_( rval.rosetta_scripts_bypass_fragments_),
	rosetta_scripts_match_rt_limit_( rval.rosetta_scripts_match_rt_limit_),
	rosetta_scripts_min_disulfides_( rval.rosetta_scripts_min_disulfides_),
	rosetta_scripts_max_disulfides_( rval.rosetta_scripts_max_disulfides_),
	rosetta_scripts_include_current_ds_( rval.rosetta_scripts_include_current_ds_),
	rosetta_scripts_keep_current_ds_( rval.rosetta_scripts_keep_current_ds_),
	rosetta_scripts_min_loop_( rval.rosetta_scripts_min_loop_),
	rosetta_scripts_( rval.rosetta_scripts_),
	last_input_pose_(rval.last_input_pose_),
	accumulator_(rval.accumulator_),
	relax_bb_for_disulf_(rval.relax_bb_for_disulf_),
	use_match_rt_(rval.use_match_rt_),
	use_disulf_fa_score_(rval.use_disulf_fa_score_),
	disulf_fa_max_(rval.disulf_fa_max_)


{
	if ( rval.vlb_.get() ) {
		vlb_ = VarLengthBuildOP( new VarLengthBuild( *rval.vlb_ ) );
	}
}

///
/// @brief
/// Default destructor. Does this need to free the VarLengthBuild memory?
///
RemodelMover::~RemodelMover() {}

///
/// @brief
/// Checks for presence of any score term weight override options and calls set_weight on the centroid scorefunction.
/// Note: options only get applied to centroid scorefunction - fullatom scorefunction left as is.
///
void RemodelMover::register_user_options() {

	// set optional weights
	if ( option[OptionKeys::remodel::vdw].user() ) {
		centroid_sfx_->set_weight( vdw,option[OptionKeys::remodel::vdw] );
		TR << "USER OVERWRITE VDW: " <<option[OptionKeys::remodel::vdw] << std::endl;
	}

	if ( option[OptionKeys::remodel::rama].user() ) {
		centroid_sfx_->set_weight( rama,option[OptionKeys::remodel::rama] );
		TR << "USER OVERWRITE RAMA: " <<option[OptionKeys::remodel::rama] << std::endl;
	}

	if ( option[OptionKeys::remodel::cbeta].user() ) {
		centroid_sfx_->set_weight( cbeta,option[OptionKeys::remodel::cbeta] );
		TR << "USER OVERWRITE CBETA: " <<option[OptionKeys::remodel::cbeta] << std::endl;
	}

	if ( option[OptionKeys::remodel::cenpack].user() ) {
		centroid_sfx_->set_weight( cenpack,option[OptionKeys::remodel::cenpack] );
		TR << "USER OVERWRITE CENPACK: " <<option[OptionKeys::remodel::cenpack] << std::endl;
	}

	if ( option[OptionKeys::remodel::rg_local].user() ) {
		centroid_sfx_->set_weight( rg_local,option[OptionKeys::remodel::rg_local] );
		TR << "USER OVERWRITE RG_LOCAL: " <<option[OptionKeys::remodel::rg_local] << std::endl;
	}
	if ( option[OptionKeys::remodel::hb_lrbb].user() ) {
		centroid_sfx_->set_weight( hbond_lr_bb,option[OptionKeys::remodel::hb_lrbb] );
		TR << "USER OVERWRITE HB_LRBB: " <<option[OptionKeys::remodel::hb_lrbb] << std::endl;
	}

	if ( option[OptionKeys::remodel::hb_srbb].user() ) {
		centroid_sfx_->set_weight( hbond_sr_bb,option[OptionKeys::remodel::hb_srbb] );
		TR << "USER OVERWRITE HB_SRBB: " <<option[OptionKeys::remodel::hb_srbb] << std::endl;
	}
	if ( option[OptionKeys::remodel::rg].user() ) {
		centroid_sfx_->set_weight( rg,option[OptionKeys::remodel::rg] );
		TR << "USER OVERWRITE RG: " <<option[OptionKeys::remodel::rg] << std::endl;
	}

	if ( option[OptionKeys::remodel::rsigma].user() ) {
		centroid_sfx_->set_weight( rsigma,option[OptionKeys::remodel::rsigma] );
		TR << "USER OVERWRITE RSIGMA: " <<option[OptionKeys::remodel::rsigma] << std::endl;
	}
	if ( option[OptionKeys::remodel::ss_pair].user() ) {
		centroid_sfx_->set_weight( ss_pair,option[OptionKeys::remodel::ss_pair] );
		TR << "USER OVERWRITE SSPAIR: " <<option[OptionKeys::remodel::ss_pair] << std::endl;
	}

}


/// @brief clone for parser
RemodelMover::MoverOP RemodelMover::clone() const {
	return RemodelMover::MoverOP( new RemodelMover( *this ) );
}


/// @brief fresh instance for parser
RemodelMover::MoverOP RemodelMover::fresh_instance() const {
	return RemodelMover::MoverOP( new RemodelMover() );
}


/// @brief clone this object
RemodelMover::MoverOP RemodelMover::clone() {
	return RemodelMover::MoverOP( new RemodelMover( *this ) );
}


/// @brief create this type of object
RemodelMover::MoverOP RemodelMover::fresh_instance() {
	return RemodelMover::MoverOP( new RemodelMover() );
}


/// @brief the centroid level score function, default "remodel_cen"
ScoreFunction const & RemodelMover::centroid_scorefunction() const {
	return *centroid_sfx_;
}


/// @brief the full-atom level score function
ScoreFunction const & RemodelMover::fullatom_scorefunction() const {
	return *fullatom_sfx_;
}


/// @brief set the centroid level score function
void RemodelMover::centroid_scorefunction( ScoreFunction const & sfx ) {
	centroid_sfx_ = sfx.clone();
}


/// @brief set the centroid level score function
void RemodelMover::centroid_scorefunction( ScoreFunctionOP const & sfx ) {
	centroid_sfx_ = sfx->clone();
}


/// @brief set the full-atom level score function
void RemodelMover::fullatom_scorefunction( ScoreFunction const & sfx ) {
	fullatom_sfx_ = sfx.clone();
}


/// @brief set the full-atom level score function
void RemodelMover::fullatom_scorefunction( ScoreFunctionOP const & sfx ) {
	fullatom_sfx_ = sfx->clone();
}

///Function for recursively creating multiple disulfides
utility::vector1< utility::vector1< std::pair<Size,Size> > > recursive_multiple_disulfide_former (
	utility::vector1< std::pair<Size,Size> >  & disulfides_formed,
	utility::vector1< std::pair<Size,Size> >  & disulfides_possible,
	Size const & max_disulfides) {

	utility::vector1< utility::vector1< std::pair<Size,Size> > > final_configurations;

	if ( disulfides_formed.size() < max_disulfides ) {

		//select one primary new disulfide to be added
		for ( utility::vector1< std::pair< Size, Size > >::iterator new_disulfide = disulfides_possible.begin(), end = disulfides_possible.end();
				new_disulfide != end;
				++new_disulfide ) {

			//add the configuration with the new disulfide
			utility::vector1< std::pair<Size,Size> > new_disulfides_formed = disulfides_formed;
			new_disulfides_formed.push_back(*new_disulfide);
			final_configurations.push_back(new_disulfides_formed);

			//storage for possible disulfides which do not clash with the new one we have chosen
			utility::vector1< std::pair<Size,Size> >  disulfides_to_be_added;

			//identify new secondary disulfides which do not clash with the primary
			utility::vector1< std::pair< Size, Size > >::iterator potential_further_disulfide = new_disulfide;
			for ( ++potential_further_disulfide; potential_further_disulfide != end;
					++potential_further_disulfide ) {

				if ( (*potential_further_disulfide).first != (*new_disulfide).first &&
						(*potential_further_disulfide).second != (*new_disulfide).first &&
						(*potential_further_disulfide).first != (*new_disulfide).second &&
						(*potential_further_disulfide).second != (*new_disulfide).second ) {

					disulfides_to_be_added.push_back(*potential_further_disulfide);
				}

			} //end identification of possible secondary disulfides

			if ( disulfides_to_be_added.size() > 0 ) { //Add new disulfides if new ones can be added

				utility::vector1< utility::vector1< std::pair<Size,Size> > > new_disulfide_configurations =
					recursive_multiple_disulfide_former(new_disulfides_formed, disulfides_to_be_added, max_disulfides);

				for ( utility::vector1< utility::vector1< std::pair<Size,Size> > >::iterator new_configuration = new_disulfide_configurations.begin(), ndcend = new_disulfide_configurations.end();
						new_configuration != ndcend;
						++new_configuration ) {
					final_configurations.push_back(*new_configuration);

				}
			} //end adding new disulfides recursively AFTER the first selected new one
		} //end addition of primary disulfide + all secondary possibilities
	} else {
		//final_configurations.push_back(disulfides_formed); //store the current configuration if we have reached max length
	}
	return final_configurations;
}

///
/// @brief
/// Apply method for Mover.
/// Checks the values of the following options
/// -remodel::checkpoint
/// -remodel::domainFusion::insert_segment_from_pdb
/// -remodel::bypass_fragments,
/// -remodel::num_trajectory,
/// -remodel::repeat_structure
/// -remodel::build_disulf
/// -remodel::quick_and_dirty
/// -remodel::use_pose_relax
/// -remodel::run_confirmation
/// -enzdes::cstfile
/// -symmetry::symmetry_definition
/// -run::show_simulation_in_pymol
///

bool RemodelMover::SamePose( Pose const & pose1, Pose const & pose2) {
	if ( pose1.total_residue() != pose2.total_residue() ) {
		return false;
	}

	for ( core::Size i = 1; i <= pose1.total_residue(); ++i ) {
		if ( std::abs(pose1.phi(i) - pose2.phi(i)) > 0.0001 ) {
			return false;
		}
		if ( std::abs(pose1.psi(i) - pose2.psi(i)) > 0.0001 ) {
			return false;
		}
		if ( std::abs(pose1.omega(i) - pose2.omega(i)) > 0.0001 ) {
			return false;
		}
	}
	return true;
}


void RemodelMover::apply( Pose & pose ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace protocols;
	using namespace chemical;

	using core::pose::metrics::CalculatorFactory;
	using core::pose::metrics::PoseMetricCalculatorOP;
	using core::pose::PDBInfoOP;
	using core::pack::task::operation::TaskOperationCOP;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_RETRY;
	using protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "Remodel" );
#endif
	if ( !last_input_pose_  || !SamePose(*last_input_pose_, pose) || accumulator_.size() == 0 ) {
		last_input_pose_ = core::pose::PoseOP( new core::pose::Pose(pose) );


		TR << "apply(): entered RemodelMover apply(). pose.total_residue(): " << pose.total_residue() << std::endl;

		// store the starting pose for KIC confirmation RMSD calculation
		native_pose_ = pose;
		// assign secondary structure
		scoring::dssp::Dssp dssp( pose );

		TR << pose.total_residue() << std::endl;
		ObjexxFCL::FArray1D_char dsspSS( pose.total_residue() );
		dssp.dssp_reduced(dsspSS);

		TR << "apply(): input PDB dssp assignment: (based on start structure)" << std::endl;
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			TR << dsspSS(i);
		}
		TR << std::endl;

		forge::remodel::RemodelData remodel_data;
		forge::remodel::RemodelWorkingSet working_model;

		// read blueprint file
		// logic: if provided at command line, use that.
		//        if not at command line but given in rosetta scripts, use that.
		//        if not given at command line or rosetta scripts, make ad hoc blueprint.
		TR << "apply(): reading blueprint file " << std::endl;
		if ( option[OptionKeys::remodel::blueprint].user() ) {
			blueprint_ =option[OptionKeys::remodel::blueprint]();
			remodel_data.getLoopsToBuildFromFile(blueprint_);
		} else {
			if ( blueprint_ != "" ) {
				remodel_data.getLoopsToBuildFromFile(blueprint_);
			} else {
				//generate blueprint on the fly for disulfide remodeling
				TR << "Generating ad hoc blueprint" << std::endl;
				std::stringstream ad_hoc_blueprint;
				//ad_hoc_blueprint = "";
				for ( Size i = 1; i <= pose.total_residue(); ++i ) {
					//number
					ad_hoc_blueprint << i;
					//residue
					ad_hoc_blueprint << "  V";
					//helix
					ad_hoc_blueprint << "   " << dsspSS(i) << std::endl;
				}
				TR << ad_hoc_blueprint.str() << std::endl;
				remodel_data.getLoopsToBuildFromBlueprint(ad_hoc_blueprint.str());
			}
		}

		remodel_data.updateWithDsspAssignment( dsspSS );
		//dssp.insert_ss_into_pose( pose ); This put the assigned sec structure into the pose, as opposed to the actual SS of the pose. Thus eliminated
		// process domain insertion option
		if ( option[OptionKeys::remodel::domainFusion::insert_segment_from_pdb].user() ) {
			TR << "apply(): INSERT SEGMENT FROM PDB" << std::endl;
			remodel_data.collectInsertionPose();
			// remodel_data will have several class member variables updated with the insertion information at this point
		}

		// create a scorefunction and score the pose first
		scoring::ScoreFunctionOP sfx = scoring::get_score_function();
		(*sfx)( pose );

		working_model.workingSetGen( pose, remodel_data );

		remodel_data_ = remodel_data; // will use the movemap for natro definition later
		working_model_ = working_model;

		//TR << "After working_model_" << std::endl;

		// test PyMol viewer
		/*
		NOTE:: If pose size changes such as in repeat_mover this can crash remodel
		if ( option[OptionKeys::run::show_simulation_in_pymol] ) {
		moves::AddPyMolObserver( pose, false, core::Real( 0.50 ) );
		}
		*/


		/*
		// DEBUG
		std::set<core::Size> up = working_model.manager.undefined_positions();
		for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++) {
		TR << *i << std::endl;
		}
		std::set<core::Size> uup = working_model.manager.union_of_intervals_containing_undefined_positions();
		for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP" <<  std::endl;
		}
		*/

		if ( option[OptionKeys::remodel::repeat_structure].user() ) {
			//for cases involve jxn, need to make pose longer so manager won't complain
			//about missing residues

			//this is pre modify, so simply extend to 2x blueprint length,
			//with extensions, the pose will go beyond the correct length.  need to fix
			//that after modify. Residues beyond first copy+ jxn doesn't really matter
			if ( pose.total_residue() < 2*remodel_data.sequence.length() ) { //just making sure it's shorter before grow, input pose can be longer
				Size len_diff = (2*remodel_data_.sequence.length()) - pose.total_residue();
				// append a tail of the same length
				String build_aa_type =option[OptionKeys::remodel::generic_aa]; //defaults to VAL
				if ( build_aa_type.size() == 1 ) {
					char build_aa_oneLetter= build_aa_type[0];
					build_aa_type = name_from_aa(aa_from_oneletter_code(build_aa_oneLetter));
				}
				for ( Size i = 1; i<= len_diff; ++i ) {
					core::chemical::ResidueTypeSetCOP rsd_set = (pose.residue(1).residue_type_set());
					core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_set->name_map(build_aa_type) ) );
					pose.conformation().safely_append_polymer_residue_after_seqpos(* new_rsd,pose.total_residue(), true);
					pose.conformation().insert_ideal_geometry_at_polymer_bond(pose.total_residue()-1);
					pose.set_omega(pose.total_residue()-1,180);
				}
			}
		}
		//
		// Pose testArc;
		// testArc = pose;
		if ( working_model.manager.size()!= 0 ) {
			if ( !option[OptionKeys::remodel::bypass_fragments] && !rosetta_scripts_bypass_fragments_ && !rosetta_scripts_fast_disulfide_ ) {
				working_model.manager.modify(pose);
			} else {
				working_model.manager.dummy_modify(pose.total_residue());

			}
			scoring::dssp::Dssp dssp( pose );
			dssp.insert_ss_into_pose( pose );

			// protocols::forge::methods::restore_residues(working_model.manager.original2modified(), testArc, pose);
			// pose.dump_pdb("testArcRestore.pdb");
			//testArc=pose;
			manager_ = working_model.manager;
			core::pose::renumber_pdbinfo_based_on_conf_chains(
				pose,
				true ,  // fix chain
				true, // start_from_existing_numbering
				false, // keep_insertion_code
				false // rotate_chain_id
			);
		}
		//finally recheck length to ensure blueprint compliance
		if ( option[OptionKeys::remodel::repeat_structure].user() ) {
			Size max_pdb_index = remodel_data_.blueprint.size()*2;
			while ( pose.total_residue() > max_pdb_index ) {
				pose.conformation().delete_residue_slow(pose.total_residue());
				pose.pdb_info()->obsolete(true); //note the previous line was having issues with the pymol observer. You may also want to add -show_simulation_in_pymol 0 to your flags.
			}
			if ( pose.total_residue() < (remodel_data_.sequence.length()*2) ) {
				Size len_diff = (2*remodel_data_.sequence.length()) - pose.total_residue();
				// append a tail of the same length
				for ( Size i = 1; i<= len_diff; ++i ) {
					core::chemical::ResidueTypeSetCOP rsd_set = (pose.residue(1).residue_type_set());
					core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_set->name_map("ALA") ) );
					pose.conformation().safely_append_polymer_residue_after_seqpos(* new_rsd,pose.total_residue(), true);
					pose.pdb_info()->obsolete(true);
					pose.conformation().insert_ideal_geometry_at_polymer_bond(pose.total_residue()-1);
					pose.set_omega(pose.total_residue()-1,180);
				}
			}

		}


		/*
		up = working_model.manager.undefined_positions();
		for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++ ) {
		TR << *i << std::endl;
		}
		uup = working_model.manager.union_of_intervals_containing_undefined_positions();
		for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP2" <<  std::endl;
		}
		*/

		//manager_.dummy_modify(testArc.n_residue());
		//core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true);
		//core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD, true);
		//protocols::forge::methods::restore_residues(manager_.original2modified(), testArc, pose);
		//pose.dump_pdb("testArcRestore2.pdb");
		//protocols::forge::methods::restore_residues(manager_.original2modified(), testArc, pose);
		//pose.update_residue_neighbors();
		//pose.dump_pdb("testArcRestore3.pdb");
		//testArc.dump_pdb("testArcRestoreSrc3.pdb");
		//protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD);
		//protocols::simple_moves::ReturnSidechainMover recover_sidechains( testArc);
		//to_all_atom.apply(pose);
		//recover_sidechains.apply(pose);
		//pose.dump_pdb("MoverREstore.pdb");

		// initialize symmetry

		// only symmetrize here if not in the repeat structure mode. for repeats, stay monomer until repeat generation.
		if ( option[OptionKeys::symmetry::symmetry_definition].user() && !option[OptionKeys::remodel::repeat_structure].user() )  {
			simple_moves::symmetry::SetupForSymmetryMover pre_mover;
			pre_mover.apply( pose );
			// Remodel assumes chain ID is ' '
			//pose::PDBInfoOP pdb_info ( pose.pdb_info() );
			//for ( Size i=1; i<= pdb_info->nres(); ++i ){
			// pdb_info->chain(i,' ');
			//
			//pose.pdb_info( pdb_info );
		}


		Size i =option[OptionKeys::remodel::num_trajectory];
		//if invoked from rosetta_scripts and num_trajectories not specified, change default to 1
		if ( !option[OptionKeys::remodel::num_trajectory].user() && rosetta_scripts_ ) {
			i = 1;
		}
		Size num_traj = i; // need this for checkpointing math
		Size prev_checkpoint = 0;

		//accumulator_ is now owned by RemodelMover and not by apply()
		//forge::remodel::RemodelAccumulator accumulator_( working_model );

		if ( option[OptionKeys::remodel::checkpoint] ) {
			prev_checkpoint = accumulator_.recover_checkpoint();
			if ( prev_checkpoint >= i ) {
				i = 0;
			} else {
				i = i - prev_checkpoint;
			}
		}

		if ( working_model.manager.size() != 0 ) {
			// setup calculators
			pose::metrics::CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
			pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc_name(),
				PoseMetricCalculatorOP( new toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() ) ) );
		}

		/*
		up = working_model.manager.undefined_positions();
		for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++){
		TR << *i << std::endl;
		}
		uup = working_model.manager.union_of_intervals_containing_undefined_positions();
		for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP2" <<  std::endl;
		}
		*/

		if ( option[OptionKeys::remodel::repeat_structure].user() ) {

			// turning on the res_type_linking constraint weight for designs
			fullatom_sfx_->set_weight( scoring::atom_pair_constraint, 1.0 );
			fullatom_sfx_->set_weight( scoring::coordinate_constraint, 1.0 );
			fullatom_sfx_->set_weight( scoring::res_type_linking_constraint, 0.3 );
			fullatom_sfx_->set_weight( scoring::res_type_constraint, 1.0 );

			//fullatom_sfx_->set_weight( core::scoring::fa_dun, 0 );
			//fullatom_sfx_->set_weight( core::scoring::fa_sol, 0 );
			//fullatom_sfx_->set_weight( core::scoring::fa_pair, 0 );
			//fullatom_sfx_->set_weight( core::scoring::hbond_sc, 0 );
			//fullatom_sfx_->set_weight( core::scoring::hbond_sc, 0 );
			//fullatom_sfx_->set_weight( core::scoring::fa_intra_rep, 0 );
			//fullatom_sfx_->set_weight( core::scoring::rama, 0 );
			//fullatom_sfx_->set_weight( core::scoring::hbond_bb_sc, 0 );
			//fullatom_sfx_->set_weight( core::scoring::p_aa_pp, 0 );
			//fullatom_sfx_->set_weight( core::scoring::ref, 10 );
		}

		// initializes a RemodelDesignMover which will be used in the loop below
		forge::remodel::RemodelDesignMover designMover( remodel_data, working_model, fullatom_sfx_ );

		Size no_attempts_at_centroid_build = 0;
		//Size no_attempts_to_make_at_centroid_build = 10;

		//bool quick_mode =option[OptionKeys::remodel::quick_and_dirty]();

		Size repeat_number =option[OptionKeys::remodel::repeat_structure];
		//rerooting tree
		if ( option[OptionKeys::remodel::repeat_structure].user() && pose.total_residue() == remodel_data.blueprint.size()*repeat_number ) {
			core::kinematics::FoldTree f = pose.fold_tree();
			f.reorder(working_model.safe_root_);
			pose.fold_tree(f);
			TR << "rerooting tree: " << pose.fold_tree() << std::endl;
		}

		while ( i > 0 ) {

			// cache the modified pose first for REPEAT
			Pose cached_modified_pose( pose );

			// do centroid build
			TR << std::endl << "apply(): BUILD CYCLE REMAINING " << i << std::endl;
			kinematics::FoldTree originalTree = pose.fold_tree();
			TR << "ORIGINAL TREE: " << pose.fold_tree() << std::endl;
			if ( working_model.manager.size() != 0 ) {
				if ( !centroid_build( pose, working_model.manager ) ) { // build failed
					no_attempts_at_centroid_build++;
					TR << "apply(): number of attempts at loop closure made: " << no_attempts_at_centroid_build << std::endl;
					set_last_move_status( FAIL_RETRY );

					/*if ( !quick_mode ) {
					continue;
					} else {
					if ( no_attempts_at_centroid_build >= no_attempts_to_make_at_centroid_build ) {
					TR << "apply(): number of attempts at loop closure exceeded limit of " << no_attempts_to_make_at_centroid_build << ". quitting." << std::endl;
					set_last_move_status( FAIL_DO_NOT_RETRY );
					return;
					} else {
					// try again. omitting this continue causes the protocol to give up after one failed iteration.
					continue;
					}
					}*/
					i--;
					continue;
					//return;
				}
			}

			if ( option[OptionKeys::remodel::repeat_structure].user() ) {
				// should fold this pose to match just the first segment of a repeat, and that will be used for next round of building
				add_lower_terminus_type_to_pose_residue(pose,1);
				for ( Size res = 1; res <= cached_modified_pose.n_residue(); res++ ) {
					cached_modified_pose.set_phi( res, pose.phi(res) );
					cached_modified_pose.set_psi( res, pose.psi(res) );
					cached_modified_pose.set_omega( res, pose.omega(res) );
					cached_modified_pose.set_secstruct(res, pose.secstruct(res));
					ResidueType const & rsd_type(pose.residue_type(res));
					replace_pose_residue_copying_existing_coordinates(cached_modified_pose,res,rsd_type);
				}
			}
			core::pose::renumber_pdbinfo_based_on_conf_chains(
				pose,
				true ,  // fix chain
				true, // start_from_existing_numbering
				false, // keep_insertion_code
				false // rotate_chain_id
			);

			//test
			//pose.dump_pdb("check.pdb");
			/*
			//extract the constraint currently in Pose for later Recycling
			ConstraintSetOP cst_set_post_built;
			if (option[OptionKeys::remodel::repeat_structure].user() ) {

			// at this stage it should hold generic cstfile and res_type_linking
			// constraints
			cst_set_post_built = new ConstraintSet( *pose.constraint_set());
			}
			*/

			designMover.set_state("stage");

			// handle constraints as soon as centroid is done.  If applying sidechain
			// constraints, replace residue to the right AA right away.
			if ( option[OptionKeys::enzdes::cstfile].user() ) {
				TR << "apply(): constraint file found on command line. updating score functions to include constraint terms." << std::endl;

				forge::remodel::RemodelEnzdesCstModuleOP cstOP( new forge::remodel::RemodelEnzdesCstModule(remodel_data) );

				// RemodelEnzdesCstModule cst(remodel_data);
				//safety
				pose.remove_constraints();
				//wipe out cst_cache
				toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> set_cst_cache( NULL );
				//wipe out observer too
				pose.observer_cache().set( pose::datacache::CacheableObserverType::ENZDES_OBSERVER, NULL , false);

				//cstOP->remove_constraints_from_pose(pose,true /*keep covalent*/, false /*fail if missing*/);

				cstOP->use_all_blocks();
				TR << "apply(): calling RemodelEnzdesCstModule apply function." << std::endl;
				cstOP->apply(pose);
				cstOP->enable_constraint_scoreterms(fullatom_sfx_);
				designMover.scorefunction(fullatom_sfx_);
			}
			if ( option[OptionKeys::constraints::cst_file].user() ) {
				//safety
				pose.remove_constraints();

				protocols::simple_moves::ConstraintSetMoverOP constraint( new protocols::simple_moves::ConstraintSetMover() );
				constraint->apply( pose );

				fullatom_sfx_->set_weight(core::scoring::atom_pair_constraint, 1.0);
				fullatom_sfx_->set_weight(core::scoring::coordinate_constraint, 1.0);
				fullatom_sfx_->set_weight(core::scoring::dihedral_constraint, 10.0);
				designMover.scorefunction(fullatom_sfx_);
			}

			if ( option[OptionKeys::remodel::build_disulf]() || rosetta_scripts_build_disulfide_ || rosetta_scripts_fast_disulfide_ ) {

				utility::vector1<std::pair <Size, Size> > disulf_partners;
				bool disulfPass = false;

				core::Energy match_rt_limit = rosetta_scripts_match_rt_limit_;
				if ( option[OptionKeys::remodel::build_disulf]() ) {
					match_rt_limit = option[OptionKeys::remodel::match_rt_limit];
				}

				disulfPass = designMover.find_disulfides_in_the_neighborhood( pose, disulf_partners, match_rt_limit, rosetta_scripts_min_loop_, rosetta_scripts_include_current_ds_, rosetta_scripts_keep_current_ds_,
					relax_bb_for_disulf_, use_match_rt_, use_disulf_fa_score_, disulf_fa_max_);
				if ( disulfPass != true ) {
					i--; //for now control disulf with num_trajectory flag, too.
					continue;
				}

				//Use the recursive multiple disulfide former
				utility::vector1< std::pair<Size,Size> > empty_disulfide_list;
				utility::vector1< utility::vector1< std::pair<Size,Size> > > disulfide_configurations =
					recursive_multiple_disulfide_former(empty_disulfide_list, disulf_partners, rosetta_scripts_max_disulfides_);

				//iterate over disulfide configurations instead of over possible disulfides
				for ( utility::vector1< utility::vector1< std::pair<Size,Size> > >::iterator ds_config = disulfide_configurations.begin();
						ds_config != disulfide_configurations.end();
						++ds_config ) {
					if ( (*ds_config).size() >= rosetta_scripts_min_disulfides_ && (*ds_config).size() <= rosetta_scripts_max_disulfides_ ) {

						//form all the disulfides in the disulfide configuration
						TR << "Building disulfide configuration ";
						for ( utility::vector1< std::pair<Size,Size> >::iterator my_ds = (*ds_config).begin(); my_ds != (*ds_config).end(); ++my_ds ) {
							TR << (*my_ds).first << "-" << (*my_ds).second << " ";
						}
						TR << std::endl;

						Pose disulf_copy_pose = pose;

						kinematics::MoveMapOP combined_mm( new kinematics::MoveMap );

						combined_mm->import( remodel_data_.natro_movemap_ );
						combined_mm->import( manager_.movemap() );

						if ( !rosetta_scripts_fast_disulfide_ ) {
							//original way of doing things, design is included through designMover.apply
							designMover.make_disulfide( disulf_copy_pose, *ds_config, combined_mm );
							designMover.apply( disulf_copy_pose );
						} else {
							//fast way of doing things, assumes design later in a rosetta script, no need to do any now.
							//designMover.apply is not called.
							//kinematics::MoveMapOP freeze_bb_combined_mm = new kinematics::MoveMap(*combined_mm);
							//(*freeze_bb_combined_mm).set_bb(false);
							designMover.make_disulfide_fast( disulf_copy_pose, *ds_config);//, freeze_bb_combined_mm );
							//designMover.apply( disulf_copy_pose );
						}

						// for now, accept all disulf build, as it is hard enough to do already.  Accept instead of cst filter?
						// accumulator_.apply(disulf_copy_pose);
						if ( option[OptionKeys::enzdes::cstfile].user() ) {
							simple_filters::ScoreTypeFilter const pose_constraint( fullatom_sfx_, atom_pair_constraint, 10 );
							bool CScore(pose_constraint.apply( disulf_copy_pose ));
							if ( !CScore ) {  // if didn't pass, rebuild
								continue;
							} else {
								accumulator_.apply(disulf_copy_pose);
							}

						} else {
							accumulator_.apply(disulf_copy_pose);
						}
					}
				}

			} else {
				// option build_disulf not specified...
				//if (option[OptionKeys::remodel::repeat_structure].user() ||option[OptionKeys::remodel::cen_minimize] ) {
				if ( option[OptionKeys::remodel::cen_minimize] ) {

					//cache current foldTree;
					kinematics::FoldTree cenFT = pose.fold_tree();
					pose.fold_tree(originalTree);

					kinematics::MoveMapOP cmmop( new kinematics::MoveMap );
					//pose.dump_pdb("pretest.pdb");

					cmmop->import( remodel_data_.natro_movemap_ );
					cmmop->import( manager_.movemap() );

					for ( Size i = 1; i<= pose.total_residue(); ++i ) {
						std::cout << "bb at " << i << " " << cmmop->get_bb(i) << std::endl;
					}

					//adding angles and bonds dof
					// cmmop->set(core::id::THETA, true);
					// cmmop->set(core::id::D, true);

					for ( Size i = 1; i <= pose.n_residue(); i++ ) {
						for ( Size j = 1; j <= pose.residue(i).nheavyatoms(); j++ ) {
							if ( cmmop->get_bb(i) == 1 ) {
								cmmop->set(core::id::DOF_ID(core::id::AtomID(j,i),core::id::THETA),true);
								cmmop->set(core::id::DOF_ID(core::id::AtomID(j,i),core::id::D),true);
							}
						}
					}


					//scoring::ScoreFunctionOP cen_min_sfxn = scoring::ScoreFunctionFactory::create_score_function("score4_smooth");

					TR << "centroid minimizing" << std::endl;
					pose::Pose archived_pose = pose;

					// flip residue type set for centroid minimize
					util::switch_to_residue_type_set( pose, chemical::CENTROID, true );

					protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;

					if ( option[OptionKeys::remodel::repeat_structure].user() ) {
						//Dihedral (NCS) Constraints, need to be updated each mutation cycle for sidechain symmetry

						Size repeat_number =option[OptionKeys::remodel::repeat_structure];
						Size segment_length = (pose.n_residue())/repeat_number;


						for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
							std::stringstream templateRangeSS;
							templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
							std::stringstream targetSS;
							targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
							TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
							setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
						}

						for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
							std::stringstream templateRangeSS;
							templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
							std::stringstream targetSS;
							targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
							TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
							setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
						}


						std::stringstream templateRangeSS;
						//take care of the terminal repeat, since the numbers are offset.
						templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
						std::stringstream targetSS;
						targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
						TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
						setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
						setup_ncs.apply(pose);

					}

					//sfx->show(TR, pose);
					//TR << std::endl;

					centroid_sfx_->set_weight( core::scoring::atom_pair_constraint, 1.0);
					centroid_sfx_->set_weight( core::scoring::coordinate_constraint, 1.0);
					centroid_sfx_->set_weight(core::scoring::dihedral_constraint, 10.0 );
					//enable cartesian bond terms
					centroid_sfx_->set_weight(core::scoring::cart_bonded_angle,  5.0 );
					centroid_sfx_->set_weight(core::scoring::cart_bonded_length,  1.0 );
					centroid_sfx_->set_weight(core::scoring::cart_bonded_torsion,  5.0 );
					centroid_sfx_->set_weight(core::scoring::omega, 0.2 );

					//only use smooth hb if either of the term is used in centroid build level
					if ( centroid_sfx_->get_weight(core::scoring::hbond_lr_bb) > 0 || centroid_sfx_->get_weight(core::scoring::hbond_sr_bb) > 0 ) {
						centroid_sfx_->set_weight( core::scoring::cen_hb, 2.0);
					}
					/*
					if (centroid_sfx_->get_weight(core::scoring::env) > 0 ){
					centroid_sfx_->set_weight( core::scoring::cen_env_smooth, centroid_sfx_->get_weight(core::scoring::env));
					centroid_sfx_->set_weight( core::scoring::env, 0.0);
					}*/

					//simple_moves::MinMoverOP minMover = new simple_moves::MinMover( cmmop , centroid_sfx_, "lbfgs_armijo", 0.01, true);
					simple_moves::MinMoverOP minMover( new simple_moves::MinMover( cmmop , centroid_sfx_, "lbfgs_armijo", 0.01, true) );
					TR << "cen_minimize pose foldtree: " << pose.fold_tree() << std::endl;
					minMover->apply(pose);

					//reset cen_hb to 0
					centroid_sfx_->set_weight( core::scoring::cen_hb, 0.0);
					//switch back the foldtree
					//pose.fold_tree(cenFT);

					// flip residue type set back, for repeat builds, currently don't do
					// restore_sidechain, as they should all be redesigned.  MAY NEED TO
					// CHANGE
					util::switch_to_residue_type_set( pose, chemical::FA_STANDARD, true );
					//forge::methods::restore_residues( manager_.original2modified(), archived_pose , pose );
					//pose.dump_pdb("test.pdb");

				}




				TR << "apply(): calling RemodelDesignMover apply function." << std::endl;

				//****Previously this loop didn't maintain the pose correctly which resulted in a difficult to track down seg fault.  Fixed, but it's questionable weather you would want to filter poses based on constraints.
				if ( option[OptionKeys::enzdes::cstfile].user() ||
						option[OptionKeys::constraints::cst_file].user()
						) {
					simple_filters::ScoreTypeFilter const  pose_constraint( fullatom_sfx_, atom_pair_constraint,option[OptionKeys::remodel::cstfilter] );
					bool CScore(pose_constraint.apply( pose ));
					if ( !CScore ) {  // if didn't pass, rebuild
						TR << "built model did not pass constraints test." << std::endl;
						if ( option[OptionKeys::remodel::repeat_structure].user() ) {
							//reset the pose to monomer
							pose = cached_modified_pose;
						} else {
							pose.fold_tree(originalTree);
						}
						i--; //decrement the number of tries, so it won't run into an infinite loop if -bypass_fragments
						continue;
					} else {
						designMover.apply(pose);
						accumulator_.apply(pose);
					}
				} else {
					designMover.apply(pose);
					accumulator_.apply(pose);
				}

				//put this ligand test here for now.
				if ( option[OptionKeys::in::file::extra_res_fa].user() && option[OptionKeys::remodel::move_ligand].user() ) {
					TR << "ligand handler" << std::endl;
					RemodelLigandHandler ligand_handler;
					ligand_handler.minimize(pose);
					pose.dump_pdb("ligand1.pdb");
				}
				//*****
			}
			if ( option[OptionKeys::remodel::checkpoint] ) {
				// debug:
				TR << "writing chkpnt at step " << num_traj-i+prev_checkpoint << std::endl;
				accumulator_.write_checkpoint(num_traj-i-prev_checkpoint);
			}
			// restore foldtree
			if ( i > 1 ) { //messes up rosetta_scripts if done on last loop.
				if ( option[OptionKeys::remodel::repeat_structure].user() ) {
					//reset the pose to monomer
					pose = cached_modified_pose;
				} else {
					pose.fold_tree(originalTree);
				}
			}
			i--; // 'i' is the number of remaining trajectories

		}
		/* DONT USE THIS FOR NOW...
		if (get_last_move_status() == FAIL_RETRY){
		return;
		}
		*/
		// take the lowest member and the cluster center
		// accumulator_.shrink_cluster();

		//move the accumulator results into results
		std::vector< pose::PoseOP > results;

		if ( accumulator_.cluster_switch() ) {
			results = accumulator_.clustered_best_poses();
			//results = accumulator_.clustered_top_poses(op_remodel_collect_clustered_top_);
		} else {
			results = accumulator_.contents_in_pose_store();
		}

		//reset the accumulator
		accumulator_.clear();

		// seriously refine the poses
		Size filecount = 1;

		TR << "clustered poses count: " << results.size() << std::endl;

		for ( std::vector< pose::PoseOP >::iterator it = results.begin(), end= results.end(); it!= end; ++it ) {

			TR << "clustered poses count: " << results.size() << std::endl;

			bool bypass_refinement = option[OptionKeys::remodel::quick_and_dirty];
			if ( rosetta_scripts_quick_and_dirty_ || rosetta_scripts_fast_disulfide_ ) {
				bypass_refinement = true;
			}

			if ( working_model.manager.size() == 0 ) {
				bypass_refinement = true;
			}

			if ( !bypass_refinement ) {

				//std::stringstream SS1;
				//SS1 << "pre-ref_" << filecount << ".pdb";
				//(*(*it)).dump_scored_pdb(SS1.str(), *fullatom_sfx_);

				TR << "aggressively refine" << std::endl;
				if ( option[OptionKeys::remodel::use_pose_relax] ) {
					if ( !design_refine_seq_relax( *(*it), designMover ) ) {
						TR << "WARNING: DESIGN REFINE SEQ RELAX FAILED!! (one should never see this)" << std::endl;
						continue;
					}
				} else if ( option[OptionKeys::remodel::use_cart_relax] ) {
					if ( !design_refine_cart_relax(*(*it), designMover) ) {
						TR << "WARNING: CARTESIAN MIN FAILED!! (one should never see this)" << std::endl;
						continue;
					}
				} else {
					if ( ! design_refine(*(*it), designMover) ) {
						TR << "WARNING: DESIGN REFINE FAILED TO CLOSE STRUCTURE!!" << std::endl;
						continue;
					}
				}
			} else { // simple design
				if ( option[OptionKeys::remodel::check_scored_centroid]() ) {
					std::stringstream SS;
					std::string prefix = option[OptionKeys::out::prefix];
					if ( !prefix.empty() ) {
						SS << prefix << "_" << filecount << "_cen.pdb";
					} else {
						SS << filecount << "_cen.pdb";
					}
					core::util::switch_to_residue_type_set( *(*it), core::chemical::CENTROID, true);
					(*(*it)).dump_scored_pdb(SS.str(), *centroid_sfx_);
				}

				designMover.set_state("finish");
				if ( !rosetta_scripts_fast_disulfide_ ) {
					designMover.apply(*(*it));
				}


			}

			if ( option[OptionKeys::remodel::run_confirmation]() ) {
				if ( !confirm_sequence(*(*it)) ) {
					TR << "WARNING: STRUCTURE DID NOT PASS KIC CONFIRMATION!!" << std::endl;
					continue;
				}
			}

			//TR << "CARMSD" << core::scoring::CA_rmsd( native_pose_, *(*it) ) << std::endl;
			//setPoseExtraScores( *(*it), "ca_rms",  core::scoring::CA_rmsd( native_pose_, *(*it) ));

			//save remarks

			pose.pdb_info( pose::PDBInfoOP( new core::pose::PDBInfo(pose) ));

			pose::PDBInfoOP temp_pdbinfo = pose.pdb_info();

			core::io::RemarkInfo remark;
			remark.value = "RMSD to starting structure: " + utility::to_string(core::scoring::CA_rmsd( native_pose_, *(*it) ));
			temp_pdbinfo->remarks().push_back( remark );

			(*(*it)).pdb_info(temp_pdbinfo);

			std::stringstream SS;
			std::string prefix = option[OptionKeys::out::prefix];
			if ( !prefix.empty() ) {
				SS << prefix << "_" << filecount << ".pdb";
			} else {
				SS << filecount << ".pdb";
			}
			// this is to make sure that the final scoring is done with SCORE12
			scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();

			bool pass_RGF_filter=true; //default to true

			if ( option[OptionKeys::remodel::repeat_structure].user() ) {
				//Experiment with RemodelGlobalFrame
				RemodelGlobalFrame RGF(remodel_data, working_model, scorefxn);
				RGF.align_segment(*(*it));
				RGF.apply(*(*it));

				Real rise = RGF.rise();
				Real radius = RGF.radius();
				Real omega = RGF.omega();

				//filter on the values; this really should be turned into a separate
				//class
				if ( option[OptionKeys::remodel::filter_rise].user() ) {
					utility::vector1<Real> values = option[OptionKeys::remodel::filter_rise];
					if ( values[1] > rise || values[2] < rise ) {
						pass_RGF_filter=false;
						TR << "failed rise filter" << std::endl;
					}
				}
				if ( option[OptionKeys::remodel::filter_radius].user() ) {
					utility::vector1<Real> values = option[OptionKeys::remodel::filter_radius];
					if ( values[1] > radius || values[2] < radius ) {
						pass_RGF_filter=false;
						TR << "failed radius filter" << std::endl;
					}
				}
				if ( option[OptionKeys::remodel::filter_omega].user() ) {
					utility::vector1<Real> values = option[OptionKeys::remodel::filter_omega];
					if ( values[1] > omega || values[2] < omega ) {
						pass_RGF_filter=false;
						TR << "failed omega filter" << std::endl;
					}
				}
			}

			//save structure, unless doing fast_disulfide.
			//with fast disulfide, structures are being passed along in rosetta scripts, so no need to output a second PDB.
			if ( !rosetta_scripts_fast_disulfide_ ) {
				if ( pass_RGF_filter ) {
					(*(*it)).pdb_info()->obsolete(true);
					(*(*it)).dump_scored_pdb(SS.str(), *scorefxn);
				}
			}


			Real score = 0.0;

			//rank poses by score, unless we are doing fast_disulfide, in which case we want to rank by entropy.
			if ( rosetta_scripts_fast_disulfide_ ) {
				simple_filters::DisulfideEntropyFilterOP DisulfideEntropy( new simple_filters::DisulfideEntropyFilter() );
				score = - DisulfideEntropy->compute_residual( *(*it) );
			} else {
				simple_filters::ScoreTypeFilter const pose_total_score( scorefxn, total_score, 100 );
				score = pose_total_score.compute( *(*it) );
			}

			//reload the pose into the new accumulator
			accumulator_.apply(*(*it), score);

			filecount++;
		}

	}
	//set the pose to the best pose in the accumulator and remove that pose.
	if ( accumulator_.size() > 0 ) {
		pose = accumulator_.pop();
	}


	TR << "Remodel poses remaining from original run: " << accumulator_.size() << std::endl;

	// setup calculators
	pose::metrics::CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc_name(),
		PoseMetricCalculatorOP( new toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() ) ) );

	/*
	// do design-refine iteration
	if ( dr_cycles_ > 0 ) {
	if ( !design_refine( pose ) ) { // design-refine failed
	set_last_move_status( FAIL_RETRY );
	return;
	}
	}
	*/

	// if we've gotten to this point, then the structure has been built properly
	set_last_move_status( MS_SUCCESS );

	// setup the PoseMetricCalculators and add them to the evaluators in the JobOutputter
	pose::metrics::CalculatorFactory::Instance().remove_calculator( loops_buns_polar_calc_name() );
	pose::metrics::CalculatorFactory::Instance().remove_calculator( neighborhood_buns_polar_calc_name() );

	pose::metrics::CalculatorFactory::Instance().register_calculator(
		loops_buns_polar_calc_name(),
		PoseMetricCalculatorOP( new toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "default", "default", manager_.union_of_intervals_containing_undefined_positions() ) )
	);

	basic::MetricValue< std::set< Size > > loops_neighborhood;
	pose.metric( neighborhood_calc_name(), "neighbors", loops_neighborhood );
	pose::metrics::CalculatorFactory::Instance().register_calculator(
		neighborhood_buns_polar_calc_name(),
		PoseMetricCalculatorOP( new toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "default", "default", loops_neighborhood.value() ) )
	);

}


/// @brief get_name function for JobDistributor
std::string RemodelMover::get_name() const {
	return "RemodelMover";
}

///
/// @brief
/// Does the same as function below, but takes in a BuildManager object.
/// Also checks the value of the option -remodel::bypass_fragments.
///
bool RemodelMover::centroid_build( Pose & pose, protocols::forge::build::BuildManager & manager ) {

	manager_ = manager;
	if ( option[OptionKeys::remodel::bypass_fragments] || rosetta_scripts_bypass_fragments_ || rosetta_scripts_fast_disulfide_ ) {
		TR << "-=BYPASSING FRAGMENT BUILD (REFINE ONLY) =-" << std::endl;

		if ( option[OptionKeys::remodel::repeat_structure].user() ) { // make sure that the repeats are created even if fragments are bypassed

			bool denovo = true;
			//find out if it's actually denovo
			for ( int i = 0, ie = (int)remodel_data_.blueprint.size(); i < ie; i++ ) {
				if ( remodel_data_.blueprint[i].sstype == "." ) { //if anywhere hits this assignment, not de novo
					denovo = false;
				}
			}
			if ( denovo ) {
				//this part really needs work....  currently doesn't allow growing a loop
				//in regional repeat building.  This section is used in de novo rebuild
				//cases where the monomer pose is extended, so to restore the sidechain,
				//the source of the sidechains has to be extended too in de novo cases.
				//but with refining an existing repeat pose, still need to extend here, unlike the other RLM instance later
				using namespace protocols::loops;
				using protocols::forge::methods::intervals_to_loops;
				std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();

				LoopsOP loops( new Loops() );

				//Temporary fix: special case for denovo type, needed artificially
				//padding residues.  in reality all loops would be padded, but for
				//archive pose extension, it is only special for de novo case
				if ( loop_intervals.size() == 1 && (*(loop_intervals.begin())).left == 1 && (*(loop_intervals.begin())).right == remodel_data_.blueprint.size() ) {
					loops->add_loop( Loop(1, remodel_data_.blueprint.size()+2, 0, 0, true) );
				} else {
					loops = LoopsOP( new Loops( intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) ));
				}

				RemodelLoopMover RLM(loops);
				RLM.set_repeat_tail_length(remodel_data_.sequence.length());
				Pose bufferPose(pose);
				//due to code change, modified_archive_pose would always be 2x length now
				RLM.repeat_generation_with_additional_residue( bufferPose, pose );
				if ( option[OptionKeys::symmetry::symmetry_definition].user() ) {
					//symmetrize if both rep+sym are used
					simple_moves::symmetry::SetupForSymmetryMover pre_mover;
					pre_mover.apply(pose);
					pose.pdb_info()->obsolete(true);
				}
			}
		}

		return true;
	}

	if ( centroid_build( pose ) ) {
		//update external manager
		//manager = manager_;
		return true;
	} else {
		return false;
	}

}

///
/// @brief
/// Runs the centroid level build state. Returns true if loop was closed, false if not.
/// Also checks the value of options -remodel:use_same_length_fragments, -remodel:use_blueprint_sequence
/// and -remodel:repeat_structure
///
bool RemodelMover::centroid_build( Pose & pose ) {

	using namespace basic::options;
	using namespace core::scoring;
	using protocols::moves::MS_SUCCESS;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::restore_residues;
	//using protocols::toolbox::pose_manipulation::construct_poly_uniq_restype_pose;
	using namespace protocols::forge::components;

	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;
	//manager_.modify( modified_archive_pose );

	// ensure modified_archive_pose is completely full-atom, otherwise mismatch
	// will occur when restoring sidechains at the end of the procedure
	bool mod_ap_is_full_atom = true;
	for ( Size i = 1, ie = modified_archive_pose.n_residue(); mod_ap_is_full_atom && i != ie; ++i ) {
		mod_ap_is_full_atom &= ( modified_archive_pose.residue( i ).residue_type_set()->name() == core::chemical::FA_STANDARD );
	}

	if ( !mod_ap_is_full_atom ) {
		core::util::switch_to_residue_type_set( modified_archive_pose, core::chemical::FA_STANDARD );
	}
	/*
	// flip to poly-ala-gly-pro-disulf pose, only in the rebuilt segment
	utility::vector1< Size > protein_residues;
	for (std::set<Size>::iterator it=rebuild.begin(), end=rebuild.end(); it != end; it++){
	//for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
	if ( pose.residue( *it ).is_protein() ) {
	protein_residues.push_back( *it );
	TR<< "turning these to ala: " << *it << std::endl;
	}
	}
	TR << "default building restype: " << "ALA" << std::endl;
	construct_poly_uniq_restype_pose( pose, protein_residues, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), true, true, true );
	*/
	// run VLB to build the new section, if no segments have been added/deleted
	// we use the same VLB so that fragment caching works properly
	if ( !vlb_.get() ) {
		vlb_ = VarLengthBuildOP( new VarLengthBuild( manager_ , remodel_data_ ) );
	}
	if ( !working_model_.abego.empty() ) {
		//the following block simply packages the string to feed to vlb
		utility::vector1<std::string> abego_vec;
		for ( Size i = 0; i < working_model_.abego.length(); i++ ) {
			std::string buffer;
			buffer.push_back(working_model_.abego[i]);
			abego_vec.push_back(buffer);
		}
		vlb_->set_abego(abego_vec);
	}

	vlb_->scorefunction( centroid_sfx_ );
	vlb_->use_fullmer(option[OptionKeys::remodel::use_same_length_fragments] );
	vlb_->max_linear_chainbreak( max_linear_chainbreak_ );
	vlb_->loop_mover_str( centroid_loop_mover_str_ );
	vlb_->restart_mode(true);
	vlb_->new_secondary_structure_override(working_model_.ss);
	if ( option[OptionKeys::remodel::use_blueprint_sequence] ) {
		if ( option[OptionKeys::remodel::repeat_structure].user() ) {
			Size copies =option[OptionKeys::remodel::repeat_structure];
			String rep_seq = remodel_data_.sequence;
			while ( copies > 1 ) {
				rep_seq.append(remodel_data_.sequence);
				copies--;
			}
			vlb_->new_sequence_override( rep_seq );
		} else vlb_->new_sequence_override( remodel_data_.sequence );
	}

	TR << "centroid_build(): calling VariableLengthBuild apply." << std::endl;
	vlb_->apply( pose );

	if ( vlb_->get_last_move_status() == MS_SUCCESS ) {

		// record the used manager w/ all mapping info
		//manager_ = vlb_->manager();

		// safety, clear all the energies before restoring full-atom residues and scoring
		pose.energies().clear();
		bool denovo = true;
		//have to loop to identify denovo case
		for ( int i = 0, ie = (int)remodel_data_.blueprint.size(); i < ie; i++ ) {
			if ( remodel_data_.blueprint[i].sstype == "." ) { //if anywhere hits this assignment, not de novo
				denovo = false;
			}
		}
		if ( denovo ) {
			modified_archive_pose = pose;
		}

		if ( option[OptionKeys::remodel::repeat_structure].user() && !denovo ) { // the previous step already swap modified_archive_pose to be full length. skip the following
			//this part really needs work....  currently doesn't allow growing a loop
			//in regional repeat building.  This section is used in de novo rebuild
			//cases where the monomer pose is extended, so to restore the sidechain,
			//the source of the sidechains has to be extended too in de novo cases.
			//but with refining an existing repeat pose, no need to extend
			// if (modified_archive_pose.total_residue() == pose.total_residue()){ //dangerous, assuming no further length change
			//do nothing.
			// } else { // if there's mismatch, assuming restoration source need extension... dangerous.
			// because of code-change for always using 2x modified
			using namespace protocols::loops;
			using protocols::forge::methods::intervals_to_loops;
			std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();

			LoopsOP loops( new Loops() );

			//Temporary fix: special case for denovo type, needed artificially
			//padding residues.  in reality all loops would be padded, but for
			//archive pose extension, it is only special for de novo case
			if ( loop_intervals.size() == 1 && (*(loop_intervals.begin())).left == 1 && (*(loop_intervals.begin())).right == remodel_data_.blueprint.size() ) {
				loops->add_loop( Loop(1, remodel_data_.blueprint.size()+2, 0, 0, true) );
			} else {
				loops = LoopsOP( new Loops( intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) ) );
			}

			RemodelLoopMover RLM(loops);
			RLM.set_repeat_tail_length(remodel_data_.sequence.length());
			Pose bufferPose(modified_archive_pose);
			//due to code change, modified_archive_pose would always be 2x length now
			RLM.repeat_generation_with_additional_residue( bufferPose, modified_archive_pose );
			if ( option[OptionKeys::symmetry::symmetry_definition].user() ) {
				//symmetrize if both rep+sym are used
				simple_moves::symmetry::SetupForSymmetryMover pre_mover;
				pre_mover.apply(modified_archive_pose);
				modified_archive_pose.pdb_info()->obsolete(true);
			}

			// }
		}

		// Swap back original sidechains.  At the moment this is a two step process
		// in case any sidechains from SegmentInsert and the like that aren't in the
		// original archive pose need to be transferred.
		restore_residues( modified_archive_pose, pose );

		// since pose is setup modified in RemodelMover, only one step will do
		//restore_residues( manager_.original2modified(), archive_pose, pose );
		// go ahead and score w/ full-atom here; we do this in case there are no
		// design-refine cycles -- it's useful to have e.g. rama in the output
		(*fullatom_sfx_)( pose );

		TR << "centroid_build(): variable length build succeeded." << std::endl;
		if ( option[OptionKeys::remodel::repeat_structure].user() ) {
			//return the modified pose to original state, otherwise it keeps growing.
			modified_archive_pose = archive_pose;
		}

		return true; // loop closed
	} else {
		TR << "centroid_build(): variable length build failed. resetting pose to archived pose." << std::endl;
		// set passed in pose reference to the pose this function was called with?  I guess that means the rebuild tried
		// to happen, failed, and so we don't change anything.
		pose = archive_pose;

	}
	TR << "centroid_build(): centroid_build unable to close loop. retrying." << std::endl;
	return false; // false if loop not closed
}

///
/// @brief
/// Sets up constraints and a modified scorefunction and run design/relax cycles.
/// Checks the value of -remodel:repeat_structure.
///
bool RemodelMover::design_refine_seq_relax( Pose & pose, RemodelDesignMover & designMover ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace core::pose::symmetry;
	using namespace protocols;

	// collect new regions/positions
	std::set< forge::build::Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	forge::build::BuildManager::Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	loops::Loops loops = forge::methods::intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

	if ( option[OptionKeys::remodel::repeat_structure].user() ||option[OptionKeys::remodel::free_relax] ) {
		//do nothing
	} else {
		protocols::forge::methods::fill_non_loop_cst_set(pose, loops);
	}

	// safety, clear the energies object
	pose.energies().clear();

	// for refinement, always use standard repulsive
	ScoreFunctionOP sfx = core::scoring::get_score_function();

	// turning on weights
	sfx->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfx->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfx->set_weight( core::scoring::angle_constraint, 1.0 );
	sfx->set_weight( core::scoring::dihedral_constraint, 10.0 ); // 1.0 originally
	sfx->set_weight( core::scoring::res_type_constraint, 1.0);
	sfx->set_weight( core::scoring::res_type_linking_constraint, 1.0);
	protocols::relax::FastRelax relaxMover(sfx);


	scoring::constraints::ConstraintSetOP cst_set_post_built;
	//if (option[OptionKeys::remodel::repeat_structure].user() ) {
	// at this stage it should hold generic cstfile and res_type_linking constraints
	cst_set_post_built = scoring::constraints::ConstraintSetOP( new scoring::constraints::ConstraintSet( *pose.constraint_set() ) );
	//}

	Size asym_length;
	if ( is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		asym_length = symm_info->num_independent_residues();
	} else {
		asym_length = pose.total_residue();
	}


	protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
	if ( option[OptionKeys::remodel::repeat_structure].user() ) {

		// Dihedral (NCS) Constraints, need to be updated each mutation cycle for sidechain symmetry
		Size repeat_number =option[OptionKeys::remodel::repeat_structure];
		Size segment_length = asym_length/repeat_number;

		for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		std::stringstream templateRangeSS;
		// take care of the terminal repeat, since the numbers are offset.
		templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
		std::stringstream targetSS;
		targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
		TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
		setup_ncs.add_group(templateRangeSS.str(), targetSS.str());

	}

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		TR << "design_refine_seq_relax(): dr_cycle: " << i << std::endl;

		designMover.set_state("finish");
		designMover.apply(pose);

		//update dihedral constraint for repeat structures
		if ( option[OptionKeys::remodel::repeat_structure].user() ) {
			setup_ncs.apply(pose);

			//total hack (for now), see if the restypeset fails when initialized twice.
			if ( option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide]() ) {
				protocols::relax::cyclize_pose(pose);
			}

			sfx->show(TR, pose);
			TR << std::endl;
		}

		TR << "design_refine_seq_relax(): calling RelaxMover apply()." << std::endl;
		relaxMover.apply(pose);

		// reset constraints without NCS
		pose.constraint_set(cst_set_post_built);
		TR << "\n";
		sfx->show(TR, pose);
		TR << std::endl;
	}


	// turning off weights
	sfx->set_weight(core::scoring::coordinate_constraint, 0.0 );
	sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	sfx->set_weight(core::scoring::angle_constraint, 0.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 0.0 );
	sfx->set_weight(core::scoring::res_type_constraint, 0.0);
	sfx->set_weight(core::scoring::res_type_linking_constraint, 0.0);

	(*sfx)( pose );

	return true;
}

bool RemodelMover::design_refine_cart_relax(
	Pose & pose,
	RemodelDesignMover & designMover
)
{
	using core::kinematics::FoldTree;
	using namespace core::pose::symmetry;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::RestrictToRepacking;
	//using core::scoring::STANDARD_WTS;
	//using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::forge::build::SegmentInsert;
	using namespace protocols::loops;
	using protocols::loops::Loops;
	using protocols::loops::loop_mover::refine::LoopMover_Refine_CCD;
	using protocols::simple_moves::PackRotamersMover;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;
	using namespace core::scoring::constraints;
	using namespace basic::options;


	using core::pose::annotated_to_oneletter_sequence;
	using protocols::forge::methods::intervals_to_loops;
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::remove_cutpoint_variants;

	//typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	Loops loops = intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

	if ( option[OptionKeys::remodel::repeat_structure].user() ||option[OptionKeys::remodel::free_relax] ) {
		//do nothing
	} else {
		protocols::forge::methods::fill_non_loop_cst_set(pose, loops);
	}

	// safety, clear the energies object
	pose.energies().clear();

	Size asym_length;

	//set simple tree
	FoldTree minFT;
	if ( is_symmetric(pose) ) {
		minFT = sealed_symmetric_fold_tree( pose );
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		asym_length = symm_info->num_independent_residues();
	} else {
		minFT.simple_tree(pose.total_residue());
		asym_length = pose.total_residue();
	}

	pose.fold_tree(minFT);

	// for refinement, always use standard repulsive
	ScoreFunctionOP sfx = core::scoring::get_score_function();
	//turning on weights
	sfx->set_weight(core::scoring::coordinate_constraint, 1.0 );
	sfx->set_weight(core::scoring::atom_pair_constraint, 1.0 );
	sfx->set_weight(core::scoring::angle_constraint, 1.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 10.0 ); // 1.0 originally
	sfx->set_weight(core::scoring::res_type_constraint, 1.0);
	sfx->set_weight(core::scoring::res_type_linking_constraint, 1.0);
	//sfx->set_weight(core::scoring::cart_bonded, 0.5);
	sfx->set_weight(core::scoring::cart_bonded_angle,  5.0 );
	sfx->set_weight(core::scoring::cart_bonded_length,  1.0 );
	sfx->set_weight(core::scoring::cart_bonded_torsion,  5.0 );


	core::kinematics::MoveMapOP cmmop( new core::kinematics::MoveMap );
	//pose.dump_pdb("pretest.pdb");

	if ( option[OptionKeys::remodel::free_relax]() ) {
		for ( Size i = 1; i<= asym_length; ++i ) {
			cmmop->set_bb(i, true);
			cmmop->set_chi(i, true);
		}
	} else {
		cmmop->import(remodel_data_.natro_movemap_);
		cmmop->import( manager_.movemap() );
	}

	for ( Size i = 1; i<= asym_length; ++i ) {
		// std::cout << "bb at " << i << " " << cmmop->get_bb(i) << std::endl;
		// std::cout << "chi at " << i << " " << cmmop->get_chi(i) << std::endl;
		cmmop->set_chi(i,true);
		// std::cout << "chi at " << i << " " << cmmop->get_chi(i) << std::endl;
	}

	//for (Size i = 1; i<= asym_length; ++i){
	// std::cout << "bbM at " << i << " " << manager_.movemap().get_bb(i) << std::endl;
	// std::cout << "chiM at " << i << " " << manager_.movemap().get_chi(i) << std::endl;
	//}

	//if(option[OptionKeys::remodel::repeat_structure].user()){

	// at this stage it should hold generic cstfile and res_type_linking
	// constraints
	ConstraintSetOP cst_set_post_built( new ConstraintSet( *pose.constraint_set() ));
	//}

	protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
	if ( option[OptionKeys::remodel::repeat_structure].user() ) {
		//Dihedral (NCS) Constraints, need to be updated each mutation cycle for sidechain symmetry

		Size repeat_number =option[OptionKeys::remodel::repeat_structure];
		Size segment_length = asym_length/repeat_number;


		//handle movemap
		for ( Size rep = 1; rep < repeat_number; rep++ ) { // from 1 since first segment don't need self-linking
			for ( Size residue = 1; residue <= segment_length ; residue++ ) {
				if ( cmmop->get_bb(residue) ) {
					cmmop->set_bb(residue+(segment_length*rep),true);
				}
			}
		}

		for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}


		std::stringstream templateRangeSS;
		//take care of the terminal repeat, since the numbers are offset.
		templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
		std::stringstream targetSS;
		targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
		TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
		setup_ncs.add_group(templateRangeSS.str(), targetSS.str());

	}

	for ( Size i = 1; i<= asym_length; ++i ) {
		std::cout << "bb at " << i << " " << cmmop->get_bb(i) << std::endl;
		std::cout << "chi at " << i << " " << cmmop->get_chi(i) << std::endl;
	}

	simple_moves::MinMoverOP minMover( new simple_moves::MinMover( cmmop , sfx , "lbfgs_armijo", 0.01, true ));
	minMover->cartesian(true);

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {


		designMover.set_state("finish");
		designMover.apply(pose);

		//update dihedral constraint for repeat structures
		if ( option[OptionKeys::remodel::repeat_structure].user() ) {
			setup_ncs.apply(pose);

			//total hack (for now), see if the restypeset fails when initialized
			//twice.
			if ( option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide]() ) {
				protocols::relax::cyclize_pose(pose);
			}

			//  sfx->show(TR, pose);
			//  TR << std::endl;
		}

		minMover->apply(pose);

		//reset constraints without NCS
		pose.constraint_set(cst_set_post_built);
		sfx->show(TR, pose);
		TR << std::endl;
	}


	//turning off weights
	sfx->set_weight(core::scoring::coordinate_constraint, 0.0 );
	sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	sfx->set_weight(core::scoring::angle_constraint, 0.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 0.0 );
	sfx->set_weight(core::scoring::res_type_constraint, 0.0);
	sfx->set_weight(core::scoring::res_type_linking_constraint, 0.0);

	(*sfx)( pose );

	return true;
}

///
/// @brief
/// Run the design-refine stage.
/// Checks the value of -remodel:repeat_structure and -remodel:swap_refine_confirm_protocols
/// NOTE: CURRENTLY ALWAYS RETURNS TRUE regardless of if chain breaks test passes or fails
///
bool RemodelMover::design_refine( Pose & pose, RemodelDesignMover & designMover ) {

	using namespace core;
	using namespace protocols;
	using namespace protocols::toolbox::task_operations;
	using core::pack::task::operation::TaskOperationCOP;
	using core::pack::task::operation::ResLvlTaskOperationCOP;

	//using core::pack::task::operation::RestrictResidueToRepacking;
	//using core::pack::task::operation::RestrictResidueToRepackingOP;
	//using core::pack::task::operation::RestrictToRepacking;
	//using protocols::forge::build::SegmentInsert;
	//using protocols::loops::Loops;
	//using protocols::loops::LoopMover_Refine_CCD;
	//using protocols::simple_moves::PackRotamersMover;
	//using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;

	//using core::pose::annotated_to_oneletter_sequence;
	//using protocols::forge::methods::intervals_to_loops;
	//using protocols::forge::methods::linear_chainbreak;
	//using protocols::loops::remove_cutpoint_variants;

	//typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< forge::build::Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	forge::build::BuildManager::Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	loops::LoopsOP loops( new loops::Loops( forge::methods::intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) ) );

	// refine Mover used doesn't setup a fold tree, so do it here
	//FoldTree loop_ft = protocols::forge::methods::fold_tree_from_loops( pose, loops );
	kinematics::FoldTree loop_ft;
	loops::fold_tree_from_loops( pose, *loops, loop_ft, true /*term cut*/);

	// save original fold tree
	kinematics::FoldTree original_ft = pose.fold_tree();

	// define the score function
	//ScoreFunctionOP sfx = fullatom_sfx_->clone();
	//for refinement always use hard repulsive
	scoring::ScoreFunctionOP sfx = core::scoring::get_score_function();

	//turning on weights, for paranoya
	sfx->set_weight(core::scoring::atom_pair_constraint, 1.0 );
	sfx->set_weight(core::scoring::coordinate_constraint, 1.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 10.0 ); // 1.0 originally

	// setup the refine TaskFactory
	pack::task::TaskFactoryOP refine_tf = generic_taskfactory();
	refine_tf->push_back( TaskOperationCOP( new toolbox::task_operations::RestrictToNeighborhoodOperation( neighborhood_calc_name() ) ) );
	refine_tf->push_back( TaskOperationCOP( new pack::task::operation::RestrictToRepacking() ) );

	// safety, clear the energies object
	pose.energies().clear();

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		// design the new section
		//PackRotamersMover design( sfx );
		//design.task_factory( design_tf );
		//design.apply( pose );
		designMover.set_state("finish");
		designMover.apply( pose );

		// set loop topology
		pose.fold_tree( loop_ft );

		if ( !option[OptionKeys::remodel::swap_refine_confirm_protocols]() ) {
			// refine the new section
			loops::loop_mover::refine::LoopMover_Refine_CCDOP refine( new loops::loop_mover::refine::LoopMover_Refine_CCD(loops, sfx) );
			kinematics::MoveMapOP combined_mm( new kinematics::MoveMap() );

			////// fix dna
			for ( Size i=1; i<=pose.total_residue() ; ++i ) {
				if ( pose.residue(i).is_DNA() ) {
					TR << "NATRO movemap setup: turning off DNA bb and chi move for refinement stage" << std::endl;
					remodel_data_.natro_movemap_.set_bb( i, false );
					remodel_data_.natro_movemap_.set_chi( i, false );
				}
			}
			////// end fix dna

			combined_mm->import( remodel_data_.natro_movemap_ );
			combined_mm->import( manager_.movemap() );

			//remodel_data_.natro_movemap_.show(pose.total_residue());
			//manager_.movemap().show(pose.total_residue());
			//combined_mm->show(pose.total_residue());
			//modify task to accept NATRO definition
			utility::vector1<core::Size> natroPositions;
			for ( Size i = 1; i<= pose.total_residue(); i++ ) {
				if ( remodel_data_.natro_movemap_.get_chi(i) == 0 ) {
					natroPositions.push_back(i);
				}
			}

			pack::task::operation::OperateOnCertainResiduesOP natroRes( new pack::task::operation::OperateOnCertainResidues );
			natroRes->residue_indices( natroPositions );
			natroRes->op( ResLvlTaskOperationCOP( new pack::task::operation::PreventRepackingRLT ) );
			refine_tf->push_back( natroRes );

			refine->false_movemap( combined_mm );
			refine->set_task_factory( refine_tf );
			refine->apply( pose );

		} else {
			loops::loop_mover::refine::LoopMover_Refine_KICOP KIC( new loops::loop_mover::refine::LoopMover_Refine_KIC(loops) );
			KIC->apply(pose);
		}

		// remove cutpoint variants -- shouldn't this happen at the end
		// of the refine Mover?
		loops::remove_cutpoint_variants( pose );


#ifdef FILE_DEBUG
		std::stringstream SS;
		SS << "RefineStage" << i << ".pdb";
		pose.dump_pdb( SS.str() );
#endif
	}

	//turning off weights, for paranoya
	sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	sfx->set_weight(core::scoring::coordinate_constraint, 0.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 0.0 ); // 1.0 originally

	// must score one last time since we've removed variants and set
	// new topology, otherwise component energies not correct for
	// e.g. structure output
	(*sfx)( pose );

	// evaluate all chainbreaks using linear chainbreak
	bool cbreaks_pass = true;
	for ( loops::Loops::const_iterator l = loops->begin(), le = loops->end(); l != le && cbreaks_pass; ++l ) {
		if ( l->cut() > 0 ) {
			Real const c = forge::methods::linear_chainbreak( pose, l->cut() );
			TR << "design_refine: final chainbreak = " << c  << " at " << l->cut() << std::endl;
			cbreaks_pass = c <= max_linear_chainbreak_;
		}
	}
	// set original topology
	pose.fold_tree( original_ft );

	return cbreaks_pass;
	//return true; //FOR NOW!!  change me back!
}

///
/// @brief
/// As best as I can tell, does some loop closure and calculates RMSD to native. Returns true.
/// NOTE: CURRENTLY ALWAYS RETURNS TRUE regardless of rmsd value, because this stage is not being used as a filter
/// Checks the value of -remodel::swap_refine_confirm_protocols
///
bool RemodelMover::confirm_sequence( core::pose::Pose & pose ) {

	using namespace protocols::forge::methods;
	using protocols::forge::build::Interval;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;
	using core::pack::task::operation::TaskOperationCOP;
	using pack::task::operation::RestrictToRepacking;

	pose::Pose archive_pose = pose;  //for rmsd

	std::set< forge::build::Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	//pose.dump_pdb("pre_KICpose.pdb");

	// collect loops
	loops::LoopsOP confirmation_loops( new loops::Loops( intervals_to_confirmation_loops( loop_intervals.begin(), loop_intervals.end(), pose.total_residue() ) ) );

	// refine Mover used doesn't setup a fold tree, so do it here
	kinematics::FoldTree loop_ft;
	loops::fold_tree_from_loops( pose, *confirmation_loops, loop_ft, true );
	TR << "confirmation loops tree" << loop_ft << std::endl;

	// save original fold tree
	kinematics::FoldTree original_ft = pose.fold_tree();

	// switch over to new tree
	pose.fold_tree(loop_ft);

	//LoopMover_Refine_KIC KIC(confirmation_loops);

	TR << "fold tree entering confirmation: " << pose.fold_tree() << std::endl;

	//KIC.apply(pose);

	if ( option[OptionKeys::remodel::swap_refine_confirm_protocols]() ) {
		TR << "REFINE USING CCD" << std::endl;
		// refine the new section
		// setup the refine TaskFactory
		//
		// protocols::forge::remodel::RemodelLoopMover scramble_mover(confirmation_loops);
		// scramble_mover.randomize_stage(pose);

		TaskFactoryOP refine_tf = generic_taskfactory();
		refine_tf->push_back( TaskOperationCOP( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) ) );
		refine_tf->push_back( TaskOperationCOP( new RestrictToRepacking() ) );

		loops::loop_mover::refine::LoopMover_Refine_CCD refine( confirmation_loops, fullatom_sfx_ );
		kinematics::MoveMapOP combined_mm( new kinematics::MoveMap() );

		////// fix dna
		for ( Size i=1; i<=pose.total_residue() ; ++i ) {
			if ( pose.residue( i ).is_DNA() ) {
				TR << "NATRO movemap setup: turning off DNA bb and chi move for refinement stage" << std::endl;
				remodel_data_.natro_movemap_.set_bb( i, false );
				remodel_data_.natro_movemap_.set_chi( i, false );
			}
		}
		////// end fix dna

		combined_mm->import(remodel_data_.natro_movemap_);
		combined_mm->import( manager_.movemap() );

		refine.false_movemap( combined_mm );
		refine.set_task_factory( refine_tf );
		refine.apply( pose );

	} else {
		TR << "REFINE USING KIC" << std::endl;
		loops::loop_mover::refine::LoopMover_Refine_KIC KIC( confirmation_loops );
		KIC.apply(pose);
	}

	// reset to original foldtree
	pose.fold_tree( original_ft );

	//pose.dump_pdb("post_KICpose.pdb");

	// rmsd_calculation:

	Real sum_sd = 0;
	Real sum_sd_native = 0;
	Real sum_sd_archive2native=0;
	Size atom_count = 0;

	for ( loops::Loops::iterator it = confirmation_loops->v_begin(), end = confirmation_loops->v_end(); it!=end; ++it ) {
		for ( Size i = it->start(); i <= it->stop(); ++i ) {
			Real dist_squared = ( pose.residue(i).xyz( "CA" ) - archive_pose.residue(i).xyz( "CA" ) ).length_squared();
			Real dist_squared_native = ( pose.residue(i).xyz( "CA" ) - native_pose_.residue(i).xyz( "CA" ) ).length_squared();
			Real dist_squared_archive2native = ( archive_pose.residue(i).xyz( "CA" ) - native_pose_.residue(i).xyz( "CA" ) ).length_squared();
			sum_sd = sum_sd + dist_squared;
			sum_sd_native = sum_sd_native + dist_squared_native;
			sum_sd_archive2native = sum_sd_archive2native + dist_squared_archive2native;
			atom_count++;
#ifdef FILE_DEBUG
				std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared= " << dist_squared << std::endl;
				std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared_native= " << dist_squared_native << std::endl;
				std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared_archive2native= " << dist_squared_archive2native << std::endl;
#endif
		}
	}

	sum_sd = sum_sd / atom_count;
	sum_sd_native = sum_sd_native / atom_count;
	sum_sd_archive2native = sum_sd_archive2native / atom_count;

	Real rmsd = sqrt(sum_sd);
	Real rmsd_native = sqrt(sum_sd_native);
	Real rmsd_archive2native = sqrt(sum_sd_archive2native);

	pose::PDBInfoOP temp_pdbinfo = pose.pdb_info();

	io::RemarkInfo remark;
	remark.value = "KIC confirmation RMSD: " + utility::to_string( rmsd ) + " to native RMSD: " + utility::to_string( rmsd_native );
	temp_pdbinfo->remarks().push_back( remark );

	remark.value = " ARCHIVE2NATIVE RMSD: " + utility::to_string(rmsd_archive2native);

	temp_pdbinfo->remarks().push_back( remark );

	pose.pdb_info(temp_pdbinfo);

	TR << "RMSD of KIC conformation: " << rmsd << std::endl;
	TR << "RMSD of KIC conformation to native: " << rmsd_native << std::endl;
	TR << "RMSD of ARCHIVE to NATIVE: " << rmsd_archive2native << std::endl;

	// currently the confirmation stage is not setup as filter so always return true
	if ( rmsd <= 1 ) {
		return true;
	} else {
		return true; //for now CHANGE IT BACK!!
	}

}

///
/// @brief
/// Returns a TaskFactory useable as a starting point for either design or refinement.
/// Only adds the NoRepackDisulfides, IncludeCurrent, and Init from command line ops. ReadResfile is not included.
///
RemodelMover::TaskFactoryOP RemodelMover::generic_taskfactory() {
	using protocols::toolbox::task_operations::LimitAromaChi2Operation;

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP tf( new TaskFactory() );

	tf->push_back( TaskOperationCOP( new InitializeFromCommandline() ) ); // also inits -ex options
	tf->push_back( TaskOperationCOP( new IncludeCurrent() ) ); // enforce keeping of input sidechains
	tf->push_back( TaskOperationCOP( new NoRepackDisulfides() ) );
	if ( !option[OptionKeys::remodel::design::allow_rare_aro_chi]() ) {
		tf->push_back( TaskOperationCOP( new LimitAromaChi2Operation() ) );
	}

	// load resfile op only if requested
	/*if ( !resfile_.empty() ) {
	ReadResfileOP rrf = new ReadResfile();
	rrf->filename( resfile_ );
	tf->push_back( rrf );
	}
	*/
	return tf;
}

///
/// @brief
/// process a continuous design string, adding appropriate operations to the TaskFactory
///
void RemodelMover::process_continuous_design_string( Interval const & original_interval, String const & design_str,
	Original2Modified const & original2modified_interval_endpoints, TaskFactoryOP design_tf ) {

	using namespace core;
	using core::pack::task::operation::TaskOperationCOP;

	Size const offset = original2modified_interval_endpoints.find( original_interval.left )->second;
	for ( Size i = 0, ie = design_str.length(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		switch ( design_str.at( i ) ) {
		case 's' : // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();
			break;
		case '.' : // protocol default design
			continue;
		default : // regular case, single aa type
			allowed_aa_types[ chemical::aa_from_oneletter_code( design_str.at( i ) ) ] = true;
			break;
		}

		design_tf->push_back( TaskOperationCOP( new pack::task::operation::RestrictAbsentCanonicalAAS( i + offset, allowed_aa_types ) ) );
	}
}

///
/// @brief
/// process a design string containing an insert, adding appropriate operations to the TaskFactory
///
void RemodelMover::process_insert_design_string( Interval const & original_interval, String const & design_str,
	Original2Modified const & original2modified_interval_endpoints, TaskFactoryOP design_tf ) {

	using namespace core;
	using core::pack::task::operation::TaskOperationCOP;

	char const insert_char = forge::build::SegmentInsert::insertion_char();

	// Figure out the number of residues in each section.
	forge::build::Interval const interval(
		original2modified_interval_endpoints.find( original_interval.left )->second,
		original2modified_interval_endpoints.find( original_interval.right )->second
	);

	Size const insert_char_idx = design_str.find( insert_char );
	Size const left_nres = insert_char_idx;
	Size const right_nres = design_str.size() - left_nres - 1;
	Size const insert_nres = interval.length() - left_nres - right_nres;

	// Make setup easy by building a new design string to expand the
	// insertion character into a series of the insertion character
	// the size of the insert.
	String aa = design_str;
	aa.replace( insert_char_idx, 1, insert_nres, insert_char );

	// setup TaskOperations
	pack::task::operation::RestrictResidueToRepackingOP repack_op( new pack::task::operation::RestrictResidueToRepacking() );

	Size const left_offset = interval.left;
	for ( Size i = 0, ie = aa.size(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		if ( aa.at( i ) == insert_char ) { // repack only
			repack_op->include_residue( i + left_offset );
			continue;

		} else if ( aa.at( i ) == 's' ) { // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();

		} else if ( aa.at( i ) == '.' ) { // protocol default design
			continue;

		} else { // regular case, single aa type
			allowed_aa_types[ chemical::aa_from_oneletter_code( aa.at( i ) ) ] = true;
		}

		design_tf->push_back( TaskOperationCOP( new pack::task::operation::RestrictAbsentCanonicalAAS( i + left_offset, allowed_aa_types ) ) );
	}

	design_tf->push_back( repack_op );
}

///
/// @brief
/// return a boolean vector specifying allowed a.a. when designing on the surface
///
utility::vector1< bool > const & RemodelMover::allowed_surface_aa() {
	using core::chemical::aa_from_oneletter_code;

	static String surface_aa = "ADEGHIKLMNPQRSTV";
	static utility::vector1< bool > v( 20, false );

	for ( Size i = 0, ie = surface_aa.length(); i < ie; ++i ) {
		v[ aa_from_oneletter_code( surface_aa.at( i ) ) ] = true;
	}

	return v;
}

/// @brief parse xml
void
RemodelMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/ )

{
	//note that we are being invoked from rosetta scripts
	rosetta_scripts_ = true;

	if ( tag->hasOption("blueprint") ) {
		blueprint_ = tag->getOption<std::string>( "blueprint" );
	} else {
		TR << "Skipping blueprint" << std::endl;
		blueprint_ = "";
	}

	if ( tag->hasOption("build_disulf") ) {
		rosetta_scripts_build_disulfide_ = tag->getOption< bool >( "build_disulf", false );
		if ( rosetta_scripts_build_disulfide_ ) {
			TR << "Setting build_disulfide to true" << std::endl;
		}
	} else {
		rosetta_scripts_build_disulfide_ = false;
	}

	if ( tag->hasOption("fast_disulf") ) {
		rosetta_scripts_fast_disulfide_ = tag->getOption< bool >( "fast_disulf", false );
		if ( rosetta_scripts_fast_disulfide_ ) {
			TR << "Setting fast_disulfide to true" << std::endl;
		}
	} else {
		rosetta_scripts_fast_disulfide_ = false;
	}

	if ( tag->hasOption("quick_and_dirty") ) {
		rosetta_scripts_quick_and_dirty_ = tag->getOption< bool >( "quick_and_dirty", false );
		if ( rosetta_scripts_quick_and_dirty_ ) {
			TR << "Setting quick_and_dirty to true" << std::endl;
		}
	} else {
		rosetta_scripts_quick_and_dirty_ = false;
	}

	if ( tag->hasOption("bypass_fragments") ) {
		rosetta_scripts_bypass_fragments_ = tag->getOption< bool >( "bypass_fragments", false );
		if ( rosetta_scripts_bypass_fragments_ ) {
			TR << "Setting bypass_fragments to true" << std::endl;
		}
	} else {
		rosetta_scripts_bypass_fragments_ = false;
	}

	if ( tag->hasOption("match_rt_limit") ) {
		rosetta_scripts_match_rt_limit_ = tag->getOption< core::Real >( "match_rt_limit", 1.0 );
	} else {
		rosetta_scripts_match_rt_limit_ = 1.0;
	}
	TR << "Setting match_rt_limit " << rosetta_scripts_match_rt_limit_ << std::endl;

	if ( tag->hasOption("min_disulfides") ) {
		rosetta_scripts_min_disulfides_ = tag->getOption< core::Real >( "min_disulfides", 1 );
	} else {
		rosetta_scripts_min_disulfides_ = 1;
	}
	TR << "Setting min_disulfides " << rosetta_scripts_min_disulfides_ << std::endl;

	if ( tag->hasOption("max_disulfides") ) {
		rosetta_scripts_max_disulfides_ = tag->getOption< core::Real >( "max_disulfides", 1 );
	} else {
		rosetta_scripts_max_disulfides_ = 1;
	}
	TR << "Setting max_disulfides " << rosetta_scripts_max_disulfides_ << std::endl;

	if ( tag->hasOption("min_loop") ) {
		rosetta_scripts_min_loop_ = tag->getOption< core::Real >( "min_loop", 1 );
	} else {
		rosetta_scripts_min_loop_ = 1;
	}
	TR << "Setting min_loop " << rosetta_scripts_min_loop_ << std::endl;

	if ( tag->hasOption("keep_current_disulfides") ) {
		rosetta_scripts_keep_current_ds_ = tag->getOption< bool>( "keep_current_disulfides", false);
	} else {
		rosetta_scripts_keep_current_ds_ = false;
	}
	TR << "Setting keep_current_disulfides " << rosetta_scripts_keep_current_ds_ << std::endl;

	if ( tag->hasOption("include_current_disulfides") ) {
		rosetta_scripts_include_current_ds_ = tag->getOption< bool>( "include_current_disulfides", false );
	} else {
		rosetta_scripts_include_current_ds_ = false;
	}
	TR << "Setting include_current_disulfides " << rosetta_scripts_include_current_ds_ << std::endl;

	if ( tag->hasOption("relax_bb_for_disulf") ) {
		relax_bb_for_disulf_ = tag->getOption< bool>( "relax_bb_for_disulf", false );
	} else {
		relax_bb_for_disulf_ = false;
	}
	TR << "Setting relax_bb_for_disulf " << relax_bb_for_disulf_ << std::endl;

	if ( tag->hasOption("use_match_rt") ) {
		use_match_rt_ = tag->getOption< bool>( "use_match_rt", false );
	} else {
		use_match_rt_ = true;
	}
	TR << "Setting use_match_rt " << use_match_rt_ << std::endl;

	if ( tag->hasOption("use_disulf_fa_score") ) {
		use_disulf_fa_score_ = tag->getOption< bool>( "use_disulf_fa_score", false );
	} else {
		use_disulf_fa_score_ = false;
	}
	TR << "Setting use_disulf_fa_score " << use_disulf_fa_score_ << std::endl;

	if ( tag->hasOption("disulf_fa_max") ) {
		disulf_fa_max_ = tag->getOption< core::Real>( "disulf_fa_max", -0.25 );
	} else {
		disulf_fa_max_ = -0.25;
	}
	TR << "Setting disulf_fa_max " << disulf_fa_max_ << std::endl;
}


} // namespace remodel
} // namespace forge
} // namespace protocols
