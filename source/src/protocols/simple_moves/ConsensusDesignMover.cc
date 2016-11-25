// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/ConsensusDesignMover.cc
/// @brief cc file for ConsensusDesignMover
/// @author Florian Richter (floric@u.washington.edu), april 2011

// Unit headers
#include <protocols/simple_moves/ConsensusDesignMover.hh>
#include <protocols/simple_moves/ConsensusDesignMoverCreator.hh>

// Project Headers
#include <basic/options/option.hh>

#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.hh>

//option key includes
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

//utility includes
#include <utility>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.ConsensusDesignMover" );

// c++ headerss
#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace simple_moves {

// XRW TEMP std::string
// XRW TEMP ConsensusDesignMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ConsensusDesignMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ConsensusDesignMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ConsensusDesignMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ConsensusDesignMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ConsensusDesignMover";
// XRW TEMP }


ConsensusDesignMover::ConsensusDesignMover()
: ptask_(/* NULL */), task_factory_(nullptr),
	sfxn_(/* NULL */), invert_task_(false),
	use_seqprof_constraints_(false), sasa_cutoff_(0.0),
	seqprof_(/* NULL */), ignore_pose_profile_length_mismatch_(false)
{}

ConsensusDesignMover::ConsensusDesignMover(
	core::pack::task::PackerTaskCOP ptask,
	core::scoring::ScoreFunctionCOP sfxn
)
: ptask_(std::move(ptask)), task_factory_(/* NULL */),
	sfxn_(std::move(sfxn)), invert_task_(false),
	use_seqprof_constraints_(false), sasa_cutoff_(0.0),
	seqprof_(/* NULL */), ignore_pose_profile_length_mismatch_(false)
{}

ConsensusDesignMover::~ConsensusDesignMover()= default;

protocols::moves::MoverOP
ConsensusDesignMover::clone() const
{
	return protocols::moves::MoverOP( new ConsensusDesignMover( *this ) );
}

protocols::moves::MoverOP
ConsensusDesignMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ConsensusDesignMover() );
}

/// @details this mover is allowed to touch all residues specified
/// as designable in the passed in task, resp. if the invert_task_
/// variable is set to true, all residues specified as non-packable
/// in the task
/// if no task is passed in, all residues will be considered legit
void
ConsensusDesignMover::apply( core::pose::Pose & pose )
{
	//first two safeguards
	if ( !sfxn_ ) sfxn_ = core::scoring::get_score_function();
	if ( use_seqprof_constraints_ && sfxn_->has_zero_weight( core::scoring::res_type_constraint ) ) {
		core::scoring::ScoreFunctionOP newsfxn = sfxn_->clone();
		newsfxn->set_weight( core::scoring::res_type_constraint, 1.0 );
		sfxn_ = newsfxn;
	}

	core::pack::task::PackerTaskCOP task = create_consensus_design_task( pose );

	core::scoring::constraints::ConstraintCOPs seqprof_constraints;

	if ( use_seqprof_constraints_ ) {
		seqprof_constraints = pose.add_constraints( create_sequence_profile_constraints( pose, *task ) );
	}

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		protocols::simple_moves::symmetry::SymPackRotamersMover packer( sfxn_, task );
		packer.apply( pose );
	} else {
		protocols::simple_moves::PackRotamersMover packer( sfxn_, task );
		packer.apply( pose );
	}

	if ( use_seqprof_constraints_ ) {
		if ( !pose.remove_constraints( seqprof_constraints ) ) utility_exit_with_message("Couldn't remove sequence profile constraints after ConsensusDesignMover packing step.");
	}

	(*sfxn_)(pose);

} //apply

/// @details
/// at every position that this mover is allowed to touch,
/// the task will be modified according to task operation SeqprofConsensusOperation
core::pack::task::PackerTaskCOP
ConsensusDesignMover::create_consensus_design_task(
	core::pose::Pose const & pose
)
{

	if ( !ptask_ ) {
		if ( task_factory_ ) ptask_ = task_factory_->create_task_and_apply_taskoperations( pose );
		else {
			ptask_ = core::pack::task::PackerTaskCOP( core::pack::task::PackerTaskOP( new core::pack::task::PackerTask_( pose ) ) );
			if ( invert_task_ ) utility_exit_with_message("invert_task_ set to true even though no task or task_factory was passed in. something probably unclean somewhere.");
		}
	}

	core::pack::task::PackerTaskOP consensus_task( new core::pack::task::PackerTask_( pose ) );
	consensus_task->initialize_from_command_line();
	toolbox::task_operations::SeqprofConsensusOperation seqprof_to;
	seqprof_to.set_ignore_pose_profile_length_mismatch( ignore_pose_profile_length_mismatch_);
	seqprof_to.apply( pose, *consensus_task );
	if ( use_seqprof_constraints_ ) seqprof_ = seqprof_to.seqprof();

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		consensus_task = core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, consensus_task );
		ptask_ = core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, ptask_ );
	}

	utility::vector1< core::Real > residue_sasa;
	bool use_sasa( sasa_cutoff_ > 0.0 );
	if ( use_sasa ) {
		core::id::AtomID_Map< core::Real > dummy;
		core::scoring::calc_per_atom_sasa( pose, dummy, residue_sasa, basic::options::option[ basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius]);
	}

	std::string touched_residues;
	core::Size num_design_residues(0);

	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		bool this_residue_allowed( invert_task_ ? !ptask_->residue_task(i).being_packed() : ptask_->residue_task(i).being_designed() );
		if ( !pose.residue_type( i ).is_protein() ) this_residue_allowed = false;
		if ( use_sasa && (residue_sasa[i] < sasa_cutoff_ ) ) this_residue_allowed = false;

		if ( !this_residue_allowed )  consensus_task->nonconst_residue_task(i).restrict_to_repacking();
		else {
			if ( consensus_task->residue_task( i ).being_designed() ) {
				touched_residues = touched_residues + utility::to_string( i ) + "+";
				num_design_residues++;
			}
		}

	} // loop over pose residues
	TR << num_design_residues << "residues (out of a total of " << pose.size() << ") for consensus design are " << touched_residues << std::endl;

	return consensus_task;
}

core::scoring::constraints::ConstraintCOPs
ConsensusDesignMover::create_sequence_profile_constraints(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & task
) const
{
	core::scoring::constraints::ConstraintCOPs csts;
	core::sequence::SequenceProfileOP temp_sp( new core::sequence::SequenceProfile(*seqprof_) ); //dumb nonconstness of seqprofile in SequenceProfileConstraint makes this necessary :(
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type(i).is_protein() && task.residue_task(i).being_designed() ) {
			csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( pose, i, temp_sp ) ) );
		}
	}
	return csts;
}

// XRW TEMP std::string
// XRW TEMP ConsensusDesignMover::get_name() const {
// XRW TEMP  return "ConsensusDesignMover";
// XRW TEMP }

void
ConsensusDesignMover::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data_map );
	if ( tag->hasOption("invert_task") ) invert_task_ = tag->getOption< bool >("invert_task", true);
	if ( tag->hasOption("use_seqprof_constraints") ) use_seqprof_constraints_ = tag->getOption< bool >("use_seqprof_constraints", true);
	if ( tag->hasOption("sasa_cutoff") ) sasa_cutoff_ = tag->getOption< core::Real >("sasa_cutoff", 1.0);
	//if( tag->hasOption("scorefxn") ) sfxn_ = new core::scoring::ScoreFunction( *data_map.get< core::scoring::ScoreFunction * >("scorefxns", tag->getOption< std::string >("scorefxn")) );
	if ( tag->hasOption("scorefxn") ) sfxn_ = protocols::rosetta_scripts::parse_score_function( tag, data_map );

	if ( tag->hasOption("ignore_pose_profile_length_mismatch") ) ignore_pose_profile_length_mismatch_ = tag->getOption< bool >("ignore_pose_profile_length_mismatch");
}

std::string ConsensusDesignMover::get_name() const {
	return mover_name();
}

std::string ConsensusDesignMover::mover_name() {
	return "ConsensusDesignMover";
}

void ConsensusDesignMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_task_operations(attlist);
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn" );
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "invert_task", xsct_rosetta_bool, "Operate on the residues specified as non-packable in the PackerTask", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_seqprof_constraints", xsct_rosetta_bool, "use sequence profile constraints", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "sasa_cutoff", xsct_real, "skip designing residues with SASA lower than this", "1.0" )
		+ XMLSchemaAttribute( "ignore_pose_profile_length_mismatch", xsct_rosetta_bool, "If true, and the pose/profile mismatch, excess pose residues are marked repackable by SeqprofConsensusOperation.  If false, and they mismatch, it crashes.  If they match, no problems!" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This mover will modify a given task according to a sequence profile and then call the PackRotamersMover. At every position that is designable in the task, AAs that have a probability greater than min_aa_probability_ and higher than the native in the sequence profile will be allowed", attlist );
}

std::string ConsensusDesignMoverCreator::keyname() const {
	return ConsensusDesignMover::mover_name();
}

protocols::moves::MoverOP
ConsensusDesignMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ConsensusDesignMover );
}

void ConsensusDesignMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConsensusDesignMover::provide_xml_schema( xsd );
}


}  // namespace simple_moves
}  // namespace protocols
