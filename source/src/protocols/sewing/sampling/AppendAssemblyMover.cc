// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AppendAssemblyMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/AppendAssemblyMover.hh>
#include <protocols/sewing/sampling/AppendAssemblyMoverCreator.hh>
#include <protocols/sewing/conformation/DisembodiedAssembly.hh>

//Package headers
#include <protocols/sewing/util/io.hh>
#include <protocols/sewing/util/util.hh>
#include <protocols/sewing/conformation/Model.hh>
#include <protocols/sewing/conformation/Assembly.hh>
#include <protocols/sewing/sampling/requirements/RequirementSet.hh>
#include <protocols/sewing/sampling/requirements/ResidueRetentionRequirement.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationCreators.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/util.hh>

#include <core/scoring/AtomVDW.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoringManager.hh>

//#include <protocols/relax/AtomCoordinateCstMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/relax/cst_util.hh>
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/toolbox/task_operations/RestrictResiduesToRepackingOperation.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/keys/relax.OptionKeys.gen.hh>

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.AppendAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP
AppendAssemblyMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AppendAssemblyMover );
}

std::string
AppendAssemblyMoverCreator::keyname() const
{
	return AppendAssemblyMoverCreator::mover_name();
}

std::string
AppendAssemblyMoverCreator::mover_name()
{
	return "AppendAssemblyMover";
}

protocols::moves::MoverOP
AppendAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new AppendAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
AppendAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AppendAssemblyMover );
}

std::string
AppendAssemblyMover::get_name() const {
	return "AppendAssemblyMover";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  AppendAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


AppendAssemblyMover::AppendAssemblyMover():
	parent(),
	initialized_(false)
{}


void
AppendAssemblyMover::add_starting_model(
	AssemblyOP assembly
) const {
	std::set<core::Size> model_nodes = graph_->get_node_indices_from_model_id(-1);
	if ( model_nodes.size() == 0 ) {
		utility_exit_with_message("No pose segments found in the graph!");
	}
	ModelNode const * const node = graph_->get_model_node(*model_nodes.begin());
	assembly->add_model(graph_, node->model());

	//If we have a partner then add it to the Assembly, but don't allow
	//other things to be built off of it
	if ( partner_pose_ ) {
		assembly->set_partner(partner_pose_);
	}
}


void
AppendAssemblyMover::apply(
	core::pose::Pose & pose
) {
	using namespace core;
	using namespace basic::options;

	if ( !initialized_ ) {
		utility::vector1<int> segment_starts;
		utility::vector1<int> segment_ends;
		if ( option[OptionKeys::sewing::pose_segment_starts].user() ) {
			segment_starts = option[OptionKeys::sewing::pose_segment_starts].value();
		} else {
			TR.Warning << "No specified segment starts for input pose, creating a single segment starting at residue 1" << std::endl;
			segment_starts.push_back(1);
		}

		if ( option[OptionKeys::sewing::pose_segment_ends].user() ) {
			segment_ends = option[OptionKeys::sewing::pose_segment_ends].value();
		} else {
			TR.Warning << "No specified segment ends for input pose, creating a single segment ending at the final residue" << std::endl;
			segment_ends.push_back((int)pose.total_residue());
		}

		if ( segment_starts.size() != segment_ends.size() ) {
			utility_exit_with_message("You must supply the same number of segment starts and ends");
		}

		utility::vector1< std::pair<core::Size, core::Size> > segments;
		for ( core::Size i=1; i<=segment_starts.size(); ++i ) {
			if ( segment_ends[i] > (int)pose.total_residue() ) {
				utility_exit_with_message("Specified segment end lies outside PDB.");
			}
			segments.push_back(std::make_pair(segment_starts[i], segment_ends[i]));
		}


		////**** Generate the Model from the PDB, and setup necessary requirements ****////
		int append_pdb_model_id = -1;
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);
		Model pdb_model = create_model_from_pose(pose, segments, append_pdb_model_id);
		TR << "PDB Model: " << std::endl;
		for ( core::Size i=1; i<=pdb_model.segments_.size(); ++i ) {
			TR << "\tSegment " << i << " (" << pdb_model.segments_[i].segment_id_ << ")" << ": "
				<< pdb_model.segments_[i].residues_.front().resnum_ << " - " << pdb_model.segments_[i].residues_.back().resnum_
				<< " (" << pdb_model.segments_[i].residues_.size() << " residues)" << std::endl;
		}
		models_.insert(std::make_pair(pdb_model.model_id_, pdb_model));

		if ( option[OptionKeys::sewing::keep_model_residues].user() ) {
			utility::vector1<int> pdb_residues_to_repack = option[OptionKeys::sewing::keep_model_residues].value();
			sampling::requirements::ResidueRetentionRequirementOP res_retention (
				new sampling::requirements::ResidueRetentionRequirement(append_pdb_model_id) );
			for ( core::Size i = 1; i <= pdb_residues_to_repack.size(); ++i ) {
				if ( pdb_residues_to_repack[i] > (int)pose.total_residue() ) {
					utility_exit_with_message("Residue not present in starting node!");
				}
				res_retention->add_resnum(pdb_residues_to_repack[i]);
			}
			requirement_set_->add_requirement(res_retention);
		}


		////**** Initialize the partner pose ****////
		if ( option[OptionKeys::sewing::partner_pdb].user() ) {
			partner_pose_ = core::import_pose::pose_from_file(option[OptionKeys::sewing::partner_pdb].value(), core::import_pose::PDB_file);
		}

		hash_pdb_model(pdb_model);

	}

	parent::apply(pose);
}

void
AppendAssemblyMover::hash_pdb_model(
	Model const & pdb_model
) {

	using namespace core;
	using namespace basic::options;

	////Figure out which segments from the starting PDB we need to match
	std::set<core::Size> segments_to_match;
	if ( option[OptionKeys::sewing::match_segments].user() ) {
		utility::vector1<int> temp = option[OptionKeys::sewing::match_segments];
		for ( core::Size i=1; i<=temp.size(); ++i ) {
			segments_to_match.insert(core::Size(temp[i]));
		}
	} else {
		utility::vector1<SewSegment> all_segments = pdb_model.segments_;
		for ( core::Size i=1; i<=all_segments.size(); ++i ) {
			segments_to_match.insert(all_segments[i].segment_id_);
		}
	}


	//Now, score the pdb model against all others
	core::Size num_segments_to_match = option[OptionKeys::sewing::num_segments_to_match].value();
	core::Size min_hash_score = option[OptionKeys::sewing::min_hash_score].value();
	core::Size max_clash_score = option[OptionKeys::sewing::max_clash_score].value();

	//init the SewGraph
	std::map<int,Model>::const_iterator it = models_.begin();

	//Breakup the models into groups of 1000 for scoring. This make the memory demand *much* less
	//with a minor hit to scoring speed.
	ScoreResults all_scores;
	core::Size group = 1;
	core::Size starttime = time(NULL);
	while ( true ) {
		++group;
		core::Size counter = 0;
		Hasher hasher;
		for ( ; it != models_.end(); ++it ) {
			hasher.insert(it->second);
			++counter;
			if ( counter > 100 ) {
				break;
			}
		}
		ScoreResults group_scores;

		group_scores = hasher.score(pdb_model, num_segments_to_match, min_hash_score, max_clash_score, segments_to_match, false, 3);
		// SewingAppend just assumes box_length=3
		//group_scores = hasher.score(pdb_model, num_segments_to_match, min_hash_score, max_clash_score, segments_to_match, false);

		all_scores.insert(group_scores.begin(), group_scores.end());
		if ( it == models_.end() ) {
			break;
		}
	}
	core::Size endtime = time(NULL);
	TR << "Scoring complete, found " << all_scores.size() << " unique alignments to input PDB. (" << endtime - starttime << " seconds)" << std::endl;

	//Combine the PDB scores with the score file and run the AssemblyMover
	models_.insert(std::make_pair(pdb_model.model_id_, pdb_model));
	if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "discontinuous" ) {
		graph_.reset( new SewGraph(models_, 2) );
	} else if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "continuous" ) {
		graph_.reset( new SewGraph(models_, 1) );
	}
	graph_->set_special_edges(all_scores);

	initialized_ = true;
}

core::pose::Pose
AppendAssemblyMover::refine_assembly(
	AssemblyOP & assembly
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	/************** Full-atom Stage ****************/
	//Before relaxing and packing
	core::pose::Pose pose = get_fullatom_pose(assembly);

	//Prepeare the TaskFactory for packing from both the command line and from
	//Assemblie's native retention
	core::pack::task::TaskFactoryOP task_factory ( new core::pack::task::TaskFactory() );

	core::pack::task::operation::InitializeFromCommandlineOP command_line(
		new core::pack::task::operation::InitializeFromCommandline() );
	task_factory->push_back(command_line);

	//Restrict the partner pose from designing
	protocols::toolbox::task_operations::RestrictChainToRepackingOperationOP restrict_chain (
		new protocols::toolbox::task_operations::RestrictChainToRepackingOperation(2) );
	task_factory->push_back(restrict_chain);

	//Keep any residues specified by the user to be kept
	if ( option[basic::options::OptionKeys::sewing::keep_model_residues].user() ) {
		utility::vector1<int> pdb_residues_to_repack = option[OptionKeys::sewing::keep_model_residues].value();
		protocols::toolbox::task_operations::RestrictResiduesToRepackingOperationOP restrict_residues(
			new protocols::toolbox::task_operations::RestrictResiduesToRepackingOperation() );

		utility::vector1<core::Size> assembly_residues_to_repack;
		for ( core::Size i=1; i<=pdb_residues_to_repack.size(); ++i ) {
			assembly_residues_to_repack.push_back( assembly->pose_num( -1, pdb_residues_to_repack[i] ) );
		}
		restrict_residues->set_residues(assembly_residues_to_repack);
		task_factory->push_back(restrict_residues);
	}

	assembly->prepare_for_packing(pose, task_factory, base_native_bonus_, neighbor_cutoff_);

	//setup a move map that allows only chain 1 to move
	core::kinematics::MoveMapOP move_map(new core::kinematics::MoveMap());
	move_map->set_jump(false);
	move_map->set_bb(false);
	move_map->set_chi(false);
	move_map->set_bb_true_range(1, pose.conformation().chain_endings().front());
	move_map->set_chi_true_range(1, pose.conformation().chain_endings().front());
	TR << "Move map in refinment stage: " << std::endl;
	move_map->show(TR);

	//FastDesign
	core::Size repeats = option[OptionKeys::relax::default_repeats].value();
	protocols::relax::FastRelaxOP fast_design( new protocols::relax::FastRelax(repeats) );
	fast_design->set_movemap(move_map);
	fast_design->set_scorefxn(fa_scorefxn_);
	fast_design->cartesian(true);
	fast_design->min_type("lbfgs_armijo_nonmonotone");
	fast_design->set_task_factory(task_factory);
	fast_design->constrain_relax_to_start_coords(true);
	fast_design->apply(pose);

	protocols::relax::delete_virtual_residues( pose );

	return pose;
}

void
AppendAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	using namespace basic::options;

	//This is silly, and causes reading the model file (which is actually quite large, twice). An alternative is to
	//modify Random assembly mover to store the map of model ids to models
	std::string model_file;
	if ( tag->hasOption("model_file") ) {
		model_file = tag->getOption<std::string>("model_file");
		models_ = read_model_file(model_file);
	} else if ( option[basic::options::OptionKeys::sewing::model_file_name].user() ) {
		model_file = option[basic::options::OptionKeys::sewing::model_file_name].value();
		models_ = read_model_file(model_file);
	} else {
		utility_exit_with_message("You must give a model file through options or tags");
	}

	parent::parse_my_tag(tag, data, filters, movers, pose);
}


} //sewing
} //protocols
