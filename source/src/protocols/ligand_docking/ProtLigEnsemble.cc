// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ProtLigEnsemble.cc
/// @brief  implementation of docking with protein and ligand ensembles
/// @author Darwin Fu (darwinyfu@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/ProtLigEnsemble.hh>
#include <protocols/ligand_docking/ProtLigEnsembleCreator.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <protocols/ligand_docking/MinimizeLigand.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/TetherLigand.hh>
#include <protocols/ligand_docking/HighResEnsemble.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/chemical/ResidueTypeFinder.hh>
//Options headers
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/ligand_docking/ligand_scores.hh>
// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ResidueTypeConstraint.hh>

// Project Headers
#include <core/scoring/Energies.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <algorithm>
#include <sstream>
#include <utility/string_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <core/kinematics/Edge.hh>
#include <core/pack/task/PackerTask.hh>

// Scripter Headers
#include <utility/tag/Tag.hh>
//#include <protocols/moves/DataMap.hh>

// Boost Headers

//STL headers
#include <string>

#include <set>

//Auto Headers
#include <protocols/ligand_docking/LigandArea.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace ligand_docking {

static basic::Tracer TR("protocols.ligand_docking.ligand_options.Protocol");

///@brief
ProtLigEnsemble::ProtLigEnsemble():
	Mover("ProtLigEnsemble"),
	pose_infos_(),
	ignore_correlation_until_(0),
	ignore_correlation_after_(0),
	distance_(0),
	num_cycles_(0),
	repack_every_Nth_(0),
	score_fxn_(nullptr),
	final_score_fxn_(nullptr),
	correlation_weight_(0),
	exp_ranks_(),
	rosetta_lowest_scores_(),
	rosetta_lowest_poses_()
{
}

ProtLigEnsemble::ProtLigEnsemble(ProtLigEnsemble const & /*that*/) = default;

ProtLigEnsemble::~ProtLigEnsemble() = default;

protocols::moves::MoverOP ProtLigEnsemble::clone() const {
	return protocols::moves::MoverOP( new ProtLigEnsemble( *this ) );
}

bool ProtLigEnsemble::reinitialize_for_each_job() const { return true; }

protocols::moves::MoverOP ProtLigEnsemble::fresh_instance() const {
	return protocols::moves::MoverOP(new ProtLigEnsemble );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
ProtLigEnsemble::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "ProtLigEnsemble" ) {
		utility_exit_with_message("This should be impossible");
	}
	num_cycles_= tag->getOption<core::Size>("cycles");
	repack_every_Nth_= tag->getOption<core::Size>("repack_every_Nth");

	/// Score Function ///
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_name);

	if ( tag->hasOption("qsar_file") ) {
		read_qsar_file(tag->getOption<std::string>("qsar_file"));
	}

	/// Use Final Minimizer as well, which needs final_score and final_move options
	if ( tag->hasOption("final_score") ) {
		std::string scorefxn_name= tag->getOption<std::string>("final_score");
		final_score_fxn_= datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_name);
	}

	/// Correlation Weight ///
	if ( ! basic::options::option[ basic::options::OptionKeys::docking::ligand::ligand_ensemble ].user() ) utility_exit_with_message("'ProtLigEnsemble' requires the docking:ligand:ligand_ensemble option to run");
	correlation_weight_ = basic::options::option[ basic::options::OptionKeys::docking::ligand::ligand_ensemble ]();

	if ( tag->hasOption("ignore_correlation") ) {
		ignore_correlation_until_ = tag->getOption<core::Size>("ignore_correlation", 0);
	}

	if ( tag->hasOption("ignore_correlation_after") ) {
		ignore_correlation_after_ = tag->getOption<core::Size>("ignore_correlation_after", 0);
	} else {
		ignore_correlation_after_ = pose_infos_.size();

	}

	if ( tag->hasOption("distance") ) {
		distance_ = tag->getOption<core::Real>("distance", 7.0);
	}

}

void ProtLigEnsemble::prepare_single_ligand_pose(core::pose::Pose & pose, ProtLigPair_info & info, core::Size count)
{

	if ( !info.wild_type ) {
		//Replace protein residue by mutation if not wild-type
		core::chemical::AA my_aa = core::chemical::aa_from_oneletter_code( info.mut_target );
		core::chemical::ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose(core::chemical::FULL_ATOM_t) );
		core::chemical::ResidueType const & rsd_type( *( residue_set->get_representative_type_aa( my_aa ) ) );
		core::pose::replace_pose_residue_copying_existing_coordinates(pose, info.mut_resid, rsd_type);
	}
	//Save one ligand
	core::conformation::Residue single_ligand = pose.residue(pose.conformation().chain_begin(core::pose::get_chain_id_from_chain(info.lig_chain, pose)));

	//Delete all ligands
	core::pose::remove_nonprotein_residues(pose);
	core::pose::remove_ligand_canonical_residues(pose);

	//Add ligand of interest back
	pose.append_residue_by_jump(single_ligand, pose.size(), "", "", true);

	core::pose::PDBInfo & pdbinfo = *pose.pdb_info();
	pdbinfo.chain( pose.size(), info.lig_chain);

	//core::pose::fix_pdbinfo_damaged_by_insertion(pose);

	//Add initial record of lowest score and pose
	exp_ranks_.push_back(std::make_pair(count,info.bind_data));
	rosetta_lowest_scores_.push_back(std::make_pair(count,(*score_fxn_)(pose)));
	rosetta_lowest_poses_.push_back(pose);

}

void ProtLigEnsemble::read_qsar_file(std::string filename)
{
	//RESID MUTATION LIGAND-CHAIN BIND
	//117 V B 5.35

	utility::vector1<std::string> lines = utility::io::get_lines_from_file_data(filename);

	for ( core::Size i = 1; i <= lines.size(); ++i ) {
		pose_infos_.push_back(process_line(lines[i]));
	}
	//convert bind data to ranks, order doesn't matter if pose files haven't been created yet
	std::sort(pose_infos_.begin(), pose_infos_.end(), sort_by_binding);  //Sorted by binding data
	info_to_rank(pose_infos_);

	//Set ignore_correlation_after equal to number of binding data available
	for ( core::Size i = 1; i <= pose_infos_.size(); ++i ) {
		if ( pose_infos_[i].has_bind == false ) {
			ignore_correlation_after_ = i-1;
			break;
		}
	}
}

ProtLigPair_info process_line(std::string & line)
{
	ProtLigPair_info pair_info;
	utility::vector1<std::string> fields= utility::string_split(line);

	if ( fields.size() < 2 || fields.size() > 4 ) {
		utility_exit_with_message("ProtLigEnsemble: Incorrect format. Either WT LIG-CHAIN BIND-DATA for wild type or RESID MUTATION LIG-CHAIN BIND-DATA for mutants");
	}

	if ( fields[1] == "WT" ) {
		//WT X
		//WT X 1.0
		pair_info.lig_chain = fields[2][0];
		pair_info.wild_type = true;

		//Check to see if no binding data
		if ( fields.size() == 2 ) {
			pair_info.has_bind = false;
			pair_info.bind_data = 0.0; //Just set it to something
		} else {
			TR << "Read bind data" << fields[3]<< std::endl;
			pair_info.has_bind = true;
			pair_info.bind_data = std::stod(fields[3]);
		}
	} else {
		//100 A X
		//100 A X 1.0

		if ( fields[2].size() != 1 || fields[3].size() != 1 ) {
			utility_exit_with_message("ProtLigEnsemble: Protein residue and ligand chain should be a single letter!");
		}

		pair_info.wild_type=false;
		pair_info.mut_resid = std::stoi(fields[1]);
		pair_info.mut_target = fields[2][0];
		pair_info.lig_chain = fields[3][0];

		//Check to see if no binding data
		if ( fields.size() == 3 ) {
			pair_info.has_bind = false;
			pair_info.bind_data = 0.0; //Just set it to something
		} else {
			pair_info.has_bind = true;
			pair_info.bind_data = std::stod(fields[4]);
		}

	}

	return pair_info;
}

void
ProtLigEnsemble::apply(core::pose::Pose & pose) {
	debug_assert(num_cycles_ > 0);

	rosetta_lowest_poses_.clear();
	rosetta_lowest_scores_.clear();
	exp_ranks_.clear();

	//original_pose has all the ligands. Work with pose
	core::pose::Pose original_pose = pose;

	// Create poses by deleting residues as appropriate based on pose_infos_ being passed. Store in Rosetta poses

	//Dock poses in order of binding affinity, starting with the strongest (set pose_infos_ to be in order)
	//Optimize all cycles for strongest first, then move on to next
	for ( core::Size pose_count = 1; pose_count <= pose_infos_.size(); ++pose_count ) {
		//create single pose for docking by modifying original pose
		//rosetta_old_poses_.push_back(pose);

		TR << "Start setup for pose " << pose_count << std::endl;

		pose = original_pose;
		prepare_single_ligand_pose(pose, pose_infos_[pose_count], pose_count);

		//Generate residue selector for docking
		core::select::residue_selector::NeighborhoodResidueSelector interface_selector = core::select::residue_selector::NeighborhoodResidueSelector();
		score_fxn_->score(pose);

		//Set all to false except for ligand
		utility::vector1<bool> focus(pose.total_residue(), false);
		//pose get residue id of ligand from mut info
		core::Size const begin(pose.conformation().chain_begin(core::pose::get_chain_id_from_chain(pose_infos_[pose_count].lig_chain, pose)));
		focus[begin] = true;

		interface_selector.set_focus(focus);
		interface_selector.set_distance (distance_);
		interface_selector.set_include_focus_in_subset(true);

		//Default angstroms/degrees based on old ligand area definitions
		protocols::moves::MoverOP rigid_body_mover( new protocols::rigid::RigidBodyPerturbMover( core::pose::get_jump_id_from_chain(pose_infos_[pose_count].lig_chain, pose), 2.8648, 0.1));

		protocols::moves::MonteCarloOP monte_carlo ( new protocols::moves::MonteCarlo(*score_fxn_, 2.0) );

		//Pose and score storing terms
		core::pose::Pose last_accepted_pose = pose;
		monte_carlo->set_last_accepted(rosetta_lowest_scores_[pose_count].second);
		monte_carlo->set_lowest(rosetta_lowest_scores_[pose_count].second);

		core::Real correlation_before = 0;
		core::Real correlation_after = 0;
		core::Real current_score = 0;
		bool mc_result;
		core::Real score_delta;
		core::Real low_score_so_far; //Temporary storage variable while considering new score
		core::Real score_being_considered; //Current score adjusted with rank correlation

		TR << "Start docking for pose " << pose_count << std::endl;

		for ( core::Size cycle = 1; cycle <= num_cycles_; ++cycle ) {

			utility::vector1< bool > interface_residues = interface_selector.apply(pose);

			TR << "Created interface residues" << std::endl;

			core::pack::task::PackerTaskOP packer_task = make_packer_task(pose, interface_residues);

			TR << "Created packer task" << std::endl;

			protocols::moves::MoverOP pack_mover;

			if ( cycle % repack_every_Nth_ == 1 ) {
				TR << "making PackRotamersMover" << std::endl;
				//Use minPackMOver instead of PackRotamers to do pack and minimization in one go.
				pack_mover = protocols::moves::MoverOP(new protocols::minimization_packing::PackRotamersMover(score_fxn_, packer_task));

				// pack_mover= moves::MoverOP( new protocols::minimization_packing::PackRotamersMover(score_fxn_, packer_task) );
			} else {
				TR << "making RotamerTrialsMover" << std::endl;
				pack_mover = protocols::moves::MoverOP( new protocols::minimization_packing::RotamerTrialsMover(score_fxn_, *packer_task) );
			}

			TR << "Making rigid body moves" << std::endl;
			//Apply movers and update scores
			rigid_body_mover->apply(pose);

			TR << "Making packing moves" << std::endl;
			pack_mover->apply(pose);
			current_score = (*score_fxn_)(pose);

			TR << "Move made and scored" << std::endl;

			//Store the existing low score aside in case we need it later
			low_score_so_far = rosetta_lowest_scores_[pose_count].second;

			//Calculate the previous correlation for comparison
			correlation_before = qsar_correlation();

			//Put current score into low so we can do the QSAR correlation calculation
			rosetta_lowest_scores_[pose_count].second = current_score;

			//Calculate QSAR Correlation
			correlation_after = qsar_correlation();

			//Restore original low score
			rosetta_lowest_scores_[pose_count].second = low_score_so_far;

			TR << "Correlation To Use for Adjustment in Cycle " << cycle << " is " << correlation_after << std::endl;

			//Adjust Score
			score_being_considered = current_score - (correlation_after - correlation_before)*correlation_weight_; //This formula assumes that positive correlation is preferred, change sign if this is not the case (with -log K for example)
			TR << "Adjusted Current Score is " << score_being_considered << "\n";

			//Find difference
			score_delta = score_being_considered - monte_carlo->last_accepted_score();
			TR << "Score Delta is " << score_delta << "\n";

			//MC Accept or Reject
			mc_result = monte_carlo->boltzmann(score_being_considered);

			//If accepted, copy current pose into old pose and update score. Check if it beats low and update if appropriate
			if ( mc_result ) {
				TR << "Monte Carlo Accepted" << "\n";
				last_accepted_pose = pose;
				monte_carlo->set_last_accepted(current_score);

				if ( current_score < rosetta_lowest_scores_[pose_count].second ) {
					TR << "Beat lowest score" << "\n";
					rosetta_lowest_scores_[pose_count].second = current_score;
					rosetta_lowest_poses_[pose_count] = pose;
					monte_carlo->set_lowest(rosetta_lowest_scores_[pose_count].second);
				}
			} else {
				//if reject, restore the old pose/score into current
				TR << "rejected" << "\n";
				pose = last_accepted_pose;
			}

			correlation_after = qsar_correlation();
			TR << "Correlation For Cycle " << cycle << " is " << correlation_after << std::endl;

		}

		//End cycle of one pose

		//Final minimized if being used (Just operates on the lowest pose)
		if ( final_score_fxn_ ) {

			TR << "Performing final minimization for pose " << pose_count << std::endl;

			core::kinematics::MoveMapOP movemap(new core::kinematics::MoveMap());
			utility::vector1< bool > interface_residues = interface_selector.apply(rosetta_lowest_poses_[pose_count]);

			//Use residue subset interface to build movemap
			for ( core::Size i=1; i <= interface_residues.size(); ++i ) {
				if ( interface_residues[i] == true ) {
					movemap->set_chi(i, true);
					movemap->set_bb(i, true);
				}
			}

			//From FinalMinimizer
			protocols::minimization_packing::MinMoverOP minmover(new protocols::minimization_packing::MinMover(movemap, final_score_fxn_, "lbfgs_armijo_nonmonotone_atol", 0.02, true));
			minmover->min_options()->nblist_auto_update(true);
			minmover->apply(rosetta_lowest_poses_[pose_count]);
			rosetta_lowest_scores_[pose_count].second=(*final_score_fxn_)(rosetta_lowest_poses_[pose_count]);


		}

		//End cycle of all poses
	}

	//Calculate final correlation
	core::Real final_correlation = qsar_correlation();

	//Output individual poses, scores and correlation

	protocols::jd2::add_string_real_pair_to_current_job("spearman", final_correlation);

	for ( core::Size pose_count = 1; pose_count <= pose_infos_.size(); ++pose_count ) {

		std::string tag( pose_infos_[pose_count].print_string() + "_" + utility::to_string(protocols::jd2::current_nstruct_index()) + ".pdb" );

		//output individual poses
		rosetta_lowest_poses_[pose_count].dump_scored_pdb(tag, *final_score_fxn_);

		if ( basic::options::option[ basic::options::OptionKeys::out::pdb_gz ]() ) {
			utility::file::gzip( tag, true );
		}

		//Add interface delta scores
		std::map< std::string, core::Real > ligand_scores( get_interface_deltas( pose_infos_[pose_count].lig_chain, rosetta_lowest_poses_[pose_count], final_score_fxn_, pose_infos_[pose_count].print_mut_string() ) );

		for ( auto const & entry : ligand_scores ) {
			protocols::jd2::add_string_real_pair_to_current_job( entry.first, entry.second );
		}

	}

}

core::Real ProtLigEnsemble::qsar_correlation()
{
	//If number of scores smaller than ignore, than return 0 to prevent too many poses from being rejected with small subset.
	//ignore_correlation_after set to ignore correlations once there is no binding data available
	if ( rosetta_lowest_scores_.size() < ignore_correlation_until_ || rosetta_lowest_scores_.size() > ignore_correlation_after_ ) {

		TR << "Ignoring correlation due to dataset size " << rosetta_lowest_scores_.size() << std::endl;
		return 0;
	}

	//Current run rosetta scores in rosetta_sorted_scores and experimental values in qsar_sorted_scores
	utility::vector1<std::pair<core::Size, core::Real> > rosetta_ranks = rosetta_lowest_scores_;

	//Convert rosetta_ranks from scores to rank
	std::sort(rosetta_ranks.begin(), rosetta_ranks.end(), sort_by_second);  //Sorted by assay now
	vector_to_rank(rosetta_ranks);
	std::sort(rosetta_ranks.begin(), rosetta_ranks.end());  //Restore to chain ID order to match exp_ranks_

	return spearman(exp_ranks_, rosetta_ranks);

}

core::pack::task::PackerTaskOP ProtLigEnsemble::make_packer_task(core::pose::Pose const & pose, utility::vector1<bool> const & interface_residues)
{
	//check packing for 2nd pose and so on
	static bool pose_already_packed=false;

	core::pack::task::PackerTaskOP pack_task( core::pack::task::TaskFactory::create_packer_task(pose) );
	pack_task->initialize_from_command_line();

	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot_( new core::pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot_->initialize_from_command_line();
	pack_task->append_rotamerset_operation( unboundrot_ );

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		/// If several params files have the same name, allow switching among them
		/// This was previously only enabled with mutate_same_name3.  Now default.
		if ( ! pose.residue(i).is_ligand() ) continue;
		enable_ligand_rotamer_packing(pose, i, pack_task);
	}

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! pose.residue(i).is_ligand() ) {
			pack_task->nonconst_residue_task( i ).restrict_to_repacking();
		}
	}


	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( interface_residues[i] == false  ) {
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	// We always want the option (after the initial unbiased pack)
	// of sticking with our current nicely minimized conformation.
	if ( pose_already_packed ) {
		pack_task->or_include_current(true);
	} else {
		pose_already_packed=true;
	}

	return pack_task;

}

void
ProtLigEnsemble::enable_ligand_rotamer_packing(
	core::pose::Pose const & pose,
	core::Size const ligand_residue_id,
	core::pack::task::PackerTaskOP & pack_task
) const{
	core::conformation::Residue const & this_residue= pose.residue(ligand_residue_id);

	core::chemical::ResidueTypeSetCOP rsd_type_set = pose.residue_type_set_for_pose( this_residue.type().mode() );
	core::chemical::ResidueTypeCOPs allowed_types = core::chemical::ResidueTypeFinder( *rsd_type_set ).name3( this_residue.name3() ).get_all_possible_residue_types(); // a vector1

	debug_assert(allowed_types.size() > 0);
	/// TODO consider removing this so resfiles can specify ligand mutations to allow
	if ( allowed_types.size() == 1 ) {
		pack_task->nonconst_residue_task( ligand_residue_id ).restrict_to_repacking();
		return;
	}
	// else
	for ( core::Size j = 1; j <= allowed_types.size(); ++j ) {
		if ( allowed_types[j]->name() == this_residue.name() ) continue; // already in the task's list
		///TODO figure out why this is nonconst.  Perhaps it could be const
		pack_task->nonconst_residue_task( ligand_residue_id ).allow_noncanonical_aa( allowed_types[j]->name() );
	}
}

std::string ProtLigEnsemble::get_name() const {
	return mover_name();
}

void ProtLigEnsemble::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("scorefxn", xs_string, "Score function to be used during docking")
		+ XMLSchemaAttribute::required_attribute("qsar_file", xs_string, "File containing protein mutations and ligand piars")
		+ XMLSchemaAttribute::required_attribute("final_score", xs_string, "Score function to be used during minimizing")
		+ XMLSchemaAttribute("distance", xsct_real, "Distance around ligand to repack")
		+ XMLSchemaAttribute("ignore_correlation", xsct_non_negative_integer, "Ignore correlation weight for first few ligands")
		+ XMLSchemaAttribute("ignore_correlation_fter", xsct_non_negative_integer, "Ignore correlation weight for at the end")
		+ XMLSchemaAttribute("cycles", xsct_non_negative_integer, "Number of cycles to run.")
		+ XMLSchemaAttribute("repack_every_Nth", xsct_non_negative_integer, "Perform side chain repacking every Nth cycle.");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "High resolution docking of a cross ensemble of protein and ligands", attlist );
}


std::string
ProtLigEnsembleCreator::keyname() const
{
	return ProtLigEnsemble::mover_name();
}

protocols::moves::MoverOP
ProtLigEnsembleCreator::create_mover() const {
	return protocols::moves::MoverOP( new ProtLigEnsemble );
}

std::string
ProtLigEnsemble::mover_name()
{
	return "ProtLigEnsemble";
}

void ProtLigEnsembleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProtLigEnsemble::provide_xml_schema( xsd );
}

// Non-member functions

//Convert info with bind_data to rank
void info_to_rank(utility::vector1<ProtLigPair_info> & vector)
{
	core::Real tied_rank_sum = 0;
	core::Real tied_rank_avg = 0;
	core::Real tied_start = 0;
	core::Real i,j;
	for ( i=1; i<vector.size(); ++i ) {
		if ( vector[i].bind_data != vector[i+1].bind_data ) {
			vector[i].bind_data = i;
		} else {
			tied_start = i;
			do
			{
				++i;
				tied_rank_sum += i;

			}
			while(i < vector.size() && vector[i].bind_data == vector[i+1].bind_data);

			tied_rank_sum += tied_start;
			tied_rank_avg = tied_rank_sum/(i-tied_start+1);
			for ( j=tied_start; j<=i; ++j ) {
				vector[j].bind_data=tied_rank_avg;
			}
			tied_rank_sum = 0;
		}
	}
	if ( i <= vector.size() ) {
		if ( vector[i].bind_data != vector[i-1].bind_data ) vector[i].bind_data = i;
	}

}


//Compare vector by binding data stored in PDB pair info
bool sort_by_binding(const ProtLigPair_info & left_pair, const ProtLigPair_info & right_pair)
{

	//Put no binding data to the larger side so it's at end of vector

	if ( left_pair.has_bind == true && right_pair.has_bind == false ) {
		return true;
	} else if ( left_pair.has_bind == false && right_pair.has_bind == true ) {
		return false;
	} else if ( left_pair.has_bind == false && right_pair.has_bind == false ) {
		//Neither has binding data. Just say it's in right order
		return true;
	} else {
		//Both has data. Return true if right pair is greater
		return left_pair.bind_data < right_pair.bind_data;
	}
}


} //namespace ligand_docking
} //namespace protocols
