// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@unc.edu) and Andrew Leaver-Fay

// Unit Headers
#include <protocols/ligand_docking/HighResEnsemble.hh>
#include <protocols/ligand_docking/HighResEnsembleCreator.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <core/pose/chains_util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>

//Options headers
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <protocols/ligand_docking/ligand_scores.hh>
// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Project Headers
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <algorithm>
#include <sstream>
#include <utility/string_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Scripter Headers
#include <utility/tag/Tag.hh>
//#include <protocols/moves/DataMap.hh>

// Boost Headers

//STL headers
#include <string>


//Auto Headers
#include <utility/vector1.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

#include <protocols/ligand_docking/FinalMinimizer.hh> // AUTO IWYU For FinalMinimizer
#include <protocols/ligand_docking/HighResDocker.hh> // AUTO IWYU For HighResDocker
#include <basic/datacache/DataMap.hh> // AUTO IWYU For DataMap
#include <utility/file/gzip_util.hh> // AUTO IWYU For gzip

namespace protocols {
namespace ligand_docking {

static basic::Tracer high_res_docker_tracer("protocols.ligand_docking.ligand_options.Protocol", basic::t_debug);

///@brief
HighResEnsemble::HighResEnsemble():
	Mover("HighResEnsemble"),
	num_cycles_(0),
	repack_every_Nth_(0),
	chains_(),
	score_fxn_(nullptr),
	movemap_builder_(nullptr),
	resfile_(""),
	final_score_fxn_(nullptr),
	final_movemap_builder_(nullptr),
	correlation_weight_(0),
	exp_ranks_(),
	rosetta_current_scores_(),
	rosetta_lowest_scores_(),
	rosetta_old_scores_(),
	rosetta_current_poses_(),
	rosetta_old_poses_(),
	rosetta_lowest_poses_(),
	rosetta_names_(),
	ligand_chains_()
{
	resfile_.clear();
	// Now use cycles and repack_every_Nth to replicate these options...
	//meiler2006: 50, 8;
	//abbreviated: 5, 4;
	//abbrev2: 6, 3;
}

HighResEnsemble::HighResEnsemble(HighResEnsemble const & /*that*/) = default;

HighResEnsemble::~HighResEnsemble() = default;

protocols::moves::MoverOP HighResEnsemble::clone() const {
	return utility::pointer::make_shared< HighResEnsemble >( *this );
}

bool HighResEnsemble::reinitialize_for_each_job() const { return true; }

protocols::moves::MoverOP HighResEnsemble::fresh_instance() const {
	return utility::pointer::make_shared< HighResEnsemble >();
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
HighResEnsemble::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->getName() != "HighResEnsemble" ) {
		utility_exit_with_message("This should be impossible");
	}

	// cycles and repack_every_Nth
	if ( ! tag->hasOption("cycles") ) utility_exit_with_message("'HighResEnsemble' mover requires cycles tag");
	if ( ! tag->hasOption("repack_every_Nth") ) utility_exit_with_message("'HighResEnsemble' mover requires repack_every_Nth tag");
	num_cycles_= tag->getOption<core::Size>("cycles");
	repack_every_Nth_= tag->getOption<core::Size>("repack_every_Nth");

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) utility_exit_with_message("'HighResEnsemble' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_name);

	/// MoveMapBuilder///
	if ( ! tag->hasOption("movemap_builder") ) utility_exit_with_message("'HighResEnsemble' requires 'movemap_builder' tag");
	std::string movemap_builder_name= tag->getOption<std::string>("movemap_builder");
	movemap_builder_= datamap.get_ptr< MoveMapBuilder >( "movemap_builders", movemap_builder_name);

	/// Resfile ///
	if ( tag->hasOption("resfile") ) {
		resfile_= tag->getOption<std::string>("resfile");
	}

	/// Use Final Minimizer as well, which needs final_score and final_move options
	if ( tag->hasOption("final_score") && tag->hasOption("final_move") ) {
		std::string scorefxn_name= tag->getOption<std::string>("final_score");
		final_score_fxn_= datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_name);

		std::string movemap_builder_name= tag->getOption<std::string>("final_move");
		final_movemap_builder_= datamap.get_ptr< protocols::ligand_docking::MoveMapBuilder >( "movemap_builders", movemap_builder_name);
	} else if ( tag->hasOption("final_score") || tag->hasOption("final_move") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "FinalMinimzer Step requires final_score and final_map tag");
	}

	/// Correlation Weight ///
	if ( ! basic::options::option[ basic::options::OptionKeys::docking::ligand::ligand_ensemble ].user() ) utility_exit_with_message("'HighResEnsemble' requires the docking:ligand:ligand_ensemble option to run");
	correlation_weight_ = basic::options::option[ basic::options::OptionKeys::docking::ligand::ligand_ensemble ]();

	use_rosetta_ranks_ = tag->getOption<bool>("rosetta", false);

	/// Chains for run ///
	if ( ! tag->hasOption("chains") ) utility_exit_with_message("'HighResEnsemble' requires 'chains' tag");
	std::string all_ligands = tag->getOption<std::string>("chains");
	utility::vector1<std::string> ligands_strs = utility::string_split(all_ligands, ',');

	//Get chain ID of ligands/jumps in order
	for ( std::string ligand : ligands_strs ) {
		ligand_chains_.push_back( ligand );
	}

}

void HighResEnsemble::prepare_single_ligand_pose(core::pose::Pose pose, core::Size chain_to_keep)
{

	utility::vector1<core::Size> chain_ids_to_delete;
	core::conformation::ResidueCOPs residues_to_delete;

	for ( std::pair<core::Size, core::Real> ligand : exp_ranks_ ) {
		if ( ligand.first != chain_to_keep ) chain_ids_to_delete.push_back(ligand.first);
	}

	residues_to_delete = core::pose::get_residues_from_chains(pose, chain_ids_to_delete);

	//Reverse iterate since deleting residue will change sequence positions of all residues above delete point, assume ligand residues at end
	std::reverse(residues_to_delete.begin(), residues_to_delete.end());

	for ( core::conformation::ResidueCOP residue_to_delete : residues_to_delete ) {

		// pose.conformation().residues_delete(residue_to_delete->seqpos());
		pose.conformation().delete_residue_slow(residue_to_delete->seqpos());
	}

	//Initial setup of all score vectors
	rosetta_old_scores_.push_back(std::make_pair(chain_to_keep,(*score_fxn_)( pose )));
	rosetta_current_scores_.push_back(rosetta_old_scores_.back());
	rosetta_lowest_scores_.push_back(rosetta_old_scores_.back());

	//Initial setup of all pose vecotrs
	rosetta_old_poses_.push_back(pose);
	rosetta_current_poses_.push_back(pose);
	rosetta_lowest_poses_.push_back(pose);

}

void
HighResEnsemble::apply(core::pose::Pose & pose) {
	debug_assert(num_cycles_ > 0);

	rosetta_old_poses_.clear();
	rosetta_old_scores_.clear();
	rosetta_current_poses_.clear();
	rosetta_current_scores_.clear();
	rosetta_lowest_poses_.clear();
	rosetta_lowest_scores_.clear();

	// Figure out exp ranks:
	for ( std::string const & chain: ligand_chains_ ) {
		core::Size chain_id = core::pose::get_chain_id_from_chain(chain, pose);
		core::conformation::ResidueCOP current_residue = core::pose::get_chain_residues(pose, chain_id)[1];
		if ( use_rosetta_ranks_ ) {
			exp_ranks_.push_back(std::make_pair(chain_id,current_residue->type().get_numeric_property("ROSETTA")));
		} else {
			exp_ranks_.push_back(std::make_pair(chain_id,current_residue->type().get_numeric_property("AFFINITY")));
		}
	}
	// //Convert exp_ranks_ to ranks
	std::sort(exp_ranks_.begin(), exp_ranks_.end(), sort_by_second);  //Sorted by assay now
	vector_to_rank(exp_ranks_);
	std::sort(exp_ranks_.begin(), exp_ranks_.end());  //Restored to ligand in order, necessary for proper deleting


	// Create HighResDocker poses by deleting all ligands except for one in each pose and passing in the parameters, also scores the poses and place into scores
	for ( std::pair<core::Size, core::Real> ligand : exp_ranks_ ) {
		prepare_single_ligand_pose(pose, ligand.first);
	}

	HighResDockerOP high_res_docker ( new HighResDocker(num_cycles_, repack_every_Nth_, score_fxn_, movemap_builder_, resfile_) );
	protocols::moves::MonteCarloOP monte_carlo ( new protocols::moves::MonteCarlo(*score_fxn_, 2.0) );/* temperature, from RosettaLigand paper */

	core::Real correlation_before = 0;
	core::Real correlation_after = 0;
	bool mc_result;
	core::Real score_delta;
	core::Real low_score_so_far; //Temporary storage variable while considering new score
	core::Real score_being_considered; //Current score adjusted with rank correlation

	for ( core::Size cycle = 1; cycle <= num_cycles_; ++cycle ) {

		for ( core::Size i=1; i <= rosetta_current_poses_.size(); ++i ) {

			//Set the MC score (last accepted and lowest) to the score of original poses
			monte_carlo->set_last_accepted(rosetta_old_scores_[i].second);
			monte_carlo->set_lowest(rosetta_lowest_scores_[i].second);

			//Call multiple ligand version of HighResDocker apply....The pose and score for that ligand is now updated in current
			high_res_docker->apply(rosetta_current_poses_[i], rosetta_current_scores_[i].second, ligand_chains_[i], cycle);

			//Store the existing low score aside in case we need it later
			low_score_so_far = rosetta_lowest_scores_[i].second;

			//Calculate the previous correlation for comparison
			correlation_before = qsar_correlation();


			//Put current score into low so we can do the QSAR correlation calculation
			rosetta_lowest_scores_[i].second = rosetta_current_scores_[i].second;

			//Calculate QSAR Correlation
			correlation_after = qsar_correlation();

			//Restore original low score
			rosetta_lowest_scores_[i].second = low_score_so_far;

			//Adjust Score
			score_being_considered = rosetta_current_scores_[i].second - (correlation_after - correlation_before)*correlation_weight_; //This formula assumes that positive correlation is preferred, change sign if this is not the case (with -log K for example)
			std::cout << "Adjusted Current Score is" << score_being_considered << "\n";


			//Find difference (CHECK SUBTRACTION)
			score_delta = score_being_considered - rosetta_old_scores_[i].second;
			std::cout << "Score Delta is" << score_delta << "\n";

			//MC Accept or Reject
			mc_result = monte_carlo->boltzmann(score_being_considered);

			//If accepted, copy current pose into old pose and update score. Check if it beats low and update if appropriate
			if ( mc_result ) {
				std::cout << "Monte Carlo Accepted" << "\n";
				rosetta_old_poses_[i] = rosetta_current_poses_[i];
				rosetta_old_scores_[i].second = rosetta_current_scores_[i].second;

				if ( rosetta_current_scores_[i].second < rosetta_lowest_scores_[i].second ) {
					std::cout << "Beat lowest score" << "\n";
					rosetta_lowest_scores_[i].second = rosetta_current_scores_[i].second;
					rosetta_lowest_poses_[i] = rosetta_current_poses_[i];
				}
			} else {
				//if reject, restore the old pose/score into current
				std::cout << "rejected" << "\n";
				rosetta_current_poses_[i] = rosetta_old_poses_[i];
				rosetta_current_scores_[i].second = rosetta_old_scores_[i].second;
			}

		}

		correlation_after = qsar_correlation();
		std::cout << "Correlation For Cycle " << cycle << " is " << correlation_after;
	}

	//recover low pose, scores, and correlation
	rosetta_current_poses_ = rosetta_lowest_poses_;
	rosetta_current_scores_ = rosetta_lowest_scores_;
	correlation_after = qsar_correlation();

	//Final minimized if being used (Just operates on lowest poses to save time)

	if ( final_score_fxn_ && final_movemap_builder_ ) {
		FinalMinimizerOP final_min( new FinalMinimizer(final_score_fxn_, final_movemap_builder_) );

		for ( core::Size i=1; i<=rosetta_lowest_poses_.size(); ++i ) {
			final_min->apply(rosetta_lowest_poses_[i]);
			rosetta_lowest_scores_[i].second=(*final_score_fxn_)(rosetta_lowest_poses_[i]);
		}

		//update scores and/or correlation after final minimization....Rosetta Lowest score and Rosetta Current scores/poses should be the same
		rosetta_current_poses_ = rosetta_lowest_poses_;
		rosetta_current_scores_ = rosetta_lowest_scores_;
		correlation_after = qsar_correlation();
	}

	//Output individual poses, scores and correlation

	//Just uses lowest poses to save time

	protocols::jd2::add_string_real_pair_to_current_job("spearman", correlation_after);

	for ( core::Size i=1; i<=rosetta_lowest_poses_.size(); ++i ) {

		std::string chain = ligand_chains_[i];

		// core::Size jump = core::pose::get_jump_id_from_chain(chain, rosetta_lowest_poses_[i]);

		core::Size nstruct = protocols::jd2::current_nstruct_index();
		std::string tag( chain + "_" + utility::to_string(nstruct) + ".pdb" );

		//output individual poses
		rosetta_lowest_poses_[i].dump_scored_pdb(tag, *final_score_fxn_);

		if ( basic::options::option[ basic::options::OptionKeys::out::pdb_gz ]() ) {
			utility::file::gzip( tag, true );
		}

		//Add interface delta scores, hopefully to overall pose output
		std::map< std::string, core::Real > ligand_scores( get_interface_deltas( chain, rosetta_lowest_poses_[i], final_score_fxn_, "" ) );
		ligand_scores[ chain ] = rosetta_lowest_scores_[i].second;

		for ( auto const & entry : ligand_scores ) {
			protocols::jd2::add_string_real_pair_to_current_job( entry.first, entry.second );
		}
	}

	std::map<std::string, core::Real> job_outputs = protocols::jd2::get_string_real_pairs_from_current_job();
	core::Real mean_interface_score = 0;

	//Calculate average interface energy from job output
	for ( core::Size i=1; i<=ligand_chains_.size(); ++i ) {
		std::stringstream ss;
		std::string name_of_term;
		std::string chain = ligand_chains_[i];

		ss << "interface_delta_" << chain;
		ss >> name_of_term;

		mean_interface_score = mean_interface_score + job_outputs[name_of_term];
	}

	mean_interface_score = (mean_interface_score / (core::Real)(ligand_chains_.size()));

	protocols::jd2::add_string_real_pair_to_current_job("mean_interface", mean_interface_score);

	//Set pose equal to first pose (temporary benchmark purposes)
	pose = rosetta_lowest_poses_[1];


}

core::Real HighResEnsemble::qsar_correlation()
{
	//If only one, then qsar = 0
	if ( rosetta_lowest_scores_.size() == 1 ) {
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

//void HighResEnsemble::get_ranks_from_file(core::pose::Pose const pose, std::string filename)
//{
// utility::io::izstream infile;
// infile.open(filename.c_str(),std::ifstream::in);
// utility::json_spirit::mValue qsar_data;
// utility::json_spirit::read(infile,qsar_data);
// infile.close();
//
// //Input File Format:
// /*
// [
//     {
//         "name" : "compound name string",
//         "chain" : "X",
//         "assay" : 5
//     }
// ]
//  */
//
// utility::json_spirit::mArray qsar_ranks = qsar_data.get_array();
// for(utility::json_spirit::mArray::iterator counter = qsar_ranks.begin(); counter != qsar_ranks.end();++counter)
// {
//  utility::json_spirit::mObject file_data(counter->get_obj());
//
//  std::string name = file_data["name"].get_str();
//  std::string chain = file_data["chain"].get_str();
//  core::Real assay = file_data["assay"].get_real();
//
//  rosetta_names_.push_back(name);
//
//  core::Size chain_id = core::pose::get_chain_id_from_chain(chain, pose);
////  qsar_jumps_.push_back(core::pose::get_jump_id_from_chain_id(chain_id, pose));
//
//
//  core::conformation::ResidueCOP current_residue = (core::pose::get_chain_residues(pose, chain_id))[1];
//  exp_ranks_.push_back(std::make_pair(chain_id,assay));
//
//
// // qsar_chars_.push_back(core::pose::get_chain_from_chain_id(qsar_chains_.back(), pose));
//
//  std::cout << "\n compound name is: " << rosetta_names_.back();
//  std::cout << "\n assay value is: " << exp_ranks_.back().second;
//  std::cout << "\n chain ID is: " << exp_ranks_.back().first;
// }
//
//}

std::string HighResEnsemble::get_name() const {
	return mover_name();
}

void HighResEnsemble::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("scorefxn", xs_string, "Score function to be used during docking")
		+ XMLSchemaAttribute::required_attribute("chains", xs_string, "Ligand chains, specified as the PDB chain IDs")
		+ XMLSchemaAttribute::required_attribute("movemap_builder", xs_string, "Name of a previously defined MoveMapBuilder for Docking phase.")
		+ XMLSchemaAttribute::required_attribute("final_move", xs_string, "Name of a previously defined MoveMapBuilder for Minimization phase.")
		+ XMLSchemaAttribute::required_attribute("final_score", xs_string, "Score function to be used during minimizing")
		+ XMLSchemaAttribute("resfile", xs_string, "Name (path to) the resfile.")
		+ XMLSchemaAttribute("cycles", xsct_non_negative_integer, "Number of cycles to run.")
		+ XMLSchemaAttribute::attribute_w_default("rosetta", xsct_rosetta_bool,
		"Use Rosetta scores as stand-in for experimental affinity (testing purposes)", "false")
		+ XMLSchemaAttribute("repack_every_Nth", xsct_non_negative_integer, "Perform side chain repacking every Nth cycle.");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Independent high resolution optimization of each protein-ligand interface in ensemble.", attlist );
}


std::string
HighResEnsembleCreator::keyname() const
{
	return HighResEnsemble::mover_name();
}

protocols::moves::MoverOP
HighResEnsembleCreator::create_mover() const {
	return utility::pointer::make_shared< HighResEnsemble >();
}

std::string
HighResEnsemble::mover_name()
{
	return "HighResEnsemble";
}

void HighResEnsembleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HighResEnsemble::provide_xml_schema( xsd );
}

// Non-member functions

//Convert experimental vector<pair(chain_id, affinity)) to rank
void vector_to_rank(utility::vector1<std::pair<core::Size, core::Real> > & vector)
{
	core::Real tied_rank_sum = 0;
	core::Real tied_rank_avg = 0;
	core::Real tied_start = 0;
	core::Real i,j;
	for ( i=1; i<vector.size(); ++i ) {
		if ( vector[i].second != vector[i+1].second ) {
			vector[i].second = i;
		} else {
			tied_start = i;
			do
			{
				++i;
				tied_rank_sum += i;

			}
			while(i < vector.size() && vector[i].second == vector[i+1].second);

			tied_rank_sum += tied_start;
			tied_rank_avg = tied_rank_sum/(i-tied_start+1);
			for ( j=tied_start; j<=i; ++j ) {
				vector[j].second=tied_rank_avg;
			}
			tied_rank_sum = 0;
		}
	}
	if ( i <= vector.size() ) {
		if ( vector[i].second != vector[i-1].second ) vector[i].second = i;
	}

}


//Calculate the spearman between two vectors
core::Real spearman(
	utility::vector1<std::pair<core::Size, core::Real> > vector_exp,
	utility::vector1<std::pair<core::Size, core::Real> > vector_rosetta
)
{

	// utility::vector1<std::pair<core::Size, core::Real> > vector_exp;
	//
	// if (create_subset){
	//  //Create subset for a partial spearman
	//  for(core::Size i=1; i<= vector_rosetta.size(); ++i)
	//  {
	//   vector_exp.push_back(vector_exp_whole[i]);
	//
	//  }
	// }
	//check two vectors have same size
	if ( vector_exp.size() != vector_rosetta.size() ) {
		std::cout << "\n exp vector:" << vector_exp.size();
		std::cout << "\n rosetta vector:" << vector_rosetta.size();
		utility_exit_with_message("'Number of ligands in exp file is not the same as ligand chains");
	}

	core::Real average = 0.5 * (1 + vector_rosetta.size());
	core::Real num = 0, denomE = 0, denomQ = 0;
	core::Size i;


	//Calculate Correlation coefficient between vector_rosetta and vector_qsar
	for ( i=1; i<= vector_rosetta.size(); ++i ) {
		vector_rosetta[i].second = vector_rosetta[i].second - average;
		vector_exp[i].second = vector_exp[i].second - average;
	}

	for ( i=1; i<=vector_rosetta.size(); ++i ) {
		num = num + (vector_rosetta[i].second * vector_exp[i].second);
		vector_rosetta[i].second = (vector_rosetta[i].second * vector_rosetta[i].second);
		vector_exp[i].second = (vector_exp[i].second * vector_exp[i].second);
		denomE = denomE + vector_rosetta[i].second;
		denomQ = denomQ + vector_exp[i].second;
	}
	return num / (sqrt (denomE * denomQ));
}

//Compare vector by second element (affinity)
bool sort_by_second(std::pair<core::Size, core::Real> left_lig, std::pair<core::Size, core::Real> right_lig)
{
	return left_lig.second < right_lig.second;
}

/// @brief Provide the citation.
void
HighResEnsemble::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationCollectionOP cc(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		"HighResEnsemble", basic::citation_manager::CitedModuleType::Mover
		)
	);

	cc->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "10.1021/acsomega.7b02059" ) );

	citations.add( cc );
}


} //namespace ligand_docking
} //namespace protocols
