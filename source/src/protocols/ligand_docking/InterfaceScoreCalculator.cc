// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/InterfaceScoreCalculator.cc
/// @brief
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculatorCreator.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <protocols/ligand_docking/ligand_scores.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer InterfaceScoreCalculator_tracer( "protocols.ligand_docking.ligand_options.InterfaceScoreCalculator", basic::t_debug );

std::string
InterfaceScoreCalculatorCreator::keyname() const
{
	return InterfaceScoreCalculatorCreator::mover_name();
}

protocols::moves::MoverOP
InterfaceScoreCalculatorCreator::create_mover() const {
	return protocols::moves::MoverOP( new InterfaceScoreCalculator );
}

std::string
InterfaceScoreCalculatorCreator::mover_name()
{
	return "InterfaceScoreCalculator";
}

/// @brief
InterfaceScoreCalculator::InterfaceScoreCalculator():
		Mover("InterfaceScoreCalculator"),
		chains_(),
		native_(/* NULL */),
		score_fxn_(/* NULL */),
		normalization_function_(/* NULL */),
		compute_grid_scores_(true),
		prefix_("")
{}

InterfaceScoreCalculator::InterfaceScoreCalculator(InterfaceScoreCalculator const & that):
	    //utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chains_(that.chains_),
		native_(that.native_),
		score_fxn_(that.score_fxn_),
		normalization_function_(that.normalization_function_),
		compute_grid_scores_(that.compute_grid_scores_),
		prefix_(that.prefix_)
{}

InterfaceScoreCalculator::~InterfaceScoreCalculator() {}

protocols::moves::MoverOP InterfaceScoreCalculator::clone() const {
	return protocols::moves::MoverOP( new InterfaceScoreCalculator( *this ) );
}

protocols::moves::MoverOP InterfaceScoreCalculator::fresh_instance() const {
	return protocols::moves::MoverOP( new InterfaceScoreCalculator );
}

std::string InterfaceScoreCalculator::get_name() const{
	return "InterfaceScoreCalculator";
}

void InterfaceScoreCalculator::chains(std::vector<std::string> const & chains)
{
	chains_ = chains;
}


void InterfaceScoreCalculator::score_fxn(core::scoring::ScoreFunctionOP const & score_fxn)
{
	score_fxn_ = score_fxn;
}


/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
InterfaceScoreCalculator::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "InterfaceScoreCalculator" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chains") ) throw utility::excn::EXCN_RosettaScriptsOption("'InterfaceScoreCalculator' requires 'chains' tag (comma separated chains to dock)");

	std::string const chains_str = tag->getOption<std::string>("chains");
	chains_= utility::string_split(chains_str, ',');

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) throw utility::excn::EXCN_RosettaScriptsOption("'HighResDocker' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_name);
	assert(score_fxn_);

	if (tag->hasOption("native") ){
		std::string const & native_str= tag->getOption<std::string>("native");
		utility::vector1<std::string> natives_strs= utility::string_split(native_str, ',');
		std::string natives_str = utility::join(natives_strs, " ");

		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*native_, natives_str);
	}
	else if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user()){
		std::string const & native_str= basic::options::option[ basic::options::OptionKeys::in::file::native ]().name();
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*native_, native_str);
	}
	if(tag->hasOption("normalize")) {
		std::string const & normalization_mode = tag->getOption<std::string>("normalize");
		normalization_function_ = protocols::qsar::scoring_grid::get_score_normalization_function(normalization_mode);
	}
	compute_grid_scores_ = tag->getOption<bool>("compute_grid_scores", true);

	prefix_ = tag->getOption<std::string>("prefix","");


}

void InterfaceScoreCalculator::apply(core::pose::Pose & pose) {
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	protocols::jd2::Job::StringStringPairs string_string_pairs(job->get_string_string_pairs());
	if(string_string_pairs.find("native_path") != string_string_pairs.end())
	{
		std::string native_string(string_string_pairs.find("native_path")->second);
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*native_,native_string);
	}
	add_scores_to_job(pose, job);
	append_ligand_docking_scores(pose, job);
}

void InterfaceScoreCalculator::add_scores_to_job(
	core::pose::Pose & pose,
	protocols::jd2::JobOP job
) const
{
	assert(score_fxn_);
	using namespace core::scoring;

	core::Real const tot_score = score_fxn_->score( pose );

	// Which score terms to use
	typedef utility::vector1<ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for(int i = 1; i <= n_score_types; ++i) {
		ScoreType ii = ScoreType(i);
		if ( score_fxn_->has_nonzero_weight(ii) ) score_types.push_back(ii);
	}

	BOOST_FOREACH(ScoreType score_type, score_types){
		std::string const score_term = name_from_score_type(score_type);
		core::Real const weight = score_fxn_->get_weight(score_type);
		if(prefix_ == "")
		{
			job->add_string_real_pair(score_term,  weight * pose.energies().total_energies()[ score_type ]);
		}else
		{
			job->add_string_real_pair(prefix_+ "_" + score_term,  weight * pose.energies().total_energies()[ score_type ]);
		}
	}
	if(prefix_ == "")
	{
		job->add_string_real_pair(name_from_score_type(core::scoring::total_score), tot_score);
	}
	else
	{
		job->add_string_real_pair(prefix_+"_"+name_from_score_type(core::scoring::total_score), tot_score);
	}

}


/// @brief For multiple ligands, append ligand docking scores for each ligand
void
InterfaceScoreCalculator::append_ligand_docking_scores(
		core::pose::Pose const & after,
		protocols::jd2::JobOP job
) const
{
	BOOST_FOREACH(std::string chain, chains_){
		InterfaceScoreCalculator_tracer.Debug << "appending ligand: "<< chain << std::endl;
		assert( core::pose::has_chain(chain, after));
		if(native_)
		{
			if( !core::pose::has_chain(chain, *native_))
			{
				utility_exit_with_message("The native pose passed to InterfaceScoreCalculator does not have chain " +chain);
			}
		}

		utility::vector1<core::Size> jump_ids= core::pose::get_jump_ids_from_chain(chain, after);
		BOOST_FOREACH(core::Size jump_id, jump_ids){
			if(normalization_function_)
			{
				append_interface_deltas(jump_id,job,after,score_fxn_,prefix_,normalization_function_);
			}else
			{
				append_interface_deltas(jump_id, job, after, score_fxn_,prefix_);
			}
			append_ligand_docking_scores(jump_id, after, job);
		}
	}
}

/// @brief Scores to be output that aren't normal scorefunction terms.
void
InterfaceScoreCalculator::append_ligand_docking_scores(
	core::Size jump_id,
	core::pose::Pose const & after,
	protocols::jd2::JobOP job
) const {
	if( jump_id == 0 || jump_id > after.num_jump() ) {
				utility_exit_with_message("The pose does not have jump number " + utility::to_string( jump_id ) );
	}

	if(native_)
	{
		if( jump_id > native_->num_jump() ) {
					utility_exit_with_message("The native pose does not have jump number " + utility::to_string( jump_id ) );
		}

		append_ligand_travel(jump_id, job, *native_, after,prefix_);
		append_radius_of_gyration(jump_id, job, *native_,prefix_);
		append_ligand_RMSD(jump_id, job, *native_, after,prefix_);
	}

	if(compute_grid_scores_)
	{
		if(normalization_function_ && !protocols::qsar::scoring_grid::GridManager::get_instance()->is_normalization_enabled())
		{
			append_ligand_grid_scores(jump_id,job,after,prefix_,normalization_function_);
		}else
		{
			append_ligand_grid_scores(jump_id,job,after,prefix_);
		}
	}


}


} //namespace ligand_docking
} //namespace protocols
