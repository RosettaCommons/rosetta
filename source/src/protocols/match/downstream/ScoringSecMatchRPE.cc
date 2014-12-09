// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/ScoringSecMatchRPE.cc
/// @brief
/// @author Kui K. Chan, kuichan@u.washington.edu, Oct, 2009

// Unit headers
#include <protocols/match/downstream/ScoringSecMatchRPE.hh>

// Package headers
// AUTO-REMOVED #include <stdio.h>

// Project headers
// AUTO-REMOVED #include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <basic/basic.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// for debug usage
// AUTO-REMOVED #include <protocols/match/Matcher.hh>
#include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <protocols/match/Hit.hh> // REQUIRED FOR WINDOWS

// Numeric headers
// AUTO-REMOVED #include <numeric/constants.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>

// C++ headers
#include <list>

#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <time.h>
// AUTO-REMOVED #include <stdlib.h>
// AUTO-REMOVED #include <core/io/pdb/file_data.hh>
#include <basic/Tracer.hh>

#if (defined _WIN32) && (!defined WIN_PYROSETTA)
// AUTO-REMOVED #include <windows.h>

#include <utility/vector1.hh>

#endif

namespace protocols {
namespace match {
namespace downstream {
static thread_local basic::Tracer TR( "core.protocols.match.downstream" );
ScoringSecMatchRPE::~ScoringSecMatchRPE() {}


void
ScoringSecMatchRPE::setPose(
core::pose::Pose const & ref_pose )
{
	ref_pose_ = core::pose::PoseOP( new core::pose::Pose(ref_pose) );
}

/*
ScoringSecMatchRPE::ScoringSecMatchRPE(
std::string const & s_in, std::string const & pdb_file)
{
	core::pose::Pose pose(pdb_file);
	this->(s_in, pose);
}
*/

/**
	* 1) I have not implement long range two bodies constraint.
	* 2) I will not check "weight" parameter for two bodies constraint.
	* 3) I check two bodies constraint for cutoff parameter.
  * 4) All inputs constraints are handled in the constructor.
	**/
ScoringSecMatchRPE::ScoringSecMatchRPE(
std::string const & s_in, core::pose::Pose const & ref_pose)
{
	//cd_2b_pose_ = new core::pose::Pose ();
	cd_2b_pose_ = core::pose::PoseOP( new core::pose::Pose(ref_pose) );

	//scoringEval_counter = 0;
	secmatch_scotypes_cutoff_.clear();
	secmatch_value_cutoff_.clear();
	ref_pose_ = core::pose::PoseOP( new core::pose::Pose(ref_pose) );
	sfxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
	cutoff_flag_ = false;
	cutoff_scoreType_flag_ = false;
	longRange_ = false;
	shortRange_ = false;

	//parse input into line
	utility::vector1<std::string> vLineString = utility::string_split( s_in, '\n' );
	TR << s_in << std::endl;
	//keyword in the CONSTRAINT file
	std::string const keywordSCORING ("SCORING_SECMATCH::");
	std::string const keywordFilename ("weights_file:");
	std::string const keywordWeight ("weights:");
	std::string const keywordCutOff ("cutoff:");
	std::string const keywordTotal_score ("total_score");

	for(utility::vector1< std::string >::iterator it_line=vLineString.begin(), end_line =vLineString.end();
			it_line != end_line; ++it_line ) {

		//if line contain CONSTRINAT SCORING_SECMATCH keyword
		if ((*it_line).find(keywordSCORING) != std::string::npos){

			//parse each line into words
			utility::vector1<std::string> vString = utility::string_split( *it_line );
			for (utility::vector1< std::string >::iterator it=vString.begin(), end = vString.end(); it != end; ++it){

				//weigh file ex. standard.wts
				if ((*it).find(keywordFilename) != std::string::npos){
					++it;
					//TR << "Filename:" + *it +"!" << std::endl;
					sfxn_->add_weights_from_file(*it);

				//Weight for each of the score type
				} else if ((*it).find(keywordWeight)!= std::string::npos){
					while ( it != end ){
						if (core::scoring::ScoreTypeManager::is_score_type(*it)){

							core::scoring::ScoreType sType = core::scoring::ScoreTypeManager::score_type_from_name (*it);
							//secmatch_scotypes_.push_back( sType );
							//TR << "Weight Input from ScoringSecMatchRPE:" + *it +"!" << std::endl;
							++it;
							core::Real weight = std::atof((*it).c_str());
							//TR << "Weight Input from ScoringSecMatchRPE:" + *it +"!" << std::endl;
							sfxn_->set_weight( sType, weight );
						}
						++it;
					}

				//CutOff
				} else if ((*it).find(keywordCutOff)!= std::string::npos){

					while ( it != end ){

						//CutOff: total_score
						if ((*it).find(keywordTotal_score) != std::string::npos){

							++it;
							total_score_cutoff_ = std::atof((*it).c_str());
							//TR << "Total_score constraint:" + *it +"!" << std::endl;
							cutoff_flag_ = true;
							TR << "TotalScore is the sum of context independent two bodies and context dependent two bodies" << std::endl;

						//CutOff: scoretype
						} else if (core::scoring::ScoreTypeManager::is_score_type(*it)){
							core::scoring::ScoreType sType = core::scoring::ScoreTypeManager::score_type_from_name (*it);

							//if scoretype is not 2b type, we skip and continue
							//if (sType > core::scoring::n_shortranged_2b_score_types)
							if (!ScoringSecMatchRPE::check2bsc( sType, sfxn_->core::scoring::ScoreFunction::get_weight( sType ))){
									std::string errorStr = "This constraint, " + *it +", is not a two bodies scoretype.";
									utility_exit_with_message (errorStr);
							}

							//TR << "CutOff Input from ScoringSecMatchRPE sType:" + *it +"!" << std::endl;
							++it;
							core::Real sc_cutoff ( std::atof((*it).c_str()) );
							//TR << "CutOff Input from ScoringSecMatchRPE cutoff:" << sc_cutoff << "!" << std::endl;
							secmatch_scotypes_cutoff_.push_back( sType );
							secmatch_value_cutoff_.push_back( sc_cutoff );
							cutoff_scoreType_flag_ = true;
						}
						++it;
					}
				}
				if (it == end) break;
			}
		}
  }


	/**check specified weight for corresponding cutoff value
		*if there is cutoff value, but no weight specified, then error
		**/
	utility::vector1< core::scoring::ScoreType > secmatch_scotypes_temp = sfxn_->get_nonzero_weighted_scoretypes();
	for( utility::vector1< core::scoring::ScoreType>::const_iterator cut_it = secmatch_scotypes_cutoff_.begin();
		cut_it != secmatch_scotypes_cutoff_.end(); ++cut_it ){
		bool noCutOff = true;

		for( utility::vector1< core::scoring::ScoreType>::const_iterator sco_it = secmatch_scotypes_temp.begin();
			sco_it != secmatch_scotypes_temp.end(); ++sco_it ){

			//TR << *sco_it << ":" << *cut_it << std::endl;
			if (*sco_it == *cut_it){
				noCutOff = false;
			}
		}
		if (noCutOff){
			TR << "ERROR: There is no specific weight for the following scoretype " << *cut_it << std::endl;
			exit(1);
		}
	}

}

/**
	*Check incoming scoreType is two bodies scoreType.
	*We also setup the long range and short range scoreType flag.
	**/
bool
ScoringSecMatchRPE::check2bsc(
	core::scoring::ScoreType sType,
	core::Real wts
)
{
	using namespace core::scoring::methods;

	core::scoring::methods::EnergyMethodOptionsOP energy_method_options( new core::scoring::methods::EnergyMethodOptions );
	utility::vector1< core::Real > wtsV;
	wtsV.push_back(wts);
	energy_method_options->set_method_weights( sType , wtsV );
	core::scoring::methods::EnergyMethodOP method (
		core::scoring::ScoringManager::get_instance()->energy_method( sType, *energy_method_options));

	bool twobodiesterm = false;
	//The switch codes are copied from core::scoring::ScoreFunction.cc (approximately line 1740)
  // PHIL replace these with utility::down_cast when it's working
  switch ( method->method_type() ) {
  case ci_2b:
		twobodiesterm = true;
		shortRange_ = true;
		TR << "ci_2b:" << sType << " weight:" << wts << std::endl;
    break;

  case cd_2b:
		twobodiesterm = true;
		shortRange_ = true;
		TR << "cd_2b:" << sType << " weight:" << wts << std::endl;
    break;

  case ci_1b:
		twobodiesterm = false;
		TR << "ci_1b:" << sType << " weight:" << wts << std::endl;
    break;

  case cd_1b:
		twobodiesterm = false;
		TR << "cd_1b:" << sType << " weight:" << wts << std::endl;
    break;

  case ci_lr_2b:
		twobodiesterm = true;
		longRange_ = true;
		TR << "ci_lr_2b:" << sType << " weight:" << wts << std::endl;
    break;

  case cd_lr_2b:
		twobodiesterm = true;
		longRange_ = true;
		TR << "cd_lr_2b:" << sType << " weight:" << wts << std::endl;
    break;

  case ws:
		twobodiesterm = false;
		TR << "ws:" << sType << " weight:" << wts << std::endl;
    break;

  default:
    utility_exit_with_message( "unrecognized scoreType:"
				+ core::scoring::ScoreTypeManager::name_from_score_type (sType));

  } // switch

	return twobodiesterm;
}

bool
ScoringSecMatchRPE::evaluate_residues(
  core::conformation::Residue const & match_res,
  core::conformation::Residue const & target_res
) const
{
	bool returnVal = false;
	//utility_exit_with_message("I am in evaluate_residue, ScoringSecMatchRPE");
	//evaluate two body short range terms
	if (shortRange_){
		returnVal = ScoringSecMatchRPE::eval_cd_2b_residues(match_res,target_res);
				//utility_exit_with_message ("ScoringSecMatchRPE evaluate_residue");

		//for debug output 101509
		//Output the PDB coordinates of match residue and target residue
		if (returnVal){
		//TR << "ScoringSecMatchRPE: match_res:" << match_res << " target_res:" << target_res
		//			<< " returnVal:" << returnVal << std::endl;


		//debug usage
		std::string randNum;
		std::stringstream out;
		out << numeric::random::random_range(0,INT_MAX);
		randNum = out.str();
#ifdef _WIN32
		#ifndef WIN_PYROSETTA
			Sleep( 10000 );
		#endif
#else
		sleep (10);
#endif
		utility::io::ozstream outstream;
		std::string timeStr("EV_" + randNum + ".run");
		outstream.open( timeStr, std::ios::app );

		//core::pose::Pose pose1 (*ref_pose_);
		//pose1.append_residue_by_jump( match_res, 1 );
		//pose1.append_residue_by_jump( target_res, 1 );
		//pose1.dump_pdb( outstream, "1");

		core::Size temp_counter (1);
		core::io::pdb::dump_pdb_residue( match_res, temp_counter, outstream);
		core::io::pdb::dump_pdb_residue( target_res, temp_counter, outstream);
		outstream.close();
		//debug
		}
	}
	if( !returnVal ) return false;

	//need to implement two body long range terms
	//currently, we are return false and exist program
	if (longRange_)
		returnVal = ScoringSecMatchRPE::eval_longRange_2b_residue(match_res,target_res);
	return returnVal;
}

bool
ScoringSecMatchRPE::eval_longRange_2b_residue(
  core::conformation::Residue const & /*match_res*/,
  core::conformation::Residue const & /*target_res*/
) const
{
	TR << "CAUTION: There is no long range implementation" << std::endl;
	utility_exit_with_message("CAUTION: There is no long range implementation");
	return false;
}

bool
ScoringSecMatchRPE::eval_cd_2b_residues(
  core::conformation::Residue const & match_res,
  core::conformation::Residue const & target_res
) const
{
	core::scoring::EnergyMap emap;


	//debug
	//emap.print();

	//TotalScore
	if (cutoff_flag_){
			core::Real value( 0.0 );
			//debug
			//TR <<"TotalScore emap size:"<< sizeof(emap) << std::endl;
			//TR <<"TatalScore sfxn_->weight():"<< sizeof(sfxn_->weights()) << std::endl;

			//context independent two bodies energy term
			sfxn_->eval_ci_2b( match_res, target_res, *ref_pose_, emap);
			//context dependent two bodies energy term code
			TR << "Caution::The context dependent two bodies energy term is computational expensive and it will slow down matcher." << std::endl;
			ref_pose_->replace_residue(match_res.seqpos(),match_res,false);
			ref_pose_->replace_residue(target_res.seqpos(),target_res,false);
			sfxn_->setup_for_scoring( *ref_pose_ );
			sfxn_->eval_cd_2b (match_res, target_res, *ref_pose_, emap);
			//end context dependent code

			//core::scoring::ScoreFunction::eval_ci_2b ( match_res, target_res, ref_pose_, emap);
			//sfxn_->eval_ci_2b ( match_res, target_res, *ref_pose_, emap);
			value =  emap.dot( sfxn_->weights() );
			TR << "Calculated TotalScore:" << value << std::endl;
			if (value > total_score_cutoff_) return false;
	}

	//debug
	//TR << "match_res:" << match_res << std::endl;
	//TR << "target_res:" << target_res << std::endl;
	//ref_pose_->dump_pdb("ref_pose_.txt","1");


	//Check each score type
	if (cutoff_scoreType_flag_){

			utility::vector1< core::scoring::ScoreType >::const_iterator sco_it = secmatch_scotypes_cutoff_.begin();

			for( utility::vector1< core::Real >::const_iterator real_it = secmatch_value_cutoff_.begin();
						real_it != secmatch_value_cutoff_.end(); ++real_it ){

				//if scoretype is not 2b type, we skip and continue
				//We will skip scoretype that is not 2b type to calculate the total score
				if (*sco_it > core::scoring::n_shortranged_2b_score_types){
					 continue;
				//context independent two bodies energy term
				} else if (*sco_it < core::scoring::n_ci_2b_score_types){
					sfxn_->eval_ci_2b( match_res, target_res, *ref_pose_, emap);

				//context dependent two bodies energy term code
				//*sco_it > core::scoring::n_ci_2b_score_types && *sco_it < core::scoring::n_shortranged_2b_score_types
				} else {

					//context dependent two bodies energy term code
					//This is computational expensive.
					//Need optimaization
					ref_pose_->replace_residue(match_res.seqpos(),match_res,false);
					ref_pose_->replace_residue(target_res.seqpos(),target_res,false);
					sfxn_->setup_for_scoring( *ref_pose_ );
					sfxn_->eval_cd_2b (match_res, target_res, *ref_pose_, emap);
				}
				//end context dependent code

				//debug
				//TR << "ScoreType cutoff:" << *sco_it << ":" << tmp_cutoff << std::endl;
				//TR << "emap[ *sco_it ]:" << emap[ *sco_it ] << " sfxn_->weights()[ *sco_it ]:" << sfxn_->weights()[ *sco_it ] <<std::endl;
				//TR << "sco_it:" << *sco_it << std::endl;

				//core::Real tmp_cutoff ( *real_it );
				core::Real checkCutoff = emap[ *sco_it ] * sfxn_->weights()[ *sco_it ];
				TR << "sco_it:" << *sco_it << " checkCutoff:" << checkCutoff << std::endl;
				if (checkCutoff > *real_it) return false ;
				++sco_it;

			}

	}

	return true;
}



core::scoring::ScoreFunctionCOP
ScoringSecMatchRPE::get_score_function() const
{
  return sfxn_;
}

core::scoring::ScoreFunctionOP
ScoringSecMatchRPE::get_score_function()
{
  return sfxn_;
}

core::Real
ScoringSecMatchRPE::getCutoff() const
{
	return total_score_cutoff_;
}

bool
ScoringSecMatchRPE::require_all_target_residue_atom_coordinates() const
{
	return true;
}

bool
ScoringSecMatchRPE::require_target_atom_coordinate( Size /*target_atom_id*/ ) const
{
	return true;
}

}
}
}

