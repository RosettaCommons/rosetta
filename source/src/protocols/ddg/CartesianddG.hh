// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be

/// @file
/// @brief
/// @author Brandon Frenz
/// @author Frank DiMaio
/// @author Hahnbeom Park

#ifndef _INCLUDED_protocols_ddg_CartesianddG_hh_
#define _INCLUDED_protocols_ddg_CartesianddG_hh_


#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <protocols/hybridization/util.hh>
#include <protocols/hybridization/CartesianSampler.fwd.hh>

#include <utility/io/ozstream.hh>
#include <utility/json_utilities.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace ddg {
namespace CartesianddG {

#ifdef _NLOHMANN_JSON_ENABLED_

class MutationSet{
public:
	MutationSet(utility::vector1<core::Size> resnums, utility::vector1<core::chemical::AA> mutations, core::Size iterations){
		assert(resnums.size() == mutations.size());
		resnums_ = resnums;
		mutations_ = mutations;
		iterations_ = iterations;
	}

	//return a list of the residues being mutated in this MutationSet
	utility::vector1<core::Size>
	get_resnums(){ return resnums_; }

	//return a list of the amino acids which the structure is being mutated to.
	utility::vector1<core::chemical::AA>
	get_mutations(){ return mutations_; }

	//add fragset to the list of fragments, resnum should be the residue
	//that corresponds to the start of the fragment
	void
	add_fragments(core::Size resnum, core::fragment::FragSetOP fragset){
		fragments_[resnum] = fragset;
	}

	//return the fragments that match the residue number for this mutationset
	std::map<core::Size,core::fragment::FragSetOP>
	get_fragments(){ return fragments_; }

	//Returns true if any mutated aas are proline
	bool
	has_proline(){
		return mutations_.has_value(core::chemical::aa_pro);
	}

	//Return the residue number that matches the index
	core::Size
	get_resnum(core::Size index){
		return resnums_[index];
	}

	//Return the aa
	core::chemical::AA
	get_aa(core::Size index){
		return mutations_[index];
	}

	//Get the number of remaining iterations
	core::Size
	iterations(){ return iterations_; };

	//Subtract n_sub from the remaining iterations
	void
	subtract_iterations(core::Size n_sub ){
		iterations_ = iterations_-n_sub;
	}

	//Check if the scores have converged
	bool
	is_converged( const core::Size n_to_check, const core::Real cutoff);

	//Return the index of all the prolines.
	utility::vector1<core::Size>
	get_prolines();

	//Returns a sorted vector of pairs with residue numbers and 3 letter codes
	utility::vector1<std::pair<core::Size,std::string>>
	get_mutation_pairs();

	//Generate the tag string based on which residues are being mutated.
	std::string
	generate_tag();

	//Add a score to the score vector
	void
	add_score(core::Real score){
		scores_.push_back(score);
	}

    bool
    is_same(MutationSet otherset){
        if(get_mutation_pairs() == otherset.get_mutation_pairs()){
            return true;
        }else{
            return false;
        }
    }

    void
    add_wildtypes(core::pose::Pose & pose){
        for( core::Size res : resnums_ ){
            wild_types_.push_back(pose.residue(res).name1());
        }
    }

    //@brief returns a json object with the data stored in this mutation set.
    nlohmann::json
    to_json(core::pose::Pose & pose);

private:
	core::Size iterations_;
	utility::vector1<core::Size> resnums_;
	utility::vector1<core::chemical::AA> mutations_;
    utility::vector1<char> wild_types_; //One letter code of the matching wild type residues;
	std::map<core::Size,core::fragment::FragSetOP> fragments_;//These are the fragments used to optimize proline
	utility::vector1<core::Real> scores_; //Scores for this mutation set
};

utility::vector1<core::Size>
find_neighbors(
	MutationSet mutations,
	core::pose::Pose const & pose,
	core::Real const heavyatom_distance_threshold
);

utility::vector1<core::Size>
find_neighbors_directional(
	MutationSet mutations,
	core::pose::Pose const & pose,
	core::Real const K
);


/// @brief Read the mutfile and parse generate the mutset vector
utility::vector1<MutationSet>
read_in_mutations(const std::string filename, core::pose::Pose & pose);

void
sample_fragments(core::pose::Pose & pose,
	MutationSet & mutations,
	core::scoring::ScoreFunctionOP sf,
	const core::Size bbnbrs,
	const core::Size ncycles
);

void
optimize_structure(
	MutationSet mutations,
	core::scoring::ScoreFunctionOP fa_scorefxn,
	core::pose::Pose & pose,
	utility::vector1<core::Size> neighbors,
	const bool flex_bb,
	const bool cartesian = true,
	const core::Size bbnbrs=0
);

/// @brief Optimizes the native structure around all the mutated sites in the mutation set.
void
optimize_native(
	utility::vector1<MutationSet> mutationsets,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP fa_scorefxn,
	const core::Size bbnbrs=0,
	const bool cartesian=true
);

/// @brief Cartesian Minimizer the residues in the set
void
min_pack_min_element(
	core::pose::Pose & pose,
	utility::vector1<core::Size> min_resis,
	core::scoring::ScoreFunctionOP sfxn
);

/// @brief Mutate the one letter sequence to match the mutation set
void
mutate_sequence(std::string & sequence, MutationSet mutations);

/// @brief Return a list of any residues in the mutation set that will involve proline.
utility::vector1<core::Size>
involves_prolines( core::pose::Pose & pose, MutationSet mutations);


/// @brief Mutate the pose to the residues specified in the mutation set.
void
mutate_pose(core::pose::Pose & pose, MutationSet mutations, core::scoring::ScoreFunctionOP score_fxn);


/// @brief Pick fragments for mutationsets involving prolines and store them in the mutation set.
void
pick_fragments(core::pose::Pose & pose, utility::vector1<MutationSet> & mutationsets, const core::Size frag_nbrs);

/// @brief Subtract an iteration from the mutationset that has been completed based on the resnum and aa pair
void
subtract_iterations(utility::vector1<MutationSet> & mutationsets, utility::vector1<std::pair<core::Size,std::string>> finished_mutations, const core::Real score);

/// @brief Trim the two input poses to match the context of the structure and the peptide around the target mutation respectively.
void
extract_element(core::pose::Pose& context_pose, core::pose::Pose &peptide_pose, MutationSet mutations, const core::Size neighbors_to_extract);

/// @brief Extract the neighbors around the mutation(s) and calculate the score ignoring the local electrostatic and hydrogen bonding.
core::Real
extracted_score(core::pose::Pose & pose, MutationSet mutations, core::scoring::ScoreFunctionOP score_fxn, const core::Size n_nbrs_to_extract);

/// @brief Finds the jump that corresponds to the interface
core::Size
find_interface_jump(core::pose::Pose & pose, const core::Size interface_ddg );

/// @brief read the existing files in json format
void
read_existing_json(utility::vector1<MutationSet> & existing_mutsets, const std::string filename, nlohmann::json & results_json, const core::Size iters);

/// @brief Reads the existing results file and decreases the number of iterations in the mutationset
// accordingly.
void
read_existing(std::string filename, utility::vector1<MutationSet> & mutationsets);

// @brief returns a json object of all the score types and their values as found in the pose energy graph.
nlohmann::json
get_scores_as_json(
    core::pose::Pose & pose,
    core::scoring::ScoreFunctionOP score_fxn,
    core::Real total_score);

// @brief write all the results stored in the json object to a file.
void
write_json(const std::string filename, nlohmann::json results_json);

void
run(core::pose::Pose & pose);

#endif //_NLOHMANN_JSON_ENABLED_
}// CartesianDDG
}//ddg
}//protocols

#endif // _INCLUDED_protocols_ddg_CartesianddG_hh_
