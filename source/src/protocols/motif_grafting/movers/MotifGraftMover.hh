// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/movers/MotifGraftMover.hh
/// @brief  Declaration of the MoverCreator class for the MotifGraftCreator
/// @author Daniel-Adriano Silva (dadriano@uw.edu) and Alex Ford (fordas@uw.edu)

//Include macro
#ifndef INCLUDED_protocols_motif_grafting_movers_MotifGraftMover_hh
#define INCLUDED_protocols_motif_grafting_movers_MotifGraftMover_hh

// Project Headers
#include <string>
#include <queue>

//Include Rosetta utilities
#include <utility/vector1.hh>

//Include Rosetta Core Stuff
#include <core/types.hh>

//Include Rosetta numeric
#include <numeric/xyz.functions.hh>

//Include Rosetta XML tag reader
#include <utility/tag/Tag.fwd.hh>

//Include Rosetta Scoring functions
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

//Include Rosetta protocols
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

//Include ObjexxFCL
#include <ObjexxFCL/FArray.all.hh>

// Unit headers

namespace protocols {
namespace motif_grafting {
namespace movers {

/**@brief structure that contains the data of corresponding fragments in the motif and scaffold**/
struct motif2scaffold_indexes
{
	core::Size motifLow;
	core::Size motifHigh;
	core::Size scaffoldLow;
	core::Size scaffoldHigh;
};

/**@brief structure that contains the motif2scaffold_indexes data in a vector and adds fields for fragment matching information**/
struct motif2scaffold_data
{
	utility::vector1< motif2scaffold_indexes > v_indexes;
	numeric::xyzMatrix< core::Real > RotM;
	numeric::xyzVector< core::Real > TvecA;
	numeric::xyzVector< core::Real > TvecB;
	core::Real RMSD;
	core::Real motif_fragments_RMSD;
	core::Size clash_score;
	utility::vector1 < utility::vector1< core::Size > > vvr_hotspots;
	bool b_allow_independent_alignment_per_fragment;
	bool b_full_motif_bb_alignment;
	bool b_graft_only_hotspots_by_sidechain_replacement;
};

//Function used to sort the motif2scaffold_indexes by > of the scaffold indexes
inline bool compare_motif2scaffold_data_by_scaffold_low2high(motif2scaffold_indexes const & a, motif2scaffold_indexes const & b)
{
	//Means that best RMSD (lower) goes first
	return (a.scaffoldHigh < b.scaffoldHigh);
}

// @brief Internal class to store generated motif match results
// Implements support to copy-by-value
// Implements '<' operator to sort results by return priority
class MotifMatch {
public:
	//default constructor
	MotifMatch() {};
	//constructor
	MotifMatch(motif2scaffold_data  data)
	{
		scaffold_fragment_data = data;
		RMSD = data.RMSD;
		motif_fragments_RMSD = data.motif_fragments_RMSD;
		clash_score = data.clash_score;
		v_indexes = data.v_indexes;
		b_allow_independent_alignment_per_fragment = data.b_allow_independent_alignment_per_fragment;
		b_full_motif_bb_alignment = data.b_full_motif_bb_alignment;
		b_graft_only_hotspots_by_sidechain_replacement = data.b_graft_only_hotspots_by_sidechain_replacement;
	};
	//Overwrites <operator
	bool operator <(const MotifMatch& other) const
	{
		if ( b_allow_independent_alignment_per_fragment ) {
			return motif_fragments_RMSD > other.motif_fragments_RMSD;
		}
		return RMSD > other.RMSD;
	};
	//accessor methods
	motif2scaffold_data get_scaffold_fragment_data() const { return scaffold_fragment_data; }
	core::Real get_RMSD() const { return RMSD; }
	core::Real get_motif_fragments_RMSD() const { return motif_fragments_RMSD; }
	core::Real get_clash_score() const { return clash_score; }
	bool b_is_full_motif_bb_alignment() const { return b_full_motif_bb_alignment; }
	std::string get_allow_independent_alignment_per_fragment_mode() const {
		std::stringstream converter;
		converter << b_allow_independent_alignment_per_fragment;
		return converter.str();
	};
	std::string get_full_motif_bb_alignment_mode() const {
		std::stringstream converter;
		converter << b_full_motif_bb_alignment;
		return converter.str();
	};
	std::string get_motif_ranges(core::Size const & ndx_shift) const {
		std::string s_out="";
		for ( core::Size i=1; i<= v_indexes.size(); ++i ) {
			s_out += " " + utility::to_string(ndx_shift + v_indexes[i].motifLow) + "," + utility::to_string(ndx_shift + v_indexes[i].motifHigh);
		}
		return s_out;
	};
	std::string get_scaffold_ranges(core::Size const & ndx_shift) const {
		std::string s_out="";
		for ( core::Size i=1; i<= v_indexes.size(); ++i ) {
			s_out += " " + utility::to_string(ndx_shift + v_indexes[i].scaffoldLow) + "," + utility::to_string(ndx_shift + v_indexes[i].scaffoldHigh);
		}
		return s_out;
	};
	std::string get_scaffold2motif_size_change() const {
		std::string s_out="";
		for ( core::Size i=1; i<= v_indexes.size(); ++i ) {
			//can be negative
			long int value= (v_indexes[i].motifHigh-v_indexes[i].motifLow) - (v_indexes[i].scaffoldHigh-v_indexes[i].scaffoldLow) ;
			s_out += " " + utility::to_string(value);
		}
		return s_out;
	};

private:
	motif2scaffold_data scaffold_fragment_data;
	core::Real RMSD;
	core::Real motif_fragments_RMSD;
	core::Real clash_score;
	utility::vector1< motif2scaffold_indexes > v_indexes;
	bool b_allow_independent_alignment_per_fragment;
	bool b_full_motif_bb_alignment;
	bool b_graft_only_hotspots_by_sidechain_replacement;
};//END class MotifMatch

class MotifGraftMover : public protocols::moves::Mover
{
public:
	/**@brief MotifGraftMover Creator**/
	MotifGraftMover();

	/**@brief MotifGraftMover parameters and options initializer**/
	void init_parameters(
		std::string const & s_contextStructure,
		std::string const & s_motif,
		core::Real  const & r_RMSD_tolerance,
		core::Real  const & r_NC_points_RMSD_tolerance,
		core::Size  const & i_clash_score_cutoff,
		std::string const & s_combinatory_fragment_size_delta,
		std::string const & s_max_fragment_replacement_size_delta,
		std::string const & s_clash_test_residue,
		std::string const & s_hotspots,
		bool              & b_full_motif_bb_alignment,
		bool        const & b_allow_independent_alignment_per_fragment,
		bool        const & b_graft_only_hotspots_by_sidechain_replacement,
		bool        const & b_only_allow_if_N_point_match_aa_identity,
		bool        const & b_only_allow_if_C_point_match_aa_identity,
		bool        const & b_revert_graft_to_native_sequence,
		bool        const & b_allow_repeat_same_graft_output);

	/**@brief MotifGraftMover Destructor**/
	~MotifGraftMover();

	/**@brief Apply mover function**/
	void apply( Pose & ) override;

	/**@brief Iterate over the results to get additional matches in the queue**/
	core::pose::PoseOP get_additional_output() override;

	/**@brief Function used by roseta to create clones of movers**/
	protocols::moves::MoverOP clone() const override;

	/**@brief Header only mover get_name**/
	// XRW TEMP  virtual std::string get_name() const
	// XRW TEMP  {
	// XRW TEMP   return "MotifGraft";
	// XRW TEMP  }

	/** @brief As the name suggests in generates all the permutations of a vector of vectors of pairs (Alex: we should templatize this! Maybe alrready there?)**/
	void permutate_n_vv_of_pairs(
		utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > const & vv_of_pairs,
		utility::vector1< std::pair< core::Size, core::Size > > & buff_combVec,
		core::Size start_index,
		utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > & vv_resulting_permutations);

	/**@brief Generate all the combination of different legths of the motif fragment as requested in combinatory_fragment_size_delta
	** Uses permutate_n_vv_of_pairs to generate the permutations**/
	void generate_combinations_of_motif_fragments_by_delta_variation(
		core::pose::PoseOP const & p_motif_,
		utility::vector1 < std::pair< long int, long int > > const & combinatory_fragment_size_delta,
		utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > & vv_resulting_permutations);

	/**@brief Fuction to parse RosettaScripts XML options**/
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &) override;

	/**@brief Identify all potential matches for the given target scaffold (this is where the motif grafting code is called)**/
	std::priority_queue<MotifMatch> generate_scaffold_matches(
		core::pose::Pose & target_scaffold,
		core::pose::PoseOP & target_motif_,
		core::pose::PoseOP & target_contextStructure_);

	/** @brief Generate pose corresponding to the given match **/
	void generate_match_pose(
		core::pose::Pose & target_pose,
		core::pose::Pose const & contextStructure,
		bool b_revert_graft_to_native_sequence,
		MotifMatch motif_match);

	/**@brief Return a priority queue with the sucessful epigrafts **/
	void get_matching_fragments(
		core::pose::Pose const & target_scaffold,
		core::pose::PoseOP const & target_motif_,
		core::pose::PoseOP const & target_contextStructure_,
		core::Real const & RMSD_tol,
		core::Real const & NC_points_RMSD_tol,
		core::Size const & clash_cutoff,
		std::string const & clash_test_residue,
		utility::vector1 < std::pair< long int, long int > > const & max_fragment_replacement_size_delta_,
		utility::vector1 < std::pair< core::Size, core::Size > > const & combinatory_fragment_size_delta,
		utility::vector1 < utility::vector1< core::Size > > const & vvr_hotspots,
		bool const & b_full_motif_bb_alignment,
		bool const & b_allow_independent_alignment_per_fragment,
		bool const & b_graft_only_hotspots_by_sidechain_replacement,
		bool const & b_only_allow_if_N_point_match_aa_identity,
		bool const & b_only_allow_if_C_point_match_aa_identity,
		std::priority_queue<MotifMatch> & pq);

	/**@brief Functions that takes the scaffold, motif, contextStructure and superposition transform data. Deletes from the supperposition data
	** those transformations that can't pass the clash score**/
	void test_epigraft_and_contextStructure_clashes(
		core::pose::Pose const & p_scaffold,
		core::pose::Pose const & p_motif_,
		core::pose::Pose const & p_contextStructure_,
		core::Size const & clash_cutoff,
		utility::vector1< motif2scaffold_data > & v_m2s_data);

	/**@brief Count the Number of Clashes between two poses*/
	core::Size count_clashes_between_two_poses(
		core::pose::Pose const & p_A,
		core::pose::Pose const & p_B,
		core::Size clash_cutoff);

	/**@brief Function that returns by reference a rotated copy of the pose */
	core::pose::Pose get_rotated_and_translated_pose(
		core::pose::Pose const & p_scaffold,
		numeric::xyzMatrix< core::Real > const & RotM,
		numeric::xyzVector< core::Real > const & TvecA,
		numeric::xyzVector< core::Real > const & TvecB);

	/**@brief returns a pose with two input poses merged (with a jump in-between) and with the PDB info corrected*/
	core::pose::Pose join_two_poses_by_jump(
		core::pose::Pose const & p_A,
		core::pose::Pose const & p_B);

	/**@brief Helper function to stich (epigraft) two poses given a set of indices in pose A and B stored in a motif2scaffold_data structure**/
	core::Real get_clash_score_from_pose(
		core::pose::Pose & p_input,
		core::scoring::ScoreFunctionOP const & scorefxn_);

	/**@brief Helper function to stich (epigraft) two poses given a set of indices in pose A and B stored in a motif2scaffold_data structure**/
	core::pose::Pose stich_motif_in_scaffold_by_indexes_rotation_and_translation(
		core::pose::Pose const & p_scaffold,
		core::pose::Pose const & p_motif_,
		motif2scaffold_data & m2s_dat,
		bool const & skip_motif_extremes);


	/**@brief Performs alignment of the protein BB on the selected aminoacids.
	**Returns the RMSD,
	**Returns by reference the rotation Matrix and Translation Vector,
	**Will fail if both poses are not protein <-This can be fixed by adding a list of the atoms to align to the function, but I am not doing it now.
	**Will fail if the number of residues to align is not the same in the two poses. **/
	core::Real get_bb_alignment_and_transformation(
		core::pose::Pose const & poseA,
		utility::vector1< core::Size > const & positions_to_alignA,
		core::pose::Pose const & poseB,
		utility::vector1< core::Size > const & positions_to_alignB,
		numeric::xyzMatrix< core::Real > & RotM,
		numeric::xyzVector< core::Real > & TvecA,
		numeric::xyzVector< core::Real > & TvecB);

	/**@brief Performs alignment of the protein BB on the selected aminoacids.
	**Returns the RMSD,
	**Returns by reference the rotation Matrix and Translation Vector,
	**Will fail if both poses are not protein <-This can be fixed by adding a list of the atoms to align to the function, but I am not doing it now.
	**Will fail if the number of residues to align is not the same in the two poses. **/
	core::Real get_bb_alignment_and_transformation_wTipsExtraInfo(
		utility::vector1< bool > const & v_isNorC,
		core::pose::Pose const & poseA,
		utility::vector1< core::Size > const & positions_to_alignA,
		core::pose::Pose const & poseB,
		utility::vector1< core::Size > const & positions_to_alignB,
		numeric::xyzMatrix< core::Real > & RotM,
		numeric::xyzVector< core::Real > & TvecA,
		numeric::xyzVector< core::Real > & TvecB,
		utility::vector1< core::Real > & RMSD_tip_elements);

	/**@brief Returns the BB distance of two poses respect to indexes**/
	core::Real get_bb_distance(
		core::pose::Pose const & poseA,
		utility::vector1< core::Size > const & positions_to_alignA,
		core::pose::Pose const & poseB,
		utility::vector1< core::Size > const & positions_to_alignB);


	/** @brief Helper Fortran wrapper to get the aligment of two matrixes as well as the corresponding transform**/
	void superposition_transform(
		core::Size natoms,
		ObjexxFCL::FArray1_double const& weights,
		ObjexxFCL::FArray2_double& ref_coords,
		ObjexxFCL::FArray2_double& coords,
		numeric::xyzMatrix< core::Real > &RotM,
		numeric::xyzVector< core::Real > &TvecA,
		numeric::xyzVector< core::Real > &TvecB);

	/** @brief performs soperposition based on the motif and returns fragments within the RMSD_tol and
	** also returns the Rotation and translation superposition information in the first [1] vector of each fragment**/
	void get_motif_scaffold_superposition_and_RMSD(
		core::pose::Pose const & p_scaffold,
		core::pose::PoseOP const & p_motif_,
		core::Real const & RMSD_tol,
		core::Real const & tip_RMSD_tol,
		utility::vector1 < utility::vector1< core::Size > > const & vvr_hotspots,
		bool const & b_full_motif_bb_alignment,
		bool const & b_allow_independent_alignment_per_fragment,
		bool const & b_graft_only_hotspots_by_sidechain_replacement,
		utility::vector1< motif2scaffold_data > & v_m2s_data);

	/** @brief returns a copy of the pose that replaces all the aminoacids for a single selected aminoacid**/
	core::pose::Pose get_mono_aa_pose_copy(
		core::pose::Pose const & p_input,
		std::string const & aminoacid_code);

	/** @brief returns by reference two vectors of indexes (vv_scaffold_fragments_indexes, v_motif_fragments_indexes)
	** that hold the lower and upper bounds of the fragments. Indeed the corresponding to the scaffold one is a vector
	**of vectors, since each pose_scaffold can have many matches**/
	bool get_fragments_by_CA_distances_and_NCpoints_restrains(
		core::pose::Pose const & p_scaffold,
		core::pose::PoseOP const & p_motif_,
		utility::vector1< utility::vector1 < std::pair< core::Size, core::Size > > > & vv_scaffold_fragments_indexes,
		utility::vector1< std::pair< core::Size, core::Size > > & v_motif_fragments_indexes,
		core::Real const & RMSD_tol,
		utility::vector1 < std::pair< long int, long int > > const & max_fragment_replacement_size_delta,
		utility::vector1< std::pair< core::Size, core::Size > > const & v_motif_fragments_permutation,
		bool const & b_only_allow_if_N_point_match_aa_identity,
		bool const & b_only_allow_if_C_point_match_aa_identity,
		bool const & b_N_point_can_replace_proline,
		bool const & b_C_point_can_replace_proline);

	/** @brief Generates all the discontinuous fragments combinations that are within the tol restriction.
	** Reduce the combinations by matching intra chains distances by pairs.
	** The method is/can be exahustive but fast (i.e. iteratively it test the restrains (tree unfolding) and skips to the next combination once one fails (branch removal) ).
	** CAUTION, Uses self recursion, so use it wisely **/
	void fragments_permutation_test_by_CA_distances(
		core::pose::Pose const & p_scaffold,
		core::pose::PoseOP const & p_motif_,
		utility::vector1< utility::vector1 < std::pair< core::Size, core::Size > > >  const & vv_scaffold_fragments_indexes,
		utility::vector1< std::pair< core::Size, core::Size > > const & v_motif_fragments_indexes,
		core::Real const & RMSD_tol,
		core::Size const & start_motif_num,
		utility::vector1< std::pair< core::Size, core::Size > > & buff_combVec,
		utility::vector1< motif2scaffold_data > & v_m2s_data);

	void parse_my_string_arguments_and_cast_to_globalPrivateSpaceVariables(
		std::string const & s_contextStructure,
		std::string const & s_motif,
		core::Real  const & r_RMSD_tolerance,
		core::Real  const & r_NC_points_RMSD_tolerance,
		core::Size  const & i_clash_score_cutoff,
		std::string const & s_combinatory_fragment_size_delta,
		std::string const & s_max_fragment_replacement_size_delta,
		std::string const & s_clash_test_residue,
		std::string const & s_hotspots,
		bool              & b_full_motif_bb_alignment,
		bool        const & b_allow_independent_alignment_per_fragment,
		bool        const & b_graft_only_hotspots_by_sidechain_replacement,
		bool        const & b_only_allow_if_N_point_match_aa_identity,
		bool        const & b_only_allow_if_C_point_match_aa_identity,
		bool        const & b_revert_graft_to_native_sequence,
		bool        const & b_allow_repeat_same_graft_output);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	core::pose::PoseOP gp_p_contextStructure_;
	core::pose::PoseOP gp_p_motif_;
	core::Real         gp_r_RMSD_tolerance_;
	core::Real         gp_r_NC_points_RMSD_tolerance_;
	core::Size         gp_i_clash_score_cutoff_;
	utility::vector1 < std::pair< core::Size, core::Size > >  gp_vp_combinatory_fragment_size_delta_;
	utility::vector1 < std::pair< long int, long int > >      gp_vp_max_fragment_replacement_size_delta_;
	std::string        gp_s_clash_test_residue_;
	utility::vector1 < utility::vector1< core::Size > >       gp_vvr_hotspots_;
	bool               gp_b_full_motif_bb_alignment_;
	bool               gp_b_allow_independent_alignment_per_fragment_;
	bool               gp_b_allow_repeat_same_graft_output_;
	bool               gp_b_is_first_run_;
	bool               gp_b_graft_only_hotspots_by_sidechain_replacement_;
	bool               gp_b_only_allow_if_N_point_match_aa_identity_;
	bool               gp_b_only_allow_if_C_point_match_aa_identity_;
	bool               gp_b_revert_graft_to_native_sequence_;

	core::pose::PoseOP gp_p_target_pose_; //Swap space for our input pose
	std::priority_queue<MotifMatch> motif_match_results_;
}; //END class MotifGraftMover

}//END namespace movers
}//END namespace motif_grafting
}//END namespace protocols

#endif
