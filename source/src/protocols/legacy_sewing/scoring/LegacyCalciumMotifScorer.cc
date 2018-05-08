// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LegacyCalciumMotifScorer.cc
///
/// @brief
/// @author Sharon Guffy

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyCalciumMotifScorer.hh>

//Core headers
#include <core/types.hh>

//see if necessary
#include <core/chemical/ChemicalManager.hh>
//probably necessary
#include <core/pose/motif/reference_frames.hh>
//secondary structure
#include <core/scoring/dssp/Dssp.hh>
//Probably necessary--see what functions are in here, might be useful
#include <core/scoring/motif/util.hh>

//Utility headers
#include <basic/Tracer.hh>
//maybe necessary
#include <numeric/xyzTransform.hh>
//output zstream--see what we're using it for
#include <utility/io/ozstream.hh>

namespace protocols {
namespace legacy_sewing {
namespace scoring {

static basic::Tracer TR("protocols.legacy_sewing.scoring.LegacyCalciumMotifScorer");

LegacyCalciumMotifScorer::LegacyCalciumMotifScorer():
	legacy_sewing::scoring::LegacyMotifScorer()
{
	//set any private variables (possibly finding calcium motif?
}

///@details Use the negative normalized motif score
core::Real
LegacyCalciumMotifScorer::score(
	legacy_sewing::AssemblyCOP assembly
) {
	return -1.0 * full_motif_score(assembly);
}






//METHODS TO CHECK
//Assembly methods
//assembly.starting_node() tells first node in assembly (returns a ModelNode)--first ModelNode in the first entry of edge_path_
//assembly.ending_node() tells the last node in assembly (returns a ModelNode)--last ModelNode in the last entry of edge_path_
//assembly.segments() gives a vector1 of SewSegment objects
//assembly.get_model_segments(int model_id) takes a model_id (an integer) and returns a utility::vector1 of SewSegment objects (maybe useful to get both segments from our starting node)
//assembly.model_ids() returns a std::Set of core::Size (model ids in the assembly)
//assembly.disconnected_segments() returns a vector1 of core::Size (segment ids) that do not have another segment immediately after them
//assembly.path() gives a vector1 of std::pairs of core::Size
//Iterators::
//assembly_begin(), assembly_end()

//Data:
//edge_path_ is a vector1 of tuples (each has 2 raw ModelNode pointers and a raw HashEdge pointer)
//path_ is a vector1 of std::pair of core::Size (path followed to make the assembly, each core::Size is a model_id)
//all_segments_ is a vector1 of vector1 of SewSegment--complete history of the assembly
//segments_ is a SegmentGraph, tells what segments are in the assembly
//partner_pose_ (I won't use this)
//std::set<core::Size> available_nodes_ tells nodes in the assembly that haven't been built off of

//ModelNode methods
//segment_ids() gives a std::set of core::Size

//SegmentGraph--wrapper around a vector1 of SewSegments
//has_next() returns forward_connections_.find( index )->second;
//so if not has_next, then this segment does not have another segment that comes sequentially after it in the assembly



//****model_id for the appended model is -1, so we want all segments with model_id of -1
//AppendAssemblyMover stores a map  models_ that has a std::pair of the model_id (-1) and the Model object
//Then it calls hash_pdb_model on the Model object (which knows its ID and its segment ids) and calls the apply() from the parent MonteCarloAssemblyMover (passing in the pose--presumably it also knows all the private data for this class (models_ is a std::map of model_id to Model object)
//***Look for the generate_assembly() function--see what it knows

//**Ask tim about adding another add_scorer statement in the Assembly Mover so this score is taken into account (or at least reported)





///@details use Will's Motif score to calculate the motif score for interactions between
///a given segment and segments from other models. Divide by total number of segments
core::Real
LegacyCalciumMotifScorer::full_motif_score(
	legacy_sewing::AssemblyCOP assembly
) {

	//First, we need to get just the segments with a model_id of -1
	utility::vector1< legacy_sewing::SewSegment > calcium_site_segments = assembly->get_model_segments( legacy_sewing::POSE_SEGMENT );
	//We also want all segments that are chimeras and have a parent that comes from model_id -1
	utility::vector1< legacy_sewing::SewSegment > chimera_segments = assembly->get_model_segments( legacy_sewing::CHIMERA_SEGMENT );
	//Loop through the chimera segments. If that segment has a parent with model_id of -1, add it to the calcium_site_segments.
	for ( core::Size chimera_index = 1; chimera_index <= chimera_segments.size(); ++chimera_index ) {
		//Now we need to loop through chimera_segments[ chimera_index ].parent_segments_
		for ( core::Size parent_index = 1; parent_index <= chimera_segments[ chimera_index ].parent_segments_.size(); ++parent_index ) {
			if ( chimera_segments[ chimera_index ].parent_segments_[ parent_index ].first == legacy_sewing::POSE_SEGMENT ) {
				//We now know this segment is a chimera, and we can stop looking at this particular segment. Otherwise do nothing.
				calcium_site_segments.push_back( chimera_segments[ chimera_index ] );
				break;
			}
		}
	}

	//If this vector1 is empty, then the final score should just be 0

	core::Real score = 0.0;
	//Do I need to store the scores individually to promote balance? Or is a total good enough?
	//There should be some sort of a penalty for a helical segment having no interactions

	utility::vector1<legacy_sewing::SewSegment> const segments = assembly->segments();
	//We'll want to evaluate interaction of all segments in calcium_site_segments with all segments that aren't in calcium_site_segments
	for ( core::Size cal_seg = 1; cal_seg <= calcium_site_segments.size(); ++cal_seg ) {
		//Score this segment's direct interactions
		core::Real cal_seg_score = 0;
		core::Size cal_seg_counter = 0;
		//What are we counting? The number of residue pairs that are being scored for this segment
		for ( core::Size res_cal = 1; res_cal <= calcium_site_segments[ cal_seg ].residues_.size(); ++res_cal ) {
			numeric::xyzTransform< core::Real > cal_stub = get_stub( calcium_site_segments, cal_seg, res_cal );
			char cal_ss = calcium_site_segments[ cal_seg ].dssp_;
			char cal_aa = res_type_set_->name_map( calcium_site_segments[ cal_seg ].residues_[ res_cal ].residue_type_ ).name1();
			//Loop through all segments
			for ( core::Size seg_num = 1; seg_num <= segments.size(); ++seg_num ) {

				//If this segment (or one matching it) is in calcium_site_segments, continue
				//== operator is defined, so we're in good shape in that regard--can use the has_value method of vector1

				//NOTE: THIS WON'T WORK WELL WITH MULTIPLE CALCIUM SITES
				if ( calcium_site_segments.has_value( segments[ seg_num ]  ) ) {
					continue;
				}
				//Otherwise, evaluate the score residue by residue
				for ( core::Size res_seg = 1; res_seg <= segments[ seg_num ].residues_.size(); ++res_seg ) {
					numeric::xyzTransform<core::Real> seg_stub = get_stub(segments, seg_num, res_seg);
					char seg_ss = segments[ seg_num ].dssp_;
					char seg_aa = res_type_set_->name_map(segments[ seg_num ].residues_[res_seg].residue_type_).name1();

					cal_seg_score += get_score(cal_stub, cal_ss, cal_aa, seg_stub, seg_ss, seg_aa);
					++cal_seg_counter;
					//What should we normalize by? The number of interactions that we're scoring (using counter, as here) seems reasonable
					//But I feel like there should be some penalty if the final score for any helical calcium segment is zero.
					//Or maybe instead a sizeable bonus for each one that is nonzero.
				}//end res_seg
			}//end seg_num
		}//end res_cal


		//Now is where we'll want to add some normalized version of cal_seg_score to score
		//Things to consider:
		//Number of interactions scored for this segment (cal_seg_counter)
		//Is the score for this segment nonzero? Add a bonus (+1 to score after normalization)
		if ( cal_seg_counter == 0 ) {
			//In this case, the score will still be zero
			continue;
		}
		//Otherwise, normalize the current score by the number of interactions and add a bonus if appropriate
		score += cal_seg_score / cal_seg_counter;
		if ( cal_seg_score > 0 ) {
			score += 1;
		}

	}//end cal_seg

	// if(counter == 0) { return score; }
	return score;
}

} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace
