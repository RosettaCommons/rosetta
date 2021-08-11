// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/util.cc
/// @brief Utility functions for Pose Sewing.
/// @author Frank Teets (frank.teets@proteininnovation.org)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/pose_sewing/util.hh>

#include <protocols/sewing/scoring/MotifScorer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/pose/datacache/CacheableObserver.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/motif/util.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzTransform.hh>
#include <core/conformation/Atom.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <numeric/HomogeneousTransform.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

#include <cmath>

static basic::Tracer TR( "protocols.pose_sewing.util" );


namespace protocols {
namespace pose_sewing {

using core::Size;
using namespace core::select::residue_selector;


void
identify_ss_blocks(std::set<core::select::residue_selector::ResidueSubset> & outset, core::select::residue_selector::ResidueSubset const & selection ){
	//std::set<core::select::residue_selector::ResidueSubset> block_selections;
	core::select::residue_selector::ResidueSubset current_selection( selection.size(), false );

	//Note - does not take into account new chains!!
	bool last_res = false;
	for ( core::Size count_resnum = 1; count_resnum <= selection.size(); ++count_resnum ) {
		if ( selection[count_resnum] ) {
			if ( !last_res ) {
				current_selection = *(new core::select::residue_selector::ResidueSubset( selection.size(), false ));
			}
			current_selection[count_resnum] = true;
		} else {
			if ( last_res ) {
				outset.insert(current_selection);
			}
		}
		last_res = selection[count_resnum];
	}
	if ( last_res ) {
		outset.insert(current_selection);
	}

	//return block_selections;
}

void
identify_ss_blocks_vec(utility::vector1<core::select::residue_selector::ResidueSubset> & out_vec, core::select::residue_selector::ResidueSubset const & selection ){
	utility::vector1<core::select::residue_selector::ResidueSubset> block_selections;
	core::select::residue_selector::ResidueSubset current_selection( selection.size(), false );

	//Note - does not take into account new chains!
	bool last_res = false;
	for ( core::Size count_resnum = 1; count_resnum <= selection.size(); ++count_resnum ) {
		if ( selection[count_resnum] ) {
			if ( !last_res ) {
				current_selection = core::select::residue_selector::ResidueSubset( selection.size(), false );
			}
			current_selection[count_resnum] = true;
		} else {
			if ( last_res ) {
				out_vec.push_back(current_selection);
				//block_selections.push_back(current_selection);
			}
		}
		last_res = selection[count_resnum];
	}
	if ( last_res ) {
		out_vec.push_back(current_selection);
		//block_selections.push_back(current_selection);
	}
	//return block_selections;
}

void
calculate_blocks(std::map< core::Size, core::Size > & outmap, core::pose::Pose const & pose ){
	//Helix and Sheet residues
	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		//TR << i<< " " << pose.secstruct(i) << std::endl;
		if ( pose.secstruct(i) != 'L' && !pose.residue(i).is_virtual_residue() ) {
			selection[i] = true;
		}
	}
	utility::vector1<core::select::residue_selector::ResidueSubset> subsets;
	identify_ss_blocks_vec(subsets, selection);
	core::Size block_num = 0;
	//std::map< core::Size, core::Size > res_blocks;
	//TR < "Total Blocks " << subsets.size() << std::endl;
	for ( auto const & subset: subsets ) {
		block_num+=1;
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( subset[i] ) {
				outmap[i] = block_num;
				//res_blocks[i] = block_num;
			}
		}
	}
	//return res_blocks;
}

void
calculate_helices(std::map< core::Size, core::Size > & outmap, core::pose::Pose const & pose, core::Size min_length){

	std::pair<core::Size,core::Size> current_pair;
	char last_dssp = 'L';

	utility::vector1<std::pair<core::Size,core::Size>> helices;

	for ( core::Size current_residue = 1; current_residue <= pose.size(); ++current_residue ) {
		if ( pose.secstruct(current_residue) == 'H' ) {
			if ( last_dssp != 'H' ) {
				current_pair.first = current_residue;
			}
			current_pair.second = current_residue;
		} else if ( last_dssp == 'H' ) {
			helices.push_back(current_pair);
		}

		last_dssp = pose.secstruct(current_residue);
	}

	if ( last_dssp == 'H' ) {
		helices.push_back(current_pair);
	}
	utility::vector1<std::pair<core::Size,core::Size>> good_helices;
	for ( auto working_pair : helices ) {
		if ( min_length == 0 || (working_pair.second - working_pair.first) >= (min_length-1) ) {
			good_helices.push_back(working_pair);
		}
	}


	//std::map< core::Size, core::Size > blocks;
	core::Size block_number = 0;

	for ( auto working_pair : good_helices ) {
		block_number+=1;
		for ( core::Size current_residue = working_pair.first; current_residue <= working_pair.second; ++current_residue ) {
			outmap[current_residue] = block_number;
			//blocks[current_residue] = block_number;
		}
	}
	//TR << "helices done" << std::endl;

	//return blocks;

}

void
calculate_blocks_from_subset(std::map< core::Size, core::Size > & outmap, utility::vector1< bool > const & selection){


	//std::map< core::Size, core::Size > block_assignments;
	bool last_res = false;
	core::Size block_count = 1;
	for ( core::Size count_resnum = 1; count_resnum <= selection.size(); ++count_resnum ) {
		if ( selection[count_resnum] ) {
			outmap[count_resnum] = block_count;
			//block_assignments[count_resnum] = block_count;
		} else {
			outmap[count_resnum] = 0;
			//block_assignments[count_resnum] = 0;
			if ( last_res ) {
				++block_count;
			}
		}
		last_res = selection[count_resnum];
	}
	//if(!last_res){
	// block_count--;
	//}

	//return block_assignments;
}

bool
all_L_dssp(core::pose::Pose const & pose){
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.secstruct(i) != 'L' ) return false;
	}
	return true;
}

core::Real
calculate_motif_score_bt_residues(
	sewing::scoring::MotifScorer const & scorer,
	core::conformation::Residue const & N_res,
	core::conformation::Residue const & C_res,
	char ss1,
	char ss2
){
	utility::vector1< core::conformation::Atom > const & N_atoms = N_res.atoms();
	numeric::xyzTransform< core::Real > const stub1 = core::pose::motif::get_backbone_reference_frame(N_atoms[ 1 ].xyz(), N_atoms[ 2 ].xyz(), N_atoms[ 3 ].xyz() );
	char aa1 = N_res.name1();

	utility::vector1< core::conformation::Atom > const & C_atoms = C_res.atoms();
	numeric::xyzTransform< core::Real > const stub2 = core::pose::motif::get_backbone_reference_frame(C_atoms[ 1 ].xyz(), C_atoms[ 2 ].xyz(), C_atoms[ 3 ].xyz() );
	char aa2 = C_res.name1();
	return scorer.get_score(stub1, ss1, aa1, stub2, ss2, aa2);
}

core::Real
calculate_distance_score_bt_residues(
	core::conformation::Residue const & N_res,
	core::conformation::Residue const & C_res,
	core::Real min_score,
	core::Real min_distance,
	core::Real dist_mult
){

	return -1*(min_score + (dist_mult * std::abs(N_res.atom(2).xyz().distance(C_res.atom(2).xyz())-min_distance)));
}
void
calculate_bw_pose_compat_motifs(
	std::map< std::string, core::Real > & outmap,
	core::pose::Pose const & pose,
	std::set<core::select::residue_selector::ResidueSubset> const & block_selections,
	bool drop_best /*true */,
	bool normalize_by_residues /*false */,
	core::Real max_pair_score /*-2.0*/,
	bool use_motifs)
{
	using namespace core::select::residue_selector;

	protocols::sewing::scoring::MotifScorer scorer = protocols::sewing::scoring::MotifScorer();

	//std::map< std::string, core::Real > out;
	core::Real worst_score = -99999.0;
	utility::vector1<core::Real> blockwise_scores;
	core::Size block_counter = 0;

	for ( ResidueSubset const & selection_1 : block_selections ) {
		TR << "Block Counter " << block_counter << std::endl;
		blockwise_scores.clear();
		++block_counter;
		for ( ResidueSubset const & selection_2 : block_selections ) {
			if ( selection_1 == selection_2 ) {
				continue;
			}

			core::Size N_size = 0;
			core::Size C_size = 0;

			for ( core::Size N_resnum = 1; N_resnum <= pose.size(); ++N_resnum ) {
				if ( selection_1[N_resnum] && !pose.residue(N_resnum).is_virtual_residue() ) {
					++N_size;
				}
			}
			for ( core::Size C_resnum = 1; C_resnum <= pose.size(); ++C_resnum ) {
				if ( selection_2[C_resnum] && !pose.residue(C_resnum).is_virtual_residue() ) {
					++C_size;
				}
			}
			if ( N_size == 0 || C_size == 0 ) {
				continue;
			}

			core::Real score = 0.0;
			core::Real running_score = 0.0;
			// begin
			core::Real global_max=0.0;
			for ( core::Size N_resnum = 1; N_resnum <= pose.size(); ++N_resnum ) {
				if ( !pose.residue(N_resnum).is_virtual_residue() && selection_1[N_resnum] ) {
					char ss1 = pose.secstruct(N_resnum);
					core::conformation::Residue const & N_res = pose.residue(N_resnum);
					core::Real max_score = 0.0;
					for ( core::Size C_resnum = 1; C_resnum <= pose.size(); ++C_resnum ) {
						if ( !pose.residue(C_resnum).is_virtual_residue() && selection_2[C_resnum] ) {
							char ss2 = pose.secstruct(C_resnum);
							core::conformation::Residue const & C_res = pose.residue(C_resnum);

							if ( use_motifs ) {
								running_score = calculate_motif_score_bt_residues(scorer, N_res, C_res, ss1, ss2);
							} else {
								running_score = calculate_distance_score_bt_residues(N_res,C_res,-2.0, 4.0, 0.1);
							}

							if ( max_pair_score != 0 ) {
								if ( running_score > -1 * max_pair_score ) {
									running_score = max_pair_score;
								}
							}
							if ( running_score > max_score ) {
								max_score = running_score;
							}
							if ( running_score > global_max ) {
								global_max = running_score;
							}
							//score -= running_score;
						}
					}
					score -= max_score;
				}
			}
			if ( drop_best ) {
				score += global_max;
			}
			if ( normalize_by_residues ) {
				if ( N_size < C_size ) {
					score = score/N_size;
				} else {
					score = score/C_size;
				}
			}

			blockwise_scores.push_back(score);
		}

		// do per element processing here
		//TR << "Per element processing" << std::endl;
		if ( blockwise_scores.size() == 0 ) {
			//return out;
			return;
		}
		std::sort(blockwise_scores.begin(),blockwise_scores.end());



		if ( blockwise_scores.size() < 2 ) {
			if ( blockwise_scores[1] > worst_score ) {
				worst_score = blockwise_scores[1];
			}
		} else {
			if ( blockwise_scores[2] > worst_score ) {
				worst_score = blockwise_scores[2];
			}
		}
		//TR << "Accessing out" << std::endl;
		outmap[std::to_string(block_counter)] = worst_score;
		//
	}
	//return out;
}

void
calculate_bw_window_motifs(
	std::map< std::string, core::Real > & outmap,
	core::pose::Pose const & pose,
	std::set< core::select::residue_selector::ResidueSubset> const & block_selections,
	core::Size window_width /* 3 */,
	bool use_motifs
){

	using namespace core::select::residue_selector;

	protocols::sewing::scoring::MotifScorer scorer = protocols::sewing::scoring::MotifScorer();

	//std::map< std::string, core::Real > out;

	core::Size window_counter = 0;
	core::Real maximum = 0;
	std::map<core::Size,core::Real> window_scores;
	for ( ResidueSubset const & selection_1 : block_selections ) {
		window_scores.clear();

		core::Size N_size = 0;
		for ( core::Size N_resnum = 1; N_resnum <= pose.size(); ++N_resnum ) {
			if ( selection_1[N_resnum] ) {
				window_scores[N_resnum] = 0;
				++N_size;
			}
		}

		for ( ResidueSubset const & selection_2 : block_selections ) {
			if ( selection_1 == selection_2 ) {
				continue;
			}
			core::Size C_size = 0;
			for ( core::Size C_resnum = 1; C_resnum <= pose.size(); ++C_resnum ) {
				if ( selection_2[C_resnum] ) {
					++C_size;
				}
			}
			if ( N_size == 0 || C_size == 0 ) {
				continue;
			}

			core::Real running_score = 0.0;
			for ( core::Size N_resnum = 1; N_resnum <= pose.size(); ++N_resnum ) {
				if ( !pose.residue(N_resnum).is_virtual_residue() && selection_1[N_resnum] ) {
					char ss1 = pose.secstruct(N_resnum);
					core::conformation::Residue const & N_res = pose.residue(N_resnum);

					core::Real max_score = 0.0;
					for ( core::Size C_resnum = 1; C_resnum <= pose.size(); ++C_resnum ) {
						if ( !pose.residue(C_resnum).is_virtual_residue() && selection_2[C_resnum] ) {
							char ss2 = pose.secstruct(C_resnum);
							core::conformation::Residue const & C_res = pose.residue(C_resnum);
							if ( use_motifs ) {
								running_score = calculate_motif_score_bt_residues(scorer, N_res, C_res, ss1, ss2);
							} else {
								running_score = calculate_distance_score_bt_residues(N_res,C_res,-2.0, 4.0, 0.1);
							}
							if ( running_score > max_score ) {
								max_score = running_score ;
							}
						}
					}
					if ( max_score > window_scores[N_resnum] ) {
						window_scores[N_resnum] = max_score;
					}
				}
			}

		}
		std::list<core::Real> score_window;
		for ( auto window_score : window_scores ) {
			score_window.push_back(window_score.second);
			if ( score_window.size() > window_width ) {
				score_window.pop_front();
			}
			maximum = *std::max_element(score_window.begin(),score_window.end());
			if ( score_window.size() == window_width ) {
				++window_counter;
				outmap[std::to_string(window_counter)] = -1 * maximum;
			}
		}
		if ( score_window.size() < window_width ) {
			++window_counter;
			outmap[std::to_string(window_counter)] = -1 * maximum;
		}
	}
	//return out;
}

} //pose_sewing
} //protocols


