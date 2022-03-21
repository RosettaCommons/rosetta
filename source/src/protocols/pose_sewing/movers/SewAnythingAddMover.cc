// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/SewAnythingAddMover.cc
/// @brief perform SEWING-like chimerization-based addition of supersecondary structural segments into previously virtual residues
/// @author frankdt (frankdt@email.unc.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/pose_sewing/movers/SewAnythingAddMover.hh>
#include <protocols/pose_sewing/movers/SewAnythingAddMoverCreator.hh>
#include <protocols/pose_sewing/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/random/WeightedSampler.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/datacache/DataCache.hh>
#include <core/scoring/methods/util.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/conformation/Atom.hh>

static basic::Tracer TR( "protocols.pose_sewing.movers.SewAnythingAddMover" );

namespace protocols {
namespace pose_sewing {
namespace movers {
using namespace protocols::filters;
using namespace protocols::pose_sewing::data_storage;
using namespace core::simple_metrics;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SewAnythingAddMover::SewAnythingAddMover():
	protocols::moves::Mover( SewAnythingAddMover::mover_name() )
{
	pose_vector_ = data_storage::TerminalDSSPSortedPoseVectorOP(new data_storage::TerminalDSSPSortedPoseVector());

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SewAnythingAddMover::SewAnythingAddMover( SewAnythingAddMover const &  ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SewAnythingAddMover::~SewAnythingAddMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

std::pair< std::string, std::string >
SewAnythingAddMover::segment_vector_dssp(char mod_terminus, char mod_terminus_dssp, std::string const & permissible_segment_ends) const
{

	if ( mod_terminus == 'N' ) {
		return std::make_pair( permissible_segment_ends,std::string(1,mod_terminus_dssp));
	} else {
		return std::make_pair(std::string(1,mod_terminus_dssp),permissible_segment_ends);
	}

}

/// @brief Apply the mover
void
SewAnythingAddMover::apply( core::pose::Pose& pose){

	protocols::simple_moves::MutateResidueOP mutator = protocols::simple_moves::MutateResidueOP(new protocols::simple_moves::MutateResidue());
	bool needs_dssp = true;
	for ( core::Size resnum = 1; resnum<= pose.size(); ++resnum ) {
		if ( pose.secstruct(resnum)!='L' ) {
			needs_dssp = false;
			break;
		}
	}
	if ( needs_dssp && allow_DSSP_insertion_ ) {
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);
	}
	char mod_terminus = 'X';
	char mod_terminus_dssp = 'X';
	core::Size attempts = 0;
	core::Size to_add_N_basis_residue;
	core::Size to_add_C_basis_residue;
	core::Size local_basis_residue = 1; //check
	core::Size incoming_basis_residue = 1; //check

	core::Size first_residue = 1; //check
	core::Size last_residue = 1; //check

	if ( pose.pdb_info()==nullptr ) {
		utility_exit_with_message("The pose does not have a PDBInfo object. Run an AddVirtualResiduesMover on it prior to chimerization.");
	}
	std::pair<core::Size,core::Size> current_residue_pair = find_sewable_region(pose);
	first_residue = current_residue_pair.first;
	last_residue = current_residue_pair.second;

	if ( pose.chain(first_residue)!=pose.chain(last_residue) ) {
		utility_exit_with_message("Sewable region spans chains!");
	}

	//first determine moddable terminus

	//Figure out matches between permissible
	//
	mod_terminus = find_mod_terminus(pose,current_residue_pair);
	//now mod it
	//new
	core::Size terminal_basis_residue = 1; //check
	core::Size inner_basis_residue = 1; //check
	int terminus_direction = 0;

	// end new
	//bool has_been_sewn = false;
	if ( mod_terminus == 'N' ) {
		terminal_basis_residue=first_residue;
		terminus_direction = 1;
	} else {
		terminal_basis_residue=last_residue;
		terminus_direction = -1;
	}
	while ( pose.residue(terminal_basis_residue).is_virtual_residue() ) {
		terminal_basis_residue = terminal_basis_residue+terminus_direction;
	}
	if ( terminal_basis_residue < 1 || terminal_basis_residue > pose.size() ) {
		utility_exit_with_message("Bad terminal basis residue!");
	}
	if ( trim_terminal_loops_ ) {
		while ( pose.secstruct(terminal_basis_residue) == 'L' ) {
			if ( pose.secstruct(terminal_basis_residue+terminus_direction) != 'L' ) {
				pose.set_secstruct(terminal_basis_residue,pose.secstruct(terminal_basis_residue+terminus_direction));
			} else {
				terminal_basis_residue = terminal_basis_residue+terminus_direction;
			}
		}
	}
	mod_terminus_dssp = pose.secstruct(terminal_basis_residue);

	//If we have settings for the mod terminus DSSP, use that, otherwise use the defaults.
	utility::vector1<core::Size> vitals;
	for ( core::Size check_residue = first_residue; check_residue <= last_residue; ++check_residue ) {
		if ( pose.pdb_info()->res_haslabel(check_residue,"VITAL_RESIDUE") ) {
			vitals.push_back(check_residue);
		}
	}
	for ( auto current: vitals ) {
		if ( (mod_terminus == 'N' && current < terminal_basis_residue) || (mod_terminus == 'C' && current > terminal_basis_residue) ) {
			utility_exit_with_message("vital residues incompatible with selected terminus, retry");
			return;
		}
	}


	// constant section
	if ( permissible_segment_ends_.find(mod_terminus_dssp) == std::string::npos ) {
		TR<< "Permissible seg ends set to " << permissible_segment_ends_ <<" but found " <<mod_terminus_dssp << std::endl;
		TR << "selected terminus DSSP is outside of permissible set. Please either expand the permissible set of terminus DSSP codes, exclude this terminus, or truncate to a permissible DSSP." << std::endl;
		set_last_move_status(protocols::moves::MS_FAIL);
		return;
	}
	if ( pose.secstruct(terminal_basis_residue+(window_width_*terminus_direction))!=mod_terminus_dssp ) {
		TR << "secstruct of " << terminal_basis_residue+(window_width_*terminus_direction) << " is " << pose.secstruct(terminal_basis_residue+(window_width_*terminus_direction)) << std::endl;
		TR << "terminal secondary structure is too short. ABORTING CHIMERIZATION!" << std::endl;
		set_last_move_status(protocols::moves::MS_FAIL);
		return;
	}
	inner_basis_residue = terminal_basis_residue;
	while ( pose.secstruct(inner_basis_residue) == mod_terminus_dssp && inner_basis_residue >= first_residue && inner_basis_residue <= last_residue ) {
		inner_basis_residue = inner_basis_residue+terminus_direction;
	}

	core::Size abs_inner_basis_residue = inner_basis_residue;
	core::Size terminal_length = 0;
	terminal_length = (inner_basis_residue-terminal_basis_residue)*terminus_direction;
	// constant section
	if ( inner_basis_residue >= pose.size() ) {
		--inner_basis_residue;
	}
	if ( inner_basis_residue == 0 ) {
		++inner_basis_residue;
	}

	//
	if ( mod_terminus == 'N' ) { //was originally N
		for ( core::Size inner_basis_check_residue = inner_basis_residue; inner_basis_check_residue >= terminal_basis_residue; --inner_basis_check_residue ) {
			if ( inner_basis_check_residue < 1 || inner_basis_check_residue > pose.size() ) {
				TR << "INNER BASIS CHECK RESIDUE ERROR: " << inner_basis_check_residue << std::endl;
				utility_exit_with_message("Bad inner basis residue: "+ utility::to_string(inner_basis_check_residue));
			}
			if ( pose.pdb_info()->res_haslabel(inner_basis_check_residue,"VITAL_RESIDUE") ) {
				inner_basis_residue = inner_basis_check_residue;
			}
		}
	} else {
		for ( core::Size inner_basis_check_residue = inner_basis_residue; inner_basis_check_residue <= terminal_basis_residue; ++inner_basis_check_residue ) {
			if ( inner_basis_check_residue < 1 || inner_basis_check_residue > pose.size() ) {
				TR << "INNER BASIS CHECK RESIDUE ERROR: " << inner_basis_check_residue << std::endl;
				utility_exit_with_message("Bad inner basis residue: " + utility::to_string(inner_basis_check_residue));
			}
			if ( pose.pdb_info()->res_haslabel(inner_basis_check_residue,"VITAL_RESIDUE") ) {
				inner_basis_residue = inner_basis_check_residue;
			}
		}
	}
	inner_basis_residue = inner_basis_residue - terminus_direction;
	if ( inner_basis_residue*terminus_direction < terminal_basis_residue*terminus_direction ) {
		TR << "set of permissible residues is empty; check vital residues relative to window width. ABORTING CHIMERIZATION!" << std::endl;
		set_last_move_status(protocols::moves::MS_FAIL);
		return;
	}
	protocols::pose_sewing::data_storage::PoseWithTerminalSegmentsOfKnownDSSPOP segments_to_add;
	core::pose::PoseOP pose_to_add;
	core::Size innermost_contiguous_segment_residue;
	bool found_good_match = false;

	utility::vector1<data_storage::PoseWithTerminalSegmentsOfKnownDSSPOP> segment_vector;

	//Attempt even sampling, but only on segments that match terminus dssp - or we just fail for no reason.
	std::pair < std::string, std::string > allowed_seg_dssp = segment_vector_dssp(mod_terminus, mod_terminus_dssp, permissible_segment_ends_);

	//Deal with segment filters.
	core::Size n_segments_to_read = max_attempts_;

	core::Size seg_filter_attempts = 0;
	if ( pose_vector_ == nullptr || pose_vector_->get_total_size() == 0 ) {

		if ( segment_file_paths_.size() > 1 && even_sampling_ ) {
			utility::vector1< core::Size > acceptable_indexes;
			std::map< core::Size, utility::vector1<data_storage::PoseWithTerminalSegmentsOfKnownDSSPOP > > segment_vectors;

			//First, populate list of vectors that are actually compatable.
			for ( core::Size index = 1; index <= segment_file_paths_.size(); ++index ) {
				//Once we populate the segment vector, we don't use the pose_vector_.
				//Figure out some way to simplify this down the line.
				segment_vector = pose_vector_->populate_from_segment_file_and_get_random_vector_set(segment_file_paths_.at(index), allowed_seg_dssp.first, allowed_seg_dssp.second, n_segments_to_read);

				if ( segment_vector.size() ) {
					segment_vectors[index] = segment_vector;
					acceptable_indexes.push_back(index);
				}
			}
			//Now select from acceptable segment vectors.
			if ( acceptable_indexes.size() ) {
				core::Size index_lookup = numeric::random::rg().random_range(1, acceptable_indexes.size());
				core::Size index = acceptable_indexes[index_lookup];

				segment_vector = segment_vectors[index];
				segment_vectors.clear();
			}
		} else {
			segment_vector = pose_vector_->populate_from_segment_file_and_get_random_vector_set(segment_file_paths_, allowed_seg_dssp.first, allowed_seg_dssp.second, n_segments_to_read);
		}
	} else {
		segment_vector = pose_vector_->get_vector(allowed_seg_dssp.first, allowed_seg_dssp.second);
	}


	numeric::random::random_permutation(segment_vector);


	//If a post_filter or sort is set, we try all basis pairs.
	bool try_all_basis_pairs = false;
	if ( sort_metric_ != nullptr || all_post_filters_.size() != 0 ) try_all_basis_pairs = true;
	utility::vector1< std::pair< core::Size, core::Size > > non_triaged_basis_pairs;

	std::map < std::pair< core::Size, core::Size >, core::pose::Pose > aligned_pieces;

	for ( core::Size outer_try = 1; outer_try <= segment_vector.size(); ++outer_try ) {
		if ( found_good_match || attempts > max_attempts_ || seg_filter_attempts > max_filter_attempts_ ) {
			//outer_try = segment_vector.size()+1;
			break;
		}
		segments_to_add = segment_vector[outer_try];

		if ( (segments_to_add->get_N_term_length() > max_nter_len_) ||
				segments_to_add->get_C_term_length() > max_cter_len_ ) {
			continue; // Continue, but don't count it as an attempt
		}

		pose_to_add = segments_to_add->get_source_pose_op( false /*clone_if_new*/);

		if ( mod_terminus == 'N' ) {
			innermost_contiguous_segment_residue = pose_to_add->size() - segments_to_add->get_C_term_length();
			if ( pose_to_add->size()-innermost_contiguous_segment_residue < window_width_ ) {
				TR << "breaking loop. size-last is " << pose_to_add->size()-innermost_contiguous_segment_residue << " versus window width " << window_width_ << std::endl;
				continue;
			}
			to_add_N_basis_residue = innermost_contiguous_segment_residue;
			to_add_C_basis_residue = pose_to_add->size()-window_width_;
		} else {
			innermost_contiguous_segment_residue = segments_to_add->get_N_term_length();
			if ( innermost_contiguous_segment_residue < window_width_ ) {
				TR << "breaking loop. first is " << segments_to_add->get_N_term_length() << " versus window width " << window_width_ << std::endl;
				continue;
			}
			to_add_N_basis_residue = 1;
			to_add_C_basis_residue = innermost_contiguous_segment_residue - window_width_;
		}

		bool passed_filters = true;

		if ( pose_to_add && (all_seg_filters_.size() || seg_filters_.count(segments_to_add->get_segfile_path())) ) {

			for ( FilterCOP filter : all_seg_filters_ ) {
				if ( ! filter->apply(*pose_to_add) ) {
					passed_filters = false;
					break;
				}
			}
			if ( passed_filters && seg_filters_.count(segments_to_add->get_segfile_path()) ) {
				for ( FilterCOP filter : seg_filters_.at(segments_to_add->get_segfile_path()) ) {
					if ( ! filter->apply(*pose_to_add) ) {
						passed_filters = false;
						break;
					}
				}
			}
		}
		if ( ! passed_filters ) {
			seg_filter_attempts+=1;
			continue; //Go to next segment!
		}

		//get the actual basis residues to use
		//size check here
		// BEGIN MOD
		utility::vector1<std::pair<core::Size,core::Size>> basis_pairs;
		std::pair<core::Size,core::Size> basis_pair;
		for ( core::Size local_basis_residue = terminal_basis_residue; local_basis_residue*terminus_direction < inner_basis_residue*terminus_direction; local_basis_residue = local_basis_residue + terminus_direction ) {
			for ( core::Size incoming_basis_residue = to_add_N_basis_residue; incoming_basis_residue != to_add_C_basis_residue; ++incoming_basis_residue ) {
				bool good_basis_residues=false;
				if ( mod_terminus=='N' ) {
					if ( !use_absolute_sizes_ || (((abs_inner_basis_residue-local_basis_residue)+(incoming_basis_residue-to_add_N_basis_residue)) <= hashable_element_max_size_ && ((abs_inner_basis_residue-local_basis_residue)+(incoming_basis_residue-to_add_N_basis_residue)) >= hashable_element_min_size_) ) {
						if ( !use_relative_sizes_ || (((abs_inner_basis_residue-local_basis_residue)+(incoming_basis_residue-to_add_N_basis_residue)) <= (terminal_length+hashable_element_relative_max_size_) && ((abs_inner_basis_residue-local_basis_residue)+(incoming_basis_residue-to_add_N_basis_residue)) >= (terminal_length-hashable_element_relative_min_size_)) ) {
							good_basis_residues = true;
						}
					}
				} else {
					if ( !use_absolute_sizes_ || (((local_basis_residue-abs_inner_basis_residue)+(to_add_C_basis_residue-incoming_basis_residue)) <= hashable_element_max_size_ && ((local_basis_residue-abs_inner_basis_residue)+(to_add_C_basis_residue-incoming_basis_residue)) >= hashable_element_min_size_) ) {
						if ( !use_relative_sizes_ || (((local_basis_residue-abs_inner_basis_residue)+(to_add_C_basis_residue-incoming_basis_residue)) <= (terminal_length+hashable_element_relative_max_size_) && ((local_basis_residue-abs_inner_basis_residue)+(to_add_C_basis_residue-incoming_basis_residue)) >= (terminal_length-hashable_element_relative_min_size_)) ) {
							good_basis_residues = true;
						}
					}
				}

				if ( good_basis_residues ) {
					basis_pair.first = local_basis_residue;
					basis_pair.second = incoming_basis_residue;
					basis_pairs.push_back(basis_pair);
				}
			}
		}
		// basis_pairs is now complete
		numeric::random::random_permutation(basis_pairs);


		//core::Size tries = 0;
		non_triaged_basis_pairs.clear();
		aligned_pieces.clear();
		for ( core::Size current_try = 1; current_try <= basis_pairs.size(); ++current_try ) {
			//tries+=1;
			local_basis_residue = basis_pairs[current_try].first;
			incoming_basis_residue = basis_pairs[current_try].second;

			if ( this->align(pose,local_basis_residue,*pose_to_add,incoming_basis_residue,window_width_,alignment_max_distance_) ) {
				if ( mod_terminus == 'N' ) {

					if ( !(this->clashes(*pose_to_add,incoming_basis_residue,pose_to_add->size(),pose,first_residue,local_basis_residue+window_width_,clash_radius_,check_all_backbone_atoms_)) ) {

						found_good_match=true;
						non_triaged_basis_pairs.push_back(basis_pairs[current_try]);
						if ( ! try_all_basis_pairs ) break;
					}
				} else {
					if ( !(this->clashes(pose,local_basis_residue,last_residue,*pose_to_add,1,incoming_basis_residue+window_width_,clash_radius_,check_all_backbone_atoms_)) ) {

						found_good_match=true;
						non_triaged_basis_pairs.push_back(basis_pairs[current_try]);

						if ( ! try_all_basis_pairs ) break;
					}
				}

			} else { }
		}
		++attempts;
	}

	core::pose::Pose best_pose;
	core::Real best_pose_score = 0;
	bool filters_passed = false;
	if ( found_good_match ) {
		core::pose::Pose original_pose = pose;
		//core::Size basis_pair_num = 0;
		for ( auto const & local_basis_pair : non_triaged_basis_pairs ) {
			//basis_pair_num+=1;
			pose = original_pose;
			local_basis_residue = local_basis_pair.first;
			//core::Size const local_basis_residue_const = local_basis_pair.first;

			incoming_basis_residue = local_basis_pair.second;
			//core::Size const incoming_basis_residue_const = local_basis_pair.second;

			//Realign
			this->align(pose,local_basis_residue,*pose_to_add,incoming_basis_residue,window_width_,alignment_max_distance_,  false /* check_alignment */);
			if ( mod_terminus == 'N' ) {

				while ( local_basis_residue > (first_residue-1) && incoming_basis_residue > 0 && local_basis_residue > 0 && local_basis_residue <= pose.size() && incoming_basis_residue <= pose_to_add->size() ) {
					if ( pose.pdb_info()->res_haslabel(local_basis_residue,"VITAL_RESIDUE") ) {
						utility_exit_with_message("chimerizing over a vital residue!");
					}
					if ( pose.pdb_info()->res_haslabel(local_basis_residue,"IMMUTABLE_RESIDUE") ) {
						mutator->set_target(incoming_basis_residue);
						mutator->set_res_name(pose_to_add->residue(local_basis_residue).name3());
						mutator->apply(*pose_to_add);
						for ( core::Size current_chi = 1; current_chi <= pose_to_add->residue(incoming_basis_residue).nchi(); ++current_chi ) {
							pose_to_add->set_chi(current_chi,current_residue_pair.first,pose.chi(current_chi,local_basis_residue));
						}
					}
					pose.replace_residue( local_basis_residue, pose_to_add->residue(incoming_basis_residue), false);
					pose.set_secstruct(local_basis_residue, pose_to_add->secstruct(incoming_basis_residue));
					pose.pdb_info()->add_reslabel(local_basis_residue,segments_to_add->get_filename());
					if ( ! pose.pdb_info()->res_haslabel(local_basis_residue, "SEWN") ) {
						pose.pdb_info()->add_reslabel(local_basis_residue, "SEWN");
					}
					local_basis_residue--;
					incoming_basis_residue--;
				}
				for ( core::Size res_to_virtualize = local_basis_residue; res_to_virtualize >= first_residue; res_to_virtualize-- ) {
					pose.real_to_virtual( res_to_virtualize );
				}
				//and we're done!
			} else {
				while ( local_basis_residue < last_residue && incoming_basis_residue <= pose_to_add->size() && incoming_basis_residue > 0 && local_basis_residue > 0 && local_basis_residue <= pose.size() ) {
					if ( pose.pdb_info()->res_haslabel(local_basis_residue,"VITAL_RESIDUE") ) {
						utility_exit_with_message("chimerizing over a vital residue!");
					}
					if ( pose.pdb_info()->res_haslabel(local_basis_residue,"IMMUTABLE_RESIDUE") ) {
						mutator->set_target(local_basis_residue);
						mutator->set_res_name(pose.residue(local_basis_residue).name3());
						mutator->apply(*pose_to_add);
						for ( core::Size current_chi = 1; current_chi <= pose_to_add->residue(incoming_basis_residue).nchi(); ++current_chi ) {
							pose_to_add->set_chi(current_chi,current_residue_pair.first,pose.chi(current_chi,local_basis_residue));
						}
					}
					pose.replace_residue( local_basis_residue, pose_to_add->residue(incoming_basis_residue), false);
					pose.set_secstruct(local_basis_residue, pose_to_add->secstruct(incoming_basis_residue));
					pose.pdb_info()->add_reslabel(local_basis_residue,segments_to_add->get_filename());
					if ( ! pose.pdb_info()->res_haslabel(local_basis_residue, "SEWN") ) {
						pose.pdb_info()->add_reslabel(local_basis_residue, "SEWN");
					}
					++local_basis_residue;
					++incoming_basis_residue;
				}
				for ( core::Size res_to_virtualize = local_basis_residue; res_to_virtualize <= last_residue; ++res_to_virtualize ) {
					pose.real_to_virtual( res_to_virtualize );
				}
			}
			core::pose::set_reasonable_fold_tree(pose);
			//Here we sort or pass and break

			//1. Filter but no sort metric.  We stop at the first passing.
			//2. Sort metric.  We just add it to the list and return later.
			//3. Sort metric and filter.  We apply the filters and sort what we have left.
			bool filter_passed = false;

			//Score if present.
			core::Real score = 0;
			if ( sort_metric_ ) score = sort_metric_->calculate(pose);
			if ( positive_scores_are_better_ ) score = -1*score;


			if ( all_post_filters_.size() && sort_metric_ == nullptr ) {
				filter_passed = true;
				for ( FilterCOP filter : all_post_filters_ ) {
					if ( !filter->apply(pose) ) {
						filter_passed = false;
						break;
					}
				}
				if ( filter_passed ) {
					found_good_match = true;
					filters_passed = true;
					break; //Break out of the main loop.
				} else {
					found_good_match = false; //Set to false as this may be the last basis pair and we still didn't break out.
					filters_passed = false;

				}

			} else if ( all_post_filters_.size() && sort_metric_ != nullptr ) {
				bool local_filters_passed = true;
				for ( FilterCOP filter : all_post_filters_ ) {

					if ( !filter->apply(pose) ) {
						local_filters_passed = false;
						break;
					}
				}

				//If no local filters passed, continue onto the next chimerization to test.
				if ( local_filters_passed ) {
					filters_passed = true;
				} else {
					continue;
				}

				//If the filters passed, and this is the first, save it
				// Else, set the best pose.
				if ( best_pose_score == 0 ) {
					best_pose = pose;
					best_pose_score = score;
				} else if ( best_pose_score != 0 && score <= best_pose_score ) {

					best_pose_score = score;
					best_pose = pose;

				} else { }
			} else if ( sort_metric_ != nullptr ) {
				if ( best_pose_score != 0 ) {
					if ( score <= best_pose_score ) {

						best_pose = pose;
						best_pose_score = score;
					}

				} else {
					best_pose = pose;
					best_pose_score = score;
				}
			}
		} //End for each basis pair
		//Set to no good match if all filters have failed.
		if ( all_post_filters_.size() && filters_passed == false ) {
			found_good_match = false;
		}
	}

	if ( (! found_good_match) && fail_on_no_match_ ) {
		set_last_move_status(protocols::moves::MS_FAIL);
		utility_exit_with_message("NO MATCHES FOUND in SewAnything! ");
	}

	//Replace best pose
	if ( sort_metric_ != nullptr ) {
		pose = best_pose;
	}

}
std::pair<core::Size,core::Size>
SewAnythingAddMover::find_sewable_region(core::pose::Pose& pose){
	std::pair<core::Size,core::Size> current_residue_pair;
	utility::vector1<std::pair<core::Size,core::Size>> residue_pairs;
	for ( core::Size resnum = 1; resnum<= pose.size(); ++resnum ) {
		if ( pose.pdb_info()->res_haslabel(resnum,"PRE_POSE_SEWING_START") ) {
			current_residue_pair.first = resnum;
		}
		if ( pose.pdb_info()->res_haslabel(resnum,"PRE_POSE_SEWING_END") ) {
			current_residue_pair.second = resnum;
			residue_pairs.push_back(current_residue_pair);
		}
	}
	return residue_pairs.at(numeric::random::random_range(1,residue_pairs.size()));

}
char
SewAnythingAddMover::find_mod_terminus(core::pose::Pose& pose,std::pair<core::Size,core::Size> current_residue_pair){
	if ( permissible_termini_ == "N" ) {
		return 'N';
	} else if ( permissible_termini_ == "C" ) {
		return 'C';
	} else {
		if ( !(pose.residue(current_residue_pair.first).is_virtual_residue()) ) {
			return 'C';
		} else if ( !(pose.residue(current_residue_pair.second).is_virtual_residue()) ) {
			return 'N';
		} else if ( numeric::random::uniform()<0.5 ) {
			return 'C';
		} else {
			return 'N';
		}
	}
}

bool
SewAnythingAddMover::align(core::pose::Pose const & stationary_pose, core::Size stationary_residue, core::pose::Pose& mobile_pose, core::Size mobile_residue, core::Size window_width, core::Real tolerance,  bool check_alignment /*true*/) const {

	if ( stationary_residue < 1 || stationary_residue > stationary_pose.size() || mobile_residue < 1 || mobile_residue > mobile_pose.size() ) {
		TR << "Mobile or stationary basis pair residues out of range.  Skipping alignment." << std::endl;
		return false;
	}


	//Create HT b/t local_basis and incoming_basis
	core::conformation::Residue const & stationary_basis_residue = stationary_pose.residue(stationary_residue);
	core::conformation::Residue const & mobile_basis_residue = mobile_pose.residue(mobile_residue);

	utility::vector1< core::conformation::Atom > const & stationary_basis_atoms = stationary_basis_residue.atoms();
	utility::vector1< core::conformation::Atom > const & mobile_basis_atoms = mobile_basis_residue.atoms();

	numeric::HomogeneousTransform< core::Real > const stationary_ht( stationary_basis_atoms[ 3 ].xyz(), stationary_basis_atoms[ 1 ].xyz(), stationary_basis_atoms[ 2 ].xyz() );
	numeric::HomogeneousTransform< core::Real > const mobile_ht( mobile_basis_atoms[ 3 ].xyz(), mobile_basis_atoms[ 1 ].xyz(), mobile_basis_atoms[ 2 ].xyz() );

	numeric::HomogeneousTransform< core::Real > const inverse_mobile_ht = mobile_ht.inverse();
	numeric::HomogeneousTransform< core::Real > const mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;


	//Align mobile pose to stationary pose using a HT at the local basis residue of the stationary pose and the incoming basis residue of the mobile pose.
	mobile_pose.apply_transform_Rx_plus_v(mobile_to_stationary_ht.rotation_matrix(),mobile_to_stationary_ht.point());


	//check alignment at basis_residue+window_width for good alignment
	if ( check_alignment ) {
		numeric::xyzVector<core::Real> const & stationary_N = stationary_pose.residue(stationary_residue+window_width).xyz(1);
		numeric::xyzVector<core::Real> const & mobile_N = mobile_pose.residue(mobile_residue+window_width).xyz(1);
		numeric::xyzVector<core::Real> const & stationary_CA = stationary_pose.residue(stationary_residue+window_width).xyz(2);
		numeric::xyzVector<core::Real> const & mobile_CA = mobile_pose.residue(mobile_residue+window_width).xyz(2);
		numeric::xyzVector<core::Real> const & stationary_CB = stationary_pose.residue(stationary_residue+window_width).xyz(3);
		numeric::xyzVector<core::Real> const & mobile_CB = mobile_pose.residue(mobile_residue+window_width).xyz(3);
		if ( stationary_N.distance(mobile_N)>tolerance || stationary_CA.distance(mobile_CA)>tolerance || stationary_CB.distance(mobile_CB)>tolerance ) {
			TR.Debug << "Alignment failed! " << tolerance << " exceeded! N-N: " << stationary_N.distance(mobile_N) << " CA-CA: " << stationary_CA.distance(mobile_CA) << " CB-CB: " << stationary_CB.distance(mobile_CB) << std::endl;
			return false;
		}
	}
	return true;
}

bool
SewAnythingAddMover::clashes(core::pose::Pose const & nterm_pose, core::Size first_nterm_res, core::Size last_nterm_res, core::pose::Pose const & cterm_pose, core::Size first_cterm_res, core::Size last_cterm_res, core::Real clash_radius, bool check_all_backbone)const{


	for ( core::Size nterm_resnum = 1; nterm_resnum <= nterm_pose.size(); ++nterm_resnum ) {
		if ( (!nterm_pose.residue(nterm_resnum).is_virtual_residue()) && (nterm_resnum < first_nterm_res || nterm_resnum > last_nterm_res) ) {
			numeric::xyzVector<core::Real> const & nterm_CA_ = nterm_pose.residue(nterm_resnum).xyz(2);
			char const nter_ss = nterm_pose.secstruct(nterm_resnum);

			for ( core::Size cterm_resnum = 1; cterm_resnum <= cterm_pose.size(); ++cterm_resnum ) {
				if ( (!cterm_pose.residue(cterm_resnum).is_virtual_residue()) && (cterm_resnum < first_cterm_res || cterm_resnum > last_cterm_res) ) {

					numeric::xyzVector<core::Real> const & cterm_CA = cterm_pose.residue(cterm_resnum).xyz(2);
					char const cter_ss = cterm_pose.secstruct(cterm_resnum);
					if ( (nter_ss == 'E') && (cter_ss == 'E') && enable_clash_check_v2_ ) {}
					else if ( nterm_CA_.distance(cterm_CA)<clash_radius ) {
						return true;
					} else if ( check_all_backbone ) {
						numeric::xyzVector<core::Real> const & nterm_N = nterm_pose.residue(nterm_resnum).xyz(1);
						numeric::xyzVector<core::Real> const & nterm_C = nterm_pose.residue(nterm_resnum).xyz(3);
						numeric::xyzVector<core::Real> const & nterm_CB = nterm_pose.residue(nterm_resnum).xyz(4);
						numeric::xyzVector<core::Real> const & cterm_N = cterm_pose.residue(cterm_resnum).xyz(1);
						numeric::xyzVector<core::Real> const & cterm_C = cterm_pose.residue(cterm_resnum).xyz(3);
						numeric::xyzVector<core::Real> const & cterm_CB = cterm_pose.residue(cterm_resnum).xyz(4);
						if ( nterm_CA_.distance(cterm_N) < clash_radius || nterm_CA_.distance(cterm_C) < clash_radius || nterm_CA_.distance(cterm_CB) < clash_radius || nterm_C.distance(cterm_N) < clash_radius || nterm_C.distance(cterm_C) < clash_radius || nterm_C.distance(cterm_CB) < clash_radius || nterm_C.distance(cterm_CA) < clash_radius || nterm_N.distance(cterm_N) < clash_radius || nterm_N.distance(cterm_C) < clash_radius || nterm_N.distance(cterm_CB) < clash_radius ||  nterm_N.distance(cterm_CA) < clash_radius || nterm_CB.distance(cterm_N) < clash_radius || nterm_CB.distance(cterm_C) < clash_radius || nterm_CB.distance(cterm_CB) < clash_radius ||  nterm_CB.distance(cterm_CA) < clash_radius ) {
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

data_storage::TerminalDSSPSortedPoseVectorOP
SewAnythingAddMover::get_pose_vector() const{
	return pose_vector_;
}

core::Real
SewAnythingAddMover::get_clash_radius() const{
	return clash_radius_;
}

utility::vector1< std::string >
SewAnythingAddMover::get_segment_file_paths() const{
	return segment_file_paths_;
}

bool
SewAnythingAddMover::get_trim_terminal_loops() const{
	return trim_terminal_loops_;
}

std::string
SewAnythingAddMover::get_permissible_segment_ends() const{
	return permissible_segment_ends_;
}

void
SewAnythingAddMover::set_permissible_segment_ends( std::string in){
	permissible_segment_ends_ = in;
}
std::string
SewAnythingAddMover::get_permissible_termini() const{
	return permissible_termini_;
}
core::Size
SewAnythingAddMover::get_max_attempts() const{
	return max_attempts_;
}

void
SewAnythingAddMover::set_segment_file_paths(utility::vector1<std::string> segment_file_paths){

	segment_file_paths_ = segment_file_paths;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
SewAnythingAddMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
SewAnythingAddMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data_map)
{

	//Re-initialize
	pose_vector_ = data_storage::TerminalDSSPSortedPoseVectorOP(new data_storage::TerminalDSSPSortedPoseVector());

	seg_filters_.clear();
	all_seg_filters_.clear();
	all_post_filters_.clear();
	sort_metric_ = nullptr;


	//SimpleMetricCOP metric = get_metric_from_datamap_and_subtags(tag, data_map, "sort_metric", false /*fail_on_missing_from_factory*/);
	if ( tag->hasOption("sort_metric") ) {
		std::string gen = tag->getOption< std::string >( "sort_metric" );
		if ( data_map.has( "SimpleMetric", gen ) ) {
			SimpleMetricCOP metric( data_map.get_ptr< SimpleMetric >( "SimpleMetric", gen ) );
			if ( metric != nullptr ) {
				if ( metric->simple_metric_type() != "RealMetric" ) {
					utility_exit_with_message("Only real metrics can be used for sorting in SewAnythingAddMover!");
				}
				sort_metric_ = utility::pointer::dynamic_pointer_cast< RealMetric const >( metric );
			}
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "SimpleMetric " + gen + " not found in basic::datacache::DataMap.");
		}
	}

	window_width_ = tag->getOption< core::Size >( "window_width", 4 ); //This will only be used in hashless sewing
	alignment_max_distance_ = tag->getOption< core::Real >( "alignment_max_distance", 0.5 ); //This will only be used in hashless sewing

	if ( tag->hasOption("clash_radius") ) {
		clash_radius_ = tag->getOption< core::Real >( "clash_radius", 4.0 ); //This will only be used in hashless sewing
	}

	if ( tag->hasOption("segment_file_path") ) {
		std::string segment_file_path = tag->getOption<std::string>("segment_file_path");
		segment_file_paths_ = utility::string_split(segment_file_path, ',');
	}

	if ( tag->hasOption("clash_check_all_backbone_atoms") ) {
		check_all_backbone_atoms_ = tag->getOption<bool>("clash_check_all_backbone_atoms");
	}
	if ( tag->hasOption("trim_terminal_loops") ) {
		trim_terminal_loops_ = tag->getOption<bool>("trim_terminal_loops");
	}
	if ( tag->hasOption("allow_DSSP_insertion") ) {
		allow_DSSP_insertion_ = tag->getOption<bool>("allow_DSSP_insertion");
	}
	if ( tag->hasOption("max_attempts") ) {
		max_attempts_ = tag->getOption<core::Size>("max_attempts");
	}
	if ( tag->hasOption("seg_attempts") ) {
		max_filter_attempts_ = tag->getOption<core::Size>("seg_attempts");
	}
	if ( tag->hasOption("hashable_element_max_size") ) {
		hashable_element_max_size_ = tag->getOption<core::Size>("hashable_element_max_size");
	}
	if ( tag->hasOption("hashable_element_min_size") ) {
		hashable_element_min_size_ = tag->getOption<core::Size>("hashable_element_min_size");
	}
	if ( tag->hasOption("hashable_element_relative_max_size") ) {
		hashable_element_relative_max_size_ = tag->getOption<core::Size>("hashable_element_relative_max_size");
	}
	if ( tag->hasOption("hashable_element_relative_min_size") ) {
		hashable_element_relative_min_size_ = tag->getOption<core::Size>("hashable_element_relative_min_size");
	}
	if ( tag->hasOption("permissible_segment_ends") ) {
		permissible_segment_ends_ = tag->getOption<std::string>("permissible_segment_ends");
	}


	fail_on_no_match_ = tag->getOption("fail_on_no_match", fail_on_no_match_);

	if ( tag->hasOption("permissible_termini") ) {
		permissible_termini_ = tag->getOption<std::string>("permissible_termini");
	}

	use_absolute_sizes_ = tag->getOption<bool>("use_absolute_sizes", use_absolute_sizes_);
	use_relative_sizes_ = tag->getOption<bool>("use_relative_sizes", use_relative_sizes_);

	if ( use_absolute_sizes_ && use_relative_sizes_ ) {
		utility_exit_with_message("Cannot use relative and absolute sizes.  One must be false. ");
	}

	//Add filters pre-defined or given as subelement
	if ( tag->hasOption("seg_filters") ) {
		utility::vector1< std::string > filter_names = utility::string_split(tag->getOption< std::string >("seg_filters"), ',');
		for ( std::string const & fname : filter_names ) {
			all_seg_filters_.push_back(protocols::rosetta_scripts::parse_filter( fname, data_map ) );
			TR<<"Added Segment filter "<<fname<<std::endl;
		}
	}
	if ( tag->hasOption("post_filters") ) {
		utility::vector1< std::string > filter_names = utility::string_split(tag->getOption< std::string >("post_filters"), ',');
		for ( std::string const & fname : filter_names ) {
			all_post_filters_.push_back(protocols::rosetta_scripts::parse_filter( fname, data_map ) );
			TR<<"Added Post-Chimerization filter "<<fname<<std::endl;
		}
	}
	enable_clash_check_v2_ = tag->getOption("clash_check_v2", enable_clash_check_v2_);

	utility::vector1< std::string > ss_elements;
	ss_elements.push_back("E");
	ss_elements.push_back("H");
	ss_elements.push_back("L");

	even_sampling_ = tag->getOption< bool >("even_sampling", even_sampling_);
	positive_scores_are_better_ = tag->getOption< bool >("positive_metric_is_better", positive_scores_are_better_);

	//Now we set per-SS settings

	if ( segment_file_paths_.empty() ) {
		utility_exit_with_message("Segment file paths must be given in either the main tag or subblocks.");
	}

}

void SewAnythingAddMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	///Attributes that can also be set per-SS!
	XMLSchemaAttribute min_seg_size = XMLSchemaAttribute::attribute_w_default( "hashable_element_min_size", xsct_positive_integer, "how short can the final chimeric segments be","4" );

	XMLSchemaAttribute max_seg_size = XMLSchemaAttribute::attribute_w_default( "hashable_element_max_size", xsct_positive_integer, "how long can the final chimeric segments be","40" );

	XMLSchemaAttribute relative_min_seg_size = XMLSchemaAttribute::attribute_w_default( "hashable_element_relative_min_size", xsct_positive_integer, "how short can the final chimeric segments be relative to the initial segment","100" );

	XMLSchemaAttribute relative_max_seg_size = XMLSchemaAttribute::attribute_w_default( "hashable_element_relative_max_size", xsct_positive_integer, "how long can the final chimeric segments be relative to the initial segment","100" );

	XMLSchemaAttribute max_nter_len = XMLSchemaAttribute::attribute_w_default( "max_nter_len", xsct_positive_integer, "Maximum length of nter element","100" );

	XMLSchemaAttribute max_cter_len = XMLSchemaAttribute::attribute_w_default( "max_cter_len", xsct_positive_integer, "Maximum length of nter element","100" );

	XMLSchemaAttribute window_width = XMLSchemaAttribute::attribute_w_default( "window_width", xsct_positive_integer, "Required number of overlapping residues for two segments to be considered a match. Used in hashless SEWING only (for hashed SEWING, this is determined by the hasher settings used when generating the edge file).", "4" );

	XMLSchemaAttribute alignment_max = XMLSchemaAttribute::attribute_w_default( "alignment_max_distance", xsct_real, "maximum distance in Angstroms any atom of residues at the far side of a window can be", "0.5" );

	attlist

		+ XMLSchemaAttribute::attribute_w_default( "use_mt", xs_boolean, "Add MoltenTopology for segments and final pose?  Useful when working with sheets to reduce re-calculations","false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_absolute_sizes", xs_boolean, "consider the absolute sizes of the hashable elements","true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_relative_sizes", xs_boolean, "consider the relative sizes of the hashable elements","false" )
		+ XMLSchemaAttribute( "segment_file_path", xs_string, "Path to directory containing segment file.  Multiple paths can be given now, comma-delimited. Can also use the Segment sub block. Required if not using subblocks." )
		+ XMLSchemaAttribute::attribute_w_default( "even_sampling", xs_boolean, "Perform even sampling when using multiple segment files?  If false, segments with larger members will dominate.", "false")
		+ XMLSchemaAttribute::attribute_w_default( "permissible_segment_ends", xs_string, "DSSP codes of segment ends we are allowed to hybridize","EHL" )
		+ XMLSchemaAttribute::attribute_w_default( "permissible_termini", xs_string, "DSSP codes of segment ends we are allowed to hybridize","NC" )
		+ XMLSchemaAttribute::attribute_w_default( "max_attempts", xsct_positive_integer, "how many attempts can we make per add?","1000" )
		+ XMLSchemaAttribute::attribute_w_default( "seg_attempts", xsct_positive_integer, "How many segment attempts can make if using filters?","500" )
		+ XMLSchemaAttribute( "seg_filters", xs_string, "Optional Comma-delimited list of previously defined filter(s) that will be run on each chosen segment before chimerization.")
		+ XMLSchemaAttribute( "post_filters", xs_string, "Optional Comma-delimited list of previously defined filter(s) that will be run on after chimerization. Will try all basis pairs and stop when filter has passed.")
		+ XMLSchemaAttribute( "sort_metric", xs_string, "Optional RealMetric to be used for sorting.  Will sort on all basis pairs.  If combined with a filter, will only sort if the filter has passed. ")
		+ XMLSchemaAttribute::attribute_w_default( "positive_metric_is_better", xs_boolean, "If we are sorting, are positive scores better?","false" )
		+ min_seg_size
		+ max_seg_size
		+ relative_min_seg_size
		+ relative_max_seg_size
		+ max_cter_len
		+ max_nter_len
		+ XMLSchemaAttribute::attribute_w_default( "clash_radius", xsct_real, "how close two atoms can get before a clash is registered","4.0" )
		+ XMLSchemaAttribute::attribute_w_default( "trim_terminal_loops", xs_boolean, "do we make a segment_file?" , "true")
		+ window_width
		+ alignment_max
		+ XMLSchemaAttribute::attribute_w_default( "clash_check_all_backbone_atoms", xs_boolean, "should we check all backbone atoms of all residues during clash check", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "fail_on_no_match", xs_boolean, "Should we utility_exit if no good matches are found at the end of max_attempts?" , "false")
		+ XMLSchemaAttribute::attribute_w_default( "clash_check_v2", xs_boolean, "Enable Clash Check v2?  This reduces E-E clash radius to 3.5." , "true")
		+ XMLSchemaAttribute::attribute_w_default( "allow_DSSP_insertion", xs_boolean, "let sewing add DSSP if it needs" , "true");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "add_only one-step SEWING. Subelements of E,L or H can be given. Filter subelements can be added to filter segments.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SewAnythingAddMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SewAnythingAddMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SewAnythingAddMover::clone() const
{
	return protocols::moves::MoverOP( new SewAnythingAddMover( *this ) );
}

std::string SewAnythingAddMover::get_name() const {
	return mover_name();
}

std::string SewAnythingAddMover::mover_name() {
	return "SewAnythingAddMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
SewAnythingAddMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SewAnythingAddMover );
}

std::string
SewAnythingAddMoverCreator::keyname() const
{
	return SewAnythingAddMover::mover_name();
}

void SewAnythingAddMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SewAnythingAddMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, SewAnythingAddMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //pose_sewing
} //movers
