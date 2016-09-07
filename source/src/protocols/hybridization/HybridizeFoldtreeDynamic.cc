// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief A mover to dynamically manipulate fold tree during the template hybridization sampling
/// @details based on Chris Miles's start tree builder
/// @file protocols/hybridization/HybridizeFoldtreeDynamic.cc
/// @author Yifan Song, David Kim

// Unit headers
#include <protocols/hybridization/HybridizeFoldtreeDynamic.hh>

// C/C++ headers
#include <string>
#include <utility>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/DistributionSampler.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/chemical/VariantType.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/util/kinematics_util.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.hybridization.HybridizeFoldtreeDynamic" );

namespace protocols {
namespace hybridization {

HybridizeFoldtreeDynamic::HybridizeFoldtreeDynamic() : virtual_res_(-1) {}

utility::vector1 < core::Size > HybridizeFoldtreeDynamic::decide_cuts(core::pose::Pose & pose, core::Size n_residues) {
	// complete the chunks to cover the whole protein and customize cutpoints
	// cutpoints in the middle of the loop
	// cut number is the residue before the cut
	std::string cut_point_decision = "middle_L";

	using utility::vector1;
	vector1<core::Size> cuts_orig;
	vector1<std::pair<core::Size, core::Size> > jumps_orig;
	jumps_and_cuts_from_foldtree(initial_asymm_foldtree_, jumps_orig, cuts_orig);

	utility::vector1<bool> cut_options(n_residues, true);

	// not allowing cuts in the core_chunks, currectly is set by the secondary structure of the initial template
	for ( core::Size ichunk = 1; ichunk<=core_chunks_.num_loop(); ++ichunk ) {
		for ( core::Size ires = core_chunks_[ichunk].start(); ires <= core_chunks_[ichunk].stop()-1; ++ires ) {
			cut_options[ires] = false;
		}
	}

	utility::vector1 < core::Size > cut_positions;
	for ( core::Size i=2; i<=core_chunks_.num_loop(); ++i ) {
		core::Size loop_start = core_chunks_[i-1].stop();
		core::Size loop_end = core_chunks_[i].start() - 1;
		core::Size cut;

		if ( loop_start >= loop_end ) {
			cut = loop_start;
		} else {
			// if there is a cut in the original foldtree (in multi-chain cases), use that cut
			bool cut_found = false;
			for ( core::Size i_cut_old = 1; i_cut_old <= cuts_orig.size(); ++i_cut_old ) {
				if ( cuts_orig[i_cut_old] >= loop_start && cuts_orig[i_cut_old] <= loop_end ) {
					cut = cuts_orig[i_cut_old];
					TR.Debug << "Respect the input cut position: " << cut << std::endl;
					cut_found = true;
					break;
				}
			}
			if ( ! cut_found ) {
				if ( cut_point_decision == "middle" ) {
					cut = (loop_start + loop_end ) /2;
				} else if ( cut_point_decision == "middle_L" ) {
					cut = (loop_start + loop_end ) /2;
					if ( pose.secstruct(cut) != 'L' ) {
						utility::vector1< core::Size > candidates;
						for ( core::Size ic = cut-1; ic>=loop_start; --ic ) {
							if ( pose.secstruct(ic) == 'L' ) {
								candidates.push_back(ic);
								break;
							}
						}
						for ( core::Size ic = cut+1; ic<=loop_end; ++ic ) {
							if ( pose.secstruct(ic) == 'L' ) {
								candidates.push_back(ic);
								break;
							}
						}
						if ( candidates.size() ) {
							cut = candidates[numeric::random::rg().random_range(1,candidates.size())];
						}
					}
				} else if ( cut_point_decision == "combine" ) {
					utility::vector1 < core::Size > cut_residues;
					for ( core::Size ires = loop_start; ires <= loop_end; ++ires ) {
						if ( cut_options[ires] ) {
							cut_residues.push_back(ires);
						}
					}
					if ( cut_residues.size() > 0 ) {
						using boost::math::normal;
						using core::Size;
						using numeric::random::DistributionSampler;

						double mu = cut_residues.size() / 2.0;
						double sigma = cut_residues.size() / 5.0;
						normal distribution(mu, sigma);
						DistributionSampler<normal> sampler(distribution);

						// Clamp insertion position to closed interval [start, stop]
						Size position = static_cast<Size>(sampler.sample());
						position = numeric::clamp<core::Size>(position, 1, cut_residues.size());

						cut = cut_residues[position];
					} else {
						cut = (loop_start + loop_end ) /2;
					}
				} else {
					utility_exit_with_message("do not know how to make cut points");
				}
			}
		}
		cut_positions.push_back(cut);
		TR << "Cutpoint " << i << " at " << cut << std::endl;
		TR << "  Loop SS: ";
		std::string jnk;
		for ( core::Size ic = loop_start; ic<=loop_end; ++ic ) {
			TR << pose.secstruct(ic);
			jnk += (ic == cut) ? "*" : "-";
		}
		TR << std::endl;
		TR << "           " << jnk << std::endl;
	}
	return cut_positions;
}

void HybridizeFoldtreeDynamic::make_complete_chunks(
	utility::vector1 < core::Size > cut_positions,
	core::Size n_residues
)
{
	assert(cut_positions.size() == core_chunks_.size() - 1);

	chunks_ = core_chunks_;
	for ( core::Size i=1; i<=chunks_.num_loop(); ++i ) {
		if ( i==1 ) {
			chunks_[i].set_start(1);
		} else {
			core::Size new_start = cut_positions[i-1] + 1;
			chunks_[i].set_start( new_start );
		}

		if ( i==chunks_.num_loop() ) {
			chunks_[i].set_stop(n_residues);
		} else {
			core::Size new_stop = cut_positions[i];
			chunks_[i].set_stop( new_stop );
		}
	}
}

void HybridizeFoldtreeDynamic::choose_anchors() {
	anchor_positions_.clear();
	for ( core::Size i=1; i<=core_chunks_.num_loop(); ++i ) {
		anchor_positions_.push_back( choose_anchor_position(core_chunks_[i]) );
	}
}

// from cmiles:
//   mu : midpoint of the chunk
//   sigma : linear function of chunk length
core::Size HybridizeFoldtreeDynamic::choose_anchor_position(const protocols::loops::Loop & chunk) const {
	using boost::math::normal;
	using core::Size;
	using numeric::random::DistributionSampler;

	double mu = chunk.start() + chunk.length() / 2.0;
	double sigma = chunk.length() / 5.0;
	normal distribution(mu, sigma);
	DistributionSampler<normal> sampler(distribution);

	// Clamp insertion position to closed interval [start, stop]
	Size position = static_cast<Size>(sampler.sample());
	return numeric::clamp<core::Size>(position, chunk.start(), chunk.stop());
}


void HybridizeFoldtreeDynamic::initialize(
	core::pose::Pose & pose,
	protocols::loops::Loops const & core_chunks, // template chunks (may include single residue strand pairings as chunks)
	utility::vector1< std::pair< core::Size, core::Size > > const & strand_pairs,  // strand pair residue positions
	std::set<core::Size> const & strand_pair_library_positions // strand pair positions that are from the pairing library (not from a pdb template)
) {
	using core::Size;
	using core::Real;
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	using utility::vector1;

	num_nonvirt_residues_ = pose.total_residue();
	num_protein_residues_ = pose.total_residue();
	saved_n_residue_ = pose.total_residue();
	saved_ft_ = pose.conformation().fold_tree();

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		initial_asymm_foldtree_ = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
	} else {
		initial_asymm_foldtree_ = pose.conformation().fold_tree();
	}

	// strand pairings
	// initialize pairings data
	strand_pairs_ = strand_pairs;
	strand_pair_library_positions_ = strand_pair_library_positions;
	// identify core_chunks that are not from the strand pair library (i.e. these should be from a template)
	std::set<core::Size>::iterator set_iter;
	for ( core::Size i=1; i<=core_chunks.num_loop(); ++i ) {
		bool pair_chunk = false;
		for ( set_iter = strand_pair_library_positions_.begin(); set_iter != strand_pair_library_positions_.end(); ++set_iter ) {
			if ( core_chunks[i].start() <= *set_iter && core_chunks[i].stop() >= *set_iter ) {
				pair_chunk = true;
				break;
			}
		}
		if ( !pair_chunk ) {
			TR << "Identified template core chunk with index: " << i << std::endl;
			TR << core_chunks[i] << std::endl;
			template_core_chunk_indices_.insert(i);
		}
	} // strand pairings

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info=nullptr;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		num_nonvirt_residues_ = symm_info->num_independent_residues();
		num_protein_residues_ = symm_info->num_independent_residues();
	}
	while ( pose.residue_type(num_nonvirt_residues_).aa() == core::chemical::aa_vrt ) num_nonvirt_residues_--;
	while ( !pose.residue(num_protein_residues_).is_protein() ) num_protein_residues_--;

	// special case for ab initio
	if ( core_chunks.size() == 0 ) {
		protocols::loops::Loops pseudochunk;
		pseudochunk.push_back( protocols::loops::Loop(1,1) );
		set_core_chunks(pseudochunk);
	} else {
		set_core_chunks(core_chunks);
	}


	// set new
	choose_anchors();
	utility::vector1 < core::Size > cut_positions = decide_cuts(pose, num_protein_residues_);
	make_complete_chunks(cut_positions, num_protein_residues_);

	assert(chunks_.num_loop());

	// Add a virtual residue at the center of mass (updates the fold tree)
	if ( !symm_info ) {  // pose is asymm
		core::pose::addVirtualResAsRoot(pose);
	}

	virtual_res_ = pose.total_residue();

	// do the actual foldtree updates
	update(pose);

	protocols::loops::add_cutpoint_variants( pose );
}


// strand pairings
void HybridizeFoldtreeDynamic::get_pair_core_chunks(
	core::Size const chunk_index,
	utility::vector1<core::Size> & pair_chunks,
	utility::vector1<std::pair<core::Size, core::Size> > & pair_chunks_pairs
) {
	utility::vector1<core::Size> pair_chunks_tmp;
	utility::vector1<std::pair<core::Size, core::Size> > pair_chunks_pairs_tmp;
	for ( core::Size i=1; i<=strand_pairs_.size(); ++i ) {
		core::Size idx = 0;
		if ( core_chunks_[chunk_index].start() <= strand_pairs_[i].first && core_chunks_[chunk_index].stop() >= strand_pairs_[i].first ) {
			get_core_chunk_index_from_position( strand_pairs_[i].second, idx );
			if ( idx ) pair_chunks_pairs_tmp.push_back(std::pair<core::Size, core::Size>( strand_pairs_[i].first, strand_pairs_[i].second ));
		} else if ( core_chunks_[chunk_index].start() <= strand_pairs_[i].second && core_chunks_[chunk_index].stop() >= strand_pairs_[i].second ) {
			get_core_chunk_index_from_position( strand_pairs_[i].first, idx );
			if ( idx ) pair_chunks_pairs_tmp.push_back(std::pair<core::Size, core::Size>( strand_pairs_[i].second, strand_pairs_[i].first ));
		}
		if ( idx ) pair_chunks_tmp.push_back(idx);
	}
	assert( pair_chunks_tmp.size() == pair_chunks_pairs_tmp.size() );
	pair_chunks = pair_chunks_tmp;
	pair_chunks_pairs = pair_chunks_pairs_tmp;
} // strand pairings

void HybridizeFoldtreeDynamic::get_core_chunk_index_from_position( core::Size const position, core::Size & index ) {
	for ( core::Size i = 1; i<=core_chunks_.num_loop(); ++i ) {
		if ( core_chunks_[i].start() <= position && core_chunks_[i].stop() >= position ) {
			index = i;
			break;
		}
	}
}


void HybridizeFoldtreeDynamic::reset( core::pose::Pose & pose ) {
	if ( pose.total_residue() > saved_n_residue_ ) {
		pose.conformation().delete_residue_range_slow(saved_n_residue_+1, pose.total_residue());
	}
	pose.conformation().fold_tree( saved_ft_ );

	//protocols::loops::remove_cutpoint_variants( pose );
	core::pose::Pose init_pose = pose;
	bool pose_changed = false;
	using namespace core::chemical;
	for ( core::Size ir=1; ir< pose.total_residue() ; ++ir ) {
		if ( pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {

			bool is_cut = false;
			for ( int ic=1; ic<=pose.fold_tree().num_cutpoint() ; ++ic ) {
				core::Size cutpoint = pose.fold_tree().cutpoint(ic);
				if ( ir == cutpoint ) {
					is_cut = true;
					break;
				}
			}
			if ( !is_cut ) {
				pose_changed = true;
				core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, ir );

				if ( pose.residue(ir+1).has_variant_type(CUTPOINT_UPPER) ) {
					core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, ir+1 );
				}
			}
		}
	}

	if ( pose_changed == true ) pose.constraint_set( init_pose.constraint_set()->remapped_clone( init_pose, pose ) );

}

// stolen from protocols::forge::methods::jumps_and_cuts_from_pose
void HybridizeFoldtreeDynamic::jumps_and_cuts_from_foldtree( core::kinematics::FoldTree & ft, utility::vector1< std::pair< core::Size, core::Size > > & jumps, utility::vector1< core::Size > & cuts) {

	for ( core::Size i = 1; i<= ft.num_jump(); ++i ) {
		core::Size down ( ft.downstream_jump_residue(i) );
		core::Size up ( ft.upstream_jump_residue(i) );
		jumps.push_back( std::pair<core::Size,core::Size>( down, up ) );
	}
	cuts =  ft.cutpoints();
}

void HybridizeFoldtreeDynamic::update(core::pose::Pose & pose) {
	using core::Size;
	using core::Real;
	using protocols::loops::Loop;
	using protocols::loops::Loops;
	using utility::vector1;

	assert(chunks_.num_loop());

	bool use_symm = core::pose::symmetry::is_symmetric(pose);

	// Define jumps, cuts
	vector1<core::Size> cuts;
	vector1<std::pair<core::Size, core::Size> > jumps;

	// "symmetry-safe" version
	core::kinematics::FoldTree tree = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
	Size jump_root = num_nonvirt_residues_+1;
	if ( use_symm ) jump_root = anchor_positions_[1];

	// keep a copy of cuts and jumps, if they are in the region outside of the chunk definition
	vector1<Size> cuts_old;
	vector1<std::pair<Size, Size> > jumps_old;
	core::Size last_chunk_residue(chunks_[chunks_.num_loop()].stop());

	//if (!use_symm) {
	//if ( ! initial_asymm_foldtree_.is_simple_tree() ) {

	jumps_and_cuts_from_foldtree(initial_asymm_foldtree_, jumps_old, cuts_old);

	for ( Size i = 1; i <= jumps_old.size(); ++i ) {
		//if (jumps_old[i].first == jump_root) continue; // no need to skip this any more
		//if (jumps_old[i].second == jump_root) continue;

		if ( jumps_old[i].first > last_chunk_residue && jumps_old[i].first <= num_nonvirt_residues_ ) {
			jumps.push_back(std::make_pair(jump_root, jumps_old[i].first));
			TR.Debug << "Adding additional jump: " << jump_root << jumps_old[i].first << std::endl;
		} else if ( jumps_old[i].second > last_chunk_residue && jumps_old[i].second <= num_nonvirt_residues_ ) {
			jumps.push_back(std::make_pair(jump_root, jumps_old[i].second));
			TR.Debug << "Adding additional jump: " << jump_root << jumps_old[i].second << std::endl;
		}
	}
	for ( Size i = 1; i <= cuts_old.size(); ++i ) {
		if ( cuts_old[i] > last_chunk_residue && cuts_old[i] <= num_nonvirt_residues_ ) {
			cuts.push_back(cuts_old[i]);
			TR.Debug << "Adding additional cut: " << cuts_old[i] << std::endl;
		}
	}
	if ( !use_symm ) {
		if ( initial_asymm_foldtree_.nres() > last_chunk_residue && initial_asymm_foldtree_.nres() <= num_nonvirt_residues_ ) {
			cuts.push_back(initial_asymm_foldtree_.nres());
			TR.Debug << "Adding a cut on the last residue: " << initial_asymm_foldtree_.nres() << std::endl;
		}
	}
	//}

	// add root jumps and cuts for chunks that should be rooted to the star fold tree
	//   these should be chunks from a template or if none exist, the first strand pair chunk
	std::set<core::Size> root_chunk_indices = template_core_chunk_indices_;
	if ( root_chunk_indices.empty() ) { //size()) {
		TR << "Template core chunks do not exist so using first chunk" << std::endl;
		root_chunk_indices.insert(1);
	}
	std::set<core::Size>::iterator set_iter;
	for ( Size i = 1; i <= chunks_.num_loop(); ++i ) {
		const Loop& chunk = chunks_[i];
		const Size cut_point  = chunk.stop();
		const Size jump_point = anchor_positions_[i];

		Size j_root = num_nonvirt_residues_+1;
		if ( use_symm && i>1 ) j_root = anchor_positions_[1];

		TR.Debug << "Adding chunk cut: " << cut_point << std::endl;
		cuts.push_back(cut_point);

		if ( root_chunk_indices.count(i) ) {
			rooted_chunk_indices_.insert(i);
			TR.Debug << "Adding root jump: " << j_root << " " << jump_point << std::endl;
			jumps.push_back(std::make_pair(j_root, jump_point));
		}
	}

	// strand pairings

	if ( !strand_pairs_.empty() ) { //size()) {
		// add jumps and cuts to strand pair chunks that are from a rooted branch
		std::set<core::Size> rooted_chunk_indices; // keep track of all rooted branch chunks
		for ( set_iter = rooted_chunk_indices_.begin(); set_iter != rooted_chunk_indices_.end(); ++set_iter ) {
			add_overlapping_pair_chunks( *set_iter, cuts, jumps, rooted_chunk_indices );
		}
		// at this point all pair chunks that overlap with a rooted chunk or overlap with the branch should have cuts and jumps. These overlapping pair chunks
		// are in the global coordinate frame since they are branched from template core chunks (or the first strand pair chunk if templates do not exist)
		// i.e.  like |*-|-|-| where |* is a star tree rooted chunk, - is a jump and | are pair chunks

		// now lets add cuts to chunks that are not rooted yet (i.e. left over pairs, these are floating pairs with no reference to the global coordinate frame)
		// i.e.      |*
		//           |-| <- this is a pair that is not rooted yet and adjacent in sequence to the rooted chunk |*, so cuts have to be made to these pair chunks
		//                  and if a cut already exists between the rooted chunk and the pair chunk, the cut has to be removed
		std::set<core::Size> floating_pair_chunks;
		for ( core::Size i=1; i<=strand_pairs_.size(); ++i ) {
			core::Size i_index = 0;
			core::Size j_index = 0;
			get_core_chunk_index_from_position( strand_pairs_[i].first, i_index );
			get_core_chunk_index_from_position( strand_pairs_[i].second, j_index );
			assert( i_index && j_index );
			core::Size paircnt = 0;
			if ( !rooted_chunk_indices.count(i_index) && !rooted_chunk_indices_.count(i_index) ) {
				floating_pair_chunks.insert(i_index);
				paircnt++;
			}
			if ( !rooted_chunk_indices.count(j_index) && !rooted_chunk_indices_.count(j_index) ) {
				floating_pair_chunks.insert(j_index);
				paircnt++;
			}
			assert( paircnt == 0 || paircnt == 2 );
		}

		// first check if adjacent chunks are from the original rooted chunks and if they are add them to the rooted branch
		// i.e. check if adjacent chunks are from a template and then add the floating chunk and overlapping floating chunks
		for ( set_iter = rooted_chunk_indices_.begin(); set_iter != rooted_chunk_indices_.end(); ++set_iter ) {
			// check upstream
			core::Size up_index = *set_iter+1;
			if ( floating_pair_chunks.count(up_index) && !rooted_chunk_indices.count(up_index) ) {
				// remove cut of rooted chunk
				TR << "Removing root cut: " << chunks_[*set_iter].stop() << std::endl;
				remove_cut( chunks_[*set_iter].stop(), cuts );
				TR.Debug << "Adding floating pair chunk cut: " << chunks_[up_index].stop() << std::endl;
				cuts.push_back(chunks_[up_index].stop());
				rooted_chunk_indices.insert(up_index);
				add_overlapping_pair_chunks( up_index, cuts, jumps, rooted_chunk_indices );
			}
			// check downstream
			if ( *set_iter > 1 ) {
				core::Size down_index = *set_iter-1;
				if ( floating_pair_chunks.count(down_index) && !rooted_chunk_indices.count(down_index) ) {
					rooted_chunk_indices.insert(down_index);
					add_overlapping_pair_chunks( down_index, cuts, jumps, rooted_chunk_indices );
				}
			}
		}

		// remove rooted floating pairs from set
		for ( set_iter = rooted_chunk_indices.begin(); set_iter != rooted_chunk_indices.end(); ++set_iter ) {
			floating_pair_chunks.erase(*set_iter);
		}

		// second check if adjacent chunks are from a built up rooted branch chunk and if they are add them to the rooted branch
		// i.e. check if adjacent chunks are from a rooted branch and then add the floating chunk and overlapping floating chunks
		while ( !floating_pair_chunks.empty() ) { //size()) {
			for ( set_iter = rooted_chunk_indices.begin(); set_iter != rooted_chunk_indices.end(); ++set_iter ) {
				// check upstream
				core::Size up_index = *set_iter+1;
				if ( floating_pair_chunks.count(up_index) && !rooted_chunk_indices.count(up_index) ) {
					// remove cut of rooted chunk
					TR << "Removing root cut: " << chunks_[*set_iter].stop() << std::endl;
					remove_cut( chunks_[*set_iter].stop(), cuts );
					TR.Debug << "Adding floating pair chunk cut: " << chunks_[up_index].stop() << std::endl;
					cuts.push_back(chunks_[up_index].stop());
					rooted_chunk_indices.insert(up_index);
					add_overlapping_pair_chunks( up_index, cuts, jumps, rooted_chunk_indices );
				}
				// check downstream
				if ( *set_iter > 1 ) {
					core::Size down_index = *set_iter-1;
					if ( floating_pair_chunks.count(down_index) && !rooted_chunk_indices.count(down_index) ) {
						rooted_chunk_indices.insert(down_index);
						add_overlapping_pair_chunks( down_index, cuts, jumps, rooted_chunk_indices );
					}
				}
			}
			// remove rooted floating pairs from set
			for ( set_iter = rooted_chunk_indices.begin(); set_iter != rooted_chunk_indices.end(); ++set_iter ) {
				floating_pair_chunks.erase(*set_iter);
			}
		}

		// hopefully all chunks are covered now
		assert( chunks_.size() == rooted_chunk_indices.size() + rooted_chunk_indices_.size() );
	} // strand pairings


	// Remember to include the original cutpoint at the end of the chain
	// (before the virtual residue)
	// TR << "Adding the last cut: " << num_nonvirt_residues_ << std::endl;
	// cuts.push_back(num_nonvirt_residues_);

	TR.Debug << "jump size: " << jumps.size() << " cut size: " << cuts.size() << std::endl;

	ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
	for ( Size i = 1; i <= jumps.size(); ++i ) {
		TR.Debug << "jump " << i << " " << jumps[i].first << " " << jumps[i].second << std::endl;
		ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
		ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
	}

	ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
	for ( Size i = 1; i <= cuts.size(); ++i ) {
		TR.Debug << "cut " << i << " " << cuts[i] << std::endl;
		ft_cuts(i) = cuts[i];
	}

	bool status = tree.tree_from_jumps_and_cuts(num_nonvirt_residues_+1,   // nres_in
		jumps.size(),   // num_jump_in
		ft_jumps,    // jump_point
		ft_cuts,  // cuts
		num_nonvirt_residues_+1);  // root


	// strand pairings
	if ( strand_pairs_.size() ) {
		// must join adjacent peptide edges (cannot use FoldTree::delete_extra_vertices() because it skips jumps)
		// this might not be necessary because of the code following this
		TR << "tree_from_jumps_and_cuts: " << tree << std::endl;
		core::kinematics::FoldTree tmp_tree;
		bool join_edges = true;
		bool updated = false;
		while ( join_edges ) {
			join_edges = false;
			core::kinematics::FoldTree const const_tree( tree );
			tmp_tree.clear();
			for ( auto it = const_tree.begin(), it_end = const_tree.end(); it != it_end; ++it ) {
				auto it_next = it+1;
				if ( it_next != it_end ) {
					if ( it->label() == core::kinematics::Edge::PEPTIDE && it_next->label() == core::kinematics::Edge::PEPTIDE &&
							it->stop() == it_next->start() ) {
						join_edges = true;
						tmp_tree.add_edge( it->start(), it_next->stop(), core::kinematics::Edge::PEPTIDE );
						++it;
						updated = true;
						continue;
					}
					tmp_tree.add_edge( it->start(), it->stop(), it->label() );
				} else {
					tmp_tree.add_edge( it->start(), it->stop(), it->label() );
				}
			}
			tree = tmp_tree;
		}
		if ( updated ) TR << "joined adjacent peptide edges: " << tree << std::endl;
		// join continuous peptide edges starting from a rooted jump point
		utility::vector1< core::Size > jump_points;
		core::kinematics::FoldTree const const_tree( tree );
		for (const auto & it : const_tree) {
			if ( it.start() == (int)jump_root ) jump_points.push_back(it.stop());
		}
		for ( core::Size i=1; i<=jump_points.size(); ++i ) {
			std::set< std::pair< core::Size, core::Size > > remove; // keep track of edges to replace
			tmp_tree.clear();
			core::kinematics::FoldTree const const_new_tree( tree );
			for ( auto it = const_new_tree.begin(), it_end = const_new_tree.end(); it != it_end; ++it ) {
				std::set< std::pair< core::Size, core::Size > > remove_tmp;
				if ( it->start() == (int)jump_points[i] && it->label() == core::kinematics::Edge::PEPTIDE ) {
					core::Size start = it->start();
					core::Size stop = it->stop();
					remove_tmp.insert(std::pair< core::Size, core::Size >( start, stop ));
					for ( auto jt = it+1, jt_end = const_new_tree.end(); jt != jt_end; ++jt ) {
						if ( jt->start() == (int)stop && jt->label() == core::kinematics::Edge::PEPTIDE ) {
							stop = jt->stop();
							remove_tmp.insert(std::pair< core::Size, core::Size >( jt->start(), jt->stop() ));
						}
					}
					tmp_tree.add_edge( start, stop, it->label() );
					if ( remove_tmp.size() > 1 ) remove.insert(remove_tmp.begin(),remove_tmp.end());
				} else {
					tmp_tree.add_edge( it->start(), it->stop(), it->label() );
				}
			}
			core::kinematics::FoldTree const const_tmp_tree( tmp_tree );
			core::kinematics::FoldTree new_tree;
			for (const auto & it : const_tmp_tree) {
				if ( !remove.count(std::pair< core::Size, core::Size >( it.start(), it.stop() )) ) {
					new_tree.add_edge( it.start(), it.stop(), it.label() );
				} else {
					updated = true;
				}
			}
			tree = new_tree;
			if ( updated ) TR << "joined continuous rooted peptide edges: " << tree << std::endl;
		}
	} // strand pairings

	if ( !status ) {
		utility_exit_with_message("HybridizeFoldtreeDynamic: failed to build fold tree from cuts and jumps");
	}

	// Update the pose's fold tree
	core::pose::symmetry::set_asymm_unit_fold_tree( pose , tree );
}


// strand pairings
void HybridizeFoldtreeDynamic::add_overlapping_pair_chunks(
	core::Size const index,
	utility::vector1<core::Size> & cuts,
	utility::vector1<std::pair<core::Size, core::Size> > & jumps,
	std::set<core::Size> & rooted_chunk_indices
) {
	utility::vector1<core::Size> pair_chunks; // pair chunks that pair with root chunk (index)
	utility::vector1<std::pair<core::Size, core::Size> > pair_chunks_pairs; // positions of pair chunks, first is from the root chunk
	get_pair_core_chunks( index, pair_chunks, pair_chunks_pairs );
	utility::vector1<std::pair< core::Size, utility::vector1<core::Size> > > start_pair_chunks;
	utility::vector1<std::pair< core::Size, utility::vector1<std::pair<core::Size, core::Size> > > > start_pair_chunks_pairs;
	if ( pair_chunks.size() ) {
		start_pair_chunks.push_back(std::pair< core::Size, utility::vector1<core::Size> >( index, pair_chunks ));
		start_pair_chunks_pairs.push_back(std::pair< core::Size, utility::vector1<std::pair<core::Size, core::Size> > >( index, pair_chunks_pairs ));
	}
	while ( start_pair_chunks.size() ) {
		utility::vector1< std::pair< core::Size, utility::vector1<core::Size> > > start_pair_chunks_tmp;
		utility::vector1< std::pair< core::Size, utility::vector1<std::pair<core::Size, core::Size> > > > start_pair_chunks_pairs_tmp;
		for ( core::Size i=1; i<=start_pair_chunks.size(); ++i ) {
			//core::Size root = start_pair_chunks[i].first;
			utility::vector1<core::Size>this_pair_chunks = start_pair_chunks[i].second;
			utility::vector1<std::pair<core::Size, core::Size> >this_pair_chunks_pairs = start_pair_chunks_pairs[i].second;
			assert( this_pair_chunks.size() == this_pair_chunks_pairs.size() );
			for ( core::Size j=1; j<=this_pair_chunks.size(); ++j ) {
				if ( rooted_chunk_indices_.count(this_pair_chunks[j]) ) {
					//TR.Debug << "WARNING! Strand pairings between two template chunks overlap. Chunk " << index << " and " << this_pair_chunks[j] << ". Geometry may be funky!" << std::endl;
					continue;
				}
				if ( !rooted_chunk_indices.count(this_pair_chunks[j]) ) {
					// add jump if not rooted already
					TR << "Adding strand pair root " << index << " tree jump: " << this_pair_chunks_pairs[j].first << " " << this_pair_chunks_pairs[j].second << std::endl;
					jumps.push_back(this_pair_chunks_pairs[j]);
					TR << "Adding strand pair root " << index << " tree cut: " << chunks_[this_pair_chunks[j]].stop() << std::endl;
					cuts.push_back(chunks_[this_pair_chunks[j]].stop());
					rooted_chunk_indices.insert(this_pair_chunks[j]);
					utility::vector1<core::Size> new_pair_chunks;
					utility::vector1<std::pair<core::Size, core::Size> > new_pair_chunks_pairs;
					get_pair_core_chunks( this_pair_chunks[j], new_pair_chunks, new_pair_chunks_pairs );
					if ( new_pair_chunks.size() ) {
						start_pair_chunks_tmp.push_back(std::pair< core::Size, utility::vector1<core::Size> >( this_pair_chunks[j], new_pair_chunks ));
						start_pair_chunks_pairs_tmp.push_back(std::pair< core::Size, utility::vector1<std::pair<core::Size, core::Size> > >( this_pair_chunks[j], new_pair_chunks_pairs ));
					}
				}
			}
		}
		start_pair_chunks = start_pair_chunks_tmp;
		start_pair_chunks_pairs = start_pair_chunks_pairs_tmp;
	}
}

void HybridizeFoldtreeDynamic::remove_cut( core::Size const cut, utility::vector1<core::Size> & cuts ) {
	utility::vector1<core::Size> cuts_tmp;
	for ( core::Size i=1; i<=cuts.size(); ++i ) {
		if ( cuts[i] != cut ) cuts_tmp.push_back(cuts[i]);
	}
	cuts = cuts_tmp;
}

void HybridizeFoldtreeDynamic::set_core_chunks(const protocols::loops::Loops & chunks) {
	core_chunks_ = chunks;
}

}  //  namespace hybridization
}  //  namespace protocols
