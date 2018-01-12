// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/evolution/AlignmentGapInserterFilter.cc
/// @brief
/// @author Christoffer Norn


//Unit Headers
#include <protocols/evolution/AlignmentGapInserterFilter.hh>
#include <protocols/evolution/AlignmentGapInserterFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <map>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <string>


//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <protocols/evolution/AlignmentCleanerTools.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/graph/Graph.hh>



#include <utility/vector1.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace evolution {

using bools = utility::vector1<bool>;
using namespace std;
using namespace core::scoring;

static basic::Tracer TR( "protocols.evolution.AlignmentGapInserter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP AlignmentGapInserterFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AlignmentGapInserter ); }

// XRW TEMP std::string
// XRW TEMP AlignmentGapInserterFilterCreator::keyname() const { return "AlignmentGapInserter"; }

//default ctor
AlignmentGapInserter::AlignmentGapInserter() :
	protocols::filters::Filter( "AlignmentGapInserter" ),
	alignment_file_(""),
	available_AAs_file_(""),
	cleaned_alignment_file_(""),
	nbr_e_threshold_( 0.1 ),
	indel_motif_radius_( 2 ),
	only_clean_seq_num_( 1 ),
	scorefxn_( /* NULL */ )
{
	max_score_diffs_.clear();
	loop_seqid_thresholds_.clear();

}

AlignmentGapInserter::~AlignmentGapInserter() = default;

void
AlignmentGapInserter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const &  )
{
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );

	max_score_diffs( utility::string_split( tag->getOption< std::string >( "max_score_diffs" ), ',', core::Real() ) );
	alignment_file( tag->getOption< std::string >( "alignment_file", "" ) );
	available_AAs_file( tag->getOption< std::string >( "available_AAs_file", "" ) );
	cleaned_alignment_file( tag->getOption< std::string >( "cleaned_alignment_file", "" ) );
	nbr_e_threshold( tag->getOption< core::Real >( "nbr_e_threshold", 0.1 ) );
	loop_seqid_thresholds( utility::string_split( tag->getOption< std::string >( "loop_seqid_thresholds" ), ',', core::Real() ) );
	indel_motif_radius( tag->getOption< core::Size >( "indel_motif_radius", 2 ) );
	only_clean_seq_num( tag->getOption< core::Size >( "only_clean_seq_num", 1 ) );

	if ( cleaned_alignment_file() == "" ) {
		utility_exit_with_message( "Outfile pathname not found" );;
	}

	TR << "Neighboring residues defined by res-res energy of at least: nbr_e_threshold = " << nbr_e_threshold() << std::endl;
	TR << "Energy threshold(s) to discriminate different chemical environments is: exclude_site_threshold = " << max_score_diffs() << std::endl;
	TR << "Seq id threshold(s) by which loops are excluded due to expected other conformation = " << loop_seqid_thresholds() << std::endl;

}

utility::vector0< core::Size >
AlignmentGapInserter::find_char_location_in_string(std::string const & string, char const findIt) const
{
	utility::vector1<core::Size> characterLocations;
	for ( core::Size i = 0; i < string.size(); i++ ) {
		if ( string[i] == findIt ) {
			characterLocations.push_back( i );
		}
	}
	return characterLocations;
}

bool
AlignmentGapInserter::apply( core::pose::Pose const & p ) const {
	//=======================================================================
	///// Clean-up the alignment, so it only contains chemical-environment
	///// relevant identities
	//=======================================================================
	// Method for accepting a position as being chemically equivalent to
	// position in the sequence environment of the target sequence:
	// 1) No nbr position* should have a neighboring +/-2 indel
	//    motif (xxPxx) in that differs from that of the target seq.
	// 2) All available amino acids** to that site in the aligment
	//    should score the same***, as in the target sequence.
	//
	// Definitions:
	// *   nbr position: Any residues within 5 AA.
	// **  Available amino acids: Residues with an energy of more than
	//     10 REU are not considered compatible with the backbone. The
	//     logic is that, if it doesn't fit in the structure, we might
	//     as well exclude from the alignment.
	// *** Score the same: max(e_vec) < 1 REU.


	using namespace core::sequence;
	using namespace protocols::evolution;

	// Create a copy of pose, which we can work with
	core::pose::Pose pose( p );

	// load in alignment
	utility::vector1< core::sequence::SequenceOP > aln_all = core::sequence::read_fasta_file( alignment_file() );

	// Check that the pose sequence matches the first sequence in alignment
	std::string const target_sequence = aln_all[1]->sequence();
	std::string target_sequence_no_gaps;
	for ( char const& c : target_sequence ) {
		if ( c != '-' ) {
			target_sequence_no_gaps += c;
		}
	}
	std::string const pose_seq = pose.sequence();
	if ( pose_seq != target_sequence_no_gaps ) {
		TR << "Pose sequence " << pose_seq << std::endl;
		TR << "First sequence in alignment " << aln_all[1]->sequence() << std::endl;
		utility_exit_with_message("The sequence of the pose is different from the first sequence in the alignment");
	}

	// Determine the secondary structure of the pose
	// map this onto the alignment of the target
	// sequence
	core::scoring::dssp::Dssp dssp( pose );
	core::Size const min_loop_length = 2;
	std::string const pose_ss = AlignmentCleanerTools::short_ss_loop_filter( dssp.get_dssp_secstruct(), min_loop_length );
	std::string pose_ss_aln;
	core::Size ss_counter = 0;
	for ( char i : target_sequence ) {
		if ( i == '-' || ss_counter > pose_ss.size() ) pose_ss_aln += '-';
		else {
			pose_ss_aln += pose_ss[ss_counter];
			++ss_counter;
		}
	}

	// Slim down alignment with two criteria:
	// (1)  User-specified focus on one alignment sequence (only_clean_seq_num)
	// (2)  We don't want to add sequences that cannot add
	//      information. For example
	//      target_seq   HWIFSAH--IFFA
	//      aln_seq      HWI--AHAFIFFA

	utility::vector1< core::sequence::SequenceOP > aln;
	aln.push_back( aln_all[1] ); // we always want to target sequence in there
	for ( core::Size aln_num = 2; aln_num<=aln_all.size(); ++aln_num ) {
		if ( only_clean_seq_num() != 1 && only_clean_seq_num() != aln_num ) {
			continue; // for large alignments it is useful to be able to run through sequences individually
		}
		std::string const aln_seq = aln_all[aln_num]->sequence();
		std::string aln_seq_no_relative_gaps;
		std::string target_seq_no_relative_gaps;
		for ( core::Size i=0; i < target_sequence.length(); ++i ) {
			if ( target_sequence[i] != '-' && aln_seq[i] != '-' ) {
				aln_seq_no_relative_gaps += aln_seq[i];
				target_seq_no_relative_gaps += target_sequence[i];
			}
		}
		bool is_adding_info = ( aln_seq_no_relative_gaps != target_seq_no_relative_gaps );
		if ( is_adding_info ) aln.push_back( aln_all[aln_num] );
		else {
			TR << "Sequence ID: " << aln_all[aln_num]->id() << " is not adding information." << std::endl;
		}
	}

	core::Size const n_removed_seqs = aln_all.size() - aln.size();
	if ( n_removed_seqs > 0 ) {
		TR << "Removed " << n_removed_seqs << " sequences from input alignment (size=" << aln_all.size() << "), as they couldn't contribute information. " << std::endl;
	}
	if ( aln.size() == 1 ) {
		utility_exit_with_message("There are no sequences in the alignment, which can add information");
	}

	// Read the available amino acids file
	utility::vector1< std::string > available_aa_identities_vector;
	std::ifstream infile( available_AAs_file() );
	core::Size resi_num;
	std::string avail_aas;
	while ( infile >> resi_num >> avail_aas ) {
		available_aa_identities_vector.push_back( avail_aas );
	}
	if ( available_aa_identities_vector.size() != pose_seq.size() ) {
		TR << "Pose sequence is " << pose_seq.size() << " aa long " << std::endl;
		TR << "The avail aas file is " << available_aa_identities_vector.size() << " aa long " << std::endl;
		utility_exit_with_message("The avail aas file must have the same length as the pose");
	}
	TR << "The per site available amino acids are" << std::endl;
	for ( core::Size pose_resi = 1; pose_resi < pose_seq.length() + 1; ++pose_resi ) {
		TR << pose_resi << " " << available_aa_identities_vector[pose_resi] << std::endl;
	}

	// Make a copy of the alignment, which will be our cleaned alignment
	// Since we might use a range of energy cutoff and loop seq id
	// values for how we are defining the chemical neighboring environment
	// we need multiple aln_cleans
	std::map< core::Real, std::map< core::Real, utility::vector1< core::sequence::SequenceOP > > > aln_clean_all;
	for ( core::Real &e: max_score_diffs() ) {
		for ( core::Real &s: loop_seqid_thresholds() ) {
			utility::vector1< core::sequence::SequenceOP > aln_clean;
			for ( core::sequence::SequenceOP& seq: aln ) {
				aln_clean.push_back(seq->clone());
			}
			aln_clean_all[e][s] = aln_clean;
		}
	}

	// Setup map between pose positions and alignment positions
	std::map< core::Size, core::Size > aln2pose;
	std::map< core::Size, core::Size > pose2aln;
	core::Size pos_gapped = 0; // aln is 0 numbered
	core::Size pos_pose = 1; // pose is 1 numbered
	for ( char const& c: target_sequence ) {
		if ( c != '-' ) {
			aln2pose[pos_gapped] = pos_pose;
			pose2aln[pos_pose] = pos_gapped;
			++pos_pose;
		}
		++pos_gapped;
	}

	// Next we go over all positions in the target sequence and filter the columns in the alignment
	for ( core::Size aln_resi = 0; aln_resi < target_sequence.length(); ++aln_resi ) {
		//TR << "" << std::endl;
		//TR << "Cleaning alignment position " << aln_resi << std::endl;

		// If the target sequence has a gap, there is no information here. Insert a gap also in the aln seq
		char const pose_resn = target_sequence[aln_resi];
		if ( pose_resn == '-' ) {
			for ( core::Size aln_num = 1; aln_num<=aln.size(); ++aln_num ) {
				for ( core::Real& e: max_score_diffs() ) {
					for ( core::Real &s: loop_seqid_thresholds() ) {
						aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
					}
				}
			}
			TR << "Excluding position " << aln_resi << " for all aligned sequences bc target sequence is gapped" << std::endl;
			continue;
		}

		// Next go over all sequences in the alignment exept the target sequence
		// for which there is nothing to exclude...
		utility::vector1< core::Real > target_energy_profile; // This is only pushed to for the first sequence in the alignment
		core::Size const pose_resi = aln2pose[aln_resi];
		std::string const avail_aas = available_aa_identities_vector[pose_resi];
		for ( core::Size aln_num = 1; aln_num <= aln.size(); ++aln_num ) {
			std::string const aln_seq = aln[aln_num]->sequence();

			// If there is already a gap in the aln_seq, there is nothing to exclude..
			char const aln_seq_resi_type = aln_seq[aln_resi];
			bool const is_gap = ( aln_seq_resi_type == '-');
			if ( is_gap ) {
				TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = " <<  pose_resi << pose_resn << ") in aln " << aln_num << " bc already gap" << std::endl;
				continue;
			}

			// If the aln_seq resi is not allowed at the position, we insert a gap and continue
			bool const is_allowed = ( std::find(avail_aas.begin(), avail_aas.end(), aln_seq_resi_type ) != avail_aas.end() );
			if ( !is_allowed ) {
				TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = " <<  pose_resi << pose_resn << ") in aln " << aln_num << " bc disallowed aa type" << std::endl;
				for ( core::Real& e: max_score_diffs() ) {
					for ( core::Real &s: loop_seqid_thresholds() ) {
						aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
					}
				}
				continue;
			}

			// If the indel motif is wrong in the primary seq
			// we insert a gap and continue
			std::string aa_indel_motif_target, binary_indel_motif_target;
			std::string aa_indel_motif_aln_seq, binary_indel_motif_aln_seq;
			std::tie(aa_indel_motif_target, binary_indel_motif_target) = AlignmentCleanerTools::indel_motif(target_sequence, indel_motif_radius(), aln_resi, pose_ss_aln );
			std::tie(aa_indel_motif_aln_seq, binary_indel_motif_aln_seq) = AlignmentCleanerTools::indel_motif(aln_seq, indel_motif_radius(), aln_resi, pose_ss_aln);
			bool is_non_matching_indel_motif = ( binary_indel_motif_target != binary_indel_motif_aln_seq );
			if ( is_non_matching_indel_motif ) {
				TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = " <<  pose_resi << pose_resn << ") in aln " << aln_num << " bc non-matching indel motif in primary seq" << std::endl;
				for ( core::Real& e: max_score_diffs() ) {
					for ( core::Real &s: loop_seqid_thresholds() ) {
						aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
					}
				}
				continue;
			}

			// If the indel motif is matching and aln_resi is
			// in a loop we also filter according to
			// sequence identity
			// (One optimization here, could be that we only filter
			// based on the loop residues...)
			bool const is_aln_resi_loop = ( pose_ss[aln_resi] == 'L' );
			if ( is_aln_resi_loop ) {
				core::Real seqid_of_loop = AlignmentCleanerTools::indel_motif_seq_id(aa_indel_motif_target, aa_indel_motif_aln_seq);
				core::Size gap_insert_counter = 0;
				for ( core::Real &s: loop_seqid_thresholds() ) {
					if ( seqid_of_loop < s ) {
						++gap_insert_counter;
						TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = " <<  pose_resi << pose_resn << ") in aln " << aln_num << " bc loop seq id is " << seqid_of_loop << " < " << s << std::endl;
						for ( core::Real &e: max_score_diffs() ) {
							aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
						}
					}
				}
				bool const gap_inserted_for_all = ( gap_insert_counter ==  loop_seqid_thresholds().size());
				if ( gap_inserted_for_all ) continue;
			}

			// Evaluate in detail the chemical environment at site
			utility::vector1< core::Real > candidate_energy_profile;
			bool gap_inserted = false;
			for ( char const& avail_resn: avail_aas ) {
				if ( gap_inserted ) break;

				// Prepare the aln_seq for threading: Three steps:
				// (1) Remove all the gaps in the aln_seq, where there are also
				//     gaps in the target_sequence
				std::string thread_seq;
				for ( core::Size i=0; i < target_sequence.length(); ++i ) {
					if ( target_sequence[i] != '-' ) thread_seq += aln_seq[i];
				}
				// (2) Insert the avail aa, for which we want to test the chemical
				//     environment
				thread_seq[pose_resi - 1] = avail_resn;

				// (3) Mutate remaining gaps in aln_seq to glycines.
				utility::vector0<core::Size> gap_positions = find_char_location_in_string(thread_seq, '-');
				for ( core::Size& gap_pos: gap_positions ) {
					thread_seq[ gap_pos ] = 'G'; // mutate to glycine
				}
				// Thread on the sequence
				core::pose::Pose pose( p ); // reset the pose to the input pose
				AlignmentCleanerTools::thread_sequence_on_pose( pose, thread_seq, scorefxn() );

				//std::stringstream ss;
				//std::string outname;
				//ss << "dumps/pose_resi_" << pose_resi << "_avail_resn_" << avail_resn << ".pdb";
				//ss >> outname;
				//pose.dump_scored_pdb(outname, *scorefxn());

				// Find the neighbors and insert gap if:
				// (1) Neighbor is gapped residue
				// (2) Neighbor residues has different indel motif
				// (3) Neighbor is in a low loop seqid environment
				utility::vector1< core::Size > e_nbrs;
				for ( utility::graph::EdgeListConstIterator egraph_it = pose.energies().energy_graph().get_node( pose_resi )->edge_list_begin();
						egraph_it != pose.energies().energy_graph().get_node( pose_resi )->edge_list_end(); ++egraph_it ) {
					core::Size const other_pose_res = (*egraph_it)->get_other_ind( pose_resi );
					core::Size const other_res_aln = pose2aln[other_pose_res];
					auto const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);
					EnergyMap const cur_weights = pose.energies().weights();
					core::Real const resresE( Eedge->dot_abs( cur_weights ) );
					bool const is_not_res_nbr = ( resresE < nbr_e_threshold() );
					if ( is_not_res_nbr ) continue;

					std::tie(aa_indel_motif_target, binary_indel_motif_target) = AlignmentCleanerTools::indel_motif(target_sequence, indel_motif_radius(), other_res_aln, pose_ss_aln );
					std::tie(aa_indel_motif_aln_seq, binary_indel_motif_aln_seq) = AlignmentCleanerTools::indel_motif(aln_seq, indel_motif_radius(), other_res_aln, pose_ss_aln);
					bool const is_non_matching_indel_motif = ( binary_indel_motif_target != binary_indel_motif_aln_seq );
					if ( is_non_matching_indel_motif ) {
						TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = "
							<< pose_resi << pose_resn << " ) in aln " << aln_num << " bc nbr-residue "
							<< other_res_aln << " (pose resi = " << other_pose_res << " ) has non-matching indel motif"
							<< std::endl;
						for ( core::Real& e: max_score_diffs() ) {
							for ( core::Real &s: loop_seqid_thresholds() ) {
								aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
							}
						}
						gap_inserted = true;
						break;
					}
					bool const nbr_is_gapped = ( aln_seq[other_res_aln] == '-' );
					if ( nbr_is_gapped && !is_non_matching_indel_motif ) {
						TR << "For debugging:" << std::endl;
						TR << "aa_indel_motif_target " << aa_indel_motif_target << std::endl;
						TR << "binary_indel_motif_target " << binary_indel_motif_target << std::endl;
						TR << "aa_indel_motif_aln_seq " << aa_indel_motif_aln_seq << std::endl;
						TR << "binary_indel_motif_aln_seq " << binary_indel_motif_aln_seq << std::endl;
						utility_exit_with_message( "This should be impossible, I think... If not we need to fix the code a bit" );
					}

					bool const is_other_res_loop = ( pose_ss[other_res_aln] == 'L' );
					core::Size gap_insert_counter = 0;
					core::Real const seqid = AlignmentCleanerTools::indel_motif_seq_id(aa_indel_motif_target, aa_indel_motif_aln_seq);
					for ( core::Real &s: loop_seqid_thresholds() ) {
						bool const is_low_seq_loop = ( is_other_res_loop && seqid < s ); //loop_seqid_thresholds());
						if ( is_low_seq_loop ) {
							++gap_insert_counter;

							// When we have been here once, we know that for that loop seq id threshold, at least
							// one of the avail aas is in a bad environment. Then we don't need to print it anymore.
							// However we can't break out of the seqid-for-loop, as another higher seq id threshold
							// might allow the environment. Moving the seq id for loop outside the avail aas is
							// possible but would prolong runtime, as we would need to thread more. Instead we
							// check whether this is the first time we are inserting a gap here. If so, we print
							// exclusion information -- otherwise not.
							bool first_gap_insertion = false;
							for ( core::Real& e: max_score_diffs() ) {
								if ( aln_clean_all[e][s][aln_num]->at(aln_resi + 1) != '-' ) first_gap_insertion = true;
								aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
							}

							if ( first_gap_insertion ) {
								TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = "
									<< pose_resi << pose_resn << " ) in aln " << aln_num << " bc nbr-residue "
									<< other_res_aln << " (pose resi = " << other_pose_res << " ) is in a low seq id loop. Seqid: "
									<< seqid << " < " << s << " based on avail aa: " << avail_resn << std::endl;
							}
						}
					}
					bool const gap_inserted_for_all = ( gap_insert_counter ==  loop_seqid_thresholds().size());
					if ( gap_inserted_for_all ) gap_inserted = true;
				}
				if ( gap_inserted ) break;

				// Measure the energy vec for the site.
				core::Real resi_score = pose.energies().residue_total_energy(pose_resi);
				if ( aln_num == 1 ) target_energy_profile.push_back(resi_score);
				else candidate_energy_profile.push_back(resi_score);
			} // avail_resn

			if ( aln_num > 1 && !gap_inserted ) {
				//TR << "For resi " << pose_resi << "    target energy vector is " << target_energy_profile << std::endl;
				//TR << "For resi " << pose_resi << " candidate energy vector is " << candidate_energy_profile << std::endl;
				std::vector< core::Real > diffs;
				for ( core::Size i=1; i<=target_energy_profile.size(); ++i ) {
					diffs.push_back(sqrt((target_energy_profile[i]-candidate_energy_profile[i])*(target_energy_profile[i]-candidate_energy_profile[i])));
				}
				core::Real max_diff = *std::max_element(diffs.begin(), diffs.end());
				for ( core::Real& e: max_score_diffs() ) {
					if ( max_diff > e ) {
						TR << "Excluding position " << aln_resi << aln_seq_resi_type << " (pose_resi = " <<  pose_resi << pose_resn << ") in aln " << aln_num << " bc different chemical environment. Site delta residue energy vector = " << diffs << " threshold = " << e << std::endl;
						for ( core::Real &s: loop_seqid_thresholds() ) {
							aln_clean_all[e][s][aln_num]->replace_char(aln_resi + 1, '-'); // the sequence object is 1 numbered
						}
					} else {
						TR << "Including position " << aln_resi << aln_seq_resi_type << " (pose_resi = " <<  pose_resi << pose_resn << ") in aln " << aln_num << " bc similar chemical environment. Site delta residue energy vector = " << diffs << " threshold = " << e << std::endl;
					}
				}
			}
		} // aln number
	} // pose_resi

	// Write the cleaned alignment with gaps
	for ( core::Real &e: max_score_diffs() ) {
		for ( core::Real &s: loop_seqid_thresholds() ) {
			TR << "Writing the cleaned alignment with gaps for same chemical environment threshold of " << e << " and loop seqid threshold of " << s << std::endl;
			std::ofstream cleanedAlignFile;
			std::stringstream ss;
			ss << cleaned_alignment_file() << "_eThreshold_" << e << "_sThreshold_" << s;
			cleanedAlignFile.open( ss.str().c_str() );
			for ( core::sequence::SequenceOP& seq: aln_clean_all[e][s] ) { // this is also printing the first element... perhaps we don't want this?
				cleanedAlignFile << ">" << seq->id() << "\n" << seq->sequence() << "\n"; //
			}
			cleanedAlignFile.close();
		}
	}


	// Write the cleaned alignment for non-gapped positions in target sequence
	for ( core::Real &e: max_score_diffs() ) {
		for ( core::Real &s: loop_seqid_thresholds() ) {
			TR << "Writing the cleaned alignment with target seq gaps to file" << std::endl;
			std::ofstream cleanedAlignFileNoGaps;
			std::stringstream ss_no_gaps;
			ss_no_gaps << cleaned_alignment_file() << "_eThreshold_" << e << "_sThreshold_" << s << "_no_gaps";
			cleanedAlignFileNoGaps.open( ss_no_gaps.str().c_str() );
			for ( core::sequence::SequenceOP& seq: aln_clean_all[e][s] ) {
				std::string out_seq;
				std::string gapped_seq = seq->sequence();
				for ( core::Size i=0; i < target_sequence.length(); ++i ) {
					if ( target_sequence[i] != '-' ) out_seq += gapped_seq[i];
				}
				cleanedAlignFileNoGaps << ">" << seq->id() << "\n" << out_seq << "\n";
				TR << ">" << seq->id() << std::endl;
				TR << out_seq << std::endl;
			}
			cleanedAlignFileNoGaps.close();
		}
	}
	return( true );
}

void
AlignmentGapInserter::report( std::ostream &, core::pose::Pose const & ) const
{
}

core::Real
AlignmentGapInserter::report_sm( core::pose::Pose const & ) const {
	return( 1 );
}

std::string AlignmentGapInserter::name() const {
	return class_name();
}

std::string AlignmentGapInserter::class_name() {
	return "AlignmentGapInserter";
}

void AlignmentGapInserter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute( "max_score_diffs", xsct_real_cslist, "Comma-separated list of exclude site thresholds" )
		+ XMLSchemaAttribute::attribute_w_default( "alignment_file", xs_string, "alignment_file to be cleaned", "null" )
		+ XMLSchemaAttribute::attribute_w_default( "only_clean_seq_num", xsct_positive_integer, "Sequence number in the alignment to filter. If unspecifed all will be filtered", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "cleaned_alignment_file", xs_string, "output file", "null" )
		+ XMLSchemaAttribute::attribute_w_default( "available_AAs_file", xs_string, "input file", "null" )
		+ XMLSchemaAttribute::attribute_w_default( "nbr_e_threshold", xsct_real, "res res energy for neighbor classification", "0.1")
		+ XMLSchemaAttribute::attribute_w_default( "loop_seqid_thresholds", xsct_real_cslist, "Sequence identity thresholds of indel motif matching loops", "0.65")
		+ XMLSchemaAttribute::attribute_w_default( "indel_motif_radius", xsct_non_negative_integer, "If indel motif radius is 2 then 2 residues will be included up and downstream from target position", "2");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Cleans an alignment, such that all amino acids are measured in the same chemical environment",
		attlist );
}

std::string AlignmentGapInserterFilterCreator::keyname() const {
	return AlignmentGapInserter::class_name();
}

protocols::filters::FilterOP
AlignmentGapInserterFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new AlignmentGapInserter );
}

void AlignmentGapInserterFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignmentGapInserter::provide_xml_schema( xsd );
}


}
}
