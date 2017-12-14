// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/evolution/AlignmentAAFinderFilter.cc
/// @brief
/// @author Christoffer Norn


//Unit Headers
#include <protocols/evolution/AlignmentAAFinderFilter.hh>
#include <protocols/evolution/AlignmentAAFinderFilterCreator.hh>
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
#include <protocols/evolution/AlignmentCleanerTools.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
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

static basic::Tracer TR( "protocols.evolution.AlignmentAAFinder" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP AlignmentAAFinderFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AlignmentAAFinder ); }

// XRW TEMP std::string
// XRW TEMP AlignmentAAFinderFilterCreator::keyname() const { return "AlignmentAAFinder"; }

//default ctor
AlignmentAAFinder::AlignmentAAFinder() :
	protocols::filters::Filter( "AlignmentAAFinder" ),
	exclude_AA_threshold_( 10.0 ),
	alignment_file_(""),
	available_AAs_file_(""),
	indel_motif_radius_( 2 ),
	loop_seqid_threshold_( 0.50 ),
	scorefxn_( /* NULL */ ),
	relax_mover_( /* NULL */ )
{
}

AlignmentAAFinder::~AlignmentAAFinder() = default;

void
AlignmentAAFinder::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &, moves::Movers_map const &movers, core::pose::Pose const & )
{
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	exclude_AA_threshold( tag->getOption< core::Real >( "exclude_AA_threshold", 10.0 ) );
	loop_seqid_threshold( tag->getOption< core::Real >( "loop_seqid_threshold", 0.50 ) );
	alignment_file( tag->getOption< std::string >( "alignment_file", "" ) );
	available_AAs_file( tag->getOption< std::string >( "available_AAs_file", "" ) );
	indel_motif_radius( tag->getOption< core::Size >( "indel_motif_radius", 2 ) );
	if ( available_AAs_file() == "" ) {
		utility_exit_with_message( "Outfile pathname not found" );;
	}
	std::string const relax_mover_name( tag->getOption< std::string >( "relax_mover", "null" ) );
	auto mover_it( movers.find( relax_mover_name ) );
	if ( mover_it == movers.end() ) {
		utility_exit_with_message( "Relax mover "+relax_mover_name+" not found" );
	}
	relax_mover( mover_it->second );

}


bool
AlignmentAAFinder::apply( core::pose::Pose const & p ) const {
	//=======================================================================
	// This is a coarse-grained filter to determine which amino acids that
	// are realitically available in the chemical environment given by the
	// target sequence. The criteria are (1) they have to have the same
	// indel-motif* as the target sequence in the alignment, and (2) they
	// have to be at least somewhat capable of fitting in the target
	// sequence background**.
	//
	// *indel motif:    If the primary sequence is
	//                  "--GIHJTS---RJG--QNFOO-----"
	//                  the position in question is no. 12
	//                  the indel motif with radius 2 is:
	//                  X---XXX--X (S---RJG--G)
	//
	// **fit:           If the aln residue is impossible
	//                  to fit in the target sequence
	//                  background (extensive steric overlap)
	//                  it must come from another chemical
	//                  background. Here it is important to
	//                  use a slacked threshold, as the
	//                  structure is not relaxed for the
	//                  residue in question. It might be
	//                  wise to use a soft rep force field
	//                  to take this into account. If one
	//                  has the time, one would also allow
	//                  for a bit of backbone minimization.

	using namespace core::sequence;
	using namespace protocols::evolution;


	// Create a copy of pose, which we can work with
	core::pose::Pose pose( p );

	// load in alignment
	utility::vector1< core::sequence::SequenceOP > aln = core::sequence::read_fasta_file( alignment_file() );

	// Check that the pose sequence matches the first sequence in alignment
	std::string const target_sequence = aln[1]->sequence();
	std::string target_sequence_no_gaps;
	for ( char const& c : target_sequence ) {
		if ( c != '-' ) target_sequence_no_gaps += c;
	}
	std::string const pose_seq = pose.sequence();
	if ( pose_seq != target_sequence_no_gaps ) {
		TR << "Pose sequence " << pose_seq << std::endl;
		TR << "First sequence in alignment " << aln[1]->sequence() << std::endl;
		utility_exit_with_message("The sequence of the pose is different from the first sequence in the alignment");
	}

	// Test that each residue in the pose is sufficently relaxed
	AlignmentCleanerTools::thread_sequence_on_pose( pose, pose_seq, scorefxn() );
	for ( core::Size pose_resi = 1; pose_resi < pose.size() + 1; ++pose_resi ) {
		core::Real const resi_score = pose.energies().residue_total_energy(pose_resi);
		TR << "The energy of pose resi " << pose.sequence()[pose_resi-1] << pose_resi << " in the target sequence is " << resi_score << std::endl;
		if ( resi_score > exclude_AA_threshold() ) {
			utility_exit_with_message("One native residue has a very high score in the target structure. Consider relaxing your input structure or increasing the exclude_AA_threshold");
		}
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

	// Setup map between pose positions and alignment positions
	std::map< core::Size, core::Size > aln2pose;
	int pos_gapped = 0; // aln is 0 numbered
	int pos_pose = 1; // pose is 1 numbered
	for ( char const& c: target_sequence ) {
		if ( c != '-' ) {
			aln2pose[pos_gapped] = pos_pose;
			++pos_pose;
		}
		++pos_gapped;
	}

	// First we find the available amino acids choices to each site
	utility::vector0< std::string > available_aa_identities_vector;
	for ( core::Size aln_resi = 0; aln_resi < target_sequence.length(); ++aln_resi ) {
		core::Size const pose_resi = aln2pose[aln_resi];
		bool const target_position_is_gapped = ( target_sequence[ aln_resi ] == '-');
		if ( target_position_is_gapped ) continue; // nothing to be done here...

		std::string available_aa_identities;
		std::string tested_aa_identities;
		std::string excluded_due_to_indel_motif_diff;
		std::string excluded_due_to_low_loop_seqid;
		core::Size aln_num = 1;
		for ( core::sequence::SequenceOP& aln_seq: aln ) {
			char const aln_seq_resi_type = aln_seq->sequence()[aln_resi];

			// We don't spend time if
			// 1) The indel motif of the aln seq at the aln_resi is different
			// 2) The aln_resi is a loop position and the sequence id is too low.
			// 3) If the we already tested that amino acid at that position
			// 4) If the aln_seq_resi_type is a gap
			std::string aa_indel_motif_target, binary_indel_motif_target;
			std::string aa_indel_motif_aln_seq, binary_indel_motif_aln_seq;
			std::tie(aa_indel_motif_target, binary_indel_motif_target) = AlignmentCleanerTools::indel_motif(target_sequence, indel_motif_radius(), aln_resi, pose_ss_aln );
			std::tie(aa_indel_motif_aln_seq, binary_indel_motif_aln_seq) = AlignmentCleanerTools::indel_motif(aln_seq->sequence(), indel_motif_radius(), aln_resi, pose_ss_aln );
			bool const is_non_matching_indel_motif = ( binary_indel_motif_target != binary_indel_motif_aln_seq );
			bool const is_tested = ( std::find(tested_aa_identities.begin(), tested_aa_identities.end(), aln_seq_resi_type ) != tested_aa_identities.end() );
			bool const is_gap = ( aln_seq_resi_type == '-' );

			if ( is_tested || is_gap ) continue;
			if ( is_non_matching_indel_motif ) {
				bool already_excluded_bc_indel_motif = ( std::find(excluded_due_to_indel_motif_diff.begin(), excluded_due_to_indel_motif_diff.end(), aln_seq_resi_type ) != excluded_due_to_indel_motif_diff.end() );
				if ( !already_excluded_bc_indel_motif ) {
					excluded_due_to_indel_motif_diff += aln_seq_resi_type;
				}
				continue;
			}

			// If the indel motif is matching and aln_resi is
			// in a loop we also filter according to
			// sequence identity
			bool const is_aln_resi_loop = ( pose_ss[aln_resi] == 'L' );
			if ( is_aln_resi_loop ) {
				core::Real seqid = AlignmentCleanerTools::indel_motif_seq_id(aa_indel_motif_target, aa_indel_motif_aln_seq);
				if ( seqid < loop_seqid_threshold() ) {
					bool already_excluded_bc_low_loop_seqid = ( std::find(excluded_due_to_low_loop_seqid.begin(), excluded_due_to_low_loop_seqid.end(), aln_seq_resi_type ) != excluded_due_to_low_loop_seqid.end() );
					if ( !already_excluded_bc_low_loop_seqid ) {
						excluded_due_to_low_loop_seqid += aln_seq_resi_type;
					}
					continue;
				}
			}

			// Test if the amino acid can actually fit in the current sequence background.
			core::pose::Pose pose( p ); // reset pose to the input pose
			std::string thread_seq_gapped = target_sequence;
			thread_seq_gapped[aln_resi] = aln_seq_resi_type;
			std::string thread_seq;
			for ( char const& c: thread_seq_gapped ) {
				if ( c != '-' ) thread_seq += c;
			}
			AlignmentCleanerTools::thread_sequence_on_pose(pose, thread_seq, scorefxn() );
			core::Real const resi_score = pose.energies().residue_total_energy(pose_resi);
			tested_aa_identities += aln_seq_resi_type;

			bool const use_stability_in_target_seq_as_cutoff = true;
			if ( use_stability_in_target_seq_as_cutoff ) {
				if ( resi_score < exclude_AA_threshold() ) {
					TR << "Including AA " << aln_seq_resi_type << " at resi " << pose_resi << " because " << resi_score << " < " << exclude_AA_threshold() << " = exclude_AA_threshold" << std::endl;
					available_aa_identities += aln_seq_resi_type;
				} else {
					TR << "Excluding AA " << aln_seq_resi_type << " at resi " << pose_resi << " because " << resi_score << " > " << exclude_AA_threshold() << " = exclude_AA_threshold" << std::endl;
				}
			} else {
				available_aa_identities += aln_seq_resi_type;
			}
			++aln_num;
		} // aln_seq
		available_aa_identities_vector.push_back( available_aa_identities );

		for ( char const& c: excluded_due_to_indel_motif_diff ) {
			bool const is_included_for_other_seqs = ( std::find(available_aa_identities.begin(), available_aa_identities.end(), c ) != available_aa_identities.end() );
			if ( !is_included_for_other_seqs ) {
				TR << "Excluding AA " << c << " at resi " << pose_resi << " because non-matching indel-motif" << std::endl;
			}
		}
		for ( char const& c: excluded_due_to_low_loop_seqid ) {
			bool const is_included_for_other_seqs = ( std::find(available_aa_identities.begin(), available_aa_identities.end(), c ) != available_aa_identities.end() );
			if ( !is_included_for_other_seqs ) {
				TR << "Excluding AA " << c << " at resi " << pose_resi << " because low seq-id of containing loop" << std::endl;
			}
		}
		//TR << "available_aa_identities_vector is now " << available_aa_identities_vector.size() << "long " << std::endl;
	} // aln_resi

	// Write out info
	//    std::ifstream infile( available_AAs_file().c_str() );
	//    infile.close();
	std::ofstream outfile;
	outfile.open( available_AAs_file().c_str() );
	TR << "The per site available amino acids are" << std::endl;
	core::Size pose_resi = 1;
	for ( std::string& avail_ids: available_aa_identities_vector ) {
		TR << pose_resi << " " << avail_ids << std::endl;
		outfile << pose_resi << " " << avail_ids << '\n';
		++pose_resi;
	}
	outfile.close();

	return( true );
}

void
AlignmentAAFinder::report( std::ostream &, core::pose::Pose const & ) const
{
}

core::Real
AlignmentAAFinder::report_sm( core::pose::Pose const & ) const {
	return( 1 );
}

std::string AlignmentAAFinder::name() const {
	return class_name();
}

std::string AlignmentAAFinder::class_name() {
	return "AlignmentAAFinder";
}

void AlignmentAAFinder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "exclude_AA_threshold", xsct_real, "How large may the score be of the an amino acid in the target sequence", "10")
		+ XMLSchemaAttribute::attribute_w_default( "alignment_file", xs_string, "alignment_file to be cleaned", "null" )
		+ XMLSchemaAttribute::attribute_w_default( "relax_mover", xs_string, "Relax mover", "null" )
		+ XMLSchemaAttribute::attribute_w_default( "loop_seqid_threshold", xsct_real, "Sequence identity threshold of indel motif matching loops", "0.50")
		+ XMLSchemaAttribute::attribute_w_default( "available_AAs_file", xs_string, "output path", "null" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Cleans an alignment, such that all amino acids are measured in the same chemical environment",
		attlist );
}

std::string AlignmentAAFinderFilterCreator::keyname() const {
	return AlignmentAAFinder::class_name();
}

protocols::filters::FilterOP
AlignmentAAFinderFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new AlignmentAAFinder );
}

void AlignmentAAFinderFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignmentAAFinder::provide_xml_schema( xsd );
}


}
}
