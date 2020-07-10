// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pdbinfo_manipulations/AlignPDBInfoToSequences.cc
/// @brief Realign poses to sequences after losing pdb_info
/// @author Dan Farrell (danpf@uw.edu)

// Unit headers
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequences.hh>
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesCreator.hh>
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequencesUtil.hh>

// Core headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SWAligner.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>
#include <utility/io/izstream.hh>
#include <utility/json_utilities.hh>
#include <utility/excn/Exceptions.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <numeric>
#include <sstream>


static basic::Tracer TR( "protocols.pdbinfo_manipulations.AlignPDBInfoToSequences" );

namespace protocols {
namespace pdbinfo_manipulations {

AlignPDBInfoToSequencesMode
get_alignpdbinfotosequencesmode_from_string(std::string const & str) {
	if ( str == "single" ) return AlignPDBInfoToSequencesMode::single;
	else if ( str == "multiple" ) return AlignPDBInfoToSequencesMode::multiple;
	else {
		std::stringstream ss;
		ss << "Bad get_alignpdbinfotosequencesmode_from_string, allowed: ['single', 'multiple'] found " << str
			<< std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput, ss.str() );
	}
}

std::string
alignpdbinfotosequencesmode_to_string(AlignPDBInfoToSequencesMode const & mode) {
	switch (mode) {
	case AlignPDBInfoToSequencesMode::single :
		return "single";
	case AlignPDBInfoToSequencesMode::multiple :
		return "multiple";
	default :
		std::stringstream ss;
		ss << "Bad alignpdbinfotosequencesmode_to_string, did not find 'single' or 'multiple'"
			<< std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput, ss.str() );
	}
}


///@brief make sure none of the non-sequence elements are larger than
///       the sequence_.
void
SequenceSpecification::check_single_format() const {
	// 1. check for bad settings
	if ( chains_.size() > sequence_.size() ) {
		std::stringstream ss;
		ss << "SequenceSpecification chains size was larger than the sequence c: "
			<< chains_.size() << " > s: " << sequence_.size() << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}

	if ( segmentIDs_.size() > sequence_.size() ) {
		std::stringstream ss;
		ss << "SequenceSpecification segmentIDs size was larger than the sequence seg: "
			<< segmentIDs_.size() << " > s: " << sequence_.size() << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}

	if ( insCodes_.size() > sequence_.size() ) {
		std::stringstream ss;
		ss << "SequenceSpecification insCodes size was larger than the sequence ins: "
			<< insCodes_.size() << " > s: " << sequence_.size() << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}

	if ( residue_numbers_.size() > sequence_.size() ) {
		std::stringstream ss;
		ss << "SequenceSpecification residue_numbers size was larger than the sequence ins: "
			<< residue_numbers_.size() << " > s: " << sequence_.size() << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}
}

///@brief make sure none of the non-sequence elements are larger than
///       the sequence_.
void
SequenceSpecification::update_multiple_format_residue_numbering() {
	if ( residue_numbers_.empty() ) {
		for ( int i = 1; i <= int(sequence_.size()); ++i ) {
			residue_numbers_.push_back(i);
		}
	} else if ( residue_numbers_.size() == core::Size(1) ) {
		// skip seed number
		for ( int i = 1; i < int(sequence_.size()); ++i ) {
			residue_numbers_.push_back(i+residue_numbers_[1]);
		}
	} else if ( residue_numbers_.size() != sequence_.size() ) {
		std::stringstream ss;
		ss << "SequenceSpecification residue_numbers set and found to not have the same size"
			<< " as the sequence (they must have the same size if set): "
			<< residue_numbers_.size() << " > s: " << sequence_.size() << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}
}


///@brief duplicate the last element (or default last element if not set) to
///       pad chains_, segmentIDs_, and insCodes_ to be the same length as
///       the sequence_.
void
SequenceSpecification::set_1_to_1_format(int const starting_number) {
	while ( chains_.size() <= sequence_.size() ) chains_.push_back(chains_.back());
	if ( segmentIDs_.size() == 0 ) segmentIDs_.push_back("    ");
	while ( segmentIDs_.size() <= sequence_.size() ) segmentIDs_.push_back(segmentIDs_.back());
	if ( insCodes_.size() == 0 ) insCodes_.push_back(" ");
	while ( insCodes_.size() <= sequence_.size() ) insCodes_.push_back(insCodes_.back());
	if ( residue_numbers_.size() == 0 ) {
		for ( int i = 0; i < int(sequence_.size()); ++i ) {
			residue_numbers_.push_back(starting_number + i);
		}
	}
}

bool
SequenceSpecification::operator==(SequenceSpecification const & alt_ss) const {
	return std::tie(sequence_, chains_, segmentIDs_, insCodes_, residue_numbers_) == std::tie(alt_ss.sequence_, alt_ss.chains_, alt_ss.segmentIDs_, alt_ss.insCodes_, alt_ss.residue_numbers_);
}




std::ostream &
operator<<( std::ostream & os, SequenceSpecification const & ss ) {
	os << "SequenceSpecification:\n"
		<< "sequence: " << ss.sequence_.size() << " : " << ss.sequence_ << "\n"
		<< "chains: " << ss.chains_.size() << " : " << ss.chains_ << "\n"
		<< "insCodes_: " << ss.insCodes_.size() << " : " << ss.insCodes_ << "\n"
		<< "segmentIDs_: " << ss.segmentIDs_.size() << " : " << ss.segmentIDs_ << "\n"
		<< "residue_numbers_: " << ss.residue_numbers_.size() << " : " << ss.residue_numbers_ << "\n"
		<< "current_idx: " << ss.current_idx_ << std::endl;
	return os;
}


///@brief serialization (json) of SequenceSpecification
void to_json(nlohmann::json& j, SequenceSpecification const & ss) {
	j = nlohmann::json{
		{"sequence", ss.sequence_},
		{"chains", ss.chains_.vector()},
		{"insCodes", ss.insCodes_.vector()},
		{"segmentIDs", ss.segmentIDs_.vector()},
		{"residue_numbers", ss.residue_numbers_.vector()},
		{"current_idx", ss.current_idx_}
		};
}

void from_json(const json& j, SequenceSpecification & ss) {
	j.at("sequence").get_to(ss.sequence_);
	j.at("chains").get_to(ss.chains_);
	if ( j.find("insCodes") != j.end() ) j.at("insCodes").get_to(ss.insCodes_);
	if ( j.find("segmentIDs") != j.end() ) j.at("segmentIDs").get_to(ss.segmentIDs_);
	if ( j.find("current_idx") != j.end() ) j.at("current_idx").get_to(ss.current_idx_);
	if ( j.find("residue_numbers") != j.end() ) j.at("residue_numbers").get_to(ss.residue_numbers_);
}



/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
AlignPDBInfoToSequences::AlignPDBInfoToSequences():
	protocols::moves::Mover( AlignPDBInfoToSequences::mover_name() ),
	mode_(AlignPDBInfoToSequencesMode::unset), throw_on_fail_(false)
{
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
AlignPDBInfoToSequences::AlignPDBInfoToSequences( AlignPDBInfoToSequences const & src ):
	protocols::moves::Mover( src ), target_sequences_(src.target_sequences_),
	mode_(src.mode_), throw_on_fail_(src.throw_on_fail_)
{}

/// @brief = operator
AlignPDBInfoToSequences&
AlignPDBInfoToSequences::operator=( AlignPDBInfoToSequences const & src ) {
	this->set_mode(src.get_mode());
	this->set_target_sequences(src.get_target_sequences());
	this->set_throw_on_fail(src.get_throw_on_fail());
	return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
AlignPDBInfoToSequences::~AlignPDBInfoToSequences(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////


void
AlignPDBInfoToSequences::apply_single_sequence( core::pose::Pose& pose ) {
	core::sequence::SWAligner sw_align;
	core::sequence::ScoringSchemeOP ss(  utility::pointer::make_shared< core::sequence::SimpleScoringScheme >(120, -999, -100, 0) );

	core::sequence::SequenceOP const pose_seq(
		utility::pointer::make_shared< core::sequence::Sequence >(pose.sequence(), "query", 1));

	std::string target_combo_str("");
	utility::vector1<int> linear_residue_numbers;
	utility::vector1<std::string> linear_chains;
	utility::vector1<std::string> linear_insCodes;
	utility::vector1<std::string> linear_segmentIDs;
	{
		int starting_number(1);
		for ( core::Size j=1; j<=target_sequences_.size(); ++j ) {
			target_sequences_[j].check_single_format();
			target_sequences_[j].set_1_to_1_format(starting_number);
			starting_number += target_sequences_[j].sequence_.size();
			target_combo_str += target_sequences_[j].sequence_;
			linear_chains.insert(linear_chains.end(), target_sequences_[j].chains_.begin(), target_sequences_[j].chains_.end());

			linear_insCodes.insert(linear_insCodes.end(), target_sequences_[j].insCodes_.begin(), target_sequences_[j].insCodes_.end());
			linear_segmentIDs.insert(linear_segmentIDs.end(), target_sequences_[j].segmentIDs_.begin(), target_sequences_[j].segmentIDs_.end());
			linear_residue_numbers.insert(linear_residue_numbers.end(), target_sequences_[j].residue_numbers_.begin(), target_sequences_[j].residue_numbers_.end());
		}
	}

	core::sequence::SequenceOP const target_seq(
		utility::pointer::make_shared< core::sequence::Sequence >( target_combo_str, "target", 1 )
	);


	core::sequence::SequenceAlignment const pose_seq_to_target(
		sw_align.align( pose_seq, target_seq, ss )
	);

	core::sequence::SequenceOP const pose_aln(pose_seq_to_target.sequence(1));
	core::sequence::SequenceOP const target_aln(pose_seq_to_target.sequence(2));
	pad_sequences(pose_seq->sequence(), pose_aln, target_seq->sequence(), target_aln);

	// There should be no gaps in target_aln
	// ez TODO to change this to string when chains are strings
	for ( core::Size i=1, pose_loc=1; i<=target_aln->sequence().size(); ++i ) {
		if ( pose_aln->sequence()[i-1] != '-' ) {
			pose.pdb_info()->set_resinfo(
				pose_loc,
				linear_chains[i].at(0),
				linear_residue_numbers[i],
				linear_insCodes[i].at(0),
				linear_segmentIDs[i]
			);
			++pose_loc;
		}
	}
	TR.Debug << "Pose to target alignment:" << std::endl;
	TR.Debug << *pose_aln << std::endl;
	TR.Debug << *target_aln << std::endl;
}


// utilities
utility::vector1< core::sequence::SequenceOP >
get_pose_sequences( core::pose::Pose const & pose ) {
	utility::vector1< core::sequence::SequenceOP > chain_sequences;
	for ( core::Size i=1, current_start=1; i<=pose.num_chains(); ++i ) {
		chain_sequences.push_back(
			utility::pointer::make_shared< core::sequence::Sequence >( pose.chain_sequence(i), std::to_string(i), current_start )
		);
		current_start += chain_sequences.back()->sequence().size();
	}
	return chain_sequences;
}

///@brief take something like ----HELLO----MY---BOY---
//        and make it ["HELLO", "MY", "BOY"]
utility::vector1< std::string >
split_string_based_on_inserts(std::string const & str_in) {
	utility::vector1<std::string> str_split( utility::string_split(str_in, '-' ) );
	str_split.erase(
		std::remove_if(
		str_split.begin(),
		str_split.end(),
		[](std::string const & s) { return s.size() == 0; }),
		str_split.end());
	return str_split;
}

void
AlignPDBInfoToSequences::apply_multi_sequence( core::pose::Pose& pose ) {
	// 1. collect all pose sequences
	core::sequence::SWAligner sw_align;
	core::sequence::ScoringSchemeOP ss(  utility::pointer::make_shared< core::sequence::SimpleScoringScheme >(120, -999, -100, 0) );
	utility::vector1< core::sequence::SequenceOP > const chain_sequences(get_pose_sequences(pose));
	// id == <chain, altcode, segmentID>
	std::map< std::tuple< std::string, char, std::string >, std::set<core::Size> > seen_resIDs;
	// If target sequence is 'CBAABC', and we have pose with seqs: ['ABC', 'CBA']
	// both chains will have the same chainID even though they are separate.
	// we should keep a tally of the last seen residue number and if we go behind
	// it we should increase the index.
	std::map< std::tuple< std::string, char, std::string >, core::Size > last_resnum_from_prev;

	// 2. Try aligning all pose sequences to given sequences
	for ( core::Size i=1; i<=chain_sequences.size(); ++i ) {
		TR.Debug << "current sequence: " << *chain_sequences[i] << std::endl;
		core::sequence::SequenceOP const & current_pose_seq( chain_sequences[i] );

		core::Size correct_seq_idx(0);
		core::sequence::SequenceOP current_aln(nullptr), target_aln(nullptr);
		utility::vector1<core::sequence::SequenceAlignment> seen_alignments;

		// current_pass 0 == perfect matches, current_pass 1 == matches with 1 jump, etc
		// try all sequences -- if no alignment found try all sequences again with +1 gap allowed
		for ( core::Size j=1, current_pass=0; j<=target_sequences_.size() && correct_seq_idx == 0 && current_pass <= sequence_alignment_cut_max_; ++j ) {
			core::sequence::SequenceOP const current_target_seq(
				utility::pointer::make_shared< core::sequence::Sequence >( target_sequences_[j].sequence_, std::to_string(1000+j), 1 )
			);

			// Since sw_align can throw an error if you're making a really bad alignment
			// (which this function is basically doing -- ie trying alignments until they work)
			// we can use a lambda to properly init chainSeq_to_target and keep track if
			// it failed or not inside a try/catch.
			bool failed_alignment(false);
			core::sequence::SequenceAlignment const chainSeq_to_target = [&](){
				try {
					return sw_align.align( current_pose_seq, current_target_seq, ss );
				} catch ( utility::excn::NullPointerError const & ) {
					TR.Debug << "caught thrown align excn: " << *current_pose_seq  << " : " << *current_target_seq << std::endl;
					failed_alignment = true;
					return core::sequence::SequenceAlignment();
				}
			}();
			if ( failed_alignment ) {
				if ( j == target_sequences_.size() ) {
					j=0;  // will get incremented to 1 on loop end
					++current_pass;
				}
				continue;
			}

			seen_alignments.push_back(chainSeq_to_target);
			current_aln = chainSeq_to_target.sequence(1);
			// Start must be 1 for pad sequences to work.
			current_aln->start(1);
			target_aln = chainSeq_to_target.sequence(2);
			utility::vector1< std::string > const split_and_trimmed_pose_aln(
				split_string_based_on_inserts(current_aln->sequence()));

			// check that target sequence has NO gaps AND all the current chain is included in the alignment AND that pose seq has NO gaps
			if ( (target_aln->sequence().find("-") == std::string::npos)
					&& (current_aln->sequence().size() >= current_pose_seq->sequence().size() )
					&& (current_pass == split_and_trimmed_pose_aln.size() ) ) {
				correct_seq_idx = j;
				break; // redundant with loop check but easy to follow
			}

			if ( j == target_sequences_.size() ) {
				j=0;  // will get incremented to 1 on loop end
				++current_pass;
			}
		}

		if ( correct_seq_idx == 0 ) {
			std::stringstream ss;
			ss << "Failed to align sequence:" << *current_pose_seq << std::endl;
			if ( throw_on_fail_ ) {
				throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
			} else {
				TR.Error << ss.str() << std::endl;
				continue;
			}
		}

		TR.Debug << "pre-padded alignments (current/target)" << std::endl;
		TR.Debug << *current_aln << std::endl;
		TR.Debug << *target_aln << std::endl;
		pad_sequences(current_pose_seq->sequence(), current_aln, target_sequences_[correct_seq_idx].sequence_, target_aln);
		TR.Debug << "padded alignments (current/target)" << std::endl;
		TR.Debug << *current_aln << std::endl;
		TR.Debug << *target_aln << std::endl;
		TR.Debug << "current idx: " << target_sequences_[correct_seq_idx].current_idx_ << std::endl;

		// 1. determine if any of the current residues have already been placed with the
		// current chainID
		bool should_increment_current_idx(false);
		{
			// 1.1 check for overflow of current_idx_
			SequenceSpecification const & current_seqSpec(target_sequences_[correct_seq_idx]);
			core::Size const & current_idx(target_sequences_[correct_seq_idx].current_idx_);
			if ( current_idx > current_seqSpec.chains_.size() ) {
				std::stringstream ss;
				ss << "Too many of chain (precheck): " << current_seqSpec.sequence_ << " current idx: " << current_idx
					<< " chains size: " << current_seqSpec.chains_.size() << std::endl;
				if ( throw_on_fail_ ) {
					throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
				} else {
					TR.Error << ss.str() << std::endl;
					continue;
				}
			}
			// 1.2 check for seen resIDs
			std::string const current_chain(current_seqSpec.chains_[current_idx]);
			char const & insCode(current_seqSpec.insCodes_.size() >= current_idx ? current_seqSpec.insCodes_[current_idx].at(0) : ' ');
			std::string const & seg_id(current_seqSpec.segmentIDs_.size() >= current_idx ? current_seqSpec.segmentIDs_[current_idx] : "    ");
			std::tuple< std::string, char, std::string > const current_id(std::tie(current_chain, insCode, seg_id));
			for ( core::Size j=0; j<current_aln->length(); ++j ) {
				if ( current_aln->sequence()[j] == '-' ) continue;
				if ( seen_resIDs.count(current_id) == 0 ) {
					seen_resIDs[current_id]; // init empty
					break;
				}
				if ( seen_resIDs[current_id].count(j+1) ||
						(last_resnum_from_prev.count(current_id) && last_resnum_from_prev[current_id] > j+1 ) ) {
					TR.Debug << "incrementing seq: " << current_seqSpec.sequence_ << " from: " << target_sequences_[correct_seq_idx].current_idx_
						<< " seen: " << seen_resIDs[current_id].count(j+1) << " or prev: " << last_resnum_from_prev.count(current_id)
						<< " last resnum: " << last_resnum_from_prev[current_id] << " vs: " << j+1 << " bool: "
						<< (last_resnum_from_prev.count(current_id) && last_resnum_from_prev[current_id] < j+1 )  << std::endl;
					should_increment_current_idx = true;
					break;
				}
			}
		}
		// 2. set the new chainIDs/residue_numbers based on current_idx_ (increment if has been found)
		if ( should_increment_current_idx ) ++target_sequences_[correct_seq_idx].current_idx_;
		{
			core::Size const & current_idx(target_sequences_[correct_seq_idx].current_idx_);
			TR.Debug << "post current idx: " << current_idx << std::endl;
			SequenceSpecification const & current_seqSpec(target_sequences_[correct_seq_idx]);
			// 2.1 check for overflow of current_idx_
			if ( current_idx > current_seqSpec.chains_.size() ) {
				std::stringstream ss;
				ss << "Too many of chain: " << current_seqSpec << " current idx: " << current_idx
					<< " chains size: " << current_seqSpec.chains_.size() << " you must either supply more chainIDs"
					<< "or check your protocol's output";
				if ( throw_on_fail_ ) {
					throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
				} else {
					TR.Error << ss.str() << std::endl;
					continue;
				}
			}

			std::string const current_chain(current_seqSpec.chains_[current_idx]);
			char const & insCode(current_seqSpec.insCodes_.size() >= current_idx ? current_seqSpec.insCodes_[current_idx].at(0) : ' ');
			std::string const & seg_id(current_seqSpec.segmentIDs_.size() >= current_idx ? current_seqSpec.segmentIDs_[current_idx] : "    ");
			std::tuple< std::string, char, std::string > const current_id(std::tie(current_chain, insCode, seg_id));
			utility::vector1<int> const & residue_numbers(current_seqSpec.residue_numbers_);
			for ( core::Size j=0, found_count=0; j<current_aln->length(); ++j ) {
				if ( current_aln->sequence()[j] == '-' ) continue;
				if ( seen_resIDs.count(current_id) == 0 ) {
					seen_resIDs[current_id]; // init empty
				}
				seen_resIDs[current_id].insert(j+1);
				last_resnum_from_prev[current_id] = j+1;

				pose.pdb_info()->set_resinfo(
					current_pose_seq->start() + found_count,
					current_chain.at(0),
					residue_numbers[j+1],
					insCode,
					seg_id
				);
				++found_count;
			}
		}
	}
}

/// @brief Apply the mover
void
AlignPDBInfoToSequences::apply( core::pose::Pose& pose ) {
	for ( auto & seq : this->get_target_sequences() ) seq.current_idx_ = 1;

	if ( !pose.pdb_info() ) {
		pose.pdb_info(utility::pointer::make_shared< core::pose::PDBInfo >( pose ));
	}

	switch (this->get_mode()) {
	case AlignPDBInfoToSequencesMode::single :
		for ( auto & seq : this->get_target_sequences() ) seq.check_single_format();
		apply_single_sequence(pose);
		break;
	case AlignPDBInfoToSequencesMode::multiple :
		for ( auto & seq : this->get_target_sequences() ) seq.update_multiple_format_residue_numbering();
		apply_multi_sequence(pose);
		break;
		// This will happen if mode_ is unset or if someone adds another mode and forgets to update this
	default :
		std::stringstream ss;
		ss << "AlignPDBInfoToSequences Did not run! mode was not 'single' or 'multiple'";
		if ( throw_on_fail_ ) throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
		TR.Error << ss.str() << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
AlignPDBInfoToSequences::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

void
AlignPDBInfoToSequences::load_json_file( std::string const & json_fn ) {
	utility::io::izstream instream(json_fn);
	nlohmann::json json_result;
	instream >> json_result;
	utility::vector1<SequenceSpecification> const tmp_in(json_result.get<utility::vector1<SequenceSpecification>>());
	this->get_target_sequences().insert(this->get_target_sequences().end(),  tmp_in.begin(), tmp_in.end());
}

SequenceSpecification
AlignPDBInfoToSequences::parse_target_tag( utility::tag::TagCOP const & tag ) const {
	std::string const sequence(tag->getOption< std::string >("sequence"));
	std::string const chains_raw(tag->getOption< std::string >("chains"));
	std::string const segmentIDs_raw(tag->getOption< std::string >("segmentIDs", ""));
	std::string const insCodes_raw(tag->getOption< std::string >("insCodes", ""));
	std::string const residue_numbers_raw(tag->getOption< std::string >("residue_numbers", ""));
	auto set_str_vec = [](std::string const & str) {
		utility::vector1<std::string> ret;
		if ( !str.empty() ) ret = utility::string_split(str, ',');
		return ret;
	};
	auto set_int_vec = [](std::string const & str) {
		nlohmann::json json = nlohmann::json::parse("[" + str + "]");
		utility::vector1<int> const ret(json.get<utility::vector1<int>>());
		return ret;
	};
	utility::vector1<std::string> const chains(set_str_vec(chains_raw));
	utility::vector1<std::string> const segmentIDs(set_str_vec(segmentIDs_raw));
	utility::vector1<std::string> const insCodes(set_str_vec(insCodes_raw));
	utility::vector1<int> const residue_numbers(set_int_vec(residue_numbers_raw));
	return SequenceSpecification(sequence, chains, segmentIDs, insCodes, residue_numbers);
}


/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
AlignPDBInfoToSequences::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&)
{
	// 1. load sequences
	if ( tag->hasOption("json_fns") ) {  // 1.a load jsons if available
		std::string const & json_fns = tag->getOption<std::string>("json_fns");
		for ( std::string const & json_fn : utility::string_split(json_fns, ',') ) {
			this->load_json_file(json_fn);
		}
	}
	{  // 1.b load from tags
		utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
		for ( auto tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
			if ( (*tag_it)->getName() == "Target" ) {
				this->get_target_sequences().push_back(parse_target_tag(*tag_it));
			}
		}
	}
	this->set_mode(get_alignpdbinfotosequencesmode_from_string(tag->getOption<std::string>("mode")));

	// ok to throw early!
	if ( this->get_mode() == AlignPDBInfoToSequencesMode::single ) {
		for ( auto seq : this->get_target_sequences() ) seq.check_single_format();
	} else if ( this->get_mode() == AlignPDBInfoToSequencesMode::multiple ) {
		for ( auto seq : this->get_target_sequences() ) seq.update_multiple_format_residue_numbering();
	}

	if ( tag->hasOption( "throw_on_fail" ) ) {
		this->set_throw_on_fail(tag->getOption<bool>("throw_on_fail"));
	}

	if ( tag->hasOption( "sequence_alignment_cut_max" ) ) {
		this->set_sequence_alignment_cut_max(tag->getOption<core::Size>("sequence_alignment_cut_max"));
	}

}

std::string
AlignPDBInfoToSequences_namer( std::string const & name ) {
	return "AlignPDBInfoToSequences_subelement_" + name + "_type";
}

void AlignPDBInfoToSequences::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "mode", xs_string, "Which mode to run in. options: ['single', 'multiple']")
		+ XMLSchemaAttribute::attribute_w_default( "json_fns", xs_string, "The name of the json sequence file(s) (separated by ',')", "")
		+ XMLSchemaAttribute::attribute_w_default( "throw_on_fail", xsct_rosetta_bool, "throw on failure.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "sequence_alignment_cut_max", xsct_positive_integer, "Maximum number of cuts to allow when doing the sequence alignment.", "10" );

	AttributeList target_subelement_attributes;
	target_subelement_attributes
		+ XMLSchemaAttribute::attribute_w_default( "name", xs_string, "unused but useful for keeping track of what protein you're looking at", "")
		+ XMLSchemaAttribute::required_attribute( "sequence", xs_string, "sequence of current chain/protein")
		+ XMLSchemaAttribute::required_attribute( "chains", xs_string, "chains to set with current protein")
		+ XMLSchemaAttribute::attribute_w_default( "segmentIDs", xs_string, "segmentIDs to set for the current protein", "")
		+ XMLSchemaAttribute::attribute_w_default( "insCodes", xs_string, "insertion codes to set for the current protein", "")
		+ XMLSchemaAttribute::attribute_w_default( "residue_numbers", xs_string, "Residue numbers to set for current chain (must be of size 1 (used as a starting number), or the size of 'sequence', or empty)", "");

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & AlignPDBInfoToSequences_namer );
	subelements.add_simple_subelement( "Target", target_subelement_attributes, "Targets and settings of the sequences to align to.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Align the PDBInfo of a pose to some given sequences. This mover does not alter"
		" geometry or sequence of the pose. The only thing that is altered is the PDBInfo."
		" The goal of this mover is to re-number/re-chain a pose so that its PDBInfo is in"
		" sync with some reference.\n"
		" There are two modes. In 'single' mode all the target sequences"
		" are appended in order and the entire pose is aligned to that sequence."
		" In 'multiple' mode we attempt to align the sequence of each chain of the pose"
		" to each of the target sequences.  If we experience a perfect match, we will"
		" renumber based on that target sequence.\n"
		" For more information, please reference the online documentation.\n"
		" WARNING 1: if you have a very long sequence, (1000s of residues) this protocol"
		" can use a lot of memory due to the SmithWaterman alignment.  Be cautious of this.\n"
		" WARNING 2: This will most likely not give you the correct results if your pose"
		" or target sequence contain carbohydrate or other residues that are represented"
		" using [ ] or ( ).",
		attlist, subelements );

	//basic direct attributes

}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AlignPDBInfoToSequences::fresh_instance() const
{
	return utility::pointer::make_shared< AlignPDBInfoToSequences >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AlignPDBInfoToSequences::clone() const
{
	return utility::pointer::make_shared< AlignPDBInfoToSequences >( *this );
}

std::string AlignPDBInfoToSequences::get_name() const {
	return mover_name();
}

std::string AlignPDBInfoToSequences::mover_name() {
	return "AlignPDBInfoToSequences";
}


void
AlignPDBInfoToSequences::set_mode_from_string(std::string const & mode) {
	this->set_mode(get_alignpdbinfotosequencesmode_from_string(mode));
}

std::string
AlignPDBInfoToSequences::get_mode_as_string() const {
	return alignpdbinfotosequencesmode_to_string(this->get_mode());
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
AlignPDBInfoToSequencesCreator::create_mover() const
{
	return utility::pointer::make_shared< AlignPDBInfoToSequences >();
}

std::string
AlignPDBInfoToSequencesCreator::keyname() const
{
	return AlignPDBInfoToSequences::mover_name();
}

void AlignPDBInfoToSequencesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignPDBInfoToSequences::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, AlignPDBInfoToSequences const & mover )
{
	mover.show(os);
	return os;
}

} //pdbinfo_manipulations
} //protocols
