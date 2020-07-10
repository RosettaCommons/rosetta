// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pdbinfo_manipulations/AlignPDBInfoToSequences.hh
/// @brief Realign poses to sequences after losing pdb_info
/// @author Dan Farrell (danpf@uw.edu)

#ifndef INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequences_HH
#define INCLUDED_protocols_pdbinfo_manipulations_AlignPDBInfoToSequences_HH

// Unit headers
#include <protocols/pdbinfo_manipulations/AlignPDBInfoToSequences.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/json_utilities.hh>


namespace protocols {
namespace pdbinfo_manipulations {

enum class AlignPDBInfoToSequencesMode {
	unset,
	single,
	multiple
};

AlignPDBInfoToSequencesMode
get_alignpdbinfotosequencesmode_from_string(std::string const & str);

std::string
alignpdbinfotosequencesmode_to_string(AlignPDBInfoToSequencesMode const & mode);

///@brief Rosetta defines the pdbinfo of an object via possible chain IDs, segment IDs,
///       and insertion codes.
struct SequenceSpecification {

	SequenceSpecification() : current_idx_(1) {;}

	SequenceSpecification(
		std::string const & sequence,
		utility::vector1<std::string> const & chains,
		utility::vector1<std::string> const & segmentIDs = utility::vector1<std::string>(),
		utility::vector1<std::string> const & insCodes = utility::vector1<std::string>(),
		utility::vector1<int> const & residue_numbers = utility::vector1<int>()) :
		sequence_(sequence), chains_(chains), segmentIDs_(segmentIDs), insCodes_(insCodes), residue_numbers_(residue_numbers), current_idx_(1) {;}

	///@brief make sure none of the non-sequence elements are larger than
	///       the sequence_.
	void
	check_single_format() const;

	///@brief make sure residue numbering is set to be the same size as the sequence_
	//        or set the residue_numbering to start from 1 -> sequence.size()
	void
	update_multiple_format_residue_numbering();

	///@brief duplicate the last element (or default last element if not set) to
	///       pad chains_, segmentIDs_, insCodes_, residue_numbers_,  to be the same length as
	///       the sequence_.
	void
	set_1_to_1_format(int const starting_number);

	///@brief This ignores current_idx_ because that is just a way to track how many times
	///       This has been used.
	bool
	operator==(SequenceSpecification const & alt_ss) const;

	std::string sequence_;
	utility::vector1<std::string> chains_;
	utility::vector1<std::string> segmentIDs_;
	utility::vector1<std::string> insCodes_;
	utility::vector1<int> residue_numbers_;
	core::Size current_idx_;
};

std::ostream &
operator<<( std::ostream & os, SequenceSpecification const & ss );

///@brief serialization (json) of SequenceSpecification
void to_json(nlohmann::json& j, SequenceSpecification const & ss);

void from_json(const json& j, SequenceSpecification & ss);



std::string
AlignPDBInfoToSequences_namer( std::string const & name );


///@brief Realign poses to sequences after losing pdb_info
class AlignPDBInfoToSequences : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	AlignPDBInfoToSequences();

	/// @brief Copy constructor (not needed unless you need deep copies)
	AlignPDBInfoToSequences( AlignPDBInfoToSequences const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~AlignPDBInfoToSequences() override;

	/// @brief = operator
	AlignPDBInfoToSequences& operator=( AlignPDBInfoToSequences const & src );

public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	void
	set_target_sequences(utility::vector1<SequenceSpecification> const & target_sequences) {
		target_sequences_ = target_sequences;
	}

	utility::vector1<SequenceSpecification> const &
	get_target_sequences() const { return target_sequences_; }

	utility::vector1<SequenceSpecification>&
	get_target_sequences() { return target_sequences_; }

	void
	set_mode(AlignPDBInfoToSequencesMode const & mode) {
		mode_ = mode;
	}

	void
	set_mode_from_string(std::string const & mode);

	AlignPDBInfoToSequencesMode const &
	get_mode() const { return mode_; }

	std::string
	get_mode_as_string() const;

	void
	set_throw_on_fail(bool const throw_on_fail) {
		throw_on_fail_ = throw_on_fail;
	}

	bool
	get_throw_on_fail() const { return throw_on_fail_; }

	void
	set_sequence_alignment_cut_max(core::Size const sequence_alignment_cut_max) {
		sequence_alignment_cut_max_ = sequence_alignment_cut_max;
	}

	core::Size
	get_sequence_alignment_cut_max() const { return sequence_alignment_cut_max_; }

public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////
	SequenceSpecification
	parse_target_tag(utility::tag::TagCOP const & tag) const;

	void
	load_json_file( std::string const & json_fn );

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // methods

	///@brief This method tries to align every sequence of the pose to every target sequence
	///       by applying a SWAligner to them, and if we find a perfect alignment we use that
	///       target sequence as the reference.
	void apply_multi_sequence( core::pose::Pose& pose );

	///@brief This method tries to align the pose to a target sequence that is generated by
	///       appending all target_sequences_ together. This allows for more residue-level
	///       control of the residue numbers/chainIDs but uses more memory (larger SWAligner).
	///       Most of the time you should probably use apply_multi_sequence, but sometimes
	///       you just crave that control!!!
	void apply_single_sequence( core::pose::Pose& pose );

private: // data
	utility::vector1<SequenceSpecification> target_sequences_;
	AlignPDBInfoToSequencesMode mode_;
	bool throw_on_fail_ = false;
	core::Size sequence_alignment_cut_max_ = 10;
};

std::ostream &
operator<<( std::ostream & os, AlignPDBInfoToSequences const & mover );

} //pdbinfo_manipulations
} //protocols

#endif //protocols_pdbinfo_manipulations_AlignPDBInfoToSequences_HH
