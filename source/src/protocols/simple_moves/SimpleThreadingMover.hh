// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SimpleThreadingMover.cc
/// @brief Very Simple class for threading a regional sequence onto a structure
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @modified NCAA support added by Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_simple_moves_SimpleThreadingMover_hh
#define INCLUDED_protocols_simple_moves_SimpleThreadingMover_hh

#include <protocols/simple_moves/SimpleThreadingMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/simple_metrics/metrics/SequenceMetric.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <basic/datacache/DataMap.fwd.hh>

#include <map>

// Forward
namespace protocols {
namespace simple_moves {


///@brief
/// This mover functions to thread the sequence of a region onto the given pose.  Nothing fancy here.
/// For more fancy things see protocols/comparative_modeling.
///
///@details
/// It does the threading by allowing the task to only enable these residues and then does a repacking. Optionally repack neighbors so we save one more step.
/// A sequence is just a string, additional '-' charactors denote to skip this position in the thread.
/// Default is 5 rounds of packing.
///
class SimpleThreadingMover : public protocols::moves::Mover {

public :

	/// @brief Default constructor.
	SimpleThreadingMover() = default;

	/// @brief Initialization constructor.
	SimpleThreadingMover(std::string thread_sequence, core::Size start_position);

	/// @brief Default copy constructor.
	SimpleThreadingMover(SimpleThreadingMover const & ) = default;

	/// @brief Destructor.
	~SimpleThreadingMover() override;


	///@brief Set the sequence to thread onto the structure used in apply and where to start.
	///@details Can have '-' charactors in sequence to denote a gap in the threaded sequence.
	void
	set_sequence(std::string thread_sequence, core::Size start_position);

	///@brief Pack the neighbor residues?
	void
	set_pack_neighbors(bool pack_neighbors);

	///@brief Set the packing distance for neighbor pack.
	void
	set_neighbor_distance(core::Real neighbor_dis);

	/// @brief Set the sequence mode, by string.
	/// @details This determines how the sequence is interpreted (one-letter codes, three-letter codes, etc.).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	set_sequence_mode ( std::string const & mode_string_in );

	/// @brief Set the sequence mode, by enum.
	/// @details This determines how the sequence is interpreted (one-letter codes, three-letter codes, etc.).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	set_sequence_mode ( core::simple_metrics::metrics::SequenceMetricMode const mode_in );

	///////////////

	bool
	get_pack_neighbors() const;

	core::Real
	get_neighbor_distance() const;

	/// @brief Get the sequence mode.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	inline core::simple_metrics::metrics::SequenceMetricMode sequence_mode() const { return sequence_mode_; }

	///@brief Set the scorefunction used for packing.
	void
	set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn);

	///@brief Set the number of pack rounds.
	void
	set_pack_rounds(core::Size pack_rounds);

	///@brief Set whether we skip unknown residue types, or throw errors.
	inline void set_skip_unknown_mutant( bool const setting ) { skip_unknown_mutant_ = setting; }

	///@brief Get whether we skip unknown residue types, or throw errors.
	inline bool skip_unknown_mutant() const { return skip_unknown_mutant_; }

	void
	apply(core::pose::Pose & pose) override;


public:


	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	protocols::moves::MoverOP
	clone() const override;

	//SimpleThreadingMover & operator=( SimpleThreadingMover const & src);

	moves::MoverOP fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public:  // Citation Management

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private: //Functions

	/// @brief Given the sequence, the interpretation mode, and the start position, fill a map of position->mutation name.
	std::map< core::Size, std::string >
	determine_mutations( core::Size const start_position, std::string const & sequence, core::simple_metrics::metrics::SequenceMetricMode const mode, core::pose::Pose const & pose ) const;

	/// @brief Given the sequence as one-letter codes and the start position, fill a map of position->mutation name.
	/// @details Supports possibility of positions of the form X[NCAA_NAME].
	std::map< core::Size, std::string >
	determine_mutations_oneletter( core::Size const start_position, std::string const & sequence ) const;

	/// @brief Given the sequence as a comma-separated list of either three-letter codes, base names, or full names, plus the start position,
	/// fill a map of position->mutation name.
	std::map< core::Size, std::string >
	determine_mutations_comma_separated( core::Size const start_position, std::string const & sequence, core::simple_metrics::metrics::SequenceMetricMode const mode, core::pose::Pose const & pose ) const;

private: //Data

	core::Size start_position_ = 1;
	std::string thread_sequence_ = "0"; //Note that "0" is the signal that the sequence has not been set.  (Changed from "NA" because ASN-ALA is a valid sequence.)

	bool pack_neighbors_ = false;
	core::Real neighbor_dis_ = 6.0;

	core::scoring::ScoreFunctionCOP scorefxn_ = nullptr;

	std::string parsed_position_ = "NA"; //Enables pose-length changes after construction of mover.  Note that "NA" is the signal that this has not been set.
	bool skip_unknown_mutant_ = false;

	core::Size pack_rounds_ = 5;

	/// @brief How should the sequence be interpreted?  Default is as a string of one-letter codes.
	core::simple_metrics::metrics::SequenceMetricMode sequence_mode_ = core::simple_metrics::metrics::SequenceMetricMode::ONELETTER_CODE;

};

}//simple_moves
}//protocols



#endif //INCLUDED_protocols_simple_moves_SimpleThreadingMover_hh







