// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceMetric.hh
/// @brief A SimpleMetric to output the single-letter OR three-letter sequence of a protein or subset of positions/regions using a ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added support for writing full residue type names or basenames.

#ifndef INCLUDED_core_simple_metrics_metrics_SequenceMetric_HH
#define INCLUDED_core_simple_metrics_metrics_SequenceMetric_HH

#include <core/simple_metrics/metrics/SequenceMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief The mode for this metric.  If you add to this list, be sure to update the map associating this enum
/// with corresponding strings in the .cc file.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
enum SequenceMetricMode {
	SMM_ONELETTER_CODE=1,
	SMM_THREELETTER_CODE,
	SMM_BASE_NAME,
	SMM_FULL_NAME,
	SMM_INVALID_MODE, //Keep this second-to-last.
	SMM_END_OF_LIST = SMM_INVALID_MODE //Keep this last
};

///@brief A SimpleMetric to output the single-letter OR three-letter sequence of a protein or subset of positions/regions using a ResidueSelector.
class SequenceMetric : public core::simple_metrics::StringMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SequenceMetric();

	SequenceMetric( select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Copy constructor (not needed unless you need deep copies)
	SequenceMetric( SequenceMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SequenceMetric() override;

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in RealMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate the metric.
	std::string
	calculate( core::pose::Pose const & pose ) const override;

public:

	void
	set_residue_selector( select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Set the output mode -- one-letter code (e.g. Y), three-letter code (e.g. DTY), residue base name (e.g. DTYR), or
	/// full residue name (e.g. DTYR:CtermProteinFull).
	/// @details Throws an error if invalid enum provided.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	void set_output_mode( SequenceMetricMode const mode_in );

	/// @brief Set the output mode using the string corresponding to the output mode.
	/// @details Throws an error if string is invalid.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	void set_output_mode( std::string const & mode_in );

public:

	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	///@brief Name of the metric
	std::string
	metric() const override;

	/// @brief Returns all allowed output modes.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	static std::string allowed_output_modes();

	/// @brief Returns all allowed output modes as a vector of strings.
	/// @details Note that this is a bit inefficient.  It generates the vector each time, and returns it by copy.
	/// Not intended for repeated calls.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	static utility::vector1< std::string > allowed_output_modes_as_vector();

	/// @brief Given an output mode enum, get its string representation.
	/// @details Returns "INVALID" if invalid.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	static std::string const & mode_name_from_enum( SequenceMetricMode const mode_enum );


	/// @brief Given an output mode string, get its enum.
	/// @details Returns SMM_INVALID_MODE if invalid.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	static SequenceMetricMode mode_enum_from_name( std::string const & mode_string );

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

private:

	select::residue_selector::ResidueSelectorCOP selector_ = nullptr;

	/// @brief The output mode -- one-letter code (e.g. Y), three-letter code (e.g. DTY), residue base name (e.g. DTYR), or
	/// full residue name (e.g. DTYR:CtermProteinFull).
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	SequenceMetricMode output_mode_;

};

} //core
} //simple_metrics
} //metrics



#endif //core_simple_metrics_metrics_SequenceMetric_HH





