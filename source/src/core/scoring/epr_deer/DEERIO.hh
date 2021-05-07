// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/DEERIO.hh
/// @brief  IO class for data obtained with double electron-electron resonance (DEER)
/// @details This class is only called once at the beginning of the application by the
///      energy method. It then parses the input file and the electron coordinate file.
///      The input file is then used to create a "map" (essentially a fake graph) which
///      then becomes the basis of the DEERData object.

/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_DEERIO_hh
#define INCLUDED_core_scoring_epr_deer_DEERIO_hh

#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/DEERIO.fwd.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzVector.hh>
#include <utility/fixedsizearray1.hh>
#include <iosfwd>
#include <string>

namespace core {
namespace scoring {
namespace epr_deer {

static std::map< Size, metrics::DEERDataOP > epr_deer_input_data_;

using SplitLine = utility::vector1< std::string >;
using SplitLines = utility::vector1< SplitLine >;
using SplitLineMap = std::map< std::string, SplitLines >;

/// @brief Alias to asve space when defining residues and spin labels
using PairSizeString = std::pair< Size, std::string >;

class DEERIO {

public:

	/// @brief    Constructor
	DEERIO() = default;

	/// @brief    Copy Constructor
	DEERIO(
		DEERIO const & io
	) = default;

	/// @brief    Destructor
	~DEERIO() = default;

	/// @brief Split an input line into vector of strings
	/// @param  line: Line of interest
	/// @return Vector of strings
	/// @detail Ignores lines starting with "#"
	SplitLine
	split_line(
		std::string const & line
	) const;

	/// @brief Generates the data for use. Core function
	/// @param  pose: Input pose (used for residue numbering)
	/// @return Map of input data
	/// @detail Will only compute this once! Otherwise, the data is saved
	/// @detail  in a static/global object (epr_deer_input_data_) defined in
	/// @detail  header file
	std::map< Size, metrics::DEERDataOP >
	generate_data(
		core::pose::Pose const & pose
	) const;

	/// @brief  Parses file
	/// @param  pose: Pose (for residue numbering)
	/// @return Map of strings (grouped by data index)
	std::map< Size, SplitLines >
	parse_input() const;

	/// @brief Parse the input lines into actual data containers
	/// @param  pose: Pose (used for residue numbering)
	/// @param  idx: Index of datum (used for error messages)
	/// @param  lines: Lines associated with this dataset
	/// @return Proper data container that can be used by Rosetta
	metrics::DEERDataOP
	parse_data(
		pose::Pose const & pose,
		SplitLines const & lines
	) const;

	/// @brief  Parse inputs lines into DEERDecayDataOP
	/// @param  pose: Pose (used for residue renumbering)
	/// @param  sorted_data: Input lines sorted by linetype
	/// @return Proper data container that can be used by Rosetta
	metrics::DEERDataOP
	parse_decay_data(
		pose::Pose const & pose,
		std::map< std::string, SplitLines > const & sorted_data
	) const;

	/// @brief  Parse inputs lines into DEERDistanceBoundsOP
	/// @param  pose: Pose (used for residue renumbering)
	/// @param  sorted_data: Input lines sorted by linetype
	/// @return Proper data container that can be used by Rosetta
	metrics::DEERDataOP
	parse_bounds_data(
		pose::Pose const & pose,
		std::map< std::string, SplitLines > const & sorted_data
	) const;

	/// @brief  Parse inputs lines into DEERDistanceDistributionOP
	/// @param  pose: Pose (used for residue renumbering)
	/// @param  sorted_data: Input lines sorted by linetype
	/// @return Proper data container that can be used by Rosetta
	metrics::DEERDataOP
	parse_dist_data(
		pose::Pose const & pose,
		std::map< std::string, SplitLines > const & sorted_data
	) const;

	/// @brief Parses lines labeled "INFO" into information usable by Rosetta
	/// @param  info_lines: Lines with the information (as strings)
	/// @return Proper data container that can be used by Rosetta
	/// @detail Note that this container, upon output, is incomplete!
	metrics::DEERDataOP
	parse_dist_info_lines(
		SplitLines const & info_lines
	) const;

	/// @brief Parses lines labeled "DIST" into information usable by Rosetta
	/// @param  dist_lines: Lines with the distance distributions
	/// @param  data: DEERData object (passed by reference and modified here)
	void
	parse_dist_lines(
		SplitLines const & dist_lines,
		metrics::DEERDataOP data
	) const;

	/// @brief  Determines function used by DEERDistanceDistribution
	/// @param  info_lines: Lines labeled "INFO" in input file
	/// @return Proper data container that can be used by Rosetta
	/// @detail Note that this container, upon output, is incomplete!
	metrics::DEERDataOP
	parse_dist_datatype(
		SplitLines const & info_lines
	) const;

	/// @brief  Parses data for DEERDistanceBounds object
	/// @param  bounds_lines: Lines labeled "BOUNDS" in input file
	/// @param  data: Proper data container that can be used by Rosetta
	/// @detail Note that only the last BOUNDS line will be used
	void
	parse_bounds_lines(
		SplitLines const & bounds_lines,
		metrics::DEERDataOP bounds_data
	) const;

	/// @brief  Parses information data for DEERDecayData object
	/// @param  info_lines: Lines labeled "INFO" in input file
	/// @param  data: Proper data container that can be used by Rosetta
	void
	parse_decay_info_lines(
		SplitLines const & info_lines,
		metrics::DEERDataOP data
	) const;

	/// @brief  Parses decay data for DEERDecayData object
	/// @param  decay_lines: Lines labeled "DECAY" in input file
	/// @param  data: Proper data container that can be used by Rosetta
	void
	parse_decay_lines(
		SplitLines const & decay_lines,
		metrics::DEERDataOP & data
	) const;

	/// @brief  Parse input lines for information on residues involved in data
	/// @param  pose: Pose (used for residue renumbering)
	/// @param  desc_lines: Input lines containing this information
	/// @return List of residues describing both the ID# and spin label type
	utility::vector1< PairSizeString >
	parse_desc_lines(
		pose::Pose const & pose,
		SplitLines const & desc_lines
	) const;

	/// @brief  Parse coordinate files for custom spin labels
	/// @param  pose: Pose (used for residue renumbering)
	/// @return List of spin labels and corresponding weights
	utility::vector1< std::pair< EPRSpinLabel, Real > >
	pull_coords(
		pose::Pose const & pose
	) const;

private:

	/// @brief Defines the default max distance (for DEERDecayData objects)
	Size const ANGSTROM_LIMIT_   = 100;

	/// @brief Defines the default bins per angstrom in distributions
	Size const BINS_PER_A_ = 2;

	/// @brief Sets the first line that will replace "PAIR" lines
	utility::vector1< std::string > const pairline1_ =
	{ "DESC", "", "DEFAULT", "", "DEFAULT", "" };

	/// @brief Sets the second line that will replace "PAIR" lines
	utility::vector1< std::string > const pairline2_ =
	{ "BOUNDS", "", "", "", "1.0" };

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
