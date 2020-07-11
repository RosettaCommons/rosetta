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

#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/DEERIO.fwd.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>
#include <utility/fixedsizearray1.hh>
#include <iosfwd>
#include <string>

namespace core {
namespace scoring {
namespace epr_deer {

class DEERIO {

public:

	/// @brief    Constructor
	DEERIO();

	/// @brief    Destructor
	~DEERIO();

	/// @brief Generates the data for use. Core function
	std::map< Size, DEERDataOP >
	generate_data(
		core::pose::Pose const & pose
	);

	/// @brief Reads the input file(s) and makes an unsorted vector of their whitespace-separated contents
	utility::vector1< utility::vector1< std::string > >
	get_splitlines();

	/// @brief Reads lines that start with DESC, which contain residue and spin label information
	void
	read_desc_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & splitlines,
		pose::Pose const & pose
	);

	/// @brief Reads lines that start with INFO, which gives additional options specific for each bit of data
	void
	read_info_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & splitlines
	);

	/// @brief Read lines for decay data (raw data)
	void
	read_decay_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & data_lines
	);

	/// @brief Read lines for bounded distance restraints
	void
	read_bounds_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & splitlines
	);

	/// @brief Read lines for gaussian-distributed distance data
	void
	read_gauss_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & splitlines
	);

	/// @brief Read lines for non-gaussian-distributed distance data
	void
	read_dist_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & splitlines
	);

	/// @brief Read the lines starting with PAIR, which are for easy input
	void
	read_pair_lines(
		std::map< Size, DEERDataOP > & output,
		utility::vector1< utility::vector1< std::string > > const & splitlines,
		pose::Pose const & pose
	);

	/// @brief Normalize distance distributions so their integrals equal 1
	Real
	normalize_distribution(
		std::map< Size, Real > & in_map
	);

	/// @brief Get custom spin label coordinates
	utility::vector1< std::pair< EPRSpinLabel, Real > >
	pull_coords() const;

private:

	Size ANGSTROM_LIMIT_   = 100;
	Size BINS_PER_ANGSTROM_ = 2;

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
