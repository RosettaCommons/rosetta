// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/util.hh
/// @brief   Util functions for Input and Output of NMR data.
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_util_hh
#define INCLUDED_core_io_nmr_util_hh

// Unit headers
#include <core/io/nmr/AtomSelection.fwd.hh>

// Project headers
#include <core/types.hh>
// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.fwd.hh>

// C++ headers
#include <string>
#include <iostream>
#include <map>

namespace core {
namespace io {
namespace nmr {

/// @brief reads residue number, atom name and chain ID for one spin
///        stores it in an AtomSelection object
void
read_atom_selection_from_string(
	std::string const & str,
	AtomSelection & atom
) /* throw (utility::excn::BadInput) */ ;

/// @brief read residue number, atom name and chain ID  for multiple
///        spins within a selection and store them in a vector
///        of AtomSelection objects
Size
read_selection_group_from_string(
	std::string const & str,
	utility::vector1< AtomSelection > & group,
	Size offset
) /* throw (utility::excn::BadInput) */ ;

/// @brief reads complete pcs data file
/// @details each line contains an AtomSelection for a protein spin, the measured pcs value and its error
void
read_pcs_datafile(std::string const & filename,
	utility::vector1< utility::vector1< AtomSelection > > & spins,
	utility::vector1< Real > & values,
	utility::vector1< Real > & errors
);

/// @brief reads complete rdc data file
/// @details each line contains two AtomSelections of the interacting protein spins,
///          the measured rdc value and its error
void
read_rdc_datafile(std::string const & filename,
	utility::vector1< utility::vector1< AtomSelection > > & spinsA,
	utility::vector1< utility::vector1< AtomSelection > > & spinsB,
	utility::vector1< Real > & values,
	utility::vector1< Real > & errors
);

/// @brief reads complete pre data file
/// @details each line contains an AtomSelection for a protein spin, the measured pre value and its error
void
read_pre_datafile(std::string const & filename,
	utility::vector1< utility::vector1< AtomSelection > > & spins,
	utility::vector1< Real > & values,
	utility::vector1< Real > & errors
);

/// @brief utility function to read in key-value pairs from the NMR data main input file
void
read_key_value_pair_from_line(std::string const & line,
	std::string const & key,
	std::map< std::string, std::string > & key_value_map,
	Size line_number
);

/// @brief strip brackets from string of dataset list
std::string
strip_brackets(std::string const & str);

/// @brief convert PCS dataset string and get vector of items
///        [ file, lanthanide, weight, single value weighting, averaging type,
///          computation type, xM, yM, zM, Xax, Xrh, alpha, beta, gamma ]
utility::vector1<std::string>
read_pcs_dataset_params_list(std::string const & str);

/// @brief convert RDC dataset string and get vector of items
///        [ file, weight, single value weighting ]
utility::vector1<std::string>
read_rdc_dataset_params_list(std::string const & str);

/// @brief convert PRE dataset string and get vector of items
///        [ file, weight, single value weighting, rate type, B0 ]
utility::vector1<std::string>
read_pre_dataset_params_list(std::string const & str);

/// @brief convert RDC tensor values from string
utility::vector1<Real>
read_rdc_tensor_values_from_string(std::string const & str);

/// @brief read gridsearch values from string
///        [ atom1, atom2, distance center-atom1, stepsize, inner radius, outer radius ]
utility::vector1<std::string>
read_gridsearch_values_from_string(std::string const & str);

} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_util_hh
