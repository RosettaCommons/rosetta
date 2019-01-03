// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/util.cc
/// @brief   Util functions for Input and Output of NMR data.
/// @details last Modified: 06/08/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/io/nmr/AtomSelection.hh>

// Package headers

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

namespace core {
namespace io {
namespace nmr {

static basic::Tracer TR( "core.io.nmr.util" );

/// @brief reads residue number, atom name and chain ID for one spin
///        stores it in an AtomSelection object
void
read_atom_selection_from_string(
	std::string const & str,
	AtomSelection & atom
) /* throw (utility::excn::BadInput) */
{
	Size rsd(0);
	std::string name("");
	char chain('^');
	std::istringstream iss(str);
	if ( !(iss >> rsd >> name >> chain) ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "AtomSelection reading error.");
	} else {
		runtime_assert_msg( rsd > 0, "ERROR: Residue number must be positive.");
		atom.set_rsd(rsd);
		atom.set_atom(name);
		atom.set_chain(chain);
	}
}

/// @brief read residue number, atom name and chain ID  for multiple
///        spins within a selection and store them in a vector
///        of AtomSelection objects
Size
read_selection_group_from_string(
	std::string const & str,
	utility::vector1< AtomSelection > & group,
	Size offset
) /* throw (utility::excn::BadInput) */
{
	AtomSelection atom;
	Size start = str.find('(', offset);
	Size end = str.find(')', offset);
	if ( start != std::string::npos && end != std::string::npos ) {
		std::string substring = str.substr(start+1, end-start-1);
		Size start_atomsel(0);
		Size found_or(0);
		do {
			found_or = substring.find("or", start_atomsel);
			Size const end_atomsel = ( found_or == std::string::npos ) ? substring.size() : found_or;
			try {
				read_atom_selection_from_string(substring.substr(start_atomsel,end_atomsel-start_atomsel), atom);
			} catch (utility::excn::BadInput & excn) {
				TR.Warning << "ERROR: " << excn.msg() << std::endl;
				group.clear();
				throw;
			}
			group.push_back(atom);
			// Make sure that the AtomSelections to be grouped have the same atom name (e.g. all "H")
			if ( group.size() > 1 ) {
				if ( atom.get_atom() != group.back().get_atom() ) {
					group.clear();
					throw CREATE_EXCEPTION(utility::excn::BadInput, "AtomSelections do not have the same atom name.");
				}
			}
			start_atomsel = found_or+2;
		} while (found_or != std::string::npos);
		return end;
	} else {
		return std::string::npos;
	}
}


/// @brief reads complete pcs data file
/// @details each line contains an AtomSelection or a vector of AtomSelections for equivalent spins
///          on equivalent subunits, the measured pcs value and its error
void
read_pcs_datafile(
	std::string const & filename,
	utility::vector1< utility::vector1< AtomSelection > > & spins,
	utility::vector1< Real > & values,
	utility::vector1< Real > & errors
)
{
	std::ifstream infile;
	std::string line;
	Real nmr_value(0);
	Real error(0);
	Size line_number(0);
	utility::vector1< AtomSelection > vec;
	Size found(0);

	TR.Info << "Opening file '" << filename <<"' " << std::endl;
	infile.open(filename.c_str(), std::ios::in);
	if ( !infile.is_open() ) {
		utility_exit_with_message( "Unable to open PCS/PRE data file: " + filename );
	}
	while ( std::getline(infile, line) ) {
		++line_number;
		vec.clear();

		// Remove leading and trailing whitespace and skip comments and blank lines.
		utility::trim( line, " \t\n" );
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) {
			continue;
		}

		// Read the first vector of AtomSelection objects
		try {
			found = read_selection_group_from_string(line, vec, 0);
		} catch (utility::excn::BadInput  & excn) {
			TR.Warning << "ERROR: " << excn.msg() << std::endl;
			TR.Warning << "Skip line " << line_number << std::endl;
			continue;
		}
		if ( found == std::string::npos ) {
			TR.Info << "ERROR: No AtomSelection found. Did you forget to place it into brackets \"( ... )\" ?" << std::endl;
			TR.Info << "Skip line " << line_number << std::endl;
			continue;
		}
		// Read the pcs value and error
		std::istringstream iss(line.substr(found+1, line.size()-found-1));
		if ( !(iss >> nmr_value >> error) ) {
			TR.Info << "ERROR: Could not read in PCS/PRE value and error." << std::endl;
			TR.Info << "Skip line " << line_number << std::endl;
			continue;
		}

		spins.push_back(vec);
		values.push_back(nmr_value);
		errors.push_back(error);
	}
	runtime_assert( (spins.size() == values.size()) && (spins.size() == errors.size()) );
	infile.close();
}

/// @brief reads complete rdc data file
/// @details each line contains two AtomSelections of the interacting protein spins,
///          the measured rdc value and its error
void
read_rdc_datafile(
	std::string const & filename,
	utility::vector1< utility::vector1< AtomSelection > > & spinsA,
	utility::vector1< utility::vector1< AtomSelection > > & spinsB,
	utility::vector1< Real > & values,
	utility::vector1< Real > & errors
)
{
	std::ifstream infile;
	std::string line;
	Real nmr_value(0);
	Real error(0);
	Size line_number(0);
	utility::vector1< AtomSelection > vec1, vec2;
	Size found(0);

	TR.Info << "Opening file '" << filename <<"' " << std::endl;
	infile.open(filename.c_str(), std::ios::in);
	if ( !infile.is_open() ) {
		utility_exit_with_message( "Unable to open RDC data file: " + filename );
	}
	while ( std::getline(infile, line) ) {
		++line_number;
		vec1.clear();
		vec2.clear();

		// Remove leading and trailing whitespace and skip comments and blank lines.
		utility::trim( line, " \t\n" );
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) {
			continue;
		}

		// Read the first group of AtomSelection objects that contains the first spin(s) A
		try {
			found = read_selection_group_from_string(line, vec1, 0);
		} catch (utility::excn::BadInput  & excn) {
			TR.Warning << "ERROR: " << excn.msg() << std::endl;
			TR.Warning << "Skip line " << line_number << std::endl;
			continue;
		}
		if ( found == std::string::npos ) {
			TR.Info << "ERROR: No first AtomSelection found. Did you forget to place it into brackets \"( ... )\" ?" << std::endl;
			TR.Info << "Skip line " << line_number << std::endl;
			continue;
		}

		// Read the second group of AtomSelection objects that contains the second spin(s) B
		try {
			found = read_selection_group_from_string(line, vec2, ++found);
		} catch (utility::excn::BadInput  & excn) {
			TR.Warning << "ERROR: " << excn.msg() << std::endl;
			TR.Warning << "Skip line " << line_number << std::endl;
			continue;
		}
		if ( found == std::string::npos ) {
			TR.Info << "ERROR: No second AtomSelection found. Did you forget to place it into brackets \"( ... )\" ?" << std::endl;
			TR.Info << "Skip line " << line_number << std::endl;
			continue;
		}

		// Make sure that both AtomSelection groups vec1 and vec2 are of equal size.
		if ( vec1.size() != vec2.size() ) {
			TR.Info << "ERROR: The two RDC AtomSelection groups have not the same length." << std::endl;
			TR.Info << "Skip line " << line_number << std::endl;
			continue;
		}

		// Read the nmr value and error
		std::istringstream iss(line.substr(found+1, line.size()-found-1));
		if ( !(iss >> nmr_value >> error) ) {
			TR.Info << "ERROR: Could not read in RDC value and error." << std::endl;
			TR.Info << "Skip line " << line_number << std::endl;
			continue;
		}

		spinsA.push_back(vec1);
		spinsB.push_back(vec2);
		values.push_back(nmr_value);
		errors.push_back(error);
	}
	runtime_assert( (spinsA.size() == spinsB.size()) && (spinsA.size() == values.size()) && (spinsA.size() == errors.size()) );
	infile.close();
}

/// @brief reads complete pre data file
/// @details each line contains an AtomSelection for a protein spin, the measured pre value and its error
void
read_pre_datafile(
	std::string const & filename,
	utility::vector1< utility::vector1< AtomSelection > > & spins,
	utility::vector1< Real > & values,
	utility::vector1< Real > & errors
)
{
	// We can use the same function as for pcs here
	read_pcs_datafile(filename, spins, values, errors);
}

/// @brief utility function to read in key-value pairs from the NMR data main input file
void
read_key_value_pair_from_line(
	std::string const & line,
	std::string const & key,
	std::map< std::string, std::string > & key_value_map,
	Size line_number
)
{
	Size idx_key = line.find(key);
	Size idx_equal_sign = line.find("=", idx_key + key.size());
	if ( idx_equal_sign == std::string::npos ) {
		utility_exit_with_message("ERROR: No equal sign found after parameter \"" + key + "\" in NMR input file. Please provide input parameter \" key = value \".");
	} else {
		std::string value = utility::trim(line.substr(idx_equal_sign+1));
		if ( key_value_map.find(key) != key_value_map.end() ) {
			TR.Warning << "Found a second key-value pair for parameter \"" << key << "\". Only one " << key << " possible per NMR MULTISET." << std::endl;
			TR.Warning << "Skip line " << line_number << std::endl;
		} else {
			key_value_map.insert(std::make_pair(key, value));
		}
	}
}

/// @brief strip brackets from string of dataset list
std::string
strip_brackets(std::string const & str)
{
	std::string::const_iterator begin = str.begin();
	std::string::const_iterator end = str.end();
	for ( std::string::const_iterator p = str.begin(); p!=str.end(); ++p ) {
		if ( *p == '[' ) begin = p+1;
		else break;
	}
	for ( std::string::const_iterator p = str.end(); p!=begin; --p ) {
		if ( *(p-1) == ']' ) end = p-1;
		else break;
	}
	return std::string(begin,end);
}

/// @brief convert PCS dataset string and get vector of items
///        [ file, lanthanide, weight, single value weighting, averaging type,
///          computation type, xM, yM, zM, Xax, Xrh, alpha, beta, gamma ]
utility::vector1<std::string>
read_pcs_dataset_params_list(std::string const & str)
{
	std::string stripped_list = strip_brackets(str);
	utility::vector1<std::string> items = utility::string_split(stripped_list, ',');
	for ( std::string & i : items ) {
		utility::trim(i, " ");
		utility::trim(i, "\t");
	}
	if ( items.size() != 14 ) {
		utility_exit_with_message("ERROR while trying to read PCS dataset parameter list. List must contain exactly 14 items in the following order: "
			"[ filename, lanthanide, weight, single value weighting, averaging type, computation type, xM, yM, zM, Xax, Xrh, alpha, beta, gamma ]");
	}
	return items;
}

/// @brief convert RDC dataset string and get vector of items
///        [ file, weight, single value weighting ]
utility::vector1<std::string>
read_rdc_dataset_params_list(std::string const & str)
{
	std::string stripped_list = strip_brackets(str);
	utility::vector1<std::string> items = utility::string_split(stripped_list, ',');
	for ( std::string & i : items ) {
		utility::trim(i, " ");
		utility::trim(i, "\t");
	}
	if ( items.size() != 3 ) {
		utility_exit_with_message("ERROR while trying to read RDC dataset parameter list. List must contain exactly 3 items in the following order: "
			"[ filename, weight, single value weighting ]");
	}
	return items;
}

/// @brief convert PRE dataset string and get vector of items
///        [ file, weight, single value weighting, rate type, B0 ]
utility::vector1<std::string>
read_pre_dataset_params_list(std::string const & str)
{
	std::string stripped_list = strip_brackets(str);
	utility::vector1<std::string> items = utility::string_split(stripped_list, ',');
	for ( std::string & i : items ) {
		utility::trim(i, " ");
		utility::trim(i, "\t");
	}
	if ( items.size() != 5 ) {
		utility_exit_with_message("ERROR while trying to read PRE dataset parameter list. List must contain exactly 5 items in the following order: "
			"[ filename, weight, single value weighting, PRE rate type, field strength ]");
	}
	return items;
}

/// @brief convert RDC tensor values from string
utility::vector1<Real>
read_rdc_tensor_values_from_string(std::string const & str)
{
	using utility::string2Real;

	std::string stripped_list = strip_brackets(str);
	utility::vector1<std::string> items = utility::string_split(stripped_list, ',');
	for ( std::string & i : items ) {
		utility::trim(i, " ");
		utility::trim(i, "\t");
	}
	if ( items.size() != 5 ) {
		utility_exit_with_message("ERROR while trying to read RDC tensor values. List must contain exactly 5 values in the following order: "
			"[ Da, R, alpha, beta, gamma ]");
	}
	utility::vector1<Real> tensor_values(5);
	tensor_values[1] = string2Real(items[3]); // alpha
	tensor_values[2] = string2Real(items[4]); // beta
	tensor_values[3] = string2Real(items[5]); // gamma
	tensor_values[4] = string2Real(items[1]); // Da
	tensor_values[5] = string2Real(items[2]); // R
	return tensor_values;
}

/// @brief read gridsearch values from string
///        [ atom1, atom2, distance center-atom1, stepsize, inner radius, outer radius ]
utility::vector1<std::string>
read_gridsearch_values_from_string(std::string const & str)
{
	std::string stripped_list = strip_brackets(str);
	utility::vector1<std::string> items = utility::string_split(stripped_list, ',');
	for ( std::string & i : items ) {
		utility::trim(i, " ");
		utility::trim(i, "\t");
	}
	if ( items.size() != 6 ) {
		utility_exit_with_message("ERROR while trying to read gridsearch parameter list. List must contain exactly 6 values in the following order: "
			"[ atom1, atom2, distance center-atom1, stepsize, inner radius, outer radius ]");
	}
	return items;
}

} // namespace nmr
} // namespace io
} // namespace core
