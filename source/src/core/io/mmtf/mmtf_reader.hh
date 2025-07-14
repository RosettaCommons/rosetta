// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmtf/mmtf_reader.hh
/// @author Daniel Farrell (danpf@uw.edu)


#ifndef INCLUDED_core_io_mmtf_mmtf_reader_HH
#define INCLUDED_core_io_mmtf_mmtf_reader_HH

// Unit headers
#include <core/io/StructFileReaderOptions.hh>

// Package headers

// When you move PDBReader and PoseUnbuilder, take these.

#include <core/io/AtomInformation.hh>
#include <core/io/ResidueInformation.hh>
#include <core/io/StructFileRep.fwd.hh>


#include <core/io/mmtf/util.hh>

// Project headers
#include <core/types.hh>

// Basic headers

// Numeric headers

// Utility headers
#include <utility/vector1.hh>

// Numeric headers

// External headers
#include <mmtf.hpp>

// C++ headers

namespace core {
namespace io {
namespace mmtf {

// @brief This actually adds the SSBondInformation/LinkInformation
//        templated to work only with LinkInformation and SSBondInformation
template < typename T >
inline void
add_xbond_information(
	std::map< ResID, utility::vector1< T > >& xbond_map,
	core::io::AtomInformation const & atm_1,
	core::io::AtomInformation const & atm_2);

/// @brief This distributes what type of bonds to make based on the mmtf bonds via add_xbond_information
void
add_link_and_ss_information(
	::mmtf::StructureData const & sd,
	core::io::StructFileRep & sfr,
	std::vector< core::io::AtomInformation > const & all_AIs,
	std::vector< core::Size > const & ai_to_model,
	core::Size const & current_model);


/// @brief This adds actuals bonds from the mmtf object to connected_indicies/connected_orders
void
add_bond_information(::mmtf::StructureData const & sd,
	std::vector< core::io::AtomInformation > & all_AIs,
	std::map<core::Size, sd_index> const & atom_num_to_sd_map);

/// @brief This adds TERs based on what is bound to the C of each AA.
//         if something isn't bound, add a TER!
void
add_ters_via_bonds(std::vector< core::io::AtomInformation > & all_AIs,
	std::vector< core::Size > const & ai_to_model);


/// @brief set_model_index_if_not_empty does what the name says :p
/// useful for loading various datas from modelProperties
template < typename T >
void
set_model_index_if_not_empty(core::Size const model_index,
	std::string const & info_tag,
	std::vector< T > const & all_model_data, T & target_data);


void
read_extra_data(utility::vector1< core::io::StructFileRepOP > & sfrs,
	::mmtf::StructureData const & sd);

/// @brief load heterogen data: heterogen_names && residue_type_base_names
//  @note this is stored in modelProperties["rosetta::heterogen_names"] and
//        modelProperties["rosetta::residue_type_base_names"]
void
load_heterogen_info(
	::mmtf::MapDecoder const & md,
	core::io::StructFileRep & sfr,
	core::Size const model_index);


/// @brief makes a linear vector of all AtomInformation from the mmtf data
std::vector< core::io::AtomInformation >
make_all_atom_information(::mmtf::StructureData const & sd,
	std::vector< core::Size > & model_indexes,
	std::map< core::Size, core::Size > & model_index_to_starting_index,
	StructFileReaderOptions const & options );

/// @brief makes a single AtomInformation from the mmtf data
core::io::AtomInformation
make_atom_information(
	::mmtf::StructureData const &sd,
	::mmtf::GroupType const & group,
	int const groupAtomIndex,
	core::Size const atomIndex,
	core::Size const atomSerial,
	core::Size const groupIndex,
	core::Size const chainIndex,
	utility::vector1<char> & known_chainIDs,
	core::io::StructFileReaderOptions const & options );

template < typename T >
void
set_model_index_if_not_empty(core::Size const model_index,
	std::string const & info_tag,
	std::vector< T > const & all_model_data, T & target_data);

/// @brief create sfr from create_sfrs_from_mmtf_filename
core::io::StructFileRepOP
create_sfr_from_mmtf_filename(
	std::string const & mmtf_filename,
	core::io::StructFileReaderOptions const & options);

/// @brief master sfr creator
/// @note: mmtf format requires us to iterate completely because order
///        doesn't matter.  So it makes more sense to iterate over all
///        and then trim back your selection at the end.
utility::vector1< core::io::StructFileRepOP >
create_sfrs_from_mmtf_filename(
	std::string const & mmtf_filename,
	core::io::StructFileReaderOptions const & options,
	utility::vector1< core::Size > const & model_indexes);

///@brief Read extra data from StructureData into the SFR
///@details SD non-const as the decode method of mmTF is NOT CONST!
void
read_extra_data(
	::mmtf::StructureData & sd,
	core::io::StructFileRep & sfr);

} // core
} // io
} // mmtf

#endif  // INCLUDED_core_io_mmtf_mmtf_reader_HH
