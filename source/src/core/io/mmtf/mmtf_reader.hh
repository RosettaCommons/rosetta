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
#include <core/io/pdb/pdb_reader.hh>  // TODO: Pull out pseudo-duplicated code and move to sfr_storage.cc.

// When you move PDBReader and PoseUnbuilder, take these.
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueConnection.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>

#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/StructFileRep.hh>

#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>

#include <core/io/Remarks.hh>
#include <core/io/mmtf/util.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/string_constants.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <mmtf.hpp>

// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>

namespace core {
namespace io {
namespace mmtf {

// @brief This actually adds the SSBondInformation/LinkInformation
//        templated to work only with LinkInformation and SSBondInformation
template < typename T >
inline void
add_xbond_information(
	std::map< std::string, utility::vector1< T > >& xbond_map,
	core::io::AtomInformation const & atm_1,
	core::io::AtomInformation const & atm_2);

/// @brief This distributes what type of bonds to make based on the mmtf bonds via add_xbond_information
void
add_link_and_ss_information(
	::mmtf::StructureData const & sd,
	core::io::StructFileRep & sfr,
	std::vector< core::io::AtomInformation > const & all_AIs,
	core::Size const atmSerial);


/// @brief This adds actuals bonds from the mmtf object to connected_indicies/connected_orders
void
add_bond_information(::mmtf::StructureData const & sd,
	std::vector< core::io::AtomInformation > & all_AIs,
	std::map<core::Size, sd_index> const & atom_num_to_sd_map,
	core::Size const atomSerialMax);

/// @brief This adds TERs based on what is bound to the C of each AA.
//         if something isn't bound, add a TER!
void
add_ters_via_bonds(std::vector< core::io::AtomInformation > & all_AIs);

/// @brief load heterogen data: heterogen_names && residue_type_base_names
//  @note this is stored in extraProperties["rosetta::heterogen_names"] and
//        extraProperties["rosetta::residue_type_base_names"]
void
load_heterogen_info(
	::mmtf::MapDecoder const & md,
	core::io::StructFileRep & sfr);

/// @brief makes a linear vector of all AtomInformation from the mmtf data
std::vector< core::io::AtomInformation >
make_all_atom_information(::mmtf::StructureData const & sd,
	core::Size & atomSerial,
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

/// @brief master create sfr
core::io::StructFileRepOP
create_sfr_from_mmtf_filename(
	std::string const & filename,
	core::io::StructFileReaderOptions const & options);

} // core
} // io
} // mmtf
#endif  // INCLUDED_core_io_mmtf_mmtf_reader_HH
