// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/mmtf_write.hh
/// @brief Functions for MMTF writing.
/// @author Daniel Farrell (danpf@uw.edu)


#ifndef INCLUDED_core_io_mmtf_mmtf_writer_HH
#define INCLUDED_core_io_mmtf_mmtf_writer_HH

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/StructFileRepOptions.fwd.hh>

// Project header
#include <core/pose/Pose.hh>
#include <core/io/mmtf/util.hh>

#include <utility/io/ozstream.hh>

#include <mmtf.hpp>

namespace core {
namespace io {
namespace mmtf {

//////////////////////////////////////////////////////////////
// dump_mmtf setup functions
//////////////////////////////////////////////////////////////

/// @brief Dump a MMTF from a pose to a file.
///  Return success or failure.
bool
dump_mmtf( core::pose::Pose const & pose,
	std::string const & file_name,
	StructFileRepOptionsCOP options =  StructFileRepOptionsCOP( new StructFileRepOptions ) );


///@brief Dump mmtf to an outstream, optionally passing options and returning a StructFileRep for further processing
StructFileRepOP
dump_mmtf( core::pose::Pose const & pose,
	std::ostream & out,
	StructFileRepOptionsCOP options =  StructFileRepOptionsCOP( new StructFileRepOptions )
);


/// @brief Dump a MMTF from a pose, optionally extracting extra info.
/// NOTE: DEPRECATED.  LEGACY JD2 function to use JOB data, which we no longer use in JD3
///
StructFileRepOP
dump_mmtf(
	core::pose::Pose const & pose,
	std::string const &jd2_job_data,
	utility::io::ozstream & out);


/// @brief dump_mmtf function.
///  Create the sfr from pose using the PoseToStructFileRepConverter class.
///  Return success or failure.
/// @param  file_name  The file name to write the mmtf to
/// @param  sfr  The StructFileRepOP to build the mmtf from
/// @param  options  The StructFileRepOptions used to build the sfr
///
bool
dump_mmtf(
	std::string const & file_name,
	StructFileRepOP sfr,
	StructFileRepOptions const & options );


//////////////////////////////////////////////////////////////
// dump_mmtf helper functions
//////////////////////////////////////////////////////////////



/// @brief add_if_not_empty -- adds to a sd {string: msgpack::object}
///  templated to allow us to add a variety of containers
///
/// @param  given_name  The name to use in the data map if not empty
/// @param  content     The content to add if not empty
/// @param  target_map  Where to add the content if it's not empty
/// @param  zone        Required zone to keep msgpack happy
template< typename T >
void
add_if_not_empty(
	std::string const & given_name,
	T const & content,
	std::map<std::string, msgpack::object> & target_map,
	msgpack::zone & zone);


/// @brief resize_add_if_not_empty -- adds T to std::vector< T > if T isn't empty
///  templated to allow us to add a variety of containers
///
/// @param  data        The data to add to destination, if !data.empty()
/// @param  destination The destination of the content if !empty
/// @param  data_index  Where to add the content if !empty (will resize destination)
template< typename T >
void
resize_and_add_if_not_empty(
	T const & data,
	std::vector< T > & destination,
	core::Size const & data_index,
	core::Size const & max_size);


/// @brief Set StructFileRepOptions defauls for mmtf
///
/// @param core::io::StructFileRepOptions& options to modify for mmtf
///
void
set_mmtf_default_options(core::io::StructFileRepOptions & options);

/// @brief group chain_atoms into groups based on chain_id and resseq
aiChain
make_chain(utility::vector0<AtomInformation> const & chain_atoms);

///
/// @brief Convert linear sfr to vec0[vec0[vec0[AtomInformation]]]
///
/// @param core::io::StructFileRep const& to build aiPodel from
///
aiPose
aiPose_from_sfr(core::io::StructFileRep const & sfr);

///
/// @brief Convert linear sfr to vec0[vec0[vec0[vec0[AtomInformation]]]]
///
/// @param core::io::StructFileRep const& to build aiModels from
aiModels
aiModels_from_sfrs(utility::vector1< core::io::StructFileRepOP > const & sfrs);

///
/// @brief Add extra info to the structureData if we are asked to add it
///        includes: heterogen info
///
/// @param utility::vector1< core::io::StructFileRepOP > const & to get infos from
/// @param ::mmtf::StructureData const& the StructureData to add to
///
void
add_extra_data(
	::mmtf::StructureData & sd,
	utility::vector1< core::io::StructFileRepOP > const & sfrs,
	core::io::StructFileRepOptions const & options);

//////////////////////////////////////////////////////////////
// real dump_mmtf function
//////////////////////////////////////////////////////////////


::mmtf::StructureData
sfrs_to_sd(
	utility::vector1< core::io::StructFileRepOP > sfrs,
	core::io::StructFileRepOptions const & options);

::mmtf::StructureData
sfr_to_sd(
	core::io::StructFileRepOP sfrs,
	core::io::StructFileRepOptions const & options);


/// @brief Main dump_mmtf function.
///  Create the sfr from pose using the PoseToStructFileRepConverter class.
void
dump_mmtf(
	std::ostream & out,
	utility::vector1< core::io::StructFileRepOP > sfrs,
	core::io::StructFileRepOptions const & options);


/// @brief dump single mmtf function.
/// @note  We sadly have to use this to call the function
///        that dumps a vector of SFRs due to the mmtf file format.
///        It only makes sense because of the compression scheme used.
void
dump_mmtf(
	std::ostream & out,
	StructFileRepOP sfr,
	StructFileRepOptions const & options );


}  // namespace mmtf
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_mmtf_mmtf_writer_HH
