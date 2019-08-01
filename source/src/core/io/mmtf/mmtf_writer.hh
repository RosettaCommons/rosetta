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
/// @param core::io::StructFileRep const& to build aiPose from
///
aiPose
aiPose_from_sfr(core::io::StructFileRep const & sfr);

///
/// @brief Add heterogen info to the structureData if we are asked to add it
///
/// @param core::io::StructFileRep const& to get heterogens from
/// @param ::mmtf::StructureData const& the StructureData to add to
///
void
add_heterogen_info_to_sd(
	::mmtf::StructureData & sd,
	core::io::StructFileRep const & sfr,
	core::io::StructFileRepOptions const & options);

//////////////////////////////////////////////////////////////
// real dump_mmtf function
//////////////////////////////////////////////////////////////

/// @brief Main dump_mmtf function.
///  Create the sfr from pose using the PoseToStructFileRepConverter class.
void
dump_mmtf(
	std::ostream & out,
	StructFileRepOP sfr,
	StructFileRepOptions const & options );


}  // namespace mmtf
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_mmtf_mmtf_writer_HH
