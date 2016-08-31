// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmcif/util.hh
/// @brief Functions for MMCIF writing.
/// @author Andy Watkins (andy.watkins2@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_io_mmcif_cif_writer_HH
#define INCLUDED_core_io_mmcif_cif_writer_HH

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/StructFileRepOptions.fwd.hh>

// Project header
#include <core/pose/Pose.hh>

#include <utility/io/ozstream.hh>

namespace core {
namespace io {
namespace mmcif {

/// @brief Dump an mmCIF from a pose to a file.
///  Return success or failure.
bool
dump_cif( core::pose::Pose const & pose,
	      std::string const & file_name,
		  StructFileRepOptionsCOP options =  StructFileRepOptionsCOP( new StructFileRepOptions ) );

///@brief Dump cif to an outstream, optionally passing options and returning a StructFileRep for further processing
StructFileRepOP
dump_cif( core::pose::Pose const & pose,
		  std::ostream & out,
	      StructFileRepOptionsCOP options =  StructFileRepOptionsCOP( new StructFileRepOptions )
);

///@brief Return an mmCIF-format string from a pose with defaults.
std::string
dump_cif( core::pose::Pose const & pose );

	
/// @brief Dump an mmCIF from a pose, optionally extracting extra info.
/// NOTE: DEPRECATED.  LEGACY JD2 function to use JOB data, which we no longer use in JD3
///
StructFileRepOP
dump_cif(
	core::pose::Pose const & pose,
	std::string const &jd2_job_data,
	utility::io::ozstream & out);




//////////////////////////////////////////////////////////////

/// @brief Main dump_cif function.
///  Return a string with contents of a CIF file, extracted from a Pose to a StructFileRep via PoseToStructFileRepConverter
std::string
dump_cif(
	StructFileRepOP sfr,
	StructFileRepOptions const & options );
	
	
/// @brief Main dump_cif function.
///  Create the sfr from pose using the PoseToStructFileRepConverter class.
///  Return success or failure.
bool
dump_cif(
	std::string const & file_name,
	StructFileRepOP sfr,
	StructFileRepOptions const & options );



/// @brief Main dump_cif function.
///  Create the sfr from pose using the PoseToStructFileRepConverter class.
void
dump_cif(
	std::ostream & out,
	StructFileRepOP sfr,
	StructFileRepOptions const & options );


}  // namespace mmcif
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_mmcif_cif_writer_HH


