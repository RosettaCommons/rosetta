// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pdb/file_data.cc
/// @brief  Method definitions for StructFileRep and related classes.
/// @author Sergey Lyskov

// Note: AVOID ACCESSING THE OPTIONS SYSTEM DIRECTLY IN THIS FILE, ESPECIALLY FOR PDB INPUT!
// Doing so will mean the Resource Manager may not work properly.
// Instead, modify StructFileRepOptions to include the option.


// Unit headers
#include <core/io/pdb/build_pose_as_is.hh>

// Package headers

#include <core/io/StructFileRep.hh>
//#include <core/io/pdb/file_data_fixup.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/StructFileReaderOptions.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>



// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// External headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>

#include <core/pose/Pose.hh> // AUTO IWYU For Pose


namespace core {
namespace io {
namespace pdb {

using core::Size;
using core::SSize;

using namespace ObjexxFCL::format;

using std::string;
using std::iostream;

// Tracer instance for this file
static basic::Tracer TR( "core.io.pdb.file_data" );


// StructFileRep to Pose ///////////////////////////////////////////////////////////

void
build_pose_from_pdb_as_is( pose::Pose & pose, std::string const & filename )
{
	StructFileReaderOptions options;
	build_pose_from_pdb_as_is( pose, filename, options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename,
	StructFileReaderOptions const & pdr_options
)
{
	using namespace chemical;
	ResidueTypeSetCOP restype_set( pose.residue_type_set_for_pose( FULL_ATOM_t ) );
	if ( restype_set == nullptr ) {
		restype_set = ChemicalManager::get_instance()->residue_type_set( FULL_ATOM_t );
	}
	build_pose_from_pdb_as_is( pose, *restype_set, filename, pdr_options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
)
{
	StructFileReaderOptions options;
	build_pose_from_pdb_as_is( pose, residue_set, filename, options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	StructFileReaderOptions const & pdr_options
)
{
	utility::io::izstream file( filename );
	if ( !file ) {
		TR.Error << "File:" << filename << " not found!" << std::endl;
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		TR.Debug << "read file: " << filename << std::endl;
	}
	build_pose_from_pdb_as_is( pose, residue_set, filename, file, pdr_options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	std::istream & file_contents,
	StructFileReaderOptions const & pdr_options
)
{
	std::string all_lines;
	utility::slurp( file_contents, all_lines );

	StructFileRep sfr = create_sfr_from_pdb_file_contents( all_lines, pdr_options );
	if ( sfr.filename() == "" ) {
		sfr.filename() = filename;
	}
	id::AtomID_Mask missing( false );
	//build_pose_as_is1( sfr, pose, residue_set, missing, pdr_options);
	pose_from_sfr::PoseFromSFRBuilder builder( residue_set.get_self_ptr(), pdr_options );
	builder.build_pose( sfr, pose );
}



} // namespace pdb
} // namespace io
} // namespace core
