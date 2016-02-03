// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data.cc
/// @brief  Method definitions for StructFileRep and related classes.
/// @author Sergey Lyskov

// Note: AVOID ACCESSING THE OPTIONS SYSTEM DIRECTLY IN THIS FILE, ESPECIALLY FOR PDB INPUT!
// Doing so will mean the Resource Manager may not work properly.
// Instead, modify StructFileRepOptions to include the option.


// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>

// Package headers

#include <core/io/StructFileRepOptions.hh>
#include <core/io/StructFileRep.hh>
//#include <core/io/pdb/file_data_fixup.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/io/NomenclatureManager.hh>
#include <core/io/util.hh>
#include <core/io/NomenclatureManager.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedAtomID_Map.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/cryst/util.hh>

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
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/lexical_cast.hpp>

// C++ headers
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <utility>


namespace core {
namespace io {
namespace pdb {

using core::Size;
using core::SSize;

using core::chemical::chr_chains;

using basic::T;
using basic::Error;
using basic::Warning;

using ObjexxFCL::strip_whitespace;
using ObjexxFCL::stripped_whitespace;
using ObjexxFCL::rstripped_whitespace;
using namespace ObjexxFCL::format;

using std::string;
using std::iostream;

// Tracer instance for this file
static THREAD_LOCAL basic::Tracer TR( "core.io.pdb.file_data" );


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
	build_pose_from_pdb_as_is( pose, * ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ), filename, pdr_options );
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
