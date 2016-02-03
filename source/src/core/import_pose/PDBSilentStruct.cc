// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/PDBSilentStruct.cc
///
/// @brief Representation of PDB files in a silent-file format.
/// @author James Thompson

// C++ Headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <sstream>

// mini headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/import_pose/PDBSilentStruct.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/Pose.hh>


#include <core/conformation/Residue.fwd.hh>

// ObjexxFCL
#include <ObjexxFCL/FArray2D.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <basic/options/option.hh>

namespace core {
namespace import_pose {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer tr( "core.io.silent.PDBSilentStruct" );

PDBSilentStruct::PDBSilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	fill_struct( pose, tag );
} // PDBSilentStruct

void PDBSilentStruct::print_header( std::ostream & out ) {
	print_score_header( out );
}

void PDBSilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	decoy_tag( tag );
	if ( tag == "empty_tag" ) set_tag_from_pose( pose );

	energies_from_pose( pose );

	io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
	pose_to_sfr.init_from_pose( pose );
	sfr_ = pose_to_sfr.sfr();

	sequence( pose.sequence() );
}

bool PDBSilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	core::io::silent::SilentFileData & container
) {
	bool success( false );
	using std::string;
	using utility::vector1;
	using namespace core::io::silent;

	vector1< std::string > energy_names_;
	vector1< std::string >::const_iterator iter = lines.begin();
	if ( iter->substr(0,9) == "SEQUENCE:" ) ++iter; // ignore sequence for now
	if ( iter->substr(0,6) != "SCORE:" ) {
		// get sequence and scorename data from the silent-file data object, because I don't have it!
		EnergyNamesOP enames = EnergyNamesOP(
			utility::pointer::static_pointer_cast< core::io::silent::EnergyNames > ( container.get_shared_silent_data( energynames ) )
		);

		energy_names_ = enames->energy_names();
	} else {
		// get scorename data from the first two lines provided, put into container
		// for further use by other SilentStruct objects.

		EnergyNamesOP enames( new EnergyNames( *iter ) );
		container.set_shared_silent_data( energynames, enames  );
		energy_names_ = enames->energy_names();
	} // get header information

	std::string concatenated_pdb_info; // concatenated pdb information
	for ( vector1< string >::const_iterator end = lines.end(); iter != end; ++iter ) {
		string tag;
		std::istringstream line_stream( *iter );

		if ( iter->substr(0,7) == "SCORE: " ) { // SCORE: line with values from this structure.
			std::string tag;
			line_stream >> tag;
			if ( line_stream.fail() || tag != "SCORE:" ) {
				tr.Error << "bad format in first score line of silent file" << std::endl;
				tr.Error << "line = " << *iter << std::endl;
				tr.Error << "tag = " << tag << std::endl;
			}

			vector1< string >::const_iterator energy_iter;
			for ( energy_iter = energy_names_.begin();
					energy_iter != energy_names_.end(); ++energy_iter
					) {
				line_stream >> tag;
				if ( *energy_iter != "description" ) { // currently the only text-based field, might change in future.
					Real score_val = (Real) float_of( tag );
					add_energy( *energy_iter, score_val );
				} else {
					line_stream >> tag;
				}
			} // for ( energy_iter ... )
			decoy_tag( tag ); // decoy_tag should be last column of this line.
		} else if ( iter->substr(0,6) == "REMARK" ) {
			// do nothing!
		} else {
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[ out::file::silent_preserve_H ]() || iter->substr(13,1) != "H" ) {
				// FileData needs \n's
				concatenated_pdb_info += iter->substr( 0, 79 ) + '\n';
				string temp_decoy_tag  = iter->substr( 80 );
			}
		} // else
	} // for ( iter ... )

	pdb_lines_ = concatenated_pdb_info;
	sfr_ = core::io::pdb::create_sfr_from_pdb_file_contents( concatenated_pdb_info ).clone();

	success = true;
	return success;
} // init_from_lines

void PDBSilentStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) const {
	// TODO: TEMP: for Andy to fix.
	io::StructFileRepOP sfr_op = sfr_->clone();
	core::import_pose::build_pose( sfr_op, pose, residue_set );
	core::import_pose::read_additional_pdb_data( pdb_lines_, pose, sfr_op );
} // fill_pose

void PDBSilentStruct::fill_pose(
	core::pose::Pose & pose
) const {

	using namespace core::chemical;
	ResidueTypeSetCOP residue_set;
	residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	fill_pose( pose, *residue_set );
	finish_pose( pose );
} // fill_pose


void PDBSilentStruct::print_conformation( std::ostream & output ) const {
	using std::string;
	// TODO: TEMP: for Andy to fix.
	io::StructFileRepOP sfr_op = sfr_->clone();
	string data = core::io::pdb::create_pdb_contents_from_sfr(*sfr_);
	output.write( data.c_str(), data.size() );
} // print_conformation

Real PDBSilentStruct::get_debug_rmsd() {
	tr.Error << "get_debug_rmsd stubbed out!" << std::endl;
	return 0.0;
}

ObjexxFCL::FArray2D< Real >
PDBSilentStruct::get_CA_xyz() const {
	tr.Error << "PDBSilentStruct::get_CA_xyz" << std::endl;
	return FArray2D< Real > ( 3, 1 );
}

PDBSilentStruct & PDBSilentStruct::operator= (
	PDBSilentStruct const &
)
{
	utility_exit_with_message( "called ProteinSilentStruct::operator=)" );
	exit(0); //  just to keep the compiler happy
}

} // namespace import_pose
} // namespace core
