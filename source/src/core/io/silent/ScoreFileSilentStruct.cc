// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/ScoreFileSilentStruct.cc
///
/// @brief Representation of PDB files in a silent-file format.
/// @author James Thompson

// C++ Headers
#include <utility>
#include <string>
#include <sstream>

// mini headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Residue.fwd.hh>

#include <ObjexxFCL/FArray2D.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/io/mpistream.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

using namespace ObjexxFCL;

static thread_local basic::Tracer tr( "core.io.silent.ScoreFileSilentStruct" );

ScoreFileSilentStruct::ScoreFileSilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	fill_struct( pose, tag );
	decoy_tag( tag );
} // ScoreFileSilentStruct

void ScoreFileSilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	basic::ProfileThis doit( basic::SILENT_FILL_STRUCT );
	energies_from_pose( pose );
	decoy_tag( tag );
	if ( tag == "empty_tag" ) set_tag_from_pose( pose );
	sequence( pose.sequence() );
}

void
ScoreFileSilentStruct::print_header( std::ostream & out ) const {
	print_score_header( out );
}

void ScoreFileSilentStruct::print_conformation( std::ostream & /* out */ ) const {
	// do nothing!
}


bool ScoreFileSilentStruct::init_from_lines(
	utility::vector1< std::string > const & lines,
	SilentFileData & container
) {
	bool success( false );
	using std::string;
	using utility::vector1;

	vector1< std::string > energy_names_;
	vector1< std::string >::const_iterator iter = lines.begin();
	if ( iter->substr(0,9) == "SEQUENCE:" ) iter++; // ignore sequence for now
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

	success = true;
	return success;
} // init_from_lines

void ScoreFileSilentStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & /* residue_set */
) const {
	basic::ProfileThis doit( basic::SILENT_FILL_POSE );

	// put energies into pose!
	energies_into_pose( pose );
} // fill_pose

Real ScoreFileSilentStruct::get_debug_rmsd() {
	tr.Error << "get_debug_rmsd stubbed out!" << std::endl;
	return 0.0;
}

ObjexxFCL::FArray2D< Real >
ScoreFileSilentStruct::get_CA_xyz() const {
	tr.Error << "get_CA_xyz() stubbed out!" << std::endl;
	return ObjexxFCL::FArray2D< Real > ( 3, 1 );
}

ScoreFileSilentStruct & ScoreFileSilentStruct::operator= (
	ScoreFileSilentStruct const &
)
{
	utility_exit_with_message( "called ScoreFileSilentStruct::operator=)" );
	exit(0); //  just to keep the compiler happy
}

} // namespace silent
} // namespace io
} // namespace core
