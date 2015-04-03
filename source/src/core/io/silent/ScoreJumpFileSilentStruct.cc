// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/ScoreJumpFileSilentStruct.cc
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
#include <core/io/silent/ScoreJumpFileSilentStruct.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.fwd.hh>

#include <ObjexxFCL/FArray2D.hh>

#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

using namespace ObjexxFCL;

static thread_local basic::Tracer tr( "core.io.silent.ScoreJumpFileSilentStruct" );

ScoreJumpFileSilentStruct::ScoreJumpFileSilentStruct(
	core::pose::Pose const & pose,
	std::string tag
) {
	fill_struct( pose, tag );
	decoy_tag( tag );
} // ScoreJumpFileSilentStruct

void ScoreJumpFileSilentStruct::fill_struct(
	core::pose::Pose const & pose,
	std::string tag
) {
	basic::ProfileThis doit( basic::SILENT_FILL_STRUCT );
	energies_from_pose( pose );
	decoy_tag( tag );
	if ( tag == "empty_tag" ) set_tag_from_pose( pose );
	sequence( pose.sequence() );

  fold_tree_ = pose.fold_tree();
  jumps_.clear();
  for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
    add_jump( pose.jump(nr) );
  }

}

void
ScoreJumpFileSilentStruct::print_header( std::ostream & out ) const {
  SilentStruct::print_header( out );
}

void ScoreJumpFileSilentStruct::print_conformation( std::ostream & output  ) const {
  output << "REMARK SOCREJUMP SILENTFILE\n";
if ( fold_tree().size() > 1 || fold_tree().num_jump() > 0 ) {
    output << "FOLD_TREE ";
    for ( kinematics::FoldTree::const_iterator
        it = fold_tree().begin(), it_end = fold_tree().end();
        it != it_end; ++it
    ) {
      output << *it;
    }
    //    output << fold_tree(); this produces a new-line --- wrong behaviour
    //    of fold_tree but I don't want to fix 1000 u-tracer unit-tests!
    output << ' ' << decoy_tag() << "\n";
  }
  for ( Size i = 1; i <= fold_tree().num_jump(); i++ ) {
    output << jump( i ) << ' ' << decoy_tag() << "\n";
  }
}


bool ScoreJumpFileSilentStruct::init_from_lines(
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

  for ( utility::vector1< std::string >::const_iterator end = lines.end(); iter != end; ++iter ) {
		std::string tag;
    std::istringstream line_stream( *iter );

	if ( iter->substr(0,6) == "REMARK" ) {
      std::string tag;
      std::string comment;
      std::string value;
      line_stream >> tag >> comment >> value;
      add_comment( comment, value );
      continue;  // skip comments
    }

	if ( iter->substr(0,7) == "SCORE: " ) {
      std::string tag;
      line_stream >> tag;
      if ( line_stream.fail() || tag != "SCORE:" ) {
        tr.Error << "bad format in first score line of silent file" << std::endl;
        tr.Error << "line = " << *iter << std::endl;
        tr.Error << "tag = " << tag << std::endl;
      }

      parse_energies( line_stream, energy_names_ );

    } else { // conformation lines

		if ( iter->substr(0,10) == "FOLD_TREE " ) {
        kinematics::FoldTree f;
        line_stream >> f;
        fold_tree( f ); // add fold-tree to this SilentStruct
        tr.Debug << "read fold-tree " << f; //"\n" is in fold-tree output
        tr.Debug << "reading " << f.num_jump() << " jumps " << std::endl;
        continue;
      } else if ( iter->substr(0,2) == "RT" ) {
        kinematics::Jump jump;
        line_stream >> jump;
        tr.Debug << "read jump " << jump << std::endl;
        add_jump( jump );
        continue;
      } else {
        tr.Error << "bad format in score_jump silent format and contains " << *iter << std::endl;
			}
		}//end of reading jump and RT
	}

	success = true;
	return success;
} // init_from_lines

void ScoreJumpFileSilentStruct::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & /* residue_set */
) const {
	basic::ProfileThis doit( basic::SILENT_FILL_POSE );

	// put energies into pose!
	energies_into_pose( pose );

  tr.Debug << "FOLD TREE: " << fold_tree();
  pose.fold_tree( fold_tree() );

	for ( Size nr = 1; nr <= fold_tree().num_jump(); nr++)  {
		      pose.set_jump( nr, jump( nr ) );
	}

} // fill_pose

Real ScoreJumpFileSilentStruct::get_debug_rmsd() {
	tr.Error << "get_debug_rmsd stubbed out!" << std::endl;
	return 0.0;
}

ObjexxFCL::FArray2D< Real >
ScoreJumpFileSilentStruct::get_CA_xyz() const {
	tr.Error << "get_CA_xyz() stubbed out!" << std::endl;
	return ObjexxFCL::FArray2D< Real > ( 3, 1 );
}

ScoreJumpFileSilentStruct & ScoreJumpFileSilentStruct::operator= (
	ScoreJumpFileSilentStruct const &
)
{
	utility_exit_with_message( "called ScoreJumpFileSilentStruct::operator=)" );
	exit(0); //  just to keep the compiler happy
}

} // namespace silent
} // namespace io
} // namespace core
