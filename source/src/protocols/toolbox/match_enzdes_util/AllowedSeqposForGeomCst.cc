// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012


/// unit headeers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>

///project headers
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

//utility headers
#include <utility/io/izstream.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <sstream>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.match_enzdes_util.AllowedSeqposForGeomCst" );

AllowedSeqposForGeomCst::AllowedSeqposForGeomCst(
	utility::vector1< utility::vector1< Size > > const & seqpos_for_geomcst )
: seqpos_for_geomcst_( seqpos_for_geomcst )
{}

AllowedSeqposForGeomCst::AllowedSeqposForGeomCst(){}

AllowedSeqposForGeomCst::~AllowedSeqposForGeomCst(){}


/// @details this function supports two behaviours
/// 1. the same list for every geomcst. this means that if
/// the size of seqpos_for_geomcst_ == 1 and a bigger number
/// than 1 is passed in, it will return the first element
/// 2. different lists for every geomcst. this means that
/// the size of seqpos_for_geomcst_ is > 1, but if a larger
/// number than the vector size is passed in, it's unclear
/// what to do, so we exit
utility::vector1< core::Size > const &
AllowedSeqposForGeomCst::seqpos_for_geomcst( Size geomcst ) const {

	if ( seqpos_for_geomcst_.size() == 1 ) return seqpos_for_geomcst_[1];

	else if ( geomcst > seqpos_for_geomcst_.size() ) utility_exit_with_message("Asking for seqpos list of geomcst "+utility::to_string( geomcst )+", but lists only exist for "+utility::to_string( seqpos_for_geomcst_.size() )+" geomcst.");

	return seqpos_for_geomcst_[ geomcst ];
}


/// @details This function reads one of two files from the command line
/// depending on which options the user has provided.  These files define
/// the set of residues on the scaffold to consider as launch-points
/// for the geometric constraints (e.g. to consider as part of the
/// the catalytic core in an enzyme). The file meanings and formats are below.
///
/// 1. A list of residue id's to consider for all of the geometric constraints.
/// Such a file can be generated for a scaffold and then used along side any
/// enzyme-design .cst file.  It is scaffold dependent and constraint-file
/// independent.  The file should list the residue indexes for the scaffold
/// on one or more lines.  The file format does not support comments.
/// Residue id's start counting at 1; the input pdb resids are ignored.
/// (It's best to renumber your scaffold resids starting from 1 to avoid confusion)
///
/// Example.
/// <begin file>
/// 104 106 108 109 117 118 137 139 143 144 36 6 85 87 88 89 91 92 97
/// <end file>
///
/// 2. A list for each geometric constraint of the residues to consider.
/// Such a file allows the user to focus on particular residues for certain
/// geometric constraints for a particular scaffold.  Such a file should be
/// depends on both the scaffold and the match constraint file and cannot
/// be generalized across either multiple scaffolds or multiple constraint files.
/// The first line of the file begins with N_CST, followed by the number of geometric
/// contraints.  This must match the number of geometric constraints in the .cst file.
/// On each subsequent line, the geometric constraint ID is given, followed by
/// a colon and then followed by all of the residue ID's that should be considered
/// for that geometric constraint.  Each geometric constraint must appear on one
/// line in the file, though they may be listed in any order.  The file format does
/// not support comments.
/// flo jan 2010: it is also possible to specify that all positions in the scaffold
/// can be used for a certain constraint. see the example for cst 4 below
/// Example.
/// <begin file>
/// N_CST 3
/// 1: 9
/// 3: 9
/// 2: 6 7 9 11 12 14 15 17 18 21 22 23 25 26 38 40 43 46 47 49 53 54 57 60 61
/// 4: all
/// <end file>
///
void
AllowedSeqposForGeomCst::initialize_from_command_line( core::pose::PoseCOP pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::match;

	if ( option[ scaffold_active_site_residues ].user() && option[ scaffold_active_site_residues_for_geomcsts ].user()  ) {
		utility_exit_with_message( "Conflicting scaffold build point defintion: both "
			"the -match::scaffold_active_site_residues <filename> flag\n"
			"and the -match::scaffold_active_site_residues_for_geomcsts <filename> flag were found on the command line" );
	}

	seqpos_for_geomcst_.clear();

	if ( option[ scaffold_active_site_residues ].user() ) {
		//share_build_points_for_geomcsts_ = true;
		seqpos_for_geomcst_.resize( 1 );

		std::string filename = option[ scaffold_active_site_residues ]();
		utility::io::izstream istr( filename.c_str() );
		std::list< Size > upstream_build_resids;
		TR << "Reading match::scaffold_active_stie_residues " << filename << ":";
		while ( istr ) {
			Size resid( 0 );
			istr >> resid;
			if ( ! istr.bad() && resid != 0 ) {
				TR << " " << resid;
				upstream_build_resids.push_back( resid );
			}

		}
		if ( upstream_build_resids.empty() ) {
			utility_exit_with_message( "Failed to read any scaffold active site residues from file " + filename + " specified by the flag match::scaffold_active_stie_residues" );
		}

		TR << std::endl;
		//generic_pose_build_resids_.resize( upstream_build_resids.size() );
		seqpos_for_geomcst_[1].resize(  upstream_build_resids.size() );
		std::copy( upstream_build_resids.begin(), upstream_build_resids.end(),  seqpos_for_geomcst_[1].begin() );
	} else if ( option[ scaffold_active_site_residues_for_geomcsts ].user() ) {
		//share_build_points_for_geomcsts_ = false;

		std::string filename = option[ scaffold_active_site_residues_for_geomcsts ]();
		utility::io::izstream istr( filename.c_str() );
		//std::list< Size > upstream_build_resids;
		TR << "Reading match::scaffold_active_site_residues_for_geomcsts " << filename << std::endl;
		std::string ncsts_string; Size ncsts;
		if ( ! istr.good() ) {
			utility_exit_with_message( "Could not read first line from match::scaffold_active_site_residues_for_geomcsts " + filename );
		}
		istr >> ncsts_string;
		if ( istr.bad() ) {
			utility_exit_with_message( "Failed to read N_CST field in first line from match::scaffold_active_site_residues_for_geomcsts " + filename );
		}
		if ( ! istr.good() ) {
			utility_exit_with_message( "Unexpected end of file after reading N_CST field in first line from match::scaffold_active_site_residues_for_geomcsts " + filename );
		}
		if ( ncsts_string != "N_CST" ) {
			utility_exit_with_message( "Failed to read N_CST field in first line from match::scaffold_active_site_residues_for_geomcsts " + filename );
		}
		istr >> ncsts;
		if ( istr.bad() ) {
			utility_exit_with_message( "Failed to read the number of geometric constraints in first line from match::scaffold_active_site_residues_for_geomcsts " + filename );
		}
		//if ( ncsts != enz_input_data_->mcfi_lists_size() ) {
		//  utility_exit_with_message( "#geometric constraints disagreement between "
		// "match::scaffold_active_site_residues_for_geomcsts " + filename +
		// " (" + utility::to_string(ncsts) + ") and Enz .cst file: " + option[ geometric_constraint_file ]()() +
		// " (" + utility::to_string(enz_input_data_->mcfi_lists_size()) + ")" );
		//}
		//per_cst_pose_build_resids_.resize( ncsts );
		seqpos_for_geomcst_.resize( ncsts );

		std::string finish_the_line;
		istr.getline( finish_the_line );
		utility::vector1< Size > data_read_for_cst( ncsts, 0 );
		Size linenum = 2;
		while ( istr ) {
			Size geomcst_id( 0 );
			istr >> geomcst_id;
			if ( istr.eof() && geomcst_id == 0 ) break;

			if ( istr.bad() ) {
				utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Failed to read the geometric constraint id at the beginning of the line." );
			}
			if ( geomcst_id > ncsts ) {
				utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Geometric constraint id at the beginning of the line is out of range: "
					+ utility::to_string( geomcst_id ) + " > " +  utility::to_string( ncsts ) );
			}
			if ( data_read_for_cst[ geomcst_id ] != 0 ) {
				utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Residue list for geometric constraint #"
					+ utility::to_string( geomcst_id ) + " appeared already on line " +  utility::to_string( data_read_for_cst[ geomcst_id ] ) );
			}
			data_read_for_cst[ geomcst_id ] = linenum;
			std::string colon;
			istr >> colon;
			if ( colon != ":" ) {
				utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Failed to read the colon (:) separating the geometric constraint id from the rest of the resids." );
			}
			TR << std::endl << geomcst_id << " :";
			std::string first_token;
			istr >> first_token;
			if ( istr.bad() ) utility_exit_with_message( "Apparently there are no residues listed for geom cst " + utility::to_string( geomcst_id ) + " in file " + filename +".");
			if ( ( first_token == "ALL" ) || (first_token == "all" ) ) {
				if ( !pose ) utility_exit_with_message("AllowedSeqposForGeomCst requested to use all build pos for a certain constraint, but no pose passed into function, can't generate list.");
				TR << "All pose positions requested, using ";
				for ( core::Size seqpos = 1; seqpos <= pose->total_residue(); ++seqpos ) {
					if ( pose->residue(seqpos).is_protein() ) {
						seqpos_for_geomcst_[geomcst_id ].push_back( seqpos );
						TR << " " << seqpos;
					}
				}
				istr.getline( finish_the_line );
			} else {
				Size first_resid(0);
				std::istringstream firststr( first_token );
				firststr >> first_resid;
				if ( first_resid != 0 ) { //&& (first_resid <= upstream_pose_->total_residue() )){
					seqpos_for_geomcst_[geomcst_id].push_back( first_resid );
					TR << " " << first_resid;
				} else {
					utility_exit_with_message("Bad first residue listed for geomcst " + utility::to_string( geomcst_id ) + " in file " + filename +": " + first_token);
				}
				istr.getline( finish_the_line );
				if ( finish_the_line != "" ) {
					std::istringstream isstr( finish_the_line );
					while ( isstr.good() ) {
						Size resid( 0 );
						isstr >> resid;
						if ( isstr.eof() && resid == 0 ) break;
						if ( ! isstr.bad() ) {
							//if ( resid > 0 && resid <= upstream_pose_->total_residue() ) {
							TR << " " << resid;
							seqpos_for_geomcst_[ geomcst_id ].push_back( resid );
							//}
							//else if ( resid > upstream_pose_->total_residue() ){
							//std::cerr << std::endl << "ERROR parsing line fragment: " << finish_the_line << std::endl;
							//utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Requested upstream build resid of " + utility::to_string(resid) + " exceeds the number of residues in the pose ( " + utility::to_string( upstream_pose_->total_residue() ) + ")"  );
							//} else {
							//std::cerr << std::endl << "ERROR parsing line fragment: " << finish_the_line << std::endl;
							//utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Failed to read an integer."  );
							//}
						} else {
							std::cerr << std::endl << "ERROR parsing line fragment: " << finish_the_line << std::endl;
							utility_exit_with_message( "Reading line " + utility::to_string( linenum ) + " of " + filename + ". Only integers may be included." );
						}
					} //while loop over line
				}//if finish_the_line has stuff in it
			}//if all pos else
			++linenum;
		} //loop over lines
		TR << std::endl;
		bool any_absent( false );
		for ( Size ii = 1; ii <= ncsts; ++ii ) {
			if ( data_read_for_cst[ ii ] == 0 ) {
				std::cerr << "ERROR reading " << filename << ": did not find residue list for constraint # " << ii << std::endl;
				any_absent = true;
			}
		}
		if ( any_absent ) {
			utility_exit_with_message( "Failed to read a residue list for one or more constraints" );
		}
	}
}


}
}
}
