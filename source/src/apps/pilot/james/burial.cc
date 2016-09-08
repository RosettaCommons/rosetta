// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>


#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <apps/pilot/james/james_util.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option key includes


//Auto Headers
#include <basic/options/keys/OptionKeys.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );

	using core::Size;
	using core::Real;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;

	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	MetaPoseInputStream input = streams_from_cmd_line();

	utility::vector1< Real > cutoffs;
	//Size const lower(  5 );
	//Size const upper( 12 );
	//for ( Size ii = lower; ii <= upper; ++ii ) {
	//	cutoffs.push_back( static_cast< Real > ( ii ) );
	//}
	cutoffs.push_back(  6 );
	cutoffs.push_back( 10 );
	cutoffs.push_back( 12 );

	while( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );

		std::ostream & output( std::cout );

		output 	<< A( 10, "resi_idx" )
						<< A(  6, "resi"     );

		typedef vector1< Real >::const_iterator iter;
		for ( iter it = cutoffs.begin(), end = cutoffs.end(); it != end; ++it ) {
			std::string const column_title( "cen" + string_of(*it) );
			output << A( 10, column_title );
		}
		output << std::endl;

		// calculate burial
		vector1< vector1< int > > burial;
		for ( Size ii = 1; ii <= cutoffs.size(); ++ii ) {
			vector1< int > this_burial = calculate_burial( pose, cutoffs[ii] );
			burial.push_back( this_burial );
		}

		for ( unsigned int ii = 1; ii <= pose.size(); ++ii ) {
			core::conformation::Residue resi = pose.residue(ii);
			output << I( 10, ii ) << A(  6, resi.name1() );
			for ( Size jj = 1; jj <= cutoffs.size(); ++jj ) {
				output << I( 10, burial[jj][ii] );
			}
			output << std::endl;
		}	// for ( unsigned int i = 1; i <= pose.size(); ++i )
	} // 	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter )

	} catch ( utility::excn::EXCN_Base const & e ) {
		return -1;
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // int main( int argc, char * argv [] )
