// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/id/AtomID.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/build_pose_as_is.hh>

#include <core/pose/Pose.hh>

#include <utility/string_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


OPT_1GRP_KEY( String, measure, D1 )
OPT_1GRP_KEY( String, measure, D2 )
OPT_1GRP_KEY( String, measure, D3 )

OPT_1GRP_KEY( String, measure, U1 )
OPT_1GRP_KEY( String, measure, U2 )
OPT_1GRP_KEY( String, measure, U3 )

OPT_1GRP_KEY( Integer, measure, upstream_res )
OPT_1GRP_KEY( Integer, measure, downstream_res )


int main( int argc, char * argv [] )
{
	try {

		using namespace core;
		using namespace core::chemical;
		using namespace core::id;
		using namespace core::io::pdb;
		using namespace core::pose;

		using namespace basic::options;
		using namespace basic::options::OptionKeys::measure;

		NEW_OPT( measure::D1, "d1 atom name", "" );
		NEW_OPT( measure::D2, "d2 atom name", "" );
		NEW_OPT( measure::D3, "d3 atom name", "" );
		NEW_OPT( measure::U1, "u1 atom name", "" );
		NEW_OPT( measure::U2, "u2 atom name", "" );
		NEW_OPT( measure::U3, "u3 atom name", "" );
		NEW_OPT( measure::upstream_res,   "upstream residue #", 1 );
		NEW_OPT( measure::downstream_res, "downstream residue #", 1 );

		devel::init( argc, argv );

		utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

		if ( input_jobs.size() != 1 ) {
			utility_exit_with_message( "Expected exactly one pdb to be specified from the -s or -l flags" );
		}

		pose::Pose pose;
		core::import_pose::pose_from_file( pose, input_jobs[ 1 ]->input_tag() , core::import_pose::PDB_file);


		Size upres, downres;
		AtomID u1, u2, u3, d1, d2, d3;

		if ( option[ upstream_res ] <= 0 ) {
			utility_exit_with_message( "Negative upstream resid is illegal: read " + utility::to_string( option[ upstream_res ]() ) + " for measure::upstream_res"  );
		} else if ( Size( option[ upstream_res ] ) > pose.total_residue() ) {
			utility_exit_with_message( "Upstream resid exceeds the number of residues in the input pose: read " + utility::to_string( option[ upstream_res ]() ) + " for measure::upstream_res" );
		} else {
			upres = option[ upstream_res ];
		}

		if ( option[ downstream_res ] <= 0 ) {
			utility_exit_with_message( "Negative downstream resid is illegal: read " + utility::to_string( option[ downstream_res ]() ) + " for measure::downstream_res"  );
		} else if ( Size( option[ downstream_res ] ) > pose.total_residue() ) {
			utility_exit_with_message( "Upstream resid exceeds the number of residues in the input pose: read " + utility::to_string( option[ downstream_res ]() ) + " for measure::downstream_res" );
		} else {
			downres = option[ downstream_res ];
		}

		if ( ! pose.residue( upres ).has( option[ U1 ] ) ) {
			utility_exit_with_message( "Upstream residue does not contain an atom named '" + option[ U1 ]() + "', resid " + utility::to_string( upres ) + " is " + pose.residue( upres ).name() );
		} else {
			u1 = AtomID( pose.residue( upres ).atom_index( option[ U1 ] ), upres );
		}

		if ( ! pose.residue( upres ).has( option[ U2 ] ) ) {
			utility_exit_with_message( "Upstream residue does not contain an atom named '" + option[ U2 ]() + "', resid " + utility::to_string( upres ) + " is " + pose.residue( upres ).name() );
		} else {
			u2 = AtomID( pose.residue( upres ).atom_index( option[ U2 ] ), upres );
		}

		if ( ! pose.residue( upres ).has( option[ U3 ] ) ) {
			utility_exit_with_message( "Upstream residue does not contain an atom named '" + option[ U3 ]() + "', resid " + utility::to_string( upres ) + " is " + pose.residue( upres ).name() );
		} else {
			u3 = AtomID( pose.residue( upres ).atom_index( option[ U3 ] ), upres );
		}

		if ( ! pose.residue( downres ).has( option[ D1 ] ) ) {
			utility_exit_with_message( "Dpstream residue does not contain an atom named '" + option[ D1 ]() + "', resid " + utility::to_string( downres ) + " is " + pose.residue( downres ).name() );
		} else {
			d1 = AtomID( pose.residue( downres ).atom_index( option[ D1 ] ), downres );
		}

		if ( ! pose.residue( downres ).has( option[ D2 ] ) ) {
			utility_exit_with_message( "Dpstream residue does not contain an atom named '" + option[ D2 ]() + "', resid " + utility::to_string( downres ) + " is " + pose.residue( downres ).name() );
		} else {
			d2 = AtomID( pose.residue( downres ).atom_index( option[ D2 ] ), downres );
		}

		if ( ! pose.residue( downres ).has( option[ D3 ] ) ) {
			utility_exit_with_message( "Dpstream residue does not contain an atom named '" + option[ D3 ]() + "', resid " + utility::to_string( downres ) + " is " + pose.residue( downres ).name() );
		} else {
			d3 = AtomID( pose.residue( downres ).atom_index( option[ D3 ] ), downres );
		}


		std::cout << "Upstream chi: ";
		for ( Size ii = 1; ii <= pose.residue( upres ).nchi(); ++ii ) {
			std::cout << pose.residue( upres ).chi( ii ) << " ";
		}
		std::cout << std::endl;

		Vector vu1( pose.xyz( u1 ));
		Vector vu2( pose.xyz( u2 ));
		Vector vu3( pose.xyz( u3 ));
		Vector vd1( pose.xyz( d1 ));
		Vector vd2( pose.xyz( d2 ));
		Vector vd3( pose.xyz( d3 ));

		std::cout << "TOR U3D1 " << numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( vu3, vu2, vu1, vd1 ) << std::endl;
		std::cout << "TOR U2D2 " << numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( vu2, vu1, vd1, vd2 ) << std::endl;
		std::cout << "TOR U1D3 " << numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( vu1, vd1, vd2, vd3 ) << std::endl;
		std::cout << "ANG U2D1 " << numeric::constants::d::radians_to_degrees * numeric::angle_radians( vu2, vu1, vd1 ) << std::endl;
		std::cout << "ANG U1D2 " << numeric::constants::d::radians_to_degrees * numeric::angle_radians( vu1, vd1, vd2 ) << std::endl;
		std::cout << "DIS U1D1 " <<  vu1.distance( vd1 ) << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

