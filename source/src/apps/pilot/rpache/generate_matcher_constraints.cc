// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief produces matcher constraints for an upstream residue against a downstream residue
/// @author Andrew Leaver-Fay
/// @author Daniel J. Mandell
/// @author Roland A. Pache, Ph.D.

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <core/id/AtomID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pose_io.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/file_data.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/string_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/constants.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


OPT_1GRP_KEY( String, measure, D1 )
OPT_1GRP_KEY( String, measure, D2 )
OPT_1GRP_KEY( String, measure, D3 )

OPT_1GRP_KEY( String, measure, U1 )
OPT_1GRP_KEY( String, measure, U2 )
OPT_1GRP_KEY( String, measure, U3 )

OPT_1GRP_KEY( Boolean, measure, T_OH )
OPT_1GRP_KEY( Boolean, measure, Y_aroC )

OPT_1GRP_KEY( Integer, measure, upstream_res )
OPT_1GRP_KEY( Integer, measure, downstream_res )

static THREAD_LOCAL basic::Tracer TR( "generate_matcher_constraints" );

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
		NEW_OPT( measure::T_OH, "use hydroxyl oxygen as main reference atom for THR", "" );
		NEW_OPT( measure::Y_aroC, "use aromatic carbon as main reference atom for TYR", "" );
		NEW_OPT( measure::upstream_res,   "upstream residue #", 1 );
		NEW_OPT( measure::downstream_res, "downstream residue #", 1 );

		devel::init( argc, argv );

		utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

		if ( input_jobs.size() != 1 ) {
			utility_exit_with_message( "Expected exactly one pdb to be specified from the -s or -l flags" );
		}

		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, input_jobs[ 1 ]->input_tag() );


		Size upres, downres;
		AtomID u1, u2, u3, d1, d2, d3;

		if ( option[ upstream_res ] <= 0 ) {
			utility_exit_with_message( "Negative upstream resid is illegal: read " + utility::to_string( option[ upstream_res ]() ) + " for measure::upstream_res"  );
		} else if ( ( unsigned int ) option[ upstream_res ] > pose.total_residue() ) {
			utility_exit_with_message( "Upstream resid exceeds the number of residues in the input pose: read " + utility::to_string( option[ upstream_res ]() ) + " for measure::upstream_res" );
		} else {
			upres = option[ upstream_res ];
		}

		if ( option[ downstream_res ] <= 0 ) {
			utility_exit_with_message( "Negative downstream resid is illegal: read " + utility::to_string( option[ downstream_res ]() ) + " for measure::downstream_res"  );
		} else if ( ( unsigned int ) option[ downstream_res ] > pose.total_residue() ) {
			utility_exit_with_message( "Upstream resid exceeds the number of residues in the input pose: read " + utility::to_string( option[ downstream_res ]() ) + " for measure::downstream_res" );
		} else {
			downres = option[ downstream_res ];
		}

		if ( ! ( option[ U1 ].user() && option[ U2 ].user() && option[ U3 ].user()) ) { // auto-select the upstream atoms
			TR << "All upstream atoms not specified. Auto-selecting from residue " << upres << std::endl;
			//utility::vector1<core::Size> testvec;
			//testvec = protocols::enzdes::EnzCstTemplateRes::atom_inds_for_restype(1, &(pose.residue( upres ).type()));
			//for (int zz=1; zz <= testvec.size(); zz++) {
			// TR << zz << std::endl;
			//}
			//TR << pose.residue( upres ).atom_name(pose.residue( upres ).atom_end()->type()) << std::endl;

			//ResidueType uptype = pose.residue( upres ).type();

			core::Size u1ind = pose.residue( upres ).nheavyatoms();
			// SPECIAL CASES
			if ( pose.residue(upres).name3()=="THR" && option[T_OH].user() ) {
				// Go back one atom to get the OH for Thr
				--u1ind;
			} else if ( pose.residue(upres).name3()=="TYR" && option[Y_aroC].user() ) {
				// Go back one atom to get an aroC for Tyr
				--u1ind;
			}

			std::string u1name = pose.residue( upres ).atom_name( u1ind );
			u1 = AtomID( pose.residue( upres ).atom_index( u1name ), upres);

			core::Size u2ind = pose.residue( upres ).atom_base( u1ind );
			std::string u2name = pose.residue( upres ).atom_name( u2ind );
			u2 = AtomID( pose.residue( upres ).atom_index( u2name ), upres );

			core::Size u3ind = pose.residue( upres ).atom_base( u2ind );
			std::string u3name = pose.residue( upres ).atom_name( u3ind );
			u3 = AtomID( pose.residue( upres ).atom_index( u3name ), upres );
		} else { // else get the upstream atoms from the command line
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
		}

		TR << "U1: " << pose.residue( upres ).atom_name( u1.atomno() ) << ", U2: " << pose.residue( upres ).atom_name( u2.atomno() ) << ", U3: " << pose.residue( upres ).atom_name( u3.atomno() ) << std::endl;
		std::string u1type = (pose.residue( upres ).atom_type( u1.atomno() )).name();
		TR << "U1 type is " << u1type << std::endl;

		if ( ! ( option[ D1 ].user() && option[ D2 ].user() && option[ D3 ].user()) ) { // auto-seslect the downstream
			TR << "All downstream atoms not specified. Auto-selecting from residue " << downres << std::endl;
			// Use the 3 middle-most heavy atoms
			/*
			Size d1ind = (pose.residue( downres ).nheavyatoms() / 2) - 1;
			Size d2ind = d1ind + 1;
			Size d3ind = d2ind + 1;
			std::string d1name = pose.residue( downres ).atom_name( d1ind );
			std::string d2name = pose.residue( downres ).atom_name( d2ind );
			std::string d3name = pose.residue( downres ).atom_name( d3ind );
			d1 = AtomID( d1ind, downres );
			d2 = AtomID( d2ind, downres );
			d3 = AtomID( d3ind, downres );
			*/
			/*
			//// Use the 3 closest atoms to the downstream residue's cartesian centroid
			// Find the centroid
			Vector sum ( 0.0 );
			for ( Size i=1; i<=pose.residue( downres ).nheavyatoms(); ++i ) {
			AtomID curid = AtomID( i, downres );
			sum += pose.xyz( curid );
			}
			Vector centroid = sum / pose.residue( downres ).nheavyatoms();

			// Find the 3 downstream atoms closest to the centroid
			Real low1 = 99999., low2 = 99999., low3 = 99999.;
			for ( Size i=1; i<=pose.residue( downres ).nheavyatoms(); ++i ) {
			AtomID curid = AtomID( i, downres );
			//Real curdist = upavg.distance( pose.xyz( curid ) );
			Real curdist = centroid.distance( pose.xyz( curid ) );
			if ( curdist < low1 ) {
			low3 = low2;
			low2 = low1;
			low1 = curdist;
			d3 = d2;
			d2 = d1;
			d1 = curid;
			}
			else if ( curdist < low2 ) {
			low3 = low2;
			low2 = curdist;
			d3 = d2;
			d2 = curid;
			}
			else if ( curdist < low3 ) {
			low3 = curdist;
			d3 = curid;
			}
			}
			core::Size d1ind = d1.atomno();
			std::string d1name = pose.residue( downres ).atom_name( d1ind );

			core::Size d2ind = d2.atomno();
			std::string d2name = pose.residue( downres ).atom_name( d2ind );

			core::Size d3ind = d3.atomno();
			std::string d3name = pose.residue( downres ).atom_name( d3ind );
			*/


			// find the 3 closest atoms to the upstream residue atoms
			//Vector upavg = ( pose.xyz( u1 ) + pose.xyz( u2 ) + pose.xyz( u3 ) ) / 3;

			// DJM: 1/13/10 -- now taking 3 closest Dstream atoms to U1
			Vector u1pos = pose.xyz( u1 );

			Real low1 = 99999., low2 = 99999., low3 = 99999.;
			for ( Size i=1; i<=pose.residue( downres ).nheavyatoms(); ++i ) {
				AtomID curid = AtomID( i, downres );
				//Real curdist = upavg.distance( pose.xyz( curid ) );
				Real curdist = u1pos.distance( pose.xyz( curid ) );
				if ( curdist < low1 ) {
					low3 = low2;
					low2 = low1;
					low1 = curdist;
					d3 = d2;
					d2 = d1;
					d1 = curid;
				} else if ( curdist < low2 ) {
					low3 = low2;
					low2 = curdist;
					d3 = d2;
					d2 = curid;
				} else if ( curdist < low3 ) {
					low3 = curdist;
					d3 = curid;
				}
			}
			core::Size d1ind = d1.atomno();
			std::string d1name = pose.residue( downres ).atom_name( d1ind );

			core::Size d2ind = d2.atomno();
			std::string d2name = pose.residue( downres ).atom_name( d2ind );

			core::Size d3ind = d3.atomno();
			std::string d3name = pose.residue( downres ).atom_name( d3ind );

		} else { // else get the downstream atoms from the command line
			if ( ! pose.residue( downres ).has( option[ D1 ] ) ) {
				utility_exit_with_message( "Downstream residue does not contain an atom named '" + option[ D1 ]() + "', resid " + utility::to_string( downres ) + " is " + pose.residue( downres ).name() );
			} else {
				d1 = AtomID( pose.residue( downres ).atom_index( option[ D1 ] ), downres );
			}

			if ( ! pose.residue( downres ).has( option[ D2 ] ) ) {
				utility_exit_with_message( "Downstream residue does not contain an atom named '" + option[ D2 ]() + "', resid " + utility::to_string( downres ) + " is " + pose.residue( downres ).name() );
			} else {
				d2 = AtomID( pose.residue( downres ).atom_index( option[ D2 ] ), downres );
			}

			if ( ! pose.residue( downres ).has( option[ D3 ] ) ) {
				utility_exit_with_message( "Downstream residue does not contain an atom named '" + option[ D3 ]() + "', resid " + utility::to_string( downres ) + " is " + pose.residue( downres ).name() );
			} else {
				d3 = AtomID( pose.residue( downres ).atom_index( option[ D3 ] ), downres );
			}
		}

		TR << "D1: " << pose.residue( downres ).atom_name( d1.atomno() ) << ", D2: " << pose.residue( downres ).atom_name( d2.atomno() ) << ", D3: " << pose.residue( downres ).atom_name( d3.atomno() ) << std::endl;
		std::string d1type = (pose.residue( downres ).atom_type( d1.atomno() )).name();
		TR << "D1 type is " << d1type << std::endl;

		TR << "Upstream chi: ";
		for ( Size ii = 1; ii <= pose.residue( upres ).nchi(); ++ii ) {
			TR << pose.residue( upres ).chi( ii ) << " ";
		}
		TR << std::endl;

		Vector vu1( pose.xyz( u1 ));
		Vector vu2( pose.xyz( u2 ));
		Vector vu3( pose.xyz( u3 ));
		Vector vd1( pose.xyz( d1 ));
		Vector vd2( pose.xyz( d2 ));
		Vector vd3( pose.xyz( d3 ));

		TR << "UPRES_NAME1 " << pose.residue( upres ).name1() << std::endl;
		TR << "UPRES_NAME3 " << pose.residue( upres ).name3() << std::endl;
		TR << "UPRES_CHAIN " << pose.pdb_info()->chain( upres ) << std::endl;
		TR << "UPRES_ATOMTYPE " << u1type << std::endl;
		TR << "DOWNRES_NAME3 " << pose.residue( downres ).name3() << std::endl;
		TR << "DOWNRES_ATOMS " << pose.residue( downres ).atom_name( d1.atomno() )
			<< " " << pose.residue( downres ).atom_name( d2.atomno() )
			<< " " << pose.residue( downres ).atom_name( d3.atomno() ) << std::endl;

		TR << "distanceAB " <<  vu1.distance( vd1 ) << std::endl;
		TR << "angle_A " << numeric::constants::d::radians_to_degrees * numeric::angle_radians( vu1, vd1, vd2 ) << std::endl;
		TR << "angle_B " << numeric::constants::d::radians_to_degrees * numeric::angle_radians( vu2, vu1, vd1 ) << std::endl;
		TR << "torsion_A " << numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( vu1, vd1, vd2, vd3 ) << std::endl;
		TR << "torsion_AB " << numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( vu2, vu1, vd1, vd2 ) << std::endl;
		TR << "torsion_B " << numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( vu3, vu2, vu1, vd1 ) << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
