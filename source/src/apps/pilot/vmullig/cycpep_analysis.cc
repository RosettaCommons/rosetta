// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/cycpep_analysis.cc
/// @brief Simple analysis of N-to-C cyclic peptides.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)


//General includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <numeric/angle.functions.hh>

#include <utility/vector1.hh>
#include <stdio.h>

//Tracer:
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.vmullig.cycpep_analysis" );

/*************************/
/* FUNCTION PROTOTYPES:  */
/*************************/

void ensure_cyclic( core::pose::PoseOP pose, core::Size const count );

bool has_cis_bonds( core::pose::PoseCOP pose );

void count_bins(
	core::pose::PoseCOP pose,
	core::Size &Acount,
	core::Size &Xcount,
	core::Size &Bcount,
	core::Size &Ycount,
	core::Size &Ocount,
	core::Size &Zcount
);

void do_cycpep_analysis();

/*************************/
/* FUNCTION DEFINITIONS: */
/*************************/

void
ensure_cyclic(
	core::pose::PoseOP pose,
	core::Size const count
) {
	core::Size const nres( pose->total_residue() );
	if( pose->residue(1).has_lower_connect() && pose->residue(nres).has_upper_connect() && pose->residue(1).is_polymer_bonded(nres) ) {
		TR << "Pose " << count << " is already cyclic." << std::endl;
	} else {
		protocols::cyclic_peptide::DeclareBond cyclize;
		cyclize.set( 1, "N", nres, "C", false );
		cyclize.apply( *pose );
		TR << "Cyclizing pose " << count << "." << std::endl;
	}
}

bool
has_cis_bonds(
	core::pose::PoseCOP pose
) {
	for(core::Size i=1, imax=pose->total_residue(); i<=imax; ++i) {
		core::Real const omval( numeric::principal_angle_degrees( pose->omega(i) ) );
		if( omval > -90 && omval <= 90 ) return true;
	}
	return false;
}

void
count_bins(
	core::pose::PoseCOP pose,
	core::Size &Acount,
	core::Size &Xcount,
	core::Size &Bcount,
	core::Size &Ycount,
	core::Size &Ocount,
	core::Size &Zcount
) {
	for(core::Size i=1, imax=pose->total_residue(); i<=imax; ++i) {
		core::Real const phival( numeric::principal_angle_degrees( pose->phi(i) ) );
		core::Real const psival( numeric::principal_angle_degrees( pose->psi(i) ) );
		core::Real const omegaval( numeric::principal_angle_degrees( pose->omega(i) ) );
		
		if( omegaval > -90 && omegaval <= 90 ) { //Cis peptide bond.
			if( phival <= 0 ) { ++Ocount; }
			else { ++Zcount; }
		} else { //Trans peptide bond.
			if( phival <= 0 ) { //Negative phi
				if( psival <= -125 || psival > 50 ) {
					++Bcount;
				} else {
					++Acount;
				}
			} else { //Positive phi
				if( psival <= -50 || psival > 125 ) {
					++Ycount;
				} else {
					++Xcount;
				}
			}
		}
		
	}
}

void
do_cycpep_analysis()
{
	core::import_pose::pose_stream::MetaPoseInputStream input( core::import_pose::pose_stream::streams_from_cmd_line() );
	core::Size total_count(0), valid_count(0), Acount(0), Xcount(0), Bcount(0), Ycount(0), Ocount(0), Zcount(0);
	
	while( input.has_another_pose() ) {
		++total_count;

		core::pose::PoseOP pose( new core::pose::Pose ); //Create storage for the pose.
		input.fill_pose( *pose ); //Import the pose
		ensure_cyclic( pose, total_count );
		
		//Skip poses with cis peptide bonds
		if( has_cis_bonds( pose ) ) {
			TR << "Pose " << total_count << " contains cis bonds.  Skipping." << std::endl;
			continue;
		}
		
		++valid_count;
		
		count_bins( pose, Acount, Xcount, Bcount, Ycount, Ocount, Zcount );
	}
	
	TR << "\nSUMMARY:\nA:\t" << Acount << "\nB:\t" << Bcount << "\nX:\t" << Xcount << "\nY:\t" << Ycount << "\nO:\t" << Ocount << "\nZ:\t" << Zcount << "\nTotal poses processed:\t" << total_count << "\nTotal poses valid:\t" << valid_count << "\n" << std::endl;
}


/*************************/
/* MAIN: */
/*************************/
int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);
		TR << "Starting cycpep_analysis.cc.\nPilot app created by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory, on 13 October 2016." << std::endl;
		do_cycpep_analysis();
	} catch ( utility::excn::EXCN_Base& excn ) {
		TR.Error << "Exception caught: " << std::endl;
		excn.show( TR.Error );
		TR.Error.flush();
		return -1;
	}

	if ( TR.visible() ) {
		TR << "Finished cycpep_analysis.cc.  Exiting." << std::endl;
		TR.flush();
	}

	return 0;
}
