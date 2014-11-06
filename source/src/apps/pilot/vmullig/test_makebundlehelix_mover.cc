// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file test_makebundlehelix_mover.cc
/// @brief A test of the mover that makes a helix in a helical bundle.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>

//Application-specific includes
#include <protocols/helical_bundle/MakeBundleHelix.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Tracer:
static basic::Tracer TR( "apps.pilot.vmullig.test_makebundlehelix_mover" );

//Options (ugh -- global variables):
OPT_KEY (Real, r0)
OPT_KEY (Real, omega0)
OPT_KEY (Real, delta_omega0)
OPT_KEY (Integer, residue_repeats)
OPT_KEY (String, residue_type)
OPT_KEY (String, tail_residue_type)
OPT_KEY (RealVector, r1)
OPT_KEY (Real, omega1)
OPT_KEY (Real, z1)
OPT_KEY (Real, delta_omega1)
OPT_KEY (RealVector, delta_omega1_per_atom)
OPT_KEY (RealVector, delta_z1_per_atom)
OPT_KEY (Boolean, invert_helix)
OPT_KEY (Real, delta_t)

///
/// @brief Set up the options for this pilot app.
void register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1 < core::Real > alpha_helix_r1;
	alpha_helix_r1.push_back( 1.5243286 );
	alpha_helix_r1.push_back( 2.2819007 );
	alpha_helix_r1.push_back( 1.7156595 );

	utility::vector1 < core::Real > alpha_helix_delta_omega1;
	alpha_helix_delta_omega1.push_back( -0.46047311 );
	alpha_helix_delta_omega1.push_back( 0.0 );
	alpha_helix_delta_omega1.push_back( 0.47947732 );

	utility::vector1 < core::Real > alpha_helix_delta_z1;
	alpha_helix_delta_z1.push_back( -0.91007351 );
	alpha_helix_delta_z1.push_back( 0.0 );
	alpha_helix_delta_z1.push_back( 1.0570712 );

	NEW_OPT( r0, "The radius of the major helix, in Angstroms.  Default 8.0.", 8.0 );
	NEW_OPT( omega0, "The turn per residue of the major helix, in radians.  Default 0.2.", 0.2 );
	NEW_OPT( delta_omega0, "The offset of the major helix, in radians.  Default 0.0.", 0.0 );
	NEW_OPT( residue_repeats, "How many residues are in the helix.  Default 20", 20 );
	NEW_OPT( residue_type, "The residue type from which the helix will be constructed.  (This needs to be the full Rosetta name, not just the 3-letter code).  Default \"ALA\"", "ALA" );
	NEW_OPT( tail_residue_type, "The residue type that will cap the helix.  (This needs to be the full Rosetta name, not just the 3-letter code).  No caps if not specified.", "" );
	NEW_OPT( r1, "The r1 values for the minor helix.  One value must be specified for each mainchain atom in the residue type used.  Default paramaters are for an alpha helix.", alpha_helix_r1 );
	NEW_OPT( omega1, "The omega1 (turn per residue, in radians) value for the minor helix.  Default parameter is for an alpha helix.", 1.7277092);
	NEW_OPT( z1, "The z1 value (rise per residue, in Angstroms) for the minor helix.  Default parameter is for an alpha helix.", 1.5513078);
	NEW_OPT( delta_omega1, "The offset of the minor helix, in radians.  Default 0.0.", 0.0 );
	NEW_OPT( delta_omega1_per_atom, "The per-atom omega1 offset of the minor helix, in radians.  Default parameters are for an alpha helix.", alpha_helix_delta_omega1 );
	NEW_OPT( delta_z1_per_atom, "The per-atom z1 offset of the minor helix, in Angstroms.  Default parameters are for an alpha helix.", alpha_helix_delta_z1 );
	NEW_OPT( invert_helix, "If this flag is added, the helix runs in the opposite direction (but with the same chirality).  Not inverted by default.", false);
	NEW_OPT( delta_t, "An offset for the value of t (the residue index), used to shift the register of the helix up or down.  Default 0.", 0.0);

	return;
}


int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try {
		register_options();
		devel::init(argc, argv);

		if(TR.visible()) {
			TR << "Starting test_makebundlehelix_mover.cc" << std::endl;
			TR << "Pilot app created 28 October 2014 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
		}

		core::pose::Pose pose; //Make the empty pose

		protocols::helical_bundle::MakeBundleHelix makehelix;

		makehelix.set_reset_pose(true);
		makehelix.set_helix_length( static_cast<core::Size>(option[residue_repeats]()) );
		makehelix.set_residue_name( option[residue_type]() );
		makehelix.set_tail_residue_name( option[tail_residue_type]() );
		makehelix.set_major_helix_params ( option[r0](), option[omega0](), option[delta_omega0]());
		makehelix.set_invert_helix( option[invert_helix].user() );

		utility::vector1 < core::Real > r1vals = option[r1]();

		utility::vector1 < core::Real > delta_omega1vals = option[delta_omega1_per_atom]();
		for(core::Size i=1, imax=delta_omega1vals.size(); i<=imax; ++i) delta_omega1vals[i] += option[delta_omega1]();

		utility::vector1 < core::Real > delta_z1vals = option[delta_z1_per_atom]();

		makehelix.set_delta_t( option[delta_t]() );

		makehelix.set_minor_helix_params( r1vals, option[omega1](), option[z1](), delta_omega1vals, delta_z1vals );

		makehelix.apply(pose);

		pose.dump_pdb("output.pdb");

		if(TR.visible()) {
			TR << "Finished test_makebundlehelix_mover.cc.  Exiting." << std::endl;
			TR.flush();
		}

	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}
