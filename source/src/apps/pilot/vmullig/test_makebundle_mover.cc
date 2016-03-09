// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file test_makebundle_mover.cc
/// @brief A test of the mover that makes a helix in a helical bundle.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

//C++ includes
#include <stdio.h>

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
#include <numeric/constants.hh>

//Application-specific includes
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/MakeBundleHelix.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Tracer:
static basic::Tracer TR( "apps.pilot.vmullig.test_makebundle_mover" );

//Options (ugh -- global variables):
OPT_KEY (Integer, animation_frames)
OPT_KEY (Boolean, animation_cyclic)
OPT_KEY (Integer, radial_symmetry)
OPT_KEY (Integer, helix_count)
OPT_KEY (RealVector, r0)
OPT_KEY (RealVector, r0_end)
OPT_KEY (RealVector, omega0)
OPT_KEY (RealVector, omega0_end)
OPT_KEY (RealVector, delta_omega0)
OPT_KEY (RealVector, delta_omega0_end)
OPT_KEY (IntegerVector, residue_repeats)
OPT_KEY (String, residue_type)
OPT_KEY (String, tail_residue_type)
OPT_KEY (RealVector, r1)
OPT_KEY (Real, omega1)
OPT_KEY (Real, z1)
OPT_KEY (RealVector, delta_omega1)
OPT_KEY (RealVector, delta_omega1_end)
OPT_KEY (RealVector, delta_omega1_per_atom)
OPT_KEY (RealVector, delta_z1_per_atom)
OPT_KEY (BooleanVector, invert_helix)
OPT_KEY (RealVector, delta_t)
OPT_KEY (RealVector, delta_t_end)
OPT_KEY (RealVector, z0_offset)
OPT_KEY (RealVector, z0_offset_end)
OPT_KEY (RealVector, z1_offset)
OPT_KEY (RealVector, z1_offset_end)
OPT_KEY (String, minor_helix_params)
OPT_KEY (Boolean, set_bondlengths)
OPT_KEY (Boolean, set_bondangles)
OPT_KEY (Boolean, set_dihedrals)


/// @brief Set up the options for this pilot app.
void register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1 <core::Real> empty_vector;

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

	utility::vector1 < core::Size > repeats_default; repeats_default.push_back(20);
	utility::vector1 < core::Real > r0_default; r0_default.push_back(8.0);
	utility::vector1 < core::Real > omega0_default; omega0_default.push_back(-0.05);
	utility::vector1 < core::Real > delta_omega0_default; delta_omega0_default.push_back(0.0);
	utility::vector1 < core::Real > delta_omega1_default; delta_omega1_default.push_back(0.0);
	utility::vector1 < bool > invert_default; invert_default.push_back(false);
	utility::vector1 < core::Real > delta_t_default; delta_t_default.push_back(0.0);
	utility::vector1 < core::Real > z0_offset_default; z0_offset_default.push_back(0.0);
	utility::vector1 < core::Real > z1_offset_default; z1_offset_default.push_back(0.0);

	NEW_OPT( animation_frames, "The number of frames of animation to produce.  If not specified (or set to 0 or 1), then no animation is generated.", 0 );
	NEW_OPT( animation_cyclic, "Should the animation loop (by rocking sinusoidally between the start and end states), or just be linear?  Default false (linear).", false);
	NEW_OPT( radial_symmetry, "The radial symmetry of the bundle. Values of 0 or 1 mean no symmetry, 2 means 2-fold symetry, 3 means 3-fold symmetry, etc.  Default 3.", 3 );
	NEW_OPT( helix_count, "The number of helices that will be defined.  The final helix count is this value multiplied by the symmetry.  Default 1", 1 );
	NEW_OPT( r0, "The radius of the major helix, in Angstroms.  Default 8.0.  One value for each helix.", r0_default);
	NEW_OPT( r0_end, "The radius of the major helix, in Angstroms, at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( omega0, "The turn per residue of the major helix, in radians.  Default -0.05.  One value for each helix.", omega0_default );
	NEW_OPT( omega0_end, "The turn per residue of the major helix, in radians, at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( delta_omega0, "The offset of the major helix, in radians.  Default 0.0.  One value for each helix.",  delta_omega0_default );
	NEW_OPT( delta_omega0_end, "The offset of the major helix, in radians, at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( residue_repeats, "How many residues are in the helix.  Default 20.  One value for each helix.", repeats_default );
	NEW_OPT( residue_type, "The residue type from which the helix will be constructed.  (This needs to be the full Rosetta name, not just the 3-letter code).  Default \"ALA\"", "ALA" );
	NEW_OPT( tail_residue_type, "The residue type that will cap the helix.  (This needs to be the full Rosetta name, not just the 3-letter code).  No caps if not specified.", "" );
	NEW_OPT( r1, "The r1 values for the minor helix.  One value must be specified for each mainchain atom in the residue type used.  Default paramaters are for an alpha helix.", alpha_helix_r1 );
	NEW_OPT( omega1, "The omega1 (turn per residue, in radians) value for the minor helix.  Default parameter is for an alpha helix.", 1.7277092);
	NEW_OPT( z1, "The z1 value (rise per residue, in Angstroms) for the minor helix.  Default parameter is for an alpha helix.", 1.5513078);
	NEW_OPT( delta_omega1, "The offset of the minor helix, in radians.  Default 0.0.  One value for each helix.", delta_omega1_default );
	NEW_OPT( delta_omega1_end, "The offset of the minor helix, in radians, at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( delta_omega1_per_atom, "The per-atom omega1 offset of the minor helix, in radians.  Default parameters are for an alpha helix.", alpha_helix_delta_omega1 );
	NEW_OPT( delta_z1_per_atom, "The per-atom z1 offset of the minor helix, in Angstroms.  Default parameters are for an alpha helix.", alpha_helix_delta_z1 );
	NEW_OPT( invert_helix, "If this flag is added, the helix runs in the opposite direction (but with the same chirality).  Not inverted by default.  One value for each helix.", invert_default);
	NEW_OPT( delta_t, "An offset for the value of t (the residue index), used to shift the register of the helix up or down.  Default 0.  One value for each helix.", delta_t_default);
	NEW_OPT( delta_t_end, "The offset for the value of t (the residue index) at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( z0_offset, "An offset in the z-direction along the major helix axis.  One value per helix.  Default 0.", z0_offset_default);
	NEW_OPT( z0_offset_end, "The offset in the z-direction along the major helix axis at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( z1_offset, "An offset in the z-direction along the minor helix axis.  One value per helix.  Default 0.", z1_offset_default);
	NEW_OPT( z1_offset_end, "The offset in the z-direction along the minor helix axis at the end of the animation.  Default unused.  One value for each helix.  If not specified, this property is not animated.", empty_vector);
	NEW_OPT( minor_helix_params, "The .crick_params file specifying the minor helix params for the type of backbone and secondary structure that one wants to use.  If this option is used, all minor helix params set with other flags are ignored.  Default unused.", "");
	NEW_OPT( set_bondlengths, "Should the generator set mainchain bond lengths?  If true, non-ideal geometry results, but a more perfect helix of helices is possible.  Default false.", false);
	NEW_OPT( set_bondangles, "Should the generator set mainchain bond angles?  If true, non-ideal geometry results, but a more perfect helix of helices is possible.  Default false.", false);
	NEW_OPT( set_dihedrals, "Should the generator set mainchain dihedrals?  Default true.", true);

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
			TR << "Starting test_makebundle_mover.cc" << std::endl;
			TR << "Pilot app created 30 October 2014 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
		}

		core::Size const helixcount = static_cast<core::Size>(option[helix_count]());

		runtime_assert_string_msg( option[residue_repeats]().size() >= helixcount, "A value must be provided with the -residue_repeats flag for EACH helix." );
		runtime_assert_string_msg( option[r0]().size() >= helixcount, "A value must be provided with the -r0 flag for EACH helix." );
		if(option[r0_end].user()) runtime_assert_string_msg( option[r0_end]().size() >= helixcount, "A value must be provided with the -r0_end flag for EACH helix if this flag is used." );
		if(option[delta_omega0_end].user()) runtime_assert_string_msg( option[delta_omega0_end]().size() >= helixcount, "A value must be provided with the -delta_omega0_end flag for EACH helix if this flag is used." );
		if(option[delta_omega1_end].user()) runtime_assert_string_msg( option[delta_omega1_end]().size() >= helixcount, "A value must be provided with the -delta_omega1_end flag for EACH helix if this flag is used." );
		if(option[delta_t_end].user()) runtime_assert_string_msg( option[delta_t_end]().size() >= helixcount, "A value must be provided with the -delta_t_end flag for EACH helix if this flag is used." );
		if(option[z0_offset_end].user()) runtime_assert_string_msg( option[z0_offset_end]().size() >= helixcount, "A value must be provided with the -z0_offset_end flag for EACH helix if this flag is used." );
		if(option[z1_offset_end].user()) runtime_assert_string_msg( option[z1_offset_end]().size() >= helixcount, "A value must be provided with the -z1_offset_end flag for EACH helix if this flag is used." );
		runtime_assert_string_msg( option[omega0]().size() >= helixcount, "A value must be provided with the -omega0 flag for EACH helix." );
		if(option[omega0_end].user()) runtime_assert_string_msg( option[omega0_end]().size() >= helixcount, "A value must be provided with the -omega0_end flag for EACH helix if this flag is used." );
		runtime_assert_string_msg( option[delta_omega0]().size() >= helixcount, "A value must be provided with the -delta_omega0 flag for EACH helix." );
		runtime_assert_string_msg( option[invert_helix]().size() >= helixcount, "A value must be provided with the -invert_helix flag for EACH helix." );
		runtime_assert_string_msg( option[delta_t]().size() >= helixcount, "A value must be provided with the -delta_t flag for EACH helix." );
		runtime_assert_string_msg( option[delta_omega1]().size() >= helixcount, "A value must be provided with the -delta_omega1 flag for EACH helix." );
		runtime_assert_string_msg( option[z0_offset]().size() >= helixcount, "A value must be provided with the -z0_offset flag for EACH helix." );
		runtime_assert_string_msg( option[z1_offset]().size() >= helixcount, "A value must be provided with the -z1_offset flag for EACH helix." );
		if(option[animation_cyclic].user()) runtime_assert_string_msg( option[animation_frames]() >= 2, "If the -animation_cyclic flag is used, the animation must be at least 2 frames long." );

		core::Size const nframes = (option[animation_frames]() > 2 ? static_cast<core::Size>(option[animation_frames]()) : 1);
		core::Real cur_frame_time(0.0);

		for(core::Size iframe=1; iframe<=nframes; ++iframe) {
			cur_frame_time = static_cast<core::Real>(iframe-1) / static_cast<core::Real>(nframes-1); //Will run from 0 to 1
			if(option[animation_cyclic].user()) {
				cur_frame_time = static_cast<core::Real>(iframe-1) / static_cast<core::Real>(nframes); //Will run from 0 to 1-(1/nframes)
				cur_frame_time = -0.5*cos(cur_frame_time*numeric::constants::d::pi_2)+0.5; //Will sinusoidally oscillate from 0 to 1 and back to (almost) 0
			}

			core::pose::Pose pose; //Make the empty pose
			protocols::helical_bundle::MakeBundle makebundle;
			makebundle.set_symmetry( static_cast<core::Size>(option[radial_symmetry]()) );
			makebundle.set_default_allow_bondlengths( option[set_bondlengths]() );
			makebundle.set_default_allow_bondangles( option[set_bondangles]() );
			makebundle.set_default_allow_dihedrals( option[set_dihedrals]() );

			for(core::Size ihelix=1; ihelix<=helixcount; ++ihelix) {
				core::Real r0val = option[r0]()[ihelix];
				if(option[animation_frames]()>2 && option[r0_end].user()) r0val = (1.0-cur_frame_time)*option[r0]()[ihelix] + cur_frame_time*option[r0_end]()[ihelix]; //Linearly interpolate start and end values.
				core::Real omega0val = option[omega0]()[ihelix];
				if(option[animation_frames]()>2 && option[omega0_end].user()) omega0val = (1.0-cur_frame_time)*option[omega0]()[ihelix] + cur_frame_time*option[omega0_end]()[ihelix]; //Linearly interpolate start and end values.
				core::Real deltaomega0val = option[delta_omega0]()[ihelix];
				if(option[animation_frames]()>2 && option[delta_omega0_end].user()) deltaomega0val = (1.0-cur_frame_time)*option[delta_omega0]()[ihelix] + cur_frame_time*option[delta_omega0_end]()[ihelix]; //Linearly interpolate start and end values.
				core::Real deltaomega1val = option[delta_omega1]()[ihelix];
				if(option[animation_frames]()>2 && option[delta_omega1_end].user()) deltaomega1val = (1.0-cur_frame_time)*option[delta_omega1]()[ihelix] + cur_frame_time*option[delta_omega1_end]()[ihelix]; //Linearly interpolate start and end values.
				core::Real deltatval = option[delta_t]()[ihelix];
				if(option[animation_frames]()>2 && option[delta_t_end].user()) deltatval = (1.0-cur_frame_time)*option[delta_t]()[ihelix] + cur_frame_time*option[delta_t_end]()[ihelix]; //Linearly interpolate start and end values.
				core::Real z0offsetval = option[z0_offset]()[ihelix];
				if(option[animation_frames]()>2 && option[z0_offset_end].user()) z0offsetval = (1.0-cur_frame_time)*option[z0_offset]()[ihelix] + cur_frame_time*option[z0_offset_end]()[ihelix]; //Linearly interpolate start and end values.
				core::Real z1offsetval = option[z1_offset]()[ihelix];
				if(option[animation_frames]()>2 && option[z1_offset_end].user()) z1offsetval = (1.0-cur_frame_time)*option[z1_offset]()[ihelix] + cur_frame_time*option[z1_offset_end]()[ihelix]; //Linearly interpolate start and end values.

				makebundle.add_helix();

				makebundle.helix(ihelix)->set_helix_length( static_cast<core::Size>(option[residue_repeats]()[ihelix]) );
				utility::vector1< std::string > restypevect; restypevect.push_back( option[residue_type]() );
				makebundle.helix(ihelix)->set_residue_name( restypevect );
				makebundle.helix(ihelix)->set_tail_residue_name( option[tail_residue_type]() );
				makebundle.helix(ihelix)->set_major_helix_params ( r0val, omega0val, deltaomega0val);
				makebundle.helix(ihelix)->set_invert_helix( option[invert_helix]()[ihelix] );
				makebundle.helix(ihelix)->set_delta_t( deltatval );
				makebundle.helix(ihelix)->set_z0_offset( z0offsetval );
				makebundle.helix(ihelix)->set_z1_offset( z1offsetval );

				utility::vector1 < core::Real > r1vals = option[r1]();

				utility::vector1 < core::Real > delta_omega1vals = option[delta_omega1_per_atom]();
				makebundle.helix(ihelix)->set_delta_omega1_all( deltaomega1val );

				utility::vector1 < core::Real > delta_z1vals = option[delta_z1_per_atom]();

				if( option[minor_helix_params].user()) {
					makebundle.helix(ihelix)->set_minor_helix_params_from_file( option[minor_helix_params]() );
				} else {
					makebundle.helix(ihelix)->set_minor_helix_params( r1vals, option[omega1](), option[z1](), delta_omega1vals, delta_z1vals );
				}
			}

			makebundle.apply(pose);

			char charbuffer[256];
			sprintf(charbuffer, "output_%04lu.pdb", iframe);

			pose.dump_pdb(std::string(charbuffer));
			TR << "Wrote " << charbuffer << "." << std::endl;
		}

		if(TR.visible()) {
			TR << "Finished test_makebundle_mover.cc.  Exiting." << std::endl;
			TR.flush();
		}

	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}
