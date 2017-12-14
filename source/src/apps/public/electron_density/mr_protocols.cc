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


#include <protocols/hybridization/MRMover.hh>

#include <protocols/viewer/viewers.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


#include <iostream>
#include <string>
#include <sstream>


OPT_1GRP_KEY(Integer, MR, max_gaplength_to_model)
OPT_1GRP_KEY(Real, MR, cen_dens_wt)
OPT_1GRP_KEY(Real, MR, fa_dens_wt)
OPT_1GRP_KEY(Real, MR, cst_wt)
OPT_1GRP_KEY(Real, MR, relax_cst_wt)
OPT_1GRP_KEY(Boolean, MR, fast)
OPT_1GRP_KEY(Real, MR, censcale)
OPT_1GRP_KEY(StringVector, MR, disulf)
OPT_1GRP_KEY(String, MR, mode)

static basic::Tracer TR( "rosetta_MR" );

void*
my_main( void* ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves;
	using namespace core::fragment;


	if ( option[ OptionKeys::MR::mode ].user() ) {
		TR << "The flag -MR::mode is no longer used.  Ignoring." << std::endl;
	}

	protocols::hybridization::MRMoverOP do_MR( new protocols::hybridization::MRMover );

	do_MR->set_max_gaplength_to_model( option[ OptionKeys::MR::max_gaplength_to_model ]() );
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() ) {
		do_MR->set_symmdef_file( option[ OptionKeys::symmetry::symmetry_definition ]() );
	}
	if ( option[ OptionKeys::MR::disulf ].user() ) {
		do_MR->set_disulf( option[ OptionKeys::MR::disulf ]() );
	}

	do_MR->set_disulf( option[ OptionKeys::MR::disulf ]() );  // force disulfides
	do_MR->set_censcale( option[ OptionKeys::MR::censcale ]() );  // scale # centroid cycles
	do_MR->set_cen_cst_weight( option[ OptionKeys::MR::cst_wt ]() );
	do_MR->set_fa_cst_weight( option[ OptionKeys::MR::relax_cst_wt ]() );

	// fragment files: BACKWARDS COMPATABILITY with -loop options
	if ( option[ OptionKeys::loops::frag_files ].user() ) {
		FileVectorOption frag_files( option[ OptionKeys::loops::frag_files ] );
		for ( core::Size i=1; i<=frag_files.size(); ++i ) {
			if ( frag_files[i] == std::string("none") ) continue;

			FragSetOP frag_lib_op = (FragmentIO().read_data( frag_files[i] ));
			if ( frag_lib_op->max_frag_length() >= 7 ) {
				do_MR->set_big_fragments( frag_lib_op );
			} else {
				do_MR->set_small_fragments( frag_lib_op );
			}
		}
	}

	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		do_MR->set_centroid_density_weight( option[ OptionKeys::MR::cen_dens_wt ]() );
		do_MR->set_fullatom_density_weight( option[ OptionKeys::MR::fa_dens_wt ](), option[ OptionKeys::MR::fast ]() );
	}

	// run
	protocols::jd2::JobDistributor::get_instance()->go( do_MR );

	return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	try {
		NEW_OPT(MR::max_gaplength_to_model, "max gaplength to rebuild", 8);
		NEW_OPT(MR::cen_dens_wt, "centroid density weight", 4.0);
		NEW_OPT(MR::fa_dens_wt, "fullatom density weight", 1.0);
		NEW_OPT(MR::cst_wt, "add constraints to centroid stage", false);
		NEW_OPT(MR::relax_cst_wt, "add constraints to fullatom stage", false);
		NEW_OPT(MR::fast, "fast mode", false);
		NEW_OPT(MR::censcale, "scale # of centroid cycles", 1.0);
		NEW_OPT(MR::disulf, "force a disulfide patterning", utility::vector1<std::string>());
		NEW_OPT(MR::mode, "legacy flag; unused", "X");

		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	}
catch (utility::excn::Exception const & e) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
