// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to relax
/// @detailed
/// @author Mike Tyka, Monica Berrondo


#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/MiniRelax.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/scoring/MembraneTopology.hh> //pba
#include <basic/datacache/BasicDataCache.hh> //pba
#include <core/pose/datacache/CacheableDataType.hh> //pba

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/SuperimposeMover.hh>
#include <protocols/moves/ConstraintSetMover.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/electron_density/util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

#include <basic/Tracer.hh>

//Auto Headers
#include <protocols/relax/util.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer tr("protocols.relax");

using namespace core;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {
////////////////////////////////////////////////////////////////////////////////////////////////////

RelaxProtocolBaseOP
generate_relax_from_cmd( bool NULL_if_no_flag ) {
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::getScoreFunction();
	if ( option[ in::file::fullatom ]() )
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn );
	else
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn );

	// now add density scores
	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
	}

	RelaxProtocolBaseOP protocol;
	if ( option[ OptionKeys::relax::sequence_file ].user() ) {
		protocol = new FastRelax( scorefxn, option[ OptionKeys::relax::sequence_file ]() );
	} else if ( option[ OptionKeys::relax::script ].user() ) {
		protocol = new FastRelax( scorefxn, option[ OptionKeys::relax::script ]() );
	} else if ( option[ OptionKeys::relax::quick ]() ){
		protocol = new FastRelax( scorefxn, option[ OptionKeys::relax::default_repeats ]() );
	} else if ( option[ OptionKeys::relax::thorough ]() ){
		protocol = new FastRelax( scorefxn, 15 );
	} else if ( option[ OptionKeys::relax::fast ]() ) {
		protocol = new FastRelax( scorefxn, option[ OptionKeys::relax::default_repeats ]() );
	} else if ( option[ OptionKeys::relax::classic ]() ) {
		protocol = new ClassicRelax ( scorefxn );
	} else if ( option[ OptionKeys::relax::mini ]() ) {
		protocol = new MiniRelax( scorefxn );
	} else {
		// default relax should be a quick sequence relax
		if ( NULL_if_no_flag ){
			tr.Debug << "no relax protocol specified at command line" << std::endl;
			return NULL;
		}
		protocol = new FastRelax( scorefxn );
	}

	return protocol;
}

void setup_membrane_topology( pose::Pose & pose, std::string spanfile ) {
  //using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
  core::scoring::MembraneTopologyOP topologyOP = new core::scoring::MembraneTopology;
  pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
  core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY )() ));
  topology.initialize(spanfile);
}





void relax_pose( pose::Pose& pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag ) {
	using namespace basic::options;
	RelaxProtocolBaseOP protocol( generate_relax_from_cmd() );
	protocol->set_current_tag( tag );
	protocol->set_scorefxn( scorefxn );
	protocol->apply( pose );
}

}
}
