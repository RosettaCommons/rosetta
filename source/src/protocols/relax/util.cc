// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/relax/util.cc
/// @brief initialization protocols for relax and utility functions
/// @details
/// @author Mike Tyka, Monica Berrondo
/// @author Roland A. Pache

//Project Headers
#include <protocols/relax/util.hh>

//Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh> //pba
#include <core/scoring/MembraneTopology.hh> //pba
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>

//Protocol Headers
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/MiniRelax.hh>
#include <protocols/relax/CentroidRelax.hh>

//Basic Headers
#include <basic/datacache/BasicDataCache.hh> //pba
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.relax" );

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
	scorefxn = core::scoring::get_score_function();
	if ( option[ in::file::fullatom ]() || option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn );
	}

	// now add density scores
	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
	}

	RelaxProtocolBaseOP protocol;
	if ( option[ OptionKeys::relax::sequence_file ].user() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, option[ OptionKeys::relax::sequence_file ]() ) );
	} else if ( option[ OptionKeys::relax::script ].user() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, option[ OptionKeys::relax::script ]() ) );
	} else if ( option[ OptionKeys::relax::quick ]() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, option[ OptionKeys::relax::default_repeats ]() /*default 5*/) );
	} else if ( option[ OptionKeys::relax::thorough ]() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, 15 ) );
	} else if ( option[ OptionKeys::relax::fast ]() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, option[ OptionKeys::relax::default_repeats ]() /*default 5*/) );
	} else if ( option[ OptionKeys::relax::classic ]() ) {
		protocol = RelaxProtocolBaseOP( new ClassicRelax ( scorefxn ) );
	} else if ( option[ OptionKeys::relax::mini ]() ) {
		protocol = RelaxProtocolBaseOP( new MiniRelax( scorefxn ) );
	} else if ( option[ OptionKeys::relax::centroid_mode ]() ) {
		protocol = RelaxProtocolBaseOP( new CentroidRelax() );
	} else {
		// default relax should be a quick sequence relax
		if ( NULL_if_no_flag ) {
			TR.Debug << "no relax protocol specified at command line" << std::endl;
			return NULL;
		}
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn ) );
	}

	return protocol;
}

// RosettaMembrane from 2006
void setup_membrane_topology( pose::Pose & pose, std::string spanfile ) {
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	core::scoring::MembraneTopologyOP topologyOP( new core::scoring::MembraneTopology );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
	core::scoring::MembraneTopology & topology=*( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
	topology.initialize(spanfile);
}

void make_dna_rigid( pose::Pose & pose, core::kinematics::MoveMap & mm){
	using namespace core::conformation;
	//if DNA present set so it doesn't move
	for ( core::Size i=1; i<=pose.total_residue() ; ++i )      {
		if ( pose.residue(i).is_DNA() ) {
			TR << "turning off DNA bb and chi move" << std::endl;
			mm.set_bb( i, false );
			mm.set_chi( i, false );
		}
	}
}

void setup_for_dna( core::scoring::ScoreFunction & scorefxn) {

	scoring::methods::EnergyMethodOptions options( scorefxn.energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn.set_energy_method_options( options );

}

void relax_pose( pose::Pose& pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag ) {
	using namespace basic::options;
	RelaxProtocolBaseOP protocol( generate_relax_from_cmd() );
	protocol->set_current_tag( tag );
	protocol->set_scorefxn( scorefxn );
	protocol->apply( pose );
}

void fixH(core::pose::Pose & pose) {
	for ( Size i = 1; i <= pose.n_residue(); ++i ) {
		pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(i);
	}
}


void
cyclize_pose(core::pose::Pose & pose) {
	using namespace core;
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace chemical;
	using core::id::AtomID;
	Size N = pose.n_residue();
	for ( Size i = 1; i <= N; ++i ) {
		if ( pose.residue(i).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
		if ( pose.residue(i).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
		if ( pose.residue(i).has_variant_type(CUTPOINT_UPPER) ) core::pose::remove_variant_type_from_pose_residue(pose,CUTPOINT_UPPER,i);
		if ( pose.residue(i).has_variant_type(CUTPOINT_LOWER) ) core::pose::remove_variant_type_from_pose_residue(pose,CUTPOINT_LOWER,i);
	}
	if ( !pose.residue(1).has_variant_type(CUTPOINT_UPPER) ) {
		core::pose::add_variant_type_to_pose_residue(pose,CUTPOINT_UPPER,1);
	}
	if ( !pose.residue(N).has_variant_type(CUTPOINT_LOWER) ) {
		core::pose::add_variant_type_to_pose_residue(pose,CUTPOINT_LOWER,N);
	}
	pose.conformation().declare_chemical_bond( 1, "N", N, "C" );
	fixH(pose);

	pose.conformation().update_polymeric_connection(1);

	using namespace core::scoring::constraints;
	AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.n_residue()).atom_index("OVL1"), pose.n_residue() );
	AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.n_residue()).atom_index("OVL2"), pose.n_residue() );
	AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.n_residue()).atom_index(   "C"), pose.n_residue() );
	core::scoring::func::FuncOP fx1( new core::scoring::func::HarmonicFunc(0.0,0.1) );
	pose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(a1,a2,fx1) ) ));
	core::scoring::func::FuncOP fx2( new core::scoring::func::HarmonicFunc(0.0,0.1) );
	pose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(b1,b2,fx2) ) ));
	core::scoring::func::FuncOP fx3( new core::scoring::func::HarmonicFunc(0.0,0.1) );
	pose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(c1,c2,fx3) ) ));
}

}
}
