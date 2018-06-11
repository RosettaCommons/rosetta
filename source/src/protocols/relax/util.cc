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
#include <core/pose/variant_util.hh>
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
#include <utility/options/keys/OptionKeyList.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.relax" );

using namespace core;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {
////////////////////////////////////////////////////////////////////////////////////////////////////

RelaxProtocolBaseOP
generate_relax_from_cmd( bool NULL_if_no_flag ) {
	return generate_relax_from_cmd( basic::options::option, NULL_if_no_flag );
}

RelaxProtocolBaseOP
generate_relax_from_cmd(
	utility::options::OptionCollection const & options,
	bool NULL_if_no_flag
) {
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;

	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::get_score_function();
	if ( options[ OptionKeys::in::file::fullatom ]() || options[ OptionKeys::constraints::cst_fa_file ].user() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scorefxn );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn );
	}

	// now add density scores
	if ( options[ OptionKeys::edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
	}

	RelaxProtocolBaseOP protocol;
	if ( options[ OptionKeys::relax::sequence_file ].user() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, options[ OptionKeys::relax::sequence_file ]() ) );
	} else if ( options[ OptionKeys::relax::script ].user() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, options[ OptionKeys::relax::script ]() ) );
	} else if ( options[ OptionKeys::relax::quick ]() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, options[ OptionKeys::relax::default_repeats ]() /*default 5*/) );
	} else if ( options[ OptionKeys::relax::thorough ]() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, 15 ) );
	} else if ( options[ OptionKeys::relax::fast ]() ) {
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn, options[ OptionKeys::relax::default_repeats ]() /*default 5*/) );
	} else if ( options[ OptionKeys::relax::classic ]() ) {
		protocol = RelaxProtocolBaseOP( new ClassicRelax ( scorefxn ) );
	} else if ( options[ OptionKeys::relax::mini ]() ) {
		protocol = RelaxProtocolBaseOP( new MiniRelax( scorefxn ) );
	} else if ( options[ OptionKeys::relax::centroid_mode ]() ) {
		protocol = RelaxProtocolBaseOP( new CentroidRelax() );
	} else {
		// default relax should be a quick sequence relax
		if ( NULL_if_no_flag ) {
			TR.Debug << "no relax protocol specified at command line" << std::endl;
			return nullptr;
		}
		protocol = RelaxProtocolBaseOP( new FastRelax( scorefxn ) );
	}

	return protocol;
}

void
options_for_generate_relax_from_cmd(
	utility::options::OptionKeyList & opts
)
{
	using namespace basic::options;

	opts
		+ OptionKeys::in::file::fullatom
		+ OptionKeys::constraints::cst_fa_file
		+ OptionKeys::edensity::mapfile
		+ OptionKeys::relax::sequence_file
		+ OptionKeys::relax::script
		+ OptionKeys::relax::quick
		+ OptionKeys::relax::thorough
		+ OptionKeys::relax::fast
		+ OptionKeys::relax::classic
		+ OptionKeys::relax::mini
		+ OptionKeys::relax::centroid_mode;
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
	for ( core::Size i=1; i<=pose.size() ; ++i )      {
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
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
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
	core::Size N = pose.size();
	for ( core::Size i = 1; i <= N; ++i ) {
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
	AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.size()).atom_index("OVL1"), pose.size() );
	AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.size()).atom_index("OVL2"), pose.size() );
	AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.size()).atom_index(   "C"), pose.size() );
	core::scoring::func::FuncOP fx1( new core::scoring::func::HarmonicFunc(0.0,0.1) );
	pose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(a1,a2,fx1) ) ));
	core::scoring::func::FuncOP fx2( new core::scoring::func::HarmonicFunc(0.0,0.1) );
	pose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(b1,b2,fx2) ) ));
	core::scoring::func::FuncOP fx3( new core::scoring::func::HarmonicFunc(0.0,0.1) );
	pose.add_constraint(scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(c1,c2,fx3) ) ));
}
core::scoring::constraints::ConstraintCOPs
add_coordinate_constraint_func_atoms( core::pose::Pose const & pose, core::Size const resnum, core::conformation::Residue const & rsd_i, core::scoring::func::HarmonicFuncOP coord_cst_func, core::id::AtomID const & anchor_atom )

{
	using namespace core::scoring::constraints;

	ConstraintCOPs cst;


	{ // defining constraints scope
		using namespace core::conformation;
		using namespace core::chemical;

		ResidueType const & rsd_type( rsd_i.type() );
		AA const aa = rsd_i.aa();
		switch ( aa ){
		using namespace core::id;
		case( aa_tyr ) : //Keep placement of OH group in place
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CZ" ), resnum ), anchor_atom, rsd_i.xyz( "CZ" ),coord_cst_func) ));
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OH" ), resnum ), anchor_atom, rsd_i.xyz( "OH" ),coord_cst_func ) ));
			break;
		case( aa_trp ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NE1" ), resnum ), anchor_atom, rsd_i.xyz( "NE1" ),coord_cst_func) ));
			break;
		case( aa_gln ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OE1" ), resnum ), anchor_atom, rsd_i.xyz( "OE1" ),coord_cst_func ) ));
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NE2" ), resnum ), anchor_atom, rsd_i.xyz( "NE2" ),coord_cst_func ) ));
			break;
		case( aa_glu ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OE1" ), resnum ), anchor_atom, rsd_i.xyz( "OE1" ),coord_cst_func ) ) );
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OE2" ), resnum ), anchor_atom, rsd_i.xyz( "OE2" ),coord_cst_func ) ) );
			break;
		case( aa_arg ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NH1" ), resnum ), anchor_atom, rsd_i.xyz( "NH1" ),coord_cst_func ) ) );
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NH2" ), resnum ) , anchor_atom, rsd_i.xyz( "NH2" ),coord_cst_func  ) ) );
			break;
		case( aa_lys ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NZ" ), resnum ), anchor_atom, rsd_i.xyz( "NZ" ),coord_cst_func ) ) );
			break;
		case( aa_asn ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OD1" ), resnum ), anchor_atom, rsd_i.xyz( "OD1" ),coord_cst_func ) ) );
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "ND2" ), resnum ), anchor_atom, rsd_i.xyz( "ND2" ),coord_cst_func ) ) );
			break;
		case( aa_his ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "ND1" ), resnum ) , anchor_atom, rsd_i.xyz( "ND1" ),coord_cst_func ) ) );
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NE2" ), resnum ), anchor_atom, rsd_i.xyz( "NE2" ),coord_cst_func ) ) );
			break;
		case( aa_asp ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OD1" ), resnum ), anchor_atom, rsd_i.xyz( "OD1" ),coord_cst_func ) ) );
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OD2" ), resnum ), anchor_atom, rsd_i.xyz( "OD2" ),coord_cst_func ) ) );
			break;
		case( aa_met ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "SD" ), resnum ), anchor_atom, rsd_i.xyz( "SD" ),coord_cst_func ) ) );
			break;
		case( aa_cys ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "SG" ), resnum ), anchor_atom, rsd_i.xyz( "SG" ), coord_cst_func ) ) );
			// Preserve the CB, since the chi angle is important in disulfides
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CB" ), resnum ), anchor_atom, rsd_i.xyz( "CB" ), coord_cst_func ) ) );
			break;
		case( aa_thr ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OG1" ), resnum ), anchor_atom, rsd_i.xyz( "OG1" ), coord_cst_func ) ) );
			break;
		case( aa_ser ) :
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OG" ), resnum ), anchor_atom, rsd_i.xyz( "OG" ), coord_cst_func ) ) );
			break;
		default :
			TR.Warning << "Cannot add functional group coordinate constraints to residue type. Ignoring." << std::endl;
		}
	}//defining constraints scope
	TR<<"Constraining residue "<<pose.residue( resnum ).name()<<resnum<<std::endl;
	return( cst );
}

}
}
