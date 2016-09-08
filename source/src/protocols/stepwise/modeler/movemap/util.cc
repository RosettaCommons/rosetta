// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/movemap/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/stream_util.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.movemap.util" );

/////////////////////////////////////////////////////////////////////////////////////////////
//
// All movemap setting for stepwise assembly/monte carlo is now tucked into here.
//  Might consider making this a class.
//  Also might move "core" figure_out_stepwise_movemap() routine that goes from
//    AtomLevelDomainMap to MoveMap into AtomLevelDomainMap. Its pretty general and good.
//
//   -- rhiju, 2014
//
/////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace movemap {

/////////////////////////////////////////////////////////////////////////////////////////////
// used by legacy StepWiseProteinMinimizer.cc
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
	core::pose::Pose const & pose,
	utility::vector1< Size > const & working_minimize_res,
	bool const move_takeoff_torsions /* = true */ ){
	toolbox::AtomLevelDomainMapOP atom_level_domain_map( new toolbox::AtomLevelDomainMap( pose ) );
	figure_out_stepwise_movemap( mm, atom_level_domain_map, pose, working_minimize_res, move_takeoff_torsions );
}

/////////////////////////////////////////////////////////////////////////////////////////////
// used by legacy StepWiseRNA_Minimizer.cc
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
	toolbox::AtomLevelDomainMapOP & atom_level_domain_map,
	core::pose::Pose const & pose,
	utility::vector1< Size > const & working_fixed_res,
	utility::vector1< Size > const & working_extra_minimize_res,
	bool const move_takeoff_torsions /* = true */ ){

	utility::vector1< core::Size > working_minimize_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !working_fixed_res.has_value( n ) || working_extra_minimize_res.has_value( n ) ) working_minimize_res.push_back( n );
	}
	atom_level_domain_map = toolbox::AtomLevelDomainMapOP( new toolbox::AtomLevelDomainMap( pose ) );
	figure_out_stepwise_movemap( mm, atom_level_domain_map, pose, working_minimize_res, move_takeoff_torsions );
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
	toolbox::AtomLevelDomainMapOP atom_level_domain_map,
	core::pose::Pose const & pose,
	utility::vector1< Size > const & working_minimize_res,
	bool const move_takeoff_torsions /* = true */ ) {
	using namespace core::id;
	Size const nres( pose.size() );
	for ( Size n = 1; n <= nres; n++ ) {
		if ( working_minimize_res.has_value( n ) )  atom_level_domain_map->set( n, true );
		if ( !working_minimize_res.has_value( n ) ) atom_level_domain_map->set_fixed_if_moving( n );
		if ( pose.residue_type( n ).has_variant_type( chemical::CUTPOINT_LOWER ) ) {
			atom_level_domain_map->set( NamedAtomID( "OVL1",  n), pose, true );
			atom_level_domain_map->set( NamedAtomID( "OVL2",  n), pose, true );
		}
		if ( pose.residue_type( n ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			atom_level_domain_map->set( NamedAtomID( "OVU1", n ), pose, true );
		}
		if ( pose.residue_type( n ).is_RNA() ) {
			// there's a problem with looking at VIRTUAL PHOSPHATE -- it can change if 5' packing phosphate is 'packed' in.
			// instead let it move -- but rely on the fact that there are no score terms computed (so derivative will be zero)
			//   if ( pose.residue_type(n).has_variant_type( chemical::VIRTUAL_PHOSPHATE ) ) atom_level_domain_map->set_phosphate( n, pose, false );

			// following are 'conventions' that need to be hard-coded -- for RNA want
			// entire suites (not just connection torsions ) to move between residues in different domains
			if ( pose.residue_type(n).has_variant_type( chemical::CUTPOINT_UPPER    ) ) atom_level_domain_map->set_phosphate( n, pose, true );

			if ( n > 1 ) {
				if ( move_takeoff_torsions && pose.residue_type( n-1 ).is_RNA() ) {
					if ( working_minimize_res.has_value( n-1 ) ) atom_level_domain_map->set_phosphate( n, pose, true );
					if ( atom_level_domain_map->get_domain( NamedAtomID( " O3'", n-1 ), pose ) !=
							atom_level_domain_map->get_domain( NamedAtomID( " P  ", n )  , pose ) ) atom_level_domain_map->set_phosphate( n, pose, true );
				}
			}

			if ( core::pose::full_model_info::check_sample_sugar_in_full_model_info( pose, n ) ) atom_level_domain_map->set_sugar( n, pose, true );
		}
		if ( pose.residue_type(n).is_protein() ) {
			if ( move_takeoff_torsions ) {
				if ( n > 1 && pose.residue_type( n-1 ).is_protein() && working_minimize_res.has_value( n-1 ) &&
						pose.residue_type( n ).has( " H  " ) /* not there for Pro*/ ) {
					atom_level_domain_map->set( NamedAtomID(" H  ", n ), pose, true );
				}
				if ( n < nres && pose.residue_type( n+1 ).is_protein() && working_minimize_res.has_value( n+1 ) ) {
					atom_level_domain_map->set( NamedAtomID(" O  ", n ), pose, true );
				}
			} else { // argh, proteins. This is to match old KIC runs. Perhaps should just deprecate.
				if ( n > 1 && pose.residue_type( n-1 ).is_protein() &&
						!working_minimize_res.has_value( n-1 ) ) {
					atom_level_domain_map->set_domain( NamedAtomID(" N  ", n ), pose,
						atom_level_domain_map->get_domain( NamedAtomID( " C  ", n-1 ), pose ) );
					atom_level_domain_map->set_domain( NamedAtomID(" CA ", n ), pose,
						atom_level_domain_map->get_domain( NamedAtomID( " C  ", n-1 ), pose ) );
				}
				if ( n < nres && pose.residue_type( n+1 ).is_protein() &&
						!working_minimize_res.has_value( n+1 ) ) {
					atom_level_domain_map->set_domain( NamedAtomID(" C  ", n ), pose,
						atom_level_domain_map->get_domain( NamedAtomID( " N  ", n+1 ), pose ) );
				}
			}
		}
	}

	// atom_level_domain_map->show();
	atom_level_domain_map->setup_movemap( mm, pose );
	// output_movemap( mm, pose, TR );
}


} //movemap
} //modeler
} //stepwise
} //protocols
