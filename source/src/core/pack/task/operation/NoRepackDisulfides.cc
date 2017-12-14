// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/NoRepackDisulfides.cc
/// @brief  prevent disulfides from being repacked; assumes disulfide info in
///         Pose is up-to-date
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/NoRepackDisulfidesCreator.hh>

// package headers
#include <core/pack/task/operation/task_op_schemas.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


// C++ headers


namespace core {
namespace pack {
namespace task {
namespace operation {


// static
static basic::Tracer TR( "core.pack.task.operation.NoRepackDisulfides" );


/// @brief default constructor
NoRepackDisulfides::NoRepackDisulfides() :
	Super()
{}


/// @brief copy constructor
NoRepackDisulfides::NoRepackDisulfides( NoRepackDisulfides const & /*rval*/ ) = default;


/// @brief default destructor
NoRepackDisulfides::~NoRepackDisulfides() = default;

/// @brief clone this object
NoRepackDisulfides::TaskOperationOP NoRepackDisulfides::clone() const {
	return NoRepackDisulfides::TaskOperationOP( new NoRepackDisulfides( *this ) );
}

/// @brief apply operations to PackerTask
void NoRepackDisulfides::apply( Pose const & pose, PackerTask & task ) const {
	using core::Size;
	using core::chemical::aa_cys;
	using core::chemical::DISULFIDE;
	using core::conformation::Residue;

	core::Size nres = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}


	for ( Size i = 1, ie = nres; i <= ie; ++i ) {
		Residue const & res = pose.residue( i );

		// residue appears to be a disulfide...?
		// amw removing the check for aa_cys
		// because we want this to be able to apply to NCAA disulfides
		if ( ! res.has_variant_type( DISULFIDE ) ) continue;

		bool no_repack = true;

		// try and double-check if possible
		Size sg_index = res.atom_index( res.type().get_disulfide_atom_name() );
		Size const sg_conn_index = res.type().residue_connection_id_for_atom( sg_index );
		Residue const & partner_res = pose.residue( res.residue_connection_partner( sg_conn_index ) );

		if ( !partner_res.has_variant_type( DISULFIDE ) ) {
			no_repack = false;
		}

		// set repack status
		if ( no_repack ) {
			TR.Debug << "found disulfide residue " << i << ", preventing repack at this position" << std::endl;
			task.nonconst_residue_task( i ).prevent_repacking();
		} else {
			TR.Warning << "residue " << i << " marked as disulfide but has no partner, allowing repack at this position" << std::endl;
		}

	} // foreach residue

}


std::string NoRepackDisulfides::keyname() { return "NoRepackDisulfides"; }

void NoRepackDisulfides::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty(
		xsd, keyname(),
		"Do not allow disulfides to repack.");
}

TaskOperationOP NoRepackDisulfidesCreator::create_task_operation() const
{
	return TaskOperationOP( new NoRepackDisulfides );
}

std::string NoRepackDisulfidesCreator::keyname() const { return NoRepackDisulfides::keyname(); }

void NoRepackDisulfidesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NoRepackDisulfides::provide_xml_schema( xsd );
}

} // namespace operation
} // namespace task
} // namespace pack
} // namespace core

