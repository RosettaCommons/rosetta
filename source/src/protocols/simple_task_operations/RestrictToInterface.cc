// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RestrictToInterface.cc
/// @brief When passed to a PackerTask, pack/design is limited to the protein interface
/// @author ashworth

#include <protocols/simple_task_operations/RestrictToInterface.hh>
#include <protocols/simple_task_operations/RestrictToInterfaceCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <protocols/scoring/Interface.hh>

#include <utility/vector0.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.toolbox.task_operations.RestrictToInterface" );

namespace protocols {
namespace simple_task_operations {

using namespace core;
using namespace scoring;
using namespace core::pack::task::operation;
using namespace utility::tag;

//RestrictTaskForDocking::RestrictTaskForDocking()
// : parent(),
//  scorefxn_( 0 ),
//  rb_jump_( 0 ),
//  include_current_( true ),
//  distance_( 0 )
//{}
//
//RestrictTaskForDocking::RestrictTaskForDocking(
// ScoreFunctionCOP scorefxn,
// core::Size rb_jump,
// bool include_current,
// core::Real distance
//) : parent(),
//  scorefxn_( scorefxn ),
//  rb_jump_( rb_jump ),
//  include_current_( include_current ),
//  distance_( distance )
//{}
//
//RestrictTaskForDocking::~RestrictTaskForDocking(){}
//
//
//task::operation::TaskOperationOP RestrictTaskForDocking::clone() const
//{
// return new RestrictTaskForDocking( *this );
//}
//
//void
//RestrictTaskForDocking::apply(
// pose::Pose const & pose,
// task::PackerTask & task
//) const
//{
// task.initialize_from_command_line().restrict_to_repacking().or_include_current( include_current_ );
//
// runtime_assert( scorefxn_ != 0 );
// runtime_assert( rb_jump_ );
// runtime_assert( distance_ );
// // (existing comment) /// why is this still necessary???
//// (*scorefxn_)(pose);
//// scorefxn_->accumulate_residue_total_energies( pose );
//
// core::conformation::Interface interface( rb_jump_ );
// interface.distance( distance_ );
// interface.calculate( pose );
// interface.set_pack( pose, &task );
//}

DockingNoRepack1::DockingNoRepack1()
: parent(), rb_jump_(1)
{}

DockingNoRepack1::DockingNoRepack1( int rb_jump_in )
: parent(), rb_jump_(rb_jump_in)
{}

DockingNoRepack1::~DockingNoRepack1()= default;

core::pack::task::operation::TaskOperationOP DockingNoRepack1::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new DockingNoRepack1( *this ) );
}

void DockingNoRepack1::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

/* AMW: No Creator
void DockingNoRepack1Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
DockingNoRepack1::provide_xml_schema( xsd );
} */

void
DockingNoRepack1::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const
{
	using namespace core;
	Size cutpoint ( pose.fold_tree().cutpoint_by_jump( rb_jump_ ) );
	for ( Size ii = 1 ; ii <= cutpoint; ++ii ) {
		task.nonconst_residue_task( ii ).prevent_repacking();
	}
}

DockingNoRepack2::DockingNoRepack2()
: parent(), rb_jump_(1)
{}

DockingNoRepack2::DockingNoRepack2( int rb_jump_in )
: parent(), rb_jump_(rb_jump_in)
{}

DockingNoRepack2::~DockingNoRepack2()= default;

void DockingNoRepack2::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

/* AMW: No Creator
void DockingNoRepack2Creator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
DockingNoRepack2::provide_xml_schema( xsd );
}*/

core::pack::task::operation::TaskOperationOP DockingNoRepack2::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new DockingNoRepack2( *this ) );
}

void
DockingNoRepack2::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const
{
	Size cutpoint ( pose.fold_tree().cutpoint_by_jump( rb_jump_ ) );
	for ( Size ii = cutpoint+1 ; ii <= pose.size(); ++ii ) {
		task.nonconst_residue_task( ii ).prevent_repacking();
	}
}

RestrictToInterface::~RestrictToInterface()= default;

RestrictToInterface::RestrictToInterface(utility::vector1<bool> loop_residues):
	parent(), distance_( 8 ), loopy_interface_( true ), include_all_water_( false ) {
	set_movable_jumps( utility::vector1_int() );
	ObjexxFCL::FArray1D_bool hack_loop_residues( loop_residues.size(), false );
	for ( core::Size ii = 1; ii <= loop_residues.size(); ii++ ) {
		if ( loop_residues[ii] ) {
			hack_loop_residues(ii) = true;
		}
	}
	loop_residues_ = hack_loop_residues;
}

core::pack::task::operation::TaskOperationOP RestrictToInterfaceCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictToInterface );
}

void RestrictToInterfaceCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToInterface::provide_xml_schema( xsd );
}

std::string RestrictToInterfaceCreator::keyname() const
{
	return RestrictToInterface::keyname();
}

core::pack::task::operation::TaskOperationOP RestrictToInterface::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictToInterface( *this ) );
}

void
RestrictToInterface::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const
{
	using core::Size;
	utility::vector1<bool> is_interface( pose.size(), false );

	core::Size num_jump_ = movable_jumps().size();
	for ( Size jj=1; jj<=num_jump_; jj++ ) {
		//protocols::scoring::Interface interface( movable_jumps()[jj] );
		//interface.distance( distance_ );
		//interface.calculate( pose );
		//for ( Size ii=1; ii<=pose.size(); ++ii ) {
		// if ( interface.is_interface(ii) ) {
		//  is_interface[ii] = true;
		// }
		//}

		//fd logic above is very problematic:
		// * it uses only the stored energy graph, which may be computed with less distance than distance_
		// * it depends on interaction radius of last scorefunction used, which may vary, giving diff results
		// * Instead we replace this logic.  Similar to it, though, we use neighboratoms for all but ligands.
		ObjexxFCL::FArray1D_bool partition( pose.size(),false );
		pose.fold_tree().partition_by_jump( movable_jumps()[jj], partition );

		for ( Size ires=1; ires<=pose.size(); ++ires ) {
			if ( pose.residue_type(ires).is_water() ) continue;
			bool i_is_lig = pose.residue_type(ires).is_ligand();
			core::Size niatm = i_is_lig ? pose.residue_type(ires).nheavyatoms() : 1;
			for ( Size iatm=1; iatm<=niatm; ++iatm ) {
				core::Vector xyz_i = i_is_lig ?
					pose.residue(ires).xyz( iatm ) : pose.residue(ires).xyz( pose.residue(ires).nbr_atom() );

				for ( Size jres=ires+1; jres<=pose.size(); ++jres ) {
					if ( partition(ires) == partition(jres) ) continue;
					if ( pose.residue_type(jres).is_water() ) continue;
					bool j_is_lig = pose.residue_type(jres).is_ligand();
					core::Size njatm = j_is_lig ? pose.residue_type(jres).nheavyatoms() : 1;
					for ( Size jatm=1; jatm<=njatm; ++jatm ) {
						core::Vector xyz_j = j_is_lig ?
							pose.residue(jres).xyz( jatm ) : pose.residue(jres).xyz( pose.residue(jres).nbr_atom() );

						if ( (xyz_i-xyz_j).length_squared() <= distance_*distance_ ) {
							is_interface[ires] = is_interface[jres] = true;
						}
					} // jatm
				} // jj

			} // iatm
		} // ii
	}

	if ( loopy_interface_ ) {
		for ( Size ii=1; ii<=pose.size(); ++ii ) {
			if ( loop_residues_(ii) ) {
				is_interface[ii] = true;
			}
		}
	}

	// set all waters to part of the interface
	if ( include_all_water_ ) {
		for ( Size ii=1; ii<=pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).is_water() ) {
				is_interface[ii] = true;
			}
		}
	}

	for ( Size ii=1; ii<=pose.size(); ++ii ) {
		if ( !is_interface[ii] ) { //|| pose.residue(ii).is_ligand() )
			task.nonconst_residue_task( ii ).prevent_repacking();
		}
	}
	if ( core::pose::symmetry::is_symmetric( pose ) ) symmetric_task( pose, task );
}

void RestrictToInterface::symmetric_task(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const {
	using namespace conformation::symmetry;

	auto const & SymmConf (
		dynamic_cast< SymmetricConformation const &> ( pose.conformation()) );

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !SymmConf.Symmetry_Info()->chi_is_independent(i) ) {
			task.nonconst_residue_task( i ).prevent_repacking();
		}
	}
}

void RestrictToInterface::rb_jump( int jump_in ) {

	add_movable_jump( jump_in );
}

/*
void RestrictToInterface::set_movable_jumps( utility::vector1_int const movable_jumps ) {
rb_jump_ = movable_jumps;
}
*/

void RestrictToInterface::distance( core::Real const distance_in ) {
	distance_ = distance_in;
}

void
RestrictToInterface::parse_tag( TagCOP tag , DataMap & )
{
	add_movable_jump( ( tag->getOption< core::Size >( "jump", 1 ) ) );
	distance_ = tag->getOption< core::Real >( "distance", 8 )  ;
	include_all_water_ = tag->getOption<bool>( "include_all_water", 0 );
}

void RestrictToInterface::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "jump", xsct_non_negative_integer, "XRW TO DO",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "distance", xsct_real, "XRW TO DO",  "8"  )
		+ XMLSchemaAttribute::attribute_w_default(  "include_all_water", xsct_rosetta_bool, "add all water to interface",  "false"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO" );
}


}
}
