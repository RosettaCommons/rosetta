// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace core::kinematics;
using core::pose::Pose;
using core::pack::task::PackerTask;

RestrictToMoveMapChiOperation::RestrictToMoveMapChiOperation():
	parent()
{
	movemap_set_ = false;
	init();
}

RestrictToMoveMapChiOperation::RestrictToMoveMapChiOperation(MoveMapCOP movemap): parent() {
	movemap_ = movemap->clone();;
	movemap_set_ = true;
	init();
}

void
RestrictToMoveMapChiOperation::init() {

	set_design(false);
	set_include_neighbors(false);
	set_cutoff_distance(10.0);

}

void
RestrictToMoveMapChiOperation::init_for_equal_operator_and_copy_constructor(RestrictToMoveMapChiOperation & lhs, const RestrictToMoveMapChiOperation& rhs){
	if ( rhs.movemap_ ) {

		lhs.movemap_ = MoveMapOP( new MoveMap( *rhs.movemap_));
	}

	lhs.design_ = rhs.design_;
	lhs.include_neighbors_ = rhs.include_neighbors_;
	lhs.movemap_set_ = rhs.movemap_set_; //Since I can't assign a COP to null it seems.
	lhs.cutoff_ = rhs.cutoff_;
}

/*void
RestrictToMoveMapChiOperation::parse_tag( TagCOP tag, basic::datacache::DataMap & data){
set_cutoff(tag->getOption< core::Real >("cutoff", 10.0));
set_design(tag->getOption< bool > "design", false));
set_include_neighbors(tag->getOption< bool > "include_neighbors", false));

protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_, data, false);
}

*/

// AMW: no parse_tag or Creator.
void RestrictToMoveMapChiOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

RestrictToMoveMapChiOperation::RestrictToMoveMapChiOperation(const RestrictToMoveMapChiOperation& src):
	core::pack::task::operation::TaskOperation(src)
{
	init_for_equal_operator_and_copy_constructor( *this, src);
}

TaskOperationOP
RestrictToMoveMapChiOperation::clone() const {
	return TaskOperationOP( new RestrictToMoveMapChiOperation( *this ) );
}

RestrictToMoveMapChiOperation::~RestrictToMoveMapChiOperation()= default;

/* AMW: No Creator (had been commented out to save PyRosetta build)
void RestrictToMoveMapChiOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
RestrictToMoveMapChiOperation::provide_xml_schema( xsd );
}*/

void
RestrictToMoveMapChiOperation::set_movemap(core::kinematics::MoveMapCOP movemap) {
	movemap_ = movemap->clone();
	movemap_set_=true;
}

void
RestrictToMoveMapChiOperation::set_design(bool setting) {
	design_ = setting;
}

void
RestrictToMoveMapChiOperation::set_include_neighbors(bool setting) {
	include_neighbors_ = setting;
}

void
RestrictToMoveMapChiOperation::set_cutoff_distance(core::Real cutoff) {
	cutoff_ = cutoff;
}

void
RestrictToMoveMapChiOperation::apply(Pose const & pose, core::pack::task::PackerTask & task) const{
	if ( ! movemap_set_ ) return;

	core::pack::task::operation::PreventRepacking prevent_repacking;
	core::pack::task::operation::RestrictResidueToRepacking restrict_to_repacking;

	utility::vector1<bool> is_packable( pose.size(), false );
	utility::vector1<bool> is_designable( pose.size(), false );

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( movemap_->get_chi(i) ) {

			is_packable[i] = true;

			if ( design_ ) {
				is_designable[i] = true;
			}
		}
	}



	//Now we go through allowed residues so we don't overwrite anything.

	if ( include_neighbors_ ) {
		utility::vector1<bool> original_is_packable = is_packable;

		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( original_is_packable[i] ) {

				//Get neighbor distance
				for ( core::Size n = 1; n <= pose.size(); ++n ) {
					core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( n ).xyz( pose.residue( n ).nbr_atom() )) );
					if ( distance <= cutoff_ ) {
						is_packable[n] = true;
					}
				}
			}
		}
	}

	//Add the residues to the taskop + then apply.
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( is_packable[ i ] &&  ! is_designable[ i ] ) {
			restrict_to_repacking.include_residue(i);
		} else if ( ! is_packable[ i ] ) {
			prevent_repacking.include_residue(i);
		}
	}

	restrict_to_repacking.apply(pose, task);
	prevent_repacking.apply(pose, task);
}



} //operation
} //toolbox
} //protocols
