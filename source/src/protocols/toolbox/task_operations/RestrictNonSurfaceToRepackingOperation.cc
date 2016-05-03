// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperation.cc
/// @brief
/// @author Ron Jacak ronj@unc.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperation.hh>
#include <protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <set>

#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictNonSurfaceToRepackingOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

// Creator method
core::pack::task::operation::TaskOperationOP
RestrictNonSurfaceToRepackingOperationCreator::create_task_operation() const {
	return core::pack::task::operation::TaskOperationOP( new RestrictNonSurfaceToRepackingOperation );
}

void RestrictNonSurfaceToRepackingOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictNonSurfaceToRepackingOperation::provide_xml_schema( xsd );
}

std::string RestrictNonSurfaceToRepackingOperationCreator::keyname() const
{
	return RestrictNonSurfaceToRepackingOperation::keyname();
}

// default constructor
RestrictNonSurfaceToRepackingOperation::RestrictNonSurfaceToRepackingOperation() :
	surface_exposed_nb_cutoff_( 16 ) {}


// constructor with custom parameters
RestrictNonSurfaceToRepackingOperation::RestrictNonSurfaceToRepackingOperation( core::Size nb_cutoff ) {
	surface_exposed_nb_cutoff_ = nb_cutoff;
}

// destructor
RestrictNonSurfaceToRepackingOperation::~RestrictNonSurfaceToRepackingOperation() {}

// clone method, required by TaskOperation interface
TaskOperationOP RestrictNonSurfaceToRepackingOperation::clone() const {
	return TaskOperationOP( new RestrictNonSurfaceToRepackingOperation( *this ) );
}


// setter for nb_count cutoff. this allows users to vary how surface exposed a residue must be for it to be designed.
// more specifically, if this value is very low (e.g. 10, or 5), only the most surface-exposed residues will be designed.
// very few residues will have 5 or fewer neighbors and remain designable.  the rest will be set to repack only.
// the apply() method checks to see if a given pose position has GREATER THAN this number of neighbors using the tenA
// neighbor graph method num_neighbors_counting_self().
//
void RestrictNonSurfaceToRepackingOperation::surface_exposed_nb_cutoff( core::Size const nb_count ) {
	surface_exposed_nb_cutoff_ = nb_count;
}


// apply method
// Because there's no guarantee that the pose object has been scored at this point, we have to construct a "neighbor
// graph" ourselves. I'm going to use the implementation written by Steven L and John K in the PoseMetricCalculators
// for this task. The alternative would be to copy the pose object to a new pose (which is somewhat slow), score it
// (not trivial either), and lookup the information in the tenA_nb_graph. This alternative is better because it
// completely avoids the pose copy operation.
//
void RestrictNonSurfaceToRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const {

	using core::conformation::PointGraph;
	using core::conformation::PointGraphOP;

	PointGraphOP pg( new PointGraph ); // create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); // create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, 10.0 /* ten angstrom distance */ ); // create edges

	core::Size num_neighbors_ = 0;
	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {

		// a PointGraph is a typedef of UpperEdgeGraph< PointGraphVertexData, PointGraphEdgeData >
		// so any of the method in UpperEdgeGraph should be avail. here. The UpperEdgeGraph provides access to nodes
		// via a get_vertex() method, and each vertex can report back how many nbs it has.
		// So something that before was really complicated (nb count calculation) is done in <10 lines of code.
		// the assumption we're making here is that a pose residue position ii is the same index as the point graph vertex
		// that is indeed the case if you look at what the function residue_point_graph_from_pose().
		num_neighbors_ = pg->get_vertex(ii).num_neighbors_counting_self();

		// what about non-canonicals?  as long as the point graph can handle them, they should work fine.

		if ( num_neighbors_ > surface_exposed_nb_cutoff_ ) {
			// it's not at or below our cutoff, so limit this position to repacking only
			task.nonconst_residue_task( ii ).restrict_to_repacking();
		}

		// reset count for next position, just to be extra careful
		num_neighbors_ = 0;

	}

}

// parse_tag method; I believe this needs to be implemented so that the TaskOperationFactory can read this particular
// TaskOperation out of an XML file specifying task operations to perform
void RestrictNonSurfaceToRepackingOperation::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) {
	surface_exposed_nb_cutoff_ = tag->getOption< core::Size >( "surface_exposed_nb_count_cutoff" );
}

void RestrictNonSurfaceToRepackingOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	activate_common_simple_type( xsd, "non_negative_integer" );

	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::required_attribute( "surface_exposed_nb_count_cutoff", "non_negative_integer" ));
	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
