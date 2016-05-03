// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/InteractingRotamerExplosion.cc
/// @brief
/// @author Florian Richter, florian.richter.1@hu-berlin.de, june '14

// Unit Headers
#include <protocols/toolbox/task_operations/InteractingRotamerExplosion.hh>
#include <protocols/toolbox/task_operations/InteractingRotamerExplosionCreator.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/rotamer_set_operations/AddGood2BPairEnergyRotamers.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pack/task/PackerTask.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

//#include <core/pack/task/operation/ResLvlTaskOperations.hh>
//#include <core/pack/task/operation/OperateOnCertainResidues.hh>


// C++ Headers
#include <set>


static THREAD_LOCAL basic::Tracer tr( "protocols.toolbox.TaskOperations.InteractingRotamerExplosion" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

InteractingRotamerExplosion::InteractingRotamerExplosion() :
	string_target_seqpos_(""),
	score_cutoff_( -0.5 ),
	ex_level_( 4 ),
	debug_( false ),
	exclude_radius_( 20.0 )
{}

InteractingRotamerExplosion::~InteractingRotamerExplosion() {}

core::pack::task::operation::TaskOperationOP InteractingRotamerExplosion::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new InteractingRotamerExplosion( *this ) );
}

/// @brief restricts to repacking all residues outside of design_shell_ around each residue
void
InteractingRotamerExplosion::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task::operation;
	core::Real excl_radius_sq( exclude_radius_ * exclude_radius_ );
	utility::vector1< core::Size > target_seqpos;
	std::set< core::Size > parsed_seqpos( core::pose::get_resnum_list( string_target_seqpos_, pose ) );
	for (  std::set< core::Size >::const_iterator set_it( parsed_seqpos.begin() ), set_end( parsed_seqpos.end() ); set_it != set_end; ++set_it ) target_seqpos.push_back( *set_it );

	std::string report_string("Potential rotamer explosion positions against target pdb pos");
	for ( core::Size j(1); j <= target_seqpos.size(); ++j ) {
		report_string += (" " + utility::to_string( pose.pdb_info()->chain( target_seqpos[j] ) )  + utility::to_string( pose.pdb_info()->number( target_seqpos[j] ) ) );
	}
	report_string += " are: ";

	for ( core::Size i(1); i <= pose.total_residue(); ++i ) {

		if ( !pose.residue_type(i).is_protein() ) continue;
		bool out_of_all_target_res_radius( true );
		for ( core::Size j(1); j <= target_seqpos.size(); ++j ) {
			if ( (pose.residue(i).atom( pose.residue(i).nbr_atom() ).xyz().distance_squared( pose.residue( target_seqpos[j] ).atom( pose.residue( target_seqpos[j] ).nbr_atom() ).xyz() ) ) <= excl_radius_sq ) {
				out_of_all_target_res_radius = false;
				break;
			}
		}
		if ( out_of_all_target_res_radius ) continue;

		rotamer_set_operations::AddGood2BPairEnergyRotamersOP rotsetop( new rotamer_set_operations::AddGood2BPairEnergyRotamers(i, ex_level_, target_seqpos, score_cutoff_, false ) );
		if ( debug_ ) rotsetop->set_debug( true );
		task.nonconst_residue_task(i).append_rotamerset_operation( rotsetop );
		report_string += (utility::to_string( i ) + ", ");
	} //loop over pose res i

	tr << report_string << std::endl;
}


void
InteractingRotamerExplosion::parse_tag( TagCOP tag , DataMap & )
{
	string_target_seqpos_ = tag->getOption< std::string >( "target_seqpos" );// these are kept in memory until the pose is available (at apply time)
	score_cutoff_ =  tag->getOption< core::Real >( "score_cutoff", -0.5 );
	ex_level_ = tag->getOption< core::Size >( "ex_level", 4 );
	debug_ = tag->getOption< bool >( "debug", 0 );
	exclude_radius_ =  tag->getOption< core::Real >( "exclude_radius", 20.0 );
}

void InteractingRotamerExplosion::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	activate_common_simple_type( xsd, "non_negative_integer" );

	attributes.push_back( XMLSchemaAttribute::required_attribute( "target_seqpos", xs_string ) );
	attributes.push_back( XMLSchemaAttribute( "score_cutoff", xs_decimal, "-0.5" ) );
	attributes.push_back( XMLSchemaAttribute( "ex_level", "non_negative_integer", "4" ) );
	attributes.push_back( XMLSchemaAttribute( "debug", xs_boolean, "false" ) );
	attributes.push_back( XMLSchemaAttribute( "exclude_radius", xs_decimal, "20.0" ) );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

core::pack::task::operation::TaskOperationOP
InteractingRotamerExplosionCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new InteractingRotamerExplosion );
}

void InteractingRotamerExplosionCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InteractingRotamerExplosion::provide_xml_schema( xsd );
}

std::string InteractingRotamerExplosionCreator::keyname() const
{
	return InteractingRotamerExplosion::keyname();
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
