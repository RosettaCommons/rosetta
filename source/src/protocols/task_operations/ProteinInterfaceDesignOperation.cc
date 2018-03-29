// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/ProteinInterfaceDesignOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/task_operations/ProteinInterfaceDesignOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/task_operations/DesignAroundOperation.hh>
// C++ Headers
#include <set>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.TaskOperations.ProteinInterfaceDesignOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
ProteinInterfaceDesignOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ProteinInterfaceDesignOperation );
}

void ProteinInterfaceDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinInterfaceDesignOperation::provide_xml_schema( xsd );
}

std::string ProteinInterfaceDesignOperationCreator::keyname() const
{
	return ProteinInterfaceDesignOperation::keyname();
}

ProteinInterfaceDesignOperation::ProteinInterfaceDesignOperation() :
	repack_chain1_( true ),
	repack_chain2_( true ),
	design_chain1_( false ),
	design_chain2_( true ),
	allow_all_aas_( false ),
	design_all_aas_( false ),
	interface_distance_cutoff_( 8.0 ),
	jump_( 1 ),
	modify_before_jump_( true ),
	modify_after_jump_ ( true )
{}

ProteinInterfaceDesignOperation::~ProteinInterfaceDesignOperation() = default;

core::pack::task::operation::TaskOperationOP ProteinInterfaceDesignOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ProteinInterfaceDesignOperation( *this ) );
}

/// @brief the default taskoperation for protein-interface design. Sets up which chains to repack/design
/// disable disulfides, prolines and glycines from design, restrict designable positions to all but pro/gly/cys
void
ProteinInterfaceDesignOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task::operation;
	using namespace protocols::task_operations;

	core::Size const chains( pose.conformation().num_chains() );
	if ( chains<2 ) {
		TR << "WARNING WARNING WARNING!!! I'm calling ProteinInterfaceDesignOperation but I have only chain in pose. No interface to design!!" << std::endl;
	}
	DesignAroundOperation dao1, dao2;
	dao1.design_shell( interface_distance_cutoff_ ); dao2.design_shell( interface_distance_cutoff_ );
	for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
		if ( resi <= pose.conformation().chain_end( jump() ) ) {
			dao1.include_residue( resi );
		} else {
			dao2.include_residue( resi );
		}
	}
	if ( modify_after_jump() ) {
		dao1.apply( pose, task );
	}
	if ( modify_before_jump() ) {
		dao2.apply( pose, task );
	}

	for ( core::Size chain=1; chain<=chains; ++chain ) {
		PreventChainFromRepackingOperation pcfr;
		if ( chain<=jump() && !repack_chain1_ && modify_before_jump() ) {
			pcfr.chain( chain );
			pcfr.apply( pose, task );
		}
		if ( chain > jump() && !repack_chain2_ && modify_after_jump() ) {
			pcfr.chain( chain );
			pcfr.apply( pose, task );
		}
		RestrictChainToRepackingOperation rctr;
		if ( chain<=jump() && !design_chain1_ && modify_before_jump() ) {
			rctr.chain( chain );
			rctr.apply( pose, task );
		}
		if ( chain > jump() && !design_chain2_ && modify_after_jump() ) {
			rctr.chain( chain );
			rctr.apply( pose, task );
		}

	}
	NoRepackDisulfides nrd;
	nrd.apply( pose, task );
	utility::vector1< core::Size > residues;
	using namespace core::chemical;
	utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, true );
	// check if we are allowing design of all cannonical aas
	if ( !allow_all_aas_ ) {
		allowed_aas[ aa_cys ] = false;
		allowed_aas[ aa_gly ] = false;
		allowed_aas[ aa_pro ] = false;
	}
	for ( core::Size i = 1; i<=pose.size(); ++i ) {
		utility::vector1< bool > specific_allowed_aas( allowed_aas );
		if ( pose.residue( i ).aa() == aa_cys ) specific_allowed_aas[ aa_cys ] = true; // allow a native cystein in packer task
		RestrictAbsentCanonicalAAS racaas( i, specific_allowed_aas );
		racaas.apply( pose, task );
		//check if we are designing all residues
		if ( !design_all_aas_ ) {
			if ( pose.residue( i ).aa() == aa_pro || pose.residue( i ).aa() == aa_gly ) {
				residues.push_back( i );
			}
		}
	}
	// needed because OperateOnCertainResidues will act on all residues in
	// the pose if the vector of indices is empty
	if ( residues.size() != 0 ) {
		OperateOnCertainResidues oocr;
		oocr.op( ResLvlTaskOperationCOP( new PreventRepackingRLT ) );
		oocr.residue_indices( residues );
		oocr.apply( pose, task );
	}
}

void
ProteinInterfaceDesignOperation::repack_chain1( bool const repack )
{
	repack_chain1_ = repack;
}
void
ProteinInterfaceDesignOperation::repack_chain2( bool const repack )
{
	repack_chain2_ = repack;
}
void
ProteinInterfaceDesignOperation::design_chain1( bool const design )
{
	design_chain1_ = design;
}
void
ProteinInterfaceDesignOperation::design_chain2( bool const design )
{
	design_chain2_ = design;
}

void
ProteinInterfaceDesignOperation::allow_all_aas( bool const allow )
{
	allow_all_aas_ = allow;
}
void
ProteinInterfaceDesignOperation::design_all_aas( bool const design_all )
{
	design_all_aas_ = design_all;
}
void
ProteinInterfaceDesignOperation::interface_distance_cutoff( core::Real const dist )
{
	interface_distance_cutoff_ = dist;
}
void
ProteinInterfaceDesignOperation::jump( core::Size const j )
{
	jump_ = j;
}
core::Size
ProteinInterfaceDesignOperation::jump() const
{
	return( jump_ );
}

void
ProteinInterfaceDesignOperation::parse_tag( TagCOP tag , DataMap & )
{
	repack_chain1( tag->getOption< core::Size >( "repack_chain1", 1 ) );
	repack_chain2( tag->getOption< core::Size >( "repack_chain2", 1 ) );
	design_chain1( tag->getOption< core::Size >( "design_chain1", 0 ) );
	design_chain2( tag->getOption< core::Size >( "design_chain2", 1 ) );
	allow_all_aas( tag->getOption< core::Size >( "allow_all_aas", 0 ) );
	design_all_aas( tag->getOption< core::Size >( "design_all_aas", 0 ) );
	jump( tag->getOption< core::Size >( "jump", 1 ) );
	interface_distance_cutoff( tag->getOption< core::Real >( "interface_distance_cutoff", 8.0 ) );
	modify_before_jump( tag->getOption< bool >( "modify_before_jump", true ) );
	modify_after_jump(  tag->getOption< bool >( "modify_after_jump",  true ) );
}

void ProteinInterfaceDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"repack_chain1", xsct_non_negative_integer,
		"first chain to repack",
		"1"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"repack_chain2", xsct_non_negative_integer,
		"2ndchain to repack",
		"1"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"design_chain1", xsct_non_negative_integer,
		"first chain to design",
		"0"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"design_chain2", xsct_non_negative_integer,
		"2nd chain to design",
		"1"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"allow_all_aas", xsct_non_negative_integer,
		"allow all amino acid types for design. default: do not design to G/C/P",
		"0"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"design_all_aas", xsct_non_negative_integer,
		"design all amino acids. default: do not design G/C/P",
		"0"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"jump", xsct_non_negative_integer,
		"residue with the jump",
		"1"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"interface_distance_cutoff", xsct_real,
		"Cbeta distance of chain where design/repack is allowed on other chain",
		"8.0"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"modify_before_jump", xsct_rosetta_bool,
		"change residues before the jump. For instance, if you want set repack on "
		"chain2 interfacial residues to true, and the rest of chain2 to false, "
		"and yet not change the task for chain1, "
		"then use this taskoperation with modify_before_jump=0",
		"true"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"modify_after_jump", xsct_rosetta_bool,
		"change residues after the jump",
		"true"  );

	task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Restricts to the task that is the basis for protein-interface design.");
}


} //namespace protocols
} //namespace task_operations
