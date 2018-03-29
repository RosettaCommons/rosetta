// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictIdentitiesAtAlignedPositionsOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <protocols/task_operations/RestrictIdentitiesAtAlignedPositions.hh>
#include <protocols/task_operations/RestrictIdentitiesAtAlignedPositionsCreator.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/task_operations/DesignAroundOperation.hh>
#include <core/pack/task/TaskFactory.hh>
// Project Headers
#include <core/pose/Pose.hh>
#include <utility/string_util.hh>


// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>


// C++ Headers

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.TaskOperations.RestrictIdentitiesAtAlignedPositionsOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using namespace std;

core::pack::task::operation::TaskOperationOP
RestrictIdentitiesAtAlignedPositionsOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictIdentitiesAtAlignedPositionsOperation );
}

void RestrictIdentitiesAtAlignedPositionsOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictIdentitiesAtAlignedPositionsOperation::provide_xml_schema( xsd );
}

std::string RestrictIdentitiesAtAlignedPositionsOperationCreator::keyname() const
{
	return RestrictIdentitiesAtAlignedPositionsOperation::keyname();
}

RestrictIdentitiesAtAlignedPositionsOperation::RestrictIdentitiesAtAlignedPositionsOperation() :
	RestrictOperationsBase(),
	chain_( 1 ),
	design_only_target_residues_( false ),
	prevent_repacking_( false ),
	keep_aas_( "ACDEFGHIKLMNPQRSTVWY" ),
	restrict_identities_( false )
{
	source_pose_ = core::pose::PoseOP( new core::pose::Pose );
	res_ids_.clear();
}

RestrictIdentitiesAtAlignedPositionsOperation::~RestrictIdentitiesAtAlignedPositionsOperation() = default;

core::pack::task::operation::TaskOperationOP RestrictIdentitiesAtAlignedPositionsOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictIdentitiesAtAlignedPositionsOperation( *this ) );
}

void
RestrictIdentitiesAtAlignedPositionsOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace protocols::rosetta_scripts;
	using namespace core::pack::task::operation;

	DesignAroundOperation dao;
	dao.design_shell( 0.01 );
	dao.repack_shell( 6.0 );
	for ( core::Size const resid : res_ids_ ) {
		core::Size const nearest_to_res = find_nearest_res( pose, *source_pose_, resid, chain() );
		if ( nearest_to_res == 0 ) {
			TR.Warning << "could not find a residue near to "<<resid<<std::endl;
			continue;
		}//fi
		RestrictAbsentCanonicalAASRLTOP racaas1( new RestrictAbsentCanonicalAASRLT ); /// used to determine the single residue identity taken from the source pose
		RestrictAbsentCanonicalAASRLTOP racaas2( new RestrictAbsentCanonicalAASRLT ); /// used to limit identities to those specified by the user
		PreventRepackingRLTOP pr( new PreventRepackingRLT );
		char const residue_id( source_pose_->residue( resid ).name1() );
		std::string residues_to_keep("");
		residues_to_keep += residue_id;
		racaas1->aas_to_keep( residues_to_keep );
		racaas2->aas_to_keep( keep_aas_ );
		OperateOnCertainResidues oocr;
		if ( prevent_repacking() && source_pose_->residue( resid ).name1() == pose.residue( nearest_to_res ).name1() ) { /// if the source and designed pose have the same residue identity we can additionally prevent repacking at this position
			oocr.op( pr );
		} else if ( restrict_identities() ) {
			oocr.op( racaas2 );
		} else {
			oocr.op( racaas1 );
		}
		utility::vector1< core::Size > temp_vec;
		temp_vec.clear();
		temp_vec.push_back( nearest_to_res );
		oocr.residue_indices( temp_vec );
		oocr.apply( pose, task );
		dao.include_residue( nearest_to_res );
	}//foreach resid
	if ( design_only_target_residues() ) {
		dao.apply( pose, task );
	}
}

void
RestrictIdentitiesAtAlignedPositionsOperation::source_pose( std::string const & s ) {
	core::import_pose::pose_from_file( *source_pose_, s , core::import_pose::PDB_file);
}

void
RestrictIdentitiesAtAlignedPositionsOperation::parse_tag( TagCOP tag , DataMap & )
{
	using namespace protocols::rosetta_scripts;
	utility::vector1< std::string > pdb_names, start_res, stop_res;
	source_pose( tag->getOption< std::string >( "source_pdb" ) );
	std::string const res_list( tag->getOption< std::string >( "resnums" ) );
	utility::vector1< std::string > const split_reslist( utility::string_split( res_list,',' ) );
	chain( tag->getOption< core::Size >( "chain", 1 ) );
	TR<<"source_pdb: "<<tag->getOption< std::string >( "source_pdb" )<<" restricting residues: ";
	for ( std::string const & res_str : split_reslist ) {
		res_ids_.push_back( core::pose::parse_resnum( res_str, *source_pose_ ) );
		TR<<res_str<<",";
	}
	design_only_target_residues( tag->getOption< bool >( "design_only_target_residues", false ) );
	prevent_repacking( tag->getOption< bool >( "prevent_repacking", false ) );
	keep_aas_ = tag->getOption< std::string >( "keep_aas", "ACDEFGHIKLMNPQRSTVWY" );
	restrict_identities( tag->hasOption( "keep_aas" ) );
	TR<<std::endl;
}

void RestrictIdentitiesAtAlignedPositionsOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::required_attribute( "source_pdb", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "resnums", xsct_int_cslist , "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default(  "chain", xsct_non_negative_integer, "XRW TO DO",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "design_only_target_residues", xsct_rosetta_bool, "XRW TO DO",  "false"  )
		+ XMLSchemaAttribute::attribute_w_default(  "prevent_repacking", xsct_rosetta_bool, "XRW TO DO",  "false"  )
		+ XMLSchemaAttribute::attribute_w_default(  "keep_aas", xs_string, "XRW TO DO",  "ACDEFGHIKLMNPQRSTVWY"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "XRW TO DO" );
}


} //namespace protocols
} //namespace task_operations
