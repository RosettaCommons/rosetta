// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/pHVariantTaskOperation.cc
/// @brief A TaskOp that enables sampling of protonation variants during repacking. The ionizable sidechains considered are ASP, GLU, LYS, TYR, and HIS. All other posiitons are restricted to repacking. The -pH_mode commandline option is still currently required to initially read the variants into the database.
/// @author Rebecca Alford (rfalford12@gmail.com)

#include <protocols/task_operations/pHVariantTaskOperation.hh>
#include <protocols/task_operations/pHVariantTaskOperationCreator.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

static basic::Tracer TR( "protocols.task_operations.pHVariantTaskOperation" );

namespace protocols {
namespace task_operations {

using namespace core::select::residue_selector;
using namespace core::pack::task;
using namespace core::pack::task::operation;

pHVariantTaskOperation::pHVariantTaskOperation():
	TaskOperation(),
	select_asp_( new ResidueNameSelector( "ASP" /*aspartate*/ ) ),
	select_glu_( new ResidueNameSelector( "GLU" /*glutamate*/ ) ),
	select_tyr_( new ResidueNameSelector( "TYR" /*tyrosine*/ ) ),
	select_lys_( new ResidueNameSelector( "LYS" /*lysine*/ ) ),
	select_his_( new ResidueNameSelector( "HIS" /*histidine*/ ) )
{}

pHVariantTaskOperation::~pHVariantTaskOperation() {}

TaskOperationOP
pHVariantTaskOperation::clone() const {
	return TaskOperationOP( new pHVariantTaskOperation( *this ) );
}

pHVariantTaskOperation::pHVariantTaskOperation( pHVariantTaskOperation const & src ):
	TaskOperation(src),
	select_asp_( src.select_asp_ ),
	select_glu_( src.select_glu_ ),
	select_tyr_( src.select_tyr_ ),
	select_lys_( src.select_lys_ ),
	select_his_( src.select_his_ )
{}

void
pHVariantTaskOperation::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap& )
{}

void
pHVariantTaskOperation::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task ) const {

	using namespace basic::options;

	if ( !option[ OptionKeys::pH::pH_mode ].value() ) {
		utility_exit_with_message( "pH variant task operation will not work without -pH_mode!" );
	}

	// For each residue type, create a selection and allow protonation variants during packing
	/// Aspartate (ASP, D)
	utility::vector1< std::string > allowed_asp;
	allowed_asp.push_back( "ASP" );
	allowed_asp.push_back( "ASP_P1" );
	allowed_asp.push_back( "ASP_P2" );
	utility::vector1< bool > asp_residues = select_asp_->apply( pose );

	/// Glutamate (GLU, E)
	utility::vector1< std::string > allowed_glu;
	allowed_glu.push_back( "GLU" );
	allowed_glu.push_back( "GLU_P1" );
	allowed_glu.push_back( "GLU_P2" );
	utility::vector1< bool > glu_residues = select_glu_->apply( pose );

	/// Tyrosine (TYR, Y)
	utility::vector1< std::string > allowed_tyr;
	allowed_tyr.push_back( "TYR" );
	allowed_tyr.push_back( "TYR_D" );
	utility::vector1< bool > tyr_residues = select_tyr_->apply( pose );

	/// Histidine (HIS, H)
	utility::vector1< std::string > allowed_his;
	allowed_his.push_back( "HIS" );
	allowed_his.push_back( "HIS_D" );
	allowed_his.push_back( "HIS_P" );
	utility::vector1< bool > his_residues = select_his_->apply( pose );

	/// Lysine (LYS, K)
	utility::vector1< std::string > allowed_lys;
	allowed_lys.push_back( "LYS" );
	allowed_lys.push_back( "LYS_D" );
	utility::vector1< bool > lys_residues = select_lys_->apply( pose );

	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( asp_residues[ii] ) {
			task.nonconst_residue_task(ii).restrict_restypes( allowed_asp );
		} else if ( glu_residues[ii] ) {
			task.nonconst_residue_task(ii).restrict_restypes( allowed_glu );
		} else if ( tyr_residues[ii] ) {
			task.nonconst_residue_task(ii).restrict_restypes( allowed_tyr );
		} else if ( his_residues[ii] ) {
			task.nonconst_residue_task(ii).restrict_restypes( allowed_his );
		} else if ( lys_residues[ii] ) {
			task.nonconst_residue_task(ii).restrict_restypes( allowed_lys );
		} else {
			task.nonconst_residue_task(ii).restrict_to_repacking();
		}
	}
}

std::string
pHVariantTaskOperation::keyname() {
	return "pHVariantTaskOperation";
}

void
pHVariantTaskOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//attlist + XMLSchemaAttribute("motif", xs_string, motif_str);

	task_op_schema_w_attributes( xsd, keyname(), attlist,
		"Author: Rebecca F. Alford (rfalford12@gmail.com)\n"
		"A TaskOp that enables sampling of protonation variants during repacking. The ionizable sidechains considered are ASP, GLU, LYS, TYR, and HIS. All other posiitons are restricted to repacking. The -pH_mode commandline option is still currently required to initially read the variants into the database." );
}

TaskOperationOP
pHVariantTaskOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new pHVariantTaskOperation );
}

std::string
pHVariantTaskOperationCreator::keyname() const
{
	return pHVariantTaskOperation::keyname();
}

void
pHVariantTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	pHVariantTaskOperation::provide_xml_schema( xsd );
}


} //protocols
} //task_operations
