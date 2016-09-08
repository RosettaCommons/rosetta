// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/DsspDesignOperation.cc
/// @brief  Design residues with selected amino acids depending on DSSP secondary structure
/// assignment. All functionality here is included in the LayerDesign task operation, but
/// this filter has significantly reduced overhead by avoiding slow SASA calculations.
/// @author Brian Koepnick (koepnick@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/DsspDesignOperation.hh>
#include <protocols/toolbox/task_operations/DsspDesignOperationCreator.hh>

// Core Headers
#include <core/conformation/Conformation.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/selection.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/types.hh>

// Protocol Headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/assign/list_inserter.hpp>
#include <boost/foreach.hpp>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <utility/vector1.hh>

// C++ Headers
#include <set>

using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.DsspDesignOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

DsspDesignOperation::DsspDesignOperation()
{
	set_default_sse_residues();
}

DsspDesignOperation::DsspDesignOperation( DsspDesignOperation const & rval ): parent( rval ),
	sse_residues_( rval.sse_residues_ ),
	blueprint_( rval.blueprint_ )
{
}

DsspDesignOperation::~DsspDesignOperation() {}

TaskOperationOP
DsspDesignOperation::clone() const
{
	return TaskOperationOP( new DsspDesignOperation( *this ) );
}

void
DsspDesignOperation::set_blueprint( BluePrintOP const bp )
{
	blueprint_ = bp;
}

// default residue sets are derived from LayerDesign defaults
// (e.g. 'Helix' = union of 'Helix' allowed residues from all layers)
void
DsspDesignOperation::set_default_sse_residues() {
	TR << "Initializing DSSP regions with default residues" << std::endl;
	boost::assign::insert( sse_residues_ )
		( std::string( "Loop" ),   std::string( "ACDEFGHIKLMNPQRSTVWY" ) )
		( std::string( "Strand" ),   std::string( "DEFHIKLNQRSTVWY" ) )
		( std::string( "Helix" ),   std::string( "ADEFIKLNQRSTVWY" ) )
		( std::string( "HelixStart" ),  std::string( "ADEFHIKLNPQRSTVWY" ) )
		( std::string( "HelixCapping" ), std::string( "DNST" ) )
		( std::string( "Nterm" ),    std::string( "ACDEFGHIKLMNPQRSTVWY" ) )
		( std::string( "Cterm" ),    std::string( "ACDEFGHIKLMNPQRSTVWY" ) );
}

void
DsspDesignOperation::set_restrictions_aa( std::string const & sse, std::string const & aas )
{
	// handle all SSEs
	if ( sse == "all" ) {
		for ( SecStructResidues::iterator it = sse_residues_.begin(), end = sse_residues_.end(); it != end; ++it ) {
			set_restrictions_aa( it->first, aas );
		}
		// check that SSE is valid
	} else if ( sse_residues_.find( sse ) == sse_residues_.end() ) {
		utility_exit_with_message( "Invalid SSE: " + sse + ", valid SSEs are Helix, Strand, Loop, HelixCapping, HelixStart, Nterm, Cterm, or all" );

	} else {
		sse_residues_[ sse ] = aas;
	}
}

void
DsspDesignOperation::set_restrictions_append( std::string const & sse, std::string const & aas )
{
	// handle all SSEs
	if ( sse == "all" ) {
		for ( SecStructResidues::iterator it = sse_residues_.begin(), end = sse_residues_.end(); it != end; ++it ) {
			set_restrictions_append( it->first, aas );
		}
		// check that SSE is valid
	} else if ( sse_residues_.find( sse ) == sse_residues_.end() ) {
		utility_exit_with_message( "Invalid SSE: " + sse + ", valid SSEs are Helix, Strand, Loop, HelixCapping, HelixStart, Nterm, Cterm, or all" );

	} else {
		const std::string sse_res = sse_residues_[ sse ];
		std::set< char > temp_def_res_set( sse_res.begin(), sse_res.end() );
		temp_def_res_set.insert( aas.begin(), aas.end() );
		sse_residues_[ sse ] = std::string( temp_def_res_set.begin(), temp_def_res_set.end() );
	}
}

void
DsspDesignOperation::set_restrictions_exclude( std::string const & sse, std::string const & aas )
{
	// handle all SSEs
	if ( sse == "all" ) {
		for ( SecStructResidues::iterator it = sse_residues_.begin(), end = sse_residues_.end(); it != end; ++it ) {
			set_restrictions_exclude( it->first, aas );
		}
		// check that SSE is valid
	} else if ( sse_residues_.find( sse ) == sse_residues_.end() ) {
		utility_exit_with_message( "Invalid SSE: " + sse + ", valid SSEs are Helix, Strand, Loop, HelixCapping, HelixStart, Nterm, Cterm, or all" );

	} else {
		const std::string sse_res = sse_residues_[ sse ];
		std::set< char > temp_def_res_set( sse_res.begin(), sse_res.end() );
		BOOST_FOREACH ( char aa, aas ) {  temp_def_res_set.erase( aa ); }
		sse_residues_[ sse ] = std::string( temp_def_res_set.begin(), temp_def_res_set.end() );
	}
}

utility::vector1< bool >
DsspDesignOperation::get_restrictions( std::string const & ss_type ) const
{
	utility::vector1< bool > restrict_to_aa( core::chemical::num_canonical_aas, false );
	BOOST_FOREACH ( char restype, sse_residues_.find( ss_type )->second ) {
		restrict_to_aa[ core::chemical::aa_from_oneletter_code( restype ) ] = true;
	}

	return restrict_to_aa;
}

/// @brief
void
DsspDesignOperation::apply( core::pose::Pose const & input_pose, core::pack::task::PackerTask & task ) const
{
	core::pose::Pose pose;

	// symmetry check
	if ( core::pose::symmetry::is_symmetric( input_pose ) ) {
		TR << "Symmetry detected, extracting asymmetric unit." << std::endl;
		core::pose::symmetry::extract_asymmetric_unit( input_pose, pose, false );
	} else {
		pose = input_pose;
	}

	// support input from blueprint files
	std::string secstruct;
	if ( blueprint_ ) {
		TR << "Using blueprint..." << std::endl;
		secstruct = blueprint_->secstruct();
	} else {
		core::scoring::dssp::Dssp dssp( pose );
		dssp.dssp_reduced();
		secstruct = dssp.get_dssp_secstruct();
	}
	TR << "Secondary structure: " << secstruct << std::endl;

	// find HelixStart and HelixCapping positions ( copied from LayerDesign )
	bool flag( false );
	utility::vector1< bool > helix_capping( pose.size(), false );
	utility::vector1< bool > initial_helix( pose.size(), false );
	for ( Size i=1; i<=pose.size(); ++i ) {

		// ignore non-protein residues
		if ( ! pose.residue( i ).is_protein() ) continue;

		char ss( secstruct[ i-1 ] );
		if ( ss == 'H' && !flag && i != 1 ) {
			initial_helix[ i ] = true;
			helix_capping[ i-1 ] = true;
			flag = true;
		}

		if ( ss != 'H' && flag ) {
			flag = false;
		}
	}


	// apply restrictions to PackerTask
	for ( Size i=1; i<=pose.size(); ++i ) {

		// if residue is not a protein, continue
		// make repackable only?
		if ( ! pose.residue( i ).is_protein() ) {
			//task.nonconst_residue_task( i ).restrict_to_repacking();
			continue;
		}

		char ss( secstruct[ i-1 ] );

		// handle termini because DSSP always assigns them as loop
		if ( pose.residue( i ).is_lower_terminus() ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "Nterm" ) );
		} else if ( pose.residue( i ).is_upper_terminus() ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "Cterm" ) );

		} else if ( helix_capping[ i ] == true ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "HelixCapping" ) );
		} else if ( initial_helix[ i ] == true ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "HelixStart" ) );
		} else if ( ss == 'E' ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "Strand" ) );
		} else if ( ss == 'L' ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "Loop" ) );
		} else if ( ss == 'H' ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( "Helix" ) );
		}
	}
}

void
DsspDesignOperation::parse_tag( TagCOP tag , DataMap & )
{
	// get secondary structure definition from blueprint file
	if ( tag->hasOption( "blueprint" ) ) {
		blueprint_ = BluePrintOP( new BluePrint( tag->getOption< std::string >( "blueprint" ) ) );
	}

	BOOST_FOREACH ( utility::tag::TagCOP const sse_tag, tag->getTags() ) {
		const std::string sse = sse_tag->getName(); // Helix, Strand, Loop, HelixCapping, HelixStart, Nterm, Cterm

		// check that SSE is valid
		if ( sse != "all" && ( sse_residues_.find( sse ) == sse_residues_.end() ) ) {
			utility_exit_with_message( "Invalid SSE: " + sse + ", valid SSEs are Helix, Strand, Loop, HelixCapping, HelixStart, Nterm, Cterm, or all" );
		}

		// explicitly define allowed residues
		if ( sse_tag->hasOption( "aa" ) ) {
			const std::string aas = sse_tag->getOption< std::string >( "aa" );
			TR << "Assigning residues " << aas << " to " << sse << std::endl;
			set_restrictions_aa( sse, aas );
		}

		// append to allowed residues
		if ( sse_tag->hasOption( "append" ) ) {
			const std::string aas = sse_tag->getOption< std::string >( "append" );
			TR << "Appending residues " << aas << " to " << sse << std::endl;
			set_restrictions_append( sse, aas );
		}

		// exclude from allowed residues
		if ( sse_tag->hasOption( "exclude" ) ) {
			const std::string aas = sse_tag->getOption< std::string >( "exclude" );
			TR << "Excluding residues " << aas << " from " << sse << std::endl;
			set_restrictions_exclude( sse, aas );
		}
	}
}

std::string dsspdo_subelement_ct_name( std::string const & name ) {
	return "dsspdo_subelement_" + name + "Type";
}

std::string dsspdo_group_name() { return "dsspdo_subelement"; }

void DsspDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// attributes for the subelements -- all the subelements have the same attributes
	AttributeList dsspdo_subtag_attributes;
	dsspdo_subtag_attributes
		+ XMLSchemaAttribute( "aa", xs_string )
		+ XMLSchemaAttribute( "append", xs_string )
		+ XMLSchemaAttribute( "exclude", xs_string );

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & dsspdo_subelement_ct_name );
	subelements.add_simple_subelement( "Loop",         dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "Strand",       dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "Helix",        dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "HelixStart",   dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "HelixCapping", dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "Nterm",        dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "Cterm",        dsspdo_subtag_attributes );
	subelements.add_simple_subelement( "all",          dsspdo_subtag_attributes );

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "name", xs_string )
		+ XMLSchemaAttribute( "blueprint", xs_string );

	XMLSchemaComplexTypeGenerator complex_type_generator;
	complex_type_generator
		.element_name( keyname() )
		.complex_type_naming_func( & complex_type_name_for_task_op )
		.add_attributes( attributes )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}

TaskOperationOP
DsspDesignOperationCreator::create_task_operation() const
{
	return TaskOperationOP( new DsspDesignOperation );
}

void DsspDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DsspDesignOperation::provide_xml_schema( xsd );
}

std::string DsspDesignOperationCreator::keyname() const
{
	return DsspDesignOperation::keyname();
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
