// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CopyRotamerMover.cc
/// @brief A mover to copy a rotamer (residue identity and conformation) from one position in a pose to another.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

// Unit headers
#include <protocols/simple_moves/CopyRotamerMover.hh>
#include <protocols/simple_moves/CopyRotamerMoverCreator.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.CopyRotamerMover" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CopyRotamerMover::CopyRotamerMover():
	protocols::moves::Mover( CopyRotamerMover::mover_name() ),
	template_res_index_(0),
	target_res_index_(0),
	copy_identity_(true),
	copy_torsions_(true)
	//copy_bondangles_(true),
	//copy_bondlengths_(true)
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CopyRotamerMover::CopyRotamerMover(
	CopyRotamerMover const & src
):
	protocols::moves::Mover( src ),
	template_res_index_(src.template_res_index_),
	target_res_index_(src.target_res_index_),
	copy_identity_(src.copy_identity_),
	copy_torsions_(src.copy_torsions_)
	//copy_bondangles_(src.copy_bondangles_),
	//copy_bondlengths_(src.copy_bondlengths_)
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CopyRotamerMover::~CopyRotamerMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
CopyRotamerMover::apply(
	core::pose::Pose &pose
) {
	runtime_assert_string_msg( template_res_index_ > 0 && template_res_index_ <=pose.total_residue(), "Error in protocols::simple_moves::CopyRotamerMover::apply(): The template residue index is out of range (less than one or greater than the number of residues in the pose)." );
	runtime_assert_string_msg( target_res_index_ > 0 && target_res_index_ <=pose.total_residue(), "Error in protocols::simple_moves::CopyRotamerMover::apply(): The target residue index is out of range (less than one or greater than the number of residues in the pose)." );

	//First, copy residue identity:
	if ( copy_identity_ ) {
		MutateResidue mutres( target_res_index_, pose.residue_type(template_res_index_).name() ); //A mover to change the residue identity.
		mutres.apply( pose );
	}

	//Next, set side-chain torsions
	if ( copy_torsions_ ) {
		for ( core::Size i=1, imax=pose.residue_type(template_res_index_).nchi(); i<=imax; ++i ) {
			if ( i <= pose.residue_type(target_res_index_).nchi() ) {
				pose.set_chi( i, target_res_index_, pose.chi( i, template_res_index_ ) );
			} else {
				TR.Warning << "skipping chi " << i << ", which is not present in target residue." << std::endl;
			}
		}
	}

	//TODO -- add support for copying bond angles and bond lenghts.

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
CopyRotamerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
CopyRotamerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Parsing options for CopyRotamerMover:" << std::endl;

	runtime_assert_string_msg( tag->hasOption("template_res_index"), "Error in protocols::simple_moves::CopyRotamerMover::parse_my_tag(): The user MUST specify a template residue with the \"template_res_index=(&int)\" option." );
	runtime_assert_string_msg( tag->hasOption("target_res_index"), "Error in protocols::simple_moves::CopyRotamerMover::parse_my_tag(): The user MUST specify a target residue with the \"target_res_index=(&int)\" option." );

	set_template_res_index( tag->getOption< core::Size >( "template_res_index", template_res_index_ ) );
	set_target_res_index( tag->getOption< core::Size >( "target_res_index", template_res_index_ ) );
	set_copy_identity( tag->getOption<bool>( "copy_identity", copy_identity_ ) );
	set_copy_torsions( tag->getOption<bool>( "copy_torsions", copy_torsions_ ) );
	//set_copy_bondangles( tag->getOption<bool>( "copy_bondangles", copy_bondangles_ ) );
	//set_copy_bondlengths( tag->getOption<bool>( "copy_bondlengths", copy_bondlengths_ ) );

	TR << "Finished parsing options for CopyRotamerMover." << std::endl;
	TR.flush();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
CopyRotamerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new CopyRotamerMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CopyRotamerMover::clone() const
{
	return protocols::moves::MoverOP( new CopyRotamerMover( *this ) );
}

///////////////////////////////
/// Setters          ///
///////////////////////////////

/// @brief Set pose index of the residue FROM which we're copying.
///
void
CopyRotamerMover::set_template_res_index(
	core::Size const index_in
) {
	runtime_assert_string_msg( index_in > 0, "Error in protocols::simple_moves::CopyRotamerMover::set_template_res_index(): The index of the template residue must be greater than zero." );
	template_res_index_ = index_in;
	TR << "Set template residue index to " << index_in << "." << std::endl;
}

/// @brief Set pose index of the residue TO which we're copying.
///
void
CopyRotamerMover::set_target_res_index(
	core::Size const index_in
) {
	runtime_assert_string_msg( index_in > 0, "Error in protocols::simple_moves::CopyRotamerMover::set_target_res_index(): The index of the target residue must be greater than zero." );
	target_res_index_ = index_in;
	TR << "Set target residue index to " << index_in << "." << std::endl;
}

/// @brief Set whether we're copying the residue identity.
///
void
CopyRotamerMover::set_copy_identity(
	bool const setting_in
) {
	copy_identity_ = setting_in;
	if ( setting_in ) {
		TR << "Set mover to copy residue identities." << std::endl;
	} else {
		TR << "Set mover NOT to copy residue identities.  Note that this could cause aberrant behaviour if side-chain DoFs don't match!" << std::endl;
	}
}

/// @brief Set whether we're copying the residue side-chain torsions.
///
void
CopyRotamerMover::set_copy_torsions(
	bool const setting_in
) {
	copy_torsions_ = setting_in;
	if ( setting_in ) {
		TR << "Set mover to copy residue side-chain torsions." << std::endl;
	} else {
		TR << "Set mover NOT to copy residue side-chain torsions." << std::endl;
	}
}

// / @brief Set whether we're copying the residue side-chain bond angles.
// /
/*void
CopyRotamerMover::set_copy_bondangles(
bool const setting_in
) {
copy_bondangles = setting_in;
if( setting_in ) {
TR << "Set mover to copy residue side-chain bond angles." << std::endl;
} else {
TR << "Set mover NOT to copy residue side-chain bond angles." << std::endl;
}
}*/

// / @brief Set whether we're copying the residue side-chain bond lengths.
// /
/*void
CopyRotamerMover::set_copy_bondlengths(
bool const setting_in
) {
copy_bondlengths = setting_in;
if( setting_in ) {
TR << "Set mover to copy residue side-chain bond lengths." << std::endl;
} else {
TR << "Set mover NOT to copy residue side-chain bond lengths." << std::endl;
}
}*/

/// @brief Get the name of the Mover
// XRW TEMP std::string
// XRW TEMP CopyRotamerMover::get_name() const
// XRW TEMP {
// XRW TEMP  return CopyRotamerMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CopyRotamerMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "CopyRotamer";
// XRW TEMP }

std::ostream &
operator<<( std::ostream & os, CopyRotamerMover const & mover )
{
	mover.show(os);
	return os;
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CopyRotamerMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new CopyRotamerMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CopyRotamerMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CopyRotamerMover::mover_name();
// XRW TEMP }

std::string CopyRotamerMover::get_name() const {
	return mover_name();
}

std::string CopyRotamerMover::mover_name() {
	return "CopyRotamer";
}

void CopyRotamerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute(
		"template_res_index", xsct_non_negative_integer,
		"The index, in Rosetta pose numbering, of the residue from which the "
		"side-chain will be copied. This residue is not altered by this operation." )
		+ XMLSchemaAttribute::required_attribute(
		"target_res_index", xsct_non_negative_integer,
		"The index, in Rosetta pose numbering, of the residue to which the "
		"side-chain will be copied. The identity and/or conformation of this "
		"residue's sidechain is altered by this operation." )
		+ XMLSchemaAttribute(
		"copy_identity", xsct_rosetta_bool,
		"Should the identity of the template residue by copied to the target? "
		"Default true. If false, only side-chain torsion values will be copied. "
		"This can create strange results if the template and target residues "
		"have different numbers of side-chain chi angles, or if they have "
		"significantly different side-chain structures." )
		+ XMLSchemaAttribute(
		"copy_torsions", xsct_rosetta_bool,
		"Should the side-chain dihedral values of the template residue be "
		"copied to the target? Default true." );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"This is a very simple mover that copies the residue identity and/or sidechain "
		"conformation from one residue in a pose (the template) to another (the target).",
		attlist );
}

std::string CopyRotamerMoverCreator::keyname() const {
	return CopyRotamerMover::mover_name();
}

protocols::moves::MoverOP
CopyRotamerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new CopyRotamerMover );
}

void CopyRotamerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CopyRotamerMover::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

} //protocols
} //simple_moves

