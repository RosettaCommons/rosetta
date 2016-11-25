// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @authors Vikram Mulligan (vmullig@uw.edu), Parisa Hosseinzadeh (parisah@uw.edu)

//Unit Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>
#include <protocols/cyclic_peptide/FlipChiralityMoverCreator.hh>

//Core Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/chemical/VariantType.hh>

//Protocol Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

//Basic Headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.FlipChiralityMover" );

namespace protocols {
namespace cyclic_peptide {

// XRW TEMP std::string
// XRW TEMP FlipChiralityMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FlipChiralityMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FlipChiralityMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FlipChiralityMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FlipChiralityMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "FlipChiralityMover";
// XRW TEMP }

/// @brief Constructor for FlipChirality mover.
///
FlipChiralityMover::FlipChiralityMover() :
	Mover("FlipChiralityMover"),
	selector_()
{
	center_assigned_ = false;
	normal_assigned_ = false;
	center_.zero();
	normal_.zero();
	//set_rosetta_scripts_tag( utility::tag::TagOP( new utility::tag::Tag() ) ); //do I need this?
}

/// @brief Copy constructor for FlipChiralityMover mover.
///
FlipChiralityMover::FlipChiralityMover( FlipChiralityMover const &src ) :
	Mover("FlipChiralityMover"),
	selector_( src.selector_ )

{
	center_=src.center_;
	normal_=src.normal_;
	normal_assigned_=src.normal_assigned_;
	center_assigned_=src.center_assigned_;
}

/// @brief destructor for GeneralizedKIC mover.
///
FlipChiralityMover::~FlipChiralityMover()= default;

/// @brief Clone operator to create a pointer to a fresh GeneralizedKIC object that copies this one.
///
protocols::moves::MoverOP FlipChiralityMover::clone() const {
	return protocols::moves::MoverOP( new FlipChiralityMover( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh GeneralizedKIC object that does NOT copy this one.
///
protocols::moves::MoverOP FlipChiralityMover::fresh_instance() const {
	return protocols::moves::MoverOP( new FlipChiralityMover );
}


////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void FlipChiralityMover::apply( core::pose::Pose & pose )
{
	if ( pose.size()==0 ) {
		return;
	} else {
		core::select::residue_selector::ResidueSubset subset( pose.size(), true );
		if ( selector_ ) {
			subset = selector_->apply( pose );
		}

		get_normal();
		get_center(subset, pose);

		for ( core::Size ir=1, irmax=pose.size(); ir<=irmax; ++ir ) {

			if ( subset[ir] ) {

				for ( core::Size ia=1, iamax=pose.residue(ir).natoms(); ia<=iamax; ++ia ) {
					//updating the pose coordination if it is selected as a subset
					pose.set_xyz( core::id::AtomID(ia,ir), FlipChiralityMover::calculate_reflect(pose,ia,ir));
				}

				core::conformation::ResidueOP new_res=(pose.residue(ir)).clone_flipping_chirality();
				pose.replace_residue ( ir,*new_res,false);

			}
		}
		return;
	}
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("FlipChiralityMover").
// XRW TEMP std::string FlipChiralityMover::get_name() const{
// XRW TEMP  return "FlipChiralityMover";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/scripting scheme)

void
FlipChiralityMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
) {

	std::string const name( tag->getOption<std::string>( "name" ));
	bool const user_specified_centerx( tag->hasOption("centerx"));
	bool const user_specified_centery( tag->hasOption("centery"));
	bool const user_specified_centerz( tag->hasOption("centerz"));
	if ( user_specified_centerx || user_specified_centery || user_specified_centerz ) { //If the user specified ANY...
		/// ...then he/she needs to specify ALL:
		runtime_assert_string_msg( user_specified_centerx && user_specified_centery && user_specified_centerz, "Error in parsing tag for FlipChiralityMover: if the center point is specified, then ALL of its components must be specified." );
		set_center( tag->getOption<core::Real>("centerx"),  tag->getOption<core::Real>("centery"),  tag->getOption<core::Real>("centerz")  );
	}
	bool const user_specified_normalx( tag->hasOption("normalx"));
	bool const user_specified_normaly( tag->hasOption("normaly"));
	bool const user_specified_normalz( tag->hasOption("normalz"));
	if ( user_specified_normalx || user_specified_normaly || user_specified_normalz ) { //If the user specified ANY...
		/// ...then he/she needs to specify ALL:
		runtime_assert_string_msg( user_specified_normalx && user_specified_normaly && user_specified_normalz, "Error in parsing tag for FlipChiralityMover: if the normal vector is specified, then ALL of its components must be specified." );
		set_normal( tag->getOption<core::Real>("normalx"),  tag->getOption<core::Real>("normaly"),  tag->getOption<core::Real>("normalz")  );
	}

	if ( tag->hasOption("residue_selector") ) {
		set_selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	}
}

//Private functions:

/// @brief does the actual calculation fo the reflection
//wasn't sure if can use pose.x so used this for now
numeric::xyzVector <core::Real> const FlipChiralityMover::calculate_reflect(core::pose::Pose const & pose,core::Size ia, core::Size ir){
	numeric::xyzVector <core::Real> coord_handle=pose.xyz(core::id::AtomID(ia,ir));
	numeric::xyzVector <core::Real> new_coord;
	//matrix calculation for mirror reflection of a plane. check http://mathworld.wolfram.com/Reflection.html
	core::Real cx=center_.x();
	core::Real cy=center_.y();
	core::Real cz=center_.z();
	core::Real nx=normal_.x();
	core::Real ny=normal_.y();
	core::Real nz=normal_.z();
	core::Real initx=coord_handle.x();
	core::Real inity=coord_handle.y();
	core::Real initz=coord_handle.z();


	core::Real const plane_offset=-nx*cx-ny*cy-nz*cz;

	core::Real nx2 = nx*nx;
	core::Real ny2 = ny*ny;
	core::Real nz2 = nz*nz;
	core::Real nxyz2Sum = nx2 + ny2 + nz2;
	core::Real newx=((-nx2+ny2+nz2)*initx-(2*nx*ny*inity)-(2*nx*nz*initz)-(2*nx*plane_offset))/nxyz2Sum;
	core::Real newy=(-(2*nx*ny*initx)+((nx2-ny2+nz2)*inity)-(2*ny*nz*initz)-(2*ny*plane_offset))/nxyz2Sum;
	core::Real newz=(-(2*nx*nz*initx)-(2*ny*nz*inity)+((nx2+ny2-nz2)*initz)-(2*nz*plane_offset))/nxyz2Sum;

	//TR << "plane equation is " << nx << "x + " << ny << "y + " << nz << "z + " << plane_offset << " and sqr normal is " << nxyz2Sum << std::endl;

	//TR << "these are old values " << initx << " " << inity << " " << initz << " and these are new values " << newx << " " << newy << " " <<  newz << std::endl;

	new_coord.x(newx);
	new_coord.y(newy);
	new_coord.z(newz);

	return new_coord;
}

void FlipChiralityMover::set_center( core::Real centerx, core::Real centery,core::Real centerz ){
	center_.x(centerx);
	center_.y(centery);
	center_.z(centerz);
	center_assigned_=true;

	return;
}

void FlipChiralityMover::set_normal( core::Real normalx, core::Real normaly,core::Real normalz ){
	normal_.x(normalx);
	normal_.y(normaly);
	normal_.z(normalz);
	normal_assigned_=true;

	return;
}
//Public functions definition
/// @brief sets the normal vector either from tags or calculates it.
numeric::xyzVector <core::Real> const & FlipChiralityMover::get_normal(){

	using namespace core;
	if ( !(normal_assigned_) ) {
		normal_.x(0);
		normal_.y(0);
		normal_.z(1);
	}
	return normal_;
}

/// @brief sets the normal vector either from tags or calculates it.
numeric::xyzVector <core::Real> const & FlipChiralityMover::get_center(core::select::residue_selector::ResidueSubset subset, core::pose::Pose const & pose){

	using namespace core;
	if ( !(center_assigned_) ) {
		center_=FlipChiralityMover::center_mass(subset, pose);
	}
	return center_;
}

numeric::xyzVector< core::Real >
FlipChiralityMover::center_mass(core::select::residue_selector::ResidueSubset subset, core::pose::Pose const & pose ){
	numeric::xyzVector< core::Real > massSum( 0.0 );
	Size const & nres = pose.size();
	Size nAtms=0;
	for ( Size i=1; i<= nres; ++i ) {
		if ( subset[i] ) {
			core::conformation::Residue const & rsd(pose.residue(i) );
			if ( rsd.aa() == core::chemical::aa_vrt ) continue;
			for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
				core::conformation::Atom const & atom( rsd.atom(j) );
				massSum += atom.xyz();
				nAtms++;
			}
		}
	}
	massSum /= nAtms;
	TR << "this is the center of mass " << massSum.x() << " " << massSum.y() << " " << massSum.z() << std::endl;
	return massSum;
}
/// @brief sets the selector
void FlipChiralityMover::set_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	if ( selector_in ) {
		selector_ = selector_in;
	} else {
		utility_exit_with_message("Error in protocols::cyclic_peptide::FlipChiralityMover::set_selector(): Null pointer passed to function!");
	}
	return;
}

std::string FlipChiralityMover::get_name() const {
	return mover_name();
}

std::string FlipChiralityMover::mover_name() {
	return "FlipChiralityMover";
}

void FlipChiralityMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "centerx", xsct_real, "X coordinate of center point" )
		+ XMLSchemaAttribute( "centery", xsct_real, "Y coordinate of center point" )
		+ XMLSchemaAttribute( "centerz", xsct_real, "Z coordinate of center point" )
		+ XMLSchemaAttribute( "normalx", xsct_real, "X coordinate of normal vector" )
		+ XMLSchemaAttribute( "normaly", xsct_real, "Y coordinate of normal vector" )
		+ XMLSchemaAttribute( "normalz", xsct_real, "Z coordinate of normal vector" );
	//get attributes for parse residue selector
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "Residue selector specifying residues for which chirality should be flipped" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Flips the chirality of the specified residues", attlist );
}

std::string FlipChiralityMoverCreator::keyname() const {
	return FlipChiralityMover::mover_name();
}

protocols::moves::MoverOP
FlipChiralityMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FlipChiralityMover );
}

void FlipChiralityMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FlipChiralityMover::provide_xml_schema( xsd );
}

} // moves
} // protocols
