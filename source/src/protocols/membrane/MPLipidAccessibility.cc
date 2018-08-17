// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/MPLipidAccessibility.cc
/// @brief Mover that computes which residues are lipid accessible and puts that information into the B-factors: 50 is lipid accessible, 0 is lipid INaccessible
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author Rebecca Alford (rfalford12@gmail.com) for minor refactoring

#ifndef INCLUDED_protocols_membrane_MPLipidAccessibility_cc
#define INCLUDED_protocols_membrane_MPLipidAccessibility_cc

// Unit headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MPLipidAccessibility.hh>
#include <protocols/membrane/MPLipidAccessibilityCreator.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// Package header
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>
#include <core/membrane/hull.hh>

// Basic/Utility headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <protocols/membrane/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.membrane.MPLipidAccessibility" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
MPLipidAccessibility::MPLipidAccessibility() : protocols::moves::Mover() {

	set_defaults();
	register_options();
	init_from_cmd();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
MPLipidAccessibility::~MPLipidAccessibility()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
MPLipidAccessibility::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
MPLipidAccessibility::parse_my_tag(
	utility::tag::TagCOP tag
) {
	if ( tag->hasOption( "angle_cutoff" ) ) {
		angle_cutoff_ = tag->getOption< core::Real >("angle_cutoff", 65.0);
	}

	if ( tag->hasOption( "slice_width" ) ) {
		slice_width_ = tag->getOption< core::Real >("slice_width", 10.0);
	}

	if ( tag->hasOption( "shell_radius" ) ) {
		shell_radius_ = tag->getOption< core::Real >("shell_radius", 6.0);
	}

	if ( tag->hasOption( "dist_cutoff" ) ) {
		dist_cutoff_ = tag->getOption< core::Real >("dist_cutoff", 10.0);
	}

	if ( tag->hasOption( "tm_alpha" ) ) {
		tm_alpha_ = tag->getOption< bool >("tm_alpha", true);
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MPLipidAccessibility::fresh_instance() const
{
	return protocols::moves::MoverOP( new MPLipidAccessibility );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MPLipidAccessibility::clone() const
{
	return protocols::moves::MoverOP( new MPLipidAccessibility( *this ) );
}

/////////////////////
/// Get Methods   ///
/////////////////////

/// @brief get angle cutoff
core::Real MPLipidAccessibility::get_angle_cutoff() const {
	return angle_cutoff_;
}

/// @brief get slice width
core::Real MPLipidAccessibility::get_slice_width() const {
	return slice_width_;
}

/// @brief get shell radius
core::Real MPLipidAccessibility::get_shell_radius() const {
	return shell_radius_;
}

/// @brief get dist cutoff
core::Real MPLipidAccessibility::get_dist_cutoff() const {
	return dist_cutoff_;
}

/// @brief get tm_alpha
bool MPLipidAccessibility::get_tm_alpha() const {
	return tm_alpha_;
}


/// @brief Return a vector1 of vector1 matching the pose and atom counts with
/// the structure-based per-atom lipid accessibility
utility::vector1< utility::vector1< core::Real > >
MPLipidAccessibility::get_per_atom_lipid_accessibility() const {
	return per_atom_lipid_accessibility_;
}

////////////////////////////////////////////////////////////////////////////////
/// MOVER METHODS ///
/////////////////////
void MPLipidAccessibility::apply( core::pose::Pose & pose ){

	using namespace core::id;
	using namespace core::pose;

	// finalize setup
	finalize_setup( pose );

	////////////////// FILL ACCESSIBILITY VECTOR ///////////////////////////
	per_atom_lipid_accessibility_.resize( nres_protein( pose ) );
	for ( core::Size ii = 1; ii <= nres_protein( pose ); ++ii ) {
		per_atom_lipid_accessibility_[ ii ].resize( pose.residue( ii ).natoms() );
	}

	////////////////// SET B-FACTOR TO ZERO FOR ALL ATOMS //////////////////
	for ( core::Size r = 1; r <= nres_protein( pose ); ++r ) {
		for ( core::Size a = 1; a <= pose.residue( r ).natoms(); ++a ) {
			per_atom_lipid_accessibility_[r][a] = 0.0;
		}
	}

	//////////////////////// COMPUTE HULL///// /////////////////////////////

	// compute concave hull with membrane boundary
	// current solution until implicit lipid PR is merged
	core::Real minz(0);
	core::Real maxz(0);
	minz = -pose.membrane_info()->membrane_thickness();
	maxz = pose.membrane_info()->membrane_thickness();

	// compute the outer shell
	// default parameters are a little different for smaller proteins, like GPCRs
	core::id::AtomID_Map< bool > shell = core::membrane::concave_shell( pose, minz, maxz, slice_width_, shell_radius_, dist_cutoff_ );

	//////////////////////// GO THROUGH SLICES /////////////////////////////

	// get lipid-exposed residues from the AtomID_map
	utility::vector1< core::Size > boundary;
	for ( core::Size r = 1; r <= nres_protein( pose ); ++r ) {

		// create AtomID for CA atom
		AtomID aid = AtomID( 2, r );

		// if atomID has a true in the map, add to boundary vector
		// this relies on the fact that CA is set for true if any of the atoms
		// is is lipid exposed!
		if ( shell.has( aid ) == true && shell.get( aid ) == true ) {
			boundary.push_back( r );
		}
	}

	// go through slices
	for ( core::Size s = 1; s < slice_zmin_.size(); ++s ) {

		// go through residues
		for ( core::Size r = 1; r <= resi_[ s ].size(); ++r ) {

			bool lipid_exposed( false );

			// if residues are on boundary, then lipid exposed
			if ( boundary.has_value( resi_[ s ][ r ] ) ) {
				lipid_exposed = true;
			}

			// if angle between COM-CA-CB is smaller than cutoff, then not lipid-exposed
			// compute angle between COM-CA-CB
			// if angle < angle_cutoff, then sidechain face COM
			// if angle >= angle_cutoff, then sidechain faces away from COM
			core::Real angle = numeric::angle_degrees( slice_com_[ s ], ca_coord_[ s ][ r ], cb_coord_[ s ][ r ] );

			// go through atoms in residue and set B-factor
			for ( core::Size a = 1; a <= pose.residue( resi_[ s ][ r ] ).natoms(); ++a ) {

				// is this a 2TM protein? if yes, all residues are lipid exposed
				if ( pose.membrane_info()->spanning_topology()->nspans() <= 2 ) {
					per_atom_lipid_accessibility_[ resi_[s][r] ][ a ] = 50.0;
				}

				// if lipid exposed and COM-CA-CB larger than cutoff, i.e. facing outwards
				if ( ( lipid_exposed == true && tm_alpha_ == true )
						|| ( lipid_exposed == true && tm_alpha_ == false && angle > angle_cutoff_ ) ) {
					per_atom_lipid_accessibility_[ resi_[s][r] ][ a ] = 50.0;
				}
			} // atoms
		} // residues
	} // slices

	// Transfering information from per_atom lipid accessibility to b_factor data at the end
	for ( core::Size ii = 1; ii <= nres_protein( pose ); ++ii ) {
		for ( core::Size jj = 1; jj <= pose.residue( ii ).natoms(); ++jj ) {
			pose.pdb_info()->bfactor( ii, jj, per_atom_lipid_accessibility_[ii][jj] );
		}
	}

	TR << "tm alpha? " << tm_alpha_ << " (helical proteins with <= 7 TMs will show up as 'not helical')" << std::endl;

}// apply

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

/// @brief set defaults
void MPLipidAccessibility::set_defaults() {

	angle_cutoff_ = 65.0;

	// slices from hull data
	slice_width_ = 10.0;
	shell_radius_ = 6.0;
	dist_cutoff_ = 10.0;
	tm_alpha_ = true;

}// set defaults

//////////////////////////////////////////
/// @brief Register Options from Command Line
void MPLipidAccessibility::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::lipid_acc::angle_cutoff );
	option.add_relevant( OptionKeys::mp::lipid_acc::slice_width );
	option.add_relevant( OptionKeys::mp::lipid_acc::shell_radius );
	option.add_relevant( OptionKeys::mp::lipid_acc::dist_cutoff );
	option.add_relevant( OptionKeys::mp::lipid_acc::tm_alpha );

}// register options

//////////////////////////////////////////
/// @brief Initialize from commandline
void MPLipidAccessibility::init_from_cmd() {

	using namespace basic::options;

	if ( option[ OptionKeys::mp::lipid_acc::angle_cutoff ].user() ) {
		angle_cutoff_ = option[ OptionKeys::mp::lipid_acc::angle_cutoff ]();
	}

	if ( option[ OptionKeys::mp::lipid_acc::slice_width ].user() ) {
		slice_width_ = option[ OptionKeys::mp::lipid_acc::slice_width ]();
	}

	if ( option[ OptionKeys::mp::lipid_acc::shell_radius ].user() ) {
		shell_radius_ = option[ OptionKeys::mp::lipid_acc::shell_radius ]();
	}

	if ( option[ OptionKeys::mp::lipid_acc::dist_cutoff ].user() ) {
		dist_cutoff_ = option[ OptionKeys::mp::lipid_acc::dist_cutoff ]();
	}

	if ( option[ OptionKeys::mp::lipid_acc::tm_alpha ].user() ) {
		tm_alpha_ = option[ OptionKeys::mp::lipid_acc::tm_alpha ]();
	}

}// init from cmd

//////////////////////////////////////////
/// @brief finalize setup
void MPLipidAccessibility::finalize_setup( core::pose::Pose & pose ){

	using namespace basic::options;
	using namespace protocols::membrane;

	////////// SWITCH TO FULL-ATOM, REQUIRED FOR VECTOR DEFINITIONS ////////
	using namespace protocols::simple_moves;
	SwitchResidueTypeSetMoverOP full_atom( new SwitchResidueTypeSetMover( "fa_standard" ) );
	full_atom->apply( pose );

	// add membrane
	if ( ! pose.conformation().is_membrane() ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// check whether protein is in membrane
	if ( ! protein_in_membrane( pose ) ) {
		TR.Warning << "YOUR PROTEIN DOES NOT SPAN THE MEMBRANE! EITHER YOU KNOW WHAT YOU ARE DOING OR YOUR PROTEIN IS NOT TRANSFORMED INTO MEMBRANE COORDINATES!!!" << std::endl;
		//  utility_exit_with_message( "Your protein is not transformed into membrane coordinates! Cannot compute lipid accessibility. Quitting." );
	}

	// fill up secondary structure information in the pose
	utility::vector1< char > secstruct( get_secstruct( pose ) );

	// figure out protein type: helical bundle or beta-barrel
	core::Size nbeta = 0;
	core::Size nmem = 0;

	// count the number of (helical) residues in the membrane
	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {
		if ( pose.membrane_info()->in_membrane( pose.conformation(), i ) ) {
			nmem++;
			if ( secstruct[ i ] == 'E' ) {
				nbeta++;
			}
		}
	}

	// helical or barrel?
	core::Real beta = static_cast< core::Real >( nbeta ) / static_cast< core::Real >( nmem );
	TR << "nbeta: " << nbeta << " nmem: " << nmem << " beta: " << beta << std::endl;
	if ( nmem > 0 && ! option[ OptionKeys::mp::lipid_acc::tm_alpha ].user() && beta >= 0.5 ) {
		tm_alpha_ = false;
	}

	// fill up the data in the slices
	fill_up_slices( pose );

	// set different default parameters for smaller proteins, like GPCRs
	if ( pose.membrane_info()->spanning_topology()->nspans() <= 7 ) {

		// this flag only specifies whether we want to use an angle cutoff:
		// beta-barrels have one but helical ones don't
		// this breaks down for smaller proteins though
		tm_alpha_ = false;
		if ( ! option[ OptionKeys::mp::lipid_acc::angle_cutoff ].user() ) {
			angle_cutoff_ = 45.0;
		}
	}

}// finalize setup

//////////////////////////////////////////
/// @brief check whether protein is transformed into membrane
bool MPLipidAccessibility::protein_in_membrane( core::pose::Pose & pose ){

	bool in_membrane( true );

	// get minimum and maximum z-coordinate
	// check whether structure is transformed into membrane, get CA position
	core::Real min_z( 10000 );
	core::Real max_z( -10000 );

	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

		// skip for membrane residue
		if ( pose.residue( i ).name3() == "MEM" ) {
			continue;
		}

		if ( pose.residue( i ).xyz( "CA" ).z() < min_z ) {
			min_z = pose.residue( i ).xyz( "CA" ).z();
		} else if ( pose.residue( i ).xyz( "CA" ).z() > max_z ) {
			max_z = pose.residue( i ).xyz( "CA" ).z();
		}
	}

	// crude check: if no CA atoms with positive AND negative z-coordinates,
	// then the protein doesn't span the membrane
	if ( !( min_z < 0 && max_z > 0 ) ) {
		in_membrane = false;
	}
	return in_membrane;

}// protein in membrane?

//////////////////////////////////////////
/// @brief fill up slice arrays with protein data
void MPLipidAccessibility::fill_up_slices( core::pose::Pose & pose ) {

	using namespace core::pose;

	// check for membrane protein
	if ( ! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Something went wrong, the AddMembraneMover still didn't make it a membrane protein. Quitting..." );
	}

	// Currently: use the default thickness in membrane info
	core::Real iter(0);
	core::Real thickness(0);
	iter = -pose.membrane_info()->membrane_thickness();
	thickness = pose.membrane_info()->membrane_thickness();

	// go through protein in the membrane and get slices and arrays for them
	while ( iter <= thickness ) {

		slice_zmin_.push_back( iter );

		// temp utility vectors for each slice
		utility::vector1< core::Size > slice_res;
		utility::vector1< core::Vector > slice_ca, slice_cb;

		// go through protein residues
		for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

			// skip membrane residue
			if ( pose.residue( i ).name3() == "MEM" ) {
				continue;
			}

			// get CA and CB coordinates
			core::Vector ca = pose.residue( i ).xyz( "CA" );
			core::Vector cb;

			// use 2HA instead of CB for glycine
			if ( pose.residue( i ).name3() == "GLY" ) {
				cb = pose.residue( i ).xyz( "2HA" );
			} else {
				cb = pose.residue( i ).xyz( "CB" );
			}

			// if CA-coordinate is within slice, add resi, CA and CB to arrays
			if ( ( ca.z() >= iter && ca.z() < iter + slice_width_ ) or
					( cb.z() >= iter && cb.z() < iter + slice_width_ ) ) {

				slice_res.push_back( i );
				slice_ca.push_back( ca );
				slice_cb.push_back( cb );

			}
		}

		// push back into slices
		resi_.push_back( slice_res );
		ca_coord_.push_back( slice_ca );
		cb_coord_.push_back( slice_cb );

		// go to next slice
		iter += slice_width_;

	} // get arrays of slices

	// compute slice COMs
	compute_slice_com();

}// fill up slice data

//////////////////////////////////////////
/// @brief compute slice COM
void MPLipidAccessibility::compute_slice_com(){

	// go through slices and compute COMs
	for ( core::Size s = 1; s < slice_zmin_.size(); ++s ) {

		core::Vector com( 0, 0, 0 );

		// go through residues, compute COM per slice
		for ( core::Size r = 1; r <= resi_[ s ].size(); ++r ) {
			com += ca_coord_[ s ][ r ];
		}
		com /= resi_[ s ].size();

		TR.Debug << "slice " << s << " com " << com.to_string() << std::endl;

		slice_com_.push_back( com );
	}

}// compute slice COM

std::string MPLipidAccessibility::get_name() const {
	return mover_name();
}

std::string MPLipidAccessibility::mover_name() {
	return "MPLipidAccessibility";
}

void MPLipidAccessibility::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "angle_cutoff", xsct_real, "Angle cutoff")
		+ XMLSchemaAttribute    ( "slice_width", xsct_real, "Slice width")
		+ XMLSchemaAttribute    ( "shell_radius", xsct_real, "Shell radius")
		+ XMLSchemaAttribute    ( "dist_cutoff", xsct_real, "Distance cutoff")
		+ XMLSchemaAttribute    ( "tm_alpha", xsct_rosetta_bool, "tm_alpha");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Mover that computes which residues are lipid accessible and puts that information "
		"into the B-factors: 50 is lipid accessible, 0 is lipid INaccessible",
		attlist );
}

std::string MPLipidAccessibilityCreator::keyname() const {
	return MPLipidAccessibility::mover_name();
}

protocols::moves::MoverOP
MPLipidAccessibilityCreator::create_mover() const {
	return protocols::moves::MoverOP( new MPLipidAccessibility );
}

void MPLipidAccessibilityCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MPLipidAccessibility::provide_xml_schema( xsd );
}


} //protocols
} //membrane

#endif // INCLUDED_protocols_membrane_MPLipidAccessibility_cc
