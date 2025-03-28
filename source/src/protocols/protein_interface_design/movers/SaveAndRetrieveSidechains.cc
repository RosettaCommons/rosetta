// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.cc
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechainsCreator.hh>

// Project headers
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ref_pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <protocols/calc_taskop_movers/DesignRepackMover.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/datacache/DataMapObj.hh> // AUTO IWYU For DataMapObj

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.SaveAndRetrieveSidechains" );




SaveAndRetrieveSidechains::SaveAndRetrieveSidechains() :
	calc_taskop_movers::DesignRepackMover( SaveAndRetrieveSidechains::mover_name() )
{
	allsc_ = false; // default
	jumpid_ = 1; //default
	ensure_variant_matching_ = false; //default
	multi_use_ = false;
	two_step_ = false;
	first_apply_ = utility::pointer::make_shared< basic::datacache::DataMapObj< bool > >();
	first_apply_->obj = true;
	init_pose_ = utility::pointer::make_shared< core::pose::Pose >();
}

SaveAndRetrieveSidechains::SaveAndRetrieveSidechains(
	core::pose::Pose const & pose,
	bool const allsc /*=false*/,
	bool const ensure_variant_matching /*=false*/,
	core::Size const jumpid /*=1*/
) :
	calc_taskop_movers::DesignRepackMover( SaveAndRetrieveSidechains::mover_name() ),
	allsc_( allsc ),
	ensure_variant_matching_(ensure_variant_matching),
	jumpid_( jumpid )
{
	init_pose_ = utility::pointer::make_shared< core::pose::Pose >( pose );
	multi_use_ = false;
	two_step_ = false;
	first_apply_ = utility::pointer::make_shared< basic::datacache::DataMapObj< bool > >();
	first_apply_->obj = true;
}

SaveAndRetrieveSidechains::~SaveAndRetrieveSidechains() = default;

void
SaveAndRetrieveSidechains::apply( Pose & pose )
{
	if ( two_step() && first_apply_->obj ) {
		TR<<"Saving sidechains."<<std::endl;
		init_pose_ = pose.clone();
		first_apply_->obj = false;
		return;
	}
	if ( multi_use() ) {
		first_apply_->obj = true;
	}
	TR << "Retrieving sidechains..."<<std::endl;
	core::Size nres = pose.size();
	if ( nres != init_pose_->size() && core::pose::symmetry::is_symmetric(pose) ) {
		auto & symm_conf (
			dynamic_cast<conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		nres = symm_conf.Symmetry_Info()->num_independent_residues();
	}
	runtime_assert( nres == init_pose_->size() );
	kinematics::Jump new_jump;
	core::Size const rb_jump( jumpid_ );
	if ( jumpid_ > 0 ) {
		new_jump = pose.jump( rb_jump );
	}

	for ( core::Size res=1; res<=nres; ++res ) {
		if ( allsc_ ) { // replace all sidechains
			pose.replace_residue( res, init_pose_->residue( res ), true/*orient_backbone*/ );
			continue;
		} else {
			if ( pose.residue( res ).name3() == "ALA" ) { // only replace Ala positions
				pose.replace_residue( res, init_pose_->residue( res ), true/*orient_backbone*/ );
			}
		}
	}
	if ( ensure_variant_matching_ ) {
		//make sure variants match, if not put back the initial variants
		using namespace core;
		for ( core::uint i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( ! variants_match( pose.residue_type( i ), init_pose_->residue_type( i ) ) ) {
				utility::vector1< std::string > const new_var_types(
					pose.residue_type( i ).properties().get_list_of_variants() );
				utility::vector1< std::string > const old_var_types(
					init_pose_->residue_type( i ).properties().get_list_of_variants() );
				for ( auto const & new_var_type : new_var_types ) {
					if ( ! (init_pose_->residue_type( i ).has_variant_type( new_var_type ) ) ) {
						core::pose::remove_variant_type_from_pose_residue( pose,
							core::chemical::ResidueProperties::get_variant_from_string( new_var_type ), i );
					}
				}

				for ( auto const & old_var_type : old_var_types ) {
					if ( !pose.residue_type( i ).has_variant_type( old_var_type ) ) {
						core::pose::add_variant_type_to_pose_residue( pose,
							core::chemical::ResidueProperties::get_variant_from_string( old_var_type ), i );
					}
				}
			} //if variants don't match
		}
	}
	if ( jumpid_ > 0 ) {
		pose.set_jump( rb_jump, new_jump );
	}
	TR.flush();
}


void
SaveAndRetrieveSidechains::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
{
	first_apply_->obj = true;
	allsc_ = tag->getOption<bool>( "allsc", false );
	multi_use( tag->getOption< bool >( "multi_use", false ) );
	two_step( tag->getOption< bool >( "two_step", false ) );
	if ( !two_step() ) {
		init_pose_ = protocols::rosetta_scripts::legacy_saved_pose_or_input( tag, data, mover_name(), /*use_native*/ false );
	}
	jumpid_ = tag->getOption<core::Size>( "jumpid", 1 );
}

protocols::moves::MoverOP
SaveAndRetrieveSidechains::clone() const {
	return( utility::pointer::make_shared< SaveAndRetrieveSidechains >( *this ));
}

std::string SaveAndRetrieveSidechains::get_name() const {
	return mover_name();
}

std::string SaveAndRetrieveSidechains::mover_name() {
	return "SaveAndRetrieveSidechains";
}

void SaveAndRetrieveSidechains::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "allsc", xsct_rosetta_bool, "Save and retrieve all sidechains", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "multi_use", xsct_rosetta_bool, "Set up so that we can use this multiple times", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "two_step", xsct_rosetta_bool, "Save and retrieve in two steps", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "jumpid", xsct_non_negative_integer, "Jump ID to keep track of", "1" );

	core::pose::attributes_for_saved_reference_pose_w_description( attlist, "If two_step is false, what structure to get the sidechains from." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SaveAndRetrieveSidechainsCreator::keyname() const {
	return SaveAndRetrieveSidechains::mover_name();
}

protocols::moves::MoverOP
SaveAndRetrieveSidechainsCreator::create_mover() const {
	return utility::pointer::make_shared< SaveAndRetrieveSidechains >();
}

void SaveAndRetrieveSidechainsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SaveAndRetrieveSidechains::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
