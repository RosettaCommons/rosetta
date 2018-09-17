// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Lei Shi <shilei@u.washington.edu>
/// @date 3/19/2013

// Unit headers
#include <protocols/protein_interface_design/movers/DockWithHotspotMover.hh>
#include <protocols/protein_interface_design/movers/DockWithHotspotMoverCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <iostream>
#include <iomanip>
#include <ObjexxFCL/string.functions.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Unit Headers

// C++ headers

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace core::chemical;
using namespace std;

using core::pose::Pose;
using core::conformation::Residue;

static basic::Tracer TR( "protocols.protein_interface_design.movers.DockWithHotspotMover" );

// XRW TEMP std::string
// XRW TEMP DockWithHotspotMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DockWithHotspotMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DockWithHotspotMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DockWithHotspotMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DockWithHotspotMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "DockWithHotspotMover";
// XRW TEMP }

/// @brief default ctor
DockWithHotspotMover::DockWithHotspotMover() :
	protocols::moves::Mover( DockWithHotspotMover::mover_name() )
{
}

DockWithHotspotMover::~DockWithHotspotMover() = default;

void DockWithHotspotMover::apply( Pose & pose ) {
	//convert chain B to be alanine
	Pose archive_pose=pose;
	protocols::protein_interface_design::movers::BuildAlaPose toAla( true,true,100);
	toAla.apply( pose );

	//scoring function
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::ScoreFunctionFactory::create_score_function("interchain_cen") );
	core::scoring::ScoreFunctionOP scorefxn_emp( core::scoring::ScoreFunctionFactory::create_score_function("empty") );
	core::scoring::ScoreFunctionOP scorefxn_cen( core::scoring::ScoreFunctionFactory::create_score_function("interchain_cen") );

	//load hotspot constraints
	//core::Size fixed_res(1);  // unused ~Labonte
	core::Size const chain_to_redesign = 2;
	//if ( chain_to_redesign == 1 ) fixed_res = pose.size();  // unused ~Labonte
	// unused variable below ~Labonte
	//core::id::AtomID fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );
	core::Real const worst_allowed_stub_bonus(-1.);
	bool const apply_self_energies(false);
	core::Real const bump_cutoff(10.);
	bool const apply_ambiguous_constraints(true);

	for ( Size i=1; i <= hotspot_filenames_.size(); i++ ) {
		protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_setOP( new protocols::hotspot_hashing::HotspotStubSet );
		hotspot_stub_setOP->read_data( hotspot_filenames_[i] );
		hotspot_stub_setOP->add_hotspot_constraints_to_wholepose( pose, chain_to_redesign, hotspot_stub_setOP,
			hotspot_distcb_weight_[i], worst_allowed_stub_bonus, apply_self_energies, bump_cutoff,
			apply_ambiguous_constraints );
	}

	scorefxn->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight_ );
	scorefxn_emp->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight_ );

	//switch to centroid
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
	core::Size const rb_move_jump = 1;
	protocols::docking::DockingLowResOP docking_lowres_mover( new protocols::docking::DockingLowRes( scorefxn, rb_move_jump ) );
	docking_lowres_mover->apply(pose);
	(*scorefxn_cen)(pose);
	utility::vector1< core::Size > movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	core::Real CenScore=protocols::docking::calc_interaction_energy(pose,scorefxn_cen,movable_jumps_);

	protocols::forge::methods::restore_residues( archive_pose, pose );
	core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );

	(*scorefxn_emp)(pose);
	core::Real CstScore= pose.energies().total_energies().dot( scorefxn_emp->weights() );

	protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );

	if ( CenScore > centroidscore_filter_ || CstScore > hotspotcst_filter_ ) {
		set_last_move_status( protocols::moves::FAIL_RETRY );
		TR<< protocols::jd2::current_input_tag() <<
			" Did not pass filter CentoridScore ( " << centroidscore_filter_ <<  " ): " << CenScore <<
			"HotspotCstScore ( " << hotspotcst_filter_ << " ): " << CstScore << std::endl;
		return;
	}
	TR << protocols::jd2::current_input_tag() <<
		" Succeed Centorid Interaction: "<< CenScore << " HotspotCstScore: " << CstScore <<std::endl;
	set_last_move_status( protocols::moves::MS_SUCCESS);
}

// XRW TEMP std::string
// XRW TEMP DockWithHotspotMover::get_name() const {
// XRW TEMP  return DockWithHotspotMover::mover_name();
// XRW TEMP }


/**
* @brief Reinitialize this protocols::moves::Mover with parameters from the specified tags.
* @details Parameters recognized:
*  - target_pdb_num or target_res_num. A single target residue to form disulfides to
*  - target_pdb_nums or target_res_nums. A list of possible target residues
*/
void DockWithHotspotMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & /*pose*/)
{
	hotspot_score_weight_=tag->getOption<core::Real>( "hotspot_score_weight", 10 );
	centroidscore_filter_=tag->getOption<core::Real>( "centroidscore_filter", 0 );

	TR << "Setup DockWithHotspotMover with backbone_stub_linear_constraint with weight "<< hotspot_score_weight_ <<
		" from:" << std::endl;
	// Set target to the residue specified by "target_pdb_num" or "target_res_num"
	utility::vector1< TagCOP > const branch_tags( tag->getTags() );
	for ( TagCOP curr_tag : branch_tags ) {
		if ( curr_tag->getName() != "HotspotFiles" ) {
			TR.Error << "No 'HotspotFiles' specified." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "");
		} else {
			utility::vector1< TagCOP > const branch_tags2( curr_tag->getTags() );
			for ( TagCOP curr_tag2 : branch_tags2 ) {
				std::string const file_name( curr_tag2->getOption< std::string >( "file_name" ) );
				hotspot_filenames_.push_back(file_name);
				auto cb_force( curr_tag2->getOption< core::Real >( "cb_force", 1.0 ) );
				hotspot_distcb_weight_.push_back(cb_force);
				TR << "file: " << file_name << " with cb_force set to " << cb_force << std::endl;
			}
		}
	}
	hotspotcst_filter_=tag->getOption<core::Real>( "hotspotcst_filter", hotspot_filenames_.size()*40.0 );
	TR << "centroid score filter: " << centroidscore_filter_ << std::endl;
	TR << "hotspotcst_filter: " << hotspotcst_filter_ << std::endl;
}

std::string DockWithHotspotMover::get_name() const {
	return mover_name();
}

std::string DockWithHotspotMover::mover_name() {
	return "DockWithHotspotMover";
}

std::string dock_with_hotspot_subtag_naming_func( std::string const & foo ) {
	return "dock_with_hotspot_subtag_" + foo + "_type";
}

std::string dock_with_hotspot_subsubtag_naming_func( std::string const & foo ) {
	return "dock_with_hotspot_subsubtag_" + foo + "_type";
}

void DockWithHotspotMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "hotspot_score_weight", xsct_real, "Weight for the hotpot scoreterm", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "centroidscore_filter", xsct_real, "Threshold for the centroid score filter", "0" )
		+ XMLSchemaAttribute::required_attribute( "hotspotcst_filter", xsct_real, "XRW TO DO");

	AttributeList subsubattributes;
	subsubattributes + XMLSchemaAttribute::required_attribute( "file_name", xs_string, "Filename for this hotspot file" )
		+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_real, "Weight on the CB atoms for these hotspot constraints", "1.0" );

	// AMW: we don't know what these subelements are traditionally called
	utility::tag::XMLSchemaComplexTypeGenerator hotspot_file_ct_gen;
	hotspot_file_ct_gen.complex_type_naming_func( & dock_with_hotspot_subsubtag_naming_func )
		.element_name( "HotspotFile" )
		.description( "Wrapper tag for individual hotspot files" )
		.add_attributes( subsubattributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "HotspotFile", & dock_with_hotspot_subsubtag_naming_func/*, 0*/ );

	utility::tag::XMLSchemaComplexTypeGenerator hotspot_files_ct_gen;
	hotspot_files_ct_gen.complex_type_naming_func( & dock_with_hotspot_subtag_naming_func )
		.element_name( "HotspotFiles" )
		.description( "Wrapper tag for individual hotspot files" )
		//.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_already_defined_subelement( "HotspotFiles", & dock_with_hotspot_subtag_naming_func/*, 0*/ );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string DockWithHotspotMoverCreator::keyname() const {
	return DockWithHotspotMover::mover_name();
}

protocols::moves::MoverOP
DockWithHotspotMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DockWithHotspotMover );
}

void DockWithHotspotMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DockWithHotspotMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols

