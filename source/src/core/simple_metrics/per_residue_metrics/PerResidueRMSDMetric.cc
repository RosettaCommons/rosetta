// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueRMSDMetric.cc
/// @brief A per-residue metric thtat will calculate the RMSD for each residue given in a residue selector to a reference pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueRMSDMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ref_pose.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/symmetry/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueRMSDMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::scoring;
using namespace core::pose;
using namespace core::select::residue_selector;
using namespace utility;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueRMSDMetric::PerResidueRMSDMetric():
	core::simple_metrics::PerResidueRealMetric()
{
	setup_name_mapping();
}

/// @brief Default constructor
PerResidueRMSDMetric::PerResidueRMSDMetric(core::pose::PoseCOP ref_pose):
	core::simple_metrics::PerResidueRealMetric(),
	ref_pose_(ref_pose)
{
	setup_name_mapping();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueRMSDMetric::~PerResidueRMSDMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueRMSDMetric::PerResidueRMSDMetric( PerResidueRMSDMetric const & ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueRMSDMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PerResidueRMSDMetric( *this ) );

}

void
PerResidueRMSDMetric::setup_name_mapping(){

	rmsd_atom_names_ = setup_rmsd_atom_names();
	name_mapping_ = get_rmsd_type_name_map();
}

std::string
PerResidueRMSDMetric::name() const {
	return name_static();
}

std::string
PerResidueRMSDMetric::name_static() {
	return "PerResidueRMSDMetric";

}
std::string
PerResidueRMSDMetric::metric() const {
	return "res_rmsd";
}

void
PerResidueRMSDMetric::set_residue_selector_reference(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ref_ = residue_selector;
}

void
PerResidueRMSDMetric::set_residue_mapping(std::map<core::Size, core::Size> const & rmsd_map ){
	rmsd_map_ = rmsd_map;
}

void
PerResidueRMSDMetric::set_comparison_pose(core::pose::PoseCOP ref_pose){
	ref_pose_ = ref_pose;
}

void
PerResidueRMSDMetric::set_corresponding_atoms_robust(bool robust){
	robust_ = robust;
}

void
PerResidueRMSDMetric::set_rmsd_type(scoring::rmsd_atoms rmsd_type){
	rmsd_type_ = rmsd_type;
}

void
PerResidueRMSDMetric::set_run_superimpose(bool super){
	superimpose_ = super;
}

void
PerResidueRMSDMetric::set_desymmetrize_residue_selector( bool desym){
	desymmetrize_res_selector_ = desym;
}

void
PerResidueRMSDMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	parse_per_residue_tag( tag, datamap );

	if ( tag->hasOption("residue_selector_ref") ) {
		set_residue_selector_reference(parse_residue_selector( tag, datamap, "residue_selector_ref" ));
	}

	set_run_superimpose(tag->getOption< bool >("super", superimpose_));

	//Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		ref_pose_ = saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption<bool>("use_native", false) && datamap.has_resource("native_pose") ) {
		ref_pose_ = saved_native_pose(datamap);
	} else {
		std::string msg = "A reference pose must be set. Please use the SavePoseMover (embed the RMSDMetric in RunSimpleMetrics ) or pass the native as in:file:native and set use_native to true.";
		utility_exit_with_message(msg);
	}

	if ( tag->hasOption("rmsd_type") ) {
		set_rmsd_type( name_mapping_[ tag->getOption<std::string>("rmsd_type")]);
	}

	set_desymmetrize_residue_selector( tag->getOption< bool >("desymmetrize_selector", true));

	set_corresponding_atoms_robust(tag->getOption< bool >("robust", robust_));
}

void
PerResidueRMSDMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector_ref",
		"Selector for the reference pose (input native or passed reference pose. ).  Residues selected must be same number of residues selected for the main selector." );


	attlist + XMLSchemaAttribute::attribute_w_default(
		"robust", xsct_rosetta_bool, "Set whether we are robust to atom mismatches for selected residues."
		"  By default we only match atoms that are corresponding. (True).", "true"
	);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false"
	);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"super", xsct_rosetta_bool, "Run a superposition on the residues in the residue_selector (or all) before RMSD calculation and the atoms slected for RMSD", "false"
	);

	utility::vector1< std::string > rmsd_type_names = get_rmsd_type_names();
	utility::tag::add_schema_restrictions_for_strings( xsd, "rmsd_types", rmsd_type_names);

	attlist + XMLSchemaAttribute("rmsd_type", "rmsd_types", "Type of calculation.  Current choices are: \n" + utility::to_string(rmsd_type_names) );
	attlist + XMLSchemaAttribute::attribute_w_default("desymmetrize_selector", xsct_rosetta_bool, "Should we desymmetrize the residue selector if we have a symmetric pose?", "true");

	std::string description = "\tThis is the RMSD for each residue between the input and the set comparison pose.\n"
		"  If native is set on the cmd-line and use_native is true, we will use that.\n"
		"  Default is to calculate all_heavy atoms - but this can be set\n."
		"\n"
		"  Make sure that reference pose and set pose are the same length ."
		"   We match all corresponding atoms for each residue to match."
		"   By default we do not fail and are robust - only matching what we can for each residue.";

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);

}

std::map< id::AtomID, id::AtomID >
PerResidueRMSDMetric::create_atom_id_map( core::pose::Pose const & pose, bool desymmetrize_res_selector) const {
	//Setup the main residue mapping we will use depending on what is set.
	std::map< core::Size, core::Size > residue_map;
	if ( rmsd_map_.empty() ) {
		residue_map = get_residue_mapping_from_selectors( get_selector(), residue_selector_ref_, pose, *ref_pose_, desymmetrize_res_selector);
	} else {
		residue_map = rmsd_map_;
	}

	std::map< core::id::AtomID, core::id::AtomID > atom_map;

	//, backbone_heavy_map, heavy_map;
	//setup_matching_CA_atoms(                     pose, *native_pose_op, CA_map );
	//setup_matching_protein_backbone_heavy_atoms( pose, *native_pose_op, backbone_heavy_map );
	//setup_matching_heavy_atoms(                  pose, *native_pose_op, heavy_map );

	//Match corresponding.  Skip atom if they do not match type.
	for ( auto res_pair : residue_map ) {
		utility::vector1< std::string > atom_names = rmsd_atom_names_.at( rmsd_type_);
		if ( ! override_atom_names_.empty() ) {
			atom_names = override_atom_names_;
		}

		//Match on atomnames
		if ( ! atom_names.empty() ) {
			for ( std::string const & atom_name : atom_names ) {
				if ( ! pose.residue_type( res_pair.first).has( atom_name ) ) continue;
				if ( ! ref_pose_->residue_type( res_pair.second).has( atom_name ) ) continue;

				core::Size pose_atom = pose.residue_type( res_pair.first).atom_index( atom_name );
				core::Size ref_atom = pose.residue_type( res_pair.second).atom_index( atom_name );

				atom_map[ core::id::AtomID( pose_atom, res_pair.first)] = core::id::AtomID( ref_atom, res_pair.second);
			}
		} else {
			core::Size total_pose_atoms = 0;
			//core::Size total_ref_atoms = 0;

			if ( ( rmsd_type_ == rmsd_sc_heavy ) || (rmsd_type_ == rmsd_all_heavy) ) {
				total_pose_atoms = pose.residue_type( res_pair.first ).nheavyatoms();
				//total_ref_atoms = ref_pose_->residue_type( res_pair.second).nheavyatoms();
			} else {
				total_pose_atoms = pose.residue_type( res_pair.first ).natoms();
				//total_ref_atoms = ref_pose_->residue_type( res_pair.second).natoms();
			}


			for ( core::Size i = 1; i <= total_pose_atoms; ++i ) {
				//Skip Virtuals
				if ( pose.residue_type( res_pair.first).is_virtual( i ) ) continue;

				//If we only sidechains - skip if not sidechains.
				if ( (rmsd_type_ == rmsd_sc ) || (rmsd_type_ == rmsd_sc_heavy) ) {
					if ( ! pose.residue_type( res_pair.first).atom_is_sidechain( i ) ) continue;
				}


				core::chemical::Atom const & pose_atom = pose.residue_type( res_pair.first).atom( i );

				//if ((ref_atom_index != i) && (robust_ == false)){
				// std::string msg = "RMSDMetric.  Robust set to false.  Residues do not match:"+utility::to_string( res_pair.first)+" atomno "+utility::to_string(i)+" "+pose_atom.name()+" :ref: "+utility::to_string( res_pair.second )+" atomno "+ utility::to_string( ref_atom_index )+" "+ pose.residue_type( res_pair.second).
				// utility_exit_with_message()

				if ( (!robust_) && (pose.residue_type( res_pair.first).natoms() != ref_pose_->residue_type( res_pair.second).natoms()) ) {
					std::string msg = "RMSDMetric.  Robust set to false.  Number of atoms in residues do not match. "+ to_string( res_pair.first)+" "+to_string(res_pair.second);
					utility_exit_with_message(msg);
				}
				if ( ref_pose_->residue_type( res_pair.second).has( pose_atom.name()) ) {
					core::Size ref_atom_index = pose.residue_type( res_pair.second).atom_index( pose_atom.name() );
					atom_map[ core::id::AtomID( i, res_pair.first)] = core::id::AtomID( ref_atom_index, res_pair.second);
				} else if ( (!robust_) ) {
					std::string msg = "RMSDMetric.  Robust set to false.  ref pose is missing "+pose_atom.name()+" at residue "+to_string( res_pair.second );
					utility_exit_with_message(msg);
				}

			}
		}
	}
	return atom_map;
}

std::map< core::Size, core::Real >
PerResidueRMSDMetric::calculate(const pose::Pose & pose) const {
	using namespace utility;
	if ( ! ref_pose_ ) {
		utility_exit_with_message( "Must pass in a reference pose for PerResidueRMSDMetric.  See RS XSD or use the set_comparison_pose function");
	}

	std::map< id::AtomID, id::AtomID> atom_map = create_atom_id_map(pose);
	utility::vector1< bool > mask = get_selector()->apply( pose );

	if ( core::pose::symmetry::is_symmetric( pose )  && desymmetrize_res_selector_ )  {
		TR << "De-Symmetrizing selector" << std::endl;
		mask = core::select::get_master_subunit_selection(pose, mask);
	}

	if ( superimpose_ ) {
		core::pose::Pose local_pose = pose;
		superimpose_pose(local_pose, *ref_pose_, atom_map);
		return per_res_rms_at_corresponding_atoms_no_super( local_pose, *ref_pose_, atom_map, mask);
	} else {
		return per_res_rms_at_corresponding_atoms_no_super( pose, *ref_pose_, atom_map, mask);
	}
}

void
PerResidueRMSDMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueRMSDMetric::provide_xml_schema( xsd );
}

std::string
PerResidueRMSDMetricCreator::keyname() const {
	return PerResidueRMSDMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueRMSDMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PerResidueRMSDMetric );

}

} //core
} //simple_metrics
} //per_residue_metrics






