// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/RMSDMetric.cc
/// @brief A metric to calculate the RMSD between two poses or the input and the set cmd-line native.  Can set a subset of residues to calculate via ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/analysis/simple_metrics/RMSDMetric.hh>
#include <protocols/analysis/simple_metrics/RMSDMetricCreator.hh>
#include <core/scoring/rms_util.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/AtomID.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.RMSDMetric" );


namespace protocols {
namespace analysis {
namespace simple_metrics {

using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::pose;
using namespace core::simple_metrics;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::rosetta_scripts;

/////////////////////
/// Constructors  ///
/////////////////////


std::map< rmsd_atoms, utility::vector1< std::string > >
setup_rmsd_atom_names(){

	std::map< rmsd_atoms, utility::vector1< std::string > > atom_names;
	for ( core::Size i = 1; i <= core::Size(rmsd_all); ++i ) {
		rmsd_atoms rmsd_type = static_cast<rmsd_atoms>(i);
		utility::vector1< std::string > name_list;
		atom_names[ rmsd_type ] = name_list;
	}

	//rmsd_bb_heavy
	atom_names[ rmsd_protein_bb_heavy ].push_back( " N  ");
	atom_names[ rmsd_protein_bb_heavy ].push_back( " CA ");
	atom_names[ rmsd_protein_bb_heavy ].push_back( " C  ");

	//rmsd_bb_heavy_including_0
	atom_names[ rmsd_protein_bb_heavy_including_O] = atom_names[ rmsd_protein_bb_heavy];
	atom_names[ rmsd_protein_bb_heavy_including_O].push_back( " O  ");

	//NEED TO ADD GLYCAN BB HERE

	//rmsd_bb_ca_only
	atom_names[ rmsd_protein_bb_ca ].push_back( " CA ");

	//OTHERS we will do as we go along and create the atom_id_map.

	return atom_names;
}


/// @brief Default constructor
RMSDMetric::RMSDMetric():
	core::simple_metrics::RealMetric()
{
	setup_name_mapping();
}

RMSDMetric::RMSDMetric( Pose const & ref_pose):
	core::simple_metrics::RealMetric()
{
	set_comparison_pose( ref_pose );
	setup_name_mapping();
}

RMSDMetric::RMSDMetric( Pose const & ref_pose, ResidueSelectorCOP selector ):
	core::simple_metrics::RealMetric()
{
	set_comparison_pose( ref_pose );
	set_residue_selector( selector );
	setup_name_mapping();
}

std::map< std::string, rmsd_atoms >
get_rmsd_type_name_map(){
	std::map< std::string, rmsd_atoms > name_mapping;

	name_mapping["rmsd_protein_bb_heavy"] = rmsd_protein_bb_heavy;
	name_mapping["rmsd_protein_bb_heavy_including_O"] = rmsd_protein_bb_heavy_including_O;
	name_mapping["rmsd_protein_bb_ca"] = rmsd_protein_bb_ca;
	name_mapping["rmsd_sc_heavy"] = rmsd_sc_heavy;
	name_mapping["rmsd_sc"] = rmsd_sc;
	name_mapping["rmsd_all_heavy"] = rmsd_all_heavy;
	name_mapping["rmsd_all"] = rmsd_all;
	return name_mapping;
}

utility::vector1< std::string >
get_rmsd_type_names(){
	std::map< std::string, rmsd_atoms > name_mapping = get_rmsd_type_name_map();
	utility::vector1< std::string > names;
	for ( auto name_pair : name_mapping ) {
		names.push_back(name_pair.first);
	}
	return names;
}

void
RMSDMetric::setup_name_mapping(){
	rmsd_atom_names_ = setup_rmsd_atom_names();
	name_mapping_ = get_rmsd_type_name_map();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
RMSDMetric::~RMSDMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
RMSDMetric::RMSDMetric( RMSDMetric const & src ):
	RealMetric( src ),
	rmsd_map_(src.rmsd_map_),
	rmsd_type_(src.rmsd_type_),
	override_atom_names_( src.override_atom_names_ ),
	rmsd_atom_names_( src.rmsd_atom_names_),
	robust_( src.robust_),
	name_mapping_( src.name_mapping_ )
{
	residue_selector_ = src.residue_selector_;
	residue_selector_ref_ = src.residue_selector_ref_;
	ref_pose_ = src.ref_pose_;
}


std::string
RMSDMetric::name() const {
	return name_static();
}

std::string
RMSDMetric::name_static() {
	return "RMSDMetric";

}
std::string
RMSDMetric::metric() const {
	return "rmsd";
}

void
RMSDMetric::set_rmsd_type( rmsd_atoms rmsd_type){
	rmsd_type_ = rmsd_type;
}

void
RMSDMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ = residue_selector;
}

void
RMSDMetric::set_residue_selector_reference(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ref_ = residue_selector;
}

void
RMSDMetric::set_residue_mapping(std::map<core::Size, core::Size> const & rmsd_map ){
	rmsd_map_ = rmsd_map;
}

void
RMSDMetric::set_comparison_pose(const core::pose::Pose &pose){
	ref_pose_ = PoseOP( new Pose( pose ));
}

void
RMSDMetric::set_corresponding_atoms_robust(bool robust){
	robust_ = robust;
}

core::simple_metrics::SimpleMetricOP
RMSDMetric::clone() const {
	return SimpleMetricOP(new RMSDMetric( *this ) );

}

void
RMSDMetric::load_native_pose_as_reference(){
	if ( option[ OptionKeys::in::file::native ].user() ) {
		ref_pose_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() , core::import_pose::PDB_file);
	}
}

void
RMSDMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(protocols::rosetta_scripts::parse_residue_selector( tag, datamap ));
	}

	if ( tag->hasOption("residue_selector_ref") ) {
		set_residue_selector_reference(protocols::rosetta_scripts::parse_residue_selector( tag, datamap, "residue_selector_ref" ));
	}


	//Comparison pose.
	if ( tag->hasOption("reference_name") ) {
		ref_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, datamap, "reference_name");
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<ref_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption< bool >("use_native", false) ) {
		load_native_pose_as_reference();
	} else {
		std::string msg = "A reference pose must be set. Please use the SavePoseMover (embed the RMSDMetric in RunSimpleMetrics ) or pass the native as in:file:native and set use_native to true.";
		utility_exit_with_message(msg);
	}

	if ( tag->hasOption("rmsd_type") ) {
		set_rmsd_type( name_mapping_[ tag->getOption<std::string>("rmsd_type")]);
	}

	set_corresponding_atoms_robust(tag->getOption< bool >("robust", robust_));
}



void
RMSDMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;
	using namespace protocols::rosetta_scripts;

	AttributeList attlist;
	attributes_for_saved_reference_pose( attlist );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Calculate the RMSD for these residues for both reference and main pose." );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector_ref",
		"Selector for the reference pose (input native or passed reference pose. ).  Residues selected must be same number of residues selected for the main selector." );


	attlist + XMLSchemaAttribute::attribute_w_default(
		"robust", xsct_rosetta_bool, "Set whether we are robust to atom mismatches for selected residues."
		"  By default we only match atoms that are corresponding. (True).", "true"
	);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false"
	);

	utility::vector1< std::string > rmsd_type_names = get_rmsd_type_names();
	utility::tag::add_schema_restrictions_for_strings( xsd, "rmsd_types", rmsd_type_names);

	attlist + XMLSchemaAttribute("rmsd_type", "rmsd_types", "Type of calculation.  Current choices are: \n" + utility::to_string(rmsd_type_names) );

	std::string description = "\tThis is the RMSD between the input and the set comparison pose.\n"
		"  If native is set on the cmd-line, we will use that.\n"
		"  Default is to calculate all_heavy atoms - but this can be set\n."
		"\n"
		"  Make sure that reference pose and set pose are the same length ."
		"   We match all corresponding atoms for each residue to match."
		"   By default we do not fail and are robust - only matching what we can for each residue.";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		description, attlist);
}

core::Real
RMSDMetric::calculate(const core::pose::Pose & pose) const {

	using namespace utility;
	if ( ! ref_pose_ ) {
		utility_exit_with_message( "Must pass in a reference pose for RMSDMetric.  See RS XSD or use the set_comparison_pose function");
	}

	//Setup the main residue mapping we will use depending on what is set.
	std::map< core::Size, core::Size > residue_map;
	if ( rmsd_map_.empty() ) {
		residue_map = get_residue_mapping_from_selectors( residue_selector_, residue_selector_ref_, pose, *ref_pose_);
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
	return rms_at_corresponding_atoms_no_super( pose, *ref_pose_, atom_map);
}



void
RMSDMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RMSDMetric::provide_xml_schema( xsd );
}

std::string
RMSDMetricCreator::keyname() const {
	return RMSDMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
RMSDMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new RMSDMetric );

}

} //core
} //simple_metrics
} //metrics






