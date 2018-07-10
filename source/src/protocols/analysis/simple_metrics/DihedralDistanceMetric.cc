// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/DihedralDistanceMetric.cc
/// @brief A metric to calculate the dihedral distance between two poses or the input and the set cmd-line native.  Can set a subset of residues to calculate via ResidueSelector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/analysis/simple_metrics/DihedralDistanceMetric.hh>
#include <protocols/analysis/simple_metrics/DihedralDistanceMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/import_pose/import_pose.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/options/OptionCollection.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <numeric/NumericTraits.hh>

#include <cmath>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.DihedralDistanceMetric" );


namespace protocols {
namespace analysis {
namespace simple_metrics {

using namespace core::pose;
using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::simple_metrics;
using namespace protocols::rosetta_scripts;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
DihedralDistanceMetric::DihedralDistanceMetric():
	core::simple_metrics::RealMetric()
{

}

DihedralDistanceMetric::DihedralDistanceMetric( ResidueSelectorCOP selector ):
	core::simple_metrics::RealMetric()
{
	set_residue_selector( selector);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
DihedralDistanceMetric::~DihedralDistanceMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
DihedralDistanceMetric::DihedralDistanceMetric( DihedralDistanceMetric const & src ):
	RealMetric( src ),
	include_protein_omega_(src.include_protein_omega_),
	res_map_(src.res_map_)
{
	residue_selector_ = src.residue_selector_;
	residue_selector_ref_ = src.residue_selector_ref_;
	ref_pose_ = ref_pose_;
}


core::simple_metrics::SimpleMetricOP
DihedralDistanceMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new DihedralDistanceMetric( *this ) );

}

std::string
DihedralDistanceMetric::name() const {
	return name_static();
}

std::string
DihedralDistanceMetric::name_static() {
	return "DihedralDistanceMetric";

}
std::string
DihedralDistanceMetric::metric() const {
	return "dihedral_distance";

}

void
DihedralDistanceMetric::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ = residue_selector;
}

void
DihedralDistanceMetric::set_residue_selector_reference( core::select::residue_selector::ResidueSelectorCOP residue_selector){
	residue_selector_ref_ = residue_selector;
}

void
DihedralDistanceMetric::set_residue_mapping(std::map< core::Size, core::Size> const & res_map){
	res_map_ = res_map;
}

void
DihedralDistanceMetric::set_comparison_pose(const core::pose::Pose &pose){
	ref_pose_ = PoseOP( new Pose( pose ));
}

void
DihedralDistanceMetric::set_include_protein_omega(bool include_omega){
	include_protein_omega_ = include_omega;
}

void
DihedralDistanceMetric::load_native_pose_as_reference(){
	if ( option[ OptionKeys::in::file::native ].user() ) {
		ref_pose_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() , core::import_pose::PDB_file);
	}
}

void
DihedralDistanceMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	SimpleMetric::parse_base_tag( tag );
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
		std::string msg = "A reference pose must be set. Please use the SavePoseMover (embed the RMSDMetric in RunSimpleMetrics) or pass the native as in:file:native and set use_native to true.";
		utility_exit_with_message(msg);
	}

	set_include_protein_omega(tag->getOption< bool >("include_omege", include_protein_omega_) );

}

void
DihedralDistanceMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;
	using namespace protocols::rosetta_scripts;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"include_omega", xsct_rosetta_bool, "Do we include protein omega in calculations?  Default false.", "false"
	);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "false"
	);

	attributes_for_saved_reference_pose( attlist );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Calculate the RMSD for these residues for both reference and main pose." );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector_ref",
		"Selector for the reference pose (input native or passed reference pose. ).  Residues selected must be same number of residues selected for the main selector." );

	std::string description =  "Calculate DihedralDistance metric.\n"
		"  This is the normalized metric in degrees\n"
		"  Works on protein and carbohydrate backbones\n"
		"\n"
		" Metric described in:\n"
		"   North B, Lehmann A, Dunbrack RL. A new clustering of antibody CDR loop conformations. J Mol Biol 2011;406:228-256.\n";

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(), description, attlist);
}

core::Real
DihedralDistanceMetric::calculate(const core::pose::Pose & pose) const {
	using namespace utility;
	using namespace core::pose::carbohydrates;

	if ( ! ref_pose_ ) {
		utility_exit_with_message( "Must pass in a reference pose for RMSDMetric.  See RS XSD or use the set_comparison_pose function");
	}

	//Setup the main residue mapping we will use depending on what is set.
	std::map< core::Size, core::Size > residue_map;
	if ( res_map_.empty() ) {
		residue_map = get_residue_mapping_from_selectors( residue_selector_, residue_selector_ref_, pose, *ref_pose_);
	} else {
		residue_map = res_map_;
	}

	//core::Size n_res = residue_map.size();

	utility::vector1< core::Real > pose_dihedrals;
	utility::vector1< core::Real > ref_pose_dihedrals;

	for ( auto res_pair : residue_map ) {
		if ( pose.residue_type( res_pair.first ).is_carbohydrate() && ref_pose_->residue_type( res_pair.second ).is_carbohydrate() ) {
			core::Size pose_torsions = get_n_glycosidic_torsions_in_res( pose, res_pair.first);
			core::Size ref_torsions = get_n_glycosidic_torsions_in_res( *ref_pose_ , res_pair.second);
			if ( pose_torsions != ref_torsions ) {
				std::string msg = "DihedralDistanceMetric: number of carbohydrate torsions do not match at positions: "+utility::to_string(res_pair.first)+" "+utility::to_string( res_pair.second);
				utility_exit_with_status(msg);
			}
			for ( core::Size i = 1; i <= pose_torsions; ++i ) {
				pose_dihedrals.push_back( get_glycosidic_torsion( i, pose, res_pair.first));
				ref_pose_dihedrals.push_back( get_glycosidic_torsion(i, *ref_pose_, res_pair.second));
			}

		} else if ( pose.residue_type( res_pair.first).is_protein() && ref_pose_->residue_type( res_pair.second).is_protein() ) {
			pose_dihedrals.push_back( pose.phi( res_pair.first ) );
			pose_dihedrals.push_back( pose.psi( res_pair.first ) );

			ref_pose_dihedrals.push_back( ref_pose_->phi( res_pair.second ) );
			ref_pose_dihedrals.push_back( ref_pose_->psi( res_pair.second ) );

			if ( include_protein_omega_ ) {
				pose_dihedrals.push_back(pose.omega( res_pair.first ));
				ref_pose_dihedrals.push_back( ref_pose_->omega( res_pair.second ));
			}
		} else {
			TR << "Skipping non-protein and non-carbohydrate residue: " << res_pair.first << std::endl;
		}
	}

	core::Real k_distance = 0.0;
	core::Real PI = numeric::NumericTraits<core::Real>::pi();
	//std::cout << "Pose Size: " << pose.size() << std::endl;
	//std::cout << "Dih  Size: " << pose_dihedrals.size() << std::endl;

	for ( Size i=1; i <= pose_dihedrals.size(); ++i ) {

		core::Real torsion_d = (2 * (1- cos ((pose_dihedrals[i]-ref_pose_dihedrals[i])*PI/180)));//Speed
		//std::cout << "Torsion: "<< i << " : " << torsion_d << std::endl;
		//core::Real torsion_d = numeric::model_quality::calculate_dihedral_distance( pose_dihedrals[i], ref_pose_dihedrals[i]);
		k_distance += torsion_d;

	}

	//Calculate normalized
	core::Real normalized_distance = k_distance/pose_dihedrals.size();
	core::Real normalized_in_degrees = acos(1 - (normalized_distance/2)) *(180/numeric::NumericTraits<core::Real>::pi());
	return normalized_in_degrees;
}



void
DihedralDistanceMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	DihedralDistanceMetric::provide_xml_schema( xsd );
}

std::string
DihedralDistanceMetricCreator::keyname() const {
	return DihedralDistanceMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
DihedralDistanceMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new DihedralDistanceMetric );

}

} //core
} //simple_metrics
} //metrics






