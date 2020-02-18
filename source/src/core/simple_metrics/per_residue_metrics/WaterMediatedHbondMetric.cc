// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/WaterMediatedHbondMetric.cc
/// @brief A metric to measure hydrogen bonds between a set of residues that are water-mediated.  Depth of 1 is default where one water mediates the interface between these residues.  Depth can be set.  Make sure to use the -include_waters flag to have Rosetta not ignore HOH

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/WaterMediatedHbondMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/Pose.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/select/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.WaterMediatedHbondMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::scoring::hbonds;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
WaterMediatedHbondMetric::WaterMediatedHbondMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
WaterMediatedHbondMetric::~WaterMediatedHbondMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
WaterMediatedHbondMetric::WaterMediatedHbondMetric( WaterMediatedHbondMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
WaterMediatedHbondMetric::clone() const {
	return utility::pointer::make_shared< WaterMediatedHbondMetric >( *this );
}

std::string
WaterMediatedHbondMetric::name() const {
	return name_static();
}

std::string
WaterMediatedHbondMetric::name_static() {
	return "WaterMediatedHbondMetric";

}
std::string
WaterMediatedHbondMetric::metric() const {
	return "water_mediated_hbonds";
}

void
WaterMediatedHbondMetric::set_include_virt_waters(bool include_virt_waters){
	include_virt_waters_ = include_virt_waters;
}

void
WaterMediatedHbondMetric::set_residue_selector2(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_two_ = selector;
}

void
WaterMediatedHbondMetric::set_depth(const core::Size depth){
	depth_ = depth;
}

void
WaterMediatedHbondMetric::set_include_self(bool include_self){
	include_self_ = include_self;
}

void
WaterMediatedHbondMetric::set_include_only_set_depth(bool include_only_depth){
	include_only_set_depth_ = include_only_depth;
}

void
WaterMediatedHbondMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );


	set_include_virt_waters(tag->getOption< bool >("include_virt_waters", include_virt_waters_));

	set_include_self(tag->getOption< bool >("include_self", include_self_));

	set_depth(tag->getOption< core::Size >("depth", depth_));

	set_include_only_set_depth(tag->getOption< bool >("include_only_set_depth", include_only_set_depth_));

	if ( tag->hasOption("residue_selector2") ) {
		set_residue_selector2(select::residue_selector::parse_residue_selector( tag, datamap));
	}
}

void
WaterMediatedHbondMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attributes_for_parse_residue_selector( attlist, "residue_selector2", "Optional Selector to measure water-mediated hbonds between residues in each selection, instead of between residues within a single selector.  If NO selector is given, will calculate ALL bridged hbonds between all residues.");

	attlist +XMLSchemaAttribute::attribute_w_default("depth", xsct_positive_integer, "Depth we will go to find hbond (paths). A depth of one (default) means a single water mediates the hydrogen bond(s).", "1");

	attlist + XMLSchemaAttribute::attribute_w_default("include_only_set_depth", xsct_rosetta_bool, "Only calculate hbond paths AT the specific depth instead of UP TO AND INCLUDING that depth.", "false");

	attlist + XMLSchemaAttribute::attribute_w_default("include_self", xsct_rosetta_bool, "Set to include self-self hydrogen bonds: Ex: resJ - water - resJ", "false");

	std::string virt_docs = "After (OptH) water packing, some waters become HOH_V if they don't match certain ref energy.\n"
		"Set this option to true to incude them in our calculation if they are still virts.\n"
		"Note: The option -include_vrt false will have them all be HOH during packing. Use if you are packing water oxygens from a crystal.";

	attlist + XMLSchemaAttribute::attribute_w_default("include_virt_waters", xsct_rosetta_bool, virt_docs, "false");


	std::string documentation =
		"A metric to measure hydrogen bonds between a set of residues that are water-mediated.\n"
		"\n"
		" DEPTH:\n"
		"  Depth is set to a default of 1 (IE a single water mediates the hydrogen bond).  Make sure to set the -ignore_waters flag to false in order to have Rosetta include the HOH residues.\n"
		"\n"
		" SELECTION:\n"
		"  If one residue selector is given, will calculate bridged water paths between residues of the selection and all [OTHER] residues that are not water, otherwise it will calculate bridges between one selection and another.\n"
		"\n"
		"If NO SELECTION is give, will report ALL bridged hbonds in a pose [by default without self-water-self hbonds]."
		"\n"
		" HBONDS:\n"
		"  Since these are bridged hbonds, and h-bond networks can be rather complex, the numbers\n"
		"   reported here are the unique h-bond paths from sele1 to sele2.  If you give only a single residue selector, these are bridged hbonds from sele1 to OTHER residues in the pose not including water. Pass the same selection to get internal bridges of a selection.\n"
		"   By default we do not include mediated hbonds back on itself, but this is an option.\n"
		"\n"
		" TIPS:\n"
		"  It is generally recommended to repack the waters using the OptH TaskOperation before input into this metric, especially if only Oxygens were present.  \n"
		"\n"
		"  Using the option -include_vrt false will keep all waters present in the resulting structure.  Use the option -corrections::water::wat_rot_sampling 10 to decrease the angle of sampling from a default of 30 to 10 or 15.  This will result in many more rotamers (and run-time), but will improve networks.\n"
		"\n"
		"  Hydration shells can be calculated by passing selection1 as all waters in the pose, selection2 as the protein or chain, and then reporting only a single depth (with max_depth at 0 and 1) for the first and second shell waters. Post-processing will be needed."
		"\n"
		"  Finally, it is generally recommended to use the -beta_nov16 and -genpot scorefunctions as this can result in better h-bond detection.  Both of these have been published and can be turned on together using the -beta flag.\n"
		"\n"
		" AUTHORS:\n"
		"  Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"   Citation: De-Novo Glycan Modeling in Rosetta (drafting)\n"
		"\n"
		"  Ryan Pavlovicz (Rosetta-ECO: Explicit Consideration of cOordinated water)"
		"   Citation: 'Efficient consideration of coordinated water molecules improves computational protein- protein and protein-ligand docking', biorxiv, PLOS CB in-review";

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		documentation, attlist);
}

///@brief Recursive function to find all unique HBond paths involving water and ending with
///  a residue that is not in our water list.
void
WaterMediatedHbondMetric::find_hb_paths(
	HBondSet const & local_hb_set,
	utility::vector1< core::Size > const & local_waters,
	std::map<std::string, core::Size> & paths,
	core::Size const current_res,
	core::Size const max_depth /*1*/,
	std::string const & current_path /*""*/,
	core::Size const current_depth /*0*/,
	HBondCOP prev_hb /*nullptr*/) const {

	//JAB - Moved to hbond util for general use
	core::scoring::hbonds::find_hb_paths(local_hb_set, local_waters, paths, current_res, max_depth, current_path, current_depth, prev_hb);
}


std::map< core::Size, core::Real >
WaterMediatedHbondMetric::calculate(pose::Pose const & pose) const {
	core::pose::Pose local_pose = pose;

	std::map<core::Size, core::Real> bridged_hbonds;

	//Initialize map to all selection as 0
	utility::vector1< bool > const mask = get_selector()->apply( pose );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( mask[i] ) {
			bridged_hbonds[i] = 0;
		}
	}

	utility::vector1< core::Size > const sele1 = get_residues_from_subset(mask);
	utility::vector1< core::Size > sele2;

	if ( selector_two_ ) {
		utility::vector1< bool > const mask2 = selector_two_->apply(pose);
		sele2 = get_residues_from_subset(mask2);
	} else {
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( (! sele1.contains(i)) || (sele1.size() == pose.size()) ) {
				sele2.push_back(i);
			}
		}
	}

	//Create HbondSet
	//NEED copy of pose here...
	HBondSet hb_set = HBondSet(local_pose, false);

	//Find waters we care about.
	utility::vector1< core::Size > waters;

	for ( core::Size i = 1; i<= pose.size(); ++i ) {
		if ( pose.residue_type(i).is_water() ) {
			if ( ! include_virt_waters_ && pose.residue(i).name() == "HOH_V" ) continue;
			waters.push_back(i);
		}
	}

	for ( core::Size i : sele1 ) {
		std::map< std::string, core::Size > all_paths;
		find_hb_paths(hb_set, waters, all_paths, i,  depth_);
		for ( auto const & path_pair : all_paths ) {

			utility::vector1< std::string > sp = utility::string_split(path_pair.first, '-');
			core::Size depth = sp.size() - 2;
			if ( depth != depth_ && include_only_set_depth_ ) continue;

			//core::Size s1 = utility::string2Size(sp[0]);
			core::Size s2 = utility::string2Size(sp[sp.size()]);
			TR.Debug << path_pair.first << std::endl;
			//Bridged h-bond found at depth
			if ( sele2.contains(s2) ) {

				if ( s2 == i && ! include_self_ ) {
					continue;
				}

				bridged_hbonds[i]+=1;
			}
		}
	}
	return bridged_hbonds;

} //End calculate

void
WaterMediatedHbondMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	WaterMediatedHbondMetric::provide_xml_schema( xsd );
}

std::string
WaterMediatedHbondMetricCreator::keyname() const {
	return WaterMediatedHbondMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
WaterMediatedHbondMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< WaterMediatedHbondMetric >();
}

} //per_residue_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION


template< class Archive >
void
core::simple_metrics::per_residue_metrics::WaterMediatedHbondMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( include_virt_waters_ ));
	arc( CEREAL_NVP( include_self_ ));
	arc( CEREAL_NVP( include_only_set_depth_ ));
	arc( CEREAL_NVP( depth_ ) );
	arc( CEREAL_NVP( selector_two_) );

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::WaterMediatedHbondMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( include_virt_waters_);
	arc(include_self_);
	arc(include_only_set_depth_);
	arc(depth_);
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_two_ = local_selector;


}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::WaterMediatedHbondMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::WaterMediatedHbondMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_WaterMediatedHbondMetric )
#endif // SERIALIZATION




