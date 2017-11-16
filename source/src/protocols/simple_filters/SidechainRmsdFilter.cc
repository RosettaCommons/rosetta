// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SidechainRmsdFilter.cc
/// @brief A filter based on automorphic sidechain RMSD
/// @author Noah Ollikainen

#include <protocols/simple_filters/SidechainRmsdFilter.hh>
#include <protocols/simple_filters/SidechainRmsdFilterCreator.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/rms_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer sidechain_rmsd_filter_tracer( "protocols.simple_filters.SidechainRmsdFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP SidechainRmsdFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SidechainRmsdFilter ); }

// XRW TEMP std::string
// XRW TEMP SidechainRmsdFilterCreator::keyname() const { return "SidechainRmsd"; }

SidechainRmsdFilter::SidechainRmsdFilter() : filters::Filter( "SidechainRmsd"  ) {}

SidechainRmsdFilter::SidechainRmsdFilter( std::string const & res1, std::string const & res2, core::Real const rmsd_threshold ) :
	Filter( "SidechainRmsd" ), res1_( res1 ), res2_( res2 ), rmsd_threshold_( rmsd_threshold ) {}

SidechainRmsdFilter::~SidechainRmsdFilter()= default;

void
SidechainRmsdFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	res1_ = core::pose::get_resnum_string( tag, "res1_" );
	res2_ = core::pose::get_resnum_string( tag, "res2_" );
	rmsd_threshold_ = tag->getOption<core::Real>("threshold", 1.0);
	include_backbone_ = tag->getOption<bool>("include_backbone", false);

	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map);
	} else {
		reference_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			core::import_pose::pose_from_file( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		}
	}

	//residue_distance_filter_tracer<<"ResidueDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<<res1_<<" and "<<res2_<<std::endl;
}

bool
SidechainRmsdFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ) );

	sidechain_rmsd_filter_tracer<<"Sidechain RMSD of residue "<<res1_<<" is "<<rmsd<<std::endl;
	return( rmsd<=rmsd_threshold_ );
}

void
SidechainRmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ) );

	out<<"Sidechain RMSD of residue "<<res1_<<" is "<<rmsd<<'\n';
}

core::Real
SidechainRmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ) );

	return( rmsd );
}
core::Real
SidechainRmsdFilter::compute( core::pose::Pose const & pose ) const {
	debug_assert( reference_pose_ );
	core::Size res1( core::pose::parse_resnum( res1_, pose ) );
	core::Size res2( core::pose::parse_resnum( res2_, pose ) );
	core::conformation::Residue const res_res1( pose.conformation().residue( res1 ) );
	core::conformation::Residue const res_res2( reference_pose_->conformation().residue( res2 ) );
	core::Real rmsd (0.0);

	// make sure we're comparing the same amino acid type
	runtime_assert( res_res1.aa() == res_res2.aa() );

	if ( include_backbone_ ) {
		rmsd = core::scoring::automorphic_rmsd( res_res1, res_res2, false /*superimpose*/ );
	} else {
		core::chemical::ResidueTypeSetCOP res1_set( pose.residue_type_set_for_pose( res_res1.type().mode() ) );
		core::chemical::ResidueType const & working_res1_type(
			res1_set->get_residue_type_with_variant_added( res_res1.type(), core::chemical::VIRTUAL_BB ) );
		core::conformation::ResidueOP working_res1 =
			core::conformation::ResidueFactory::create_residue( working_res1_type );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms(
			res_res1, *working_res1, pose.conformation() );

		core::chemical::ResidueTypeSetCOP res2_set( reference_pose_->residue_type_set_for_pose( res_res2.type().mode() ) );
		core::chemical::ResidueType const & working_res2_type(
			res2_set->get_residue_type_with_variant_added( res_res2.type(), core::chemical::VIRTUAL_BB ) );
		core::conformation::ResidueOP working_res2 =
			core::conformation::ResidueFactory::create_residue( working_res2_type );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms(
			res_res2, *working_res2, reference_pose_->conformation() );

		rmsd = core::scoring::automorphic_rmsd(*working_res1, *working_res2, false /*superimpose*/);
	}

	return( rmsd );
}

filters::FilterOP SidechainRmsdFilter::clone() const {
	return filters::FilterOP( new SidechainRmsdFilter( *this ) );
}

filters::FilterOP SidechainRmsdFilter::fresh_instance() const{
	return filters::FilterOP( new SidechainRmsdFilter() );
}

std::string SidechainRmsdFilter::name() const {
	return class_name();
}

std::string SidechainRmsdFilter::class_name() {
	return "SidechainRmsd";
}

void SidechainRmsdFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold"  , xsct_real , "In a truth value context, what's the maximum RMSD value which is considered to be passing." , "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "include_backbone" , xsct_rosetta_bool , "Whether to include the backbone in the RMSD calculation. It is recommended to set this to 'true' for ligands and other residues which don't have a backbone." , "false" ) ;

	core::pose::attributes_for_get_resnum_string( attlist , "res1_" ) ;
	//The residue number for the active pose. see RosettaScripts#rosettascripts-conventions_specifying-residues
	core::pose::attributes_for_get_resnum_string( attlist , "res2_" ) ;
	//The residue number for the reference pose. see RosettaScripts#rosettascripts-conventions_specifying-residues

	protocols::rosetta_scripts::attributes_for_saved_reference_pose( attlist , "reference_name" ) ;
	//The name of the reference pose as saved with the SavePoseMover . If not given, will default to the structure passed to -in:file:native if it set, or the input structure, if not.

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates the all atom RMSD for a single residue, either with or without the backbone atoms. The RMSD calculated is the automorphic RMSD, so it will compensate for symmetric rearrangments. (For example, Phe ring flips.) No superposition is performed prior to rmsd calculation.", attlist );
}

std::string SidechainRmsdFilterCreator::keyname() const {
	return SidechainRmsdFilter::class_name();
}

protocols::filters::FilterOP
SidechainRmsdFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SidechainRmsdFilter );
}

void SidechainRmsdFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SidechainRmsdFilter::provide_xml_schema( xsd );
}


}
}
