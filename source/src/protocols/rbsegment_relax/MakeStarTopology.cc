// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/rbsegment_relax/MakeStarTopology.hh>
#include <protocols/rbsegment_relax/MakeStarTopologyCreator.hh>

#include <protocols/loops/loops_main.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ref_pose.hh>

#include <protocols/rbsegment_relax/util.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>



namespace protocols {
namespace rbsegment_relax {

static basic::Tracer TR( "protocols.cryst.cryst_movers" );

using namespace protocols;
using namespace core;


//////////////////
//////////////////
/// creator






//////////////////
//////////////////
/// mover

void MakeStarTopologyMover::apply( core::pose::Pose & pose ) {
	using namespace protocols::rbsegment_relax;

	if ( restore_ ) {
		if ( ft_restore_->empty() && input_pose_ != nullptr ) {
			// First time around, grab foldtree from input.
			*ft_restore_ = input_pose_->fold_tree();
		}
		pose.fold_tree( *ft_restore_ );
		// ??? cutpoint variants
		protocols::loops::remove_cutpoint_variants( pose );
	} else {
		*ft_restore_ = pose.fold_tree();
		if ( mode_ == "disconnected" ) {
			setup_disconnected( pose );
		} else {
			//default
			setup_star_topology( pose );
		}
	}
}

void MakeStarTopologyMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	mode_ = tag->getOption<std::string>("mode", "default");


	restore_ = tag->getOption<bool>("restore", false);
	tag_ = tag->getOption<std::string>("tag", "nulltag");

	// look for tag on datamap
	if ( data.has( "foldtrees", tag_ ) ) {
		ft_restore_ = data.get_ptr<core::kinematics::FoldTree>( "foldtrees", tag_ );
		TR << "Found foldtree " << tag_ << " on datamap" << std::endl;
	} else {
		// Create a new FoldTree on the datamap
		ft_restore_ = utility::pointer::make_shared< core::kinematics::FoldTree >();
		data.add( "foldtrees", tag_, ft_restore_ );
		TR << "Adding foldtrees " << tag_ << " to datamap" << std::endl;
		// Get a reference to the input pose for later reference if we need it.
		input_pose_ = protocols::rosetta_scripts::legacy_saved_pose_or_input( tag, data, mover_name() );
	}
}

std::string MakeStarTopologyMover::get_name() const {
	return mover_name();
}

std::string MakeStarTopologyMover::mover_name() {
	return "MakeStarTopology";
}

void MakeStarTopologyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction valid_modes;
	valid_modes.name( "valid_make_star_topology_modes" );
	valid_modes.base_type( xs_string );
	valid_modes.add_restriction( xsr_enumeration, "default" );
	valid_modes.add_restriction( xsr_enumeration, "disconnected" );
	xsd.add_top_level_element( valid_modes );

	//default, disconnected
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "mode", "valid_make_star_topology_modes", "Set up default or disconnected star topology?", "default" )
		+ XMLSchemaAttribute::attribute_w_default( "restore", xsct_rosetta_bool, "Apply the fold tree set up on a previous run of this mover to the pose?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "tag", xs_string, "Name of a previously defined fold tree to be used in this mover. If the fold tree has not been defined, it will be set as the default fold tree.", "nulltag" );

	core::pose::attributes_for_saved_reference_pose_w_description( attlist, "The saved pose to take the forldtree from if restoring without a prior save." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sets up a star topology fold tree for this pose", attlist );
}

std::string MakeStarTopologyMoverCreator::keyname() const {
	return MakeStarTopologyMover::mover_name();
}

protocols::moves::MoverOP
MakeStarTopologyMoverCreator::create_mover() const {
	return utility::pointer::make_shared< MakeStarTopologyMover >();
}

void MakeStarTopologyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeStarTopologyMover::provide_xml_schema( xsd );
}


}
}
