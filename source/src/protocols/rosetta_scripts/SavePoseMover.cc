// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/SavePoseMover.cc
/// @author Florian Richter (floric@u.washington.edu)

// Unit Headers
#include <protocols/rosetta_scripts/SavePoseMover.hh>
#include <protocols/rosetta_scripts/SavePoseMoverCreator.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.SavePoseMover" );

namespace protocols {
namespace rosetta_scripts {

SavePoseMover::SavePoseMover() :
	Mover( "SavePoseMover" ),
	reference_pose_(/* NULL */),
	restore_pose_(false)
{
}

SavePoseMover::~SavePoseMover() = default;

void
SavePoseMover::apply( core::pose::Pose & pose )
{
	if ( !restore_pose_ ) *reference_pose_ = pose;
	else pose = *reference_pose_;

	//make sure this always counts as success
	this->set_last_move_status( protocols::moves::MS_SUCCESS );
}

void
SavePoseMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = saved_reference_pose(tag,data_map );
	} else throw utility::excn::EXCN_RosettaScriptsOption("Need to specify name under which to save pose.");

	if ( tag->hasOption( "pdb_file" ) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "pdb_file" ));
		core::import_pose::pose_from_file( *reference_pose_, template_pdb_fname , core::import_pose::PDB_file);
		TR <<"reading in " << template_pdb_fname << " pdb with " << reference_pose_->size() <<" residues"<<std::endl;
	}


	if ( tag->hasOption("restore_pose") ) {
		restore_pose_ = tag->getOption<bool>("restore_pose",1);
	}
}

protocols::moves::MoverOP
SavePoseMover::clone() const{
	return protocols::moves::MoverOP( new SavePoseMover( *this ) );
}

protocols::moves::MoverOP
SavePoseMover::fresh_instance() const{
	return protocols::moves::MoverOP( new SavePoseMover );
}

std::string SavePoseMover::get_name() const {
	return mover_name();
}

std::string SavePoseMover::mover_name() {
	return "SavePoseMover";
}

void SavePoseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	// I would consider a complex type that ends in .pdb, but I want users to be
	// able to provide cifs! Or name their files whatever!
	attlist + XMLSchemaAttribute("pdb_file", xs_string, "PDB file to be read in as a reference pose." )
		+ XMLSchemaAttribute::attribute_w_default("restore_pose", xsct_rosetta_bool, "Restore the specified reference pose INTO the current pose rather than saving the current pose INTO the reference pose.", "1" );

	protocols::rosetta_scripts::attributes_for_saved_reference_pose( attlist , "reference_name" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SavePoseMoverCreator::keyname() const {
	return SavePoseMover::mover_name();
}

protocols::moves::MoverOP
SavePoseMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SavePoseMover );
}

void SavePoseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SavePoseMover::provide_xml_schema( xsd );
}



} // moves
} // protocols
