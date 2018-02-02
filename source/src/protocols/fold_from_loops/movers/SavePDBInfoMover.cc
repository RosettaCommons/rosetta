// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/SavePDBInfoMover.cc
/// @brief Adds labels to residues selected through a ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/SavePDBInfoMover.hh>
#include <protocols/fold_from_loops/movers/SavePDBInfoMoverCreator.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.SavePDBInfoMover", basic::t_trace );

namespace protocols {
namespace fold_from_loops {
namespace movers {

SavePDBInfoMover::SavePDBInfoMover():
	protocols::moves::Mover( mover_name() )
{}

SavePDBInfoMover::~SavePDBInfoMover()= default;

void
SavePDBInfoMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	core::pose::PDBInfoOP refinfo(nullptr);
	std::string refinfo_name(tag->getOption<std::string>( "reference_info_name" ) );

	if ( !data.has("spm_ref_info", refinfo_name ) ) {
		refinfo = core::pose::PDBInfoOP( new core::pose::PDBInfo() );
		data.add("spm_ref_info", refinfo_name, refinfo );
	} else refinfo = data.get_ptr<core::pose::PDBInfo>("spm_ref_info", refinfo_name );

	pdb_info_     = refinfo;
	TR.Trace << TR.Green << "pdbinfo pointer address " << pdb_info_ << TR.Reset << std::endl;
	restore_info_ = tag->getOption<bool>("restore_info", true);

}

protocols::moves::MoverOP
SavePDBInfoMover::clone() const
{
	return protocols::moves::MoverOP( new SavePDBInfoMover( *this ) );
}

protocols::moves::MoverOP
SavePDBInfoMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SavePDBInfoMover );
}

void
SavePDBInfoMover::apply( core::pose::Pose & pose )
{
	if ( !restore_info_ ) {
		runtime_assert_msg( pose.pdb_info(), "There is no PDBInfo in the Pose to store." );
		TR << "Storing PDBInfo" << std::endl;
		core::pose::PDBInfoOP tmp = pose.pdb_info();
		*pdb_info_ = *tmp;
		pdb_info_->show( TR );
	} else {
		TR << "Current PDBInfo" << std::endl;
		if ( pose.pdb_info() ) {
			pose.pdb_info()->show( TR );
		} else {
			TR << "No PDBInfo available." << std::endl;
		}
		TR << "Loading PDBInfo" << std::endl;
		pdb_info_->show( TR );
		pose.pdb_info( pdb_info_ );
		pose.pdb_info()->obsolete(false);
	}
}

std::string SavePDBInfoMover::get_name() const {
	return mover_name();
}

std::string SavePDBInfoMover::mover_name() {
	return "SavePDBInfoMover";
}

void SavePDBInfoMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "reference_info_name", xs_string, "Name to assign to the given PDBInfo." )
		+ XMLSchemaAttribute::attribute_w_default( "restore_info", xsct_rosetta_bool, "If False, store the PDBInfo; if True, place it back on the pose.", std::to_string( true ) );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Saves/Restores the PDBInfo of a Pose.", attlist );
}

std::string SavePDBInfoMoverCreator::keyname() const {
	return SavePDBInfoMover::mover_name();
}

protocols::moves::MoverOP
SavePDBInfoMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SavePDBInfoMover );
}

void SavePDBInfoMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SavePDBInfoMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
