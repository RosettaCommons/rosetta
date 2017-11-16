// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyAddStartnodeFragments.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyAddStartnodeFragments.hh>
#include <protocols/legacy_sewing/sampling/LegacyAddStartnodeFragmentsCreator.hh>
#include <protocols/legacy_sewing/conformation/DisembodiedAssembly.hh>

//Package headers
#include <protocols/legacy_sewing/util/io.hh>
#include <protocols/legacy_sewing/util/util.hh>

#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/rosetta_scripts/util.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/legacy_sewing.OptionKeys.gen.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace legacy_sewing  {

static basic::Tracer TR( "protocols.legacy_sewing.LegacyAddStartnodeFragments" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LegacyAddStartnodeFragmentsCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new LegacyAddStartnodeFragments );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LegacyAddStartnodeFragmentsCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LegacyAddStartnodeFragments::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LegacyAddStartnodeFragments::mover_name()
// XRW TEMP {
// XRW TEMP  return "LegacyAddStartnodeFragments";
// XRW TEMP }

protocols::moves::MoverOP
LegacyAddStartnodeFragments::clone() const {
	return( protocols::moves::MoverOP( new LegacyAddStartnodeFragments( *this ) ) );
}
protocols::moves::MoverOP
LegacyAddStartnodeFragments::fresh_instance() const {
	return protocols::moves::MoverOP( new LegacyAddStartnodeFragments );
}

// XRW TEMP std::string
// XRW TEMP LegacyAddStartnodeFragments::get_name() const {
// XRW TEMP  return "LegacyAddStartnodeFragments";
// XRW TEMP }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  LegacyAddStartnodeFragments function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


LegacyAddStartnodeFragments::LegacyAddStartnodeFragments():
	Mover()
{}


void
LegacyAddStartnodeFragments::apply(
	core::pose::Pose & pose
) {
	using namespace core;
	using namespace core::fragment;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !option[ in::file::native ].user() || !option[in::file::frag9].user() ) {
		utility_exit_with_message("Must specify native and frag9 files for LegacyAddStartnodeFragments mover");
	}

	core::pose::PoseOP native_pose( new pose::Pose );
	core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);

	core::Size start_window = protocols::rosetta_scripts::find_nearest_res( pose, *native_pose, start_res_);
	core::Size end_window = protocols::rosetta_scripts::find_nearest_res( pose, *native_pose, end_res_);

	std::string frag9_file = option[ in::file::frag9 ]();
	core::fragment::FragmentIO frag_io(
		//option[ OptionKeys::abinitio::number_9mer_frags ](),
		200,
		option[ OptionKeys::frags::nr_large_copies ](),
		option[ OptionKeys::frags::annotate ]());
	core::fragment::FragSetOP fragset = frag_io.read_data( frag9_file );
	TR << "Stealing native fragments for windows " << start_window << " to " << end_window << std::endl;
	frag_io.write_data(frag9_file+".nonat", *fragset);
	core::fragment::steal_frag_set_from_pose( *native_pose, start_window, end_window, *fragset, core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset->max_frag_length() ) ) );
	frag_io.write_data(frag9_file+".nat", *fragset);
}

void
LegacyAddStartnodeFragments::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	using namespace basic::options;

	if ( !tag->hasOption("start_res") || !tag->hasOption("end_res") ) {
		utility_exit_with_message("Must specify start_res and end_res");
	}

	start_res_ = tag->getOption<core::Size>("start_res");
	end_res_ = tag->getOption<core::Size>("end_res");

}

std::string LegacyAddStartnodeFragments::get_name() const {
	return mover_name();
}

std::string LegacyAddStartnodeFragments::mover_name() {
	return "LegacyAddStartnodeFragments";
}

void LegacyAddStartnodeFragments::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "start_res", xsct_positive_integer, "Beginning of window for fragment stealing" )
		+ XMLSchemaAttribute::required_attribute( "end_res", xsct_positive_integer, "End of window for fragment stealing" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Steals native fragment for residues between the provided start_res and end_res and adds it to the provided frag9 file", attlist );
}

std::string LegacyAddStartnodeFragmentsCreator::keyname() const {
	return LegacyAddStartnodeFragments::mover_name();
}

protocols::moves::MoverOP
LegacyAddStartnodeFragmentsCreator::create_mover() const {
	return protocols::moves::MoverOP( new LegacyAddStartnodeFragments );
}

void LegacyAddStartnodeFragmentsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LegacyAddStartnodeFragments::provide_xml_schema( xsd );
}



} //legacy_sewing
} //protocols
