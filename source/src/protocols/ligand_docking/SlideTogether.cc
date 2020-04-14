// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/SlideTogether.hh>
#include <protocols/ligand_docking/SlideTogetherCreator.hh>

#include <protocols/docking/DockingInitialPerturbation.hh>

// Utility Headers

#include <utility>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace ligand_docking {

static basic::Tracer slide_together_tracer( "protocols.ligand_docking.ligand_options.slide_together" );




SlideTogether::SlideTogether()= default;

SlideTogether::SlideTogether(std::string const & chain):
	chains_( 1, chain ) // vector of size 1, containing one chain
{}

protocols::moves::MoverOP SlideTogether::clone() const {
	return utility::pointer::make_shared< SlideTogether >( *this );
}

protocols::moves::MoverOP SlideTogether::fresh_instance() const {
	return utility::pointer::make_shared< SlideTogether >();
}


//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
SlideTogether::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/
)
{
	if ( tag->getName() != "SlideTogether" ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	if ( ! tag->hasOption("chains") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'SlideTogether' mover requires chains tag");

	std::string const chains_str = tag->getOption<std::string>("chains");
	chains_ = utility::string_split(chains_str, ',');
}

void
SlideTogether::apply( core::pose::Pose & pose ){
	slide_together_tracer<< "Applying slide_together"<< std::endl;

	utility::vector1< core::Size > jumps;
	for ( std::string const & chain_str : chains_ ) {
		// This two-step is historical, we probably should be able to short-circuit it.
		utility::vector1<core::Size> chain_ids= get_chain_ids_from_chain(chain_str, pose);
		for ( core::Size const chain_id : chain_ids ) {
			core::Size jump_id= core::pose::get_jump_id_from_chain_id(chain_id, pose);
			jumps.push_back(jump_id);
		}
	}

	protocols::docking::FaDockingSlideIntoContact slideTogether(jumps);
	slideTogether.apply(pose);
}

std::string SlideTogether::get_name() const {
	return mover_name();
}

std::string SlideTogether::mover_name() {
	return "SlideTogether";
}

void SlideTogether::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("chains", xs_string, "Comma separated list of chain IDs.");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Move the small molecule and protein into close proximity.", attlist );
}

std::string SlideTogetherCreator::keyname() const {
	return SlideTogether::mover_name();
}

protocols::moves::MoverOP
SlideTogetherCreator::create_mover() const {
	return utility::pointer::make_shared< SlideTogether >();
}

void SlideTogetherCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SlideTogether::provide_xml_schema( xsd );
}


} //namespace ligand_docking
} //namespace protocols
