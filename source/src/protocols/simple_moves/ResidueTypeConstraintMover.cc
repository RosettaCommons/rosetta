// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueTypeConstraintMover.cc
/// @brief Assigns a ResidueTypeConstraint to a pose.
/// @author Doo Nam Kim (doonam.kim@gmail.com) (All I did is just making a simple mover using Sarel's ResidueTypeConstraint)

#include <protocols/simple_moves/ResidueTypeConstraintMover.hh>
#include <protocols/simple_moves/ResidueTypeConstraintMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.ResidueTypeConstraintMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace scoring;
using namespace constraints;

using namespace utility::tag;

// XRW TEMP std::string
// XRW TEMP ResidueTypeConstraintMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ResidueTypeConstraintMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ResidueTypeConstraintMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ResidueTypeConstraintMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ResidueTypeConstraintMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ResidueTypeConstraintMover";
// XRW TEMP }

// constructors
ResidueTypeConstraintMover::ResidueTypeConstraintMover()
: protocols::moves::Mover( ResidueTypeConstraintMover::mover_name() )
{
}

ResidueTypeConstraintMover::ResidueTypeConstraintMover( std::string const & type )
: protocols::moves::Mover(type)
{
}

// destructor
ResidueTypeConstraintMover::~ResidueTypeConstraintMover()= default;


void
ResidueTypeConstraintMover::apply( Pose & pose )
{
	int AA_name3_length_1 = std::count(AA_name3_.begin(), AA_name3_.end(), ',');
	std::string delimiter = ",";
	for ( Size resnum=1; resnum<=pose.size(); ++resnum ) {
		int substr_index = 0;
		for ( int i =0; i<=(AA_name3_length_1); i++ ) {
			size_t pos = AA_name3_.find(delimiter);

			std::string favored_res = AA_name3_.substr(substr_index, pos);
			substr_index = substr_index + pos + 1;
			// essentially just same as "substr_index = substr_index + 4";

			ResidueTypeConstraintOP res_favor( new core::scoring::constraints::ResidueTypeConstraint(pose, resnum, favored_res, favor_bonus_) );

			pose.add_constraint(res_favor);
		}
	}
}

// XRW TEMP std::string
// XRW TEMP ResidueTypeConstraintMover::get_name() const {
// XRW TEMP  return ResidueTypeConstraintMover::mover_name();
// XRW TEMP }

protocols::moves::MoverOP ResidueTypeConstraintMover::clone() const { return protocols::moves::MoverOP( new protocols::simple_moves::ResidueTypeConstraintMover( *this ) ); }
protocols::moves::MoverOP ResidueTypeConstraintMover::fresh_instance() const { return protocols::moves::MoverOP( new ResidueTypeConstraintMover ); }

void
ResidueTypeConstraintMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	if ( tag->hasOption("AA_name3") ) {
		AA_name3_ = tag->getOption<std::string>("AA_name3"); // for example: ASP,GLU
	} else {
		TR << "please specifiy AA_name3 like SER,THR" << std::endl;
	}
	favor_bonus_ = tag->getOption<Real>("favor_bonus", 0.5);
	// since this mover does not necessarily favor native sequence only, I named it "favor_bonus" instead of "native_bonus" to make it more general
	// positively higher bonus gives more favorable selection to (a) specified residue(s)

}

std::string ResidueTypeConstraintMover::get_name() const {
	return mover_name();
}

std::string ResidueTypeConstraintMover::mover_name() {
	return "ResidueTypeConstraintMover";
}

void ResidueTypeConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "AA_name3", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "favor_bonus", xsct_real, "XRW TO DO", "0.5" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string ResidueTypeConstraintMoverCreator::keyname() const {
	return ResidueTypeConstraintMover::mover_name();
}

protocols::moves::MoverOP
ResidueTypeConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ResidueTypeConstraintMover );
}

void ResidueTypeConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueTypeConstraintMover::provide_xml_schema( xsd );
}

} // moves
} // protocols
