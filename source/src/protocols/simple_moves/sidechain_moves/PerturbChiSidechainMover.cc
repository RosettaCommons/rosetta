// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMoverBase.cc
/// @brief implementation of PerturbChiSidechainMover class and functions
/// @author Oliver Lange ( oliver.lange@tum.de )


#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/PerturbChiSidechainMoverCreator.hh>

// Procols Headers
#include <basic/datacache/DataMap.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/prof.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Utility
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// C++ Headers
#include <sstream>
#include <fstream>
#include <utility/fixedsizearray1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
using namespace core;
using namespace core::pose;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.sidechain_moves.PerturbChiSidechainMover" );

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace chemical;
using namespace conformation;


// XRW TEMP std::string
// XRW TEMP PerturbChiSidechainMoverCreator::keyname() const {
// XRW TEMP  return PerturbChiSidechainMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP PerturbChiSidechainMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PerturbChiSidechainMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PerturbChiSidechainMover::mover_name() {
// XRW TEMP  return "PerturbChiSidechain";
// XRW TEMP }


PerturbChiSidechainMover::PerturbChiSidechainMover() {
	protocols::moves::Mover::type( "PerturbChiSidechain" );
	set_defaults();
}

PerturbChiSidechainMover::PerturbChiSidechainMover(
	pack::dunbrack::RotamerLibrary const & rotamer_library
) : Parent( rotamer_library ) {
	set_defaults();
}

PerturbChiSidechainMover::PerturbChiSidechainMover(
	PerturbChiSidechainMover const & mover
) : Parent ( mover ) {
	set_defaults();
}

protocols::moves::MoverOP
PerturbChiSidechainMover::clone() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::sidechain_moves::PerturbChiSidechainMover(*this) );
}

void
PerturbChiSidechainMover::set_defaults() {
	magnitude_ = 10;
	gaussian_ = false;
}

void
PerturbChiSidechainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	pose::Pose const & /*pose*/
) {
	magnitude_ = tag->getOption<Real>( "magnitude", magnitude_ );
	gaussian_ = tag->getOption<bool>( "gaussian", gaussian_ );
}

// XRW TEMP std::string
// XRW TEMP PerturbChiSidechainMover::get_name() const {
// XRW TEMP  return "PerturbChiSidechainMover";
// XRW TEMP }

void
PerturbChiSidechainMover::make_chi_move(
	conformation::Residue const&,
	ChiVector const& old_chi,
	ChiVector& new_chi
) {
	new_chi.resize( old_chi.size() );
	for ( Size i = 1; i <= old_chi.size(); i++ ) {
		if ( !gaussian_ ) {
			Real rand = numeric::random::rg().uniform();
			new_chi[ i ] = basic::periodic_range( (( 2.0*rand-1.0 )*magnitude_ + old_chi[ i ]) , 360.0 );
		} else {
			new_chi[ i ] = basic::periodic_range( old_chi[i] + numeric::random::rg().gaussian()*magnitude_, 360.0 );
		}
	}
}

std::string PerturbChiSidechainMover::get_name() const {
	return mover_name();
}

std::string PerturbChiSidechainMover::mover_name() {
	return "PerturbChiSidechain";
}

void PerturbChiSidechainMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	// AMW: "effective default" issue
	attlist + XMLSchemaAttribute::attribute_w_default("magnitude", xsct_real, "Size of chi perturbations.", "10" )
		+ XMLSchemaAttribute::attribute_w_default("gaussian", xsct_rosetta_bool, "Sample from a normal distribution (true) or uniformly (false).", "false" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sample side chain chi angles.", attlist );
}

std::string PerturbChiSidechainMoverCreator::keyname() const {
	return PerturbChiSidechainMover::mover_name();
}

protocols::moves::MoverOP
PerturbChiSidechainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PerturbChiSidechainMover );
}

void PerturbChiSidechainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PerturbChiSidechainMover::provide_xml_schema( xsd );
}


} // sidechain_moves
} // simple_moves
} // protocols
