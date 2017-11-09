// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PeptideStapleDesignMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PeptideStapleDesignMover.hh>
#include <protocols/protein_interface_design/movers/PeptideStapleDesignMoverCreator.hh>

// Package headers

// Project headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <protocols/simple_moves/PeptideStapleMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.PeptideStapleDesignMover" );

PeptideStapleDesignMover::PeptideStapleDesignMover() :
	protocols::moves::Mover( PeptideStapleDesignMover::mover_name() )
{}

PeptideStapleDesignMover::PeptideStapleDesignMover( std::string const & seqpos, core::Size const staple_gap ) :
	protocols::moves::Mover( PeptideStapleDesignMover::mover_name() ),
	staple_start_( seqpos ),
	staple_gap_( staple_gap )
{}

PeptideStapleDesignMover::~PeptideStapleDesignMover() = default;

protocols::moves::MoverOP
PeptideStapleDesignMover::clone() const {
	return( protocols::moves::MoverOP( new PeptideStapleDesignMover( *this ) ));
}

void PeptideStapleDesignMover::apply( core::pose::Pose & pose )
{
	core::Size staple_start( core::pose::parse_resnum( staple_start_, pose ) );

	protocols::simple_moves::PeptideStapleMover stapler( staple_start, staple_gap_ );
	stapler.apply( pose );
}

void
PeptideStapleDesignMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	staple_start_ = core::pose::get_resnum_string( tag );
	staple_gap_ = tag->getOption<core::Size>( "staple_gap", 4 );
}

std::string PeptideStapleDesignMover::get_name() const {
	return mover_name();
}

std::string PeptideStapleDesignMover::mover_name() {
	return "StapleMover";
}

void PeptideStapleDesignMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "resnum", xsct_refpose_enabled_residue_number, "Residue number where the staple starts, in Rosetta or PDB or reference pose numbering" )
		+ XMLSchemaAttribute( "pdb_num", xsct_refpose_enabled_residue_number, "Residue number where the staple starts, in Rosetta or PDB or reference pose numbering" )
		+ XMLSchemaAttribute::attribute_w_default( "staple_gap", xsct_non_negative_integer, "Gap from starting to ending stapled residue", "4" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Design a stapled peptide at a protein interface", attlist );
}

std::string PeptideStapleDesignMoverCreator::keyname() const {
	return PeptideStapleDesignMover::mover_name();
}

protocols::moves::MoverOP
PeptideStapleDesignMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PeptideStapleDesignMover );
}

void PeptideStapleDesignMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideStapleDesignMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols

