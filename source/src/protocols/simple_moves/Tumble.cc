// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/simple_moves/Tumble.hh>
#include <protocols/simple_moves/TumbleCreator.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/select/util.hh>
#include <core/select/residue_selector/ResidueSpanSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/pose/selection.hh>
#include <core/pose/ResidueIndexDescription.hh>
#include <core/kinematics/FoldTree.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.Tumble" );

namespace protocols {
namespace simple_moves {

Tumble::Tumble() = default;

Tumble::~Tumble() = default;

numeric::xyzVector<core::Real>
Tumble::center_of_mass(core::pose::Pose const & pose) {
	int nAtms = 0;
	numeric::xyzVector<core::Real> massSum(0.,0.,0.), CoM;

	utility::vector1< core::Size > residue_list;
	if ( residue_list_ != nullptr ) {
		residue_list = core::select::get_residues_from_subset( residue_list_->apply( pose ) );
	}

	for ( core::Size ires: residue_list ) {
		if ( pose.residue_type(ires).aa() == core::chemical::aa_vrt ) continue;

		for ( core::Size iatom = 1; iatom <= pose.residue_type(ires).nheavyatoms(); ++iatom ) {
			core::conformation::Atom const & atom( pose.residue(ires).atom(iatom) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	CoM = massSum / (core::Real)nAtms;
	return CoM;
}

void Tumble::apply( core::pose::Pose & pose ) {
	numeric::xyzVector<core::Real> CoM;
	CoM = center_of_mass(pose);

	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > coords;

	numeric::xyzVector< core::Real > axis1 = numeric::xyzVector< core::Real >( 0.0, 0.0, 1.0 );
	core::Real angle1 = 360. * numeric::random::rg().uniform();
	numeric::xyzMatrix< core::Real > rotation_matrix1( numeric::rotation_matrix_degrees(axis1, angle1 ) );

	numeric::xyzVector< core::Real > axis2 = numeric::xyzVector< core::Real >( numeric::random::rg().uniform(), numeric::random::rg().uniform(), 0.0 );
	core::Real angle2 = 180. * numeric::random::rg().uniform();
	numeric::xyzMatrix< core::Real > rotation_matrix2( numeric::rotation_matrix_degrees(axis2, angle2 ) );

	utility::vector1< core::Size > residue_list;
	if ( residue_list_ != nullptr ) {
		residue_list = core::select::get_residues_from_subset( residue_list_->apply( pose ) );
	}

	for ( core::Size ires: residue_list ) {
		for ( core::Size iatom = 1; iatom <= pose.residue_type(ires).natoms(); ++iatom ) {

			numeric::xyzVector< core::Real > translated2origin = pose.residue(ires).atom(iatom).xyz() - CoM;
			numeric::xyzVector< core::Real > rotated_at_origin = rotation_matrix2 * rotation_matrix1 * translated2origin;
			numeric::xyzVector< core::Real > new_coord = rotated_at_origin + CoM;

			ids.push_back(core::id::AtomID(iatom,ires));
			coords.push_back(new_coord);
		}
	}

	pose.batch_set_xyz(ids, coords);
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Tumble::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & ,
	Filters_map const & ,
	moves::Movers_map const & ,
	Pose const &
) {
	if ( tag->hasOption( "chain_number" ) ) {
		auto chain = tag->getOption< std::string >( "chain_number" );
		residue_list_ = core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::ChainSelector( chain ) );
	} else {
		core::pose::ResidueIndexDescriptionCOP start_res;
		core::pose::ResidueIndexDescriptionCOP stop_res;

		if ( tag->hasOption( "start_res" ) ) {
			start_res = core::pose::parse_resnum( tag->getOption< std::string >( "start_res" ) );
		} else {
			start_res = core::pose::make_rid_posenum( 1 );
		}

		if ( tag->hasOption( "stop_res" ) ) {
			stop_res = core::pose::parse_resnum( tag->getOption< std::string >( "stop_res" ) );
		} else {
			stop_res = core::pose::ResidueIndexDescriptionCOP( new core::pose::ResidueIndexDescriptionLastResidue );
		}

		residue_list_ = core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::ResidueSpanSelector( stop_res, start_res ) );
	}
}

moves::MoverOP Tumble::clone() const {
	return moves::MoverOP( new Tumble( *this ) );
}
moves::MoverOP Tumble::fresh_instance() const {
	return moves::MoverOP( new Tumble );
}

std::string Tumble::get_name() const {
	return mover_name();
}

std::string Tumble::mover_name() {
	return "Tumble";
}

void Tumble::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"chain_number", xsct_non_negative_integer,
		"The chain to tumble - not compatible with start_res/stop_res.");

	core::pose::attributes_for_parse_resnum( attlist, "start_res", "The first residue of the residue span to tumble - not compatible with chain_number." );
	core::pose::attributes_for_parse_resnum( attlist, "stop_res", "The last residue of the residue span to tumble - not compatible with chain_number." );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Tumble (randomly rotate around the center of mass) a given subsection of the protein.",
		attlist );
}

std::string TumbleCreator::keyname() const {
	return Tumble::mover_name();
}

protocols::moves::MoverOP
TumbleCreator::create_mover() const {
	return protocols::moves::MoverOP( new Tumble );
}

void TumbleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Tumble::provide_xml_schema( xsd );
}


} // moves
} // protocols
