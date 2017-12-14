// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/PertMinMover.cc
/// @brief Roberto Chica inspired random perturbation followed by minimization
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/simple_moves/PertMinMover.hh>
#include <protocols/simple_moves/PertMinMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/xyzVector.io.hh>

static basic::Tracer TR( "protocols.simple_moves.PertMinMover" );


namespace protocols {
namespace simple_moves {

PertMinMover::PertMinMover():
	protocols::moves::Mover( mover_name() ),
	pert_size_(0.001),
	uniform_(true),
	sc_only_(false)
{}

PertMinMover::~PertMinMover()= default;

PertMinMover::PertMinMover( PertMinMover const & ) = default;

std::string
PertMinMover::mover_name() {
	return "PertMinMover";
}

void PertMinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "pert_size", xsct_real, "How large a perturbation to make", "0.001" )
		+ XMLSchemaAttribute::attribute_w_default( "uniform", xsct_rosetta_bool, "Do the perturbation in a uniform sphere, rather than Gaussian", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "sc_only", xsct_rosetta_bool, "Only do the perturbations on sidechain atoms (not backbone)", "false" );

	protocols::rosetta_scripts::attributes_for_parse_residue_selector( attlist );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	utility::tag::XMLSchemaSimpleSubelementList subelements;
	protocols::rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Do a Roberto Chica inspired random atom coordinate perturbation followed by (Cartesian) minimization.",
		attlist, subelements );
}

void
PertMinMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	pert_size( tag->getOption< core::Real >( "pert_size", pert_size() ) );
	uniform( tag->getOption< bool >( "uniform", uniform() ) );
	sc_only( tag->getOption< bool >( "sc_only", sc_only() ) );
	residues( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );

	movemap_factory( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data ) );
}

protocols::moves::MoverOP
PertMinMover::clone() const{
	return protocols::moves::MoverOP( new PertMinMover( *this ) );
}

/*
PertMinMover & PertMinMoveroperator=( PertMinMover const & src){
return PertMinMover( src );
}
*/


moves::MoverOP
PertMinMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PertMinMover );
}

std::string
PertMinMover::get_name() const {
	return "PertMinMover";
}

void
PertMinMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, PertMinMover const &mover)
{
	mover.show(os);
	return os;
}


void
PertMinMover::apply( core::pose::Pose& pose){
	runtime_assert( scorefxn_ );

	// Determine which residues to use: if not specified, use all
	utility::vector1< bool > resi( pose.total_residue(), true );
	if ( residues_ ) {
		TR << "Applying ResidueSelector." << std::endl;
		resi = residues_->apply( pose );
	}

	pert(pose, resi);
	min(pose);
}

void
PertMinMover::pert( core::pose::Pose & pose, utility::vector1< bool > const & resi ) const {
	// Pert
	for ( core::Size ii(1); ii <= pose.total_residue(); ++ii ) {
		if ( ! resi[ii] ) { continue; }
		core::conformation::Residue const & res( pose.residue(ii) );

		for ( core::Size jj(1); jj <= res.natoms(); ++jj ) {
			if ( sc_only_ && res.atom_is_backbone( jj ) ) { continue; }

			// Find a random perturbation
			core::Vector pert;
			if ( uniform_ ) {
				pert = numeric::random::uniform_vector_sphere();
			} else {
				pert = numeric::random::random_vector_spherical();
			}
			pert *= pert_size_;

			// Move the atom
			core::id::AtomID atomid( jj, ii ); // atom, then res
			pose.set_xyz( atomid, pose.xyz( atomid ) + pert );
		}
	}
}

void
PertMinMover::min( core::pose::Pose & pose ) const {
	core::kinematics::MoveMapOP mm;
	if ( movemap_factory_ ) {
		mm = movemap_factory_->create_movemap_from_pose( pose );
	} else {
		TR << "Using a default MoveMap." << std::endl;
		mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}

	core::scoring::ScoreFunctionCOP sfxn( scorefxn() );
	if ( ! sfxn ) {
		TR << "Using the default ScoreFunction." << std::endl;
		sfxn = core::scoring::get_score_function();
	} // get a default (INITIALIZED!) ScoreFunction

	TR.Debug << "Score before: " << (*sfxn)(pose) << std::endl;

	TR.Debug << "MoveMap " << *mm << std::endl;

	// We *need* to use Cartesian minimization here.
	// The Cartesian pertubation will have messed up the bond lengths and angles
	core::optimization::CartesianMinimizer minimizer;
	(*sfxn)(pose);
	// Could this be default (or exposed?)
	core::optimization::MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	minimizer.run( pose, *mm, *sfxn, min_options );

	TR.Debug << "Score after: " << (*sfxn)(pose) << std::endl;
}


/////////////// Creator ///////////////

protocols::moves::MoverOP
PertMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PertMinMover );
}

std::string
PertMinMoverCreator::keyname() const {
	return PertMinMover::mover_name();
}

void PertMinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PertMinMover::provide_xml_schema( xsd );
}

} //protocols
} //simple_moves


