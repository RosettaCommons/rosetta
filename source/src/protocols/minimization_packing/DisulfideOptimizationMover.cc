// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/minimization_packing/DisulfideOptimizationMover.cc
/// @brief A Mover to jointly optimize the geometry of a pair of disulfide-bonded residues.
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/minimization_packing/DisulfideOptimizationMover.hh>
#include <protocols/minimization_packing/DisulfideOptimizationMoverCreator.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR( "protocols.minimization_packing.DisulfideOptimizationMover" );

namespace protocols {
namespace minimization_packing {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
DisulfideOptimizationMover::DisulfideOptimizationMover():
	protocols::moves::Mover( DisulfideOptimizationMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
DisulfideOptimizationMover::DisulfideOptimizationMover( DisulfideOptimizationMover const & src ):
	protocols::moves::Mover( src ),
	selector_( src.selector_->clone() ),
	sfxn_( src.sfxn_->clone() ),
	final_optimization_n_iter_( src.final_optimization_n_iter_ )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
DisulfideOptimizationMover::~DisulfideOptimizationMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////


void
DisulfideOptimizationMover::break_repack_reform( Pose & pose, utility::vector1< Size > const & cys_pos ) {
	// Mutate both residues to their corresponding thiol type.
	// In the future here we could be careful to remove the HG,
	// but they can easily be optimized away.
	core::conformation::break_disulfide( pose.conformation(), cys_pos[ 1 ], cys_pos[ 2 ] );

	// Repack.
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();
	task->restrict_to_repacking();
	core::pack::pack_rotamers( pose, *sfxn_, task );

	// Mutate back to a disulfide. The bond length will likely not be ideal.
	core::conformation::form_disulfide( pose.conformation(), cys_pos[ 1 ], cys_pos[ 2 ], false, false );
}

/// @brief Apply the mover
void
DisulfideOptimizationMover::apply( core::pose::Pose & pose ) {

	if ( !sfxn_ ) {
		sfxn_ = core::scoring::get_score_function();
	}

	if ( !selector_ ) {
		TR.Error << "Oddly, there is no residue selector defined. This won't do anything, and that is almost certainly not intended." << std::endl;
		utility_exit();
	}

	// Take the disulfide indicated by the residue selector
	auto disulfide_cys = selector_->apply( pose );
	utility::vector1< Size > cys_pos;
	for ( Size ii = 1; ii <= disulfide_cys.size(); ++ii ) {
		if ( disulfide_cys[ ii ] ) cys_pos.push_back( ii );
	}

	// Confirm that both are cysteine-type residues in a disulfide bond to each other.
	if ( !pose.residue_type( cys_pos[ 1 ] ).is_disulfide_bonded() || !pose.residue_type( cys_pos[ 1 ] ).is_disulfide_bonded() ) {
		TR.Error << "Residue " << cys_pos[ 1 ] << " name is " << pose.residue_type( cys_pos[ 1 ] ).name() << "." << std::endl;
		TR.Error << "Residue " << cys_pos[ 2 ] << " name is " << pose.residue_type( cys_pos[ 2 ] ).name() << "." << std::endl;
		TR.Error << "Residues " << cys_pos[ 1 ] << " and " << cys_pos[ 2 ] << " are not disulfide bonded at all.";
		utility_exit();
	}

	if ( core::conformation::get_disulf_partner( pose.conformation(), cys_pos[ 1 ] ) != cys_pos[ 2 ] ||
			core::conformation::get_disulf_partner( pose.conformation(), cys_pos[ 2 ] ) != cys_pos[ 1 ] ) {
		TR.Error << "Residue " << cys_pos[ 1 ] << " name is " << pose.residue_type( cys_pos[ 1 ] ).name() << "." << std::endl;
		TR.Error << "Residue " << cys_pos[ 2 ] << " name is " << pose.residue_type( cys_pos[ 2 ] ).name() << "." << std::endl;
		TR.Error << "Residues " << cys_pos[ 1 ] << " and " << cys_pos[ 2 ] << " are disulfide bonded, but are not disulfide bonded to each other.";
		utility_exit();
	}


	break_repack_reform( pose, cys_pos );

	if ( final_optimization_n_iter_ > 0 ) {
		// Minimize with that bond length turned on (on the off chance that it's
		// in the atom tree at all)
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		//movemap->set_bb( cys_pos[ 1 ], true );
		//movemap->set_bb( cys_pos[ 2 ], true );
		//movemap->set_chi( cys_pos[ 1 ], true );
		//movemap->set_chi( cys_pos[ 2 ], true );
		movemap->set_chi( true );
		//movemap.set() it's almost certainly not there, don't bother.

		protocols::minimization_packing::MinMover minmover( movemap, sfxn_, "lbfgs_armijo_nonmonotone", 0.0001, true );

		minmover.apply( pose );
		for ( Size ii = 1; ii <= final_optimization_n_iter_; ++ii ) {
			break_repack_reform( pose, cys_pos );
			minmover.apply( pose );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
DisulfideOptimizationMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
DisulfideOptimizationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	core::select::residue_selector::ResidueSelectorCOP selector = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	if ( selector ) {
		TR << "Setting selector " << selector->get_name() << std::endl;
		set_selector( selector );
	}

	if ( tag->hasOption("scorefxn") ) {
		set_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	} else {
		set_score_function( core::scoring::get_score_function() );
	}

	if ( tag->hasOption("final_optimization_n_iter") ) {
		set_final_optimization_n_iter( tag->getOption< Size >( "final_optimization_n_iter", 0 ) );
	}
}

void DisulfideOptimizationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "final_optimization", xsct_non_negative_integer,
		"Number of iterations of final optimization (minimization and three more rounds of repack-minimize)", "0" );


	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Optimize the pair of disulfide-bonded residues indicated in the residue selector provided.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
DisulfideOptimizationMover::fresh_instance() const
{
	return utility::pointer::make_shared< DisulfideOptimizationMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
DisulfideOptimizationMover::clone() const
{
	return utility::pointer::make_shared< DisulfideOptimizationMover >( *this );
}

std::string DisulfideOptimizationMover::get_name() const {
	return mover_name();
}

std::string DisulfideOptimizationMover::mover_name() {
	return "DisulfideOptimizationMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
DisulfideOptimizationMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< DisulfideOptimizationMover >();
}

std::string
DisulfideOptimizationMoverCreator::keyname() const
{
	return DisulfideOptimizationMover::mover_name();
}

void DisulfideOptimizationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisulfideOptimizationMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, DisulfideOptimizationMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //minimization_packing
