// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeMinMover.cc
/// @brief A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/carbohydrates/GlycanTreeMinMover.hh>
#include <protocols/carbohydrates/GlycanTreeMinMoverCreator.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/kinematics/util.hh>
#include <core/select/residue_selector/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <numeric/random/random.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.carbohydrates.GlycanTreeMinMover" );

namespace protocols {
namespace carbohydrates {

using namespace core::kinematics;
using namespace core::scoring;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycanTreeMinMover::GlycanTreeMinMover():
	protocols::simple_moves::MinMover()
{

}

GlycanTreeMinMover::GlycanTreeMinMover(
	core::kinematics::MoveMapOP movemap_in,
	ScoreFunctionCOP scorefxn_in,
	std::string const & min_type_in /*dfpmin_armijo_nonmonotone */,
	Real tolerance_in /* .01 */):

	protocols::simple_moves::MinMover(movemap_in, scorefxn_in, min_type_in, tolerance_in, false /*use nblist*/)

{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycanTreeMinMover::GlycanTreeMinMover( GlycanTreeMinMover const & src ):
	protocols::simple_moves::MinMover( src ),
	mm_residues_(src.mm_residues_)
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycanTreeMinMover::~GlycanTreeMinMover(){}


////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
GlycanTreeMinMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new GlycanTreeMinMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanTreeMinMover::clone() const
{
	return protocols::moves::MoverOP( new GlycanTreeMinMover( *this ) );
}

/// @brief Get the name of the Mover
std::string
GlycanTreeMinMover::get_name() const
{
	return GlycanTreeMinMover::class_name();
}

std::string
GlycanTreeMinMover::class_name()
{
	return "GlycanTreeMinMover";
}


/////////////// Creator ///////////////

protocols::moves::MoverOP
GlycanTreeMinMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GlycanTreeMinMover );
}

std::string
GlycanTreeMinMoverCreator::keyname() const
{
	return GlycanTreeMinMover::class_name();
}

void
GlycanTreeMinMover::set_movemap(core::kinematics::MoveMapCOP movemap_in){
	runtime_assert( movemap_in != nullptr );
	protocols::simple_moves::MinMover::set_movemap(core::kinematics::MoveMapOP( new MoveMap( *movemap_in ) ));
	mm_residues_.clear();

}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////


/// @brief Apply the mover
void
GlycanTreeMinMover::apply( core::pose::Pose& pose){

	core::select::residue_selector::GlycanResidueSelector selector = core::select::residue_selector::GlycanResidueSelector();
	selector.set_include_root( true ); //Since we wont really have the ASN root, and we want the possibility to sample the whole glycan.

	utility::vector1< bool > root(pose.total_residue(), false);

	MoveMapOP active_movemap = MoveMapOP( new MoveMap());

	//Set default movemap, only include carbohydrate residues.
	if ( ! movemap() ) {
		MoveMapOP full_mm = MoveMapOP( new MoveMap ) ;
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( !pose.residue(i).is_carbohydrate() ) {
				full_mm->set_bb(false, i);
				full_mm->set_chi(false, i);
			} else {
				full_mm->set_bb(true, i);
				full_mm->set_chi(true, i);
			}
			set_movemap(full_mm);
		}
	}
	if ( mm_residues_.size() == 0 ) {
		mm_residues_ = core::kinematics::get_residues_from_movemap_bb_or_chi(*movemap(), pose.total_residue());
	}


	//Randomly select a residue from the movemap.
	core::Size index = numeric::random::rg().random_range( 1, mm_residues_.size() );
	core::Size resnum = mm_residues_[ index ];

	TR << "Minimizing from carbohydrate root: " << resnum << std::endl;
	core::Size min_residue_n = 0;

	//Get downstream Tree.
	selector.set_select_from_branch_residue( resnum );

	utility::vector1< bool > subset = selector.apply(pose);

	//Setup the active movemap
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( subset[ i ] ) {
			min_residue_n+=1;
			if ( movemap()->get_bb(i) ) {
				active_movemap->set_bb(i, true);
			}
			if ( movemap()->get_chi(i) ) {
				active_movemap->set_chi(i, true);
			}
			//This needs to be total number of possible glycan torsions!
			for ( core::Size xx = 1; xx <= 10; ++xx ) {
				if ( movemap()->get_bb( i, xx) ) {
					active_movemap->set_bb(i, xx, true);
				}
			}
		} else {
			active_movemap->set_bb( i, false );
			active_movemap->set_chi( i, false );

		}
	}

	//Minimize
	TR << "Minimizing " << min_residue_n << " residues " << std::endl;
	apply_dof_tasks_to_movemap(pose, *active_movemap);

	if ( ! score_function() ) score_function(get_score_function()); // get a default (INITIALIZED!) ScoreFunction

	minimize( pose, *active_movemap );


}

void GlycanTreeMinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = MinMover::complex_type_generator_for_min_mover( xsd );
	ct_gen->element_name( class_name() )
		.description( "A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector" )
		.write_complex_type_to_schema( xsd );
}

void GlycanTreeMinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanTreeMinMover::provide_xml_schema( xsd );
}




} //protocols
} //carbohydrates

