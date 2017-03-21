// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SimpleHbondsToAtomFilter.cc
/// @brief Simple filter for detercting Hbonds to atom with energy < energy cutoff
/// @author Benjamin Basanta (basantab@uw.edu)

#include <protocols/simple_filters/SimpleHbondsToAtomFilter.hh>
#include <protocols/simple_filters/SimpleHbondsToAtomFilterCreator.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

//Project headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <core/kinematics/tree/Atom_.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.SimpleHbondsToAtomFilter" );

namespace protocols {
namespace simple_filters {

SimpleHbondsToAtomFilter::SimpleHbondsToAtomFilter():
	protocols::filters::Filter( "SimpleHbondsToAtomFilter" ),
	target_atom_name_( /* NULL */ ),
	target_residue_( /* NULL */ ),
	n_partners_( 1 ),
	scorefxn_( /* NULL */ ),
	hb_e_cutoff_( -0.5 )
{}

SimpleHbondsToAtomFilter::~SimpleHbondsToAtomFilter()
{}

// Setters and getters

core::Size
SimpleHbondsToAtomFilter::get_n_partners() const{
	return n_partners_;
}

void
SimpleHbondsToAtomFilter::set_n_partners( core::Size n_partners ){
	n_partners_ = n_partners;
}

core::Size
SimpleHbondsToAtomFilter::get_target_residue() const{
	return target_residue_;
}

void
SimpleHbondsToAtomFilter::set_target_residue( core::Size target_residue ){
	target_residue_ = target_residue;
}

std::string
SimpleHbondsToAtomFilter::get_target_atom_name() const{
	return target_atom_name_;
}

void
SimpleHbondsToAtomFilter::set_target_atom_name( std::string target_atom_name ){
	target_atom_name_ = target_atom_name;
}

core::scoring::ScoreFunctionOP
SimpleHbondsToAtomFilter::get_scorefxn() const{
	return scorefxn_;
}

void
SimpleHbondsToAtomFilter::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::Real
SimpleHbondsToAtomFilter::get_hb_e_cutoff() const{
	return hb_e_cutoff_;
}

void
SimpleHbondsToAtomFilter::set_hb_e_cutoff( core::Real hb_e_cutoff ){
	hb_e_cutoff_ = hb_e_cutoff;
}

// Setters and getters END

void
SimpleHbondsToAtomFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & pose )
{
	TR << "SimpleHbondsToAtomFilter"<<std::endl;
	set_target_residue(core::pose::get_resnum( tag, pose ));
	set_n_partners( tag->getOption< core::Size >( "n_partners", 1 ) );
	set_target_atom_name( tag->getOption< std::string >( "target_atom_name", "" ) );
	set_hb_e_cutoff(tag->getOption< core::Real >( "hb_e_cutoff", -0.5) );
	set_scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	TR << "Will look for "<<get_n_partners()<<" Hbonds to atom "<<get_target_atom_name()<<" in residue "<<get_target_residue()<<"."<<std::endl;
}

protocols::filters::FilterOP
SimpleHbondsToAtomFilter::clone() const
{
	return protocols::filters::FilterOP( new SimpleHbondsToAtomFilter( *this ) );
}


protocols::filters::FilterOP
SimpleHbondsToAtomFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new SimpleHbondsToAtomFilter );
}

bool
SimpleHbondsToAtomFilter::apply( core::pose::Pose const & pose ) const
{
	if ( compute( pose ) >= get_n_partners() ) {
		return true;
	} else {
		return false;
	}
}

core::Real
SimpleHbondsToAtomFilter::report_sm( core::pose::Pose const & pose ) const
{
	return core::Real( compute( pose ) );
}

/////////////* COMPUTE *///////////////
core::Size
SimpleHbondsToAtomFilter::compute( core::pose::Pose const & pose) const
{
	// I will count Hbonds with this counter:
	core::Size counter = 0;

	// Get cp of the pose so you can score it:
	core::pose::Pose p = pose;

	// Score pose to populate energy graph (necessary to fill HBondSet):
	get_scorefxn()->score( p );

	// Make the HBondSet object that will contain all the HBond objects:
	core::scoring::hbonds::HBondSet hb_container;

	// Populate it:
	core::scoring::hbonds::fill_hbond_set( p, false /* calculate_derivatives */, hb_container );
	TR << "Done populating HBondSet, "<<hb_container.nhbonds()<<" Hbonds found."<<std::endl;

	// Get the residue index from the pose so you can look for it in the HB container:
	core::conformation::Residue const res = p.residue( get_target_residue() );

	// Get the atom index from the pose so you can look for it in the HB container:
	core::Size const atm_index = res.atom_index( get_target_atom_name() );
	core::id::AtomID atom_id( atm_index, get_target_residue() );

	// Get vector of HBonds involving atom atom_id:
	utility::vector1< core::scoring::hbonds::HBondCOP > Hbond_list( hb_container.atom_hbonds( atom_id ) );

	// Iterate thorugh it to make sure energy is less than hb_e_cutoff
	if ( Hbond_list.size() > 0 ) {
		TR << "Done looking for Hbonds to "<< get_target_atom_name() <<", found "<<Hbond_list.size()<<" Hbonds, here's the list:"<<std::endl;
		for ( auto hbondcop : Hbond_list ) {
			hbondcop->show(TR);
			if ( ( hbondcop->energy() * hbondcop->weight() ) <= get_hb_e_cutoff() ) {
				TR << "This Hbond passes cutoff "<< get_hb_e_cutoff()<<", adding to counter." <<std::endl;
				counter += 1;
			} else {
				TR << "This Hbond doesn't pass cutoff "<< get_hb_e_cutoff()<<", skipping." <<std::endl;
			}
		}
	} else {
		TR << "Done looking for Hbonds to "<< get_target_atom_name() <<", no Hbonds found."<<std::endl;
	}
	TR << "Current number of Hbonds: "<< counter <<"."<<std::endl;

	// Now get the H childs of the target heavy atom and find Hbonds to them.
	TR << "Now looking for Hbonds to H atoms attached to heavy atom."<<std::endl;
	core::kinematics::AtomTree atmtree( p.conformation().atom_tree() );
	debug_assert( atmtree.has( atom_id ));
	core::kinematics::tree::AtomCOP my_heavyatm( atmtree.atom( atom_id ).get_self_ptr() );
	for ( core::Size i=0; i<my_heavyatm->n_nonjump_children(); ++i ) {
		core::kinematics::tree::AtomCOP child( my_heavyatm->child(i));
		core::id::AtomID child_ID( child->atom_id() );
		TR <<"Found child: "<<p.residue( child_ID.rsd() ).atom_name( child_ID.atomno() ) <<std::endl;
		if ( p.residue( child_ID.rsd() ).atom_type( child_ID.atomno() ).is_polar_hydrogen() ) {
			// Look for Hbonds
			utility::vector1< core::scoring::hbonds::HBondCOP > child_Hbond_list( hb_container.atom_hbonds( child_ID ) );
			TR <<"Found "<< child_Hbond_list.size()<<" Hbonds to this child." <<std::endl;
			// If there are Hbonds, iterate through them and check theypass the cutoff.
			if ( child_Hbond_list.size() > 0 ) {
				for ( auto hbondcop : child_Hbond_list ) {
					hbondcop->show(TR);
					if ( ( hbondcop->energy() * hbondcop->weight() ) <= get_hb_e_cutoff() ) {
						TR << "Hbond to child "<< p.residue( child_ID.rsd() ).atom_name( child_ID.atomno() ) <<" passes cutoff "<< get_hb_e_cutoff()<<", adding to counter" <<std::endl;
						counter += 1;
					} else {
						TR << "Hbond to child "<< p.residue( child_ID.rsd() ).atom_name( child_ID.atomno() ) <<" doesn't pass cutoff "<< get_hb_e_cutoff()<<", skipping" <<std::endl;
					}
				}
			} else {
				TR <<"Found no Hbonds to this child." <<std::endl;
			}
		}
	}

	TR <<"Found "<< counter <<" Hbonds, "<<get_n_partners()<<" expected."<<std::endl;
	return counter;
}

void
SimpleHbondsToAtomFilter::report( std::ostream &, core::pose::Pose const & ) const
{

}

std::string SimpleHbondsToAtomFilter::name() const {
	return class_name();
}

std::string SimpleHbondsToAtomFilter::class_name() {
	return "SimpleHbondsToAtomFilter";
}

void SimpleHbondsToAtomFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "n_partners"  , xsct_positive_integer , "Number of partners expected to be Hbonding to tatget atom." , "1" ) + XMLSchemaAttribute::attribute_w_default( "target_atom_name" , xs_string , "Atom of interest involved in Hbonding." , "" ) + XMLSchemaAttribute::attribute_w_default( "hb_e_cutoff" , xsct_real , "Energy of Hbond cutoff." , "-0.5" );

	// Residue number where the atom of interest is. This adds res_num and pbd_num to the attribute list.
	core::pose::attributes_for_get_resnum( attlist ) ;
	// Score function used to fill in the interaction graph.
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist ) ;

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"Check if atom of interest has exactly n_partners Hbonds",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
SimpleHbondsToAtomFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new SimpleHbondsToAtomFilter );
}

std::string
SimpleHbondsToAtomFilterCreator::keyname() const
{
	return SimpleHbondsToAtomFilter::class_name();
}

void SimpleHbondsToAtomFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SimpleHbondsToAtomFilter::provide_xml_schema( xsd );
}

} //protocols
} //simple_filters
