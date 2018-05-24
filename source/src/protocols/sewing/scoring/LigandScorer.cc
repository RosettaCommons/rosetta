// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LigandScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/LigandScorer.hh>
#include <protocols/sewing/scoring/LigandScorerCreator.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>

#include <core/scoring/func/AmberPeriodicFunc.hh>
#include <core/scoring/func/Func.hh>

#include <numeric/xyz.functions.hh>
#include <core/scoring/dssp/Dssp.hh>


//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.LigandScorer");

LigandScorer::LigandScorer():
	LigandAssemblyScorer()
{
	core::Real k = 0.5;
	core::Real n_period = ( 5.0 / 4.0 );
	func_ = core::scoring::func::FuncOP( new core::scoring::func::AmberPeriodicFunc( cutoff_angle_, k, n_period ) );
}

LigandScorer::LigandScorer( LigandScorer const & src ):
	LigandAssemblyScorer( src ),
	weight_( src.weight_ ),
	func_( src.func_ ),
	cutoff_distance_(src.cutoff_distance_ ),
	cutoff_angle_( src.cutoff_angle_ )
{
}
core::Real
LigandScorer::get_angle_multiplier( core::Real theta_positive ){
	// TR << "calculating multiplier" << std::endl;
	if ( theta_positive < cutoff_angle_ ) {
		return 1.0;
	} else {
		return func_->func( theta_positive );
	}
}


///@details Use the negative normalized motif score
core::Real
LigandScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	core::chemical::AtomTypeSetCAP atom_types_ = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	set_old_last_score( this->get_last_score() );

	TR << "Beginning ligand score" << std::endl;
	core::Real current_score = 0.0;
	core::Size num_ignored_glycines = 0;
	data_storage::SmartSegmentOP active_segment = assembly->get_n_terminal_segment();
	while ( active_segment != nullptr ) {
		//  TR << "Checking Segment: " << active_segment->get_segment_id() << std::endl;
		for ( core::Size active_resnum = 1; active_resnum <= active_segment->get_length(); ++active_resnum ) {
			//   TR << "#";
			data_storage::SmartSewingResidueOP active_residue = active_segment->get_residue( active_resnum );
			if ( active_residue->get_atom_vector().size() < 5 ) {
				TR << "Resiude" << active_resnum << " in segment " << active_segment->get_segment_id() << " doesn't have a C_beta" << std::endl;
				++num_ignored_glycines;
				continue;
			}
			for ( std::pair< const core::Size, data_storage::LigandResidueOP > const active_ligand : assembly->get_local_ligands() ) {
				//    TR << "*";
				core::Real current_ligand_score = 0.0;
				core::Size num_scored_atoms = 0.0;
				data_storage::LigandResidueOP active_ligand_residue = active_ligand.second;
				for ( core::conformation::Atom active_ligand_atom : active_ligand_residue->get_atom_vector() ) {
					core::Size ligand_atom_type_num = active_ligand_atom.type();
					std::string ligand_atom_type = ( *atom_types_.lock() )[ ligand_atom_type_num ].element();
					if ( ligand_atom_type == "N" || ligand_atom_type == "O" || ligand_atom_type == "H" ) {
						continue; // We don't want to favor burying these
					}
					++num_scored_atoms;
					if ( active_residue->get_atom( 2 ).xyz().distance( active_ligand_atom.xyz()) <= cutoff_distance_ ) {
						//      TR << "found possible packing interaction" << std::endl;
						core::Real theta_positive = numeric::angle_radians( active_residue->get_atom( 5 ).xyz(), active_residue->get_atom( 2 ).xyz(), active_ligand_atom.xyz() );
						current_score = current_score + ( -1.0 * get_angle_multiplier( theta_positive ) );
					}
				}
				if ( !( current_ligand_score > 0.0 ) ) {
					current_ligand_score = current_ligand_score / core::Real( num_scored_atoms );
				}
				current_score = current_score + current_ligand_score;
			}
		}
		//  TR << std::endl;
		active_segment = active_segment->get_c_terminal_neighbor();
	}
	// TR << "done iterating over segments" << std::endl;
	// if( ! (current_score > 0.0 ) ){
	//  current_score = current_score / core::Real ( assembly->get_length() - num_ignored_glycines );
	// }
	set_last_score( current_score );
	TR << "Final ligand score: " << this->get_last_score() << std::endl;
	return this->get_last_score();
}

void
LigandScorer::set_options_from_tag(
	utility::tag::TagCOP scorer,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	if ( !scorer->hasOption( "weight" ) ) {
		TR << "No weight provided for scorer " << scorer->getName() << " Using default value" << std::endl;
	}

	weight_ = scorer->getOption<core::Real>( "weight", 1.0 );
	cutoff_distance_ = scorer->getOption<core::Real>( "ligand_interaction_cutoff_distance", 5.0 );

	TR << "Created LigandScorer from tag." << std::endl;
}

std::string
LigandScorer::get_name() const {
	return "LigandScorer";
}

core::Real
LigandScorer::get_weight() const {
	return weight_;
}

void
LigandScorer::set_weight( core::Real weight ){
	weight_ = weight;
}

core::Real
LigandScorer::get_last_score() const {
	return last_score_;
}

void
LigandScorer::set_last_score( core::Real last_score ){
	last_score_ = last_score;
}

core::Real
LigandScorer::get_old_last_score() const {
	return old_last_score_;
}

void
LigandScorer::set_old_last_score( core::Real old_last_score ){
	old_last_score_ = old_last_score;
}




//Additional getters and setters
core::scoring::func::FuncCOP
LigandScorer::get_func() const{
	return core::scoring::func::FuncCOP( func_ );
}

core::Real
LigandScorer::get_cutoff_distance() const{
	return cutoff_distance_;
}

core::Real
LigandScorer::get_cutoff_angle() const{
	return cutoff_angle_;
}

void
LigandScorer::set_func( core::scoring::func::FuncOP setting ){
	func_ = setting;
}

void
LigandScorer::set_cutoff_distance( core::Real setting){
	cutoff_distance_ = setting;
}

void
LigandScorer::set_cutoff_angle( core::Real setting ){
	cutoff_angle_ = setting;
}

void
LigandScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real,  "How heavily will this term be weighted during scoring?", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "ligand_interaction_cutoff_distance", xsct_real, "The distance cutoff between ligand atom and c alpha that is considered an interaction." , "5.0" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( type_name() )
		.description( "Scores how well ligand is buried based on orientation of nearby Ca's" )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}





//////Creator methods//////////
void
LigandScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	LigandScorer::provide_xml_schema( xsd );

}

AssemblyScorerOP
LigandScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new LigandScorer() );
}

std::string
LigandScorer::type_name(){
	return "LigandScorer";
}


std::string
LigandScorerCreator::keyname() const{
	return LigandScorer::type_name();
}
/////End Creator methods///////






} //scoring namespace
} //sewing namespace
} //protocols namespace
