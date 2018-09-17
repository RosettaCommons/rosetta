// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ForceDisulfidesMover.cc
///
/// @brief
/// @author Sarel Fleishman

// unit headers
#include <protocols/calc_taskop_movers/ForceDisulfidesMover.hh>
#include <protocols/calc_taskop_movers/ForceDisulfidesMoverCreator.hh>

// type headers
#include <core/types.hh>
#include <core/id/types.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>

// utility header
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/ResidueIndexDescription.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/task_operations/DesignAroundOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/scoring/ScoreFunction.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace calc_taskop_movers {

static basic::Tracer TR( "protocols.simple_moves.ForceDisulfidesMover" );

ForceDisulfidesMover::ForceDisulfidesMover() :
	protocols::moves::Mover("ForceDisulfidesMover"),
	scorefxn_( /* NULL */ )
{ disulfides_.clear(); }

ForceDisulfidesMover::~ForceDisulfidesMover() = default;

protocols::moves::MoverOP
ForceDisulfidesMover::clone() const
{
	return protocols::moves::MoverOP( new ForceDisulfidesMover( *this ) );
}

protocols::moves::MoverOP
ForceDisulfidesMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ForceDisulfidesMover() );
}

void
ForceDisulfidesMover::apply( Pose & pose ) {
	TR<<"Fixing disulfides"<<std::endl;
	utility::vector1< std::pair< core::Size, core::Size > > parsed_disulfides( disulfides(pose) );

	pose.conformation().fix_disulfides( parsed_disulfides );

	using namespace core::pack::task;
	using namespace protocols::task_operations;
	using core::pack::task::operation::TaskOperationCOP;
	using namespace protocols::minimization_packing;

	DesignAroundOperationOP dao( new DesignAroundOperation );
	dao->design_shell( 0.0 );
	dao->repack_shell( 6.0 );
	for ( auto disulfide_pair: parsed_disulfides ) {
		dao->include_residue( disulfide_pair.first );
		dao->include_residue( disulfide_pair.second );
	}
	TaskFactoryOP tf( new TaskFactory );
	tf->push_back( dao );
	tf->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
	PackerTaskOP ptask = tf->create_task_and_apply_taskoperations( pose );
	PackRotamersMover prm( scorefxn(), ptask );
	TR<<"repacking disulfide surroundings"<<std::endl;
	prm.apply( pose );
}


void
ForceDisulfidesMover::scorefxn( core::scoring::ScoreFunctionOP sf ){
	scorefxn_ = sf;
}

core::scoring::ScoreFunctionOP
ForceDisulfidesMover::scorefxn() const{
	return scorefxn_;
}

// XRW TEMP std::string
// XRW TEMP ForceDisulfidesMover::get_name() const {
// XRW TEMP  return "ForceDisulfidesMover";
// XRW TEMP }

void
ForceDisulfidesMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	utility::vector1< std::string > const residue_pairs( utility::string_split( tag->getOption< std::string >( "disulfides" ), ',' ) );
	TR<<"Setting fix disulfides on residues: ";
	for ( std::string const & residue_pair : residue_pairs ) {
		utility::vector1< std::string > const residues( utility::string_split( residue_pair, ':' ));
		runtime_assert( residues.size() == 2);
		core::pose::ResidueIndexDescriptionCOP res1( core::pose::parse_resnum( residues[ 1 ] ) );
		core::pose::ResidueIndexDescriptionCOP res2( core::pose::parse_resnum( residues[ 2 ] ) );
		disulfides_.push_back( std::make_pair( res1, res2 ) );
		TR<<res1<<':'<<res2<<',';
	}
	TR<<std::endl;
}

void
ForceDisulfidesMover::disulfides( utility::vector1< std::pair < core::Size, core::Size >  > d ){
	disulfides_.clear();
	for ( auto size_pair: d ) {
		core::pose::ResidueIndexDescriptionCOP first( core::pose::make_rid_posenum( size_pair.first ) );
		core::pose::ResidueIndexDescriptionCOP second( core::pose::make_rid_posenum( size_pair.first ) );
		disulfides_.push_back( std::make_pair( first, second ) );
	}
}

std::string ForceDisulfidesMover::get_name() const {
	return mover_name();
}

std::string ForceDisulfidesMover::mover_name() {
	return "ForceDisulfides";
}

void ForceDisulfidesMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	XMLSchemaRestriction pair_cslist;
	pair_cslist.name( "colon_sep_size_pair_cslist" );
	pair_cslist.base_type( xs_string );
	pair_cslist.add_restriction( xsr_pattern, "[+]?[0-9]+:[+]?[0-9]+(,[+]?[0-9]+:[+]?[0-9]+)*" );
	xsd.add_top_level_element( pair_cslist );

	attlist + XMLSchemaAttribute::required_attribute( "disulfides", "colon_sep_size_pair_cslist",
		"For instance: 23A:88A,22B:91B. Can also take regular Rosetta numbering as in: 24:88,23:91." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Set a list of cysteine pairs to form disulfides and repack their surroundings."
		" Useful for cases where the disulfides aren't recognized by Rosetta.", attlist );
}

std::string ForceDisulfidesMoverCreator::keyname() const {
	return ForceDisulfidesMover::mover_name();
}

protocols::moves::MoverOP
ForceDisulfidesMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ForceDisulfidesMover );
}

void ForceDisulfidesMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ForceDisulfidesMover::provide_xml_schema( xsd );
}


utility::vector1< std::pair< core::Size, core::Size > >
ForceDisulfidesMover::disulfides( core::pose::Pose const & pose ) const {
	utility::vector1< std::pair< core::Size, core::Size > > disulfides;
	for ( auto d: disulfides_ ) {
		debug_assert( d.first != nullptr );
		debug_assert( d.second != nullptr );
		core::Size first( d.first->resolve_index( pose ) );
		core::Size second( d.second->resolve_index( pose ) );
		disulfides.push_back( std::make_pair( first, second ) );
	}
	return disulfides;
}

} // moves
} // protocols
