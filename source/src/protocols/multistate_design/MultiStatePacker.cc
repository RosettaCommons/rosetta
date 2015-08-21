// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MultiStatePacker.hh
/// @brief
/// @author ashworth

#include <protocols/multistate_design/MultiStatePacker.hh>
#include <protocols/multistate_design/PackingState.hh>
#include <protocols/multistate_design/MultiStateEntity.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/string_util.hh>
#include <utility/exit.hh>
using utility::string_split;

#include <utility/vector0.hh>
#include <utility/vector1.hh>
using utility::vector1;

#include <basic/Tracer.hh>
using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static thread_local basic::Tracer TR( "protocols.multistate_design.MultiStatePacker", t_info );

#include <ObjexxFCL/format.hh>
using namespace ObjexxFCL::format;

#include <boost/functional/hash.hpp>

#include <cmath> // std::exp
#include <vector>
#include <iostream>

//Auto Headers
#include <protocols/multistate_design/SingleStateEntityData.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
#include <utility/options/IntegerVectorOption.hh>


namespace protocols {
namespace multistate_design {

genetic_algorithm::EntityElementRegistrator< PosTypeCreator > pt_registrator;

using core::Size;
using core::Real;

PosType::PosType() : parent(), type_(core::chemical::aa_unk) {}
PosType::~PosType() {}
PosType::PosType( Size index, core::chemical::AA type ) : parent( index ), type_(type) {}
PosType::PosType( std::string word ) : parent( word ), type_( core::chemical::aa_unk )
{
	if ( word.size() != 1 ) {
		utility_exit_with_message( "PosType std::string constructor failed to read a one-character word: " + word );
	}
	type_ = core::chemical::aa_from_oneletter_code( word[ 0 ] );
}

PosType::EntityElementOP
PosType::clone() {
	return PosType::EntityElementOP( new PosType( *this ) );
}

PosType::EntityElementOP PosType::fresh_instance() {
	return PosType::EntityElementOP( new PosType );
}

Size PosType::hash() const {
	return boost::hash_value( type_ );
}

bool PosType::operator <  ( EntityElement const & rhs ) const
{
	if ( parent::operator < ( rhs ) ) {
		return true;
	} else if ( parent::operator == ( rhs ) ) {
		if ( ! dynamic_cast< PosType const * > ( &rhs ) ) {
			utility_exit_with_message( "operator < unable to compare a " + name() + " object to a " + rhs.name() + " object!" );
		}
		PosType const & pt_rhs( static_cast< PosType const & > ( rhs ) );
		return type_ < pt_rhs.type_;
	}
	return false;
}

bool PosType::operator == ( EntityElement const & rhs ) const
{
	if ( parent::operator == ( rhs ) ) {
		if ( ! dynamic_cast< PosType const * > ( &rhs ) ) {
			utility_exit_with_message( "operator < unable to compare a " + name() + " object to a " + rhs.name() + " object!" );
		}
		PosType const & pt_rhs( static_cast< PosType const & > ( rhs ) );
		if ( type_ == pt_rhs.type_ ) {
			return true;
		}
	}
	return false;
}

genetic_algorithm::EntityElement const &
PosType::operator =  ( EntityElement const & rhs )
{
	if ( this != &rhs ) {
		parent::operator = ( rhs );

		if ( ! dynamic_cast< PosType const * > ( &rhs ) ) {
			utility_exit_with_message( "operator < unable to compare a " + name() + " object to a " + rhs.name() + " object!" );
		}
		PosType const & pt_rhs( static_cast< PosType const & > ( rhs ) );

		type_ = pt_rhs.type_;
	}
	return *this;
}

std::string
PosType::to_string() const
{
	return parent::to_string() + core::chemical::oneletter_code_from_aa( type_ );
}

std::string
PosType::name() const
{
	PosTypeCreator ptc;
	return ptc.widget_name();
}

core::chemical::AA PosType::type() const { return type_; }

/// PosTypeCreator

PosTypeCreator::~PosTypeCreator() {}

std::string PosTypeCreator::widget_name() const { return "AA"; }

PosTypeCreator::EntityElementOP
PosTypeCreator::new_entity( std::string const & word )
{
	return PosTypeCreator::EntityElementOP( new PosType( word ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultiStatePacker::single_state_design( bool restrict_to_canonical /* = true */ )
{
	for ( SingleStateOPs::iterator ss( states().begin() ), end( states().end() );
			ss != end; ++ss ) {
		PackingStateOP state = utility::pointer::dynamic_pointer_cast< PackingState >( (*ss) );
		runtime_assert( state != 0 );
		utility::vector0< int > rot_to_pack;
		// this is important if alternate states are represented in the rotamer set (e.g. for DNA)
		if ( restrict_to_canonical ) restrict_to_canonical_aas( *state, rot_to_pack );
		state->run_packer( rot_to_pack );
		( *scorefxn() )( state->nonconst_pose() );
		Real const score( state->fitness_function()->evaluate(state->nonconst_pose()) );
		state->set_best_score( score );
		TR(t_info) << "Best single-state design score: " << F(8,2,score) << std::endl;
	}
}

Real
MultiStatePacker::evaluate(
	protocols::genetic_algorithm::Entity & entity,
	core::Size single_state_num
)
{
	PackingStateOP state = utility::pointer::dynamic_pointer_cast< PackingState >( states()[single_state_num] );
	runtime_assert( state != 0 );

	// Filter down to the rotamers needed for this single sequence
	utility::vector0<int> rot_to_pack;
	limit_rotamer_set( rot_to_pack, *state, entity.traits() );

	// optionally pack multiple times to find best energy
	Real E(0.), bestE(0.);
	for ( Size i(0); i < num_packs_; ++i ) {

		state->run_packer( rot_to_pack );

		( *scorefxn() )( state->nonconst_pose() );
		E = state->fitness_function()->evaluate(state->nonconst_pose());
		if ( E < bestE || i == 0 ) {
			bestE = E;
			if ( dynamic_cast< protocols::multistate_design::MultiStateEntity * >( &entity ) ) {
				protocols::multistate_design::MultiStateEntity & multi_state_entity =
					static_cast< protocols::multistate_design::MultiStateEntity & >( entity );
				multi_state_entity.single_state_entity_data()[single_state_num].fitness(E);
				for ( MetricValueGetterMap::const_iterator iter = metric_value_getters().begin();
						iter != metric_value_getters().end(); ++iter ) {
					multi_state_entity.single_state_entity_data()[single_state_num].metric_value(
						iter->first,
						iter->second.get(state->pose())
					);
				}
			}
		}
	}
	return bestE;
}

void
limit_rotamer_set(
	utility::vector0< int > & rot_to_pack,
	PackingState const & state,
	genetic_algorithm::EntityElements const & seq
)
{
	rot_to_pack.clear();

	// Allocate enough to accomodate full rotamer set
	core::pack::rotamer_set::RotamerSets const & rotsets( *state.rotamersets() );
	Size const nrotamers( rotsets.nrotamers() );

	for ( Size rot_i(1); rot_i <= nrotamers; ++rot_i ) {

		Size const rot_pos( rotsets.res_for_rotamer( rot_i ) );
		core::chemical::ResidueTypeCOP rot_type( rotsets.rotamer( rot_i )->type().get_self_ptr() );

		core::chemical::AA seq_type( core::chemical::aa_unk );
		for ( vector1< genetic_algorithm::EntityElementOP >::const_iterator
				it( seq.begin() ), end( seq.end() );
				it != end; ++it ) {
			if ( ! it->get() ) {
				utility_exit_with_message( "Null pointer in EntityElement array" );
			}
			PosTypeCOP postype( utility::pointer::dynamic_pointer_cast< PosType const > ( *it ) );
			if ( ! postype ) {
				utility_exit_with_message( "Dynamic cast to PosType failed for object of type " + (*it)->name() );
			}
			if ( postype->index() == rot_pos ) { seq_type = postype->type(); break; }
		}

		if ( seq_type == core::chemical::aa_unk ) {
			// this is not a position whose mutation is controlled by genetic algorithm,
			// and should just repack for this state--allow the rotamer if it is not a mutation
			if ( rot_type->aa() != state.pose().residue_type( rot_pos ).aa() ) continue;
			rot_to_pack.push_back( rot_i );
		} else {
			// this is not a position whose mutation is controlled by genetic algorithm:
			// accept this rotamer only if its identity matches that which is specified in seq
			if ( rot_type->aa() == seq_type ) rot_to_pack.push_back( rot_i );
		}
	}
}

void
restrict_to_canonical_aas(
	PackingState const & state,
	utility::vector0< int > & rot_to_pack
)
{
	rot_to_pack.clear();
	core::pack::rotamer_set::RotamerSets const & rotsets( *state.rotamersets() );
	Size const nrot( rotsets.nrotamers() );

	for ( Size i(1); i <= nrot; ++i ) {
		core::conformation::Residue const & rot( *rotsets.rotamer(i) );
		if ( rot.aa() > core::chemical::num_canonical_aas ) {
			// if non-canonical rotamer, restrict to existing residue type at this position
			if ( rot.type().name3() ==
					state.pose().residue( rotsets.res_for_rotamer(i) ).type().name3() ) {
				rot_to_pack.push_back(i);
			}
		} else rot_to_pack.push_back(i);
	}
}

} // namespace multistate_design
} // namespace protocols
