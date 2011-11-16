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
static basic::Tracer TR("protocols.multistate_design.MultiStatePacker",t_info);

#include <ObjexxFCL/format.hh>
using namespace ObjexxFCL::fmt;

#include <boost/functional/hash.hpp>

#include <cmath> // std::exp
#include <vector>
#include <iostream>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <protocols/genetic_algorithm/Entity.fwd.hh>
#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/genetic_algorithm/FitnessFunction.hh>
#include <protocols/multistate_design/MultiStateAggregateFunction.hh>
#include <protocols/multistate_design/MultiStateFitnessFunction.fwd.hh>
#include <protocols/multistate_design/MultiStateFitnessFunction.hh>
#include <protocols/multistate_design/MultiStatePacker.fwd.hh>
#include <protocols/multistate_design/PackingState.fwd.hh>
#include <protocols/multistate_design/SingleState.fwd.hh>
#include <protocols/multistate_design/SingleState.hh>
#include <protocols/multistate_design/SingleStateEntityData.fwd.hh>
#include <protocols/multistate_design/SingleStateEntityData.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.fwd.hh>
#include <protocols/multistate_design/SingleStateFitnessFunction.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
// Auto-header: duplicate removed #include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/factory/WidgetFactory.fwd.hh>
#include <utility/factory/WidgetFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>
// AUTO-REMOVED #include <utility/file/FileName.fwd.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/PathName.fwd.hh>
// AUTO-REMOVED #include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.hh>
// AUTO-REMOVED #include <utility/options/FileOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileOption.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.fwd.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.hh>
// AUTO-REMOVED #include <utility/options/PathOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathOption.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.hh>
// AUTO-REMOVED #include <utility/options/RealOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealOption.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.hh>
// AUTO-REMOVED #include <utility/options/StringOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringOption.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
// AUTO-REMOVED #include <utility/options/mpi_stderr.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <basic/MetricValue.fwd.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>


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
	return new PosType( *this );
}

PosType::EntityElementOP PosType::fresh_instance() {
	return new PosType;
}

Size PosType::hash() const {
	return boost::hash_value( type_ );
}

bool PosType::operator <  ( EntityElement const & rhs ) const
{
	if ( parent::operator < ( rhs ) ) {
		return true;
	} else if ( parent::operator == ( rhs ) ) {
		PosType const * pt_rhs( dynamic_cast< PosType const * > ( &rhs ) );
		if ( ! pt_rhs ) {
			utility_exit_with_message( "operator < unable to compare a " + name() + " object to a " + rhs.name() + " object!" );
		}
		return type_ < pt_rhs->type_;
	}
	return false;
}

bool PosType::operator == ( EntityElement const & rhs ) const
{
	if ( parent::operator == ( rhs ) ) {
		PosType const * pt_rhs( dynamic_cast< PosType const * > ( &rhs ) );
		if ( ! pt_rhs ) {
			utility_exit_with_message( "operator == unable to compare a " + name() + " object to a " + rhs.name() + " object!" );
		}
		if ( type_ == pt_rhs->type_ ) {
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
		PosType const * pt_rhs( dynamic_cast< PosType const * > ( &rhs ) );
		if ( ! pt_rhs ) {
			utility_exit_with_message( "operator = unable to assign to a " + name() + " object from a " + rhs.name() + " object!" );
		}
		type_ = pt_rhs->type_;
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
	return new PosType( word );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
MultiStatePacker::single_state_design( bool restrict_to_canonical /* = true */ )
{
	for ( SingleStateOPs::iterator ss( states().begin() ), end( states().end() );
	      ss != end; ++ss ) {
		PackingStateOP state = dynamic_cast< PackingState* >( (*ss)() );
		runtime_assert( state );
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
	PackingStateOP state = dynamic_cast< PackingState* >( states()[single_state_num]() );
	runtime_assert( state );
	protocols::multistate_design::MultiStateEntity* multi_state_entity =
		dynamic_cast< protocols::multistate_design::MultiStateEntity* >( &entity );

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
			if (multi_state_entity) {
				multi_state_entity->single_state_entity_data()[single_state_num].fitness(E);
				for (MetricValueGetterMap::const_iterator iter = metric_value_getters().begin();
				     iter != metric_value_getters().end(); ++iter) {
					multi_state_entity->single_state_entity_data()[single_state_num].metric_value(
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
		core::chemical::ResidueTypeCAP rot_type( rotsets.rotamer( rot_i )->type() );

		core::chemical::AA seq_type( core::chemical::aa_unk );
		for ( vector1< genetic_algorithm::EntityElementOP >::const_iterator
				it( seq.begin() ), end( seq.end() );
				it != end; ++it ) {
			if ( ! it->get() ) {
				utility_exit_with_message( "Null pointer in EntityElement array" );
			}
			PosType const * postype( dynamic_cast< PosType const * > ( it->get() ) );
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
		}
		else rot_to_pack.push_back(i);
	}
}

} // namespace multistate_design
} // namespace protocols
