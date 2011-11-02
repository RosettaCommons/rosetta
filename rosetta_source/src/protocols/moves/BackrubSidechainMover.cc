// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/BackrubSidechainMover.cc
/// @brief BackrubSidechainMover methods implemented
/// @author


// Unit Headers
#include <protocols/moves/BackrubSidechainMover.hh>
#include <protocols/moves/BackrubSidechainMoverCreator.hh>

// Package Headers

// Project Headers
#include <protocols/moves/BackrubSegment.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ThermodynamicObserver.hh> // needed for Windows build
#include <protocols/rosetta_scripts/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <numeric/MultiDimensionalHistogram.fwd.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
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
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/ResFilterFactory.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/BackrubMover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarloStatus.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/SidechainMover.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key3Vector.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
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
#include <utility/tag/Tag.fwd.hh>
#include <numeric/MultiDimensionalHistogram.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.moves.BackrubSidechainMover" );
static numeric::random::RandomGenerator RG(2772);

namespace protocols {
namespace moves {

std::string
BackrubSidechainMoverCreator::keyname() const {
	return BackrubSidechainMoverCreator::mover_name();
}

protocols::moves::MoverOP
BackrubSidechainMoverCreator::create_mover() const {
	return new BackrubSidechainMover;
}

std::string
BackrubSidechainMoverCreator::mover_name() {
	return "BackrubSidechain";
}

BackrubSidechainMover::BackrubSidechainMover(
) :
	ThermodynamicMover(),
	backrub_mover_(new protocols::moves::BackrubMover),
	sidechain_mover_(new protocols::moves::SidechainMover),
	record_statistics_(false),
	statistics_filename_("brsc_stats.txt")
{
	backrub_mover_->set_min_atoms(7);
	backrub_mover_->set_max_atoms(7);
}

BackrubSidechainMover::BackrubSidechainMover(
	BackrubSidechainMover const & mover
) :
	//utility::pointer::ReferenceCount(),
	ThermodynamicMover(mover),
	valid_segments_(mover.valid_segments_),
	last_valid_segment_index_(mover.last_valid_segment_index_),
	last_chi1_pre_(mover.last_chi1_pre_),
	last_chi1_post_(mover.last_chi1_post_),
	record_statistics_(mover.record_statistics_),
	statistics_filename_(mover.statistics_filename_),
	proposal_hists_(mover.proposal_hists_),
	accept_hists_(mover.accept_hists_)
{
	if (mover.backrub_mover_) {
		backrub_mover_ = dynamic_cast<protocols::moves::BackrubMover *>(mover.backrub_mover_->clone()());
		runtime_assert(backrub_mover_);
	}
	if (mover.sidechain_mover_) {
		sidechain_mover_ = dynamic_cast<protocols::moves::SidechainMover *>(mover.sidechain_mover_->clone()());
		runtime_assert(sidechain_mover_);
	}
}

BackrubSidechainMover::~BackrubSidechainMover(){}

MoverOP
BackrubSidechainMover::clone() const
{
	return new BackrubSidechainMover( *this );
}

MoverOP
BackrubSidechainMover::fresh_instance() const
{
	return new BackrubSidechainMover;
}

std::string
BackrubSidechainMover::get_name() const
{
	return "BackrubSidechainMover";
}

void
BackrubSidechainMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( tag->hasOption("pivot_residues") ) {
		set_pivot_residues(protocols::rosetta_scripts::get_resnum_list(tag, "pivot_residues", pose));
	}

	core::pack::task::TaskFactoryOP new_task_factory( new core::pack::task::TaskFactory );

	if ( tag->hasOption("task_operations") ) {

		std::string const t_o_val( tag->getOption<std::string>("task_operations") );
		typedef utility::vector1< std::string > StringVec;
		StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
					t_o_key != end; ++t_o_key ) {
			if ( data.has( "task_operations", *t_o_key ) ) {
				new_task_factory->push_back( data.get< core::pack::task::operation::TaskOperation* >( "task_operations", *t_o_key ) );
			} else {
				utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
			}
		}

	} else {

		new_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
	}

	set_task_factory(new_task_factory);

	set_prob_uniform( tag->getOption<core::Real>( "prob_uniform", prob_uniform() ) );
	set_prob_withinrot( tag->getOption<core::Real>( "prob_withinrot", prob_withinrot() ) );
	set_prob_random_pert_current( tag->getOption<core::Real>( "prob_random_pert_current", prob_random_pert_current() ) );
	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance() ) );
	set_require_mm_bend( tag->getOption<bool>( "require_mm_bend", require_mm_bend() ) );
	set_record_statistics( tag->getOption<bool>( "record_statistics", record_statistics() ) );
	set_statistics_filename( tag->getOption<std::string>( "statistics_filename", statistics_filename() ) );

	update_segments(pose);
}

void
BackrubSidechainMover::update_segments(
	core::pose::Pose const & pose
)
{
	if (!(backrub_mover_->get_input_pose() && backrub_mover_->get_input_pose()->fold_tree() == pose.fold_tree())) {
		backrub_mover_->set_input_pose(new core::pose::Pose(pose));
	}

	backrub_mover_->set_input_pose(new core::pose::Pose(pose));
	backrub_mover_->clear_segments();
	backrub_mover_->add_mainchain_segments();

	sidechain_mover_->init_task(pose);

	valid_segments_.clear();

	for (core::Size i = 1; i <= backrub_mover_->num_segments(); ++i) {
		BackrubSegment const & segment(backrub_mover_->segment(i));
		if (pose.residue(segment.start_atomid().rsd()).atom_name(segment.start_atomid().atomno()) == " CA " &&
		    pose.residue(segment.end_atomid().rsd()).atom_name(segment.end_atomid().atomno()) == " CA " &&
		    segment.size() == 7) {
			core::Size middle_rsd = (segment.start_atomid().rsd() + segment.end_atomid().rsd())/2;
			if (sidechain_mover_->residue_packed()[middle_rsd]) {
				//TR << "Adding segment start:" << segment.start_atomid() << "end:"
				//	 << segment.end_atomid() << "size: " << segment.size() << std::endl;
				valid_segments_.push_back(i);
			}
		}
	}

	if (record_statistics_) reset_statistics();
}

void
BackrubSidechainMover::initialize_simulation(
	core::pose::Pose & pose,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	if (!(valid_segments_.size() && backrub_mover_->get_input_pose() && backrub_mover_->get_input_pose()->fold_tree() == pose.fold_tree())) {
  	update_segments(pose);
	}

	// because the segments have already been set up, they shouldn't be set up again in the next call
	backrub_mover_->initialize_simulation(pose, metropolis_hastings_mover);
	sidechain_mover_->initialize_simulation(pose, metropolis_hastings_mover);

	if (record_statistics_) reset_statistics();
}

void
BackrubSidechainMover::apply(
	core::pose::Pose & pose
)
{
	if (!(valid_segments_.size() && backrub_mover_->get_input_pose() && backrub_mover_->get_input_pose()->fold_tree() == pose.fold_tree())) {
  	update_segments(pose);
	}

	last_valid_segment_index_ =RG.random_range(1, valid_segments_.size());
	BackrubSegment const & segment(backrub_mover_->segment(valid_segments_[last_valid_segment_index_]));
	core::Size middle_rsd = (segment.start_atomid().rsd() + segment.end_atomid().rsd())/2;

	last_chi1_pre_ = numeric::conversions::radians(pose.chi(1, middle_rsd));

	sidechain_mover_->next_resnum(middle_rsd);
	sidechain_mover_->apply(pose);

	last_chi1_post_ = numeric::conversions::radians(pose.chi(1, middle_rsd));

	backrub_mover_->set_next_segment_id(valid_segments_[last_valid_segment_index_]);
	backrub_mover_->apply(pose);

	update_type();
}

void
BackrubSidechainMover::observe_after_metropolis(
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	if (record_statistics_) record_histograms(metropolis_hastings_mover.monte_carlo()->mc_accepted());
}

void
BackrubSidechainMover::finalize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	if (record_statistics_) {

		std::ostringstream filename;
		if (metropolis_hastings_mover.output_name() != "") {
			filename << metropolis_hastings_mover.output_name() << "_";
		}
		filename << statistics_filename();

		utility::io::ozstream statistics_stream(filename.str());
		output_statistics(statistics_stream);
		statistics_stream.close();
	}
}

utility::vector1<core::Size> const &
BackrubSidechainMover::pivot_residues() const
{
	return backrub_mover_->pivot_residues();
}

void
BackrubSidechainMover::set_pivot_residues(
	utility::vector1<core::Size> const & pivot_residues
)
{
	backrub_mover_->set_pivot_residues(pivot_residues);
}

core::pack::task::TaskFactoryCOP
BackrubSidechainMover::task_factory() const
{
	return sidechain_mover_->task_factory();
}

void
BackrubSidechainMover::set_task_factory(
	core::pack::task::TaskFactoryOP task_factory
)
{
	task_factory->push_back(new core::pack::task::operation::RestrictToRepacking);
	sidechain_mover_->set_task_factory(task_factory);
}

core::Real
BackrubSidechainMover::prob_uniform() const
{
	return sidechain_mover_->prob_uniform();
}

void
BackrubSidechainMover::set_prob_uniform(
	core::Real prob_uniform
)
{
	sidechain_mover_->set_prob_uniform(prob_uniform);
}

core::Real
BackrubSidechainMover::prob_withinrot() const
{
	return sidechain_mover_->prob_withinrot();
}

void
BackrubSidechainMover::set_prob_withinrot(
	core::Real prob_withinrot
)
{
	sidechain_mover_->set_prob_withinrot(prob_withinrot);
}

core::Real
BackrubSidechainMover::prob_random_pert_current() const
{
  return sidechain_mover_->prob_random_pert_current();
}

void
BackrubSidechainMover::set_prob_random_pert_current(
	core::Real prob_pert
)
{
  sidechain_mover_->set_prob_random_pert_current(prob_pert);
}

bool
BackrubSidechainMover::preserve_detailed_balance() const
{
	return backrub_mover_->preserve_detailed_balance() && sidechain_mover_->preserve_detailed_balance();
}

void
BackrubSidechainMover::set_preserve_detailed_balance(
	bool preserve_detailed_balance
)
{
	backrub_mover_->set_preserve_detailed_balance(preserve_detailed_balance);
	sidechain_mover_->set_preserve_detailed_balance(preserve_detailed_balance);
}

bool
BackrubSidechainMover::require_mm_bend() const
{
	return backrub_mover_->require_mm_bend();
}

void
BackrubSidechainMover::set_require_mm_bend(
	bool require_mm_bend
)
{
	backrub_mover_->set_require_mm_bend(require_mm_bend);
}

utility::vector1<core::id::TorsionID_Range>
BackrubSidechainMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}

bool
BackrubSidechainMover::record_statistics() const
{
	return record_statistics_;
}

void
BackrubSidechainMover::set_record_statistics(
	bool record_statistics
)
{
	bool const needs_reset(record_statistics_ != record_statistics);
	record_statistics_ = record_statistics;
	if (needs_reset) reset_statistics();
}

std::string const &
BackrubSidechainMover::statistics_filename() const
{
	return statistics_filename_;
}

void
BackrubSidechainMover::set_statistics_filename(
	std::string const & statistics_filename
)
{
	statistics_filename_ = statistics_filename;
}

void
BackrubSidechainMover::reset_statistics()
{
	setup_histograms();
}

void
BackrubSidechainMover::output_statistics(
	std::ostream & out
)
{
	for (core::Size i = 1; i <= proposal_hists_.size(); ++i) {
		out << proposal_hists_[i] << accept_hists_[i];
	}
}

void
BackrubSidechainMover::setup_histograms()
{
	if (!(record_statistics_ && backrub_mover_->get_input_pose())) {
		proposal_hists_.resize(0);
		accept_hists_.resize(0);
		return;
	}

	core::pose::Pose const & pose(*backrub_mover_->get_input_pose());

	core::Real max_angle_disp(ceil(numeric::conversions::degrees(backrub_mover_->max_angle_disp_7())));
	core::Size const num_bins(static_cast<core::Size>(2*max_angle_disp));
	numeric::conversions::to_radians(max_angle_disp);

	proposal_hists_.resize(valid_segments_.size());
	accept_hists_.resize(valid_segments_.size());

	for (core::Size i = 1; i <= valid_segments_.size(); ++i) {

		core::Size res_num(backrub_mover_->segment(valid_segments_[i]).start_atomid().rsd()+1);

		std::ostringstream res_label_stream;
		res_label_stream << pose.residue(res_num).name3() << " "
		                 << pose.pdb_info()->chain(res_num) << " "
		                 << pose.pdb_info()->number(res_num);

		std::ostringstream proposal_label_stream;
		proposal_label_stream << res_label_stream.str() << " Proposal";
		proposal_hists_[i].label(proposal_label_stream.str());
		proposal_hists_[i].reset_counts();
		proposal_hists_[i].num_dimensions(3);
		proposal_hists_[i].set_dimension(1, num_bins, -max_angle_disp, max_angle_disp, "backrub_disp");
		proposal_hists_[i].set_dimension(2, 3, 0, numeric::constants::r::pi_2, "chi1_pre");
		proposal_hists_[i].set_dimension(3, 3, 0, numeric::constants::r::pi_2, "chi1_post");

		std::ostringstream accept_label_stream;
		accept_label_stream << res_label_stream.str() << " Accept";
		accept_hists_[i].label(accept_label_stream.str());
		accept_hists_[i].reset_counts();
		accept_hists_[i].num_dimensions(3);
		accept_hists_[i].set_dimension(1, num_bins, -max_angle_disp, max_angle_disp, "backrub_disp");
		accept_hists_[i].set_dimension(2, 3, 0, numeric::constants::r::pi_2, "chi1_pre");
		accept_hists_[i].set_dimension(3, 3, 0, numeric::constants::r::pi_2, "chi1_post");
	}
}

void
BackrubSidechainMover::record_histograms(
	bool accepted
)
{
	if (proposal_hists_.size() == 0) setup_histograms();

	//BackrubSegment const & segment(backrub_mover_->segment(backrub_mover_->last_segment_id()));
	//core::Size middle_rsd = (segment.start_atomid().rsd() + segment.end_atomid().rsd())/2;
	//TR << "Changed residue " << middle_rsd << " chi1 from "
	//   << numeric::conversions::degrees(last_chi1_pre_) << " to "
	//   << numeric::conversions::degrees(last_chi1_post_) << " backrub angle "
	//   << numeric::conversions::degrees(backrub_mover_->last_angle()) << " accepted "
	//   << accepted << std::endl;

	utility::vector1<core::Real> values(3);
	values[1] = backrub_mover_->last_angle();
	values[2] = numeric::nonnegative_principal_angle_radians(last_chi1_pre_);
	values[3] = numeric::nonnegative_principal_angle_radians(last_chi1_post_);

	proposal_hists_[last_valid_segment_index_].record(values);
	if (accepted) accept_hists_[last_valid_segment_index_].record(values);
}

void
BackrubSidechainMover::update_type()
{
	std::stringstream mt;

	char bin_letters[] = {'p', 't', 'm'};

	core::Size chi1_pre_bin(static_cast<core::Size>(floor(numeric::nonnegative_principal_angle_radians(last_chi1_pre_)/numeric::constants::r::pi_2_over_3)));
	if (chi1_pre_bin == 3) chi1_pre_bin = 2;
	core::Size chi1_post_bin(static_cast<core::Size>(floor(numeric::nonnegative_principal_angle_radians(last_chi1_post_)/numeric::constants::r::pi_2_over_3)));
	if (chi1_post_bin == 3) chi1_post_bin = 2;

	mt << "brsc_" << bin_letters[chi1_pre_bin] << bin_letters[chi1_post_bin] << "_"
	   << (sidechain_mover_->last_nchi() ? (sidechain_mover_->last_uniform() ? "unif" : (sidechain_mover_->last_withinrot() ? "withinrot" : "rot")) : "none");

	std::string const new_type(mt.str());
	type(new_type);
}

} //moves
} //protocols

