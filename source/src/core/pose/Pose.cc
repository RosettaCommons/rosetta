// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/Pose.cc
/// @brief  Pose class
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu)

// Unit headers
#include <core/pose/Pose.hh>

// Package headers
#include <core/pose/util.hh>
#include <core/pose/reference_pose/ReferencePoseSet.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>
//#include <core/pose/carbohydrates/GlycanTreeSetObserver.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/metrics/PoseMetricContainer.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/io/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/mmcif/cif_writer.hh>

// Basic headers
#include <basic/init.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <basic/prof.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector0.hh>

// C++ headers
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/unordered_map.hpp>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {

using namespace core::conformation;

typedef boost::indirect_iterator< conformation::ResidueOPs::const_iterator, Residue const > const_iterator;
typedef boost::indirect_iterator< conformation::ResidueOPs::iterator, Residue const > iterator;

iterator       Pose::begin() noexcept { return conformation().begin(); }
const_iterator Pose::begin()   const noexcept { return conformation().begin(); }
iterator       Pose::end() noexcept { return conformation().end(); }
const_iterator Pose::end()     const noexcept { return conformation().end(); }


/// @details default init function
void Pose::init(void)
{
#ifdef PYROSETTA
		// Sanity check: check if core::init was called already and abort otherwise with helpful message...
		if( !basic::was_init_called() ) utility_exit_with_message("Attempt to initialize Pose object before core::init was called detected… Have you forgot to call core::init?");
#endif

	conformation_ = ConformationOP( new Conformation() );

	// have the Pose observe it's Conformation for XYZ changes
	// we discard the Link because we own the Conformation
	conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );

	energies_ = scoring::EnergiesOP( new scoring::Energies() );
	energies_->set_owner( this );

	data_cache_ = BasicDataCacheOP( new BasicDataCache( datacache::CacheableDataType::num_cacheable_data_types ) );

	constant_cache_ = ConstDataMapOP( new ConstDataMap );

	observer_cache_ = ObserverCacheOP( new ObserverCache( datacache::CacheableObserverType::num_cacheable_data_types, *this ) );

	metrics_ = metrics::PoseMetricContainerOP( new metrics::PoseMetricContainer );
	metrics_->attach_to( *this );
}


/// @details default constructor
Pose::Pose() :
	pdb_info_( /* NULL */ ),
	constraint_set_( /* 0 */ ),
	reference_pose_set_()
{
	init();
}

/// @details destructor -- > kill data on the heap
Pose::~Pose()
{
	//std::cout << "Pose dstor" << std::endl;
	notify_destruction_obs( DestructionEvent( this ) );
	clear();
}

/// @brief copy constructor
Pose::Pose( Pose const & src ) :
	ReferenceCount ( src ),
	utility::pointer::enable_shared_from_this< Pose >()
{
	*this = src;
}

/// @brief partial copy constructor
Pose::Pose( Pose const & src, Size begin, Size const end):
	pdb_info_( /* NULL */ ),
	constraint_set_( /* 0 */ ),
	reference_pose_set_( )
{
	init();
	utility::vector1< core::Size > residue_indices;
	for ( ; begin <= end; ++begin ) {
		residue_indices.push_back(begin);
	}
	core::io::pose_from_pose(*this, src, residue_indices); //TODO copy reference poses
}

/// @brief copy assignment
Pose &
Pose::operator=( Pose const & src )
{
	PROF_START ( basic::POSE_COPY );

	// Don't send signals to the observers during the incomplete state of a Pose copy
	buffer_observers();

	ConformationOP old_conf = conformation_; // track for observer transfer

	if ( conformation_ && conformation_->same_type_as_me( *src.conformation_, true ) ) {
		(*conformation_) = (*src.conformation_);
	} else {
		conformation_ = src.conformation_->clone();
		conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );
	}
	//std::cout << "Trying to copy the energies!" << std::endl;

	if ( energies_ && energies_->same_type_as_me( *src.energies_ ) ) {
		(*energies_) = (*src.energies_);
	} else {
		energies_ = src.energies_->clone();
		energies_->set_owner( this );
	}

	//std::cout << "Done copying energies, copying data cache" << std::endl;

	// Deep copy of the data held in the non-constant cache
	data_cache_ = BasicDataCacheOP( new BasicDataCache( datacache::CacheableDataType::num_cacheable_data_types ) );
	*data_cache_ = *(src.data_cache_);

	// Shallow copy of the data held in the constant cache
	if ( ! constant_cache_ )  constant_cache_ = ConstDataMapOP( new ConstDataMap );
	*constant_cache_ = *src.constant_cache_;

	// Special observer case: Preserve our PyMOL Observer! (or lack there of)
	// We are specifically NOT cloning any PyMOLObserver objects from src,
	// as this behavior is confusing to users (whom normally create these observers
	// interactively, and are thus surprised when more than one exists because of cloning).
	core::pose::datacache::CacheableObserverOP pymol_observer = nullptr;
	bool was_pymol_observer_attached = false;
	if ( observer_cache_ ) {
		if ( observer_cache_->has( datacache::CacheableObserverType::PYMOL_OBSERVER ) ) {
			pymol_observer = observer_cache_->get_ptr( datacache::CacheableObserverType::PYMOL_OBSERVER );
		}
		was_pymol_observer_attached = observer_cache_->is_attached( datacache::CacheableObserverType::PYMOL_OBSERVER );
	}


	observer_cache_ = ObserverCacheOP( new ObserverCache( datacache::CacheableObserverType::num_cacheable_data_types, *this ) );
	*observer_cache_ = *src.observer_cache_;

	// Restore PyMOLObserver state
	observer_cache_->set( datacache::CacheableObserverType::PYMOL_OBSERVER, pymol_observer, was_pymol_observer_attached );

	this->pdb_info( src.pdb_info_ );
	//std::cout << "Done pdb info" << std::endl;

	metrics_ = metrics::PoseMetricContainerOP( new metrics::PoseMetricContainer( *src.metrics_ ) );
	metrics_->attach_to( *this );
	//std::cout << "Done metrics" << std::endl;

	// TEMP DEBUG CONSTRAINTS
	// if ( constraint_set_ ) constraint_set_->detach_from_conformation();
	// constraint_set_ = src.constraint_set_ ? src.constraint_set_->clone() : src.constraint_set_;
	// if ( constraint_set_ ) constraint_set_->attach_to_conformation( ConformationCAP( conformation_ ));

	// Copy the constraint set -- this will perform a shallow copy of the ConstraintOPs held in
	// the ConstraintSet object (and the Constraints objects that the ConstraintSet holds) if
	// the types of the two constraint set objects match.
	if ( src.constraint_set_ ) {
		if ( constraint_set_ && src.constraint_set_->same_type_as_me( *constraint_set_ ) ) {
			// shallow copy, possibly performing no reference count increments or decrements
			// where possible, and avoids new and delete where possible; esp if several
			// Poses are being copied back and forth into each other
			*constraint_set_ = *src.constraint_set_;
		} else {
			if ( constraint_set_ ) {
				constraint_set_->detach_from_conformation();
			}
			constraint_set_ = src.constraint_set_->clone();
			constraint_set_->attach_to_conformation( ConformationAP( conformation_ ));
		}
	} else {
		if ( constraint_set_ ) {
			constraint_set_->detach_from_conformation();
			constraint_set_.reset();
		}
	}

	//std::cout << "Done constraints" << std::endl;
	//Clone the reference poses:
	if ( src.reference_pose_set_ ) reference_pose_set_ = src.reference_pose_set_->clone();

	// transfer remaining observers that honor the Conformation::TRANSFER
	// event after everything else is done
	if ( old_conf ) {
		conformation_->receive_observers_from( *old_conf );
		old_conf.reset(); // force clear
	}

	unbuffer_observers();
	//std::cout << "Done ref pose and observers" << std::endl;

	PROF_STOP ( basic::POSE_COPY );

	return *this;
}

/// @details This is basically a stub for an actual implementation of detached_copy
/// which would not contain the "receive_observers_from" call that the operator =
/// assignment does.
void
Pose::detached_copy( Pose const & src ) {

	// TEMP! before Pose copying gets refactored, just use the existing operator =
	*this = src;

	// Now perform detached copying on the data members that have observer-
	// or shallow-copy behavior. These are:
	// 1. Conformation (since the AtomTree acts as an observer)
	// 2. Constraint set (since some Constraints contain mutable data)

	if ( conformation_ && conformation_->same_type_as_me( *src.conformation_, true ) ) {
		conformation_->detached_copy( *src.conformation_ );
	} else {
		conformation_ = src.conformation_->clone();
		conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );
	}

	if ( src.constraint_set_ ) {
		if ( constraint_set_ && constraint_set_->same_type_as_me( *src.constraint_set_ ) ) {
			constraint_set_->detached_copy( *src.constraint_set_ );
		} else {
			if ( constraint_set_ ) {
				constraint_set_->detach_from_conformation();
			}
			constraint_set_ = src.constraint_set_->clone();
			constraint_set_->attach_to_conformation( ConformationAP( conformation_ ) );
		}
	} else {
		if ( constraint_set_ ) {
			constraint_set_->detach_from_conformation();
			constraint_set_.reset();
		}
	}

}

/// @brief clone the pose
PoseOP
Pose::clone() const
{
	return PoseOP( new Pose( *this ) );
}

/// @brief Returns the pose Conformation pointer (const access)
ConformationCOP
Pose::conformation_ptr() const
{
	return conformation_;
}

/// @brief Returns the pose Conformation pointer (nonconst access)
ConformationOP &
Pose::conformation_ptr()
{
	return conformation_;
}


kinematics::FoldTree const &
Pose::fold_tree() const
{
	return conformation_->fold_tree();
}

void
Pose::fold_tree( kinematics::FoldTree const & fold_tree_in )
{
	conformation_->fold_tree( fold_tree_in );
}

/// @brief get the atom_tree
kinematics::AtomTree const &
Pose::atom_tree() const
{
	return conformation_->atom_tree();
}

/// @details  Note that we do not clone the input conformation -- we just take it directly. This could be unsafe (?)
///  but it's more efficient. Maybe we want to switch to cloning... Of course we already hand out nonconst refs to
///  our conformation, which is a little unsafe anyhow.
/// @warning Classes observing the Pose's old conformation will not automatically
///  be re-attached/listening to the new Conformation.  Please pay special attention
///  if you have a CacheableObserver in the ObserverCache that listens to
///  a Pose's Conformation.  The prior PDBInfo, ConstraintSet and Energies will be
///  cleared as well.
void
Pose::set_new_conformation( conformation::ConformationCOP new_conformation )
{

	/// drop stuff
	pdb_info_.reset(); // set to NULL
	constraint_set_.reset(); // set to NULL
	observer_cache_->detach();

	// initialize new
	energies_ = scoring::EnergiesOP( new scoring::Energies() );
	energies_->set_owner( this );

	data_cache_ = BasicDataCacheOP( new BasicDataCache( datacache::CacheableDataType::num_cacheable_data_types ) );

	metrics_ = metrics::PoseMetricContainerOP( new metrics::PoseMetricContainer );
	metrics_->attach_to( *this );

	/// clone and reassign the pointer

	conformation_ = new_conformation->clone();

	conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );

	/* OPERATIONS AFTER THIS POINT NEED TO HAPPEN AFTER THE CONFORMATION
	HAS BEEN SWAPPED */

}

/// @details This function allow us to attach an Energies object from a derived class. What are the proper
/// checks to do for this? This should probably be a protected function, since we do not want this function
/// to be used regularly
void
Pose::set_new_energies_object( scoring::EnergiesOP energies )
{
	energies_ = energies->clone();
	// Set the owner of the new energies object
	energies_->set_owner( this );
}


void
Pose::metric( std::string const & calculator_name, std::string const & key, basic::MetricValueBase & val ) const
{ metrics_->get(calculator_name, key, val, *this); return; }

std::string
Pose::print_metric( std::string const & calculator_name, std::string const & key ) const
{ return metrics_->print(calculator_name, key, *this); }

void
Pose::append_residue_by_jump(
	conformation::Residue const & new_rsd,
	Size const jump_anchor_residue,
	std::string const& jump_anchor_atom, // = "",
	std::string const& jump_root_atom, // = "",
	bool const start_new_chain //= false
)
{
	//PyAssert( (jump_anchor_residue>0) && (jump_anchor_residue<size()), "Pose::append_residue_by_jump( ...Size const jump_anchor_residue... ): variable jump_anchor_residue is out of range!" );    // check later: may be fixed in conformation
	energies_->clear(); // TEMPORARY
	conformation_->append_residue_by_jump( new_rsd, jump_anchor_residue, jump_anchor_atom, jump_root_atom, start_new_chain );
	//Since this is appending a new index at the end of the pose, there should be no change to residue index mappings in ReferencePose objects.
}

void
Pose::append_residue_by_bond(
	conformation::Residue const & new_rsd,
	bool const build_ideal_geometry, // = false,
	int const connection, // = 0,
	Size const anchor_residue, // = 0,
	int const anchor_connection, // = 0,
	bool const start_new_chain, // = false
	bool const lookup_bond_length // = false
)
{
	//PyAssert( (anchor_residue>0) && (anchor_residue<=size()), "Pose::append_residue_by_bond( ...Size const anchor_residue... ): variable anchor_residue is out of range!" );    // check later: may be fixed in conformation
	energies_->clear(); // TEMPORARY
	conformation_->append_residue_by_bond( new_rsd, build_ideal_geometry, connection, anchor_residue, anchor_connection, start_new_chain, lookup_bond_length);
	//Since the new residue's index is at the end of the pose, there should be no change to residue mappings in ReferencePose objects.
}

void
Pose::append_residue_by_atoms(
	conformation::Residue const & new_rsd,
	bool const build_ideal_geometry,
	std::string const & connect_atom,
	Size const anchor_rsd_seqpos,
	std::string const & anchor_connect_atom,
	bool const start_new_chain /*= false*/,
	bool const lookup_bond_length /*= false*/
) {
	energies_->clear();
	uint anchor_connection( 0 );
	uint connection( 0 );
	uint connect_atom_index, anchor_connect_atom_index;

	// First, see if the atoms exist.
	conformation::Residue const & anchor_rsd( residue( anchor_rsd_seqpos ) );
	if ( ! new_rsd.has( connect_atom ) || ! anchor_rsd.has( anchor_connect_atom ) ) {
		utility_exit_with_message( "Can't append by these atoms, "
			"since they are not found in the Residues in question!" );
	}

	// Next, convert to indices.
	connect_atom_index = new_rsd.atom_index( connect_atom );
	anchor_connect_atom_index = anchor_rsd.atom_index( anchor_connect_atom );

	// Determine both connections...
	Size n_connections( anchor_rsd.n_possible_residue_connections() );
	for ( Size ii = 1; ii <= n_connections; ++ii ) {
		if ( anchor_rsd.residue_connect_atom_index( ii ) == anchor_connect_atom_index ) {
			anchor_connection = ii;
			break;
		}
	}
	n_connections = new_rsd.n_possible_residue_connections();
	for ( Size ii = 1; ii <= n_connections; ++ii ) {
		if ( new_rsd.residue_connect_atom_index( ii ) == connect_atom_index ) {
			connection = ii;
			break;
		}
	}

	if ( ( ! anchor_connection ) || ( ! connection ) ) {
		utility_exit_with_message( "Can't append by these atoms, "
			"since they do not correspond to connections on the Residues in question!" );
	}

	conformation_->append_residue_by_bond( new_rsd, build_ideal_geometry, connection, anchor_rsd_seqpos,
		anchor_connection, start_new_chain, lookup_bond_length);
}

void
Pose::insert_residue_by_jump(
	Residue const & new_rsd_in,
	Size const seqpos, // desired seqpos of new_rsd
	Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
	std::string const& anchor_atomno, // = "",
	std::string const& root_atomno // = ""
)
{
	PyAssert( (anchor_pos<=size()&&anchor_pos!=0), "Pose::insert_residue_by_jump( ...Size anchor_pos... ): variable anchor_pos is out of range!" );    // check later:
	energies_->clear(); // TEMPORARY
	conformation_->insert_residue_by_jump( new_rsd_in, seqpos, anchor_pos, anchor_atomno, root_atomno );
	increment_reference_pose_mapping_after_seqpos( seqpos ); //All mappings in the new pose after seqpos must be incremented by 1 in all ReferencePose objects.
}

void
Pose::insert_residue_by_bond(
	Residue const & new_rsd_in,
	Size const seqpos, // desired seqpos of new_rsd
	Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
	bool const build_ideal_geometry, // = false,
	std::string const& anchor_atom, // could be ""
	std::string const& root_atom, // ditto
	bool new_chain,
	bool const lookup_bond_length // default false
)
{
	PyAssert( (anchor_pos<=size()&&anchor_pos!=0), "Pose::insert_residue_by_jump( ...Size anchor_pos... ): variable anchor_pos is out of range!" );    // check later:
	energies_->clear(); // TEMPORARY
	conformation_->insert_residue_by_bond( new_rsd_in, seqpos, anchor_pos, build_ideal_geometry, anchor_atom, root_atom, new_chain, lookup_bond_length );
	increment_reference_pose_mapping_after_seqpos( seqpos ); //All mappings in the new pose after seqpos must be incremented by 1 in all ReferencePose objects.
}

void
Pose::replace_residue(
	Size const seqpos,
	Residue const & new_rsd_in,
	bool const orient_backbone
)
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::replace_residue( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later: may become unecessary
	conformation_->replace_residue( seqpos, new_rsd_in, orient_backbone );
	//No change to residue mappings in any ReferencePoses that might exist, since we assume that the replaced residue corresponds to whatever existed previously.
}

void
Pose::replace_residue(
	int const seqpos,
	Residue const & new_rsd_in,
	utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
)
{
	PyAssert( (seqpos<=static_cast<int>(size())&&seqpos!=0), "Pose::replace_residue( ...Size const seqpos... ): "
		"variable seqpos is out of range!" );    // check later: may become unnecessary
	conformation_->replace_residue( seqpos, new_rsd_in, atom_pairs );
	//No change to residue mappings in any ReferencePoses that might exist, since we assume that the replaced residue corresponds to whatever existed previously.
}

void
Pose::append_polymer_residue_after_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry
)
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::append_polymer_residue_after_seqpos( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later: may become unecessary
	energies_->clear(); // TEMPORARY
	conformation_->safely_append_polymer_residue_after_seqpos( new_rsd, seqpos, build_ideal_geometry );
	increment_reference_pose_mapping_after_seqpos( seqpos ); //All mappings in the new pose after seqpos must be incremented by 1 in all ReferencePose objects.
}

void
Pose::prepend_polymer_residue_before_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry
)
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::prepend_polymer_residue_before_seqpos( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later:
	energies_->clear(); // TEMPORARY
	conformation_->safely_prepend_polymer_residue_before_seqpos( new_rsd, seqpos, build_ideal_geometry );
	increment_reference_pose_mapping_after_seqpos( seqpos-1 ); //All mappings in the new pose after seqpos-1 must be incremented by 1 in all ReferencePose objects.
}

void
Pose::append_pose_by_jump(
	Pose const & src,
	Size const jump_anchor_residue,
	std::string const & jump_anchor_atom,
	std::string const & jump_root_atom)
{
	using namespace pose::carbohydrates;
	using namespace conformation::carbohydrates;


	energies_->clear();

	core::Size old_size = size();

	conformation().insert_conformation_by_jump(
		src.conformation(),
		size() + 1,
		num_jump() + 1,
		jump_anchor_residue,
		0, // Set default anchor jump number
		jump_anchor_atom,
		jump_root_atom);

	if ( pdb_info().get() != nullptr && src.pdb_info().get() != nullptr ) {
		pdb_info()->copy(
			*src.pdb_info(),
			1,
			src.pdb_info()->nres(),
			old_size + 1);
	}

	if ( conformation().contains_carbohydrate_residues() && ! glycan_tree_set() ) {
		conformation_->setup_glycan_trees();
	}
	//No change to residue mappings in ReferencePose objects.
}

/// @details delete a polymer residues
///  Fires a LengthEvent::RESIDUE_DELETE signal.
void
Pose::delete_polymer_residue( Size const seqpos )
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::delete_polymer_residue( Size const seqpos ): variable seqpos is out of range!" );
	energies_->clear(); // TEMPORARY
	conformation_->delete_polymer_residue( seqpos );
	zero_reference_pose_mapping_at_seqpos( seqpos ); //All mappings in the new pose pointing to seqpose must now point to 0 in all ReferencePose objects.
	decrement_reference_pose_mapping_after_seqpos( seqpos ); //All mappings in the new pose after seqpos must be decremented by 1 in all ReferencePose objects.
}

/// @details  Delete a residue from the Conformation the slow way -- triggers a rebuild of the atomtree
///  Fires a LengthEvent::RESIDUE_DELETE signal.
/// @note  Could be upstream and/or downstream of a jump or chemical edge, or the root of the tree
/// @note  Not well-tested. Expect funny behavior in new or different situations (email pbradley@fhcrc.org)
///
/// LOGIC: uses fold_tree.delete_seqpos to handle shifting the topology around if necessary, then calls setup_atom_tree
void
Pose::delete_residue_slow( Size const seqpos )
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::delete_residue_slow( Size const seqpos ): variable seqpos is out of range!" );
	energies_->clear(); // TEMPORARY
	conformation_->delete_residue_slow( seqpos );
	zero_reference_pose_mapping_at_seqpos( seqpos ); //All mappings in the new pose pointing to seqpose must now point to 0 in all ReferencePose objects.
	decrement_reference_pose_mapping_after_seqpos( seqpos ); //All mappings in the new pose after seqpos must be decremented by 1 in all ReferencePose objects.
}

/// @brief Delete a range of residues in the pose.
/// @details Calls confromation::delete_residue_range_slow().  Also, updates
/// reference poses, if present.
void Pose::delete_residue_range_slow( Size const start, Size const end) {
	PyAssert( (start>0&&end<=size()), "Pose::delete_residue_range_slow( Size const start, Size const end ): range extends past end(s) of pose!" );
	runtime_assert_string_msg( start <= end, "Error in core::pose::Pose::delete_residue_range_slow(): start must be less than or equal to end." );
	energies_->clear();
	conformation_->delete_residue_range_slow( start, end );
	//Update reference poses, if present:
	for ( core::Size ir=end; ir>=start; --ir ) {
		decrement_reference_pose_mapping_after_seqpos( ir );
	}

	return;
}

void
Pose::copy_segment(
	Size const size,
	Pose const & src,
	Size const begin,
	Size const src_begin
)
{
	conformation_->copy_segment( size, src.conformation(), begin, src_begin );

	// now copy any other data
}


basic::datacache::ConstDataMap const &
Pose::const_data_cache() const
{
	return *constant_cache_;
}

Size
Pose::size() const
{
	return conformation_->size();
}

/// @brief Returns the total number of atoms in the pose conformation
/// example:
///   pose.total_atoms()
Size
Pose::total_atoms() const{
	Size atomno = 0;
	for ( Size ii = 1; ii <= size(); ++ii )  {
		atomno += residue_type( ii ).natoms();
	}
	return atomno;
}

bool
Pose::empty() const
{
	return conformation_->empty();
}

Size
Pose::num_jump() const
{
	return conformation_->fold_tree().num_jump();
}

chemical::AA
Pose::aa( Size const seqpos ) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::aa( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->aa( seqpos );
}

char
Pose::secstruct( Size const seqpos ) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::secstruct( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->secstruct(seqpos);
}

std::string
Pose::secstruct() const {
	std::string ss="";
	for ( Size i = 1; i <= size(); i++ ) {
		ss += secstruct( i );
	}
	return ss;
}

void
Pose::set_secstruct( Size const seqpos, char const setting )
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_secstruct( Size const seqpos , char const setting ): variable seqpos is out of range!" );
	// check variable "setting" to ensure it is logical?
	conformation_->set_secstruct( seqpos, setting );
}


// Sequence accessors
std::string
Pose::sequence() const
{
	std::string seq;
	for ( Size i=1; i<= conformation_->size(); ++i ) {
		seq += residue_type(i).name1();
	}
	return seq;
}

std::string
Pose::sequence(core::Size resnum_start, core::Size resnum_end) const
{
	PyAssert((resnum_start >= 1), "Pose:sequence(core::Size resnum_start, core::Size resnum_end): resnum_start must be greater than or equal to 1!");
	PyAssert((resnum_end <= size()), "Pose:sequence(core::Size resnum_start, core::Size resnum_end): resnum_end must be greater less than or equal to total residues!");
	PyAssert((resnum_end >= resnum_start), "Pose:sequence(core::Size resnum_start, core::Size resnum_end): resnum_end must be greater than or equal to resnum_start!");
	return sequence().substr(resnum_start - 1, resnum_end - resnum_start + 1);

}


std::string
Pose::annotated_sequence( bool show_all_variants ) const
{
	return conformation_->annotated_sequence( show_all_variants );
}

std::string
Pose::chain_sequence(core::Size const chain_in) const
{
	using namespace std;

	debug_assert(chain_in <= conformation_->num_chains());
	PyAssert((chain_in <= conformation_->num_chains()),
		"Pose::chain_sequence(core::Size const chain_in): variable chain_in is out of range!");

	stringstream seq(stringstream::out);

	Size const begin = conformation_->chain_begin(chain_in);
	Size const end = conformation_->chain_end(chain_in);

	if ( !residue_type(begin).is_carbohydrate() ) {
		for ( Size i = begin; i <= end; ++i ) {
			seq << residue_type(i).name1();
		}
	} else /*is carbohydrate*/ {
		// Carbohydrate sequences are listed in the opposite direction as they are numbered.
		for ( Size i = end; i >= begin; --i ) {
			seq << residue_type(i).carbohydrate_info()->short_name_w_linkage_notation();
			if ( i != begin ) {
				seq << "(";
				seq << residue_type(i).carbohydrate_info()->anomeric_carbon();
			}
		}
	}
	return seq.str();
}

/// @brief    Return the IUPAC sequence of a single branch of a glycan tree,
/// Beginning with residue number <start_residue>.
/// @details  This function is called recursively to generate the sequence of all side-branches.
/// @return   Return an empty string if <start_residue> does not correspond to a saccharide residue.
/// @author   Labonte <JWLabonte@jhu.edu>
std::string
Pose::glycan_tree_branch_sequence( core::uint const start_residue ) const
{
	using namespace std;
	using namespace utility;
	using namespace core::conformation::carbohydrates;
	using namespace core::chemical::carbohydrates;

	if ( ! glycan_tree_set()->is_residue_in_tree( start_residue ) ) { return ""; }

	GlycanTreeCOP tree( glycan_tree_set()->get_tree_containing_residue( start_residue ) );
	GlycanNodeCOP current_branch_node( tree->get_node( start_residue ) );

	vector1< string > residue_names;

	// The first residue is special, because it may or may not be the base of the tree.
	CarbohydrateInfoCOP start_residue_info( residue( start_residue ).carbohydrate_info() );
	string start_residue_name( start_residue_info->short_name() );
	core::uint const starting_linkage_position( current_branch_node->get_linkage_position() );
	if ( starting_linkage_position ) {
		start_residue_name +=
			"(" + to_string( start_residue_info->anomeric_carbon() ) +
			"->" + to_string( starting_linkage_position ) + ")";
	}

	residue_names.push_back( start_residue_name );

	core::uint next_main_chain_residue_num( current_branch_node->get_mainchain_child() );
	while ( next_main_chain_residue_num ) {  // If this is 0, then we are at a "tip" and have nothing more to do.
		vector1< core::uint > branch_starts( current_branch_node->get_children() );
		for ( core::uint branch_start : branch_starts ) {
			// Recursively progress along each side-branch off this branch in turn,
			// skipping the main chain of this branch.
			if ( branch_start == next_main_chain_residue_num ) { continue; }
			residue_names.push_back( "[" + glycan_tree_branch_sequence( branch_start ) + "]-" );
		}

		// Now, having dealt with each side-branch, progress along the main chain of this branch.
		current_branch_node = tree->get_node( next_main_chain_residue_num );
		CarbohydrateInfoCOP current_node_info( residue( next_main_chain_residue_num ).carbohydrate_info() );
		residue_names.push_back( current_node_info->short_name() +
			"(" + to_string( current_node_info->anomeric_carbon() ) +
			"->" + to_string( current_branch_node->get_linkage_position() ) + ")-" );

		next_main_chain_residue_num = current_branch_node->get_mainchain_child();
	}

	// Glycan sequences are given from tips to base.
	string sequence;
	for ( core::uint i( residue_names.size() ); i > 0; --i ) {
		sequence += residue_names[ i ];
	}
	return sequence;
}

/// @brief   Return the IUPAC sequence of the entire glycan tree encompassing residue number <residue_in_tree>.
/// @return  Return an empty string if <residue_in_tree> does not correspond to a saccharide residue.
/// @author  Labonte <JWLabonte@jhu.edu>
std::string
Pose::glycan_tree_sequence( core::uint const residue_in_tree ) const
{
	if ( glycan_tree_set()->is_residue_in_tree( residue_in_tree ) ) {
		conformation::carbohydrates::GlycanTreeCOP tree(
			glycan_tree_set()->get_tree_containing_residue( residue_in_tree ) );
		return glycan_tree_branch_sequence( tree->get_start() );
	}
	return "";
}


// Residue at position seqpos
Pose::Residue const &
Pose::residue(
	Size const seqpos
) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::residue( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue( seqpos );
}

basic::datacache::BasicDataCache &
Pose::residue_data(
	Size const seqpos
)
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::residue( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue_data( seqpos );
}

basic::datacache::BasicDataCache const &
Pose::residue_data(
	Size const seqpos
) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::residue( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue_data( seqpos );
}

chemical::ResidueType const &
Pose::residue_type(
	Size const seqpos
) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::residue_type( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue_type( seqpos );
}

chemical::ResidueTypeCOP
Pose::residue_type_ptr(
	Size const seqpos
) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::residue_type_ptr( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue_type_ptr( seqpos );
}

// backbone torsions
// peptides and saccharides

/// @brief Returns the value of the phi backbone dihedral angle.
/// @details  For proteins and peptoids, phi is defined as C(n-1)-N(n)-CA(n)-C(n).
/// For beta-amino acids, phi is defined as C(n-1)-N(n)-CA(n)-CM(n).
/// For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).
/// Et cetera...
Real
Pose::phi( Size const seqpos ) const
{
	using namespace id;

	debug_assert( residue_type( seqpos ).is_protein() || residue_type( seqpos ).is_peptoid() ||
		residue_type( seqpos ).is_carbohydrate() );
	PyAssert( (seqpos <= size() && seqpos != 0 ), "Pose::phi( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type( seqpos ).is_protein() || residue_type( seqpos ).is_peptoid() ||
		residue_type( seqpos ).is_carbohydrate() ),
		"Pose::phi( Size const seqpos ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type( seqpos ).is_protein() || residue_type( seqpos ).is_peptoid() ) {
		if ( residue_type( seqpos ).is_beta_aa() ) {
			return residue( seqpos ).mainchain_torsion(phi_torsion_beta_aa );
		} else if ( residue_type( seqpos ).is_oligourea() ) {
			return residue( seqpos ).mainchain_torsion(phi_torsion_oligourea );
		} else { //Default case, including peptoids and alpha-amino acids:
			return residue( seqpos ).mainchain_torsion( phi_torsion );
		}
	} else /*is carbohydrate*/ {
		return carbohydrates::get_glycosidic_torsion( phi_torsion, *this, seqpos );
	}
}

/// @brief Sets the value of the phi backbone dihedral angle.
/// @details  For proteins, phi is defined as C(n-1)-N(n)-CA(n)-C(n).
/// For beta-amino acids, phi is defined as C(n-1)-N(n)-CA(n)-CM(n).
/// For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).
/// Et cetera...
void
Pose::set_phi( Size const seqpos, Real const setting )
{
	using namespace id;

	debug_assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_phi( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate()),
		"Pose::set_phi( Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		if ( residue_type(seqpos).is_beta_aa() ) {
			conformation_->set_torsion( TorsionID( seqpos, BB, phi_torsion_beta_aa ), setting );
		} else if ( residue_type(seqpos).is_oligourea() ) {
			conformation_->set_torsion( TorsionID( seqpos, BB, phi_torsion_oligourea ), setting );
		} else { //Default case, including peptoids and alpha-amino acids:
			conformation_->set_torsion( TorsionID( seqpos, BB, phi_torsion ), setting );
		}
	} else /*is carbohydrate*/ {
		carbohydrates::set_glycosidic_torsion(phi_torsion, *this, seqpos, setting);
	}
}

/// @brief Returns the value of the psi backbone dihedral angle.
/// @details  For proteins, psi is defined as N(n)-CA(n)-C(n)-N(n+1).
/// For beta-amino acids, psi is defined as CA(n)-CM(n)-C(n)-N(n+1).
/// For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),
/// where X is the position of the glycosidic linkage.
Real
Pose::psi( Size const seqpos ) const
{
	using namespace id;

	debug_assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::psi( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate()),
		"Pose::psi( Size const seqpos ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		if ( residue_type(seqpos).is_beta_aa() ) {
			return residue(seqpos).mainchain_torsion(psi_torsion_beta_aa);
		} else if ( residue_type(seqpos).is_oligourea() ) {
			return residue(seqpos).mainchain_torsion(psi_torsion_oligourea);
		} else { //Default case, including peptoids and alpha-amino acids:
			return residue(seqpos).mainchain_torsion(psi_torsion);
		}
	} else /*is carbohydrate*/ {
		return carbohydrates::get_glycosidic_torsion(psi_torsion, *this, seqpos);
	}
}

/// @brief Sets the value of the psi backbone dihedral angle.
/// @details  For proteins, psi is defined as N(n)-CA(n)-C(n)-N(n+1).
/// For beta-amino acids, psi is defined as CA(n)-CM(n)-C(n)-N(n+1).
/// For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),
/// where X is the position of the glycosidic linkage.
void
Pose::set_psi( Size const seqpos, Real const setting )
{

	using namespace id;

	debug_assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_psi( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate()),
		"Pose::set_psi( Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		if ( residue_type(seqpos).is_beta_aa() ) {
			conformation_->set_torsion( TorsionID( seqpos, BB, psi_torsion_beta_aa ), setting);
		} else if  ( residue_type(seqpos).is_oligourea() ) {
			conformation_->set_torsion( TorsionID( seqpos, BB, psi_torsion_oligourea ), setting);
		} else { //Default case, including peptoids and alpha-amino acids:
			conformation_->set_torsion( TorsionID( seqpos, BB, psi_torsion ), setting );
		}
	} else /*is carbohydrate*/ {
		carbohydrates::set_glycosidic_torsion(psi_torsion, *this, seqpos, setting);
	}
}

/// @brief Returns the value of the omega backbone dihedral angle.
/// @details  For proteins, omega is defined as CA(n)-C(n)-N(n+1)-CA(n+1).
/// For beta-amino acids, omega is defined as CM(n)-C(n)-N(n+1)-CA(n+1).
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-CX-2(n-1),
/// where X is the position of the glycosidic linkage.  (Note that every atom
/// defining this torsion comes from the previous residue!)
Real Pose::omega( Size const seqpos ) const
{
	using namespace id;

	debug_assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::omega( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() ),
		"Pose::omega( Size const seqpos ): residue seqpos is not part of a protein,peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		if ( residue_type(seqpos).is_beta_aa() ) {
			return residue(seqpos).mainchain_torsion(omega_torsion_beta_aa);
		} else if ( residue_type(seqpos).is_oligourea() ) {
			return residue(seqpos).mainchain_torsion(omega_torsion_oligourea);
		} else { //Default case, including peptoids and alpha-amino acids:
			return residue(seqpos).mainchain_torsion(omega_torsion);
		}
	} else /*is carbohydrate*/ {
		return carbohydrates::get_glycosidic_torsion(omega_torsion, *this, seqpos);
	}
}

/// @brief Sets the value of the omega backbone dihedral angle.
/// @details  For proteins, omega is defined as CA(n)-C(n)-N(n+1)-CA(n+1).
/// For beta-amino acids, omega is defined as CM(n)-C(n)-N(n+1)-CA(n+1).
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-OX-1(n-1),
/// where X is the position of the glycosidic linkage.  (Note that every atom
/// defining this torsion comes from the previous residue!)
void
Pose::set_omega( Size const seqpos, Real const setting )
{
	using namespace id;

	debug_assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_omega( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() ),
		"Pose::set_omega( Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		if ( residue_type(seqpos).is_beta_aa() ) {
			conformation_->set_torsion( TorsionID( seqpos, BB, omega_torsion_beta_aa ),  setting);
		} else if ( residue_type(seqpos).is_oligourea() ) {
			conformation_->set_torsion( TorsionID( seqpos, BB, omega_torsion_oligourea ),  setting);
		} else { //Default case, including peptoids and alpha-amino acids:
			conformation_->set_torsion( TorsionID( seqpos, BB, omega_torsion ),  setting );
		}
	} else /*is carbohydrate*/ {
		carbohydrates::set_glycosidic_torsion(omega_torsion, *this, seqpos, setting);
	}
}

/// @brief For a beta-amino acid or oligourea, get the theta backbone dihedral angle.
/// @details Theta is defined as N(n)-CA(n)-CM(n)-C(n) for a beta-amino acid, and
/// N(n)-CA(n)-CM(n)-NU(n) for an oligourea.
Real
Pose::theta( Size const seqpos ) const
{
	using namespace id;

	debug_assert( residue_type(seqpos).is_beta_aa() || residue_type(seqpos).is_oligourea() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::theta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( residue_type(seqpos).is_beta_aa() || residue_type(seqpos).is_oligourea(), "Pose::theta( Size const seqpos ): residue seqpos is not a beta-amino acid or oligourea residue!" );
	if ( residue_type(seqpos).is_beta_aa() ) {
		return residue(seqpos).mainchain_torsion(theta_torsion_beta_aa);
	} else if ( residue_type(seqpos).is_oligourea() ) {
		return residue(seqpos).mainchain_torsion(theta_torsion_oligourea);
	} else {
		// Undefined.
		return 0.0;
	}
}

/// @brief For a beta-amino acid or oligourea residue, set the theta backbone dihedral angle.
/// @details Theta is defined as N(n)-CA(n)-CM(n)-C(n) for a beta-amino acid, and N(n)-CA(n)-CM(n)-NU(n)
/// for an oligourea.
/// @details Theta is defined as N(n)-CA(n)-CM(n)-C(n).
void
Pose::set_theta( Size const seqpos, Real const setting)
{
	using namespace id;

	debug_assert( residue_type(seqpos).is_beta_aa() || residue_type(seqpos).is_oligourea() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_theta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( residue_type(seqpos).is_beta_aa() || residue_type(seqpos).is_oligourea(), "Pose::set_theta( Size const seqpos, Real const setting ): residue seqpos is not a beta-amino acid or oligourea residue!" );
	if ( residue_type(seqpos).is_beta_aa() ) conformation_->set_torsion( TorsionID( seqpos, BB, theta_torsion_beta_aa ), setting );
	else if ( residue_type(seqpos).is_oligourea() ) conformation_->set_torsion( TorsionID(seqpos, BB, theta_torsion_oligourea), setting);
}


/// @brief Returns the mu torsion angle of oligourea residue <seqpos>.
/// @details Mu is defined as CM(n)-NU(n)-C(n)-N(n+1)
/// for an oligourea.
/// @note  Assumes residue is an oligourea.
///
/// example(s):
///     pose.mu(21)
/// See also:
///     Pose
///     Pose.set_mu
///     Pose.residue
///     Residue
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Real
Pose::mu( Size const seqpos ) const {
	debug_assert( residue_type(seqpos).is_oligourea() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::mu( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( residue_type(seqpos).is_oligourea(), "Pose::mu( Size const seqpos ): residue seqpos is not an oligourea residue!" );
	if ( residue_type(seqpos).is_oligourea() /*Silly, redundant check*/ ) {
		return residue(seqpos).mainchain_torsion(id::mu_torsion_oligourea);
	}
	return 0.0; //Undefined.
}

/// @brief Sets the mu torsion angle of oligourea residue <seqpos> to <setting>.
/// @details Mu is defined as CM(n)-NU(n)-C(n)-N(n+1) for an oligourea.
/// @note  <setting>  must be in degrees.  Assumes residue is an oligourea.
///
/// example(s):
///     pose.set_mu(21, 58.9)
/// See also:
///     Pose
///     Pose.mu
///     Pose.residue
///     Residue
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
Pose::set_mu(
	Size const seqpos,
	Real const setting
) {
	debug_assert( residue_type(seqpos).is_oligourea() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_mu( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( residue_type(seqpos).is_oligourea(), "Pose::set_mu( Size const seqpos, Real const setting ): residue seqpos is not an oligourea residue!" );
	if ( residue_type(seqpos).is_oligourea() ) conformation_->set_torsion( id::TorsionID(seqpos, id::BB, id::mu_torsion_oligourea), setting);
}

// nucleic acids

Real
Pose::alpha( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::alpha( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::alpha( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 1 ) );
}

void
Pose::set_alpha( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_alpha( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_alpha( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 1 ), setting );
}

Real
Pose::beta( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::beta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::beta( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 2 ) );
}

void
Pose::set_beta( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_beta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_beta( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 2 ), setting );
}

Real
Pose::gamma( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::gamma( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::gamma( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 3 ) );
}

void
Pose::set_gamma( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_gamma( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_gamma( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 3 ), setting );
}

Real
Pose::delta( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::delta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::delta( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 4 ) );
}

void
Pose::set_delta( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_delta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_delta( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 4 ), setting );
}

Real
Pose::epsilon( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::epsilon( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::epsilon( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 5 ) );
}

void
Pose::set_epsilon( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_epsilon( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_epsilon( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 5 ), setting );
}

Real
Pose::zeta( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::zeta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::zeta( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 6 ) );
}

void
Pose::set_zeta( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_zeta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_zeta( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 6 ), setting );
}

// sidechain torsions
// peptides and saccharides

Real
Pose::chi(
	int const chino,
	Size const seqpos
) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::chi( int const chino , Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (chino>0) && (chino<=static_cast<int>(residue_type(seqpos).nchi())),
		"Pose::chi( int const chino , Size const seqpos ): variable chino innappropriate for this residue!" );

	return residue( seqpos ).chi( chino );
}

void
Pose::set_chi(
	int const chino,
	Size const seqpos,
	Real const setting
)
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_chi( int const chino , Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein()  || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() || residue_type(seqpos).is_ligand() ),
		"Pose::set_chi( int const chino , Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, ligand or carbohydrate!" );
	PyAssert( (chino>0) && (chino<=static_cast<int>(residue_type(seqpos).nchi())),
		"Pose::set_chi( int const chino , Size const seqpos ): variable chino innappropriate for this residue!" );

	conformation_->set_torsion( TorsionID(seqpos, id::CHI, chino), setting);
}

// nucleic acids

Real
Pose::chi( Size const seqpos ) const{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::chi( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::chi( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::CHI, 1 ) );
}

void
Pose::set_chi( Size const seqpos, Real const setting )
{
	debug_assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::set_chi( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_chi( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::CHI, 1 ), setting );
}

// Set the given residue's ring conformation, if appropriate.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @remark  See core/chemical/rings/RingConformerSet.hh and .cc for more information about RingConformers.
void
Pose::set_ring_conformation(
	uint const seqpos,
	uint const ring_num,
	core::chemical::rings::RingConformer const & conformer )
{
	using namespace std;
	using namespace id;
	using namespace numeric;

	Residue const & res( residue( seqpos ) );

	debug_assert( res.type().is_cyclic() );
	PyAssert( ( seqpos <= size() && seqpos != 0 ),
		"Pose::set_ring_conformation(uint const seqpos, uint const ring_num, "
		"core::chemical::rings::RingConformer const & conformer): variable seqpos is out of range!" );
	PyAssert( ( ring_num <= res.type().n_rings() ),
		"Pose::set_ring_conformation(uint const seqpos, uint const ring_num, "
		"core::chemical::rings::RingConformer const & conformer): variable ring_num is out of range!" );

	// First, figure out which nus belong to this ring.
	Size n_nus_on_previous_rings( 0 );
	for ( uint previous_ring_num( ring_num -1 ); previous_ring_num > 0; --previous_ring_num ) {
		// ( Each ring has one fewer nus associated with it than the size of the ring. )
		n_nus_on_previous_rings += res.type().ring_atoms( previous_ring_num ).size() - 1;
	}

	// Then, set the nus, which DEFINE the ideal ring conformer.
	Size const n_nus( res.type().ring_atoms( ring_num ).size() - 1 );
	for ( uint i( 1 ); i <= n_nus; ++i ) {
		set_torsion( TorsionID( seqpos, NU, n_nus_on_previous_rings + i ), conformer.nu_angles[ i ] );
	}

	// Then, set the taus, which result from ring strain.
	Size const n_taus( n_nus + 1 );  // There will always be one more tau stored in the conformer than nus.
	for ( uint i( 1 ); i < n_taus; ++i ) {
		// The reference atoms for the bond angle can be extracted from those used for the corresponding nu angle.
		// For example, nu2 is defined as C1-C2-C3-C4, and tau 1 is defined as C1-C2-C3.
		AtomID const ref1( res.type().nu_atoms( n_nus_on_previous_rings + i )[ 1 ], seqpos );
		AtomID const ref2( res.type().nu_atoms( n_nus_on_previous_rings + i )[ 2 ], seqpos );
		AtomID const ref3( res.type().nu_atoms( n_nus_on_previous_rings + i )[ 3 ], seqpos );
		conformation_->set_bond_angle( ref1, ref2, ref3, conversions::radians( conformer.tau_angles[ i ] ) );
	}

	// Since one fewer nus are stored than taus, we need the LAST 3 reference atoms from the last nu, instead of the
	// 1st 3 atoms as we used above.
	AtomID ref1( res.type().nu_atoms( n_nus_on_previous_rings + n_nus )[ 2 ], seqpos );
	AtomID ref2( res.type().nu_atoms( n_nus_on_previous_rings + n_nus )[ 3 ], seqpos );
	AtomID ref3( res.type().nu_atoms( n_nus_on_previous_rings + n_nus )[ 4 ], seqpos );
	conformation_->set_bond_angle( ref1, ref2, ref3, conversions::radians( conformer.tau_angles[ n_taus ] ) );
}


/////////////////////////////////////////////////////////////////////////////
// generic torsion-angle access

/// @brief  get the torsion angle identified by id
Real
Pose::torsion( TorsionID const & id ) const
{
	return conformation_->torsion( id );
}

/// @brief  set the torsion angle identified by id
void
Pose::set_torsion( TorsionID const & id, Real const setting )
{
	conformation_->set_torsion( id, setting );
}


/////////////////////////////////////////////////////////////////////////////
// chains

core::Size
Pose::num_chains() const
{
	return conformation_->num_chains();
}

core::Size
Pose::chain( Size const seqpos ) const
{
	PyAssert( (seqpos<=size()&&seqpos!=0), "Pose::chain( Size const seqpos ): variable seqpos is out of range!" );
	return residue( seqpos ).chain();
}

/// @details splits the current pose into several poses containing only a single chain each.
utility::vector1<PoseOP>
Pose::split_by_chain() const
{
	utility::vector1<PoseOP> singlechain_poses;

	for ( core::Size i = 1; i <= conformation_->num_chains(); i++ ) {
		singlechain_poses.push_back( split_by_chain(i) );
	}

	return singlechain_poses;
}

/// @details Returns a pose containing only the given chain.
PoseOP
Pose::split_by_chain(Size const chain_id) const
{

	core::pose::PoseOP chain_pose( new Pose(*this) );
	Size chain_begin, chain_end, delete_begin, delete_end;

	chain_begin = chain_pose->conformation().chain_begin( chain_id );
	chain_end = chain_pose->conformation().chain_end( chain_id );

	// if there is only one chain in the pose do nothing and return a copy of the pose
	if ( (conformation_->num_chains() == 1) && (chain_id == 1) ) {}
	// if this is the first chain, delete chain_end to the end of the pose
	else if ( chain_begin == 1 ) {
		delete_begin = chain_end + 1;
		delete_end = chain_pose->size();
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
	} else if ( chain_end == chain_pose->size() ) {
		// if this is the last chain, delete the start of the pose to chain_begin
		delete_begin = 1;
		delete_end = chain_begin - 1;
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
	} else {
		// otherwise, make two deletes around the chain of interest
		delete_begin = 1;
		delete_end = chain_begin - 1;
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
		// just deleted residues --> renumbering pose, so need to reset deletion mask
		delete_begin = chain_end - chain_begin + 2;
		delete_end = chain_pose->size();
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
	}

	// disulfides
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	if ( option[ in::detect_disulf ].user() ?
			option[ in::detect_disulf ]() : // detect_disulf true
			chain_pose->is_fullatom() ) { // detect_disulf default but fa pose
		chain_pose->conformation().detect_disulfides();
	}

	// restore broken pdb_info to new pose ~ Labonte
	if ( chain_pose->pdb_info() != nullptr ) {
		chain_pose->pdb_info()->obsolete(false);
	}

	return chain_pose;
}


// TODO: Move to util.cc.
/// @details  This method updates the pose chain IDs to match the chain IDs
/// found in pdb_info().  In some applications, it is more intuitive to change
/// pdb chain ID letters than it is to change pose chain IDs.  This method
/// adds chain endings between pdb chains and re-derives the pose chain IDs.
/// Currently, disconnected segments with the same pdb chain ID character are
/// treated as separate pose chains, e.g., it is possible for pose chains 1,
/// 3, and 5 to all be chain X.  In the future, I will add a flag to force a
/// one-to-one correspondence between the two chain designations, e.g., if
/// residues 6 through 10 are chain B and residues 1 through 5 AND residues 11
/// through 15 are chain A, then the pose will be reordered to place all res-
/// idues with the same pdb chain ID into a single pose chain. (This is how it
/// works when a pose is loaded from a pdb file.)  I personally have needed
/// use of both functionalities.  I have chosen to create this as a separate
/// method, rather than a part of a change_pdb_chain_ID_of_range(), to avoid
/// multiple calls to Conformation.rederive_chain_ids().  Thus, this method
/// should be called once after all modifications to pdb_info() have been
/// made. ~ Labonte
///
/// See also:
///  PDBInfo.chain()
///  PDBInfo.set_resinfo()
///  Pose.split_by_chain()
///  Conformation.rederive_chain_IDs()
///  Conformation.rederive_chain_endings()
void
Pose::update_pose_chains_from_pdb_chains()
{
	// Declare a vector for storing new (between-residue) chain endings.
	utility::vector1<Size> new_endings;

	char last_pdb_chain = pdb_info_->chain(1);

	for ( Size i = 1; i <= conformation_->size(); ++i ) {
		char current_pdb_chain = pdb_info_->chain(i);
		if ( current_pdb_chain != last_pdb_chain ) {
			new_endings.push_back(i - 1);
			last_pdb_chain = current_pdb_chain;
		}

		// (chain_endings() includes a call to Conformer.rederive_chain_IDs.)
		conformation_->chain_endings(new_endings);
	}
}

/////////////////////////////////////////////////////////////////////////////
// jumps


void
Pose::set_jump(
	int const jump_number,
	const kinematics::Jump & new_jump
)
{
	conformation_->set_jump( jump_number, new_jump );
}


kinematics::Jump const &
Pose::jump( int const jump_number ) const
{
	return conformation_->jump( jump_number );
}


void
Pose::set_jump(
	AtomID const & id,
	const kinematics::Jump & new_jump
)
{
	conformation_->set_jump( id, new_jump );
}


kinematics::Jump const &
Pose::jump( AtomID const & id ) const
{
	return conformation_->jump( id );
}

///////////////////////////////////////////////////////////////////////////
// access atomtree dof's

/// get the value of the atomtree degree of freedom (DOF)
Real
Pose::dof( DOF_ID const & id ) const
{
	return conformation_->dof( id );
}


/// set the value of the atomtree degree of freedom (DOF)
void
Pose::set_dof( DOF_ID const & id, Real const setting )
{
	conformation_->set_dof( id, setting );
}

bool
Pose::has_dof(
	DOF_ID const & did
) const
{
	if ( Size(did.rsd()) > size() || did.rsd() < 1 ) return false;
	if ( Size(did.atomno()) > residue(did.rsd()).natoms() || did.atomno() < 1 ) return false;
	if ( id::PHI == did.type() || id::THETA == did.type() || id::D == did.type() ) {
		if ( conformation_->atom_tree().atom(did.atom_id()).is_jump() ) return false;
	}
	// TODO SHEFFLER MAKE THIS RIGHT!!!!!
	return true;
}


/// get the location of an atom
PointPosition const &
Pose::xyz( AtomID const & id ) const
{
	return conformation_->xyz( id );
}

/// get the location of an atom
PointPosition const &
Pose::xyz( NamedAtomID const & id ) const
{
	return conformation_->residue(id.rsd()).xyz(id.atom());
}

/// set the location of an atom
void
Pose::set_xyz( AtomID const & id, PointPosition const & point )
{
	conformation_->set_xyz( id, point );
}

/// set the location of an atom
void
Pose::set_xyz(
	NamedAtomID const & id,
	PointPosition const & point
)
{
	conformation_->set_xyz( named_atom_id_to_atom_id(id, *this), point );
}

/// set the locations of a vector of atoms
void
Pose::batch_set_xyz( utility::vector1< AtomID > const & ids, utility::vector1< PointPosition > const & points )
{
	conformation_->batch_set_xyz( ids, points );
}


/// get the locations of a vector of atoms
void
Pose::batch_get_xyz( utility::vector1< AtomID > const & ids, utility::vector1< PointPosition > & points ) const
{
	conformation_->batch_get_xyz( ids, points );
}

kinematics::Stub
Pose::stub_from_id(
	id::NamedStubID const& id
)
{
	return conformation_->stub_from_id( named_stub_id_to_stub_id( id, *this ) );
}

void
Pose::center()
{
	PointPosition cog(0,0,0);

	Size count=0;
	for ( Size ir = 1; ir <= size(); ir++ ) {
		for ( Size at = 1; at <= residue( ir ).natoms(); at++ ) {
			cog += xyz( AtomID( at, ir ) );
			count++;
		}
	}

	//std::cout << cog.x() << std::endl;
	cog /= (Real) count;

	//std::cout << cog.x() << std::endl;

	for ( Size ir = 1; ir <= size(); ir++ ) {
		for ( Size at = 1; at <= residue( ir ).natoms(); at++ ) {
			set_xyz(  AtomID( at, ir ),  xyz( AtomID( at, ir )) - cog );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details transfers domain map information into the Energies object, and then
/// resets the domain map information from the Conformation object
void
Pose::update_residue_neighbors()
{
	if ( size() > 0 ) {
		residue( size() ); // completely unnecessary temporary hack to force refold if nec.
	}

	if ( conformation_->structure_moved() ) {
		energies_->structure_has_moved( size() );
		conformation_->reset_structure_moved();
	} else if ( energies_->residue_neighbors_updated() ) {
		return;
	}

	// figure out what's changed since the residue neighbor update
	kinematics::DomainMap domain_map( size() );
	conformation_->update_domain_map( domain_map );

	energies_->update_residue_neighbors( domain_map, *this );
	if ( energies_->discard_conformation_domain_map() ) conformation_->reset_move_data();
}

///////////////////////////////////////////////////////////////////////////////
/// @details called by the ScoreFunction at the start of scoring.  If the score
/// function has changed since the last round of scoring, then cached energies
/// may have become invalidated -- the Energies object makes that decision.
void
Pose::scoring_begin(
	scoring::ScoreFunction const & sfxn
)
{
	// notify of any structure changes
	if ( conformation_->structure_moved() ) {
		conformation_->reset_structure_moved();
		energies_->structure_has_moved( size() );
	}

	energies_->scoring_begin( sfxn, *this );
	// figure out what's changed since the last score evaluation
	/// and update the Energies object
	update_residue_neighbors();
}

void
Pose::scoring_end( scoring::ScoreFunction const & scorefxn )
{
	// reset the data that describes what has moved
	//conformation_->reset_move_data();
	energies_->scoring_end( scorefxn );
	notify_energy_obs( EnergyEvent( this ) );
}

/// @brief called by PairEPotential
void
Pose::update_actcoords()
{
	conformation_->update_actcoords();
}

void
Pose::update_actcoord( Size resid )
{
	conformation_->update_actcoord( resid );
}

void
Pose::update_orbital_coords( Size resid )
{
	conformation_->update_orbital_coords( resid );
}


/// @brief Applies a transform of the form Rx + v, where R is a rotation
/// matrix, V is a vector, and x is the original position in xyz space.
void
Pose::apply_transform_Rx_plus_v(
	numeric::xyzMatrix< Real > const & R,
	Vector const & v
) {
	//fpd move this to conformation since behvaior is very different for symm versus non-symm conformations!
	conformation_->apply_transform_Rx_plus_v(R,v);
}

void
Pose::clear()
{
	conformation_->clear();
	energies_->clear();
	constraint_set_.reset();
	metrics_->clear();
	data_cache_->clear();
	observer_cache_->clear();
	pdb_info( nullptr ); // will check for existence, remove observers, etc
}

// @details Switch residues to VIRTUAL
/// @author Sebastian Raemisch <raemisch@scripps.edu>
/// @author Jared Adolf-Bryfogle <jadolfbr@gmail.com>
void
Pose::real_to_virtual( core::Size seqpos ){
	//std::cout << "Original residue type:" << std::endl;
	//std::cout << this->residue_type(seqpos) << std::endl;


	//work on a copy of the restype. Then add the restype to the pose
	core::chemical::ResidueTypeOP newRT( new core::chemical::ResidueType ( this->residue_type(seqpos) ) );
	newRT->real_to_virtual();
	// Should the new, virtual ResidueType be placed in the pose's RTS, so we don't make a new TR for each replacement?
	replace_pose_residue_copying_existing_coordinates(*this, seqpos, *newRT);

	// update connections
	for ( Size i_con=1; i_con<=this->residue_type(seqpos).n_possible_residue_connections(); ++i_con ) {
		if ( this->conformation().residue(seqpos).connected_residue_at_resconn(i_con) != 0 ) {
			Size connected_seqpos = this->conformation().residue(seqpos).connected_residue_at_resconn(i_con);
			Size connected_id = this->residue(seqpos).connect_map(i_con).connid();
			this->conformation().update_noncanonical_connection(seqpos, i_con, connected_seqpos, connected_id);
		}
	}
}


/// @details Switch residues from VIRTUAL to full_atom
/// @author Sebastian Raemisch <raemisch@scripps.edu>
/// @author Jared Adolf-Bryfogle <jadolfbr@gmail.com>
void
Pose::virtual_to_real( core::Size seqpos ){
	core::chemical::ResidueTypeOP oldRT( new core::chemical::ResidueType ( this->residue_type(seqpos) ) );
	std::string const base_name( residue_type_base_name( *oldRT ) );

	//Remove variant type before getting list.
	oldRT->delete_property("VIRTUAL_RESIDUE");

	utility::vector1< std::string > const & variants( oldRT->properties().get_list_of_variants() );
	chemical::ResidueTypeSetCOP base_RTset( residue_type_set_for_pose( oldRT->mode() ) );
	core::chemical::ResidueTypeCOP RT = core::chemical::ResidueTypeFinder(*base_RTset).residue_base_name( base_name ).variants( variants ).get_representative_type();
	// swap in new (old) residue type
	core::pose::replace_pose_residue_copying_existing_coordinates(*this,seqpos,*RT);
}

conformation::carbohydrates::GlycanTreeSetCOP
Pose::glycan_tree_set() const {
	return conformation_->glycan_tree_set();

}
// @brief Return a const memrbane info object
conformation::membrane::MembraneInfoOP
Pose::membrane_info() const {
	return conformation_->membrane_info();
}

// @brief Returning a nonconst membrane info object
conformation::membrane::MembraneInfoOP
Pose::membrane_info() {
	return conformation_->membrane_info();
}

/// @details Dumps an mmcif formatted file if the extension on the input file_name_string is .cif
void
Pose::dump_file(const std::string & file_name_string) const {
	utility::file::FileName fname( file_name_string );
	if ( fname.ext() == "cif" ) {
		core::io::mmcif::dump_cif( *this, file_name_string );
	} else {
		core::io::pdb::dump_pdb(*this, file_name_string );
	}
}

void
Pose::dump_cif(std::string const &file_name) const {
	core::io::mmcif::dump_cif( *this, file_name );
}


///////////////////////////////////////////////////////////////////////////////
/// @details save pose data to file with supplied file_name
bool
Pose::dump_pdb(std::string const &file_name) const
{
	return core::io::pdb::dump_pdb(*this, file_name);
}

/// @brief  Dump a pdbfile with some score info at the end.
void
Pose::dump_scored_pdb(
	std::string const &file_name,
	scoring::ScoreFunction const & scorefxn
) {
	Real const total_score( scorefxn( *this ) );
	/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( *this );
	std::ofstream out( file_name.c_str() );
	core::io::pdb::dump_pdb( *this, out);
	// verbose output
	out << "END\n";
	std::string secstruct;
	for ( Size i=1; i<= size(); ++i ) secstruct += conformation_->secstruct(i);
	out << "SS: " << secstruct << '\n';
	out << "SCORE_INFO:\n";
	out << "TOTAL_SCORE: " << total_score << '\n';
	scoring::EnergyMap const & wts( scorefxn.weights() );
	out << "WTS: " << wts.show_nonzero() << '\n';
	out << "TOTAL_WTD: " << this->energies().total_energies().weighted_string_of( wts ) << '\n';
	for ( Size i=1; i<= size(); ++i ) {
		out << "RSD_WTD: " << i << ' ' << this->energies().residue_total_energies( i ).weighted_string_of( wts ) << '\n';
	}
	scorefxn.show( out );
	out.close();
}

void Pose::dump_pdb(std::ostream & out) const
{
	return core::io::pdb::dump_pdb(*this, out);
}

/// @brief for writing a specified subset of residues in pdb format
void
Pose::dump_pdb(
	std::ostream & out,
	utility::vector1< Size > const & residue_indices
) const
{
	return core::io::pdb::dump_pdb( *this, out, residue_indices);
}


Pose::ConstraintSetCOP
Pose::constraint_set() const
{
	if ( constraint_set_ == nullptr ) {
		return Pose::ConstraintSetCOP( Pose::ConstraintSetOP( new scoring::constraints::ConstraintSet ) ); // create an empty constraint set
	}
	return constraint_set_;
}

scoring::constraints::ConstraintCOP
Pose::add_constraint( scoring::constraints::ConstraintCOP cst )
{
	energies_->clear();
	if ( constraint_set_ == nullptr ) {
		constraint_set_ = ConstraintSetOP( new scoring::constraints::ConstraintSet ); // create an empty constraint set the first time it's asked for
		constraint_set_->attach_to_conformation( ConformationAP( conformation_ ) );
	}
	scoring::constraints::ConstraintCOP new_cst( cst->clone() );
	constraint_set_->add_constraint( new_cst );
	return( new_cst );
}

scoring::constraints::ConstraintCOPs
Pose::add_constraints( scoring::constraints::ConstraintCOPs csts )
{
	energies_->clear();
	if ( constraint_set_ == nullptr ) {
		constraint_set_ = ConstraintSetOP( new scoring::constraints::ConstraintSet ); // create an empty constraint set the first time it's asked for
		constraint_set_->attach_to_conformation( ConformationAP( conformation_ ) );
	}
	using namespace scoring::constraints;
	ConstraintCOPs new_csts;
	for ( ConstraintCOPs::const_iterator cst_it = csts.begin(); cst_it != csts.end(); ++cst_it ) {
		new_csts.push_back( (*cst_it)->clone() );
	}
	constraint_set_->add_constraints( new_csts );
	return( new_csts );
}

bool
Pose::remove_constraint(
	scoring::constraints::ConstraintCOP cst,
	bool object_comparison )
{
	if ( constraint_set_ == nullptr ) return false;
	energies_->clear();
	return constraint_set_->remove_constraint(cst, object_comparison);
}

bool
Pose::remove_constraints(
	scoring::constraints::ConstraintCOPs csts,
	bool object_comparison )
{
	if ( constraint_set_ == nullptr ) return false;
	energies_->clear();
	return constraint_set_->remove_constraints(csts, object_comparison);
}

bool
Pose::remove_constraints(){
	if ( constraint_set_ == nullptr ) return false;
	energies_->clear();
	constraint_set_->clear();
	return true;
}

void
Pose::clear_sequence_constraints() {
	if ( constraint_set_ == nullptr ) return; //Do nothing if no constraint set.
	constraint_set_->clear_sequence_constraints();
}


/// @details FIX ME! This function should clone all of the constraints held
/// in the input constraint set instead of allowing the input constraint set
/// to perform a shallow copy of its constraints.  Someone may be nefariously
/// holding onto non-const pointers to the constraints in the input constraint
/// set.
void
Pose::constraint_set( ConstraintSetOP constraint_set )
{
	energies_->clear();

	if ( constraint_set_ != nullptr ) constraint_set_->detach_from_conformation();

	if ( constraint_set != nullptr ) {
		constraint_set_ = constraint_set->clone();
		constraint_set_->attach_to_conformation( ConformationAP( conformation_ ) );
	} else constraint_set_ = constraint_set;
}

void Pose::transfer_constraint_set( const pose::Pose &pose ){
	energies_->clear();

	if ( constraint_set_ != nullptr ) constraint_set_->detach_from_conformation();

	if ( pose.constraint_set_ != nullptr ) {
		constraint_set_ = pose.constraint_set_->clone();
		constraint_set_->attach_to_conformation( ConformationAP( conformation_ ) );
	} else constraint_set_ = pose.constraint_set_;
}

//////////////////////////////// ReferencePose and ReferencePoseSet methods /////////////////////////////////////

/// @brief Create a new reference pose from the current state of the pose.
/// @details If a ReferencePoseSet object does not exist, this function will create it.
void Pose::reference_pose_from_current( std::string const &ref_pose_name, bool override_current ) {
	if ( !reference_pose_set_ ) reference_pose_set_= core::pose::reference_pose::ReferencePoseSetOP(new core::pose::reference_pose::ReferencePoseSet); //Create a ReferencePoseSet if it doesn't already exist.

	reference_pose_set_->add_and_initialize_reference_pose( ref_pose_name, *this, override_current );
	return;
}

/// @brief Access the ReferencePoseSet object (non-const).
/// @details If a ReferencePoseSet object does not exist, this function will create it.
core::pose::reference_pose::ReferencePoseSetOP Pose::reference_pose_set()
{
	if ( !reference_pose_set_ ) reference_pose_set_= core::pose::reference_pose::ReferencePoseSetOP(new core::pose::reference_pose::ReferencePoseSet); //Create a ReferencePoseSet if it doesn't already exist.
	return reference_pose_set_;
}

/// @brief Const-access the ReferencePoseSet object.
/// @details If a ReferencePoseSet object does not exist, this function will throw an error.
core::pose::reference_pose::ReferencePoseSetCOP Pose::reference_pose_set_cop() const
{
	runtime_assert_string_msg(reference_pose_set_,
		"Programming error in core::pose::Pose::reference_pose_set_cop(): No ReferencePoseSet object has been created, so it's not possible to request a pointer to the object."
	);
	return core::pose::reference_pose::ReferencePoseSetCOP( reference_pose_set_ );
}

/// @brief Returns the index of a residue in this pose corresponding to a residue in a reference pose.
/// @details Throws an error if the reference pose with the given name doesn't exist, or the residue number
/// doesn't exist in that reference pose.  Returns zero if no corresponding residue exists in this pose (e.g.
/// if the residue in question has been deleted.
core::Size Pose::corresponding_residue_in_current( core::Size const ref_residue_index, std::string const &ref_pose_name ) const
{
	return reference_pose_set_cop()->corresponding_residue_in_current(ref_residue_index, ref_pose_name);
}

/// @brief Find all mappings in the new pose after seqpos in all ReferencePose objects, and increment them by 1.
/// @details If there is no ReferencePose object, do nothing.
void Pose::increment_reference_pose_mapping_after_seqpos( core::Size const seqpos )
{
	if ( !reference_pose_set_ ) return; //Do nothing if there is no ReferencePoseSet object.
	reference_pose_set_->increment_reference_pose_mapping_after_seqpos( seqpos );
	return;
}

/// @brief Find all mappings in the new pose after seqpos in all ReferencePose objects, and decrement them by 1.
/// @details If there is no ReferencePose object, do nothing.
void Pose::decrement_reference_pose_mapping_after_seqpos( core::Size const seqpos )
{
	if ( !reference_pose_set_ ) return; //Do nothing if there is no ReferencePoseSet object.
	reference_pose_set_->decrement_reference_pose_mapping_after_seqpos( seqpos );
	return;
}

/// @brief Find all mappings in the new pose to seqpos in all ReferencePose objects, and set them to point to residue 0 (deletion signal).
/// @details If there is no ReferencePose object, do nothing.
void Pose::zero_reference_pose_mapping_at_seqpos( core::Size const seqpos )
{
	if ( !reference_pose_set_ ) return; //Do nothing if there is no ReferencePoseSet object.
	reference_pose_set_->zero_reference_pose_mapping_at_seqpos( seqpos );
	return;
}

/// @brief Returns true if a pose has at least one reference pose, false otherwise.
///
bool Pose::has_reference_pose() const
{
	if ( !reference_pose_set_ ) return false;
	return ( !reference_pose_set_->empty() );
}

//////////////////////////////// PDBInfo methods /////////////////////////////////////


/// @brief get pdb info (const)
/// @return NULL if no PDBInfo instance exists, the pdb info instance otherwise
PDBInfoCOP
Pose::pdb_info() const
{
	debug_assert( pdb_info_ ? pdb_info_->nres() == size() : true );
	return pdb_info_;
}

/// @brief get pdb info
/// @return NULL if no PDBInfo instance exists, the pdb info instance otherwise
PDBInfoOP
Pose::pdb_info()
{
	debug_assert( pdb_info_ ? pdb_info_->nres() == size() : true );
	return pdb_info_;
}

/// @brief copy new pdb info into this Pose
/// @param[in] new_info  the new pdb info to copy, pass NULL if you want to zero
///  the existence of pdb info inside this Pose
/// @return the prior pdb info instance
PDBInfoOP
Pose::pdb_info( PDBInfoOP new_info )
{
	if ( pdb_info_ ) {
		pdb_info_->detach_from();
	}

	PDBInfoOP prior_pdb_info = pdb_info_;

	if ( new_info ) {
		pdb_info_ = PDBInfoOP( new PDBInfo( *new_info ) ); // make a copy
		pdb_info_->attach_to( *conformation_ );
	} else {
		pdb_info_.reset(); // set to NULL
	}

	PyAssert( pdb_info_ ? pdb_info_->nres() == size() : true, "Invalid PDBInfo, pdb_info_->nres() != size()" );

	debug_assert( pdb_info_ ? pdb_info_->nres() == size() : true );

	return prior_pdb_info;
}

bool Pose::is_fullatom() const {
	return conformation_->is_fullatom();
}

bool Pose::is_centroid() const {
	return conformation_->is_centroid();
}

core::chemical::ResidueTypeSetCOP
Pose::residue_type_set_for_pose( core::chemical::TypeSetMode mode ) const {
	return conformation().residue_type_set_for_conf( mode );
}

void Pose::store_const_data(
	std::string const & category,
	std::string const & name,
	utility::pointer::ReferenceCountCOP ptr
)
{
	constant_cache_->add( category, name, ptr );
}

/// @brief notify DestructionEvent observers
/// @remarks called only upon destruction of the Pose
void
Pose::notify_destruction_obs( DestructionEvent const & e ) {
	destruction_obs_hub_( e );
}


/// @brief notify GeneralEvent observers
/// @remarks should only be called when there are no other suitable event types
///  since specific event notifications will automatically fire a GeneralEvent signal
void
Pose::notify_general_obs( GeneralEvent const & e ) {
	general_obs_hub_( e );
}


/// @brief notify EnergyEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Pose::notify_energy_obs( EnergyEvent const & e, bool const fire_general ) {
	energy_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}

/// @brief notify ConformationEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Pose::notify_conformation_obs( ConformationEvent const & e, bool const fire_general ) {
	conformation_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}

/// @brief Temporarily turn off observer notification
/// Used for places where the Pose is in a temporarily inconsistent state
void
Pose::buffer_observers() const {
	destruction_obs_hub_.buffer();
	general_obs_hub_.buffer();
	energy_obs_hub_.buffer();
	conformation_obs_hub_.buffer();
}

/// @brief Turn back on observer notification
/// Used for places where the Pose is in a temporarily inconsistent state
void
Pose::unbuffer_observers() const {
	destruction_obs_hub_.unblock();
	general_obs_hub_.unblock();
	energy_obs_hub_.unblock();
	conformation_obs_hub_.unblock();
}

/// @brief upon receiving a conformation::signals::XYZEvent
void
Pose::on_conf_xyz_change( core::conformation::signals::XYZEvent const & ) {
	notify_conformation_obs( ConformationEvent( this ) );
}


std::ostream & operator << ( std::ostream & os, Pose const & pose)
{
	PDBInfoCOP p = pose.pdb_info();
	if ( p ) {
		os << "PDB file name: "<< p->name() << std::endl;
	}
	os << "Total residues:" << pose.size() << std::endl;
	os << "Sequence: " << pose.sequence() << std::endl;
	os << "Fold tree:" << std::endl << pose.fold_tree();

	if ( pose.conformation().is_membrane() ) {
		os << pose.conformation().membrane_info() << std::endl;
	}

	// AMW TODO: full model info specific output?
	if ( full_model_info::full_model_info_defined( pose ) ) {
		os << *(full_model_info::const_full_model_info( pose ).full_model_parameters()) << std::endl;
	}

	return os;
}

} // pose
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::Pose::save_with_options( Archive & arc, bool save_observers ) const {
	arc( CEREAL_NVP( conformation_ ) ); // ConformationOP
	arc( CEREAL_NVP( energies_ ) ); // scoring::EnergiesOP
	arc( CEREAL_NVP( metrics_ ) ); // metrics::PoseMetricContainerOP
	arc( CEREAL_NVP( data_cache_ ) ); // BasicDataCacheOP
	arc( CEREAL_NVP( constant_cache_ ) ); // ConstDataMapOP

	if ( save_observers ) arc( CEREAL_NVP( observer_cache_ ) ); // ObserverCacheOP
	// else {
	//  Pose dummy_pose;
	//  auto dummy_observer_cache = ObserverCacheOP( new ObserverCache( datacache::CacheableObserverType::num_cacheable_data_types, dummy_pose) );
	//  arc( CEREAL_NVP( dummy_observer_cache ) );
	// }

	arc( CEREAL_NVP( pdb_info_ ) ); // PDBInfoOP
	arc( CEREAL_NVP( constraint_set_ ) ); // ConstraintSetOP
	arc( CEREAL_NVP( reference_pose_set_ ) ); // core::pose::reference_pose::ReferencePoseSetOP
	// Observers are not serialized arc( CEREAL_NVP( destruction_obs_hub_ ) );
	// Observers are not serialized arc( CEREAL_NVP( general_obs_hub_ ) );
	// Observers are not serialized arc( CEREAL_NVP( energy_obs_hub_ ) );
	// Observers are not serialized arc( CEREAL_NVP( conformation_obs_hub_ ) );
	// EXEMPT destruction_obs_hub_ general_obs_hub_ energy_obs_hub_ conformation_obs_hub_
}

template< class Archive >
void core::pose::Pose::save( Archive & arc) const { save_with_options(arc, true); }


/// @Brief Automatically generated deserialization method
template< class Archive >
void
core::pose::Pose::load_with_options( Archive & arc, bool load_observers ) {
	arc( conformation_ ); // ConformationOP
	conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );

	arc( energies_ ); // scoring::EnergiesOP
	energies_->set_owner( this );

	arc( metrics_ ); // metrics::PoseMetricContainerOP
	arc( data_cache_ ); // BasicDataCacheOP
	arc( constant_cache_ ); // ConstDataMapOP

	if ( load_observers ) {
		arc( observer_cache_ ); // ObserverCacheOP
		observer_cache_->attach_pose( *this );
	}
	// else {
	//  Pose dummy_pose;
	//  auto dummy_observer_cache = ObserverCacheOP( new ObserverCache( datacache::CacheableObserverType::num_cacheable_data_types, dummy_pose) );
	//  arc( dummy_observer_cache );
	// }

	arc( pdb_info_ ); // PDBInfoOP
	if ( pdb_info_ ) {
		pdb_info_->attach_to( *conformation_ );
	}

	arc( constraint_set_ ); // ConstraintSetOP
	if ( constraint_set_ ) constraint_set_->attach_to_conformation( ConformationAP( conformation_ ) );

	arc( reference_pose_set_ ); // core::pose::reference_pose::ReferencePoseSetOP
	// Observers are not serialized arc( destruction_obs_hub_ );
	// Observers are not serialized arc( general_obs_hub_ );
	// Observers are not serialized arc( energy_obs_hub_ );
	// Observers are not serialized arc( conformation_obs_hub_ );
	// EXEMPT destruction_obs_hub_ general_obs_hub_ energy_obs_hub_ conformation_obs_hub_

}

template< class Archive >
void core::pose::Pose::load( Archive & arc) { load_with_options(arc, true); }


template void core::pose::Pose::save_with_options( typename cereal::BinaryOutputArchive &, bool) const;
template void core::pose::Pose::load_with_options( typename cereal::BinaryInputArchive &, bool);


SAVE_AND_LOAD_SERIALIZABLE( core::pose::Pose );
CEREAL_REGISTER_TYPE( core::pose::Pose )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_Pose )
#endif // SERIALIZATION
