// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/PolymerBondedEnergyContainer.cc
/// @brief  A container interface long range energies for polymer-bonded residue interactions only.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Modified 21 February 2016 so that this no longer just scores n->n+1 interactions, but includes
/// anything that's polymer-bonded (whether or not it's adjacent in linear sequence).  This includes N-to-C cyclic peptide bonds.

// Unit headers
#include <core/scoring/PolymerBondedEnergyContainer.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

PolymerBondedNeighborIterator::~PolymerBondedNeighborIterator(){}

PolymerBondedNeighborIterator::PolymerBondedNeighborIterator(
	Size const base_in,
	utility::vector1< Size > const & pos_in,
	PolymerBondedEnergyContainer & parent
):
	base_( base_in ),
	pos_( pos_in ),
	parent_( &parent )
{
	curr_idx_=1;
}

ResidueNeighborIterator & PolymerBondedNeighborIterator::operator = ( ResidueNeighborIterator const & src ) {
	debug_assert( dynamic_cast< PolymerBondedNeighborIterator const * >( &src ) );
	PolymerBondedNeighborIterator const & my_src( static_cast< PolymerBondedNeighborIterator const & >( src ) );
	base_ = my_src.base_;
	pos_ = my_src.pos_;
	parent_ = my_src.parent_;

	return *this;
}

ResidueNeighborIterator const & PolymerBondedNeighborIterator::operator ++ () {
	curr_idx_++;
	return *this;
}

bool PolymerBondedNeighborIterator::operator == ( ResidueNeighborIterator const & other ) const
{
	return (
		residue_iterated_on() == other.residue_iterated_on() &&
		neighbor_id() == other.neighbor_id() );
}

bool PolymerBondedNeighborIterator::operator != ( ResidueNeighborIterator const & other ) const {
	return !( *this == other );
}

Size PolymerBondedNeighborIterator::upper_neighbor_id() const {
	core::Size other = 0;
	if ( curr_idx_ <= pos_.size() ) other = pos_[curr_idx_];
	return other;
}

Size PolymerBondedNeighborIterator::lower_neighbor_id() const {
	return base_;
}

Size PolymerBondedNeighborIterator::residue_iterated_on() const {
	return base_;
}

Size PolymerBondedNeighborIterator::neighbor_id() const {
	core::Size other = 0;
	if ( curr_idx_ <= pos_.size() ) other = pos_[curr_idx_];
	return other;
}

void PolymerBondedNeighborIterator::save_energy( EnergyMap const & emap ) {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		Real const energy( emap[ parent_->score_types()[i] ] );
		parent_->set_energy( base_ , i , energy);
	}
}

void PolymerBondedNeighborIterator::retrieve_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] = parent_->get_energy( base_ , i );
	}
}

void PolymerBondedNeighborIterator::accumulate_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] += parent_->get_energy( base_ , i );
	}
}

void PolymerBondedNeighborIterator::mark_energy_computed() {
	parent_->set_computed( base_ , true );
}

void PolymerBondedNeighborIterator::mark_energy_uncomputed() {
	parent_->set_computed( base_ , false );
}

bool PolymerBondedNeighborIterator::energy_computed() const {
	return parent_->get_computed( base_ );
}

/////////////////////////////////////////////////////

PolymerBondedNeighborConstIterator::~PolymerBondedNeighborConstIterator(){}

PolymerBondedNeighborConstIterator::PolymerBondedNeighborConstIterator(
	Size const base_in,
	utility::vector1< Size > const & pos_in,
	PolymerBondedEnergyContainer const & parent
):
	base_( base_in ),
	pos_( pos_in ),
	parent_( &parent )
{
	curr_idx_=1;
}

ResidueNeighborConstIterator & PolymerBondedNeighborConstIterator::operator = ( ResidueNeighborConstIterator const & src ) {
	debug_assert( dynamic_cast< PolymerBondedNeighborConstIterator const * >( &src ) );
	PolymerBondedNeighborConstIterator const & my_src( static_cast< PolymerBondedNeighborConstIterator const & >( src ) );
	base_ = my_src.base_;
	pos_ = my_src.pos_;
	parent_ = my_src.parent_;
	return *this;
}

ResidueNeighborConstIterator const & PolymerBondedNeighborConstIterator::operator ++ () {
	curr_idx_++;
	return *this;
}

bool PolymerBondedNeighborConstIterator::operator == ( ResidueNeighborConstIterator const & other ) const {
	return (
		residue_iterated_on() == other.residue_iterated_on() &&
		neighbor_id() == other.neighbor_id() );
}

bool PolymerBondedNeighborConstIterator::operator != ( ResidueNeighborConstIterator const & other ) const {
	return !( *this == other );
}

Size PolymerBondedNeighborConstIterator::upper_neighbor_id() const {
	core::Size other = 0;
	if ( curr_idx_ <= pos_.size() ) other = pos_[curr_idx_];
	return other;
}

Size PolymerBondedNeighborConstIterator::lower_neighbor_id() const {
	return base_;
}

Size PolymerBondedNeighborConstIterator::residue_iterated_on() const {
	return base_;
}

Size PolymerBondedNeighborConstIterator::neighbor_id() const {
	core::Size other = 0;
	if ( curr_idx_ <= pos_.size() ) other = pos_[curr_idx_];
	return other;
}

void PolymerBondedNeighborConstIterator::retrieve_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] = parent_->get_energy( base_ , i );
	}
}

void PolymerBondedNeighborConstIterator::accumulate_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] += parent_->get_energy( base_ , i );
	}
}

bool PolymerBondedNeighborConstIterator::energy_computed() const {
	return parent_->get_computed( base_ );
}

/////////////////////////////////////////////////////////////////////////

/// @brief Destructor.
///
PolymerBondedEnergyContainer::~PolymerBondedEnergyContainer() {}

LREnergyContainerOP PolymerBondedEnergyContainer::clone() const {
	return LREnergyContainerOP( new PolymerBondedEnergyContainer( *this ) );
}

/// @brief Pose constructor.
/// @details Initializes PolymerBondedEnergyContainer from a pose, facilitating calculations involving non-canonical connections
/// (e.g. terminal peptide bonds).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
PolymerBondedEnergyContainer::PolymerBondedEnergyContainer(
	core::pose::Pose const & pose,
	utility::vector1< ScoreType > const & score_type_in
) :
	size_(0),
	score_types_( score_type_in )
{
	initialize_peptide_bonded_pair_indices( pose ); //Sets size_ and peptide_bonded_pair_indices_.
	computed_.resize( size_, false);
	tables_.resize( size_, utility::vector1< core::Real >( score_types_.size(), 0.0) );
}


bool PolymerBondedEnergyContainer::empty() const {
	return ( size_ == 0 );
}

bool
PolymerBondedEnergyContainer::any_neighbors_for_residue( int /*resid*/ ) const {
	return true;
}

bool
PolymerBondedEnergyContainer::any_upper_neighbors_for_residue( int /*resid*/ ) const {
	return true;
}

Size
PolymerBondedEnergyContainer::size() const {
	return size_;
}

ResidueNeighborConstIteratorOP
PolymerBondedEnergyContainer::const_neighbor_iterator_begin( int resid ) const {
	if ( resid <= (int)(0) || resid > static_cast<int>(size_) ) {
		//If out of range, return no interactions
		return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, utility::vector1<core::Size>(), *this ) );
	}

	utility::vector1<core::Size> neighbors;
	std::map <core::Size, core::Size >::const_iterator it = peptide_bonded_pair_indices_.find(resid);
	if ( it != peptide_bonded_pair_indices_.end() ) {
		neighbors.push_back(it->second);
	}
	std::map <core::Size, core::Size >::const_iterator itinv = inv_peptide_bonded_pair_indices_.find(resid);
	if ( itinv != inv_peptide_bonded_pair_indices_.end() ) {
		neighbors.push_back(itinv->second);
	}

	return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, neighbors, *this ) );
}

ResidueNeighborConstIteratorOP
PolymerBondedEnergyContainer::const_neighbor_iterator_end( int resid ) const {
	return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, utility::vector1<core::Size>(), *this ) );
}

ResidueNeighborConstIteratorOP
PolymerBondedEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const {
	if ( resid <= (int)(0) || resid > static_cast<int>(size_) ) {
		//If out of range, return no interactions
		return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, utility::vector1<core::Size>(), *this ) );
	}

	utility::vector1<core::Size> neighbors;
	std::map <core::Size, core::Size >::const_iterator it = peptide_bonded_pair_indices_.find(resid);
	if ( it != peptide_bonded_pair_indices_.end() ) {
		neighbors.push_back(it->second);
	}

	return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, neighbors, *this ) );
}

ResidueNeighborConstIteratorOP
PolymerBondedEnergyContainer::const_upper_neighbor_iterator_end( int resid ) const {
	return const_neighbor_iterator_end( resid );
}

//////////////////// non-const versions
ResidueNeighborIteratorOP
PolymerBondedEnergyContainer::neighbor_iterator_begin( int resid ) {
	if ( resid <= (int)(0) || resid > static_cast<int>(size_) ) {
		//If out of range, return no interactions
		return ResidueNeighborIteratorOP( new PolymerBondedNeighborIterator( resid, utility::vector1<core::Size>(), *this ) );
	}

	utility::vector1<core::Size> neighbors;
	if ( peptide_bonded_pair_indices_.find(resid) != peptide_bonded_pair_indices_.end() ) {
		neighbors.push_back(peptide_bonded_pair_indices_[resid]);
	}
	if ( inv_peptide_bonded_pair_indices_.find(resid) != inv_peptide_bonded_pair_indices_.end() ) {
		neighbors.push_back(inv_peptide_bonded_pair_indices_[resid]);
	}

	return ResidueNeighborIteratorOP( new PolymerBondedNeighborIterator( resid, neighbors, *this ) );
}

ResidueNeighborIteratorOP
PolymerBondedEnergyContainer::neighbor_iterator_end( int resid ) {
	return ResidueNeighborIteratorOP( new PolymerBondedNeighborIterator( resid, utility::vector1<core::Size>(), *this ) );
}

ResidueNeighborIteratorOP
PolymerBondedEnergyContainer::upper_neighbor_iterator_begin( int resid )
{
	if ( resid <= (int)(0) || resid > static_cast<int>(size_) ) {
		//If out of range, return no interactions
		return ResidueNeighborIteratorOP( new PolymerBondedNeighborIterator( resid, utility::vector1<core::Size>(), *this ) );
	}

	utility::vector1<core::Size> neighbors;
	if ( peptide_bonded_pair_indices_.find(resid) != peptide_bonded_pair_indices_.end() ) {
		neighbors.push_back(peptide_bonded_pair_indices_[resid]);
	}

	return ResidueNeighborIteratorOP( new PolymerBondedNeighborIterator( resid, neighbors, *this ) );
}

ResidueNeighborIteratorOP
PolymerBondedEnergyContainer::upper_neighbor_iterator_end( int resid ) {
	return neighbor_iterator_end( resid );
}

/// @brief Is this PolymerBondedEnergyContainer properly set up for the pose?
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
PolymerBondedEnergyContainer::is_valid(
	core::pose::Pose const &pose
) const {
	//First, check that the number of residues is correct:
	core::Size nres( pose.n_residue() );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres = core::pose::symmetry::symmetry_info(pose)->last_independent_residue();
	}
	if ( size_ != nres ) return false;

	//Next, check that the correct connections have been set up.
	core::Size connect_count(0);
	bool correct_connections_found(true);
	for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all residues (or all residues in the asymmetric unit, if this is a symmetric pose).
		core::Size other_res( other_res_index( pose, ir ) ); //Get the index of the residue that this one is connected to at its upper_connect (if any).
		if ( !other_res ) continue; //These are not connected by normal polymeric connection.

		//OK, these two residues are connected by a peptide bond (or, at least, a polymer bond).  Check that they're stored:
		if ( peptide_bonded_pair_indices_.count(ir) == 0 || peptide_bonded_pair_indices_.at(ir) != other_res ) {
			correct_connections_found = false;
			break;
		}

		//Increment the count of connections found.  At the end, this should match the number of elements in the peptide_bonded_pair_indices_ map.
		++connect_count;
	} //End loop through all residues (in asymmetric unit, in symmetric case).

	return ( correct_connections_found && ( connect_count == peptide_bonded_pair_indices_.size() ) );
}

/// @brief Given a pose, set up a list of pairs of peptide-bonded residues.
/// @details The first index is the lower residue (the residue contributing the carboxyl/UPPER_CONNECT) and the second is the
/// upper residue (the residue contributing the amide/LOWER_CONNECT).
/// @note Since the PolymerBondedEnergyContainer is also apparently used for RNA, I'm going to base this on UPPER_CONNECT/LOWER_CONNECT,
/// not on any peptide-specific geometry.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
PolymerBondedEnergyContainer::initialize_peptide_bonded_pair_indices(
	core::pose::Pose const &pose
) {
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		size_ = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	} else {
		size_ = pose.total_residue();
	}

	peptide_bonded_pair_indices_.clear();
	inv_peptide_bonded_pair_indices_.clear();
	for ( core::Size ir=1; ir<=size_; ++ir ) { //Loop through all resiudes (asymmetric case) or residues in asymmetric unit (symmetric case).
		core::Size other_res( other_res_index( pose, ir ) ); //Get the index of the residue that this one is connected to at its upper_connect (if any).
		if ( !other_res ) continue; //These are not connected by normal polymeric connection.

		//OK, these two residues are connected by a peptide bond (or, at least, a polymer bond).  Store them:
		debug_assert(peptide_bonded_pair_indices_.count(ir) == 0); //Should always be true.
		debug_assert(inv_peptide_bonded_pair_indices_.count(other_res) == 0); //Should always be true.
		peptide_bonded_pair_indices_[ir] =  other_res;
		inv_peptide_bonded_pair_indices_[other_res] =  ir;
	}
}

/// @brief Logic to get the residue that this residue is connected to at its upper_connect, and which is connected to this residue at its lower_connect.
/// @details Returns 0 if this residue is not polymeric, if it lacks an upper_connect, if it's not connected to anything at its upper connect, if the
/// thing that it's connected to is not polymeric, if the thing that it's connected to lacks a lower_connect, or if the thing that it's connected to is not
/// connected to it at its lower connect.  Phew, what a mouthful!
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::Size
PolymerBondedEnergyContainer::other_res_index(
	core::pose::Pose const &pose,
	core::Size const this_res_index
) const {
	//fpd if this logic is updated, the logic in CartesianBondedEnergy::eval_intrares_energy must similarly be updated

	if ( !pose.residue(this_res_index).type().is_polymer() ) return 0; //Skip this residue if it's not a polymer.
	if ( !pose.residue(this_res_index).has_upper_connect() ) return 0; //Skip this residue if it's not one that has an upper_connect.
	core::Size const other_res( pose.residue(this_res_index).residue_connection_partner( pose.residue(this_res_index).upper_connect().index() ) ); //Index of residue connected to at upper_conect.
	if ( other_res == 0 ) return 0; //Skip this residue if it's not connected to anything at its upper_connect.

	//Skip this pair if the connection is not polymer-like:
	if ( !pose.residue(other_res).type().is_polymer() ) return 0;
	if ( !pose.residue(other_res).has_lower_connect() ) return 0;
	if ( pose.residue(other_res).residue_connection_partner( pose.residue(other_res).lower_connect().index() ) != this_res_index ) return 0;

	//OK, all checks have passed at this point.  This res is connected normally to other_res.  Return the index of the other residue.
	return other_res;
}


} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::PolymerBondedEnergyContainer::PolymerBondedEnergyContainer() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::PolymerBondedEnergyContainer::save( Archive & arc ) const {
	arc( CEREAL_NVP( size_ ) ); // Size
	arc( CEREAL_NVP( score_types_ ) ); // utility::vector1<ScoreType>
	arc( CEREAL_NVP( tables_ ) ); // utility::vector1<utility::vector1<Real> >
	arc( CEREAL_NVP( computed_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( peptide_bonded_pair_indices_ ) ); // std::map<Size,Size>
	arc( CEREAL_NVP( inv_peptide_bonded_pair_indices_ ) ); // std::map<Size,Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::PolymerBondedEnergyContainer::load( Archive & arc ) {
	arc( size_ ); // Size
	arc( score_types_ ); // utility::vector1<ScoreType>
	arc( tables_ ); // utility::vector1<utility::vector1<Real> >
	arc( computed_ ); // utility::vector1<_Bool>
	arc( peptide_bonded_pair_indices_ ); // std::map<Size,Size>
	arc( inv_peptide_bonded_pair_indices_ ); // std::map<Size,Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::PolymerBondedEnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::PolymerBondedEnergyContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_PolymerBondedEnergyContainer )
#endif // SERIALIZATION
