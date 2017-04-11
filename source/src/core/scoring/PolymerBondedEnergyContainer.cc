// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

PolymerBondedNeighborIterator::~PolymerBondedNeighborIterator()= default;

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
		parent_->set_energy( base_, neighbor_id() , i , energy);
	}
}

void PolymerBondedNeighborIterator::retrieve_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] = parent_->get_energy( base_, neighbor_id() , i );
	}
}

void PolymerBondedNeighborIterator::accumulate_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] += parent_->get_energy( base_, neighbor_id() , i );
	}
}

void PolymerBondedNeighborIterator::mark_energy_computed() {
	parent_->set_computed( base_, neighbor_id() , true );
}

void PolymerBondedNeighborIterator::mark_energy_uncomputed() {
	parent_->set_computed( base_, neighbor_id(), false );
}

bool PolymerBondedNeighborIterator::energy_computed() const {
	return parent_->get_computed( base_, neighbor_id() );
}

/////////////////////////////////////////////////////

PolymerBondedNeighborConstIterator::~PolymerBondedNeighborConstIterator()= default;

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
		emap[ parent_->score_types()[i] ] = parent_->get_energy( base_, neighbor_id() , i );
	}
}

void PolymerBondedNeighborConstIterator::accumulate_energy( EnergyMap & emap ) const {
	for ( Size i=1; i<=parent_->score_types().size(); ++i ) {
		emap[ parent_->score_types()[i] ] += parent_->get_energy( base_, neighbor_id() , i );
	}
}

bool PolymerBondedNeighborConstIterator::energy_computed() const {
	return parent_->get_computed( base_, neighbor_id() );
}

/////////////////////////////////////////////////////////////////////////

/// @brief Destructor.
///
PolymerBondedEnergyContainer::~PolymerBondedEnergyContainer() = default;

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
	utility::vector1<core::Size> neighbors;
	typedef std::multimap<core::Size, core::Size>::const_iterator it;
	std::pair<it,it> range = chemical_edges_.equal_range(resid);

	for ( auto p = range.first; p != range.second; ++p ) {
		neighbors.push_back(p->second);
	}

	return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, neighbors, *this ) );
}

ResidueNeighborConstIteratorOP
PolymerBondedEnergyContainer::const_neighbor_iterator_end( int resid ) const {
	return ResidueNeighborConstIteratorOP( new PolymerBondedNeighborConstIterator( resid, utility::vector1<core::Size>(), *this ) );
}

ResidueNeighborConstIteratorOP
PolymerBondedEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const {
	utility::vector1<core::Size> neighbors;
	typedef std::multimap<core::Size, core::Size>::const_iterator it;
	std::pair<it,it> range = chemical_edges_.equal_range(resid);

	for ( auto p = range.first; p != range.second; ++p ) {
		// use residue ordering to decide lower vs upper
		if ( (int)p->second > resid ) {
			neighbors.push_back(p->second);
		}
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
	utility::vector1<core::Size> neighbors;
	typedef std::multimap<core::Size, core::Size>::const_iterator it;
	std::pair<it,it> range = chemical_edges_.equal_range(resid);

	for ( auto p = range.first; p != range.second; ++p ) {
		neighbors.push_back(p->second);
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
	utility::vector1<core::Size> neighbors;
	typedef std::multimap<core::Size, core::Size>::const_iterator it;
	std::pair<it,it> range = chemical_edges_.equal_range(resid);

	for ( auto p = range.first; p != range.second; ++p ) {
		// use residue ordering to decide lower vs upper
		if ( (int)p->second > resid ) {
			neighbors.push_back(p->second);
		}
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
	//fd (Nov 2016): Note that this container iterates over a set of edges independent of pose size
	//   So as long as the #edges stays fixed, we are fine
	//   Note that this means adding residues by jumps to a pose makes this still valid
	//   This is fine.

	// check that the correct connections have been set up.
	core::Size nres( pose.size() );
	core::Size connect_count(0);
	for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all residues (or all residues in the asymmetric unit, if this is a symmetric pose).
		if ( core::pose::symmetry::is_symmetric(pose) && !core::pose::symmetry::symmetry_info(pose)->bb_is_independent(ir) ) continue;

		core::conformation::Residue const & rsd = pose.residue(ir);
		core::Size nconnected = rsd.n_possible_residue_connections();

		typedef std::multimap<core::Size, core::Size>::const_iterator it;
		std::pair<it,it> range = chemical_edges_.equal_range(ir);

		for ( core::Size ic=1; ic<=nconnected; ++ic ) {
			if ( rsd.connected_residue_at_resconn( ic ) == 0 ) continue;
			core::Size connect_i = rsd.residue_connection_partner( ic );
			bool connect_found = false;

			// ensure ir,rsd.residue_connection_partner( ic ) is in the list
			for ( auto p = range.first; p != range.second && !connect_found; ++p ) {
				if ( p->second == connect_i ) connect_found = true;
			}

			if ( !connect_found ) {
				return false;
			}
			connect_count++;
		}
	}

	// every edge in the pose has been accounted for
	// ensure there are no additional edges in current container
	return ( connect_count == size_ );
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
	core::Size nres = pose.total_residue();
	chemical_edges_.clear();

	for ( core::Size ir=1; ir<=nres; ++ir ) {
		if ( core::pose::symmetry::is_symmetric(pose) && !core::pose::symmetry::symmetry_info(pose)->bb_is_independent(ir) ) continue;
		core::conformation::Residue const & rsd = pose.residue(ir);
		core::Size nconnected = rsd.n_possible_residue_connections();
		for ( core::Size ic=1; ic<=nconnected; ++ic ) {
			if ( rsd.connected_residue_at_resconn( ic ) == 0 ) continue;
			std::pair<core::Size,core::Size> edge = std::make_pair(ir,rsd.residue_connection_partner( ic ));
			chemical_edges_.insert(edge);
			if ( edge.first < edge.second ) {
				tables_[edge] = utility::vector1<core::Real>(score_types_.size(),0.0);
				computed_[edge] = false;
			}
		}
	}

	size_ = chemical_edges_.size();
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
	arc( CEREAL_NVP( chemical_edges_ ) ); // std::multimap<Size,Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::PolymerBondedEnergyContainer::load( Archive & arc ) {
	arc( size_ ); // Size
	arc( score_types_ ); // utility::vector1<ScoreType>
	arc( tables_ ); // utility::vector1<utility::vector1<Real> >
	arc( computed_ ); // utility::vector1<_Bool>
	arc( chemical_edges_ ); // std::multimap<Size,Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::PolymerBondedEnergyContainer );
CEREAL_REGISTER_TYPE( core::scoring::PolymerBondedEnergyContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_PolymerBondedEnergyContainer )
#endif // SERIALIZATION
