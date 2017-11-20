// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerSet/RotamerSet.hh
/// @brief  Residue Set class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSet_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSet_hh

// Unit Headers
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/Residue.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
#include <core/scoring/trie/RotamerTrieBase.hh>
#include <core/conformation/Residue.hh> // WIN32 INCLUDE
#endif

// Basic headers
#include <basic/datacache/BasicDataCache.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

typedef std::pair<core::Angle, core::Angle> RotamerBin;
typedef utility::vector1<RotamerBin> RotamerBins;


class RotamerSet : public conformation::RotamerSetBase
{

private: // typedefs
	typedef conformation::RotamerSetBase parent;

public:
	RotamerSet();
	virtual ~RotamerSet();

	void set_resid( Size resid );

	virtual
	void build_rotamers(
		pose::Pose const & the_pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	) = 0;

	virtual
	void build_dependent_rotamers(
		RotamerSets const & rotamer_sets,
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph
	) = 0;

	/// @brief Append a rotamer to the list; it will not be sorted into the same group as other rotamers
	/// of the same group or amino acid unless the last group/amino acid is already the same. Instead it
	/// will sit at the end of the list of Rotamers. There are performance implications of using this
	/// function instead of add_rotamer_into_existing_group: there are several places in the code which
	/// scale quadratically with the number of amino acid groups (but where we assume that this number is small)
	/// so if you call this function N times oscilating between ASP and ASN rotamers, you will get
	/// O(N^2) performance. If you do not need your rotamers to appear in a particular order,
	/// use add_rotamer_into_existing_group instead.
	virtual
	void
	add_rotamer(
		conformation::Residue const & rotamer
	) = 0;

	/// @brief Add a rotamer to the RotamerSet where you will group it with other residues of the same type
	/// or barring that, the same group. This will keep the total number of residue type groups down. It will
	/// not guarantee (it cannot guarantee) that the newly added rotamer will appear after existing rotamers
	/// or at the end of the list of rotamers -- if you need that kind of guarantee, use add_rotamer instead.
	virtual
	void
	add_rotamer_into_existing_group(
		conformation::Residue const & rotamer
	) = 0;

	/// @brief Return the number of different residue types; two residue types are considered
	/// different if they have a different address.
	virtual
	Size
	get_n_residue_types() const = 0;

	/// @brief Return the number of different residue groups.  Two residue types are considered
	/// to be part of the same block of residues if 1. they have the same address or 2. they have
	/// the same "name3" and the same neighbor radius.
	virtual
	Size
	get_n_residue_groups() const = 0;

	/// @brief Return the first rotamer of a particular residue type
	virtual
	Size
	get_residue_type_begin( Size which_restype ) const = 0;

	/// @brief Return the first rotamer that belongs to a particular rotamer group
	virtual
	Size
	get_residue_group_begin( Size which_resgroup ) const = 0;

	virtual
	Size
	get_n_rotamers_for_residue_type( Size which_restype ) const = 0;

	virtual
	Size
	get_n_rotamers_for_residue_group( Size which_resgroup ) const = 0;

	/// @brief Rotamers i to i+j of all the same residue type are grouped together.
	/// This function returns the index of the residue type in a contiguous block
	/// of rotamers.  E.g. rotamers 100 to 120 might all be lysine rotamers, and might
	/// be the 8th residue type, with the first 7 residue types spanning rotamers 1 to 99.
	/// If new lysine rotamers are appended to the end of the rotamer set, they are
	/// considered to be in a separate residue type block.  Lysine rotamers 200 to 210 might
	/// be block 15 while lysine rotamers 100 to 120 are still block 7.
	virtual
	Size
	get_residue_type_index_for_rotamer( Size which_rotamer ) const = 0;

	/// @brief Return the index of the rotamer group for a particular rotamer.
	virtual
	Size
	get_residue_group_index_for_rotamer( Size which_rotamer ) const = 0;

	virtual
	void
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & energies ) const = 0;

	virtual
	void
	compute_one_and_two_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & one_body_energies,
		utility::vector1< utility::vector1< core::PackerEnergy > > & two_body_energies,
		utility::vector1< core::Size > & packable_neighbors ) const = 0;

	/// for OptE
	virtual
	void
	compute_one_body_energy_maps(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		utility::vector1< scoring::EnergyMap > & energies ) const = 0;

	virtual
	Size
	num_rotamers() const = 0;

	/// @brief Return the index in the RotamerSet for the current rotamer
	virtual
	Size
	id_for_current_rotamer() const = 0;

	virtual
	Size resid() const { return resid_;}

	virtual
	conformation::ResidueCOP
	rotamer( Size rot_id ) const = 0;

	virtual
	basic::datacache::BasicDataCache &
	rotamer_data_cache( Size rot_id ) const = 0;

	// iterate over rotamers directly
	virtual Rotamers::const_iterator begin() const = 0;
	virtual Rotamers::const_iterator end() const = 0;

	/// @brief mutatable access to a single rotamer in the set.
	virtual
	conformation::ResidueOP
	nonconst_rotamer( Size rot_id ) = 0;

	virtual
	void
	store_trie( Size method_enum_id, conformation::AbstractRotamerTrieOP trie ) = 0;

	virtual
	conformation::AbstractRotamerTrieCOP
	get_trie( Size method_enum_id ) const = 0;

	/// @brief removes a single rotamer and causes a rotamer index update
	virtual
	void
	drop_rotamer( Size rot_id ) = 0;

	/// @brief rotamers_to_delete must be of size nrotmaers -- each position
	/// in the array that's "true" is removed from the set of rotamers
	virtual
	void
	drop_rotamers( utility::vector1< bool > const & rotamers_to_delete ) = 0;

	/// @brief deletes the rotamers in the list with the given indices.
	/// The indices of these rotamers is presumed to be those before any delete operation.
	/// e.g. if there are four rotamers, and rotamer_indices_to_delete includes 1 & 3,
	/// then the rotamers that will remain are the rotamers originally indexed as 2 and 4,
	/// even though their new indices will be 1 & 2.
	virtual
	void
	drop_rotamers_by_index( utility::vector1< Size > const & rotamer_indices_to_delete ) = 0;

	virtual
	void
	show( std::ostream & out ) const = 0;

	virtual
	void
	initialize_pose_for_rotset_creation(
		pose::Pose & pose
	) const = 0;

private:
	// deny use of the copy constructor (no pass-by-value)
	RotamerSet( RotamerSet const & );

	Size resid_; //which residue is this?

	/// @brief BasicDataCache indexed by enum in core/pack/rotamer_set/RotamerSetCacheableDataType.hh
	/// @warning DataCache must always be initialized with the number of cacheable
	///  data types -- see the last enum entry.
	//BasicDataCache data_cache_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

inline
std::ostream &
operator<<( std::ostream & out, RotamerSet const & rs ) {
	rs.show(out);
	return out;
}

} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSet )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSet_HH
