// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh
/// @brief  Declaration for a class to hold a set of rotamers for a flexible packing run
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu), oct 08

#ifndef INCLUDED_protocols_flexpack_rotamer_set_FlexbbRotamerSets_hh
#define INCLUDED_protocols_flexpack_rotamer_set_FlexbbRotamerSets_hh

/// Unit headers
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.fwd.hh>

/// Package headers
#ifdef WIN32
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.hh>
#endif

#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.fwd.hh>
#include <protocols/flexpack/interaction_graph/OTFFlexbbInteractionGraph.fwd.hh>
#include <protocols/flexpack/interaction_graph/PrecomputedFlexbbInteractionGraph.fwd.hh>

/// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/types.hh>

/// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace flexpack {
namespace rotamer_set {

class FlexbbRotamerSets : public core::pack::rotamer_set::RotamerSetsBase, public utility::pointer::enable_shared_from_this< FlexbbRotamerSets >
{
public:
	typedef core::Size Size;
	typedef core::PackerEnergy PackerEnergy;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::methods::LongRangeTwoBodyEnergy LongRangeTwoBodyEnergy;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;

public:
	FlexbbRotamerSets( PackerTaskCOP task );

	virtual ~FlexbbRotamerSets();

	inline FlexbbRotamerSetsCOP get_self_ptr() const { return shared_from_this(); }
	inline FlexbbRotamerSetsOP  get_self_ptr() { return shared_from_this(); }
	inline FlexbbRotamerSetsCAP get_self_weak_ptr() const { return FlexbbRotamerSetsCAP( shared_from_this() ); }
	inline FlexbbRotamerSetsAP  get_self_weak_ptr() { return FlexbbRotamerSetsAP( shared_from_this() ); }

	void
	set_frames(
		Pose const & pose,
		utility::vector1< core::fragment::FrameCOP > const & frames
	);

	utility::graph::GraphOP
	flexpack_neighbor_graph(
		Pose const & pose,
		ScoreFunction const & sfxn,
		PackerTaskCOP task
	) const;

	core::Distance
	determine_res_cb_deviation(
		Pose const & pose,
		Size resid
	) const;

	void build_rotamers(
		Pose const & pose,
		ScoreFunction const & sfxn,
		utility::graph::Graph const & flexpack_neighbor_graph
	);

	void
	dump_pdbs( core::pose::Pose const & pose, std::string const & filename_base ) const;

	void
	precompute_energies(
		Pose const & pose,
		ScoreFunction const & sfxn,
		utility::graph::GraphCOP flexpack_neighbor_graph,
		interaction_graph::FlexbbInteractionGraph & flexbb_ig
	) const;

	void precompute_all_energies(
		Pose const & pose,
		ScoreFunction const & sfxn,
		utility::graph::GraphCOP flexpack_neighbor_graph,
		interaction_graph::PrecomputedFlexbbInteractionGraph & flexbb_ig
	) const;

	/* void compute_one_body_energies_for_precomputed_ig(
	Pose const & pose,
	ScoreFunction const & sfxn,
	utility::graph::GraphCOP flexpack_neighbor_graph,
	interaction_graph::FlexbbInteractionGraph & flexbb_ig
	) const; */

	void compute_one_body_energies_for_otf_ig(
		Pose const & pose,
		ScoreFunction const & sfxn,
		utility::graph::GraphCOP flexpack_neighbor_graph,
		interaction_graph::OTFFlexbbInteractionGraph & flexbb_ig
	) const;

	/* void compute_two_body_energies(
	Pose const & pose,
	ScoreFunction const & sfxn,
	utility::graph::Graph const & flexpack_neighbor_graph,
	interaction_graph::PrecomputedFlexbbInteractionGraph & flexbb_ig
	) const; */

	//void prepare_sets_for_packing( Pose const & pose, ScoreFunction const &);

	/// Virtual functions from base class
	virtual core::uint nrotamers() const;
	virtual core::uint nrotamers_for_moltenres( core::uint ) const;

	virtual core::uint nmoltenres() const;

	virtual core::uint total_residue() const;

	virtual
	core::uint
	moltenres_2_resid( core::uint ) const;

	virtual
	core::uint
	resid_2_moltenres( core::uint ) const;

	virtual
	core::uint
	moltenres_for_rotamer( core::uint ) const;

	virtual
	core::uint
	res_for_rotamer( core::uint ) const;

	virtual
	core::conformation::ResidueCOP
	rotamer( core::uint ) const;

	virtual
	core::conformation::ResidueCOP
	rotamer_for_moltenres( core::uint moltenres_id, core::uint rotamerid ) const;

	virtual
	core::uint
	nrotamer_offset_for_moltenres( core::uint ) const;

	/// @brief convert rotid in full rotamer enumeration into rotamer id on its source residue
	virtual
	core::uint
	rotid_on_moltenresidue( core::uint rotid ) const;

	/// @brief convert moltenres rotid to id in full rotamer enumeration
	virtual
	core::uint
	moltenres_rotid_2_rotid( core::uint moltenres, core::uint moltenresrotid ) const;

	Size
	nmolten_res() const
	{ return nmoltenres_;}

	//Size
	//moltenres_2_resid( Size molten_resid ) const
	//{ return moltenres_2_resid_[ molten_resid ]; }

	//Size
	//resid_2_moltenres( Size resid ) const
	//{ return resid_2_moltenres_[ resid ]; }

	/// @brief The total number of rotamers across all residues and all backbone conformations
	//Size
	//nrotamers() const {
	// return nrotamers_;
	//}

	//Size
	//nrotamers_for_moltenres( Size molten_resid ) const
	//{ return nrotamers_for_moltenres_[ molten_resid ]; }

	Size
	nrotamers_for_res( Size resid ) const
	{ return nrotamers_for_moltenres_[ resid_2_moltenres_[ resid ] ]; }

	Size
	nbbconfs_for_moltenres( core::Size moltenres ) const;

	Size
	nbbconfs_for_flexseg( core::Size flex_segment_id ) const
	{ return nbbconfs_for_flexseg_[ flex_segment_id ]; }

	/// @brief the "total" number of backbone conformations, counting the input conformation
	/// for all flexible segments (but not counting the input conformation for residues with only a single
	/// backbone conformation)
	Size
	nbackbone_conformations() const;

	Size
	nbbconfs_for_res( core::Size resid ) const;

	utility::vector1< Size > const &
	num_states_per_backbone_for_moltenres( Size moltenres ) const {
		return nrots_for_moltenres_bbconf_[ moltenres ];
	}

	Size
	nrotamers_for_moltenres_in_bbconf( Size moltenres, Size bbconf ) const
	{ return nrots_for_moltenres_bbconf_[moltenres][bbconf]; }

	Size
	nrotamers_for_resid_in_bbconf( Size resid, Size bbconf ) const
	{ return nrots_for_moltenres_bbconf_[ resid_2_moltenres_[ resid ] ][ bbconf ]; }

	Size
	global_rotid_start_for_moltenres_in_bbconf( Size moltenres, Size bbconf ) const
	{ return nrotoffset_for_moltenres_bbconf_[ moltenres ][ bbconf ] ;}

	Size
	local_rotid_start_for_moltenres_in_bbconf( Size moltenres, Size bbconf ) const
	{ return nrotoffset_for_moltenres_bbconf_[ moltenres ][ bbconf ] - nrotoffset_for_moltenres_[ moltenres ]; }

	//virtual
	//uint
	//nrotamer_offset_for_moltenres( uint moltenresid ) const {
	// return nrotoffset_for_moltenres_[ moltenresid ];
	//}

	/// @brief Input: rotid in global enumeration, Output rotid in local enumeration
	Size
	local_rotid_for_rotamer_on_moltenres( Size moltenres, Size rotamer_id ) const
	{ return rotamer_id - nrotoffset_for_moltenres_[ moltenres ] ; }


	Size
	local_rotid_for_rotamer( Size rotamer_id ) const
	{ return rotamer_id - nrotoffset_for_moltenres_[ moltenres_for_rotamer_[ rotamer_id ] ]; }

	Size
	nflexible_segments() const
	{ return flexsegment_span_.size(); }

	Size
	flexsegment_start_moltenresid( Size flex_segment_id ) const
	{ return resid_2_moltenres_[ flexsegment_span_[ flex_segment_id ].first ]; }

	Size
	flexsegment_stop_moltenresid( Size flex_segment_id ) const
	{ return resid_2_moltenres_[ flexsegment_span_[ flex_segment_id ].second ]; }

	Size
	flexsegment_size( Size flex_segment_id ) const
	{ return flexsegment_span_[ flex_segment_id ].second - flexsegment_span_[ flex_segment_id ].first + 1; }

	Size
	flexsegment_start_resid( Size flex_segment_id ) const
	{ return flexsegment_span_[ flex_segment_id ].first; }

	Size
	flexsegment_stop_resid( Size flex_segment_id ) const
	{ return flexsegment_span_[ flex_segment_id ].second; }

	Size
	flexsegid_for_moltenres( Size molten_resid ) const
	{ return moltenres_2_flexseg_[ molten_resid ]; }

	Size
	flexsegid_for_res( Size resid ) const
	{ return moltenres_2_flexseg_[ resid_2_moltenres_[ resid ] ]; }

	bool
	moltenres_part_of_flexsegment( Size molten_resid ) const
	{ return moltenres_2_flexseg_[ molten_resid ] != 0; }

	bool
	res_part_of_flexsegment( Size resid ) const
	{ return moltenres_2_flexseg_[ resid_2_moltenres_[ resid ] ] != 0; }

	FlexbbRotamerSetCOP
	rotset_for_moltenres( Size molten_resid, Size bbconf = 1 ) const;

	FlexbbRotamerSetCOP
	rotset_for_residue( Size resid, Size bbconf = 1 ) const;

	/// @brief Rotamer indexed locally
	core::conformation::ResidueCOP
	rotamer_for_residue( Size resid, Size rotindex_on_residue ) const;

	/// @brief Rotamer indexed locally
	//core::conformation::ResidueCOP
	//rotamer_for_moltenres( Size moltenres, Size rotindex_on_residue ) const;

	/// @brief Rotamer indexed globally.
	//core::conformation::ResidueCOP
	//rotamer( Size rotindex ) const;

	core::conformation::Residue const &
	backbone_for_resid_bbconf( Size resid, Size bbconf ) const
	{ return *conformations_for_flexible_segments_[ resid_2_moltenres_[ resid ] ][ bbconf ]; }

	core::conformation::Residue const &
	backbone_for_moltenres_bbconf( Size moltenres, Size bbconf ) const
	{ return *conformations_for_flexible_segments_[ moltenres ][ bbconf ]; }


	/// @brief rotamer_id in the global enumeration of rotamers.
	//Size
	//moltenres_for_rotamer( Size rotamer_id ) const {
	// return moltenres_for_rotamer_[ rotamer_id ];
	//}

	Size
	bbconf_for_rotamer( Size rotamer_id ) const {
		Size moltenres = moltenres_for_rotamer_[ rotamer_id ];
		return bbconf_for_rotamer_of_moltenres_[ moltenres ][ local_rotid_for_rotamer_on_moltenres( moltenres, rotamer_id ) ];
	}

	Size
	bbconf_for_rotamer_on_moltenres( Size moltenres, Size rotamer_id ) const {
		return bbconf_for_rotamer_of_moltenres_[ moltenres ][ rotamer_id ];
	}

	/// @brief Do two rotamers originate from the same backbone conformation?
	/// Local enumeration for both rot1 and rot2; rot1 and rot2 must be rotamers of the
	/// same molten residue.
	bool
	rotamers_on_same_bbconf( Size moltenres, Size rot1, Size rot2 ) const {
		utility::vector1< Size > const & bbconf( bbconf_for_rotamer_of_moltenres_[ moltenres ] );
		return bbconf[ rot1 ] == bbconf[ rot2 ];
	}

	virtual
	void
	initialize_pose_for_rotsets_creation(
		core::pose::Pose & /*pose*/
	) const {}

	virtual
	void
	show( std::ostream & out ) const;

protected:

	void
	build_residue_vector_from_fragment(
		Pose & pose,
		core::fragment::FrameCOP frame,
		Size frag_num,
		utility::vector1< core::conformation::ResidueOP > & fragment_res
	);

private:
	void update_offset_data();

	void
	compute_sr_one_body_energies_for_flexsets(
		Size lowermoltenres,
		Size uppermoltenres,
		utility::vector1< utility::vector1< Size > > const & regular_representatives,
		utility::vector1< utility::vector1< Size > > const & proline_representatives,
		utility::vector1< utility::vector1< Size > > const & glycine_representatives,
		Pose const & pose,
		ScoreFunction const & sfxn,
		interaction_graph::OTFFlexbbInteractionGraph & flexbb_ig
	) const;

	void
	compute_onebody_interactions_with_background(
		Size moltenres,
		Pose const & pose,
		ScoreFunction const & sfxn,
		utility::graph::GraphCOP flexpack_neighbor_graph,
		interaction_graph::FlexbbInteractionGraph & flexbb_ig
	) const;

private:

	//not derived from fixbb RotamerSets, so we need to copy some of the basic data structs
	Size nmoltenres_;
	Size total_residue_;
	Size nbbconfs_;

	PackerTaskCOP task_;  // initialized at construction

	utility::vector1< utility::vector1< rotamer_set::FlexbbRotamerSetOP > > rotamers_;
	Size nrotamers_;
	utility::vector1< Size > nrotamers_for_moltenres_;
	utility::vector1< Size > moltenres_2_resid_;
	utility::vector1< Size > resid_2_moltenres_;

	//andrew, did you forget to put these two vectors in?
	utility::vector1< Size > moltenres_for_rotamer_;
	utility::vector1< utility::vector1< Size > > bbconf_for_rotamer_of_moltenres_;

	//why is this a vector of vector of residues?
	utility::vector1< utility::vector1< core::conformation::ResidueCOP > > conformations_for_flexible_segments_;
	utility::vector1< std::pair< Size, Size > > flexsegment_span_;
	utility::vector1< Size > nbbconfs_for_flexseg_;
	utility::vector1< Size > moltenres_2_flexseg_;
	utility::vector1< Size > nrotoffset_for_moltenres_;
	utility::vector1< utility::vector1< Size > > nrots_for_moltenres_bbconf_;
	utility::vector1< utility::vector1< Size > > nrotoffset_for_moltenres_bbconf_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	FlexbbRotamerSets();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_flexpack_rotamer_set_FlexbbRotamerSets )
#endif // SERIALIZATION


#endif
