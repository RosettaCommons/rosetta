// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author John Karanicolas, Jacob Corn (jecorn@u.washington.edu), Sarel Fleishman


#ifndef INCLUDED_protocols_hotspot_hashing_HotspotStubSet_hh
#define INCLUDED_protocols_hotspot_hashing_HotspotStubSet_hh

#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/exit.hh>

// C++ Headers
#include <map>
#include <vector>
#include <set>

#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>


#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh> // WIN32 INCLUDE
#endif

namespace protocols {
namespace hotspot_hashing {
typedef platform::Size Size;

class HotspotStubSet : public utility::pointer::ReferenceCount {
public:
	typedef std::multimap< core::Real, HotspotStubOP > Hotspots;
	typedef std::map< std::string, Hotspots > Hs_map;
	typedef std::pair< std::string, std::pair< core::Real, HotspotStubOP > > Hs_data;
	// for iterator type access of stubs
	typedef std::vector< Hs_data > Hs_vec;
	typedef Hs_vec::iterator iterator;
	typedef Hs_vec::const_iterator const_iterator;
public:

	HotspotStubSet();
	HotspotStubSet( HotspotStubSet const & init );
	virtual ~HotspotStubSet();
	void clear();

	//iterator functions
	// iterators
	inline HotspotStubSet::const_iterator begin() const {
		runtime_assert( stub_set_vec_.size() == size() );
		return stub_set_vec_.begin();
	}
	inline HotspotStubSet::const_iterator end() const { return stub_set_vec_.end(); }
	inline HotspotStubSet::iterator begin(){
		runtime_assert( stub_set_vec_.size() == size() );
		return stub_set_vec_.begin();
	}
	inline HotspotStubSet::iterator end(){ return stub_set_vec_.end(); }

	/// @brief set length of hotspot peptide used for finding stubs. Default=1, Hotspot stored in set will always be only 1 residue
	Size hotspot_length( ) const;
	void hotspot_length( Size const length );

	/// @brief returns a new stub_set with stub scores recalculated by colony energy (Xiang, Soto, and Honig PNAS 2002)
	HotspotStubSetOP colonyE( );
	/// @brief cluster stubs, returning a new stubset of cluster centers for each residue
	HotspotStubSetOP cluster( );

	/// @brief retrieve all stubs with a given residue name
	Hotspots retrieve( std::string const & residue_name3 );

	/// @brief gets a stub from the stub_set_. A version returning an OP is private
	Hotspots::const_iterator get_stub( std::string const & residue_name3, core::Real const score ) const;
	// @brief build a new stubset containing stubs with a given residue name3 and score cutoff
	HotspotStubSetOP subset( std::string const & residue_name3, core::Real const scorecut );
	/// @brief build a new stubset containing stubs that pass a score cutoff
	HotspotStubSetOP subset( core::Real const scorecut ) const;

	/// @brief gets the best energy stub in the set
	HotspotStubCOP get_best_energy_stub() const;
	/// @brief get the stub that's nearest to a residue
	HotspotStubCOP get_nearest_stub( core::conformation::ResidueCOP stub ) const;
	/// @brief finds neighbors to stub based on distance between atoms.
	std::set< std::pair< std::string, core::Real > > find_neighboring_stubs( HotspotStubCOP stub ) const;
	/// @brief removes neighbors of stub based on repulsive energy between the pair of residues
	void remove_stubs_from_set( std::set< std::pair< std::string, core::Real > > const );
	void remove_random_stubs_from_set( int const num_to_remove );
	/// @brief removes a single stub. Reports false if stub is not found
	bool remove_stub( HotspotStubCOP stub );

	/// @brief add to stubset by reading from a file
	void read_data( std::string const & filename );
	// unfortunately, this won't compile on Windows with the BOINC libraries. Ask
	// tex for more information.
	//void read( std::string const filename );

	/// @brief fill the stub set with n_stubs by Rosetta residue name
	void fill( core::pose::Pose const & reference_pose, core::scoring::ScoreFunctionCOP scorefxn_in, std::string const & residue_name3, Size const n_stubs );
	/// @brief only keep stubs within a certain distance of a residue on the target pose.
	void fill( core::pose::Pose const & reference_pose, core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const target, core::Real const distance, std::string const & residue_name3, Size const n_stubs );

	/// @brief rescore all stubs in this set based on current flags (eg - sc_only() )
	HotspotStubSetOP rescore( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn );


	/// @brief fill the stub set with n_stubs each of A, R, N, D, E, Q, H, I, L, K, M, P, F, S, T, W, Y, and V
	void autofill ( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn, Size const n_stubs );
	void autofill( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn, core::Size const target, core::Real const distance, Size const n_stubs );

	/// @brief only accept stubs that score better than this threshold (default is -1.0)
	void score_threshold( core::Real const threshold );

	/// @brief write all stubs contained in this set to filename
	void write_all( std::string const & filename ) const;
	/// @brief write one stub with a user-supplied tag
	void write_stub( utility::io::ozstream & outstream, HotspotStubCOP stub, std::string const & tag ) const;

	/// @brief associate all stubs in the set with a scaffold partner
	//SJF does it make sense to associate the entire stubset with a filter? The filter is going to change ALL the time.
	void pair_with_scaffold( core::pose::Pose const & pose, core::Size const partner, protocols::filters::FilterCOP filter  ) ;
	/// @brief set the filter to use for scaffold matching within this set
	void filter( protocols::filters::FilterCOP filter );

	/// @brief true if these stubs are sidechain-only (defaults true)
	bool sc_only() const;

	/// @brief set whether or not sidechains are included
	void sc_only( bool sc_switch );


	/// @brief how many total stubs are in the set (all residues)?
	core::Size size() const;
	/// @brief how many stubs are in the set by residue?
	core::Size size( std::string const & resname );

	/// @brief returns a random stub either from the entire set or based on residue name
	HotspotStubOP random_stub();
	HotspotStubOP random_stub( std::string const & resname );

	/// @brief Sets up constraints using a given partner (auto choose packer task and fixed reference atom)
	void add_hotspot_constraints_to_pose(
		core::pose::Pose & pose,
		core::Size const partner,
		HotspotStubSetOP hotspot_stub_set,
		core::Real const & CB_force_constant,
		core::Real const & worst_allowed_stub_bonus,
		bool const apply_self_energies,
		core::Real const & bump_cutoff,
		bool const apply_ambiguous_constraints = false
	);

	void add_hotspot_constraints_to_wholepose(
		core::pose::Pose & pose,
		core::Size const partner,
		HotspotStubSetOP hotspot_stub_set,
		core::Real const & CB_force_constant,
		core::Real const & worst_allowed_stub_bonus,
		bool const apply_self_energies,
		core::Real const & bump_cutoff,
		bool const apply_ambiguous_constraints = false
	);

	/// @brief Sets up constraints with user-supplied packer task and fixed reference atom
	void add_hotspot_constraints_to_pose(
		core::pose::Pose & pose,
		core::id::AtomID const & fixed_atom,
		core::pack::task::PackerTaskCOP const packer_task,
		HotspotStubSetOP hotspot_stub_set,
		core::Real const & CB_force_constant,
		core::Real const & worst_allowed_stub_bonus, // = 0.
		bool const apply_self_energies,
		core::Real const & bump_cutoff,
		bool const apply_ambiguous_constraints = false
	);

	/// @brief Sets up constraints with user-supplied packer task and fixed reference atom
	void add_hotspot_constraints_to_wholepose(
		core::pose::Pose & pose,
		core::id::AtomID const & fixed_atom,
		core::pack::task::PackerTaskCOP const packer_task,
		HotspotStubSetOP hotspot_stub_set,
		core::Real const & CB_force_constant,
		core::Real const & worst_allowed_stub_bonus, // = 0.
		bool const apply_self_energies,
		core::Real const & bump_cutoff,
		bool const apply_ambiguous_constraints = false
	);

	/// @brief remove all ambiguous constraints that contain backbone_stub_constraints from the supplied pose
	bool remove_all_hotspot_constraints( core::pose::Pose & pose ) const;
	void set_chain_to_design( core::Size const chain_to_design = 2 );
	// Constraint creation methods
	core::pack::task::PackerTaskOP prepare_hashing_packer_task_(
		core::pose::Pose const & pose,
		core::Size const chain_to_redesign = 2
	);

	/// @brief return bbcst's associatied with this stub set
	core::scoring::constraints::ConstraintCOPs constraints() const;

	void add_stub_( HotspotStubCOP stub );
	void add_stub_set( HotspotStubSet const & stubset );
private:
	friend class HotspotStub;
	// Stub data storage
	// This was previously a map of multisets, but multisets only allow const iteration (even if you use a normal iterator)
	// The use of a map of multimaps allows us to make changes to the stub_status, but somewhat obscures stub iteration
	// So when iterating over a retrieved stubset, be sure you get the second element
	Hs_map stub_set_;
	Hs_vec stub_set_vec_; // A mirror of the data in stub_set_ for easier iterator access
	bool sc_only_;
	core::Size target_resnum_; // resnum to use for focused hashing
	core::Real target_distance_; // radius from resnum for focused hashing

	core::Size chain_to_design_;
	core::pose::PoseCOP pose_;
	core::Real score_threshold_; // only accept stubs with this score or better
	protocols::filters::FilterCOP filter_;
	Size hotspot_length_; // length of peptide to use for hotspot searching (polyAla, except for central hotspot). only hotspot itself is scored/stored in the set.


	/// @brief clears stub_set_vec_ and inserts all the elements in stub_set_ to it.
	void handshake_stub_sets( void );

	core::scoring::constraints::ConstraintCOPs constraints_;


	// predicate for keeping the multiset sorted by bonus value
	// obsolete, since we're now using a multimap
	/* struct stubsort_pred_
	{
	bool operator () ( HotspotStub const & left, HotspotStub const & right )
	{
	return left.bonus_value() < right.bonus_value();
	}
	};
	*/
	// Stub creation methods
	//void dock_residue_lowres_ ( core::pose::Pose & pose, platform::Size const jump_number ) ;
	//void dock_residue_highres_ ( core::pose::Pose & pose, core::scoring::ScoreFunctionOP scorefxn, Size const jump_number ) ;
	void create_hotspot_after_pose ( core::pose::Pose & pose, std::string const & resname ) ;
	void setup_hotspot_foldtree_ ( core::pose::Pose & pose ) ;
	/// @brief utility function to find distance of stub from the end of the pose.
	core::Size stub_offset();

	core::Real get_residue_score_ ( core::pose::Pose const & pose, core::scoring::ScoreFunctionCOP scorefxn, Size const seqpos) ;

	core::Real evaluate_stub_bumps_(
		core::conformation::Residue const & placed_stub_residue,
		core::pose::Pose const & unbound_pose,
		core::Size const resnum,
		core::scoring::TenANeighborGraph const & unbound_neighbor_graph,
		core::scoring::ScoreFunctionCOP const & bump_scorefxn,
		core::Real const & max_bump_energy
	);

	core::Real evaluate_stub_self_energy_(
		core::conformation::Residue const & stub_residue,
		core::pose::Pose const & unbound_pose,
		core::Size const resnum,
		core::scoring::TenANeighborGraph const & unbound_neighbor_graph,
		core::scoring::ScoreFunctionCOP const & full_scorefxn
	);

	void generate_unbound_pose_( core::pose::Pose & pose );

	/// @brief gets a stub from the stub_set_
};

/// @brief utility function for deleting all backbone stub constraints from a pose.
/// Returns the removed constraints (ambiguous).
core::scoring::constraints::ConstraintCOPs
remove_hotspot_constraints_from_pose( core::pose::Pose & );

/// @brief utility function to calculate per-residue sidechain rmsd without superposition
core::Real residue_sc_rmsd_no_super( core::conformation::ResidueCOP res1, core::conformation::ResidueCOP res2, bool const fxnal_group_only=false );

/// @brief utility function to make sure stub's Cbeta is not pointed away from the target.
core::Real stub_tgt_angle( core::pose::Pose const & pose, core::conformation::ResidueCOP stub, core::Size const target_res );


} // namespace hotspot_hashing
} // namespace protocols


#endif
