// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/connection/Connection.hh
/// @brief Connection classes for building structures from components
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_connection_Connection_hh
#define INCLUDED_protocols_denovo_design_connection_Connection_hh

// Unit headers
#include <protocols/denovo_design/connection/Connection.fwd.hh>
#include <protocols/denovo_design/components/NamedMover.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>

// Package headers
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.fwd.hh>

// Core headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ headers
#include <string>
#include <set>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace denovo_design {
namespace connection {

class Connection : public components::NamedMover {
	// sub-objects
public:
	struct Motif {
		Motif() :
			len( 0 ),
			ss( "" ),
			abego( "" ) {}

		Motif( core::Size const len_val, char const ss_val, std::string const & abego_val ) :
			len( 0 ),
			ss( "" ),
			abego( "" )
		{
			add( len_val, ss_val, abego_val );
		}

		Motif( core::Size const len_val, char const ss_val, char const abego_val ) :
			len( 0 ),
			ss( "" ),
			abego( "" )
		{
			add( len_val, ss_val, abego_val );
		}

		void add( core::Size const len_val, char const ss_val, std::string const & abego_val )
		{
			len += len_val;
			for ( core::Size i=1; i<=len_val; ++i ) {
				ss += ss_val;
				abego += abego_val;
			}
		}

		void add( core::Size const len_val, char const ss_val, char const abego_val )
		{
			len += len_val;
			for ( core::Size i=1; i<=len_val; ++i ) {
				ss += ss_val;
				abego += abego_val;
			}
		}

		friend std::ostream & operator<<( std::ostream & os, Motif const & motif );
		core::Size len;
		std::string ss;
		std::string abego;
	};
	typedef utility::vector1< Motif > MotifList;

public:
	Connection();
	virtual ~Connection();

	/// @brief setup the parameters via an xml tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief performs setup and applies loop building
	/// @details Steps to apply():
	/// 1. Pulls StructureData from the pose
	/// 2. setup_permutation() stores data about this connection which is not
	///    static for every apply() call
	/// 3. apply_permutation() uses the StructureData object to build the loop
	/// 4. check() checks the built loop
	virtual void apply( core::pose::Pose & pose );

	/// @brief applies the loop building
	/// @throw EXCN_ConnectionFailed on failure to connect
	/// @details You should overload this and insert your loop-building method here.
	/// The input StructureData object contains a pose which has the loop residues
	/// built in extended conformation with a cutpoint in a user-specified location
	/// (or random location if the user hasn't specified anything)
	/// The expected result of this function is a connected loop with sealed cut.
	/// Cutpoints/fold tree can be easily sealed by calling:
	///
	/// perm.delete_jump_and_intervening_cutpoint( loop_lower(perm), loop_upper(perm) );
	///
	virtual void apply_connection( components::StructureData & perm ) const = 0;

	/// @brief uses the permutation to set up, builds, and updates permutation
	/// @details Steps are:
	///  1. Collect information about the pose/connection
	///  2. Move residues such that residues to be joined are
	///     adjacent in sequence.
	///  3. Removes terminal variants if necessary
	///  4. Adds cutpoint variants if necessary
	///  5. Declares covalent bond to the conformation if necessary
	///  6. Adds necessary residues to the pose in extended conformation with
	///     cutpoint at cut_resi( perm )
	///  7. Calls pure virtual apply_connection( perm )
	///  8. If apply_connection fails, stop and set status to
	///     protocols::moves::FAIL_RETRY
	///  9. Performs bookkeeping on the StuctureData object
	/// 10. Calls check() and sets status to protocols::moves::FAIL_RETRY
	///     if the check fails.
	/// @throws EXCN_ConnectionFailed() on connection failure
	void apply_permutation( components::StructureData & perm ) const;

	/// @brief set up the connection mover based on the information in the permutation
	/// @throw EXCN_Setup if no valid free connection points can be found
	/// @details You can store build-specified information (loop length, desired
	/// abego, etc.) in the StructureData object which will be used later when
	/// apply_permutation() is called.
	virtual void setup_permutation( components::StructureData & perm ) const;

	/// @brief sets up the connection mover based on the information in the permutation
	/// @throw EXCN_Setup if no valid free termini can be found
	virtual void
	setup_from_random( components::StructureData & perm, core::Real random ) const;

	/// @brief checks the generated StructureData object to ensure it fits with what the user wants
	/// before building
	/// @details default behavior is to always pass this check. Subclasses can implement more stringent checks.
	virtual bool check_permutation( components::StructureData const & perm ) const;

	/// @brief checks the inserted region vs. the desired ss/abego.  True if it matches, false otherwise
	/// @details One can add checks by overriding this function
	virtual bool check( components::StructureData const & perm ) const;

	/// @brief derived classes should return true if the Connection forms a polymer bond (e.g. peptide)
	/// and false if it forms another type of chemical linkage (e.g. disulfide)
	virtual bool polymer_connection() const = 0;

	/// @brief derived classes can override this to check whether the given segments are connectable
	/// the default implementation takes into account whether or not the connection performs
	/// orientation
	virtual bool are_connectable(
		components::StructureData const & sd,
		std::string const & segment1,
		std::string const & segment2,
		Motif const & motif ) const;

	virtual core::pose::PoseOP
	build_pose( components::StructureData const & perm ) const;

	/// @throw EXCN_Process on stochastic failures
	virtual void
	process_permutation( components::StructureData & perm ) const;

	///////////////////////////////////////////////////////////////////////////
	/// Methods for setting/getting data in StructureData Object
	///////////////////////////////////////////////////////////////////////////
public:
	/// @brief get "left" residue of loop region taking overlap into account
	/// @details For example, if a user specifies a loop to start at residue
	/// 20 with an overlap of three, the "left" of the build region will
	/// actually be residue 17 (or the chain terminus, whichever is larger)
	core::Size build_left( components::StructureData const & perm ) const;

	/// @brief get "right" residue of loop region taking overlap into account
	/// @details For example, if a user specifies a loop to stop at residue
	/// 20 with an overlap of three, the "right" of the build region will
	/// actually be residue 23 (or the chain terminus, whichever is smaller)
	core::Size build_right( components::StructureData const & perm ) const;

	/// @brief length of the connection to be built
	core::Size build_len( components::StructureData const & perm ) const;
	void set_build_len( components::StructureData & perm, core::Size const len_val ) const;

	/// @brief secondary structure of the connection to be built
	std::string build_ss( components::StructureData const & perm ) const;
	void set_build_ss( components::StructureData & perm, std::string const & ss_val ) const;

	/// @brief abego of the connection to be built
	std::string build_abego( components::StructureData const & perm ) const;
	void set_build_abego( components::StructureData & perm, std::string const & abego_val ) const;

	/// @brief location of the cut to be placed (within the loop)
	core::Size cut_resi( components::StructureData const & perm ) const;
	void set_cut_resi( components::StructureData & perm, core::Size const cut_val ) const;

	/// @brief get id of the segment immediately lower to the loop
	/// @details The upper terminus of this segment will become the first residue
	/// of the loop
	std::string const & lower_segment_id( components::StructureData const & perm ) const;

	/// @brief get id of the segment immediately upper to the loop
	/// @details The lower terminus of this segment will become the last residue
	/// of the loop
	std::string const & upper_segment_id( components::StructureData const & perm ) const;

	/// @brief get id of the segment at the lower terminus of the first chain
	/// being connected
	//std::string const & comp1_lower( components::StructureData const & perm ) const;

	/// @brief get id of the segment at the upper terminus of the second chain
	/// being connected
	//std::string const & comp2_upper( components::StructureData const & perm ) const;

	/// @brief get id of the lower component actually being connected
	/// @details If the loop has non-zero length, this will be the id() of
	/// the Connection object. If the loop has zero length, this will be
	/// the same as lower_segment_id()
	std::string const & loop_lower( components::StructureData const & perm ) const;
	void set_loop_lower( components::StructureData & perm, std::string const & comp ) const;

	/// @brief get id of the upper component actually being connected
	/// @details If the loop has non-zero length, this will be the id() of
	/// the Connection object plus "_1". If the loop has zero length, this
	/// will be the same as upper_segment_id()
	std::string const & loop_upper( components::StructureData const & perm ) const;
	void set_loop_upper( components::StructureData & perm, std::string const & comp ) const;

private:
	void set_lower_segment_id( components::StructureData & perm, std::string const & comp ) const;
	void set_upper_segment_id( components::StructureData & perm, std::string const & comp ) const;
	//void set_comp1_lower( components::StructureData & perm, std::string const & comp ) const;
	//void set_comp2_upper( components::StructureData & perm, std::string const & comp ) const;
	///////////////////////////////////////////////////////////////////////////
	/// Class Member Variable Accessors/Mutators
	///////////////////////////////////////////////////////////////////////////
public:
	/// @brief derived classes should return true if the Connection moves any jumps (e.g. StapleConnection)
	/// and false if it keeps rigid body dofs fixed (e.g. ConnectJumps)
	inline bool performs_orientation() const { return performs_orientation_; }
	inline void set_performs_orientation( bool const or_val ) { performs_orientation_ = or_val; }

	/// @brief the number of residues adjacent to the loop to also include in remodeling
	core::Size lower_overlap() const { return lower_overlap_; }
	core::Size upper_overlap() const { return upper_overlap_; }
	/// @brief sets overlap for both lower and upper segments
	void set_overlap( core::Size const overlap_val );
	/// @brief sets overlap for the lower segment only
	void set_lower_overlap( core::Size const overlap_val );
	/// @brief sets overlap for the upper segment only
	void set_upper_overlap( core::Size const overlap_val );

	/// @brief sets whether to allow components to connect to themselves to create cyclic peptides
	inline bool allow_cyclic() const { return allow_cyclic_; }

	/// @details tells whether we should remodel the extended-conformation
	/// residues built by the Connection class. Default is true, and this should
	/// only be set to false for unit testing purposes.
	inline bool do_remodel() const { return do_remodel_; }
	inline void set_do_remodel( bool const val ) { do_remodel_ = val; }

	// max distance for termini to be "connectable" per loop residue
	inline core::Real connecting_bond_dist() const { return connecting_bond_dist_; }

	/// @brief Returns list of possible lengths/SS/Abegos
	inline MotifList const & motifs() const { return motifs_; }
	/// @brief set possible motifs of the connection with a comma-separated string
	void set_motifs( std::string const & motif_str );
	/// @brief set possible motifs of the connection with a vector of strings
	void set_motifs( utility::vector1< std::string > const & motif_strs );

	/// @brief name of component1
	inline utility::vector1< std::string > const & comp1_ids() const { return comp1_ids_; }
	void set_comp1_ids( std::string const & id_str );
	void set_comp1_ids( utility::vector1< std::string > const & id_list );

	/// @brief name of component2
	inline utility::vector1< std::string > const & comp2_ids() const { return comp2_ids_; }
	void set_comp2_ids( std::string const & id_str );
	void set_comp2_ids( utility::vector1< std::string > const & id_list );

	/// @brief Tells whether or not to construct motifs based on Nobu/Rie/YuRu abego rules
	inline bool idealized_abego() const { return idealized_abego_; }
	inline void set_idealized_abego( bool const val ) { idealized_abego_ = val; }

	/// @brief returns user-set chain1
	inline core::Size user_chain1() const { return chain1_; }
	inline void set_user_chain1( core::Size const chain1 ) { chain1_ = chain1; }

	/// @brief returns user-set chain2
	inline core::Size user_chain2() const { return chain2_; }
	inline void set_user_chain2( core::Size const chain2 ) { chain2_ = chain2; }

	/// @brief if set to true, abego of the build insert will be checked
	inline void set_check_abego( bool const val ) { check_abego_ = val; }

	/// @brief sets number of trials
	inline void set_trials( core::Size const tr ) { trials_ = tr; }

	// @brief set possible lengths of the connection with a string
	void set_lengths( std::string const & length_str );
	/// @brief set possible lengths of the connection with a vector
	void set_lengths( utility::vector1< core::Size > const & lengths_val );

	/// @brief set possible cut residues (relative to the start of the loop)
	void set_cut_resis( std::string const & cut_str );
	/// @brief set possible lengths of the connection with a vector
	void set_cut_resis( utility::vector1< core::Size > const & cuts_val );

	/// @brief sets the given pair of segments as not allowed
	/// @details for example, if you don't want "h1" and "h2" to be joined
	/// by a loop, call disallow_pair( "h1", "h2" )
	void disallow_pair( std::string const & c1, std::string const & c2 );
	/// @brief returns true if this specific pairing of components hasn't been
	/// explicitly disallowed by the user
	bool pair_allowed( std::string const & c1, std::string const & c2 ) const;
	/// @brief erases all information on disallowed pairs of segments
	void clear_disallowed_pairs();

	/// @brief sets whether to allow components to connect to themselves to create cyclic peptides
	inline void set_allow_cyclic( bool const cyc_val ) { allow_cyclic_ = cyc_val; }

	/// @brief sets whether or not to extend SS elements to try to connect them,
	///        or to use loop residues only. Default = true
	void set_extend_ss( bool const extend_ss );

	inline void set_connecting_bond_dist( core::Real const val ) { connecting_bond_dist_ = val; }

	void add_rcg( protocols::forge::remodel::RemodelConstraintGeneratorOP rcg );
	void clear_rcgs();

public:
	/// @brief Performs pre-build setup and makes loop residues
	void setup( components::StructureData & perm ) const;

	/// @brief checks the inserted region vs. the desired ss/abego.  True if it matches, false otherwise
	bool check_insert(
		core::pose::Pose const & pose,
		std::string const & build_ss,
		utility::vector1< std::string > const & build_abego,
		core::Size const left ) const;

	/// @brief parse a motif string, return a list of paired lengths and ss+abegos
	Motif parse_motif( std::string const & motif_str ) const;

	/// @brief finds usable/available upper termini (i.e. those for comp1)
	utility::vector1< std::string >
	find_available_upper_termini( components::StructureData const & perm ) const;

	/// @brief finds usable/available upper termini (i.e. those for comp1)
	utility::vector1< std::string >
	find_available_lower_termini( components::StructureData const & perm ) const;

	/// @brief Given desired lengths, compute a set of idealized loop motifs via Nobu/Rie/YuRu rules
	MotifList calc_idealized_motifs(
		std::string const & abego1,
		std::string const & abego2,
		std::set< core::Size > const & len_set ) const;

protected:
	void move_segments( components::StructureData & perm, StringList const & desired_order ) const;
	void move_segments_cyclic( components::StructureData & perm ) const;
	void connect_lower_loop( components::StructureData & perm ) const;
	void connect_upper_loop( components::StructureData & perm ) const;

	/// @brief performs post-connection-building tasks
	/// @details Result is a closed chain with closed fold tree and covalent bond
	/// @throw EXCN_ConnectionFailed if checks fail
	void post_process_permutation( components::StructureData & perm ) const;

	void setup_structuredata(
		components::StructureData & perm,
		core::Size const len,
		std::string const & ss,
		std::string const & abego,
		std::string const & seg1,
		std::string const & seg2,
		std::string const & seg1_lower,
		std::string const & seg2_upper,
		core::Size const cut_resi_val ) const;

	/// @brief parses subtag
	void parse_subtag( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers );

	void apply_constraints( components::StructureData & sd ) const;
	void remove_constraints( components::StructureData & sd ) const;

private:
	// component on the n-terminal side
	utility::vector1< std::string > comp1_ids_;
	// component on the c-terminal side
	utility::vector1< std::string > comp2_ids_;
	// chains to connect
	core::Size chain1_;
	core::Size chain2_;
	// number of times to try to build a connection before failing
	core::Size trials_;
	// number of residues adjacent to the loop to include in modeling
	core::Size lower_overlap_;
	core::Size upper_overlap_;
	// call remodel to fold connecting loop or just build residues?
	bool do_remodel_;
	// tells whether cyclic connections should be allowed (Default=false)
	bool allow_cyclic_;
	// the length of the bond which will join the connecting atoms of the upper/lower residue
	core::Real connecting_bond_dist_;
	// possible loop motifs
	MotifList motifs_;
	// possible locations in loops for cuts
	utility::vector1< core::Size > cut_resis_;
	// explicitly disabled pairings
	std::set< std::pair< std::string, std::string > > disallowed_pairs_;
	// constraint generators
	utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP > rcgs_;
	// Tells whether or not to construct motifs based on Nobu/Rie/YuRu abego rules
	bool idealized_abego_;
	// Tells whether or not to include SS extensions in the loop length, or just build loop residues only (defualt=true)
	bool extend_ss_;
	// whether or not orientation is performed, or just connection (default=false)
	bool performs_orientation_;
	// whether or not to check abegos after loop construction
	bool check_abego_;
};

/// @brief "Dummy" connection used when we wnat to track connection information but not build anything
class GenericConnection : public Connection {
public:
	GenericConnection();
	virtual ~GenericConnection();

	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual protocols::moves::MoverOP clone() const;

	/// @brief applies the loop building, in this case does nothing
	/// @throw EXCN_ConnectionFailed on failure to connect
	virtual void apply_connection( components::StructureData & perm ) const;

	virtual bool polymer_connection() const { return true; }

	/// @brief check whether the given segments are connectable
	virtual bool are_connectable(
		components::StructureData const & sd,
		std::string const & segment1,
		std::string const & segment2,
		Motif const & motif ) const;
};

class StapleChains : public Connection {
public:
	StapleChains();
	virtual ~StapleChains();

	//overloaded virtuals
public:
	/// @brief name
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual protocols::moves::MoverOP clone() const;

	/// @brief setup the parameters via an xml tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief sets up the connection mover based on the information in the permutation
	/// @throw EXCN_Setup if no valid connection endpoints can be found
	virtual void setup_permutation( components::StructureData & perm ) const;

	/// @brief applies the loop building
	/// @throw EXCN_ConnectionFailed on failure to connect
	virtual void apply_connection( components::StructureData & perm ) const;

	/// @brief derived classes should return true if the Connection forms a polymer bond (e.g. peptide)
	/// and false if it forms another type of chemical linkage (e.g. disulfide)
	virtual bool polymer_connection() const { return true; }

	// methods
public:

	// accessor/mutators
public:

private:
};

class BridgeTomponents : public Connection {
public:
	BridgeTomponents();
	virtual ~BridgeTomponents();

	//overloaded virtuals
public:
	/// @brief name
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief setup the parameters via an xml tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief applies the loop building
	/// @throw EXCN_ConnectionFailed on failure to connect
	virtual void apply_connection( components::StructureData & perm ) const;

	/// @brief derived classes should return true if the Connection forms a polymer bond (e.g. peptide)
	/// and false if it forms another type of chemical linkage (e.g. disulfide)
	virtual bool polymer_connection() const { return true; }

	// methods
public:
	/// @brief generates a list of loop residues for the given permutation
	/// includes N- and C- terminal residues as anchors
	std::pair< utility::vector1< core::Size >, core::Size >
	compute_loop_residues( components::StructureData const & perm ) const;

	/// @brief create and return kic mover
	protocols::generalized_kinematic_closure::GeneralizedKICOP
	create_kic_mover( components::StructureData const & perm,
		utility::vector1< core::Size > const & lres,
		core::Size const pre_overlap ) const;

	/// @brief setup kic for simple closure
	virtual void
	setup_kic_closure(
		components::StructureData const & perm,
		protocols::generalized_kinematic_closure::GeneralizedKICOP kic,
		utility::vector1< core::Size > const & lres,
		core::Size const pre_overlap ) const;

	/// @brief sets kic mover
	inline void set_kic_mover( protocols::generalized_kinematic_closure::GeneralizedKICOP const kic_val ) { kic_template_ = kic_val; }
	/// @brief sets kic selector scorefunction
	void set_selector_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

private:
	protocols::generalized_kinematic_closure::GeneralizedKICCOP kic_template_;
};

/// @brief uses KIC to introduce a disulfide between the two termini
class ConnectTerminiWithDisulfide : public BridgeTomponents {
public:
	ConnectTerminiWithDisulfide();
	virtual ~ConnectTerminiWithDisulfide();
public:
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual protocols::moves::MoverOP clone() const;

	virtual void
	setup_from_random( components::StructureData & perm, core::Real random ) const;

	/// @throw EXCN_ConnectionFailed on failure to connect
	virtual void apply_connection( components::StructureData & perm ) const;
	virtual bool polymer_connection() const { return false; }

	// methods
public:
	/// @brief given a pose and list of loop residues, creates a CYD pair
	std::pair< core::Size, core::Size >
	create_cyd_pair( components::StructureData & perm, utility::vector1< core::Size > const & loop_residues ) const;

	/// @brief setup kic for simple closure
	void
	setup_disulf_kic_closure(
		protocols::generalized_kinematic_closure::GeneralizedKICOP kic,
		utility::vector1< core::Size > const & loop_residues,
		std::pair< core::Size, core::Size > const & disulf_pos ) const;

	/// @brief setup kic for disulfide closure
	virtual void
	setup_kic_closure(
		components::StructureData const & perm,
		protocols::generalized_kinematic_closure::GeneralizedKICOP kic,
		utility::vector1< core::Size > const & lres,
		core::Size const pre_overlap ) const;
};

void staple_work_function(
	components::StructureData & perm,
	std::string const & loop_lower,
	std::string const & loop_upper,
	bool const polymer_connection,
	bool const reverse );

core::Real
calc_approx_loop_length( std::string const & abego );

/// @brief compares desired insert ss and abego to pose values, returns number of mismatches
core::Size compare_insert_ss_and_abego(
	std::string const & pose_ss,
	utility::vector1< std::string > const & pose_abego,
	std::string const & build_ss,
	utility::vector1< std::string > const & build_abego,
	core::Size const left );

bool check_insert_ss_and_abego(
	std::string const & pose_ss,
	utility::vector1< std::string > const & pose_abego,
	std::string const & build_ss,
	utility::vector1< std::string > const & build_abego,
	core::Size const left );

class EXCN_ConnectionFailed : public utility::excn::EXCN_Base {
public:
	EXCN_ConnectionFailed( std::string const & msg ):
		utility::excn::EXCN_Base(),
		message_( msg ) {}
	std::string const & message() const { return message_; }
	virtual void show( std::ostream & os ) const { os << message_; }
private:
	std::string const message_;
};

class EXCN_UnknownSubtag : public utility::excn::EXCN_RosettaScriptsOption {
public:
	EXCN_UnknownSubtag( std::string const & msg ) : utility::excn::EXCN_RosettaScriptsOption( msg ) {}
};

} // connection
} // denovo_design
} // protocols

#endif
