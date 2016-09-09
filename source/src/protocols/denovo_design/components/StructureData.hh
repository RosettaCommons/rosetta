// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/components/StructureData.hh
/// @brief StructureData functions for building structures from components
/// @details
/// @author Tom Linsky

#ifndef INCLUDED_protocols_denovo_design_components_StructureData_hh
#define INCLUDED_protocols_denovo_design_components_StructureData_hh

// Unit headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.fwd.hh>
#include <protocols/denovo_design/components/StructureDataObserver.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Core headers
#include <utility/graph/Graph.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/Remarks.fwd.hh>
#include <core/types.hh>

// Basic/Numeric/Utility Headers
#include <basic/datacache/WriteableCacheableData.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

// C++ Headers
#include <map>
#include <set>
#include <string>

namespace protocols {
namespace denovo_design {
namespace components {

typedef std::map< std::string, std::string > StringMap;
typedef std::map< std::string, Segment > SegmentMap;

class Alias : public std::pair< std::string, core::Size > {
public:
	Alias():
		std::pair< std::string, core::Size >() {}

	Alias( std::string const & segment, core::Size const resid ):
		std::pair< std::string, core::Size >( segment, resid ) {}
};
typedef std::map< std::string, Alias > AliasMap;

class BondInfo : public utility::pointer::ReferenceCount {
public:
	BondInfo(): seg1( "" ), seg2( "" ), res1( 0 ), res2( 0 ), atom1( "" ), atom2( "" ) {}

	BondInfo(
		std::string const & s1,
		std::string const & s2,
		core::Size const r1,
		core::Size const r2,
		std::string const & a1,
		std::string const & a2
	): seg1( s1 ), seg2( s2 ), res1( r1 ), res2( r2 ), atom1( a1 ), atom2( a2 ) {}

	bool operator==( BondInfo const & other ) const {
		if ( res1 != other.res1 ) return false;
		if ( res1 != other.res2 ) return false;
		if ( seg1 != other.seg1 ) return false;
		if ( seg2 != other.seg2 ) return false;
		if ( atom1 != other.atom1 ) return false;
		if ( atom2 != other.atom2 ) return false;
		return true;
	}
	friend std::ostream & operator<<( std::ostream & os, BondInfo const & b ) {
		os << "\t<CovalentBond segment1=\"" << b.seg1 << "\" segment2=\"" << b.seg2
			<< "\" residue1=\"" << b.res1 << "\" residue2=\"" << b.res2
			<< "\" atom1=\"" << b.atom1 << "\" atom2=\"" << b.atom2 << "\" />";
		return os;
	}
	std::string seg1, seg2;
	core::Size res1, res2;
	std::string atom1, atom2;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
};
typedef utility::vector1< BondInfo > BondInfos;

typedef utility::vector1< core::Size > Cutpoints;
typedef utility::vector1< MovableGroup > MovableGroups;
typedef utility::vector1< SegmentPairingCOP > SegmentPairingCOPs;

// StructureData objects -- contains functionality for building components
class StructureData : public basic::datacache::WriteableCacheableData {
private:
	typedef basic::datacache::CacheableDataOP CacheableDataOP;
	typedef std::string SegmentName;
	typedef std::set< SegmentName > SegmentNameSet;

	// object construction
public:
	StructureData();

	StructureData( std::string const & id_val );

	virtual ~StructureData();

	static std::string
	class_name();

public:
	// WriteableCacheableData virtuals
	virtual CacheableDataOP
	clone() const;

	virtual void
	write( std::ostream & os ) const;

	virtual std::string
	datatype() const;

public:
	/// @brief Retrieves data from an XML tag
	void
	parse_tag( utility::tag::TagCOP tag );

	/// @brief Retrieves data from an XML subtag
	void
	parse_subtag( utility::tag::TagCOP tag );

	// constants
public:
	static char const DATA_DELIMETER;

protected:
	/// @brief overridden by derived classes if they need to do anything when the SD changes
	void
	on_change();

	// non-const methods
public:
	/// @brief sets id name
	void
	set_id( std::string const & id_str );

	/// @brief attaches a template pose to the given segment
	void
	set_template_pose(
		std::string const & segment,
		core::pose::Pose const & template_pose,
		core::Size const start_resid,
		core::Size const stop_resid );

	/// @brief renames a residue segment and updates all connections
	void
	add_prefix_to_segments( std::string const & prefix );

	/// @brief renames a residue segment and updates all connections
	void
	rename_segment( std::string const & old_name, std::string const & new_name );

	/// @brief re-arranges segments suchs that segment2 follows segment1 in sequence
	void
	move_segment( std::string const & segment1, std::string const & segment2 );

	/// @brief connects the given chains together, doesn't update anything -- don't call this on its own unless you know what you're doing
	void
	connect_segments(
		std::string const & segment1,
		std::string const & segment2 );

	/// @brief connects the given chains together, doesn't update anything -- don't call this on its own unless you know what you're doing
	void
	disconnect_segments(
		std::string const & segment1,
		std::string const & segment2 );

	/// @brief merges two segments into one that has the name new_name. They must be next to each other in sequence.
	void
	merge_segments(
		std::string const & segment1,
		std::string const & segment2,
		std::string const & new_name );

	/// @brief removes all traces of the given segment from the object
	void
	delete_segment( std::string const & segment );

	void
	delete_residue( core::Size const pose_resid );

	/// @brief declares a covalent bond between the specified atoms
	void declare_covalent_bond(
		std::string const & seg1, core::Size const res1, std::string const & atom1,
		std::string const & seg2, core::Size const res2, std::string const & atom2 );

	/// @brief declares a covalent bond using pose residues
	void declare_covalent_bond(
		core::Size const res1, std::string const & atom1,
		core::Size const res2, std::string const & atom2 );

	/// @brief returns list of all cutpoints, in N-->C order
	Cutpoints
	cutpoints() const;

	/// @brief marks the ith residue as a cutpoint
	void
	set_cutpoint( core::Size const resid );

	/// @brief marks the resi-th residue of segment as a cutpoint
	void
	set_cutpoint( std::string const & seg, SegmentResid const resi );

	/// @brief replace segment named seg_name with the given new segment
	void
	replace_segment( SegmentName const & seg_name, Segment const & segment );

	/// @brief adds a residues segment to the end of the list
	void
	add_segment( Segment const & resis );

	/// @brief adds a residues segment -- will be ordered before the given segment
	void
	add_segment( Segment const & resis, std::string const & insert_before_segment );

	/// @brief adds a residues segment -- will be inserted just before at given iterator
	void
	add_segment( Segment const & resis, SegmentNameList::iterator insert_pos );

	/// @brief Returns a StructureData containing only the given segments
	///        The resulting StructureData will contain all data from the current
	StructureData
	slice( SegmentNameSet const & segments, bool const force_padding ) const;

	/// @brief merge all data and segments from "other" into this StructureData
	void merge( StructureData const & other );
	/// @brief merge given data and segments from "other" into this StructureData
	void merge( StructureData const & other, SegmentNames const & segments );

	/// @brief merge all data and segments from "other" into this StructureData before the given position
	void merge_before( StructureData const & other, std::string const & position );
	/// @brief merge all data and given segments from "other" into this StructureData before the given position
	void merge_before( StructureData const & other, std::string const & position, SegmentNames const & segments );

	/////////////////////////////////////////////////////////////////////////////
	/// Aliases
	/////////////////////////////////////////////////////////////////////////////
public:
	/// @brief given a residue alias, returns a pose residue number
	bool has_alias( std::string const & alias ) const { return ( aliases_.find(alias) != aliases_.end() ); }

	/// @brief given a residue alias, returns a pose residue number
	core::Size
	alias( std::string const & alias ) const;

	/// @brief sets an "alias" for a particular residue inside a segment which allows for it to be easily accessed
	void
	set_alias(
		std::string const & alias_name,
		std::string const & segment_name,
		core::Size const resi );

	/// @brief sets an "alias" for a particular residue which allows for it to be easily accessed
	void
	set_alias(
		std::string const & alias_name,
		core::Size const resi );

	/////////////////////////////////////////////////////////////////////////////
	/// User-specified Data
	/////////////////////////////////////////////////////////////////////////////
public:
	/// @brief check for integer data
	bool
	has_data_int( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief check for real number data
	bool
	has_data_real( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets real number data
	bool
	has_data_str( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets integer data
	int
	get_data_int( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets real number data
	core::Real
	get_data_real( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets real number data
	std::string const &
	get_data_str( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief sets real number data
	void
	set_data_int( std::string const & segment_id, std::string const & data_name, int const val );

	/// @brief sets real number data
	void
	set_data_real( std::string const & segment_id, std::string const & data_name, core::Real const val );

	/// @brief sets real number data
	void
	set_data_str( std::string const & segment_id, std::string const & data_name, std::string const & val );

	/// @brief copies user data fields from one permutation to this one -- overwrites existing data
	void copy_data( StructureData const & perm );

private:
	/// @brief copies user data fields from one permutation to this one -- optionally overwrites
	void copy_data( StructureData const & perm, bool const overwrite );

	// const methods
public:
	/// @brief Total number of chains WARNING: This is an O(n) operation, where n is number of residue segments
	core::Size
	num_chains() const;

	/// @brief returns the actual residue number of the given name and res #
	core::Size
	pose_residue( std::string const & segment_name, core::Size const local_res ) const;

	/// @brief true if this permutation contains a residue segment named seg
	bool
	has_segment( std::string const & seg ) const;

	/// @brief returns an ordered list of segments which are all connected containing seg
	/// @param stop_at_cutpoint - if true, segments past a cutpoint will not be counted as connected
	SegmentNameList
	connected_segments( std::string const & seg, bool const stop_at_cutpoint ) const;

	/// @brief computes and returns a set of segments which are in the given movable group
	SegmentNames
	segments_in_movable_group( core::Size const group ) const;

	/// @brief computes and returns a set of movable groups
	MovableGroups const &
	movable_groups() const;

	/// @brief movable group of segment which contains residue resid
	core::Size
	movable_group( core::Size const resid ) const;

	/// @brief returns segments which have free lower termini
	SegmentNames
	available_lower_termini() const;

	/// @brief returns segments which have free upper termini
	SegmentNames
	available_upper_termini() const;

	/// @brief start of segments list
	SegmentNameList::const_iterator
	segments_begin() const;

	/// @brief end of segment list
	SegmentNameList::const_iterator
	segments_end() const;

	/// @brief start/end of covalent bonds list
	BondInfos::const_iterator
	covalent_bonds_begin() const;

	/// @brief end of covalent bonds list
	BondInfos::const_iterator
	covalent_bonds_end() const;

	/// @brief tells whether a non-polymer bond exists between the given segments
	bool
	non_polymer_bond_exists( std::string const & seg1, std::string const & seg2 ) const;

	/// @brief finds a non-peptide bond between two segments, returns end() if there isn't one
	BondInfo const &
	non_polymer_bond( std::string const & seg1, std::string const & seg2 ) const;

	/// @brief finds a non-peptide bond between two segments, returns end() if there isn't one
	BondInfos::const_iterator
	find_non_polymer_bond( std::string const & seg1, std::string const & seg2 ) const;

	////////////////////////////////////////////////////////////////////
	/// Segment retrieval
	////////////////////////////////////////////////////////////////////

	/// @brief returns segment represented by given string
	Segment const &
	segment( std::string const & id_val ) const;

	/// @brief finds a segment in the segment_order list and returns an iterator to it
	SegmentNameList::const_iterator
	find_segment_name( std::string const & segname ) const;

	/// @brief finds a segment in the segment map and returns an iterator to it
	SegmentMap::const_iterator
	find_segment( std::string const & segname ) const;

	/// @brief returns n-terminal residue of the chain represented by given string
	core::Size lower_anchor( std::string const & id_val ) const;

	/// @brief returns c-terminal residue of the chain represented by given string
	core::Size upper_anchor( std::string const & id_val ) const;

	/// @brief returns segment which includes residue number res
	std::string const & segment_name( core::Size const res ) const;

	/// @brief returns n and c terminal segments of the chain which includes seg
	std::pair< std::string, std::string > termini( std::string const & seg ) const;

	/// @brief tells if the segment given has an available lower terminus
	bool has_free_lower_terminus( std::string const & id_val ) const;

	/// @brief tells if the segment given has an available lower terminus
	bool has_free_upper_terminus( std::string const & id_val ) const;

	/// @brief returns the id of this permutation
	std::string const &
	id() const;

	/// @brief returns the length of this permutation
	core::Size
	length() const;

	/// @brief returns the total length of this permutation, including n-, c-terminal loop residues which are basically for show
	core::Size
	pose_length() const;

	/////////////////////////////////////////////////////////////////////////////
	/// Data storage/access
	/////////////////////////////////////////////////////////////////////////////
public:
	/// @brief returns true if this object has a group of segments with the given name
	bool has_segment_group( std::string const & sname ) const;

	/// @brief returns true if this object has a group of segments with the given name
	SegmentNames
	segment_group( std::string const & sname ) const;

	/// @brief gets all real number data
	inline std::map< std::string, int > const & data_int() const { return data_int_; }

	/// @brief gets all real number data
	inline std::map< std::string, core::Real > const & data_real() const { return data_real_; }

	/// @brief gets all string data
	inline std::map< std::string, std::string > const & data_str() const { return data_str_; }

	/// @brief gets all alias data
	AliasMap const & aliases() const { return aliases_; }

	/// @brief return secondary structure string
	std::string const &
	ss() const;

	/// @brief abego of residue resid
	char
	ss( core::Size const resid ) const;

	/// @brief sets secondary structure for residue resid
	void
	set_ss( core::Size const resid, char const ss_type );

	/// @brief return abego string
	std::string const &
	abego() const;

	/// @brief abego of residue resid
	char
	abego( core::Size const resid ) const;

	void
	set_abego( std::string const & segment, std::string const & abego );

	void
	set_abego( std::string const & segment, utility::vector1< std::string > const & abego );

	/// @brief given an input stream, substitute all variables
	/// @details variables are of the form: %%SEGMENTNAME#residue%%
	/// SEGMENTNAME = name of the segment
	/// residue = local residue number within the segment
	/// The substituted value will be an core::Size corresponding to the pose residue
	std::string substitute_variables( std::istream & input ) const;

	/// @brief checks consistency of the data
	/// @throws EXCN_PoseInconsistent if there is a problem
	void
	check_consistency() const;

	/// @brief checks the permutation for internal consistency vs a pose
	/// @throws EXCN_PoseInconsistent if there is a problem
	void
	check_pose_consistency( core::pose::Pose const & pose ) const;

	/// @brief for output
	friend std::ostream & operator<<( std::ostream & os, StructureData const & perm );

	/// @brief marks the given segments as covanlently connected
	void mark_connected(
		std::string const & lower_seg,
		std::string const & upper_seg );

	/// @brief unmarks the given segments as covalently connected
	void mark_disconnected(
		std::string const & seg1,
		std::string const & seg2 );

	/// @brief removes jump and cutpoint between the two segments to create a single polymer chain
	void
	delete_jump_and_intervening_cutpoint( std::string const & segment1, std::string const & segment2 );

	/// @brief deletes the residues between the segment N terminus and the N anchor point
	void delete_leading_residues( std::string const & seg );

	/// @brief deletes the residues between the segment C terminus and the C anchor point
	void delete_trailing_residues( std::string const & seg );

	/// @brief chooses a new movable group which doesn't conflict with existing ones
	core::Size choose_new_movable_group() const;

	/// @brief sets movable group of a segment
	void set_movable_group( std::string const & id, core::Size const mg );

	/// @brief renumbers movable group "oldg" to have new number "newg"
	void renumber_movable_group( core::Size const oldg, core::Size const newg );

	/// @brief updates numbering based on the saved order of Segment objects
	void
	update_numbering();

	/// @brief Saves given remarks changes enzdes residues to generic segment name/number
	void
	save_remarks( core::io::Remarks const & remarks );

	/// @brief retrieves saved remarks, makes any enzdes residues specific to the given pose
	core::io::Remarks
	retrieve_remarks( core::pose::Pose const & pose ) const;

protected:
	/// helper functions for data access
	int get_data_int( std::string const & data_name ) const;
	core::Real get_data_real( std::string const & data_name ) const;
	std::string const & get_data_str( std::string const & data_name ) const;
	bool has_data_int( std::string const & data_name ) const { return data_int_.find(data_name) != data_int_.end(); }
	bool has_data_real( std::string const & data_name ) const { return data_real_.find(data_name) != data_real_.end(); }
	bool has_data_str( std::string const & data_name ) const { return data_str_.find(data_name) != data_str_.end(); }
	void set_data_int( std::string const & data_name, int const val );
	void set_data_real( std::string const & data_name, core::Real const val );
	void set_data_str( std::string const & data_name, std::string const & val );

private:
	/// @brief blocks on_change() signals
	void
	block_signals();

	/// @brief unblocks on_change() signals and calls on_change()
	void
	unblock_signals();

	/// @brief should be called when something changes
	void
	changed();

protected:
	/// @brief returns constant list of residue ranges
	inline SegmentMap const & segments() const { return segments_; }

	/// @brief returns non-const access to residue range of the segment represented by given string
	Segment &
	segment_nonconst( std::string const & id_val );

	/// @brief finds a segment in the segment_order list and returns an iterator to it
	SegmentNameList::iterator
	find_segment_name( std::string const & segname );

	/// @brief non-const iterator to end of segment names list
	SegmentNameList::iterator
	segments_end_nonconst();

	/// @brief finds a segment in the segments map and returns a non-const iterator
	SegmentMap::iterator
	find_segment( std::string const & segname );

	/// @brief performs dfs in lower direction looking for termini
	std::string
	find_lower_terminus( std::set< std::string > & visited, std::string const & seg ) const;

	/// @brief performs dfs in upper direction looking for termini
	std::string
	find_upper_terminus( std::set< std::string > & visited, std::string const & seg ) const;

	/// @brief safely slide a jump, avoiding foldtree segmentation faults
	core::kinematics::FoldTree
	slide_jump(
		core::kinematics::FoldTree const & ft_orig,
		core::Size const jump_idx,
		core::Size const new_start,
		core::Size const new_stop ) const;

public:
	// temprarily public
	void
	add_covalent_bond(
		core::Size const res1, std::string const & atom1,
		core::Size const res2, std::string const & atom2 );

	void
	clear_pairings();

	void
	add_pairing( SegmentPairing const & pairing );

	SegmentPairingCOPs::const_iterator
	pairings_begin() const;

	SegmentPairingCOPs::const_iterator
	pairings_end() const;

protected:
	void
	add_covalent_bond(
		std::string const & seg1, core::Size const res1, std::string const & atom1,
		std::string const & seg2, core::Size const res2, std::string const & atom2 );

	void
	add_covalent_bond( BondInfo const & bi );

	/// @brief renames a residue segment and updates all connections
	void add_prefix_to_segments( std::string const & prefix, char const delimeter );

protected:
	void
	move_segments(
		std::string const & segment1,
		std::string const & segment2_lower,
		std::string const & segment2_upper );

	// internal bookkeeping functions
private:
	/// @brief checks pose vs. StructureData info
	/// @throw EXCN_PoseInconsistent if things don't match
	void check_pose( core::pose::Pose const & pose ) const;

	/// @brief checks residues in SD -- makes sure everything is sequential and accounted for
	void check_residues() const;

	/// @brief checks chain termini in pose vs SD
	/// @throws EXCN_PoseInconsistent if there is a problem
	void check_improper_termini( core::pose::Pose const & pose ) const;
	void check_chains( core::pose::Pose const & pose ) const;

	/// @brief check pose movable groups vs SD
	void check_movable_groups() const;

	/// @brief called to generalize an enzdes header to track it as pose changes
	std::string
	generic_enzdes_header( std::string const & remark_str ) const;

	/// @brief called to make an enzdes header specific to a pose
	std::string
	specific_enzdes_header(
		core::pose::Pose const & pose,
		std::string const & generic_remark_str ) const;

	// member variables
private:
	// for identification with a segment
	std::string id_;
	// secondary structure string
	std::string ss_;
	// abego string
	std::string abego_;
	// length of the pose, including n-, c-terminal loop residues
	core::Size pose_length_;
	// usable length of the pose, not including n-, c-terminal loops
	core::Size length_;
	// sub-permutations for compound objects
	std::map< std::string, int > data_int_;
	// arbitrary real number data that can be set by movers
	std::map< std::string, core::Real > data_real_;
	// arbitrary string data that can be set by movers
	StringMap data_str_;
	// map of segment to residue [start, end]
	SegmentMap segments_;
	// names given to special single residues -- value stored is segment name + intra-segment residue number
	AliasMap aliases_;
	// non-canonical covalent bonds
	utility::vector1< BondInfo > covalent_bonds_;
	// segments listed in order
	SegmentNameList segment_order_;
	// movable groups
	MovableGroups movable_groups_;
	// Segment pairings
	SegmentPairingCOPs pairings_;
	// whether on_change() is called
	bool block_signals_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class EXCN_PoseInconsistent : public utility::excn::EXCN_Base {
public:
	EXCN_PoseInconsistent( std::string const & your_msg ):
		utility::excn::EXCN_Base(), msg_( your_msg ) {}

	virtual void
	show( std::ostream & os ) const { os << msg_; }
	std::string msg_;
private:
	EXCN_PoseInconsistent();
};

} // components
} // denovo_design
} // protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_denovo_design_components_StructureData )
CEREAL_FORCE_DYNAMIC_INIT( protocols_denovo_design_components_BondInfo )
#endif // SERIALIZATION


#endif // header guard
