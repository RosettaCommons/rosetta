// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/StructureData.hh
/// @brief StructureData functions for building structures from components
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_StructureData_hh
#define INCLUDED_protocols_denovo_design_components_StructureData_hh

// Unit headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Core headers
#include <core/graph/Graph.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/types.hh>

// Basic/Numeric/Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <map>
#include <set>
#include <string>

namespace protocols {
namespace denovo_design {
namespace components {

typedef std::map< std::string, Segment > SegmentMap;
class Alias : public std::pair< std::string, core::Size > {
public:
	Alias():
		std::pair< std::string, core::Size >() {}

	Alias( std::string const & segment, core::Size const resid ):
		std::pair< std::string, core::Size >( segment, resid ) {}
};

struct BondInfo {
private:
	BondInfo(): seg1( "" ), seg2( "" ), res1( 0 ), res2( 0 ), atom1( "" ), atom2( "" ) {}
public:
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
};

// StructureData objects -- contains functionality for building components
class StructureData : public utility::pointer::ReferenceCount {

	// object construction
public:
	StructureData( std::string const & id_val );
	StructureData( StructureData const & perm );
	virtual ~StructureData();
	virtual bool is_multi() const = 0;
	virtual StructureDataOP fresh_instance() const = 0;
	virtual StructureDataOP clone() const = 0;

	/// @brief creates permutation from a given pose. Remark records in the pose are NOT used.
	static StructureDataOP infer_from_pose( core::pose::Pose const & pose, std::string const & id );

	/// @brief sets data from a given pose.  Data is stored in pose
	static StructureDataOP create_from_pose( core::pose::Pose const & pose, std::string const & id );

	/// @brief parses PDB remarks and creates a permutation from them
	static StructureDataOP create_from_remarks( core::pose::Remarks const & rem, std::string const & newid );

	/// @brief loads data from pdb remarks into this permutation
	static StructureDataOP parse_remarks( core::pose::Remarks const & rem, std::string const & newid );

	/// @brief creates a StructureData from an xml stringstream
	static StructureDataOP create_from_xml( std::istream & xmltag, std::string const & newid );

	// constants
public:
	static int const REMARK_NUM;
	static std::string const DATA_NAME;
	static char const DATA_DELIMETER;

	// non-const methods
public:
	/// @brief returns the pose which defines this permutation, building it if necessary
	void consolidate_movable_groups( utility::vector1< std::string > const & root_segments );

	/// @brief moves around jumps so they are located in the center of rigid residue groups
	void move_jumps_to_safety();

	/// @brief renames a residue segment and updates all connections
	void add_prefix_to_segments( std::string const & prefix );

	/// @brief renames a residue segment and updates all connections
	void rename_segment( std::string const & old_name, std::string const & new_name );

	/// @brief aligns the upper-terminal residue of segment1 to the start "anchor" residue of segment2
	void align_segments(
		std::string const & segment1,
		std::string const & segment2 );

	/// @brief aligns the lower-terminal residue of segment1 to the end "anchor" residue of segment2
	void align_segments_rev(
		std::string const & segment1,
		std::string const & segment2 );

	/// @brief moves jump pointing to the first segment so that it will move with the second segment
	void slide_jump(
		std::string const & child_segment,
		std::string const & parent_segment );

	/// @brief aligns the template_res'th residue of template_segment to the movable_res'th residue of movable_seg
	void align_residues(
		std::string const & template_seg,
		core::Size const template_res,
		std::string const & movable_seg,
		core::Size const movable_res );

	/// @brief re-arranges residues such that segment 2 follows segment 1 in sequence, return value is number of jump to segment2
	void move_segment(
		std::string const & segment1_n,
		std::string const & segment1_c,
		std::string const & segment2_n,
		std::string const & segment2_c );

	/// @brief connects the given chains together, doesn't update anything -- don't call this on its own unless you know what you're doing
	void connect_segments(
		std::string const & segment1,
		std::string const & segment2 );

	/// @brief connects the given chains together, doesn't update anything -- don't call this on its own unless you know what you're doing
	void disconnect_segments(
		std::string const & segment1,
		std::string const & segment2 );

	/// @brief merges two segments into one that has the name new_name. They must be next to each other in sequence.
	void merge_segments(
		std::string const & segment1,
		std::string const & segment2,
		std::string const & new_name );

	/// @brief removes all traces of the given segment from the object
	void delete_segment( std::string const & segment );

	/// @brief declares a covalent bond between the specified atoms
	void declare_covalent_bond(
		std::string const & seg1, core::Size const res1, std::string const & atom1,
		std::string const & seg2, core::Size const res2, std::string const & atom2 );

	/// @brief declares a covalent bond using pose residues
	void declare_covalent_bond(
		core::Size const res1, std::string const & atom1,
		core::Size const res2, std::string const & atom2 );

	/// @brief declares a covalent bond between the upper-terminal atom of segment 1 and the lower-terminal atom of segment 2
	void declare_polymer_bond( std::string const & segment1, std::string const & segment2 );

	/// @brief marks the resi-th residue of segment as a cutpoint
	void set_cutpoint( std::string const & seg, core::Size const resi );

	/// @brief deletes the pose
	inline void clear_pose() { pose_ = core::pose::PoseOP( NULL ); }

	/// @brief adds a residues segment to the end of the list
	void add_segment( std::string const & id_val, Segment const & resis );

	/// @brief adds a residues segment -- will be ordered at the given index
	void add_segment( std::string const & id_val, Segment const & resis, std::string const & insert_before_segment );

	/// @brief adds a residues segment -- will be ordered at the given index
	void add_segment( std::string const & id_val, Segment const & resis, StringList::iterator insert_pos );

	/// @brief adds a residues segment -- also adds pose residues
	void add_segment(
		std::string const & id_val,
		Segment const & resis,
		StringList::iterator insert_pos,
		core::pose::PoseCOP residues );

	/// @brief adds a segment of residues -- start_resid MUST be the nterm_resi() of a segment if segments exist yet
	void add_segment(
		std::string const & id_val,
		StringList::iterator insert_before_pos,
		core::Size const segment_length,
		core::Size const local_safe_residue,
		core::Size const local_cutpoint,
		core::Size const movable_group,
		bool const is_loop,
		bool const nterm_included,
		bool const cterm_included,
		std::string const & lower_conn,
		std::string const & upper_conn,
		std::string const & ss,
		utility::vector1< std::string > const & abego );

	/// @brief merge all data and segments from "other" into this StructureData
	void merge( StructureData const & other );
	/// @brief merge given data and segments from "other" into this StructureData
	void merge( StructureData const & other, StringList const & segments );

	/// @brief merge all data and segments from "other" into this StructureData before the given position
	void merge_before( StructureData const & other, std::string const & position );
	/// @brief merge all data and given segments from "other" into this StructureData before the given position
	void merge_before( StructureData const & other, std::string const & position, StringList const & segments );

	/// @brief sets real number data
	void set_data_int( std::string const & segment_id, std::string const & data_name, int const val );

	/// @brief sets real number data
	void set_data_real( std::string const & segment_id, std::string const & data_name, core::Real const val );

	/// @brief sets real number data
	void set_data_str( std::string const & segment_id, std::string const & data_name, std::string const & val );

	/// @brief sets an "alias" for a particular residue inside a segment which allows for it to be easily accessed
	void set_resnum_alias(
		std::string const & alias_name,
		std::string const & segment_name,
		core::Size const resi );

	/// @brief sets an "alias" for a particular residue which allows for it to be easily accessed
	void set_resnum_alias(
		std::string const & alias_name,
		core::Size const resi );

	/// @brief copies user data fields from one permutation to this one -- overwrites existing data
	void copy_data( StructureData const & perm );

private:
	/// @brief copies user data fields from one permutation to this one -- optionally overwrites
	void copy_data( StructureData const & perm, bool const overwrite );

	// const methods
public:
	/// @brief Total number of chains WARNING: This is an O(n) operation, where n is number of residue segments
	core::Size num_chains() const;

	/// @brief returns the pose which defines this permutation
	core::pose::PoseCOP pose() const;

	/// @brief returns the actual residue number of the given name and res #
	core::Size pose_residue( std::string const & segment_name, core::Size const local_res ) const;

	/// @brief true if this permutation contains a residue segment named seg
	bool has_segment( std::string const & seg ) const;

	/// @brief true if this permutation contains a residue segment named seg
	inline core::Size num_segments() const { return segments_.size(); }

	/// @brief finds a jump pointing directly to the segment represented by given string
	/// @details returns 0 if the foldtree is rooted at the segment, returns -1 if there is no jump pointing to the segment
	int find_jump( std::string const & seg ) const;

	/// @brief returns an ordered list of segments which are all connected containing seg
	StringList connected_segments( std::string const & seg ) const;

	/// @brief computes and returns a set of segments which are in the given movable group
	utility::vector1< std::string > segments_in_movable_group( core::Size const group ) const;

	/// @brief computes and returns a set of movable groups
	std::set< core::Size > movable_groups() const;

	/// @brief counts and returns the number of movable residue groups
	core::Size movable_group( std::string const & id ) const;

	/// @brief returns segments which have free lower termini
	utility::vector1< std::string > available_lower_termini() const;

	/// @brief returns segments which have free upper termini
	utility::vector1< std::string > available_upper_termini() const;

	/// @brief start of segments list
	StringList::const_iterator segments_begin() const;

	/// @brief end of segment list
	StringList::const_iterator segments_end() const;

	/// @brief start/end of covalent bonds list
	utility::vector1< BondInfo >::const_iterator covalent_bonds_begin() const;
	utility::vector1< BondInfo >::const_iterator covalent_bonds_end() const;

	/// @brief finds a non-peptide bond between two segments, returns end() if there isn't one
	utility::vector1< BondInfo >::const_iterator
	non_peptidic_bond( std::string const & seg1, std::string const & seg2 ) const;

	/// @brief finds a segment in the segment_order list and returns an iterator to it
	StringList::const_iterator find_segment( std::string const & segname ) const;

	/// @brief returns n-terminal residue of the chain represented by given string
	core::Size lower_anchor( std::string const & id_val ) const;

	/// @brief returns c-terminal residue of the chain represented by given string
	core::Size upper_anchor( std::string const & id_val ) const;

	/// @brief returns non-const access to residue range of the segment represented by given string
	Segment & segment_nonconst( std::string const & id_val );

	/// @brief returns residue range of the segment represented by given string
	Segment const & segment( std::string const & id_val ) const;

	/// @brief finds all segments that are loops
	utility::vector1< std::string > loops() const;

	/// @brief finds and returns whether each residue is a loop
	utility::vector1< bool > loop_residues() const;

	/// @brief returns segment which includes residue number res
	std::string const & segment_name( core::Size const res ) const;

	/// @brief returns n and c terminal segments of the chain which includes seg
	std::pair< std::string, std::string > termini( std::string const & seg ) const;

	/// @brief tells if the segment given has an available lower terminus
	bool has_free_lower_terminus( std::string const & id_val ) const;

	/// @brief tells if the segment given has an available lower terminus
	bool has_free_upper_terminus( std::string const & id_val ) const;

	/// @brief gives a quick yes-no answer as to whether it might be possible to connect these termini with an nres-residue loop
	/// pose_ MUST BE SET if use_distance=true!!
	bool are_connectable(
		std::string const & id1,
		std::string const & id2,
		core::Size const nres,
		bool const use_distance,
		bool const connection_performs_orientation,
		bool const allow_cyclic,
		core::Real const bond_dist,
		core::Real const max_dist_per_res ) const;

	/// @brief returns the id of this permutation
	inline std::string const & id() const { return id_; }

	/// @brief returns the length of this permutation
	inline core::Size length() const { return length_; }

	/// @brief returns the total length of this permutation, including n-, c-terminal loop residues which are basically for show
	inline core::Size pose_length() const { return pose_length_; }

	/////////////////////////////////////////////////////////////////////////////
	/// I/O with pose
	/////////////////////////////////////////////////////////////////////////////
public:
	/// @brief checks to see whether a string for this permutation exists in the pose's datacache
	static bool has_cached_string( core::pose::Pose const & pose );

	/// @brief retrieve a string stored in the pose's datacache
	static std::string cached_string( core::pose::Pose const & pose );

	/// @brief stores the data of this permutation into a pose remarks
	void save_into_pose( core::pose::Pose & pose ) const;

protected:
	/// @brief stores the data in the permutation
	void save_into_pose();

	/// @brief checks to see whether a string for this permutation exists in the pose's datacache
	bool has_cached_string() const;

	/// @brief retrieve a string stored in the pose's datacache
	std::string cached_string() const;

	/// @brief retrieve a string stored in the pose's datacache
	static std::string cached_string( core::pose::Pose const & pose, std::string const & data_name );

	/// @brief retrieves cached remarks from pose datacache
	core::pose::Remarks cached_remarks() const;

	/// @brief retrieves cached remarks from pose datacache
	core::pose::Remarks cached_remarks( core::pose::Pose const & pose ) const;

	/// @brief stores a string in the pose's datacache
	void set_cached_string( std::string const & ss );

	/// @brief stores a string in the pose's datacache
	static void set_cached_string( core::pose::Pose & pose, std::string const & ss );

	/// @brief stores a string in the pose's datacache
	static void set_cached_string( core::pose::Pose & pose, std::string const & ss, std::string const & data_name );

private:
	/// @brief loads data from pdb remarks into this permutation
	void load_pdb_info_old( core::pose::Remarks const & rem, std::string const & prefix );

	/////////////////////////////////////////////////////////////////////////////
	/// Data storage/access
	/////////////////////////////////////////////////////////////////////////////
public:
	/// @brief check for integer data
	bool has_data_int( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief check for real number data
	bool has_data_real( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets real number data
	bool has_data_str( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets integer data
	int get_data_int( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets real number data
	core::Real get_data_real( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief gets real number data
	std::string const & get_data_str( std::string const & segment_id, std::string const & data_name ) const;

	/// @brief given a residue alias, returns a pose residue number
	bool has_alias( std::string const & alias ) const { return ( aliases_.find(alias) != aliases_.end() ); }

	/// @brief given a residue alias, returns a pose residue number
	core::Size alias_resnum( std::string const & alias ) const;

	/// @brief returns true if this object has a group of segments with the given name
	bool has_segment_group( std::string const & sname ) const;

	/// @brief returns true if this object has a group of segments with the given name
	StringList segment_group( std::string const & sname ) const;

	/// @brief gets all real number data
	inline std::map< std::string, int > const & data_int() const { return data_int_; }

	/// @brief gets all real number data
	inline std::map< std::string, core::Real > const & data_real() const { return data_real_; }

	/// @brief gets all string data
	inline std::map< std::string, std::string > const & data_str() const { return data_str_; }

	/// @brief gets all alias data
	std::map< std::string, Alias > const & aliases() const { return aliases_; }

	/// @brief return secondary structure string
	inline std::string const & ss() const { return ss_; }

	/// @brief return secondary structure string
	inline utility::vector1< std::string > const & abego() const { return abego_; }

	/// @brief given an input stream, substitute all variables
	/// @details variables are of the form: %%SEGMENTNAME#residue%%
	/// SEGMENTNAME = name of the segment
	/// residue = local residue number within the segment
	/// The substituted value will be an core::Size corresponding to the pose residue
	std::string substitute_variables( std::istream & input ) const;

	/// @brief returns true of the last residue of segment1 contains a covalent bond to the first residue of segment2
	bool polymer_bond_exists( std::string const & segment1, std::string const & segment2 ) const;

	/// @brief checks the permutation for internal consistency
	void check_consistency() const;
	void check_consistency( core::pose::Pose const & pose ) const;

	/// @brief for output
	friend std::ostream & operator<<( std::ostream & os, StructureData const & perm );

	/// @brief sets the pose and does checks to ensure data is consistent
	void set_pose( core::pose::Pose const & new_pose );

	/// @brief sets the pose and does checks to ensure data is consistent
	void set_pose( core::pose::PoseOP new_pose );

	/// @brief sets the fold tree in the pose
	void set_fold_tree( core::kinematics::FoldTree const & ft );

	/// @brief marks the given segments as covanlently connected
	void mark_connected(
		std::string const & lower_seg,
		std::string const & upper_seg );

	/// @brief unmarks the given segments as covalently connected
	void mark_disconnected(
		std::string const & seg1,
		std::string const & seg2 );

	// pose modification methods

	/// @brief returns the chain number of the given residue. Pose MUST be set.
	core::Size chain( core::Size const resid ) const;

	/// @brief determines pose chains based on the termini in the pose
	void chains_from_termini();

	/// @brief adds upper terminal variant to a residue
	void add_upper_terminus_variant_type( core::Size const resi );

	/// @brief adds lower terminal variant to a residue
	void add_lower_terminus_variant_type( core::Size const resi );

	/// @brief removes upper terminal variant from a residue
	void remove_upper_terminus_variant_type( core::Size const resi );

	/// @brief removes lower terminal variant from a residue
	void remove_lower_terminus_variant_type( core::Size const resi );

	/// @brief creates a new jump and cutpoint at the given residue
	/// returns the number of the new jump
	int new_jump_and_cutpoint( protocols::loops::Loop const & loop, core::Size const loop_overlap );

	/// @brief removes jump and cutpoint between the two segments to create a single polymer chain
	void delete_jump_and_intervening_cutpoint( std::string const & segment1, std::string const & segment2 );

	/// @brief given a jump and a cutpoint residue, move the jump around the cut residue and delete to form a single edge
	void delete_jump_and_intervening_cutpoint( int const jnum, core::Size const cut_resi1, core::Size const cut_resi2 );

	/// @brief applies an arbitrary mover to the contained pose -- the mover must NOT change pose length or numbering!
	void apply_mover( protocols::moves::Mover & mover );

	/// @brief applies an arbitrary mover to the contained pose -- the mover must NOT change pose length or numbering!
	void apply_mover( protocols::moves::MoverOP mover );

	/// @brief switches residue type set of the contained pose
	void switch_residue_type_set( std::string const & typeset );

	/// @brief removes constraints added by the given RCG
	void remove_constraints_from_pose( protocols::forge::remodel::RemodelConstraintGeneratorOP rcg );

	/// @brief sets jump with index jumpidx to the given datastructure
	void set_jump( int const jumpidx, core::kinematics::Jump const & j );

	/// @brief sets phi torsion for a specific residue
	void set_phi( core::Size const seqpos, core::Real const phi_val );

	/// @brief sets psi torsion for a specific residue
	void set_psi( core::Size const seqpos, core::Real const psi_val );

	/// @brief sets omega for a specific residue
	void set_omega( core::Size const seqpos, core::Real const omega_val );

	/// @brief detects and sets disulfides
	void detect_disulfides( core::scoring::ScoreFunctionOP sfx );

	/// @brief sets bond length in pose
	void set_bond_length(
		std::string const & segmentname,
		core::Size const resid,
		std::string const & atom1,
		std::string const & atom2,
		core::Real const newlength );

	/// @brief expands the segment so that the trailing pad residue(s) become part of the segment
	void engulf_leading_residues( std::string const & seg );

	/// @brief expands the segment so that the trailing pad residue(s) become part of the segment
	void engulf_trailing_residues( std::string const & seg );

	/// @brief deletes the residues between the segment N terminus and the N anchor point
	void delete_leading_residues( std::string const & seg );

	/// @brief deletes the residues between the segment C terminus and the C anchor point
	void delete_trailing_residues( std::string const & seg );

	/// @brief replaces one residue with another
	void replace_residue(
		std::string const & target_segment,
		core::Size const target_res,
		core::conformation::Residue const & res_in );

	/// @brief replaces one residue with another
	void replace_residue(
		core::Size const resnum,
		core::conformation::Residue const & res_in,
		bool const orient_bb );

	/// @brief chooses a new movable group which doesn't conflict with existing ones
	core::Size choose_new_movable_group() const;

	/// @brief sets movable group of a segment
	void set_movable_group( std::string const & id, core::Size const mg );

	/// @brief renumbers movable group "oldg" to have new number "newg"
	void renumber_movable_group( core::Size const oldg, core::Size const newg );

	/// @brief updates numbering based on the saved order of Segment objects
	void update_numbering();
	void update_covalent_bonds_in_pose();

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

protected:
	/// @brief returns constant list of residue ranges
	inline SegmentMap const & segments() const { return segments_; }

	/// @brief finds a segment in the segment_order list and returns an iterator to it
	StringList::iterator find_segment( std::string const & segname );

	/// @brief performs dfs in lower direction looking for termini
	std::string find_lower_terminus( std::set< std::string > & visited, std::string const & seg ) const;

	/// @brief performs dfs in upper direction looking for termini
	std::string find_upper_terminus( std::set< std::string > & visited, std::string const & seg ) const;

	/// @brief aligns two residues so their backbones are completely superimposed
	void align_residues(
		core::Size const jump_idx,
		core::Size const align_res_target,
		core::Size const align_res_movable,
		core::Size const res_with_torsions );

	void add_covalent_bond(
		core::Size const res1, std::string const & atom1,
		core::Size const res2, std::string const & atom2 );

	void add_covalent_bond(
		std::string const & seg1, core::Size const res1, std::string const & atom1,
		std::string const & seg2, core::Size const res2, std::string const & atom2 );

	void add_covalent_bond( BondInfo const & bi );

	void declare_covalent_bond_in_pose(
		core::Size const res1, std::string const & atom1,
		core::Size const res2, std::string const & atom2 );

	/// @brief moves a segment of the pose such that the segment from start2<=res<=end2 is moved so that it starts at end1+1
	/// returns the jump number that is to be ignored in fold tree searches
	/// WARNING: all parameters are by REFERENCE
	void move_segment_in_pose(
		core::Size start1,
		core::Size end1,
		core::Size start2,
		core::Size end2 );

	/// @brief copies and inserts a segment into the permutation after the given segment
	void insert_after_residue_in_pose(
		core::Size segment1_start,
		core::Size segment1_end,
		core::Size segment2_start,
		core::Size segment2_end );

	/// @brief deletes the given residues from the pose
	void delete_residues_in_pose(
		core::Size start,
		core::Size end );

	/// @brief deletes segment in SD without touching pose -- doesn't update numbering
	void delete_segment_nopose( std::string const & seg_val, SegmentMap::iterator r );

	/// @brief add lower cutpoint to residue cut and upper cutpoint to residue cut+1
	void add_cutpoint_variants( core::Size const cut_res );

	/// @brief removes cutpoint variants from residues cut and cut+1
	void remove_cutpoint_variants( core::Size const cut_res );

	/// @brief renames a residue segment and updates all connections
	void add_prefix_to_segments( std::string const & prefix, char const delimeter );

	/// @brief adds a remark to remarks object
	void add_perm_remark( core::pose::Remarks & remarks, std::string const & rem_value ) const;

	/// @brief Saves remarks of the given pose into the pose's datacache -- changes enzdes residues to segment name/number
	void save_remarks_to_datacache( core::pose::Remarks const & remarks );

	/// @brief moves around jumps so that movable groups will all move together during folding
	void consolidate_movable_groups( core::pose::PoseOP pose, utility::vector1< std::string > const & root_segments );

	// internal bookkeeping functions
private:
	/// @brief updates movable group numbering after a deletion -- deleted mg is passed
	void update_movable_groups_after_deletion( core::Size const mg_old );

	/// @brief checks pose vs. StructureData info
	/// @throw EXCN_PoseInconsistent if things don't match
	void check_pose( core::pose::Pose const & pose ) const;

	/// @brief checks residues in SD -- makes sure everything is sequential and accounted for
	void check_residues() const;

	/// @brief checks chain termini in pose vs SD
	/// @throws EXCN_PoseInconsistent if there is a problem
	void check_improper_termini( core::pose::Pose const & pose ) const;
	void check_chain_beginnings( core::pose::Pose const & pose ) const;
	void check_chain_endings( core::pose::Pose const & pose ) const;

	/// @brief check pose movable groups vs SD
	void check_movable_groups() const;

	// member variables
private:
	// permutation itself is defined by a pose
	core::pose::PoseOP pose_;
	// for identification with a segment
	std::string id_;
	// secondary structure string
	std::string ss_;
	// abego string
	utility::vector1< std::string > abego_;
	// length of the pose, including n-, c-terminal loop residues
	core::Size pose_length_;
	// usable length of the pose, not including n-, c-terminal loops
	core::Size length_;
	// sub-permutations for compound objects
	std::map< std::string, int > data_int_;
	// arbitrary real number data that can be set by movers
	std::map< std::string, core::Real > data_real_;
	// arbitrary string data that can be set by movers
	std::map< std::string, std::string > data_str_;
	// map of segment to residue [start, end]
	SegmentMap segments_;
	// names given to special single residues -- value stored is segment name + intra-segment residue number
	std::map< std::string, Alias > aliases_;
	// non-canonical covalent bonds
	utility::vector1< BondInfo > covalent_bonds_;
	// segments listed in order
	StringList segment_order_;
};

class SingleChainStructureData : public StructureData {
public:
	SingleChainStructureData( std::string const & id_val );

	SingleChainStructureData(
		std::string const & id_val,
		core::Size const length_val,
		core::Size const pose_len_val,
		bool const is_loop,
		std::string const & ss_val,
		utility::vector1< std::string > const & abego_val );

	virtual ~SingleChainStructureData();

	virtual bool is_multi() const { return false; }

	virtual StructureDataOP clone() const { return StructureDataOP( new SingleChainStructureData( *this ) ); }

	virtual StructureDataOP fresh_instance() const { return StructureDataOP( new SingleChainStructureData( id() ) ); }

private:
};

class MultiChainStructureData : public StructureData {
public:
	MultiChainStructureData( std::string const & id_val ) :
		StructureData( id_val )
	{
	}
	virtual ~MultiChainStructureData();

	virtual bool is_multi() const { return true; }

	virtual StructureDataOP fresh_instance() const { return StructureDataOP( new MultiChainStructureData( id() ) ); }

	virtual StructureDataOP clone() const { return StructureDataOP( new MultiChainStructureData( *this ) ); }

private:
};

class EXCN_PoseInconsistent : public utility::excn::EXCN_BadInput {
public:
	EXCN_PoseInconsistent( std::string const & msg ):
		utility::excn::EXCN_BadInput( msg ) {}
};

/// dump contents of residues map
std::ostream & operator<<( std::ostream & os, SegmentMap const & resmap );
std::ostream & operator<<( std::ostream & os, Alias const & alias );

} // components
} // denovo_design
} // protocols

#endif
