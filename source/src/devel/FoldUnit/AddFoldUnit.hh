// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddFoldUnitMover.hh

#ifndef INCLUDED_protocols_splice_AddFoldUnitMover_hh
#define INCLUDED_protocols_splice_AddFoldUnitMover_hh

#include <protocols/splice/AddFoldUnit.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/splice/Splice.hh>
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map.hpp>

// C++ Headers
namespace protocols {
namespace splice {

/// @brief an object to manage the torsion/sequence data and to perform individual moves, such as introducing a new fragment
// The main method is pose_from_fragment_info which creates a pose directly from the PoseFragmentInfo
class FoldUnitUtils : public utility::pointer::ReferenceCount {
public:
	FoldUnitUtils() : fragment_dbase_( "" ), pair_dbase_( "" ), overlap_length_( 10 ), n_term_entry_( 0 ), c_term_entry_( 0 ), min_overlap_( 0 ), max_length_( 0 ), max_rmsd_( 0.0 ) { legal_bbdofs_.clear(); entry_pairs_.clear(); overlap_length_.clear(); overlap_rmsd_.clear(); };
	~FoldUnitUtils() override;
	bool read_dbase(); //return whether read is successful
	void fragment_dbase( std::string const & dbase ){ fragment_dbase_ = dbase; }
	std::string fragment_dbase() const{ return fragment_dbase_; }
	void pair_dbase( std::string const & dbase ){ pair_dbase_ = dbase; }
	std::string pair_dbase() const{ return pair_dbase_; }

	//  utility::vector1< ResidueBBDofs > & bbdofs(){ return bbdofs_; }// useful to generate bbdofs subsets
	utility::vector1< ResidueBBDofs > const & bbdofs() const { return bbdofs_; }// useful to generate bbdofs subsets
	utility::vector1< core::Size > & overlap_length() { return overlap_length_; }
	utility::vector1< core::Size > overlap_length() const { return overlap_length_; }
	utility::vector1< core::Real > & overlap_rmsd() { return overlap_rmsd_; }

	bool fragment_compatibility_check( Size const i, Size const j, utility::vector1< std::pair< core::Size, core::Size > >::const_iterator & it, core::Real const max_rmsd = 10.0 ) const;
	/// methods below add or replace fragments without testing that the fragments are compatible, so they should be used with care
	void add_fragment_to_pose( core::pose::Pose & pose, PoseFragmentInfo & pose_fragment_info, core::Size const entry, bool c_term/*or n-term*/ ) const;
	void replace_fragment_in_pose( core::pose::Pose & pose, PoseFragmentInfo & pose_fragment_info, core::Size const entry, core::Size const fragment_num );

	core::Size n_term_entry() const{ return n_term_entry_; }
	void n_term_entry( core::Size const s ){ n_term_entry_ = s ;}
	core::Size c_term_entry() const{ return c_term_entry_; }
	void c_term_entry( core::Size const s ){ c_term_entry_ = s ;}
	core::Size min_overlap() const{ return min_overlap_; }
	void min_overlap( core::Size const s ){ min_overlap_ = s ;}
	core::Size max_length() const{ return max_length_; }
	void max_length( core::Size const s ){ max_length_ = s ;}
	core::Real max_rmsd() const{ return max_rmsd_; }
	void max_rmsd( core::Real const s ){ max_rmsd_ = s ;}
	core::Size pair_index( core::Size const i, core::Size const j ) const; //return the index into the vector1 objects of pair i,j

	utility::vector1< bool > legal_bbdofs() const{ return legal_bbdofs_; }
	utility::vector1< core::Size > entry_subset() const; // find the subset of entries that match all of the criteria (n_term_entry, c_term_entry, etc.)
	utility::vector1< core::Size > entry_subset_slow() const; // uses the entry_pairs vector. Low memory but very slow
	void pose_from_fragment_info( core::pose::Pose & pose, PoseFragmentInfo const & pose_fragment_info ) const; // generate a new pose from the PFI

private:

	std::string fragment_dbase_, pair_dbase_;
	utility::vector1< ResidueBBDofs > bbdofs_;// from fragment database
	utility::vector1< bool > legal_bbdofs_; // some entries are corrupted but need to be kept as placeholders. legal_bbdofs_ tells us whether the entry is valid for use.
	utility::vector1< std::pair< core::Size/*i*/, core::Size/*j*/ > > entry_pairs_; // is this connection valid? i->j
	utility::vector1< core::Size > overlap_length_; //what is the overlap between i,j fragments?
	utility::vector1< core::Real > overlap_rmsd_; // what's the rmsd over this overlap (rmsd is computed in degrees over dihedral angles in MatLab

	/// values used to generate subsets of entries that match the criteria
	core::Size n_term_entry_, c_term_entry_;// dflt 0 0; entries on both sides to restrict selection
	core::Size min_overlap_, max_length_;
	core::Real max_rmsd_;
	boost::unordered_multimap< core::Size, core::Size > entry_pairs_quick_access_N_C_, entry_pairs_quick_access_C_N_;

};

/// @brief manage pose-fragmnet information.
class PoseFragmentInfo : public utility::pointer::ReferenceCount{
public:
	PoseFragmentInfo(){ fragment_map_.clear(); };
	~PoseFragmentInfo() override;
	void load_fragment_info_from_pose( core::pose::Pose const & pose ); // set the object according to comments in the pose
	void set_fragment_info_in_pose( core::pose::Pose & pose ) const; // set the pose comments according to the object
	void fragment_start_end( FoldUnitUtils const & fuu, core::Size const fragment, core::Size & start, core::Size & end ); // provide start/end residues for fragment #fragment
	core::Size & operator[]( core::Size const s );
	core::Size operator[] ( core::Size const s ) const;
	core::Size size() const{ return fragment_map_.size(); }
	void add_pair( core::Size const segment_num, core::Size const dbase_entry ){ fragment_map_[ segment_num ] = dbase_entry; }
	std::map< core::Size, core::Size > fragment_map() const { return fragment_map_; }
private:
	std::map< core::Size/*pose segment #*/, core::Size/*fragment_dbase entry*/ > fragment_map_;
};

class AddFoldUnitMover : public protocols::moves::Mover {
public:
	AddFoldUnitMover();
	~AddFoldUnitMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string fragment_dbase() const{ return fragment_dbase_; }
	void fragment_dbase( std::string const & s ){ fragment_dbase_ = s; }
	std::string pair_dbase() const{ return pair_dbase_; }
	void pair_dbase( std::string const & s ){ pair_dbase_ = s; }

	core::Size max_length() const{ return max_length_; }
	void max_length( core::Size const s ){ max_length_ = s; }
	char terminus() const { return terminus_; }
	void terminus( char const c ){ terminus_ = c; }
	core::Size replace_fragment_segment_number() const{ return replace_fragment_segment_number_; }
	void replace_fragment_segment_number( core::Size const s ){ replace_fragment_segment_number_ = s; }
	bool replace_fragment() const{ return replace_fragment_; }
	void replace_fragment( bool const b ){ replace_fragment_ = b; }


	FoldUnitUtilsOP fold_unit_utils() const{ return fold_unit_utils_; }
	void fold_unit_utils( FoldUnitUtilsOP f ){ fold_unit_utils_ = f; }
	void min_overlap( core::Size const m ){ min_overlap_ = m; }
	core::Size min_overlap() const{ return min_overlap_; }
	void max_rmsd( core::Real const r ){ max_rmsd_ = r; }
	core::Real max_rmsd() const{ return max_rmsd_; }
	void max_segments( core::Size const m ){ max_segments_ = m; }
	core::Size max_segments() const{ return max_segments_; }


private:
	std::string fragment_dbase_, pair_dbase_;
	core::Size max_length_, min_overlap_; //dflt 40, 5 ; largest fragment to insert
	core::Real max_rmsd_; //dflt 10o.
	char terminus_; //dflt 'b'; add to C-terminus or N-terminus stochastically? 'c' and 'n' can be used to specify one or the other
	bool replace_fragment_; // dflt false; replace a segment rather than adding one?
	core::Size replace_fragment_segment_number_; // dflt 0 ; if 0, randomize segment, else replace the specified segment
	FoldUnitUtilsOP fold_unit_utils_;
	core::Size max_segments_; //dflt 5; how many segments to add?
};

/// @brief silly mover that nukes out all information from the pose. RosettaScripts can't start without a pose...
class StartFreshMover : public protocols::moves::Mover {
public:
	StartFreshMover();
	~StartFreshMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string residue_type_set() const{ return residue_type_set_; }
	void residue_type_set( std::string const & s ){ residue_type_set_ = s; }
private:
	std::string residue_type_set_; // dflt "CENTROID"; centroid or fullatom?
};
} // simple_moves
} // protocols

#endif
