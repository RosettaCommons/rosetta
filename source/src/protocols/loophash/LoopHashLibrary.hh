// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/LoopHashLibrary.hh
/// @brief
/// @author Mike Tyka



#ifndef INCLUDED_protocols_loophash_LoopHashLibrary_hh
#define INCLUDED_protocols_loophash_LoopHashLibrary_hh

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
// AUTO-REMOVED #include <protocols/loophash/LoopHashSampler.fwd.hh>
// AUTO-REMOVED #include <protocols/loophash/LocalInserter.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/loops/Loop.hh>
#include <protocols/loophash/LoopHashMap.hh>
// AUTO-REMOVED #include <protocols/frag_picker/VallChunk.hh>
// AUTO-REMOVED #include <protocols/frag_picker/VallProvider.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <string>
#include <vector>
#include <map>

#include <protocols/loops/Loop.fwd.hh>


namespace protocols {
namespace loophash {


class LoopHashLibrary: public protocols::moves::Mover {

public:

	LoopHashLibrary( const utility::vector1< core::Size > &init_sizes = utility::vector1< core::Size >(), const core::Size num_partitions = 1, const core::Size assigned_num = 0);

	void extract_data_from_pose( core::pose::Pose& pose, core::Size nres, protocols::frag_picker::VallChunkOP chunk = NULL );

	void extract_data_from_pose( core::pose::Pose& pose );

	bool test_saving_library( core::pose::Pose, core::Size ir, bool deposit );

	void test_loop_sample( core::pose::Pose& pose, core::Size nres );

	void apply_random(
		core::pose::Pose& pose,
		core::Size &fir,
		core::Size &fjr,
		core::Real min_rms,
		core::Real max_rms
	);

	void get_all(
		core::pose::Pose& pose,
    std::vector< core::io::silent::SilentStructOP > &lib_structs,

    core::Size start_res = 1,
    core::Size stop_res  = 0,

    core::Real min_bbrms = 0.0,
		core::Real max_bbrms = 100000.0,
		core::Real min_rms   = 0.0,
		core::Real max_rms   = 100.0
	);

	virtual void apply( core::pose::Pose& pose );

  virtual protocols::moves::MoverOP clone() const {
		return new LoopHashLibrary( *this );
	}


	virtual std::string get_name() const {
		return "LoopHashLibrary";
	}

	virtual	protocols::moves::MoverOP	fresh_instance() const {
		return new LoopHashLibrary;
	}


	void setup_hash_maps();

  // simple accessors
	LoopHashMap & gethash( core::Size size );


	const std::vector < core::Size > & hash_sizes() const { return hash_sizes_; }

	const BackboneDB & backbone_database() const { return bbdb_; }

	const std::pair< core::Size, core::Size > loopdb_range() { return loopdb_range_; }

	// Only support writing to text, always include extra data
	void save_db();

	// Only support reading from text, extra data is mandatory
	// used when created merged text db
	void load_db();

	// Only support reading from merged text
	// Extra data is optional and handled by extra_
	void load_mergeddb();

	// For backwards-compability
	//void load_legacydb();

	void delete_db();

	// Takes two LoopHashLibrary and merges them
	// Throws out nonunique fragments defined as follows
	// If a new frag is in a bucket with a key also in the master library, and if the rms between
	// the new frag and the members of the master bucket is lower than the rms_cutoff, it is nonunique
	// Thus, an rms_cutoff of zero will not throw away any fragments
	void merge( LoopHashLibraryOP second_lib, utility::vector1< core::Real> rms_cutoffs );

	// Takes two bbdb and merges them
	// Just takes Data structs from one and copies to the other
	bool merge_bbdb( const BackboneDB & second_bbdb, core::Size & index_offset );

	// This function destroys the backbone_index_map reference, so sorting is only done when the hash will not longer be used
	// Required after merging, before saving so that merged db can be loaded as slices
	void sort();

	void create_db();
	void set_create_db( bool setting = true ){ create_db_ = setting; }
	void set_db_path( std::string setting ){ db_path_ =  setting; }

  void graft_loop(
    const core::pose::Pose& src_pose,
    core::pose::Pose& tgt_pose,
    protocols::loops::Loop myloop
  );

  // setup scorefunctions for
  void set_default_score_functions();


	void mem_foot_print();

	bool get_extra() const { return extra_; }

private:

  // The backbone library for this HashLibrary (the actual data)
	BackboneDB bbdb_;

  // A map of Hashmaps, each for a different loop size
	std::map    < core::Size, LoopHashMap > hash_;

  // Kindof a redundant store of sizes - remove this ?
	std::vector < core::Size > hash_sizes_;

  // Path to db
    std::string db_path_;

  // the number of partitions the database is split into
    core::Size num_partitions_;

  // which partition slice is assigned to this Library
    core::Size assigned_num_;

  // Need so we can set the db string once
    std::string assigned_string_;

  // Whether this database holds extra data in the bbdb or not
    bool extra_;

	// The proteins of the backbone db that are loaded (when loading a merged db, otherwise (0,0)
		std::pair< core::Size, core::Size > loopdb_range_;

  // Some basic flags
	bool do_sanity_check_;
	bool create_db_;

  // Used for grafting - ultimately to move into a seperate Mover class.

	core::scoring::ScoreFunctionOP scorefxn_rama_cst;
	core::scoring::ScoreFunctionOP scorefxn_cen_cst;

	core::optimization::MinimizerOptions options;
	core::optimization::MinimizerOptions options2;
};






//class LoopHashSampler {
//  public:
//
//  LoopHashSampler(
//    LoopHashLibraryOP library,
//    LocalInserterOP inserter
//  ):
//    library_(library),
//    inserter_(inserter),
//    start_res_ ( 1 ),
//    stop_res_  ( 0 ),
//    min_bbrms_ ( 0.0 ),
//		max_bbrms_ ( 100000.0 ),
//		min_rms_   ( 0.0 ),
//		max_rms_   ( 100.0 )
//  {
//
//  }
//
//  // @brief create a set of structures for a the given range of residues and other parameters
//  void build_structures(
//		core::pose::Pose& start_pose,
//    std::vector< core::io::silent::SilentStructOP > &lib_structs
//	);
//
//
//
//  private:
//    // pointer to the library used for insertion
//    LoopHashLibraryOP library_;
//
//    // pointer to the insertion functor which provides the peptide insertion facility
//    LocalInserterOP inserter_;
//
//    // parameters for insertion positions
//    core::Size start_res_;
//    core::Size stop_res_ ;
//    core::Real min_bbrms_;
//		core::Real max_bbrms_;
//		core::Real min_rms_  ;
//		core::Real max_rms_  ;
//
//};





} // namespace loops
} // namespace protocols



#endif



