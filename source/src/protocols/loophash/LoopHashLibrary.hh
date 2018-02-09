// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashLibrary.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_LoopHashLibrary_hh
#define INCLUDED_protocols_loophash_LoopHashLibrary_hh

#include <protocols/loophash/LoopHashLibrary.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loophash/LoopHashMap.hh>
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

	void extract_data_from_pose( core::pose::Pose& pose, core::Size nres, protocols::frag_picker::VallChunkOP chunk = nullptr );

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

	void apply( core::pose::Pose& pose ) override;

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new LoopHashLibrary( *this ) );
	}


	std::string get_name() const override {
		return "LoopHashLibrary";
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return protocols::moves::MoverOP( new LoopHashLibrary );
	}


	void setup_hash_maps();

	// simple accessors
	LoopHashMap const & gethash( core::Size size ) const;

	// Requires non-const access to the loop hash library -- note that this is not
	// an option if the LoopHashLibrary is being managed by the ResourceManager
	LoopHashMap & gethash( core::Size size );


  /// @brief Return a list of the loop lengths present in the library.
	std::vector< core::Size > const & hash_sizes() const { return hash_sizes_; }

	BackboneDB const & backbone_database() const { return bbdb_; }

	std::pair< core::Size, core::Size > loopdb_range() const { return loopdb_range_; }

	/// @details Only support writing to text, always include extra data
	void save_db() const;

	/// @details Only support reading from text, extra data is mandatory.  Used 
	/// when created merged text db.
	void load_db();

	/// @details Only support reading from merged text.  Extra data is optional 
	/// and handled by extra_.
	void load_mergeddb();

	// For backwards-compability
	//void load_legacydb();

	void delete_db();

	/// @brief Take two LoopHash libraries and merge them.
	/// @details Throws out non-unique fragments defined as follows: if a new 
	/// fragment is in a bucket with a key also in the master library, and if the 
	/// RMS between the new fragment and the members of the master bucket is 
	/// lower than the @p rms_cutoff, it is non-unique.  Thus, an @p rms_cutoff 
	/// of zero will not throw away any fragments.
	void merge( LoopHashLibraryOP second_lib, utility::vector1< core::Real> rms_cutoffs );

	/// @brief Take two backbone databases and merge them.
	/// @details Just takes Data structs from one and copies to the other
	bool merge_bbdb( const BackboneDB & second_bbdb, core::Size & index_offset );

	/// @details This function destroys the backbone_index_map reference, so 
	/// sorting is only done when the hash will not longer be used.  Required 
	/// after merging, before saving so that merged db can be loaded as slices.
	void sort();

	void create_db();
	void set_create_db( bool setting = true ){ create_db_ = setting; }
	void set_db_path( std::string setting ){ db_path_ = setting; }

	void graft_loop(
		const core::pose::Pose& src_pose,
		core::pose::Pose& tgt_pose,
		protocols::loops::Loop myloop
	);

	// setup scorefunctions for
	void set_default_score_functions();


	void mem_foot_print() const;

	bool get_extra() const { return extra_; }

private:

	/// @brief The backbone library for this HashLibrary (the actual data)
	BackboneDB bbdb_;

	/// @brief A map of Hashmaps, each for a different loop size.
	std::map< core::Size, LoopHashMap > hash_;

	/// @brief Kind-of a redundant store of sizes (remove this?)
	std::vector < core::Size > hash_sizes_;

	/// @brief Path to database.
	std::string db_path_;

	/// @brief The number of partitions the database is split into.
	core::Size num_partitions_;

	/// @brief Which partition slice is assigned to this library.
	core::Size assigned_num_;

	/// @brief Need so we can set the database string once.
	std::string assigned_string_;

	/// @brief Whether this database holds extra data in the backbone database or 
	/// not.
	bool extra_;

	/// @brief The proteins of the backbone database that are loaded (when 
	/// loading a merged database, otherwise (0,0))
	std::pair< core::Size, core::Size > loopdb_range_;

	// Some basic flags
	bool do_sanity_check_;
	bool create_db_;

	// Used for grafting - ultimately to move into a separate Mover class.

	core::scoring::ScoreFunctionOP scorefxn_rama_cst_;
	core::scoring::ScoreFunctionOP scorefxn_cen_cst_;

	core::optimization::MinimizerOptions options_;
	core::optimization::MinimizerOptions options2_;
};

} // namespace loops
} // namespace protocols

#endif
