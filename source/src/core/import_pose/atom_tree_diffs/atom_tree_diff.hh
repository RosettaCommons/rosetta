// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/atom_tree_diffs/atom_tree_diff.hh
///
/// @brief  Silent-file format based on "diffs" of AtomTree DOFs
/// @author Ian W. Davis

#ifndef INCLUDED_core_import_pose_atom_tree_diffs_atom_tree_diff_hh
#define INCLUDED_core_import_pose_atom_tree_diffs_atom_tree_diff_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#ifdef WIN32
#include <core/pose/Pose.hh> // WIN32 INCLUDE
#endif
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <fstream>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace atom_tree_diffs {


class AtomTreeDiff; // fwd declaration
typedef utility::pointer::shared_ptr< AtomTreeDiff > AtomTreeDiffOP;
typedef utility::pointer::shared_ptr< AtomTreeDiff const > AtomTreeDiffCOP;

typedef std::pair<std::string, core::Real> ScorePair;
typedef std::map< std::string, core::Real > Scores;
typedef std::pair< std::string, Scores > ScoresPair;

//typedef std::map< std::string, Scores > ScoresMap;
/// Just like ScoresMap, but can be sorted into some particular order *by scores*.
/// Maps can only be sorted by key, which here is just a pose tag.
typedef utility::vector1< ScoresPair > ScoresPairList;

typedef std::pair<std::string, int> RefTag;
typedef std::map<std::string, int> RefTags;

typedef std::pair<std::string, Size> TagScorePair;
typedef std::map<std::string, Size> TagScoreMap;


/// @brief An object wrapper for reading atom_tree_diff files,
/// complete with embedded reference structures.
///
/// @details Only works with uncompressed files, because we have to be able to
/// to random access (seekg) to pull out single structures in random order.
///
class AtomTreeDiff : public utility::pointer::ReferenceCount
{
public:

	AtomTreeDiff();
	AtomTreeDiff(std::string filename);
	virtual ~AtomTreeDiff();

	/// @brief returns true if a reference struct with the given tag is present
	bool has_ref_pose(std::string const & tag) const;

	/// @brief True if a (non-reference) structure with the given tag is present in the file
	bool has_tag(std::string const & tag) const;

	bool has_ref_tag(std::string const & tag) const;

	/// @brief Return list of (pose tag, score sets) pairs for all poses, in file order.
	ScoresPairList const & scores() const { return scores_; }

	/// @brief Utility function for selecting subsets of structures by score.
	static void sort_by(std::string const & score_name, ScoresPairList & scores, bool descending=false);

	/// @brief Reads the pose data from file and reconstructs the complete pose.
	void read_pose(std::string const & tag, core::pose::Pose & pose_out);

	/// @brief Reads the pose data from file and reconstructs the complete pose, using the supplied reference pose.
	void read_pose(std::string const & tag, core::pose::Pose & pose_out, core::pose::Pose const & ref_pose);

	void read_file(std::string filename);

	/// @brief Returns the default reference pose for the given tag.  Fails if none is available.
	core::pose::PoseCOP ref_pose_for(std::string const & tag);

	/// @brief Allows access to and mutation of (!) the references poses stored in this file.  Use with caution.
	/// @details This exists to allow setup on stored reference poses for properties that don't get saved/restored
	/// in PDB format, like covalent constraints for enzyme design.
	core::pose::PoseOPs const & all_ref_poses() const
	{ return unique_ref_poses_; }

	RefTags const & get_ref_tags() const {return ref_tags_;}

	TagScoreMap const & get_tag_score_map() const{return tag_idx_;}

	int count(std::string const & tag) const{
		std::map<std::string, int>::const_iterator const iter= ref_tags_.find(tag);
		return iter->second;
	}

private:
	AtomTreeDiff(AtomTreeDiff const &);

	/// The list of (tag,scores) pairs for the diffed structures, for efficient selection of subsets.
	/// (E.g. the top 5% by total score.) Memory cost is probably ~1 kb / structure.
	ScoresPairList scores_;

	/// Maps reference tags to how many associated atom_tree_diff structs exist...
	/// int is needed here because there are 0 atom tree diffs to start with
	RefTags ref_tags_;

	/// Maps tags to indices in scores_
	TagScoreMap tag_idx_;

	/// The map of (tag, file position) for getting random access to the structural data.
	std::map< std::string, long > offsets_;

	/// The map of (tag, reference pose) for decoding the diffs.  Some could be null.
	std::map< std::string, core::pose::PoseOP > ref_poses_;

	/// All references poses from the file, one copy each for convenience.
	core::pose::PoseOPs unique_ref_poses_;

	// The file being read.
	std::ifstream in_;

	mutable bool file_read_;

}; // class AtomTreeDiff

/// @brief Helper function for writing entries -- not usually called by clients.
void dump_score_line(
	std::ostream & out,
	std::string const & pose_tag,
	std::map< std::string, core::Real > const & scores
);

/// @brief Helper function for writing entries -- not usually called by clients.
/* void dump_score_line(
	std::ostream & out,
	std::string const & pose_tag
); */

/// @brief Embeds a reference pose as PDB coords + foldtree; will be used for reconstructing subsequent diffs.
void dump_reference_pose(
	std::ostream & out,
	std::string const & pose_tag,
	std::map< std::string, core::Real > const & scores,
	core::pose::Pose const & pose
);

/// @brief Encodes pose relative to ref_pose by noting which atom_tree DOFs are different.
void dump_atom_tree_diff(
	std::ostream & out,
	std::string const & pose_tag,
	std::map< std::string, core::Real > const & scores,
	core::pose::Pose const & ref_pose_in,
	core::pose::Pose const & pose,
	int bb_precision = 6,
	int sc_precision = 4,
	int bondlen_precision = 2
);

/// @brief Gets next tag and scores from the stream, or returns false if none.
/// Call this to find desired structure, then call pose_from_atom_tree_diff().
bool header_from_atom_tree_diff(
	std::istream & in,
	std::string & pose_tag_out,
	std::map< std::string, core::Real > & scores_out
);


/// @brief Sets pose = ref_pose and then starts modifying DOFs in pose to recreate a saved structure.
/// Call after header_from_atom_tree_diff().  Returns false on error.
bool pose_from_atom_tree_diff(
	std::istream & in,
	core::pose::Pose const & ref_pose,
	core::pose::Pose & pose
);


/// @brief Helper for dump_atom_tree_diff(), fills map with weighted score terms.
void map_of_weighted_scores(
	core::pose::Pose & pose, //< pose is not modified but scoring is a non-const op
	core::scoring::ScoreFunction const & sfxn,
	std::map< std::string, core::Real > & scores_out
);


/// @brief For use in deciding how many digits of precision you need when diffing an atom tree.
void rms_error_with_noise(
	core::pose::Pose const & ref_pose,
	int bb_precision = 6,
	int sc_precision = 4
);

/// @brief Test if given file is an atom_tree_diff
bool file_is_atom_tree_diff( std::string const & filename );

/// @brief Test if given stream is an atom_tree_diff
/// @details If everything goes right, after the call, the read position should be at the same place it was to start with
bool file_is_atom_tree_diff( std::istream & in );

} // silent
} // io
} // core

#endif // INCLUDED_core_import_pose_silent_atom_tree_diff_HH
