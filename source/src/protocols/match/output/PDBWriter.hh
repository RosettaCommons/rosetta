// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/PDBWriter.hh
/// @brief
/// @author Florian Richter (floric@u.washington.edu) june 2009

#ifndef INCLUDED_protocols_match_output_PDBWriter_hh
#define INCLUDED_protocols_match_output_PDBWriter_hh

// Unit headers
#include <protocols/match/output/PDBWriter.fwd.hh>

// Package Headers
#include <core/conformation/Residue.fwd.hh>
#include <protocols/match/output/MatchEvaluator.fwd.hh>
#include <protocols/match/output/MatchScoreWriter.fwd.hh>

#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/MatcherTask.fwd.hh>
#include <protocols/match/output/MatchGrouper.fwd.hh>
#include <protocols/match/output/OutputWriter.hh>
#include <protocols/match/output/UpstreamHitCacher.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>
#include <string>

#include <core/pose/Pose.fwd.hh>
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>
#include <utility/vector1.hh>
#include <iterator>


namespace protocols {
namespace match {
namespace output {

class PDBWriter : public OutputWriter
{

public:
	typedef std::pair< core::Size, core::Size > SizePair;

	typedef OutputWriter Parent;

public:

	PDBWriter();

	virtual ~PDBWriter();

	virtual
	void
	prepare_for_output_writing();

	virtual
	void
	record_match( match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer );

	/// @brief evaluator and score writer are not passed in because single-downstream-position match
	///currently have no way of being evaluated
	virtual
	void
	record_match( match_dspos1 const & m );

	//core::pose::PoseCOP
	//create_output_pose_from_dspos1_match(
	// match_dspos1 const & m );

	void
	set_coordinate_cacher( UpstreamHitCacherOP );

	void
	set_prefix( std::string const & prefix );

	void
	initialize_from_matcher_task(
		MatcherTaskCOP mtask
	);

	void
	set_downstream_builder(
		Size geomcst_id,
		downstream::DownstreamBuilderCOP dsbuilder
	);

	void
	assemble_remark_lines(
		core::pose::Pose & outpose,
		utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres,
		std::map< core::Size, core::Size > const & redundant_upstream_res,
		utility::vector1< core::Size > const & ex_geom_ids_for_upstream_res
	) const;

	core::pose::PoseCOP
	create_output_upstream_pose(
		utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres,
		std::map< core::Size, core::Size > const & redundant_upstream_res,
		utility::vector1< core::Size > const & ex_geom_ids_for_upstream_res
	);

	std::string
	signature_string(
		utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres
	) const;

	core::Size
	num_geom_cst() const;

	std::string
	scaf_name() const;

	std::string
	cstfile_name() const;

	std::string
	prefix() const;

protected:

	UpstreamHitCacherOP
	coordinate_cacher();

	utility::vector1< downstream::DownstreamBuilderCOP > const &
	dsbuilders();

	std::string
	assemble_outtag(
		utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres
	);


private:

	UpstreamHitCacherOP coordinate_cacher_;

	// whether to write the whole pose or only the matched residues
	bool write_matchres_only_;

	std::string scaf_name_, cstfile_name_, prefix_;

	core::Size num_written_;
	core::Size num_geom_cst_;

	//needed to build the downstream target for pdb writing
	utility::vector1< downstream::DownstreamBuilderCOP > dsbuilders_;


	core::pose::PoseCOP orig_upstream_pose_;
	core::pose::PoseCOP downstream_pose_from_task_;

	//map to keep track of the output names
	std::map< std::string, SizePair > signature_map_;

	//which of the cst res to write separate matches for
	utility::vector1< bool > output_dsgeom_for_geomcst_;

	//info about what the targets are for upstream only matching
	std::map< core::Size, core::Size > upstream_only_geom_cst_;

};

/// @brief an output writer that uses a grouper to group matches and
/// then writes out one pdb file per group, with the different hits
/// from the group in different MODEL sections
class CloudPDBWriter : public PDBWriter
{

public:
	typedef PDBWriter parent;
	typedef utility::vector1< std::set< upstream_hit > > UpstreamHitSets;
	typedef utility::vector1< std::set< downstream_hit > > DownstreamHitSets;

	CloudPDBWriter( MatchGrouperOP grouper );

	virtual ~CloudPDBWriter();

	virtual
	void
	prepare_for_output_writing();

	void
	end_output_writing();

	/// @brief no writing in this function, only saving the hits according
	/// to what group they belong to
	virtual
	void
	record_match( match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer );

	virtual
	void
	record_match( match_dspos1 const & m );


	/// @brief this is where the actual writing happens
	void
	write_match_groups();

	//accessors
public:

	utility::vector1< UpstreamHitSets > const & match_groups_ushits() const;
	utility::vector1< DownstreamHitSets > const &  match_groups_dshits() const;
	utility::vector1< match_dspos1 > const & representative_group_matches() const;

	utility::vector1< std::set< downstream_hit >::const_iterator > const & ds_hitset_its() const;
	utility::vector1< std::set< downstream_hit >::const_iterator > const & ds_hitset_end_its() const;
protected:

	void
	setup_hitset_iterators_for_group(
		core::Size const group );

	void
	clear_match_data();

private:

	MatchGrouperOP grouper_;

	utility::vector1< UpstreamHitSets > match_groups_ushits_;
	utility::vector1< DownstreamHitSets > match_groups_dshits_;
	utility::vector1< std::string > unique_match_names_;
	utility::vector1< match_dspos1 > representative_group_matches_;

	utility::vector1< std::set< upstream_hit >::const_iterator > us_hitset_its_;
	utility::vector1< std::set< upstream_hit >::const_iterator > us_hitset_end_its_;
	utility::vector1< std::set< downstream_hit >::const_iterator > ds_hitset_its_;
	utility::vector1< std::set< downstream_hit >::const_iterator > ds_hitset_end_its_;

};

/// @brief helper class for the MatcherMover that will put a match
/// into a supplied pose
class PoseMatchOutputWriter : public CloudPDBWriter {

public:

	typedef CloudPDBWriter parent;

	PoseMatchOutputWriter( MatchGrouperOP grouper );

	~PoseMatchOutputWriter();

	void
	insert_match_into_pose(
		core::pose::Pose & pose,
		core::Size match_group
	);

	void
	insert_match_into_pose(
		core::pose::Pose & pose
	);

	void
	end_output_writing();


};


}
}
}

#endif
