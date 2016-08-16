// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/Mather.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_Matcher_hh
#define INCLUDED_protocols_match_Matcher_hh

// Unit headers
#include <protocols/match/Matcher.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/BumpGrid.fwd.hh>
#include <protocols/match/MatchSet.fwd.hh>
#include <protocols/match/MatcherTask.fwd.hh>
#include <protocols/match/OccupiedSpaceHash.fwd.hh>

#include <protocols/match/downstream/ActiveSiteGrid.fwd.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.fwd.hh>
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.fwd.hh>

#include <protocols/match/upstream/ProteinUpstreamBuilder.fwd.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.fwd.hh>
#include <protocols/match/upstream/UpstreamBuilder.fwd.hh>

#include <protocols/match/output/MatchProcessor.fwd.hh>

// Project headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.fwd.hh>

#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/LexicographicalIterator.fwd.hh>
#include <utility/vector1_bool.hh>

// Numeric headers
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <list>
#include <map>

#include <utility/vector1.hh>


//auto headers
#ifdef WIN32
#include <core/id/AtomID.hh>
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>
#include <protocols/match/output/MatchProcessor.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>
#endif


namespace protocols {
namespace match {

/// Overview:
/// The matcher algorithm was originally concieved of within the domain of enzyme design.
/// The transition state for the desired reqction is contacted by several amino acids
/// each with a particular geometry.  The goal of the matcher is to find a set of backbone
/// positions on a given protein-backbone scaffold where those amino acids could be grafted
/// such that they would contact the ligand in the desired geometry.
///
/// Consider a case where the transition state is contacted by an asparagine, an aspartate
/// and a histadine.  The user designing an enzyme for this transition state knows the geometry
/// that describes the orientation of the transition state with respect to each of these
/// side chains; what they do not know is into what protein and at what positions they should
/// introduce these amino acids.  The user will give the matcher a description of the geometry
/// between the amind acids and the transition state. This geometry is in the form of 6 parameters:
/// 3 diherals, 2 angles, and 1 distance (more on these later).  Given the coordinates
/// of a particular side chain and the geometry describing the transition state relative to
/// the side chain, the coordinates of the transition state may be computed.  In a sense,
/// the transition state may be grown off of the end of a side chain in the desired geoemtry.
///
/// (Usually, the user will specify many different possible values for each of the 6 parameters,
/// and the matcher will then consider all combinations of those values.  E.g. the ideal
/// distance might be 2.0 A, but the user might ask the matcher to consider the values
/// 1.95 and 2.05 A  additionally.  Each assignment of values to these 6 parameters fully
/// specifies the coordinates of the transition state.)
///
/// The matcher examines each geometric constraint one at a time.  It builds rotamers
/// for one or more amino acids capable of satisfying a desired geometry (e.g. both ASP and GLU
/// if an acid group is needed) at each of several active-site positions, and for each rotamer, it
/// grows the transition state.  The matcher does a quick collision check between the atoms of the
/// transition state and the backbone of the protein, rejecting transition-state conformations
/// that collide.  If the conformation is collision-free, then the matcher measures the coordinates
/// of the transition state as a point in a 6-dimensional space.  (It is no coincidence that
/// there are 6 geometric parameters and that there are 6 dimensions in the space describing the
/// transition state's coordinates). With this 6-dimensional coordinate, the matcher can recover
/// the coordinates for the transition state -- the 6-D coordinate and the full euclidean coordinates
/// of the transition state are interconvertable.  The matcher can also bin the coordinate.  If two
/// coordinates in 6-D are close, they will be assigned to the same bin.  This is the fundamental insight of
/// the matching algorithm: the matcher will grow the transition state from different catalytic
/// residues, and when the 6-d coordinates from different catalytic residues are assigned to the
/// same bin, then the matcher has found a set of conformations of the transition state that are
/// compatible with more than one catalytic geometry.
///
/// Each collision-free placement of the transition state is called a "hit".  If there are N
/// geometric-constrains that the matcher is asked to satisfy, then a set of N hits, one per
/// constraint, that fall into the same bin are called a "match".
///
/// In the general case, the Matcher builds hits for each of several geometric constraints.  The
/// protein scaffold in the enzyme-design example generalizes to any macro-molecular polymer
/// scaffold.  The protein rotamers in the enzyme-design example generalizes to a set of conformations
/// for the  "upstream" partner. The transition state gene in the enzyme-design example generalizes to
/// a "downstream" partner, which itself may have multiple conformations.  "Upstream" and "Downstream"
/// refer to the  order in which the coordinates of the two partners are computed.  The upstream
/// coordinates are built first, the downstream coordinates second.  Changes to the coordinates of
/// the upstream partner propagate to the coordinates of the downstream partner.
/// In the enzyme-design example, the transition state is considered to be rigid; in the general case
/// the transition state may have multiple conformations.  The downstream partner could also
/// be an entire protein -- and may have it's own set of rotameric states.  E.G. one might want to
/// match a hydrogen-bond donor on the scaffold to a serine side-chain on the target (downstream) protein.
/// The downstream partner should then be able to examine many serine rotamers for each conformation of
/// the upstream rotamer.
///
/// A hit is represented in two parts: a discrete part and a continuous part.  The discrete portion consists
/// of four integers: 1. the build-point index on the scaffold, 2. the rotamer index on the upstream partner,
/// 3. the external-geometry index, and 4. the rotamer index on the downstream partner.  The continuous portion
/// consists of 6 double-precision values representing the coordinate of the downstream partner in 6D.
/// The first three values are the x,y and z coordinates of a particular atom in the downstream partner.
/// The second three values are the phi, psi, and theta values describing the coordinate frame at this atom.
/// These three "Euler angle" parameters describe three rotations: Z(psi) * X(theta) * Z(phi) * I.
/// They are described in greater detail in src/numeric/HomogeneousTransform.hh.
/// "Phi" and "psi" here have nothing to do with the protein-backbone angles.  When a hit is binned, there
/// are two sets of parameters that describe how wide the bins in each dimension should be: the Euclidean
/// bin widths are for the xyz coordinates, and the Euler bin widths are for the Euler angles.  The
/// Euclidean bin widths are in Angstroms and the Euler bin widths are in degrees.
///
/// A Matcher object should be initialized from a MatcherTask object through the
/// intialize_from_task() method. A MatcherTask will contain an EnzConstraintIO object, and the function
/// Matcher::initialize_from_file() will be invoked as the Matcher is intialied from a MatcherTask.
/// The documentation for Matcher::inialize_from_file() describes the format of extra data
/// that may be included in the enzyme-design constraint file.  This data should live within a
/// ALGORITHM_INFO:: match ... ALGORITHM::END block inside a CST::BEGIN ... CST::END block in
/// the constraint file.
///
/// find_hits() is the main worker function.  After the matcher finishes find_hits(),
/// the matches can be read by a MatchProcessor in a call to process_matches.
class Matcher : public utility::pointer::ReferenceCount {
public:
	typedef core::Real                               Real;
	typedef core::Size                               Size;
	typedef core::Vector                             Vector;
	typedef numeric::geometry::BoundingBox< Vector > BoundingBox;
	typedef std::list< Hit >                         HitList;
	typedef std::list< Hit >::iterator               HitListIterator;
	typedef std::list< Hit >::const_iterator         HitListConstIterator;

public:
	/// Construction and Destruction
	Matcher();
	virtual ~Matcher();


public:

	/// Setup
	void set_upstream_pose( core::pose::Pose const & pose );
	void set_downstream_pose(
		core::pose::Pose const & pose,
		utility::vector1< core::id::AtomID > orientation_atoms
	);

	void set_original_scaffold_build_points( utility::vector1< Size > const & resids );

	void set_original_scaffold_build_points_for_constraint(
		Size cst_id,
		utility::vector1< Size > const & resids
	);

	void set_n_geometric_constraints( Size n_constraints );

	Size n_geometric_constraints() const {
		return n_geometric_constraints_;
	}


	void add_upstream_restype_for_constraint(
		Size cst_id,
		core::chemical::ResidueTypeCOP restype
	);

	void desymmeterize_upstream_restype_for_constraint( Size cst_id );

	void set_sample_startegy_for_constraint(
		Size cst_id,
		core::chemical::ResidueTypeCOP restype,
		Size chi,
		upstream::SampleStrategyData const & strat
	);

	void
	set_fa_dun_cutoff_for_constraint(
		Size cst_id,
		core::chemical::ResidueTypeCOP restype,
		core::Real fa_dun_cutoff
	);

	void add_external_geometry_samples_for_constraint(
		Size cst_id,
		core::chemical::ResidueTypeCOP restype,
		utility::vector1< std::string >  const & upstream_launch_atoms,
		utility::vector1< core::id::AtomID > const & downstream_3atoms,
		toolbox::match_enzdes_util::ExternalGeomSampler const & exgeom,
		Size const exgeom_id,
		bool enumerate_ligand_rotamers = false,
		bool catalytic_bond = false,
		bool build_round1_hits_twice = false
	);

	void add_secondary_upstream_match_geometry_for_constraint(
		Size geom_cst_id,
		Size target_geom_cst_id,
		core::chemical::ResidueTypeCOP candidate_restype,
		core::chemical::ResidueTypeCOP target_restype,
		utility::vector1< Size > const & candidate_atids,
		utility::vector1< Size > const & target_atids,
		toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi,
		std::string SecMatchStr,
		core::pose::Pose const & upstream_pose
	);

	void add_secondary_downstream_match_geometry_for_constraint(
		Size geom_cst_id,
		core::chemical::ResidueTypeCOP candidate_restype,
		core::chemical::ResidueTypeCOP downstream_restype,
		utility::vector1< Size > const & candidate_atids,
		utility::vector1< Size > const & target_atids,
		toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi,
		std::string SecMatchStr,
		core::pose::Pose const & upstream_pose,
		bool catalytic_bond
	);

	void set_occupied_space_bounding_box( BoundingBox const & bb );
	void set_hash_euclidean_bin_width( Real width );
	void set_hash_euler_bin_width(     Real width );
	void set_hash_euclidean_bin_widths( Vector widths );
	void set_hash_euler_bin_widths(     Vector widths );

	void set_bump_tolerance( Real permitted_overlap );

	/// @brief The primary way to initialize a Matcher is through a MatcherTask.
	void
	initialize_from_task(
		MatcherTask const & task
	);

	/// @brief Intialize the geometric constraints from the EnzConstraionIO object.
	void initialize_from_file(
		toolbox::match_enzdes_util::EnzConstraintIO const & enz_data,
		MatcherTask const & task
	);

public:

	/// @brief Main worker function
	bool find_hits();

	/// @brief After find_hits completes, use this function to have the
	/// Matcher enerate the hit-combinations (matches) and send those matches
	/// to the specified match-processor.  The match processor may do what it
	/// pleases with the matches.
	void
	process_matches( output::MatchProcessor & processor ) const;

public:

	/// Data accessors
	core::pose::PoseCOP
	upstream_pose() const;

	core::pose::PoseCOP
	downstream_pose() const;

	upstream::ScaffoldBuildPointCOP
	build_point( Size index ) const;

	upstream::UpstreamBuilderCOP
	upstream_builder( Size cst_id ) const;

	//Author: Kui Chan
	//access function to pose_build_resids_
	//Reason: Use to update the SecondaryMatcherToUpstreamResidue hit.second()
	utility::vector1< Size > const &
	get_pose_build_resids() const;

	/// @brief Return const access to a representative downstream builder for a particular
	/// geometric constraint.  All downstream builders for a single geometric constraint
	/// are required to behave the same when reconstructing the coordinates of the downstream
	/// partner from a hit; therefore, a single representative is sufficient to recover
	/// hit coordinates for any hit from a particular geometric constraint.
	downstream::DownstreamBuilderCOP
	downstream_builder( Size cst_id ) const;

	std::list< downstream::DownstreamAlgorithmCOP >
	downstream_algorithms( Size cst_id ) const;

	downstream::DownstreamAlgorithmCOP
	representative_downstream_algorithm( Size cst_id ) const;

	HitList const &
	hits( Size cst_id ) const;

	OccupiedSpaceHashCOP
	occ_space_hash() const;

	utility::vector1< upstream::ScaffoldBuildPointCOP > const &
	per_constraint_build_points( Size cst_id ) const;

	/// Non-const access

	upstream::ScaffoldBuildPointOP
	build_point( Size index );

	upstream::UpstreamBuilderOP
	upstream_builder( Size cst_id );

	bool
	has_upstream_only_geomcsts() const;

	/// @brief Return non-const access to a representative downstream builder
	/// for a particular geometric constraint
	downstream::DownstreamBuilderOP
	downstream_builder( Size );

	/// @brief Return non-const access to all of the downstream builders
	/// for a particular geometric constraint
	std::list< downstream::DownstreamBuilderOP > const &
	downstream_builders( Size cst_id ) const;


	/// @brief Non-const access to the set of downstream algorithms for a
	/// particular geometric constraint -- note that the list containing
	/// these algorithms is itself const.
	std::list< downstream::DownstreamAlgorithmOP > const &
	nonconst_downstream_algorithms( Size cst_id );

	OccupiedSpaceHashOP
	occ_space_hash();

	/// @brief Return a non-constant iterator to a HitList for a particular geometric constraint.
	/// DANGER DANGER DANGER.
	/// This access is intended to allow a DownstreamAlgorithm to delete its own non-viable hits
	/// and also to allow a DownstreamAlgorithm to delete another algorithm's non-viable hits;
	/// Actual deletion requires invoking the method Matcher::erase_hit().
	/// This non-const access is not intended for any other class.
	HitListIterator
	hit_list_begin( Size geom_cst_id );

	/// @brief Return a non-constant iterator to the end position for a
	/// HitList for a particular geometric constraint. See comments for hit_list_begin()
	HitListIterator
	hit_list_end( Size geom_cst_id );

	/// @brief To be invoked by a downstream algorithm.  Downstream algorithms may prune their
	/// old, inviable hits, through this method -- they should pass themselves in as an argument
	/// -- and they may also prune this hits for other rounds.
	/// If the should prune other-round hits, then they will trigger an update to the
	/// hit_lists_with_primary_modificiations_ list, leading to an additional pass over the
	/// geometric constraints in a primary/peripheral pattern.
	void erase_hit(
		downstream::DownstreamAlgorithm const & dsalg,
		Size geom_cst_id_for_hit,
		HitListIterator const & iter
	);

private:
	bool generate_hits();
	void prepare_for_hit_generation_for_constraint( Size cst_id );
	void generate_hits_for_constraint( Size cst_id );
	void regenerate_round1_hits();
	bool finish_hit_generation_for_constraint( Size cst_id );

	bool initialize_scaffold_build_points();
	void initialize_bump_grids();
	void initialize_active_site_grid();
	void initialize_occupied_space_hash();
	void initialize_downstream_algorithms();

	downstream::DownstreamBuilderOP
	create_ds_builder(
		Size const cst_id,
		core::chemical::ResidueTypeCOP restype,
		utility::vector1< std::string >  const & upstream_launch_atoms,
		utility::vector1< core::id::AtomID > const & downstream_3atoms,
		bool enumerate_ligand_rotamers,
		bool catalytic_bond
	);

	/// @brief Selects a subset of all possible hit combinations (e.g. by clustering hits)
	/// for a particular bin.  Useful if there are too many matches found.
	void
	select_hit_representatives(
		utility::vector1< utility::vector1< Hit const * > > const & hit_vectors,
		utility::vector1< Size > & n_hits_per_geomcst,
		utility::vector1< utility::vector1< Size > > & reps
	) const;

	bool
	check_non_upstream_only_hit_incompatibility(
		match_dspos1 const & m1,
		utility::LexicographicalIterator & lex,
		output::MatchProcessor const & processor
	) const;

	bool
	check_downstream_hit_incompatibility(
		match const & m,
		utility::LexicographicalIterator & lex
	) const;


	bool
	test_upstream_only_hit_incompatibility(
		match const & m,
		utility::vector1< HitPtrListCOP > const & upstream_only_hits,
		utility::vector1< std::list< Hit const * >::const_iterator > & upstream_only_hit_iterators,
		Size & last_upstream_only_geomcst_advanced,
		output::MatchProcessor const & processor
	) const;

	bool
	test_upstream_only_hit_incompatibility(
		match_dspos1 const & m1,
		utility::vector1< HitPtrListCOP > const & upstream_only_hits,
		utility::vector1< std::list< Hit const * >::const_iterator > & upstream_only_hit_iterators,
		Size & last_upstream_only_geomcst_advanced,
		output::MatchProcessor const & processor
	) const;

	bool
	increment_upstream_only_hit_combination(
		utility::vector1< HitPtrListCOP > const & upstream_only_hits,
		Size starting_point,
		utility::vector1< std::list< Hit const * >::const_iterator > & upstream_only_hit_iterators,
		Size & last_upstream_only_geomcst_advanced
	) const;

	void process_matches_main_loop_enumerating_all_hit_combos( output::MatchProcessor & processor ) const;

	utility::vector1< std::list< Hit const * > >
	refine_grid_and_subsample_for_hit_subsets(
		Vector & good_euclidean_bin_widths,
		Vector & good_euler_bin_widths,
		utility::vector1< std::list< Hit const * > > const & neighbor_hits
	) const;

	utility::vector1< std::list< Hit const * > >
	subsample_hits(
		Vector const & euclidean_bin_widths,
		Vector const & euler_bin_widths,
		utility::vector1< std::list< Hit const * > > const & neighbor_hits
	) const;


	Size
	predict_n_matches_for_hit_subsets(
		Vector const & euclidean_bin_widths,
		Vector const & euler_bin_widths,
		utility::vector1< std::list< Hit const * > > const & neighbor_hits,
		Size accuracy_threshold
	) const;

	MatcherOutputStats
	process_matches_all_hit_combos_for_hit_subsets(
		output::MatchProcessor & processor,
		HitHasher & hit_hasher,
		utility::vector1< std::list< Hit const * > > const & neighbor_hits
	) const;


	void
	process_matches_where_one_geomcst_defines_downstream_location( output::MatchProcessor & processor ) const;

	/// @brief Note that a change has occurred for a particular hit list
	void
	note_primary_change_to_geom_csts_hitlist( Size geom_cst_id );

private:
	/// uncopyable -- unimplemented
	Matcher( Matcher const & );
	Matcher const & operator = ( Matcher const & rhs );

private:

	core::pose::PoseOP upstream_pose_;

	core::pose::PoseOP downstream_pose_;
	utility::vector1< core::id::AtomID > downstream_orientation_atoms_;

	bool same_build_resids_for_all_csts_;
	utility::vector1< Size > pose_build_resids_;
	utility::vector1< utility::vector1< Size > > per_cst_build_resids_;

	utility::vector1< upstream::ScaffoldBuildPointOP > all_build_points_;
	utility::vector1< utility::vector1< upstream::ScaffoldBuildPointCOP > > per_constraint_build_points_;


	Size n_geometric_constraints_;

	utility::vector1< HitList > hits_;

	utility::vector1< upstream::UpstreamBuilderOP >                    upstream_builders_;
	utility::vector1< std::map< std::string, Size > >                  build_set_id_for_restype_;
	utility::vector1< downstream::DownstreamAlgorithmOP >              representative_downstream_algorithm_;
	utility::vector1< std::list< downstream::DownstreamAlgorithmOP > > downstream_algorithms_;

	/// Does the downstream algorithm want the 6D coordinate stored in hit.second() to
	/// be hashed?
	utility::vector1< bool > geomcst_is_upstream_only_;

	std::list< Size > geom_csts_with_primary_hitlist_modificiations_;
	utility::vector1< bool > geom_cst_has_primary_modification_;

	utility::vector1< std::list< downstream::DownstreamBuilderOP > > downstream_builders_;
	std::list< downstream::DownstreamBuilderOP >   all_downstream_builders_;
	std::list< downstream::DownstreamAlgorithmOP > all_downstream_algorithms_;

	BumpGridOP bb_grid_;
	utility::vector1< BumpGridOP > original_scaffold_residue_bump_grids_;

	BoundingBox occ_space_bounding_box_;
	Vector euclidean_bin_widths_;
	Vector euler_bin_widths_;
	OccupiedSpaceHashOP occ_space_hash_;

	downstream::ActiveSiteGridOP active_site_grid_;

	bool read_gridlig_file_;
	std::string gridlig_fname_;
	std::list< std::pair< Size, Real > > upstream_resids_and_radii_defining_active_site_;
	std::list< core::id::AtomID > downstream_atoms_required_inside_active_site_;
	utility::vector1< core::id::AtomID > relevant_downstream_atoms_;

	bool use_input_sc_;

	bool dynamic_grid_refinement_;
	bool output_matches_as_singular_downstream_positioning_; // use match_dspos1 output pathway?
	mutable bool check_potential_dsbuilder_incompatibility_;
	utility::vector1< bool > output_match_dspos1_for_geomcst_;

	bool build_round1_hits_twice_;
};

class MatcherOutputStats
{
public:
	MatcherOutputStats() :
		num_potential_matches( 0 ),
		num_sent_to_proc( 0 ),
		num_non_up_only_incompatible( 0 ),
		num_up_only_incompatible( 0 ),
		num_considered_muliple_origins( 0 ),
		all_lex_states( 0 ),
		num_ds_hit_incompatible( 0 ),
		num_empty_uplist( 0 )
	{}

	MatcherOutputStats const & operator += ( MatcherOutputStats const & rhs )
	{
		num_potential_matches += rhs.num_potential_matches;
		num_sent_to_proc += rhs.num_sent_to_proc;
		num_non_up_only_incompatible += rhs.num_non_up_only_incompatible;
		num_up_only_incompatible += rhs.num_up_only_incompatible;
		num_considered_muliple_origins += rhs.num_considered_muliple_origins;
		all_lex_states += rhs.all_lex_states;
		num_ds_hit_incompatible += rhs.num_ds_hit_incompatible;
		num_empty_uplist += rhs.num_empty_uplist;

		return *this;
	}

	core::Size num_potential_matches;
	core::Size num_sent_to_proc;
	core::Size num_non_up_only_incompatible;
	core::Size num_up_only_incompatible;
	core::Size num_considered_muliple_origins;
	core::Size all_lex_states;
	core::Size num_ds_hit_incompatible;
	core::Size num_empty_uplist;

};

}
}

#endif
