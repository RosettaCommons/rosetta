// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/znhash/SymmZnMoversAndTaskOps.hh
/// @brief  Declaration of rosetta-scripts friendly movers and task operations for designing
///         tetrahedral zinc coordination sites at symmetric protein interfaces
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Bryan Der (bder@email.unc.edu)

#ifndef INCLUDED_devel_znhash_SymmZnMoversAndTaskOps_HH
#define INCLUDED_devel_znhash_SymmZnMoversAndTaskOps_HH

// Unit headers
#include <devel/znhash/SymmZnMoversAndTaskOps.fwd.hh>

// Package headers
#include <devel/znhash/ZnHash.hh>

// Protocols headers
#include <protocols/simple_moves/ReturnSidechainMover.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
//#include <utility/fixedsizearray1.hh>

// Numeric headers
//#include <numeric/polynomial.fwd.hh>
#include <numeric/HomogeneousTransform.hh>
//#include <numeric/xyzVector.hh>
//#include <numeric/geometry/BoundingBox.hh>
//#include <numeric/geometry/hashing/SixDHasher.hh>

// Boost headers
#include <boost/unordered_map.hpp>

namespace devel {
namespace znhash {

class InitializeZNCoordinationConstraintMover : public protocols::moves::Mover
{
public:
	typedef protocols::moves::Mover parent;
	typedef protocols::moves::MoverOP MoverOP;

public:
	InitializeZNCoordinationConstraintMover();

	void set_zn_reach( core::Real reach );
	void set_orbital_dist( core::Real dist );
	void set_orbital_reach( core::Real reach );
	void set_zn_well_depth( core::Real depth );
	void set_clash_weight( core::Real weight );
	void require_3H( bool setting );
	void set_idealize_input_virtual_atoms( bool setting );
	void set_matcher_constraint_file_name( std::string const & fname );
	void set_reference_pdb( std::string const & fname );
	void set_match_pdb_listfilename( std::string const & fname );

	virtual MoverOP clone() const;
	virtual std::string get_name() const;
	virtual void apply( core::pose::Pose & p );
	ZnCoordinationScorerOP zn_score() const;

	protocols::simple_moves::ReturnSidechainMoverOP recover_sidechains() const;

private:
	core::Real znreach_; // how far does the zn score extend out to?
	core::Real orbital_dist_; // how far from the Zn center should the orbitals be located?
	core::Real orbital_reach_; // at what point are the orbitals not within striking distance?
	core::Real znwelldepth_; // how much should zn proximity count relative to orbital proximity?  <-- scale the whole score by modifying the metalhash_constraint weight.

	// how much more does the alignment of the zn atom contribute to the score
	// than the alignment of the oribital locations
	core::Real clash_weight_;
	bool require_3H_;
	bool idealize_input_virtual_atoms_;
	std::string matcher_constraint_file_name_;
	std::string reference_pdb_;
	std::string match_pdb_listfilename_;

	ZnCoordinationScorerOP zn_score_;
	protocols::simple_moves::ReturnSidechainMoverOP recover_sidechains_;

};

class ZNCoordinationConstraintReporterMover : public protocols::moves::Mover
{
public:
	typedef protocols::moves::Mover parent;
	typedef protocols::moves::MoverOP MoverOP;

public:
	ZNCoordinationConstraintReporterMover(
		InitializeZNCoordinationConstraintMoverOP init_zn
	);

	virtual MoverOP clone() const;
	virtual std::string get_name() const;
	virtual void apply( core::pose::Pose & p );
private:
	InitializeZNCoordinationConstraintMoverOP init_zn_;

};

struct matchfile_header {
	std::string remark_string_;
	std::string remark_number_;
	std::string match_string_;
	std::string template_string_;
	std::string lig_chain_;
	std::string lig_resname_;
	std::string lig_resnum_;
	std::string matchstr2_;
	std::string motifstr_;
	std::string prot_chain_;
	std::string aa_str_;
	std::string aares_;
	std::string geocst_indstr_;
	std::string geocst_indstr2_;
};

matchfile_header
read_match_header_line_from_pdb(
	std::string const & fname,
	core::Size linenum,
	std::istringstream & matchline_stream
);

struct remark_protres_coordination_data
{
	remark_protres_coordination_data() :
		name3_("XXX"),
		resindex_( 0 ),
		exgeom_index_("0"),
		chain_('X' )
	{}
	remark_protres_coordination_data( core::Size ) :
		name3_("XXX"),
		resindex_( 0 ),
		exgeom_index_("0"),
		chain_('X' )
	{}

	std::string name3_;
	core::Size  resindex_;
	std::string exgeom_index_;
	char        chain_;
};

struct remark_znres_coordination_data
{
	remark_znres_coordination_data() :
		resind_( 0 ),
		chain_( 'X' )
	{}

	remark_znres_coordination_data( core::Size ) :
		resind_( 0 ),
		chain_( 'X' )
	{}

	core::Size resind_;
	char       chain_;
};

struct symdes_znx_coordination_data
{
	utility::fixedsizearray1< remark_znres_coordination_data, 3 > zinc_data_;
	utility::fixedsizearray1< remark_protres_coordination_data, 6 > protein_data_;
};

void
add_znx_coordination_remark_lines_to_pose(
	core::pose::Pose & p,
	symdes_znx_coordination_data const & coordination_data
);

class ZNCoordinationConstraintPlacerMover : public protocols::moves::Mover
{
public:
	typedef protocols::moves::Mover parent;
	typedef protocols::moves::MoverOP MoverOP;
	typedef numeric::HomogeneousTransform< core::Real > HTReal;

public:
	ZNCoordinationConstraintPlacerMover(
		InitializeZNCoordinationConstraintMoverOP init_zn
	);

	void set_constraint_energy_cutoff( core::Real setting );
	void set_four_residue_cst_fname( std::string const & fname );

	virtual MoverOP clone() const;
	virtual std::string get_name() const;
	virtual void apply( core::pose::Pose & p );
private:
	void mutate_the_interface_to_alanine( core::pose::Pose & p );

	void insert_zn_residues_into_pose( core::pose::Pose & p );
	void add_matcher_remark_lines_for_zn_coordination( core::pose::Pose &  p );

	core::scoring::ScoreFunctionOP
	fa_sfxn_w_cstterms() const;

	void minimize_zinc_coordination( core::pose::Pose &  p );
	void quick_hires_symdock( core::pose::Pose & p );
	void restore_alanine_interface_residues_to_wtconf(
		core::pose::Pose const & pose_w_scs,
		core::pose::Pose & p
	);

	void filter_by_constraint_score( core::pose::Pose const & p );

	devel::znhash::ZnMatchData const & m1() const { return m1_; }
	devel::znhash::ZnMatchData const & m2() const { return m2_; }

private:
	InitializeZNCoordinationConstraintMoverOP init_zn_;
	devel::znhash::ZnMatchData m1_; // best match chain A
	devel::znhash::ZnMatchData m2_; // best match chain B
	core::pack::task::PackerTaskOP mutate_interface_residues_to_alanine_task_;
	core::Real constraint_energy_cutoff_;
	std::string cst_fname_;
};


class FindZnCoordinatingResidues : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FindZnCoordinatingResidues();
	FindZnCoordinatingResidues();
	void fail_on_absent_coordinators( bool setting );
	void find_coordinating_residues( core::pose::Pose const & p );
	utility::vector1< core::Size > const & resinds() const { return resinds_; }
	utility::vector1< core::id::AtomID > const & atomids() const { return atomids_; }
private:

	std::pair< core::Real, core::Size >
	closest_distance_to_desired_vrt(
		core::Vector const & zncoord,
		core::Vector const & vcoord,
		core::conformation::Residue const & nbr,
		utility::vector1< std::string > const & coord_atnames
	) const;

private:
	bool fail_on_absent_coordinators_;
	utility::vector1< core::Size > resinds_;
	utility::vector1< core::id::AtomID > atomids_;
};

class InsertZincCoordinationRemarkLines : public protocols::moves::Mover
{
public:
	InsertZincCoordinationRemarkLines();
	~InsertZincCoordinationRemarkLines();

	virtual protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	virtual void apply( core::pose::Pose & p );

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP const,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );
};


class DisableZnCoordinationResiduesTaskOp : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation    TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP  TaskOperationOP;
	typedef core::pack::task::operation::TaskOperation    parent;
public:
	DisableZnCoordinationResiduesTaskOp();
	virtual ~DisableZnCoordinationResiduesTaskOp();

	virtual TaskOperationOP clone() const;

	virtual	void apply(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask & task
	) const;

	virtual void parse_tag( TagCOP, DataMap & );

private:
	utility::vector1< Size > resids_to_fix_;
};

class ZnCoordNumHbondCalculator : public core::pose::metrics::EnergyDependentCalculator {
public:

	ZnCoordNumHbondCalculator();

	core::pose::metrics::PoseMetricCalculatorOP clone() const;
	virtual void notify_structure_change();
	virtual void notify_energy_change();

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:
	core::Size all_Hbonds_;
	//core::Size special_region_Hbonds_;
	core::id::AtomID_Map< core::Size > atom_Hbonds_;
	utility::vector1< core::Size > residue_Hbonds_;

	protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator nhbcalc_;
	FindZnCoordinatingResidues finder_;
};


/// @brief Handles sphere-sphere overlap calculations
class LoadZnCoordNumHbondCalculatorMover : public protocols::moves::Mover {

public:
	LoadZnCoordNumHbondCalculatorMover();
	~LoadZnCoordNumHbondCalculatorMover();

	virtual protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	virtual void apply( core::pose::Pose & p );

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP const,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

};


}
}

#endif
