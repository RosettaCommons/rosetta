// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/protein_interface_design/design_utils.hh
/// @brief definition of various classes for interface design.
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_design_utils_hh
#define INCLUDED_protocols_protein_interface_design_design_utils_hh

// Project Headers
#include <utility/exit.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <list>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {

// @brief returns the weighted total energy of a given residue. assumes that the pose has been previously scored.
core::Real sum_total_residue_energy( core::pose::Pose const & pose, core::Size const resid );

/// @details Class ReportSequenceDifferences takes in two poses and provides information on the sequence
/// changes between them, including the residue energies associated with those changes.
class ReportSequenceDifferences
{
public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
public:
	ReportSequenceDifferences( core::scoring::ScoreFunctionOP scorefxn )
	{
		scorefxn_ = scorefxn;
	}
	ReportSequenceDifferences( ReportSequenceDifferences const & init ) { // copy constructor
		res_energy1_ = init.res_energy1_;
		res_energy2_ = init.res_energy2_;
		res_name1_ = init.res_name1_;
		res_name2_ = init.res_name2_;
		scorefxn_ = init.scorefxn_->clone();
	}
	void calculate( Pose const & pose1, Pose const & pose2 );
	std::map< Size, Real> const * get_res_energy( Size const num ) const {
		runtime_assert(num==1 || num==2);
		return( num==1 ? &res_energy1_ : &res_energy2_ );
	}
	void report( std::ostream & out ) const;
	std::map< Size, std::string > const & res_name1() const { return res_name1_; }
	std::map< Size, std::string > const & res_name2() const { return res_name2_; }
	virtual ~ReportSequenceDifferences() {};
private:
	std::map< Size, Real > res_energy1_;
	std::map< Size, Real > res_energy2_;
	std::map< Size, std::string > res_name1_;
	std::map< Size, std::string > res_name2_;
	core::scoring::ScoreFunctionOP scorefxn_;
};

/// @details class Revert takes in 'wt' and 'designed' poses and attempts to revert all substitutions in the
/// design to their wt identities. Each substitution is tried separately in the context of the designed protein
/// and reversions that don't adversely affect ddg are made. If the energy of the residue in the design is
/// higher than 0, but the reversion did not succeed, Revert will attempt an Ala substitution.
class Revert
{
public:
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
public:
	Revert( ScoreFunctionCOP scorefxn, core::Real const ddg_tolerance, core::Size ddg_cycles = 5 ){
		scorefxn_ = scorefxn->clone();
		ddg_tolerance_ = ddg_tolerance;
		ddg_cycles_ = ddg_cycles;
	}
	Revert( Revert const & init ) { // copy constructor
		scorefxn_ = init.scorefxn_->clone();
		ddg_tolerance_ = init.ddg_tolerance_;
		ddg_cycles_ = init.ddg_cycles_;
	}
	void apply( core::pose::Pose & pose_wt, core::pose::Pose & pose_des ) const;
	virtual ~Revert(){};
private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real ddg_tolerance_;
	core::Size ddg_cycles_;
};

/// @details class FavorNativeResidue changes a pose object so that its residue identities at the
/// initialization of FavorNativeResidue are kept in memory. If the res_type_constraint score term is set
/// to a value other than 0, an energy bonus will be assigned if the residue doesn't change. This is useful
/// e.g., in design based on a native scaffold where we want a barrier to mutation.
class FavorNativeResidue
{
public:
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
public:
	FavorNativeResidue( Pose & pose, Real const native_residue_bonus );
	FavorNativeResidue( Pose & pose, utility::vector1< Real > const native_residue_bonus );
	virtual ~FavorNativeResidue(){};
private:
	void add_residue_constraints( Pose & pose ) const;
	utility::vector1< Real > native_residue_bonus_;
};

/// @details class FavorNonNativeResidue changes a pose object so that its residue identities at the
/// initialization of FavorNonNativeResidue are kept in memory. If the res_type_constraint score term is set
/// to a value other than 0, an energy bonus will be assigned if the residue changes. This is useful
/// if we want to encourage mutation.  A negative score would discourage mutation.
class FavorNonNativeResidue
{
public:
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
public:
	FavorNonNativeResidue( Pose & pose, Real const native_residue_bonus );
	FavorNonNativeResidue( Pose & pose, utility::vector1< Real > const native_residue_bonus );
	virtual ~FavorNonNativeResidue(){};
private:
	void add_residue_constraints( Pose & pose ) const;
	utility::vector1< Real > non_native_residue_bonus_;
};


} // protein_interface_design
} // devel

/// @brief utility function for minimizing sidechain in rigid-body dof, the interface sc, and bb in the entire protein.
// The fold_tree for minimization will be set from the centre of target_residues to the closest residue on the partner.
// The packertask is used to decide which residues to minimize (those that are not set to prevent_repacking)
void MinimizeInterface( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP scorefxn, utility::vector1< bool > const min_bb, utility::vector1< bool > const min_sc, utility::vector1< bool > const min_rb, bool const optimize_foldtree, utility::vector1< core::Size > const target_residues, bool const simultaneous_minimization = false );

void SymMinimizeInterface( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP scorefxn, utility::vector1< bool > const min_bb, utility::vector1< bool > const min_sc, utility::vector1< bool > const min_rb, /*bool const optimize_foldtree, utility::vector1< core::Size > const target_residues*/ bool const simultaneous_minimization = false );

/// @brief utility function for finding hbonding partners among a list of potential binder residues to a specific target
// residue
std::list< core::Size >
hbonded( core::pose::Pose const & pose, core::Size const target_residue, std::set< core::Size > const & binders,
		 bool const bb, bool const sc, core::Real const energy_thres, bool const bb_bb = false );

/// @brief utility function for finding hbonding partners among a list of potential binder residues to a specific target
// residue and atom
std::list< core::Size >
hbonded_atom( core::pose::Pose const & pose, core::Size const target_residue, std::string target_atom, std::set< core::Size > const & binders,
		 bool const bb, bool const sc, core::Real const energy_thres, bool const bb_bb = false );

#endif /*INCLUDED_DESIGN_UTILS_H_*/


