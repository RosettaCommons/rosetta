// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/LigandAligner.hh
///
/// @brief  Align a ligand to a reference ligand or an inferred reference ligand
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_LigandAligner_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_LigandAligner_hh

#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/func/SumFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/func/ConstantFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

/// @brief
/// the atom properties important in constraint generation
class
	AtomProperties {
private:
	bool is_donor;
	bool is_acceptor;
	bool is_H;
	bool is_polarH;
	bool is_halogen;
	core::Real ambiguity_;
	core::Real score_;
	std::string tag_;

public:
	AtomProperties() : is_donor(false),is_acceptor(false),is_H(false),is_polarH(false),is_halogen(false), ambiguity_( 0.0 ), score_( 0.0 ), tag_( "" ) {}
	AtomProperties( bool donor, bool acceptor, bool H, bool polarH, bool halogen, core::Real ambiguity, core::Real score, std::string tag) :
		is_donor(donor),is_acceptor(acceptor),is_H(H),is_polarH(polarH),is_halogen(halogen),ambiguity_(ambiguity),score_(score), tag_( tag ) {}

	// property resolution: does a pair of atoms match?  What is the scale?
	// halogen_specific treats halogens separately, only matching to other halogens
	core::Real
	match( AtomProperties const &other, core::Real polar_scale, bool halogen_specific ) const;
	bool donor() const { return is_donor; }
	bool polarH() const { return is_polarH; }
	bool acceptor() const { return is_acceptor; }
	core::Real score() const { return score_; }
	void score( core::Real value ) { score_ = value; }

	core::Real ambiguity() const { return ambiguity_; } // uncertainty of position

	std::string tag() const { return tag_; }
	std::string show() const;
	bool used_for_phore() const { return (is_polarH || is_acceptor); }
};

/// @brief
/// Pharmacophore (aka binding-motif) definition
/// @details
/// Stores collection of atom properties and their clustering info

class
	Pharmacophore {

public:
	Pharmacophore();
	Pharmacophore( Pharmacophore const &p );

	Pharmacophore( utility::vector1< core::Size > atmno,
		utility::vector1< AtomProperties > const & p,
		utility::vector1< utility::vector1< core::Real > > const & Dmtrx,
		utility::vector1< numeric::xyzVector< core::Real > > const & coords );

	core::Size size() const { return atms_.size(); }
	utility::vector1< core::Size > atms() const { return atms_; } // copy
	utility::vector1< AtomProperties > props() const { return props_; } // copy
	utility::vector1< utility::vector1< core::Real > > dist() const { return dist_; } // copy
	core::Real dist( core::Size i, core::Size j) const { return dist_[i][j]; }

	core::Size atm( core::Size i ) const { return atms_[i]; }
	bool is_polarH( core::Size i ) const { return props_[i].polarH(); }
	bool is_acceptor( core::Size i ) const { return props_[i].acceptor(); }
	bool is_donor( core::Size i ) const { return props_[i].donor(); }

	bool has( core::Size i ) const  { return atms_.contains( i ); }
	std::string index() const { return index_; }

	numeric::xyzVector< core::Real > const & com() const { return com_; }
	void update_com( utility::vector1< numeric::xyzVector< core::Real > > const & coords );

	std::string show() const;

	/// @brief whether j is close to all own members (Dcol == Dmtrx[j])
	bool
	is_close( utility::vector1< core::Real > const &Dcol,
		core::Real dcut ) const;

	/// @brief find matches with proper types (e.g. acceptor - donor)
	core::Real
	find_type_matches( Pharmacophore const &other,
		bool const just_count = false ) const;

	/*
	// unused
	bool
	find_matching_permutation( Pharmacophore const &other,
	utility::vector1< core::Size > &map_index,
	core::Real &best_score ) const;
	*/

	/// @brief find matches with the partner pharmacophore provided
	core::Real
	match( Pharmacophore const &other,
		utility::vector1< core::Size > const &map_index,
		core::Real const unmatch_score ) const;

	void
	set_map_index( utility::vector1< core::Size > map_index ){ map_index_ = map_index; }

	utility::vector1< core::Size > const & map_index() const { return map_index_; }


private:
	utility::vector1< core::Size > atms_;
	utility::vector1< AtomProperties > props_;
	utility::vector1< utility::vector1< core::Real > > dist_;
	utility::vector1< core::Size > map_index_; // how ligatm maps to recatm
	numeric::xyzVector< core::Real > com_;
	std::string index_; // index within whole context
};

/// @brief
/// Constraint information about AtomProperties and/or Pharmacophores
/// @details
/// includes constrained xyzs & important properties
/// we do not use a residue type (which already has this info)
///   because in "pharmicophore mode" we generate virtual sites
///   from the protein side.
/// we want targets generated in both ways to run through the same machinery
/// therefore, we instead we have two initialization routines:
///   - from the ligand side (from restype)
///   - from the receptor side (building programmatically)
class
	ConstraintInfo {
public:
	ConstraintInfo( )
	:  limit_(4.0), weight_(3.0), polar_scale_(2.0),
		phore_dcut_(4.0), n_match_per_ligand_motif_(3), debug_(false)
	{ }

	ConstraintInfo( core::conformation::Residue const &lig,
		bool use_pharmacophore,
		bool report_phore_info=false ){
		if ( use_pharmacophore ) {
			init_from_ligand_pharmacophore(lig,report_phore_info);
		} else {
			init_from_ligand(lig);
		}
	}

	ConstraintInfo(
		core::pose::Pose const &pose,
		GridScorerOP gridscore,
		bool const simple=false
	) {
		init_from_receptor(pose, gridscore, simple );
	}

	/// @brief use a (non-identical) target ligand (for reference ligand matching)
	void
	init_from_ligand( core::conformation::Residue const &mylig);

	/// @brief defines pharmacophore partners in ligand (for pharmacophore matching)
	void
	init_from_ligand_pharmacophore( core::conformation::Residue const &mylig,
		bool report_phore_info );

	/// @details build up binding locations from a receptor
	/// gridscorer is used to define the bounding box
	void
	init_from_receptor(
		core::pose::Pose const &pose,
		GridScorerOP gridscore,
		bool const simple
	);

	core::Size natoms() const { return coords_.size(); }
	numeric::xyzVector< core::Real > const & coord(core::Size i) const { return coords_[i]; }
	AtomProperties const & properties(core::Size i) const{ return properties_[i]; }
	core::id::AtomID const & atomid(core::Size i) const { return atmids_[i]; }

	core::Real limit() { return limit_; }
	core::Real weight() { return weight_; }
	core::Real polar_scale() { return polar_scale_; }

	bool pharmacophores_defined() const { return phores_.size(); }
	utility::vector1< Pharmacophore > const& phores() const { return phores_; }

	core::Size n_phore_match() const { return phore_match_.size(); }
	//core::Size receptor_vsite_min() const { return RECEPTOR_VSITE_MIN; }
	//core::Size receptor_vsite_max() const { return RECEPTOR_VSITE_MAX; }

	// below functioned called from LigandAligner instance...

	/// @brief selects a phore_match and updates current_phore_match_
	void
	select_phore_match( core::Size const run_index );

	/// @brief try aligning to the phore defined by current_phore_match_ (defined by select_phore_match)
	void
	align_to_current_phore_match( core::pose::Pose &pose,
		ConstraintInfo const &recinfo,
		core::Size const ligid,
		utility::vector1< std::pair< core::Size, core::Size > > &marked_pairs,
		utility::vector1< core::Size > &SrcPriorIDs,
		utility::vector1< core::Size > &TgtPriorIDs
	);

	/// @brief matches ligand phores (==own) to provided receptor phores and stores at phore_match_
	void
	map_phores( utility::vector1< Pharmacophore > const &receptor_phores,
		bool const update,
		core::Size const nmax = 100 );


private:
	void set_default();

	/// @brief findes all virtual sites (for h-bonding candidates) in receptor
	void
	define_active_virtual_sites( core::pose::Pose const &pose,
		GridScorerOP gridscore,
		utility::vector1< std::pair< core::Real, core::Size > > &Vdonor_sort,
		utility::vector1< std::pair< core::Real, core::Size > > &Vacceptor_sort,
		core::Size const nneigh_cut_for_exposed_sc,
		core::Size const nneigh_cut_for_exposed_bb,
		bool const report=false
	);

	void
	append_metal_vsites( core::pose::Pose const &pose,
		core::Size const &ires,
		utility::vector1< std::pair< core::Real, core::Size > > &Vacceptor_sort,
		bool report=false
	);

	void
	build_Dmtrx();

	/// @brief defines all ligand phores
	void
	define_all_ligand_phores( core::Size const minsize,
		bool report_phore_info );

	bool
	phore_overlaps_with_existing( Pharmacophore const &phore_i,
		core::Size &jphore ) const;

	/*
	bool
	cluster_overlaps_any( utility::vector1< utility::vector1< core::Size > > const &clusters,
	utility::vector1< core::Size > const &cluster_i,
	core::Size const n,
	core::Real const fraction ) const;
	*/

	/// @brief defines receptor phores
	void
	define_receptor_phores( utility::vector1< std::pair< core::Real, core::Size > > &Vdonor_sort,
		utility::vector1< std::pair< core::Real, core::Size > > &Vacceptor_sort,
		core::Size const nmin = 10,
		bool report_phore_info = false
	);

	void
	update_ligand_coord( core::conformation::Residue const & ligand );

private:
	utility::vector1< core::id::AtomID > atmids_;
	utility::vector1< numeric::xyzVector< core::Real > > coords_;
	utility::vector1< AtomProperties > properties_;

	// protocol parameters
	core::Real limit_, weight_, polar_scale_;
	bool is_ligand_;

	// pharmacophore docking params
	core::Real phore_dcut_;
	core::Size n_match_per_ligand_motif_;
	core::Size NNEIGH_CUT_FOR_EXPOSED_SC; // nCB around to define as exposed, cutoff for sc site
	core::Size NNEIGH_CUT_FOR_EXPOSED_BB; // nCB around to define as exposed, cutoff for bb site
	core::Size RECEPTOR_PHORE_MIN; // lower-limit for receptor num. phores
	core::Size RECEPTOR_PHORE_MAX; // upper-limit for receptor num. phores

	utility::vector1< utility::vector1< core::Real > > Dmtrx_;
	utility::vector1< Pharmacophore > phores_;
	utility::vector1< std::pair< core::Size, Pharmacophore > > phore_match_;
	std::pair< core::Size, Pharmacophore > current_phore_match_;

	bool debug_;
};

/// @brief
/// Aligns ligand using defined constraint information
/// @details
/// Performs ligand alignments using constraints derived from
/// - reference ligand pose coords or
/// - receptor pharmacophore info
/// Constraint info should be properly provided in order to run this instance

class LigandAligner {
public:
	LigandAligner();
	LigandAligner( bool use_pharmacophore,
		utility::vector1< core::Size > const & movable_scs,
		bool faster );

	/// @brief main apply function
	void apply( LigandConformer & lig );

	/// @brief setup pharmacophore info from receptor
	void set_pharmacophore_reference( core::pose::Pose const &pose );

	/// @brief gets estimate of n-pharmacophore-search considering problem complexity
	core::Size estimate_nstruct_sample( core::conformation::Residue const & ligand,
		core::Size const ntotal );

	// Setters&getters
	// set target coords and scorefunction
	void set_target( ConstraintInfo const & tgt_in ) { target_ = tgt_in; }
	void set_sf( GridScorerOP sf_in ) { sf_ = sf_in; }

	// set parameters
	void set_trans_step( core::Size trans_step_in ) { trans_step_ = trans_step_in; }
	void set_rot_step( core::Size rot_step_in ) { rot_step_ = rot_step_in; }
	void set_chi_step( core::Size chi_step_in ) { chi_step_ = chi_step_in; }

	void set_use_pharmacophore( bool setting ) { use_pharmacophore_ = setting; }
	bool use_pharmacophore() const { return use_pharmacophore_; }

	void refine_input( bool setting ){ refine_input_ = setting; }

private:
	/// @brief set constraints to target
	void set_constraints( core::pose::Pose & pose, core::Size ligid,
		utility::vector1< std::pair< core::Size, core::Size > > &marked_pairs,
		core::Real const w_prior = 1.0, // default no upper limit
		utility::vector1< core::Size > const &SrcPriorIDs = utility::vector1< core::Size >(),
		utility::vector1< core::Size > const &TgtPriorIDs = utility::vector1< core::Size >()
	);

	/// @brief set stronger constraints on specific constraint set pairs
	void
	set_hard_constraint_on_marked(
		core::pose::Pose & pose, core::Size const ligid,
		utility::vector1< std::pair< core::Size, core::Size > > const &marked_pairs
	) const;

	/// @brief randomize ligand about a new center 'T'
	void randomize_lig( core::pose::Pose & pose, core::Size ligid, numeric::xyzVector<core::Real> const &T );

	/// @brief perturb ligand
	void perturb_lig( core::pose::Pose & pose, core::Size ligid );

private:
	ConstraintInfo target_; // target pose to which we are aligning
	GridScorerOP sf_;      // target scorefunction

	// perturbation parameters
	utility::vector1< core::Size > movable_scs_;
	bool refine_input_;
	core::Real trans_step_, rot_step_, chi_step_;
	utility::vector1< utility::vector1< core::Real > > weighted_score_ij_;

	// history of howmany structs generated through
	core::Size istruct_;

	// pharmacophore docking
	bool use_pharmacophore_;

	// estimation of nsamples to try
	//core::Size nstruct_sample_estimation_;

	// faster version (support for VSX mode in GAdock)
	bool faster_;

};

typedef utility::pointer::shared_ptr< LigandAligner > LigandAlignerOP;
typedef utility::pointer::shared_ptr< LigandAligner const > LigandAlignerCOP;

}
}
}

#endif

