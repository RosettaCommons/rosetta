// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg

#ifndef INCLUDED_protocols_ddg_ddGMover_hh
#define INCLUDED_protocols_ddg_ddGMover_hh

/// Where is the ddGMover.fwd.hh include?

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/scoring/constraints/ConstraintSet.fwd.hh>


#include <protocols/scoring/Interface.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

#include <core/chemical/AA.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ddg {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace scoring;

typedef std::vector<double> ddGs;

class ddGMover: public moves::Mover{
public:
	ddGMover();

	ddGMover(
		core::scoring::ScoreFunctionOP s, // for scoring the pose?
		core::scoring::ScoreFunctionOP m, // for minimizing the pose?
		utility::vector1<core::chemical::AA> res_to_mutate
	);

	virtual ~ddGMover();

	//setter function
	void neighbor_cutoff(double cutoff);
	void restrict_to_nbrs(bool truefalse);
	void score_function( core::scoring::ScoreFunctionOP s);
	void set_minimization_score_function( core::scoring::ScoreFunctionOP s);
	// Undefined, commenting out to fix PyRosetta build
	//void set_constraint_set( core::scoring::constraints::ConstraintSetOP cs);
	void num_iterations(int num);
	void dump_pdbs(bool truefalse);
	void debug_output(bool truefalse);
	void is_interface_ddg(bool truefalse);
	void wt_score_components(ObjexxFCL::FArray2D<double> wsc);
	void wt_unbound_score_components(ObjexxFCL::FArray2D<double> wusc);
	void mutant_score_components(ObjexxFCL::FArray2D<double> msc);
	void residues_to_mutate(utility::vector1<core::chemical::AA> residues);
	void set_min_cst(bool truefalse);
	void set_mean(bool truefalse);
	void set_min(bool truefalse);
	void set_num_decoys_used_in_calculations(core::Real num_lowe_used);

	//accessor function
	core::Real neighbor_cutoff();
	bool restrict_to_nbrs();
	core::scoring::ScoreFunctionOP score_function();
	core::scoring::ScoreFunctionOP minimization_score_function();
	// Undefined, commenting out to fix PyRosetta build
	//core::scoring::constraints::ConstraintSetOP constraint_set();
	int num_iterations();
	bool is_interface_ddg();
	ObjexxFCL::FArray2D<double> wt_score_components();
	ObjexxFCL::FArray2D<double> wt_unbound_score_components();
	ObjexxFCL::FArray2D<double> mutant_score_components();
	utility::vector1<core::chemical::AA> residues_to_mutate();

	utility::vector1<double> get_wt_min_score_components();
	utility::vector1<double> get_wt_averaged_score_components();
	utility::vector1<double> get_mutant_min_score_components();
	utility::vector1<double> get_mutant_averaged_score_components();
	utility::vector1<double> get_delta_energy_components();
	void get_scorefunction_header(
		core::scoring::ScoreFunctionOP sfxn,
		utility::vector1<std::string> & components
	);
	std::string mutation_label(pose::Pose const & pose) const;
	double get_wt_averaged_totals();
	double get_wt_min_totals();
	double get_mutant_averaged_totals();
	double get_mutant_min_totals();
	double ddG();
	virtual void apply(core::pose::Pose & pose);
	virtual std::string get_name() const;
	bool is_wt_calc_complete();
	bool is_mutant_calc_complete();
	bool is_properly_initialized(pose::Pose & pose);
	bool get_min_cst();
	bool use_mean();
	bool use_min();
	core::Real get_num_decoys_used_in_calculations();

private:
	core::Real num_decoys_used_in_calculations_;
	bool mean_;
	bool min_;
	bool min_cst_;
	core::Real nbr_cutoff_;
	bool restrict_to_nbrhood_;
	bool interface_ddg_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP min_cst_sfxn_; /// minimize the score with this score function
	core::scoring::ScoreFunctionOP min_cst_sfxn_no_cst_weight_; /// report the score with this score function
	core::scoring::constraints::ConstraintSetOP cst_set_;

	core::scoring::constraints::ConstraintSetOP min_cst_set_wt_;
	utility::vector1<char> min_cst_wt_types_;
	core::scoring::constraints::ConstraintSetOP min_cst_set_mut_;
	utility::vector1<char> min_cst_mut_types_;
	core::scoring::constraints::ConstraintSetOP repack_cst_set_wt_;
	utility::vector1<char> repack_wt_types_;
	core::scoring::constraints::ConstraintSetOP repack_cst_set_mut_;
	utility::vector1<char> repack_mut_types_;

	utility::vector1<core::chemical::AA> residues_to_mutate_;
	int num_iterations_;
	bool dmp_pdb_;
	bool dbg_output_;
	ObjexxFCL::FArray2D<double> wt_components_;
	ObjexxFCL::FArray2D<double> wt_unbound_components_; //ek interface ddg capabilities;
	ObjexxFCL::FArray2D<double> mutant_components_;
	ObjexxFCL::FArray2D<double> mutant_unbound_components_; //ek interface ddg capabilities;
	utility::vector1<pose::Pose> natives_;
	utility::vector1<pose::Pose> mutants_;

	utility::vector1<int>
	find_nbrs(
		core::pose::Pose & p,
		utility::vector1<int> & mutation_position,
		double radii
	);

	void setup_packertask(core::pack::task::PackerTaskOP pt, core::pose::Pose & p);
	//void refine_structure();


	bool is_complete(ObjexxFCL::FArray2D<double> to_check);

	double sum(ddGs &scores_to_sum);
	double average( utility::vector1<double> &scores_to_average);
	int store_energies(
		ObjexxFCL::FArray2D< double > &two_d_e_arrays,
		core::scoring::ScoreFunction &s,
		core::pose::Pose &p,
		int next_index,
		int size_to_expect
	);

	int average_score_components(
		ObjexxFCL::FArray2D< double > &scores_to_average,
		utility::vector1<double> &averaged_scores
	);

	//int which_iteration(ObjexxFCL::FArray2D<double> scores);//locates uninitialized values
	//and returns the iteration to continue on
	//above not used/???

	void calculate_interface_unbound_energy(core::pose::Pose & p,core::scoring::ScoreFunctionOP s,core::pack::task::PackerTaskOP pt);

	//void setup_rotamer_constraints(pose::Pose & p,ScoreFunction & s,
	//                               utility::vector1<int> & res_to_include)
	void setup_rotamer_constraints(
		core::pose::Pose & p,
		core::scoring::ScoreFunction & s,
		utility::vector1<int> & mutation_position,
		bool all_but_mutation_site,
		utility::vector1<char> & constraints_at_pos
	);

	/// @brief sets up only rotamer constraints available during repacking
	void setup_repack_constraints(
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionOP sfxn,
		bool all_but_mut_site,
		utility::vector1<char> & constraints_at_pos
	);

	void setup_constraints(
		pose::Pose & pose,
		core::scoring::ScoreFunctionOP sfxn,
		float const cst_tol,
		bool all_but_mut_site,
		utility::vector1<char> & constraints_at_pos
	);

	void minimize_with_constraints(
		pose::Pose & pose,
		core::scoring::ScoreFunctionOP s,
		core::pack::task::PackerTaskOP pt
	);

	void print_verbose_ddgs(
		utility::vector1<pose::Pose> mut,
		utility::vector1<pose::Pose> wt,
		core::scoring::ScoreFunctionOP sfxn,
		bool mean,
		bool min,
		std::string filename
	);

private:
	/// APL's attempts at refactoring this code

	/// @brief Apply the protocol to the wild type species.
	void relax_wildtype_structure(
		core::pose::Pose & pose,
		protocols::scoring::Interface & protein_interface,
		std::string const & wt_traj,
		bool output_silent
	);

	/// @brief Initialize a packer task to reflect the residues that should be repacked
	/// based on the structure of the wildtype pose.
	void
	setup_packer_task_for_mutations(
		core::pose::Pose const & pose,
		protocols::scoring::Interface const & protein_interface,
		utility::vector1<int> const & mutations,
		core::pack::task::PackerTaskOP packer_task
	);

	void
	initialize_task_level_behavior(
		core::pack::task::PackerTaskOP packer_task
	);


	bool
	is_any_pdb_empty(
		core::pose::Pose const & pose,
		std::string const & mutant_traj,
		bool const output_silent
	) const;

	utility::vector1< bool >
	neighborhood_of_mutations(
		core::pose::Pose const & pose,
		utility::vector1< int > const & mutations
	) const;

	void
	initialize_rotamer_behavior_for_residue_level_task(
		core::pack::task::ResidueLevelTask & rltask
	) const;

};

} //moves
} //protocols

#endif

