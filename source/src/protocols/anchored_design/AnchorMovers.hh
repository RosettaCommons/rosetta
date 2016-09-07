// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/anchored_design/AnchorMovers.hh
/// @brief protocol-level outermost movers for anchored design; wraps a set of other moves
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_anchored_design_AnchorMovers_hh
#define INCLUDED_protocols_anchored_design_AnchorMovers_hh

// Unit Headers
#include <protocols/anchored_design/AnchorMovers.fwd.hh>
#include <protocols/anchored_design/AnchorMoversData.fwd.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


// Utility Headers

namespace protocols {
namespace anchored_design {

class AnchoredDesignMover : public protocols::moves::Mover {

public:
	/// @brief constructor with arguments
	AnchoredDesignMover( protocols::anchored_design::AnchorMoversDataOP interface_in);

	/// @brief ctor with no arguments
	AnchoredDesignMover();

	/// @brief copy ctor
	AnchoredDesignMover( AnchoredDesignMover const & rhs );

	/// @brief assignment operator
	AnchoredDesignMover & operator=( AnchoredDesignMover const & rhs );

	~AnchoredDesignMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	bool reinitialize_for_new_input() const override;
	bool reinitialize_for_each_job() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

private:
	/// @initializes internals; must wait for a pose to initialize
	void init_on_new_input( core::pose::Pose const & pose );

public:
	/// @brief creates the anchored design fold tree and applies it to the pose
	void set_fold_tree_and_cutpoints( core::pose::Pose & pose );

	/// @brief runs varous filtering checks on the finished pose; sets MoverStatus for failure as needed
	void filter( core::pose::Pose & pose );

	/// @brief implements the "extended" field of the loop file specification - sets extended phi/psi as needed
	void forget_initial_loops( core::pose::Pose & pose );

	/// @brief calculate RMSD if desired; protected internally
	void calculate_rmsd( core::pose::Pose const & pose, core::pose::PoseCOP start_pose );

	/// @brief randomize the input loop sequence.  Useful if you have reason to believe the starting loop sequence is biasing to a particular unwanted structure in centroid mode.  Acts only on designable positions.
	void randomize_input_sequence( core::pose::Pose & pose ) const;

	/// @brief This function repacks the interface with use_input_sc forcibly off for benchmarking purposes.
	void delete_interface_native_sidechains( core::pose::Pose & pose ) const;

	/// @brief handles perturbing the initial anchor placement
	void perturb_anchor( core::pose::Pose & pose ) const;

	//option system replacement getters and setters
	/// @brief run RMSD calculations
	bool get_rmsd() const;
	/// @brief run only RMSD calculations against this input, don't do actual AnchoredDesign
	std::string const & get_RMSD_only_this() const;
	/// @brief delete the input sidechains (independently from use_input_sc in the packer) - used to prevent leakage of sidechains in benchmarking mode
	bool get_delete_interface_native_sidechains() const;
	/// @brief show_extended demonstrates that the code really forgets the input structure
	bool get_show_extended() const;
	/// @brief randomize_input_sequence to complement loop extension in forgetting the input
	bool get_randomize_input_sequence() const;
	/// @brief pick a different cutpoint than the input; useful when you want to sample cutpoints
	bool get_vary_cutpoints() const;
	/// @brief skip the perturbation step - useful when you already have a good structure
	bool get_refine_only() const;
	/// @brief filter based on total complex score
	core::Real get_filter_score() const;
	/// @brief filter based on complex SASA
	core::Real get_filter_SASA() const;
	/// @brief filter based on omega angles in the loops - filter out cis omegas
	bool get_filter_omega() const;
	/// @brief whether to automatically initialize from the options system; defaults to true
	bool get_autoinitialize() const;

	/// @brief run RMSD calculations
	void set_rmsd(bool const rmsd);
	/// @brief run only RMSD calculations against this input, don't do actual AnchoredDesign
	void set_RMSD_only_this(std::string const & RMSD_only_this);
	/// @brief delete the input sidechains (independently from use_input_sc in the packer) - used to prevent leakage of sidechains in benchmarking mode
	void set_delete_interface_native_sidechains(bool const delete_interface_native_sidechains);
	/// @brief show_extended demonstrates that the code really forsets the input structure
	void set_show_extended(bool const show_extended);
	/// @brief randomize_input_sequence to complement loop extension in forgetting the input
	void set_randomize_input_sequence(bool const randomize_input_sequence);
	/// @brief pick a different cutpoint than the input; useful when you want to sample cutpoints
	void set_vary_cutpoints(bool const vary_cutpoints);
	/// @brief skip the perturbation step - useful when you already have a good structure
	void set_refine_only(bool const refine_only);
	/// @brief filter based on total complex score
	void set_filter_score(core::Real const filter_score);
	/// @brief filter based on complex SASA
	void set_filter_SASA(core::Real const filter_SASA);
	/// @brief filter based on omega angles in the loops - filter out cis omegas
	void set_filter_omega(bool const filter_omega);
	/// @brief whether to automatically initialize from the options system; defaults to true
	void set_autoinitialize(bool const autoinitialize);

	/// @brief read in options from the options system
	void read_options();


private:

	protocols::anchored_design::AnchorMoversDataOP interface_;

	/// @details used for RMSD comparisons with RMSD_only_this mode
	core::pose::PoseCOP RMSD_only_this_pose_;

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

	//option system replacement
	//benchmarking mode options
	/// @brief run RMSD calculations
	bool rmsd_;
	/// @brief run only RMSD calculations against this input, don't do actual AnchoredDesign
	std::string RMSD_only_this_;
	/// @brief delete the input sidechains (independently from use_input_sc in the packer) - used to prevent leakage of sidechains in benchmarking mode
	bool delete_interface_native_sidechains_;
	/// @brief show_extended demonstrates that the code really forgets the input structure
	bool show_extended_;
	/// @brief randomize_input_sequence to complement loop extension in forgetting the input
	bool randomize_input_sequence_;
	//regular mode options
	/// @brief pick a different cutpoint than the input; useful when you want to sample cutpoints
	bool vary_cutpoints_;
	/// @brief skip the perturbation step - useful when you already have a good structure
	bool refine_only_;
	//filtering options
	/// @brief filter based on total complex score
	core::Real filter_score_;
	/// @brief filter based on total complex score
	bool use_filter_score_;
	/// @brief filter based on complex SASA
	core::Real filter_SASA_;
	/// @brief filter based on complex SASA
	bool use_filter_SASA_;
	/// @brief filter based on omega angles in the loops - filter out cis omegas
	bool use_filter_omega_;
	/// @brief whether to automatically initialize from the options system; defaults to true
	bool autoinitialize_;

	/// @brief used to determine the validity of the rest of the internals, especially the interface_ object
	bool init_for_input_yet_;

}; //class AnchoredDesignMover

///////////////////////////////////////////////////////////////////////////////////////////////////////////
class AnchoredPerturbMover : public protocols::moves::Mover {

public:
	//@brief constructor with arguments
	AnchoredPerturbMover( protocols::anchored_design::AnchorMoversDataOP interface_in );

	~AnchoredPerturbMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	//option system replacement
	/// @brief debugging mode activates a bunch of extra output
	bool get_debug() const;
	/// @brief do not perform CCD style closure (use KIC only)
	bool get_perturb_CCD_off() const;
	/// @brief do not perform KIC style closure (use CCD only)
	bool get_perturb_KIC_off() const;
	/// @brief use nonpivot torsion sampling for KIC?
	bool get_nonpivot_torsion_sampling() const;
	/// @brief MC temperature
	core::Real get_perturb_temp() const;
	/// @brief number of MC cycles
	core::Size get_perturb_cycles() const;
	/// @brief do not use fragments?
	bool get_no_frags() const;
	/// @brief what minimizer type to use?
	std::string const & get_min_type() const;
	/// @brief show perturb result structure?
	bool get_perturb_show() const;

	/// @brief debugging mode activates a bunch of extra output
	void set_debug(bool const debug);
	/// @brief do not perform CCD style closure (use KIC only)
	void set_perturb_CCD_off(bool const perturb_CCD_off);
	/// @brief do not perform KIC style closure (use CCD only)
	void set_perturb_KIC_off(bool const perturb_KIC_off);
	/// @brief use nonpivot torsion sampling for KIC?
	void set_nonpivot_torsion_sampling(bool const nonpivot_torsion_sampling);
	/// @brief MC temperature
	void set_perturb_temp(core::Real const perturb_temp);
	/// @brief number of MC cycles
	void set_perturb_cycles(core::Size const perturb_cycles);
	/// @brief do not use fragments?
	void set_no_frags(bool const no_frags);
	/// @brief what minimizer type to use?
	void set_min_type(std::string const & min_type);
	/// @brief show perturb result structure?
	void set_perturb_show(bool const perturb_show);

	/// @brief read in options from the options system
	void read_options();

private:

	/// @details AnchorMoversData object holds scorefunctions, etc, for the AnchoredDesign suite (shared between several movers)
	protocols::anchored_design::AnchorMoversDataOP interface_;

	//option system replacement
	/// @brief debugging mode activates a bunch of extra output
	bool debug_;
	/// @brief do not perform CCD style closure (use KIC only)
	bool perturb_CCD_off_;
	/// @brief do not perform KIC style closure (use CCD only)
	bool perturb_KIC_off_;
	/// @brief use nonpivot torsion sampling for KIC?
	bool nonpivot_torsion_sampling_;
	/// @brief MC temperature
	core::Real perturb_temp_;
	/// @brief number of MC cycles
	core::Size perturb_cycles_;
	/// @brief do not use fragments?
	bool no_frags_;
	/// @brief what minimizer type to use?
	std::string min_type_;
	/// @brief show perturb result structure?
	bool perturb_show_;
}; //class AnchoredPerturbMover

class AnchoredRefineMover : public protocols::moves::Mover {

public:
	//@brief constructor with arguments
	AnchoredRefineMover( protocols::anchored_design::AnchorMoversDataOP interface_in );

	~AnchoredRefineMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	//option system replacement
	/// @brief debugging mode activates a bunch of extra output
	bool get_debug() const;
	/// @brief do not perform CCD style closure (use KIC only)
	bool get_refine_CCD_off() const;
	/// @brief do not perform KIC style closure (use CCD only)
	bool get_refine_KIC_off() const;
	/// @brief use nonpivot torsion sampling for KIC?
	bool get_nonpivot_torsion_sampling() const;
	/// @brief KIC use vicinity sampling?
	bool get_vicinity_sampling() const;
	/// @brief KIC vicinity sampling degrees
	core::Real get_vicinity_degree() const;
	/// @brief MC temperature
	core::Real get_refine_temp() const;
	/// @brief number of MC cycles
	core::Size get_refine_cycles() const;
	/// @brief what minimizer type to use?
	std::string const & get_min_type() const;
	/// @brief how many cycles between repack/design opportunities?
	core::Size get_refine_repack_cycles() const;

	/// @brief debugging mode activates a bunch of extra output
	void set_debug(bool const debug);
	/// @brief do not perform CCD style closure (use KIC only)
	void set_refine_CCD_off(bool const refine_CCD_off);
	/// @brief do not perform KIC style closure (use CCD only)
	void set_refine_KIC_off(bool const refine_KIC_off);
	/// @brief use nonpivot torsion sampling for KIC?
	void set_nonpivot_torsion_sampling(bool const nonpivot_torsion_sampling);
	/// @brief KIC use vicinity sampling?
	void set_vicinity_sampling(bool const vicinity_sampling);
	/// @brief KIC vicinity sampling degrees
	void set_vicinity_degree(core::Size const vicinity_degree);
	/// @brief MC temperature
	void set_refine_temp(core::Real const refine_temp);
	/// @brief number of MC cycles
	void set_refine_cycles(core::Size const refine_cycles);
	/// @brief what minimizer type to use?
	void set_min_type(std::string const & min_type);
	/// @brief how many cycles between repack/design opportunities?
	void set_refine_repack_cycles(core::Size const refine_repack_cycles);

	/// @brief read in options from the options system
	void read_options();

private:

	protocols::anchored_design::AnchorMoversDataOP interface_;

	//option system replacement
	/// @brief debugging mode activates a bunch of extra output
	bool debug_;
	/// @brief do not perform CCD style closure (use KIC only)
	bool refine_CCD_off_;
	/// @brief do not perform KIC style closure (use CCD only)
	bool refine_KIC_off_;
	/// @brief use nonpivot torsion sampling for KIC?
	bool nonpivot_torsion_sampling_;
	/// @brief KIC use vicinity sampling?
	bool vicinity_sampling_;
	/// @brief KIC vicinity sampling degrees
	core::Real vicinity_degree_;
	/// @brief MC temperature
	core::Real refine_temp_;
	/// @brief number of MC cycles
	core::Size refine_cycles_;
	/// @brief what minimizer type to use?
	std::string min_type_;
	/// @brief how many cycles between repack/design opportunities?
	core::Size refine_repack_cycles_;

}; //class AnchoredRefineMover


}//AnchoredDesign
}//protocols

#endif //INCLUDED_protocols_AnchoredDesign_AnchorMovers_HH
