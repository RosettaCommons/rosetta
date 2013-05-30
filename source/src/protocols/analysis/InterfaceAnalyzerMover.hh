// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/moves/InterfaceAnalyzerMover.hh
/// @brief InterfaceAnalyzerMover examines interface structures and packages extra scores into a Job object
/// @author Steven Lewis, Bryan Der, Ben Stranges

#ifndef INCLUDED_protocols_analysis_InterfaceAnalyzerMover_HH
#define INCLUDED_protocols_analysis_InterfaceAnalyzerMover_HH

// Unit Headers
#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh> //This non-obviously needs to be the full header because the constructors have NULL default values for their ScoreFunctionCOP arguments.  If you try to use the default ctor with no arguments, it looks like you don't need ScoreFunction.hh, but gcc insists on it anyway.  Putting it here is the cleanest and most documentable solution.
#include <protocols/moves/Mover.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility Headers
// AUTO-REMOVED #include <core/init/init.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <set>

#include <core/id/AtomID.fwd.hh>



namespace protocols {
namespace analysis {

class InterfaceAnalyzerMover : public protocols::moves::Mover {

public:
	//stole these directly from InterGroupNeighborsCalculator
	typedef std::set< core::Size > one_group;
	typedef std::pair< one_group, one_group > group_pair;
	typedef utility::vector1< group_pair > group_set;

	//constructor for 2 chain poses, separates them by jump number
	InterfaceAnalyzerMover(
		core::Size interface_jump = 1,
		bool const tracer = false,
		core::scoring::ScoreFunctionCOP sf = NULL,
		bool compute_packstat = false,
		bool pack_input = false,
		bool pack_separated = false,
		bool use_jobname = true
	);
	//constructor for poses with > 1 jump, keeps defined chains together
	InterfaceAnalyzerMover(
		std::set<int> fixed_chains,
		bool const tracer = false,
		core::scoring::ScoreFunctionCOP sf = NULL,
		bool compute_packstat = false,
		bool pack_input = false,
		bool pack_separated = false,
		bool use_jobname = true
	);

	virtual ~InterfaceAnalyzerMover();

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	/// @brief Called by MoverFactory when constructing new Movers. Takes care of the specific mover's parsing.
	virtual
	void parse_my_tag(
		utility::tag::TagPtr const,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	///@brief apply function will calculate data about the input pose.  It is not intended to modify the pose itself (conformation and energies objects) although it may toss data into the DataCache or a Job object.
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void set_pose_info( core::pose::Pose const & pose );
	/* Undefined, commenting out to fix PyRosetta build   void set_scorefunction( core::scoring::ScoreFunctionCOP sf ); */
	void use_tracer( bool const tracer ) { tracer_ = tracer; }

	///@brief funtions to make interface sets needed
	virtual void make_multichain_interface_set( core::pose::Pose & Pose,
																							std::set<int> & fixed_chains );
	virtual void make_interface_set( core::pose::Pose & pose );

	///@brief assigns the complexed and separated poses for the entire mover
	virtual core::pose::Pose make_separated_pose( core::pose::Pose & pose );

	///@brief reorder the fold tree to allow multichain interfaces to be evaluated
	///       returns the jump number to use to define the interface
	virtual core::Size reorder_foldtree_find_jump( core::pose::Pose & pose,
																								 std::set<int> & fixed_chains);

	//functions to compute various parameters used by the analyzer
	//these should be private...
	virtual	void compute_separated_sasa(core::pose::Pose & complexed_pose,
																			core::pose::Pose & separated_pose);
	virtual	void compute_interface_energy(core::pose::Pose & complexed_pose,
																				core::pose::Pose & separated_pose);
	virtual	void compute_interface_packstat( core::pose::Pose & pose );
	virtual	void compute_interface_delta_hbond_unsat(core::pose::Pose & complexed_pose,
																									 core::pose::Pose & separated_pose );
	virtual	void score_separated_chains(core::pose::Pose & complexed_pose,
																			core::pose::Pose & separated_pose);
	virtual void report_data();
	///@brief setters for the various parameters used by the analyze

	void set_use_resfile(bool const use_resfile);
	bool get_use_resfile() const;

	void set_use_centroid_dG(bool const use_centroid);
	bool get_use_centroid_dG() const;

	///setters for various computations - each of these are expensive, so you can turn them off if desired
	void set_compute_packstat(bool const compute_packstat);
	void set_compute_interface_sc(bool const compute_interface_sc) {compute_interface_sc_ = compute_interface_sc;}
	void set_compute_separated_sasa(bool const compute_separated_sasa) {compute_separated_sasa_ = compute_separated_sasa;}
	void set_compute_interface_energy(bool const iface_en) {compute_interface_energy_ = iface_en;}
	void set_calc_hbond_sasaE(bool const calc_hbond_sasaE) {calc_hbond_sasaE_ = calc_hbond_sasaE;}
	void set_compute_interface_delta_hbond_unsat(bool const IDHU) {compute_interface_delta_hbond_unsat_ = IDHU;}

	void set_pack_input(bool const pack_input);
	void set_pack_separated(bool const pack_separated);

	void set_interface_jump(core::Size const interface_jump);

	void set_skip_reporting(bool const skip_reporting) {skip_reporting_ = skip_reporting;}
	void set_tracer(bool const tracer);
	void set_use_jobname(bool const use_jobname);

	//void set_calcs_ready(bool const calcs_ready); //why is this externally accessible?

	///@brief getters for the various parameters used by the analyzer
	core::Real get_total_sasa();
	core::Real get_interface_delta_sasa();
	core::Size get_num_interface_residues();
	core::Real get_complex_energy();
	core::Real get_per_residue_energy();
	core::Real get_separated_interface_energy();
	core::Real get_crossterm_interface_energy();
	core::Real get_separated_interface_energy_ratio();
	core::Real get_crossterm_interface_energy_ratio();
	core::Real get_interface_packstat();
	core::Size get_interface_delta_hbond_unsat();
	core::Real get_side1_score() { return side1_score_; }
	core::Real get_side2_score() { return side2_score_; }
	core::Size get_side1_nres() { return side1_nres_; }
	core::Size get_side2_nres() { return side2_nres_; }
	std::string get_pymol_sel_interface();
	std::string get_pymol_sel_hbond_unsat();
	std::string get_pymol_sel_packing();
	core::Real get_interface_ddG() const;
	bool get_multichain_constructor();
	std::set<int> get_fixed_chains();
	std::set<core::Size> get_interface_set();
	group_set get_chain_groups();
	bool get_pack_input();
	core::pack::task::PackerTaskOP get_packer_task();
	///@brief return the interface energy of an all glycine interface (like bb-bb energies)
	core::Real get_gly_interface_energy();
	///@brief
	core::Real get_centroid_dG();

	///@brief
	// Undefined, commenting out to fix PyRosetta build  core::Real get_interface_Hbond_sasa();

	///@brief the exposure/possible ratio avg for hbonds in the interface
	//core::Real get_Hbond_exposure_ratio();
	///@brief total hbond energy for pose
	core::Real get_total_Hbond_E();
private:
	///@brief registers the posemetric calculators
	void register_calculators();
	void register_intergroup_calculator();
	///@brief sets up the packer task within the mover
	void setup_task(core::pose::Pose & pose);
	void print_pymol_selection_of_interface_residues( core::pose::Pose const & pose, std::set< core::Size > const interface_set );
	void print_pymol_selection_of_hbond_unsat( core::pose::Pose & pose, utility::vector1< core::id::AtomID > delta_unsat_hbond_atid_vector );
	void print_pymol_selection_of_packing( core::pose::Pose const & pose,	utility::vector1< core::Real > interface_pack_scores );
	///@brief calculate the average energy per residue in the interface
	void calc_per_residue_energy( core::pose::Pose pose);
	///@brief mutate all residue in the interface to Gly and recalc the energy - not used right now
	void mut_to_gly( core::pose::Pose complex_pose, core::pose::Pose separated_pose );
	///@brief report the dG of a centroid pose with score3
	void calc_centroid_dG ( core::pose::Pose complex_pose, core::pose::Pose separated_pose );
	///@brief fill in later
	void calc_hbond_sasaE( core::pose::Pose pose );
	///@brief find the interface shape compementarity value between the chains
	void compute_interface_sc( core::Size & interface_jump, core::pose::Pose const & complexed_pose);
private:
	///@brief jump to define which interface is interesting
	core::Size interface_jump_;
	//@brief what chains are fixed in an interface
	std::set<int> fixed_chains_;
	///@brief scorefunction
	core::scoring::ScoreFunctionCOP sf_;

	utility::file::FileName posename_;
	std::string posename_base_;

	core::Size chain1_;
	core::Size chain2_;
	std::string chain1_char_;
	std::string chain2_char_;
	std::set<core::Size> upstream_chains_;
	std::set<core::Size> downstream_chains_;

	///@brief output to tracer or PDB/silent file
	bool tracer_;
	///@brief are calculators ready?
	bool calcs_ready_;

	///@brief bother with computing packstat
	bool compute_packstat_;
	///@brief bother with computing interface sc
	bool compute_interface_sc_;

	///@brief skip this expensive calculation
	bool compute_separated_sasa_;
	///@brief skip this expensive calculation
	bool compute_interface_energy_;
	///@brief skip this expensive calculation
	bool calc_hbond_sasaE_;
	///@brief skip this expensive calculation
	bool compute_interface_delta_hbond_unsat_;

	///@brief hush!  be quiet!  silence!
	bool skip_reporting_;

	//bool is_compute_hbond_unsat_;
	///@brief which constructor are we using
	bool multichain_constructor_;
	///@brief pack the input pose
	bool pack_input_;
	///@brief pack the separated poses default is false
	bool pack_separated_;
	///@brief just the jobname (true value) or the pose name (flase value)
	bool use_jobname_;
	///@brief use a resfile during the pack_input and pack_separated operations
	bool use_resfile_;
	///@brief skip the centroid_dG step, for incoming poses not centroid convertible
	bool use_centroid_;

	///@brief interface energy
	core::Real ddG_;
	core::Real total_sasa_;
	core::Real interface_delta_sasa_;
	core::Real interface_hsasa_;
	core::Real interface_polar_sasa_;
	//number of interface residues
	core::Size n_interface_res_;
	//energy of the unseparated complex (just the score)
	core::Real per_residue_energy_;
	core::Real complex_energy_;
	core::Real separated_interface_energy_;
	core::Real crossterm_interface_energy_;
	core::Real separated_interface_energy_ratio_;
	core::Real crossterm_interface_energy_ratio_;
	core::Real interface_packstat_;
	///@brief Energy of a all Gly interface
	core::Real gly_dG_;
	///@brief centroid_dG of interface
	core::Real centroid_dG_;
	///@brief number of unsat hbonds in complex
	core::Size delta_unsat_hbond_counter_;
	///@brief pymol style selections
	std::string pymol_sel_interface_;
	std::string pymol_sel_hbond_unsat_;
	std::string pymol_sel_packing_;

	core::Real side1_score_;
	core::Real side2_score_;
	core::Size nres_;
	core::Size side1_nres_;
	core::Size side2_nres_;

	///@brief avg hbond exposure ratio
	//core::Real hbond_exposure_ratio_;
	//core::Real total_hb_sasa_;
	///@brief total energy of interface Hbonds
	core::Real total_hb_E_;

	///@breif shape complementarity values
	core::Real sc_value_;

	///@brief set of residues at the interface in question
	std::set< core::Size > interface_set_;
	///@brief group of residue ids of fixed chains and mobile chains (see typedef)
	group_set chain_groups_;
	///@brief packer task used to repack pulled apart chains
	core::pack::task::PackerTaskOP task_;
	//name strings for PoseMetricCalculators
	///@brief Sasa calculator name string
	std::string Sasa_;
	///@brief InterfaceNeighborDefinition calculator name string
	std::string InterfaceNeighborDefinition_;
	///@brief InterfaceSasaDefinition calculator name string
	std::string InterfaceSasaDefinition_;
	///@brief InterfaceDeltaEnergetics calculator name string
	std::string InterfaceDeltaEnergetics_;
	///@brief NumberHBonds calculator name string
	std::string NumberHBonds_;
	///@brief BuriedUnsatisfiedPolars calculator name string
	std::string BuriedUnsatisfiedPolars_;
///@brief InterGroupNeighborsCalculator calculator name string
	std::string InterGroupNeighborsCalculator_;
}; //class InterfaceAnalyzerMover

}//analysis
}//protocols

#endif //INCLUDED_protocols_analysis_InterfaceAnalyzerMover_HH
