// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/moves/InterfaceAnalyzerMover.hh
/// @brief InterfaceAnalyzerMover examines interface structures and packages extra scores into a Job object
/// @author Steven Lewis, Bryan Der, Ben Stranges, Jared Adolf-Bryfogle


#ifndef INCLUDED_protocols_analysis_InterfaceAnalyzerMover_HH
#define INCLUDED_protocols_analysis_InterfaceAnalyzerMover_HH

// Unit Headers
#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh> //This non-obviously needs to be the full header because the constructors have NULL default values for their ScoreFunctionCOP arguments.  If you try to use the default ctor with no arguments, it looks like you don't need ScoreFunction.hh, but gcc insists on it anyway.  Putting it here is the cleanest and most documentable solution.
#include <protocols/moves/Mover.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/DockingPartners.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <set>
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

#include <core/id/AtomID.fwd.hh>

#include <map> // AUTO IWYU For map

#ifdef WIN32
#undef interface
#undef set_interface
#endif

namespace protocols {
namespace analysis {

enum InterfaceRegion {
	total = 1,
	side1,
	side2

};

/// @brief All per residue interface data and residue averages.
///  Interface Residues only
/// @details Vector1 correspond to residues in the pose.  avg vector1 sets correspond to InterfaceRegion enums: total, side1, side
struct PerResidueInterfaceData {
	utility::vector1< bool > interface_residues;

	utility::vector1< core::Real > separated_sasa;
	utility::vector1< core::Real > complexed_sasa;
	utility::vector1< core::Real > dSASA;
	utility::vector1< core::Real > dSASA_sc;
	/// @brief 'Hydrophobic' dSASA
	utility::vector1< core::Real > dhSASA;
	utility::vector1< core::Real > dhSASA_sc;


	/// @brief Relative Hydrophobic dSASA.  Calculated by: rel_dSASA = atom_dSASA*(1 - atom_charge).
	///Subtracting this by the real dSASA gives you the relative polar dSASA.
	utility::vector1< core::Real > dhSASA_rel_by_charge;

	utility::vector1< core::Real > SASA;
	utility::vector1< core::Real > dSASA_fraction; //fraction of SASA separated to SASA buried.  If all of it is buried, fraction is 1.0 (dSASA[i]/separated_sasa[i])


	//Complex/Separated energies of Interface residues only.
	utility::vector1< core::Real > separated_energy;
	utility::vector1< core::Real > complexed_energy;
	utility::vector1< core::Real > dG;


	/// @brief Average per residue change in energy in each InterfaceRegion
	utility::vector1<core::Real> regional_avg_per_residue_dG;
	/// @brief Average per residue energy of the complexed interface in each InterfaceRegion
	utility::vector1<core::Real> regional_avg_per_residue_energy_int;
	/// @brief Average per residue energy of the separated interface in each InterfaceRegion
	utility::vector1<core::Real> regional_avg_per_residue_energy_sep;

	/// @brief Average per residue change in solvent accessible surface area in each InterfaceRegion
	utility::vector1<core::Real> regional_avg_per_residue_dSASA;
	/// @brief Average per residue SASA of the complexed interface in each InterfaceRegion
	utility::vector1<core::Real> regional_avg_per_residue_SASA_sep;
	utility::vector1<core::Real> regional_avg_per_residue_SASA_int;

};

/// @brief All interface data. Unless otherwise specified, they refer specifically to the interface
/// @details All vectors correspond to InterfaceRegion enums
struct InterfaceData {

	/// @brief Number of residues in each InterfaceRegion
	utility::vector1< core::Size > interface_nres;

	/// @brief Boolean of interface residues [side1][2]
	utility::vector1< utility::vector1< bool > > interface_residues;

	/// @brief Change in Solvent Accessible Surface Area upon complexion in each InterfaceRegion
	utility::vector1< core::Real > dSASA;
	utility::vector1< core::Real > dSASA_sc;

	/// @brief 'Hydrophobic' dSASA
	utility::vector1< core::Real > dhSASA;
	utility::vector1< core::Real > dhSASA_sc;

	/// @brief Relative Hydrophobic dSASA.  Calculated by: sum(atom_dSASA*(1 - atom_charge)).
	/// Note that this includes hydrogens, and subtracting this by the real dSASA gives you the relative polar dSASA.
	utility::vector1< core::Real > dhSASA_rel_by_charge;


	core::Real complexed_SASA;
	core::Real separated_SASA;

	/// @brief Change in energy upon complexion in each InterfaceRegion
	utility::vector1<core::Real > dG;
	core::Real gly_dG; //Energy of a all Gly interface

	/// @brief Centroid_dG of interface
	core::Real centroid_dG;

	core::Real dG_dSASA_ratio;

	/// @brief Number of unsaturated hbonds in complex
	core::Size delta_unsat_hbonds;

	/// @brief Total energy of interface Hbonds
	core::Real total_hb_E;

	core::Real hbond_E_fraction;

	core::Size interface_hbonds;

	/// @brief Shape Complementarity value
	core::Real sc_value;

	/// @brief PackStat value of the interface
	core::Real packstat;

	///@brief Difference of a RealMetric between complexed and separated pose
	utility::vector1< core::Real > delta_metrics;


	/// @brief Total energy of the complex in each InterfaceRegion
	utility::vector1< core::Real > complex_total_energy;

	/// @brief Total energy of the separated pose in each InterfaceRegion
	utility::vector1< core::Real > separated_total_energy;

	core::Real crossterm_interface_energy;
	core::Real crossterm_interface_energy_dSASA_ratio;

	/// @brief pymol style selections
	std::string pymol_sel_interface;
	std::string pymol_sel_hbond_unsat;
	std::string pymol_sel_packing;

	utility::vector1< core::Real > complexed_interface_score;
	utility::vector1< core::Real > separated_interface_score;

	utility::vector1< core::Size > aromatic_nres;
	utility::vector1< core::Real > aromatic_dSASA_fraction;

	/// @brief Aromatic contribution to dG in each InterfaceRegion
	utility::vector1< core::Real > aromatic_dG_fraction;

	utility::vector1< core::Size > ss_helix_nres;
	utility::vector1< core::Size > ss_loop_nres;
	utility::vector1< core::Size > ss_sheet_nres;

	/// @brief Fraction of interface nres to total surface residues in separated pose for total, side1 and side2.
	utility::vector1< core::Real > interface_to_surface_fraction;


};


/// @brief Class for analyzing interfaces of a pose.  Many metrics are calculated and accessible after the apply method.
class InterfaceAnalyzerMover : public protocols::moves::Mover {

public:
	//stole these directly from InterGroupNeighborsCalculator
	typedef std::set< core::Size > one_group;
	typedef std::pair< one_group, one_group > group_pair;
	typedef utility::vector1< group_pair > group_set;

	// I think this ctor is also the 'default ctor' because it can take no args
	/// @brief Constructor for 2 chain poses, separates them by jump number
	/// @details pack_separated and pack_input only pack the detected interface residues.
	InterfaceAnalyzerMover(
		core::Size interface_jump = 1,
		bool const tracer = false,
		core::scoring::ScoreFunctionCOP sf = nullptr,
		bool compute_packstat = false,
		bool pack_input = false,
		bool pack_separated = false,
		bool use_jobname = true,
		bool detect_disulfide_in_separated_pose = true
	);

	/// @brief Constructor for poses with >= 1 jump, keeps defined chains together
	/// @details pack_separated and pack_input only pack the detected interface residues.
	InterfaceAnalyzerMover(
		std::set<int> const & fixed_chains,
		bool const tracer = false,
		core::scoring::ScoreFunctionCOP sf = nullptr,
		bool compute_packstat = false,
		bool pack_input = false,
		bool pack_separated = false,
		bool use_jobname = true,
		bool detect_disulfide_in_separated_pose = true
	);

//	/// @brief Constructor for any interface in a pose.  Uses string designation (ex LH_A) to keep the left chain/chains fixed ALA docking.
//	/// @details Can be used for subsets of interfaces, for example L_H in a LHA pose.
//	/// pack_separated and pack_input only pack the detected interface residues.
//	InterfaceAnalyzerMover(
//		std::string dock_chains,
//		bool const tracer = false,
//		core::scoring::ScoreFunctionCOP sf = nullptr,
//		bool compute_packstat = false,
//		bool pack_input = false,
//		bool pack_separated = false,
//		bool use_jobname = true,
//		bool detect_disulfide_in_separated_pose = true
//	);

	/// @brief Constructor for any interface in a pose.  Uses string designation (ex LH_A) to keep the left chain/chains fixed ALA docking.
	/// @details Can be used for subsets of interfaces, for example L_H in a LHA pose.
	/// pack_separated and pack_input only pack the detected interface residues.
	InterfaceAnalyzerMover(
		core::pose::DockingPartners const & dock_chains,
		bool const tracer = false,
		core::scoring::ScoreFunctionCOP sf = nullptr,
		bool compute_packstat = false,
		bool pack_input = false,
		bool pack_separated = false,
		bool use_jobname = true,
		bool detect_disulfide_in_separated_pose = true
	);

	~InterfaceAnalyzerMover() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	//virtual bool
	//reinitialize_for_new_input() const {return true;}

	/// @brief Explicitly initialize settings on apply - not at the constructor, since this can hold state, and some protocols use apply multiple times.
	void
	init_on_new_input(const core::pose::Pose & pose);

	/// @brief Called by MoverFactory when constructing new Movers. Takes care of the specific mover's parsing.

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &
	) override;

	/// @brief apply function will calculate data about the input pose.  It is not intended to modify the pose itself (conformation and energies objects).
	/// It will add data to the pose for output into a scorefile.
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Apply method for const pose.  Used by InterfaceFeatures reporter
	virtual void
	apply_const(core::pose::Pose const & pose);


	///Print data to tracer or scorefile if tracer is not set (default).
	virtual void
	report_data();

	/// @brief Add all the scoring data normally output via job to the pose instead.  Can't add this to method due to const apply.
	/// Use after you have used the apply for the class.
	///
	void
	add_score_info_to_pose(core::pose::Pose & pose);


	bool reinitialize_for_each_job() const override { return true; }

	bool reinitialize_for_new_input() const override { return true; }


	///////////////////////////////////////////////////
	//Setters for various computations and options.
	//Each computation is expensive, so you can turn them off if desired
	//

	void
	set_defaults();

	void
	set_scorefunction( core::scoring::ScoreFunctionCOP sf );

	void
	set_use_centroid_dG(bool const use_centroid);
	bool
	get_use_centroid_dG() const;

	void
	set_compute_packstat(bool const compute_packstat);

	///@brief Compute the Lawrence and Coleman Shape Complementarity score?
	void
	set_compute_interface_sc(bool const compute_interface_sc) {compute_interface_sc_ = compute_interface_sc;}

	void
	set_compute_separated_sasa(bool const compute_separated_sasa) {compute_separated_sasa_ = compute_separated_sasa;}

	void
	set_calc_dSASA( bool calc_dSASA ){ compute_separated_sasa_ = calc_dSASA; }

	void
	set_compute_interface_energy(bool const iface_en) {compute_interface_energy_ = iface_en;}

	void
	set_calc_hbond_sasaE(bool const calc_hbond_sasaE) {calc_hbond_sasaE_ = calc_hbond_sasaE;}

	void
	set_compute_interface_delta_hbond_unsat(bool const IDHU) {compute_interface_delta_hbond_unsat_ = IDHU;}

	void
	set_detect_disulfide_in_separated_pose( bool const detect_disulfide_in_separated_pose) {detect_disulfide_in_separated_pose_ = detect_disulfide_in_separated_pose;}

	/// @brief Repack the interface of the complex before separation.
	void
	set_pack_input(bool const pack_input);

	/// @brief Repack the interface of the complex after separation.
	void
	set_pack_separated(bool const pack_separated);

	///@brief If pack separated or together is set, this controls how many rounds of packing we do.
	/// Default is 1, but that certainly may not be enough...
	void
	set_pack_rounds( core::Size pack_rounds ){ pack_rounds_ = pack_rounds; }

	void
	set_interface_jump(core::Size const interface_jump);

//	void
//	set_interface( std::string const & interface );

	void
	set_interface( core::pose::DockingPartners const & interface );

	void
	set_skip_reporting(bool const skip_reporting) {skip_reporting_ = skip_reporting;}

	void
	set_use_tracer(bool const tracer);

	void
	set_use_jobname(bool const use_jobname);

	void
	set_scorefile_reporting_prefix( std::string const & prefix );

	//void set_calcs_ready(bool const calcs_ready); //why is this externally accessible?

	void set_use_resfile(bool const use_resfile); //does not do anything.
	//bool get_use_resfile() const;

	void set_metrics( utility::vector1<core::simple_metrics::RealMetricCOP> const & metrics );

	///////////////////////////////////////////////////
	//Getters for the various parameters used by the analyzer
	//

	/// @brief Get all interface data.
	InterfaceData
	get_all_data() {
		return data_;
	};

	/// @brief Get all per residue interface data.
	PerResidueInterfaceData
	get_all_per_residue_data() {
		return per_residue_data_;
	};


	core::Size
	get_num_interface_residues();

	core::Real
	get_complexed_sasa();

	core::Real
	get_interface_delta_sasa();


	core::Real
	get_complex_energy();

	/// @brief return the average per residue interface energy
	core::Real
	get_per_residue_energy();

	core::Real
	get_total_Hbond_E();

	core::Real
	get_separated_interface_energy();

	core::Real
	get_crossterm_interface_energy();

	core::Real
	get_separated_interface_energy_ratio();

	core::Real
	get_crossterm_interface_energy_ratio();


	/// @brief Return the interface energy of an all glycine interface (like bb-bb energies)
	core::Real
	get_gly_interface_energy();


	core::Real
	get_interface_packstat();

	core::Size
	get_interface_delta_hbond_unsat();


	core::Real
	get_side1_score() { return data_.complexed_interface_score[side1]; }

	core::Real
	get_side2_score() { return data_.complexed_interface_score[side2]; }

	core::Size
	get_side1_nres() { return data_.interface_nres[side1]; }

	core::Size
	get_side2_nres() { return data_.interface_nres[side2]; }


	std::string
	get_pymol_sel_interface();

	std::string
	get_pymol_sel_hbond_unsat();

	std::string
	get_pymol_sel_packing();


	core::Real
	get_interface_dG() const;

	/// @brief Return boolean if a multichain (fixedchain) constructor was used
	bool
	get_multichain_constructor();

	std::set<int>
	get_fixed_chains();

	/// @brief Get the residues at the interface in question
	std::set<core::Size>
	get_interface_set();

	group_set
	get_chain_groups();

	bool
	get_pack_input();

	/// @brief
	core::Real
	get_centroid_dG();

	utility::vector1<core::simple_metrics::RealMetricCOP>
	get_metrics();

	std::string
	get_name() const override;

	std::string
	get_scorefile_reporting_prefix() const;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/// @brief
	// Undefined, commenting out to fix PyRosetta build  core::Real get_interface_Hbond_sasa();

	/// @brief the exposure/possible ratio avg for hbonds in the interface
	//core::Real get_Hbond_exposure_ratio();
	/// @brief total hbond energy for pose


private:

	/// @brief registers the posemetric calculators
	void register_calculators();

	/// @brief sets up the packer task.  Used for both sep and tog pose
	core::pack::task::PackerTaskOP
	setup_task(core::pose::Pose & pose);

	/// @brief functions to make interface sets needed
	virtual void
	make_multichain_interface_set( core::pose::Pose & Pose, std::set<int> & fixed_chains );

	/// @brief makes the interface sets for either constructor and sets up any other basic interface info.  (chains, upstream chains, etc.)
	virtual void
	make_interface_set( core::pose::Pose & pose );

	/// @brief Creates a single pose that is separated structurally at the interface
	/// @details Uses the jump number to find downstream residues and move them in the Z axis in cartesian space according to step_size (A)
	virtual core::pose::Pose
	make_separated_pose( core::pose::Pose & pose, core::Size interface_jump, core::Size step_size = 1000 );

	/// @brief Creates a single pose that is separated structurally at the interface
	/// @details Moves the moving_chains in cartesian space a distance according to step_size (A) in the Z axis
	virtual core::pose::Pose
	make_separated_pose( core::pose::Pose & pose, std::set< core::Size > const & moving_chains, core::Size step_size = 1000 );

	/// @brief reorder the fold tree to allow multichain interfaces to be evaluated
	///       returns the jump number to use to define the interface
	virtual core::Size
	reorder_foldtree_find_jump( core::pose::Pose & pose, std::set<int> & fixed_chains);


	///////////////////////////////////////////////////
	//Interface Computations
	//

	/// @brief find the interface shape compementarity value between the chains using old fortran SC score.
	/// @details (From the wiki) Calculates the Lawrence & Coleman shape complementarity using a port of the original Fortran code from CCP4's sc.
	///
	void
	compute_interface_sc( core::Size & interface_jump, core::pose::Pose const & complexed_pose);

	virtual void
	compute_separated_sasa(core::pose::Pose & complexed_pose, core::pose::Pose & separated_pose);

	virtual void
	compute_interface_energy(core::pose::Pose & complexed_pose, core::pose::Pose & separated_pose);

	virtual void
	compute_interface_packstat( core::pose::Pose & pose );

	virtual void
	compute_interface_delta_hbond_unsat(core::pose::Pose & complexed_pose, core::pose::Pose & separated_pose );

	virtual void
	score_separated_chains(core::pose::Pose & complexed_pose, core::pose::Pose & separated_pose);

	///@brief Calculates the difference of a RealMetric in the complexed vs seperated pose
	virtual void
	compute_delta_metric(core::pose::Pose const & complexed_pose, core::pose::Pose const & separeted_pose);


	/// @brief Calculate the average energy per residue in the interface as well as other data.
	void
	calc_per_residue_and_regional_data( core::pose::Pose & complexed_pose, core::pose::Pose & separated_pose);

	/// @brief calculates dSASA for each residue
	utility::vector1< core::Real >
	calc_per_residue_dSASA(const core::pose::Pose & complexed_pose, const utility::vector1< core::Real >  & separated_sasa, const utility::vector1< core::Real >&  complexed_sasa);

	core::Real
	calc_per_residue_dSASA_general(
		const core::Size resnum,
		const core::pose::Pose & complexed_pose,
		const core::Real separated_sasa,
		const core::Real complexed_sasa,
		utility::vector1< core::Real> & regional);

	utility::vector1< core::Real >
	calc_per_residue_dG(core::pose::Pose & complexed_pose, const utility::vector1< core::Real > & separated_energy, const utility::vector1< core::Real > & complexed_energy);

	void
	calc_hbond_sasaE( core::pose::Pose & pose );

	/// @brief report the dG of a centroid pose with score3
	void
	calc_centroid_dG ( core::pose::Pose const & complex_pose, core::pose::Pose const & separated_pose );

	/// @brief Calculate the number of residues in the interface vs the surface based on SASA cutoff (40A^2 as per LayerDesign and selection operations).
	/// @details Correct way to calculate this though is through % maximal SASA buried, but calculations will need to be done later to find these values.
	void
	calc_interface_to_surface_fraction(core::pose::Pose const & separated_pose, const utility::vector1<core::Real> & separated_sasa);

	void print_pymol_selection_of_interface_residues( core::pose::Pose const & pose, std::set< core::Size > const & interface_set );
	void print_pymol_selection_of_hbond_unsat( core::pose::Pose & pose, utility::vector1< core::id::AtomID > delta_unsat_hbond_atid_vector );
	void print_pymol_selection_of_packing( core::pose::Pose const & pose, utility::vector1< core::Real > & interface_pack_scores );


	/// @brief mutate all residue in the interface to Gly and recalc the energy - not used right now
	void mut_to_gly( core::pose::Pose complex_pose, core::pose::Pose separated_pose );

	/// @brief sets up the pose information such as the name and chain ids
	void
	set_pose_info( core::pose::Pose const & pose );

	/// @brief Initialize the per residue data structure
	void
	init_per_residue_data( core::pose::Pose const & pose);

	/// @brief Initialize the data structure
	void
	init_data(const core::pose::Pose & pose);

	/// @brief Setup the scorefunction to include hbond energies in the EnergyGraph.  Yay forums and Rocco for this bug fix!
	void
	setup_scorefxn();

	/// @brief Setup for the dock_chains constructor
	void
	setup_for_dock_chains(core::pose::Pose & pose, core::pose::DockingPartners const & partners);

	std::string
	scorefile_column_name(std::string const & base_name) const;

	void
	setup_score_data();

private:
	/// @brief jump to define which interface is interesting
	core::Size interface_jump_;
	/// @brief what chains are fixed in an interface
	utility::vector1< std::string > fixed_chain_strings_; // Will be transfered to fixed_chains_ at the head of apply()
	std::set<int> fixed_chains_;

	/// @brief the ligand chain that will be moved in an interface
	///In this scheme, everything else is fixed
	std::string ligand_chain_;
	/// @brief scorefunction
	core::scoring::ScoreFunctionOP sf_;

	utility::file::FileName posename_;
	std::string posename_base_;

	///@brief a RealMetric to calculate the difference between complexed and single state for
	utility::vector1<core::simple_metrics::RealMetricCOP> delta_metrics_;

	core::Size chain1_ = 1;
	core::Size chain2_ = 2;
	core::pose::DockingPartners dock_chains_;
	std::set<core::Size> upstream_chains_; //fixed chains
	std::set<core::Size> downstream_chains_; //other chains
	std::set<core::Size> ignored_chains_; //Ignored chains from the dock_chains constructor.  These will be a part of fixed chains.
	/// @brief output to tracer or PDB/silent file
	bool tracer_;
	/// @brief are calculators ready?
	bool calcs_ready_;

	/// @brief bother with computing packstat
	bool compute_packstat_;
	/// @brief bother with computing interface sc
	bool compute_interface_sc_;

	/// @brief skip this expensive calculation
	bool compute_separated_sasa_;
	/// @brief skip this expensive calculation
	bool compute_interface_energy_;
	/// @brief skip this expensive calculation
	bool calc_hbond_sasaE_;
	/// @brief skip this expensive calculation
	bool compute_interface_delta_hbond_unsat_;

	/// @brief hush!  be quiet!  silence!
	bool skip_reporting_;

	//bool is_compute_hbond_unsat_;
	/// @brief which constructor are we using
	bool explicit_constructor_;
	/// @brief pack the input pose
	bool pack_input_;
	/// @brief pack the separated poses default is false
	bool pack_separated_;
	/// @brief just the jobname (true value) or the pose name (flase value)
	bool use_jobname_;
	/// @brief use a resfile during the pack_input and pack_separated operations
	bool use_resfile_;
	/// @brief skip the centroid_dG step, for incoming poses not centroid convertible
	bool use_centroid_;
	/// @brief detect disulfide in the separated pose
	bool detect_disulfide_in_separated_pose_;

	InterfaceData data_;
	PerResidueInterfaceData per_residue_data_;

	/// @brief avg hbond exposure ratio
	//Real hbond_exposure_ratio_;
	//Real total_hb_sasa_;

	core::Size included_nres_;
	/// @brief set of residues at the interface in question
	std::set< core::Size > interface_set_;
	utility::vector1< bool > include_residue_; //All residues in the pose minus any that are ignored from dock_chains constructor, or any other function that changes this value.

	/// @brief group of residue ids of fixed chains and mobile chains (see typedef)
	group_set chain_groups_;


	//name strings for PoseMetricCalculators
	/// @brief Sasa calculator name string
	std::string Sasa_;
	/// @brief InterfaceNeighborDefinition calculator name string
	std::string InterfaceNeighborDefinition_;
	/// @brief InterfaceDeltaEnergetics calculator name string
	std::string InterfaceDeltaEnergetics_;
	/// @brief NumberHBonds calculator name string
	std::string NumberHBonds_;
	/// @brief BuriedUnsatisfiedPolars calculator name string
	std::string BuriedUnsatisfiedPolars_;

	std::map< std::string, core::Real > score_data_;

	// The prefix to use when reporting scores to the score file
	// via the PoseExtraScores.
	std::string scorefile_reporting_prefix_;

	core::Size pack_rounds_ = 1;

}; //class InterfaceAnalyzerMover

}//analysis
}//protocols

#endif //INCLUDED_protocols_analysis_InterfaceAnalyzerMover_HH
