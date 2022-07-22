// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/protein_interface_design/movers/SecretionOptimizationMover.hh
///
/// @brief AKA the Degreaser, this Mover calculates transmembrane insertion potential (dG_ins) at the sequence level, then quasi-exhaustively attempts to increase dG_ins if it falls below a user-given threshold. By default, makes three or one mutation that most increases the lowest dG_ins stretch in the protein (while remaining under user-given score threshold).
/// @author John Wang

#ifndef INCLUDED_protocols_protein_interface_design_movers_SecretionOptimizationMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_SecretionOptimizationMover_hh

#include <protocols/protein_interface_design/movers/SecretionOptimizationMover.fwd.hh>
#include <protocols/moves/MoveMapMover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/vector1.hh>

#include <core/scoring/ScoreType.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/energy_methods/PoissonBoltzmannEnergy.hh>

#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <string>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/pack/task/TaskFactory.hh>

#include <vector>
#include <unordered_map>
#include <map>
#include <protocols/minimization_packing/PackRotamersMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {


class SecretionOptimizationMover : public protocols::moves::Mover
{
public:

	typedef core::Real Real;
	typedef core::scoring::ScoreType ScoreType;
	typedef core::pose::Pose Pose;

public :
	SecretionOptimizationMover();
	SecretionOptimizationMover( core::scoring::ScoreFunctionCOP scorefxn );
	void apply (Pose &) override;
	~SecretionOptimizationMover() override;
	protocols::moves::MoverOP fresh_instance() const override;// { return (protocols::moves::MoverOP) protocols::moves::MoverOP( new SecretionOptimizationMover ); }
	protocols::moves::MoverOP clone() const override;
	void parse_my_tag(  utility::tag::TagCOP, basic::datacache::DataMap & ) override ;


	void task_factory( core::pack::task::TaskFactoryCOP task_factory ) {
		task_factory_ = task_factory;
	}
	core::pack::task::TaskFactoryCOP task_factory() const {
		return task_factory_;
	}

	//setters and getters for parameters
	void set_max_mutations(core::Size const max_mutations) {
		max_mutations_ = max_mutations;
	}
	core::Size get_max_mutations() {
		return max_mutations_ ;
	}
	void set_aas_allowed(std::string const aas_allowed){
		keep_ = aas_allowed;
	}
	std::string get_aas_allowed(){
		return keep_;
	}
	void set_dG_ins_threshold(core::Real const dG_ins_threshold) {
		dG_ins_threshold_ = dG_ins_threshold;
	}
	core::Real get_dG_ins_threshold() {
		return dG_ins_threshold_;
	}
	void set_score_tolerance(core::Real const score_tolerance) {
		score_tolerance_ = score_tolerance;
	}
	core::Real get_score_tolerance() {
		return score_tolerance_;
	}
	void set_chain_to_degrease(std::string const chain_to_degrease){
		chain_to_degrease_ = chain_to_degrease;
	}
	std::string get_chain_to_degrease() {
		return chain_to_degrease_;
	}
	void set_dG_tolerance(core::Real const dG_tolerance){
		dG_tolerance_ = dG_tolerance;
	}
	core::Real get_dG_tolerance() {
		return dG_tolerance_;
	}
	void set_dump_multis(bool dump_multis){
		dump_multis_ = dump_multis;
	}
	bool get_dump_multis(){
		return dump_multis_;
	}
	void set_dump_singles(bool dump_singles){
		dump_singles_ = dump_singles;
	}
	bool get_dump_singles(){
		return dump_singles_;
	}
	void set_dump_dir(std::string const dump_dir){
		dump_dir_ = dump_dir;
	}
	std::string get_dump_dir() {
		return dump_dir_;
	}
	void use_largest_ddg(bool largest_ddg){
		largest_ddg_ = largest_ddg;
	}
	void use_smallest_dscore(bool smallest_dscore){
		smallest_dscore_ = smallest_dscore;
	}
	bool largest_ddg() {
		return largest_ddg_;
	}
	bool smallest_dscore() {
		return smallest_dscore_;
	}
	//end setters and getters

	void make_table() ; //sets up table of coefficients for dG_ins calculation, based on Hessa et. al 2007

	std::string get_name() const override;
	Real dG_term_for_residue_at_position(double & position, std::string & residue);  //one component of the overall equation, takes a helix position (can be non-integer
	//if window size is not 19) and a residue and returns the dG_term for that residue at that position
	Real normalize_position_in_helix_given_length(core::Size & position, double & helix_length); //this is what calculates the effective position if the window length is not 19 residues
	Real dG_ins_for_window( std::string & window ); // takes X-residue sequence string and calculates dG_ins for that string
	std::string get_window_at_index( Pose & pose, core::Size & wlen, core::Size & start ); // takes index and returns X-residue sequence string at that index
	static std::string mover_name();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	struct mutt // encodes the relevant Degreaser outputs on a per-residue level
	{
		core::Size resnum; //residue number
		std::string aa; //amino acid it's mutated to
		Real ddG_ins; //change in dG_ins that the mutation causes
		Real dscore; //change in score that the mutation causes
		bool operator< (const mutt & rhs) const
		{
			return ddG_ins < rhs.ddG_ins; //to compare mutants by their ddG_ins
		}
	};
	struct tm_region //each region of interest within a protein, includes the mutt structures so that each possible variant in a window is saved
	{
		core::Size index; //where the region begins
		std::string sequence; //the sequence of that region
		core::Real dG_ins; //the dG_ins of that region
		bool operator< (const tm_region & a) const
		{
			return dG_ins < a.dG_ins; //to compare regions by their dG_ins
		}
		core::Real score;
		utility::vector1<SecretionOptimizationMover::mutt> possible_mutants; //each region of interest can have multiple possible mutants
	} ;

	virtual void report_test( std::ostream & out, Pose & pose); // maybe should be changed to just "report", but confirms the inputs given by the user

	utility::vector1<SecretionOptimizationMover::tm_region> find_tm_regions( std::ostream & out, Pose & pose); // rosetta-independent; finds the sequence windows to design/repack

	Real repack_shell() const {
		return repack_shell_;
	} //setters for the DesignAroundApply
	void repack_shell( Real const repack_shell) {
		repack_shell_ = repack_shell;
		repack_shell_squared_ = repack_shell * repack_shell;
	} //setters for the DesignAroundApply

	void DesignAroundApply( core::pose::Pose const & pose, utility::vector1<SecretionOptimizationMover::mutt>& targets) const; // updated (as of Nov 2017) DesignAroundApply that works with symmetry. Baked into this protocol because it's dependent on it, and maybe DesignAround will have been updated before I check this in

	void packeruser( core::pose::Pose & pose) const; // just calls the Packer
	virtual void test_move(core::pose::Pose & pose) override; // for UMoverTest testing (in theory)

private :
	utility::vector1<SecretionOptimizationMover::mutt> resids_; //this is the list of residues that we can design to, specified for each position
	utility::vector1<SecretionOptimizationMover::mutt> mutids_; //this is the output list of residues (mutants)
	std::string window_ = "AAAAAAAAAAAAAAAAAAA"; //19 (or whatever) -residues to be evaluated
	core::Size window_size_ = 19; //length of evaluation segment
	core::pack::task::TaskFactoryCOP task_factory_;
	core::Size max_mutations_ = 3; //maximum number of mutations to try at once, dumbly
	core::Real dG_ins_threshold_ = 3.5; //transmembrane insertion threshold below which to detect. +3.5 covered the bases that we needed initially.
	core::Real dG_tolerance_ = 0.27; //change in dG_ins by a particular variant to accept. we used +2.7/10 to count something as "significant perturbation"
	core::Real score_tolerance_ = 15; //change in rosetta score above which to reject a variant. started with 0 but +15 REU seems to be a happy medium.
	core::Real repack_shell_ = 8.0; //distance around which to repack residues when evaluating variants
	core::Real repack_shell_squared_ = 64.0;
	std::string keep_ = "DEKRQNSTY"; //string of residues to keep (starts at "aas_allowed", can go down to "")
	std::string chain_to_degrease_ = "A"; //new-ish feature, designate which chain in a pose you want to apply to
	bool dump_multis_ = false;
	bool dump_singles_ = true;
	std::string dump_dir_ = ".";
	bool largest_ddg_ = true; //make the largest dG_ins change for the lowest dG_ins in teh sequence
	bool smallest_dscore_ = false; //make the smallest dscore change
	std::unordered_map<std::string, std::unordered_map<std::string, core::Real>> coeffs; //to initialize the coefficient table, maybe can be moved to DB
	protocols::minimization_packing::PackRotamersMoverOP pack_it_; //call the packer
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1<SecretionOptimizationMover::tm_region> found_regions_; //list of all X-residue windows within the sequence
	core::pack::task::PackerTaskOP task_;
};

} // movers
} // protein_interface_design
} // protocols

#endif

