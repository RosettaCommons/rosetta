// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/pose_creation/MakeJunctionsMover.hh
/// @brief
/// @detailed
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_pose_creation_MakeJunctionsMover_HH
#define INCLUDED_protocols_pose_creation_MakeJunctionsMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/pose_creation/MakeJunctionsMover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <stack>
#include <queue>
#include <set>
#include <map>

namespace protocols {
namespace pose_creation {

class MakeJunctionsMover: public protocols::moves::Mover {
public:
	struct Design{
		std::stack<std::string> n_term_stack;  // I decided to use stacks today cause I'm bored and want to play with a rarely used data type.
		std::stack<std::string> c_term_stack;  // & oddly this case seems like a decent use of stacks.
		std::string start_struct;
		std::string name;
		core::pose::PoseOP n_term_pose;
		core::pose::PoseOP c_term_pose;
		Design(std::string name_i,std::stack<std::string> n_term_stack_i,std::stack<std::string> c_term_stack_i, std::string start_struct_i){
			name = name_i;
			n_term_stack = n_term_stack_i;
			c_term_stack = c_term_stack_i;
			start_struct = start_struct_i;
			n_term_pose = nullptr;
			c_term_pose = nullptr;
		}
		void print(){
			for ( std::stack<std::string> dump = n_term_stack; !dump.empty(); dump.pop() ) {
				std::cout << dump.top() << std::endl;
			}
			std::cout << start_struct << std::endl;
			for ( std::stack<std::string> dump = c_term_stack; !dump.empty(); dump.pop() ) {
				std::cout << dump.top() << std::endl;
			}
		}
	};
	std::string get_name() const override;
	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	MakeJunctionsMover();

	Design line_to_design(std::string line);
	std::queue<Design> read_in_designs();
	void parse_attach_description(std::string attach_description,std::string & pdb_location,std::string & chain, core::Size & n_term_trim, core::Size & c_term_trim, core::Size & n_repeats, core::Size &n_term_attach_length, core::Size &c_term_attach_length, std::string & n_term_seq, std::string & c_term_seq);
	void generate_start_pose(core::pose::Pose & pose, core::pose::Pose & background_pose, std::string attach_description);
	bool attach_next_part(core::pose::Pose & pose, std::string attach_termini, std::string attach_description);
	void trim_pose(core::pose::Pose & pose, core::Size n_term_trim,core::Size c_term_trim);
	void assign_seq(core::pose::Pose & pose, std::string n_term_seq, std::string c_term_seq);
	core::pose::Pose get_and_cache_pdb(std::string pdb_location);
	bool check_all_junctions_good(MakeJunctionsMover::Design design);
	bool make_pose_from_design(MakeJunctionsMover::Design design,core::pose::PoseOP & return_pose);
	void apply( core::pose::Pose & pose ) override;
	core::pose::PoseOP get_additional_output()override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

private:
	std::map<std::string,core::pose::Pose> pdb_cache_;
	core::scoring::ScoreFunctionOP sfxn_;
	std::string designs_fn_;
	std::queue<Design> design_q_;
	bool pdb_cache_bool_;
	core::Size pdb_cache_max_size_;
	core::Real junction_rmsd_thresh_;
	std::multiset<std::string> junction_failure_set_;
	std::string chain_;
};

} //pose_creation
} //protocols

#endif
