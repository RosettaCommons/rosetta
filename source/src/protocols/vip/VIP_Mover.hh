// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_vip_VIP_Mover_HH
#define INCLUDED_protocols_vip_VIP_Mover_HH

#include <protocols/vip/VIP_Report.hh>

namespace protocols {
namespace vip {

class VIP_Mover
{

                core::pose::Pose initial_pose;
		core::pose::Pose cavity_pose;
		core::pose::Pose final_pose;
		core::pose::Pose final_unrelaxed_pose;
		utility::vector1<core::conformation::ResidueOP> temp_residues;
		utility::vector1<core::Size> temp_positions;
		utility::vector1<core::Real> temp_energies;
		utility::vector1<core::conformation::ResidueOP> favorable_residues;
		utility::vector1<core::Size> favorable_positions;
		utility::vector1<core::Real> favorable_energies;
		core::Size number_cavities;
		utility::vector1<core::Size> cavity_balls;
		utility::vector1<std::string> favorable_mutations;
		utility::vector1<core::Size> void_neighbors;
		utility::vector1<core::Size> void_mutatables;
		core::Real final_energy;
		core::Real energy_to_beat;

		utility::vector1<core::Size> excluded_positions;
		core::Size iteration_;
		bool use_stored_best_energy;

	public:
		VIP_Mover();
		VIP_Mover(
                        core::pose::Pose,
                        core::pose::Pose,
                        core::pose::Pose,
                        core::pose::Pose,
                        utility::vector1<core::conformation::ResidueOP>,
                        utility::vector1<core::Size>,
                        utility::vector1<core::Real>,
                        utility::vector1<core::conformation::ResidueOP>,
                        utility::vector1<core::Size>,
                        utility::vector1<core::Real>,
                        core::Size,
                        utility::vector1<core::Size>,
                        utility::vector1<std::string>,
                        utility::vector1<core::Size>,
                        utility::vector1<core::Size>,
			core::Real);
		virtual ~VIP_Mover();

		void set_iteration( core::Size it ) { iteration_ = it; }
		core::Size iteration() { return iteration_; }
		void set_initial_pose( core::pose::Pose );
		core::pose::Pose & get_unrelaxed_pose() { return final_unrelaxed_pose; };
		void minimize_conformation();
		void compute_number_cavities();
		void get_cavity_positions();
		void apply_holes();
		void dump_pdb_to_file( core::pose::Pose &, std::string );
		void get_neighbors();
		void try_point_mutants();
		void relax_favorable_poses();
		void cull_mutatable_residues();
		void sort_fill_energies();
		core::Real get_cav_approx( core::Size );
		// Undefined, commenting out to fix PyRosetta build  void print_favorable_mutations();
		void skip_relax();
		void sort_relaxed_poses();
		void print_pack_report();
		void print_relax_report();
		void nook_finder();
		void cranny_packer();

		void set_excluded_positions();

		core::Real get_final_energy(){
			return final_energy;}
		void set_energy_to_beat( core::Real in_value ){ energy_to_beat = in_value; }
		void set_use_stored_energy( bool in_value ){ use_stored_best_energy = in_value; }
		core::pose::Pose get_final_pose(){
			return final_pose;}
		void apply();


};

bool are_seqs_different( core::pose::Pose & p1, core::pose::Pose & p2 );

}}
#endif
