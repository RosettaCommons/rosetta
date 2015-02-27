// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_helix_capper_HelixNCapperMover_HH
#define INCLUDED_protocols_helix_capper_HelixNCapperMover_HH

#include <core/pose/Pose.hh>
#include <utility/libsvm/Svm_rosetta.hh>
#include <utility/libsvm/Svm_rosetta.fwd.hh>

namespace protocols {
namespace helix_capper {

class HelixNCapperMover
{

	public:
		HelixNCapperMover();
		HelixNCapperMover(
			core::pose::Pose start_pose
		);
		virtual ~HelixNCapperMover();

		//void set_initial_pose( core::pose::Pose );
		void dump_pdb_to_file( core::pose::Pose &, std::string );
		// Undefined, commenting out to fix PyRosetta build  void print_favorable_mutations();

		void set_excluded_positions();

		void setup_svm();
		void get_start_positions();
		void get_Ncap_scores();
		core::Real ncap_prob_from_svm( utility::vector1< core::Real > & );

		void apply();

	private:
		core::pose::Pose start_pose_;
		utility::libsvm::Svm_rosettaOP ncap_model_;
		utility::vector1<core::Size> excluded_positions_;
		utility::vector1<core::Size> helix_start_positions_;
		utility::vector1<core::Real> ncap_scores_;

};

}}
#endif
