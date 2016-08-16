// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/mc_convergence_checks/HPool.hh
/// @brief hierarchical pool
/// @author Yuan Liu (wendao@u.washington.edu)

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_HPool_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_HPool_hh

#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2.hh>
#include <core/pose/Pose.hh>

#include <protocols/canonical_sampling/mc_convergence_checks/HPool.fwd.hh>
#include <deque>

#include <utility/vector1.hh>


namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

class HPool_RMSD : public Pool_RMSD
{
public:
	HPool_RMSD(std::string silent_file, core::Size lv=1);

	//eval with a required radius
	core::Size evaluate(
		core::pose::Pose& pose,
		core::Real resolution,
		std::string& best_decoy,
		core::Real& best_rmsd );

	core::Size evaluate(
		core::io::silent::SilentStruct& pss,
		core::Real resolution,
		std::string& best_decoy,
		core::Real& best_rmsd );

	//load the correspond subcluster
	bool load_lib(const core::Size);
	//clear
	//void clear_lib();
	void clear_lib(const core::Size);
	//get_size
	core::Size get_size(core::Size nsubc)
	{
		return (subpools_[nsubc].get())?subpools_[nsubc]->size():0;
	}

	//void add(core::pose::Pose const &, std::string &);
	//void add(core::io::silent::SilentStruct const& pss, std::string &tag);

	void debug();

private:
	void build_pair_dis_matrix();
	core::Real get_pair_dist(core::Size, core::Size) const;
	core::Real dist_square(ObjexxFCL::FArray2_double &, ObjexxFCL::FArray2_double &);
	core::Real dist_square(core::Size, core::Size);

protected:
	//hierarchy evaluate FArray2D
	core::Size evaluate(
		ObjexxFCL::FArray2D_double& coord,
		core::Real resolution,
		std::string& best_decoy,
		core::Real& best_rmsd );

	//fast eval returns the index and tags at
	core::Size evaluate_core(
		ObjexxFCL::FArray2D_double& coord,
		std::string& best_decoy,
		core::Real& best_rmsd,
		core::Size index ) const;

	void add(ObjexxFCL::FArray2D_double &, std::string &);

private:
	utility::vector1< HPool_RMSD_OP > subpools_;
	utility::vector1< core::Real > pair_dis_;
	std::string silent_file_;
	std::string tag_prefix_;
	core::Size old_size_;
	core::Size level_;
	core::Size n_extra_;
	core::Real radius_;
	bool has_child_;

	//deque
	std::deque< core::Size > sub_ndx_deque_;
};

std::string lib_full_path(std::string tag);

}//mc_conv
}//moves
}//prot

#endif

