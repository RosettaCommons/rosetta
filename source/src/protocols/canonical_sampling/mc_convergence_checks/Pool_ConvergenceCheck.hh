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
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_Pool_ConvergenceCheck_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_Pool_ConvergenceCheck_hh


// type headers
#include <core/types.hh>

// unit headers
#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// utility headers

#include <utility/pointer/ReferenceCount.hh>
// #include "utility/basic_sys_util.h"

// C++ headers
#include <map>
#include <string>

#include <utility/vector1.hh>


// Forward declarations

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {


class EXCN_Pool_Converged : public moves::EXCN_Converged
{};

class Pool_RMSD;

typedef utility::pointer::shared_ptr< Pool_RMSD > Pool_RMSD_OP;


class Pool_RMSD : public utility::pointer::ReferenceCount {
	typedef utility::vector1< std::string > Tags;
public:
	//c'stor supply file with structures and threshold (Angstroem RMSD)
	Pool_RMSD( std::string silent_file )
	{ fill_pool( silent_file); };

	Pool_RMSD() {};

	virtual ~Pool_RMSD();

	/// @brief return position in pool for the best_decoy
	core::Size evaluate( core::pose::Pose const&, std::string& best_decoy, core::Real& best_rmsd ) const;

	core::Size evaluate( core::io::silent::SilentStruct const&, std::string& best_decoy, core::Real& best_rmsd ) const;
	virtual core::Size evaluate_and_add( core::pose::Pose const&  pose, std::string& cluster_center, core::Real& best_rmsd, core::Real transition_threshold );

	core::Size size() const { return pool_.n_decoys(); };


	void add( core::io::silent::SilentStruct const&, std::string new_tag = "keep_silent_tag" );
	void add( core::pose::Pose const&, std::string new_tag );

	virtual void finalize() { return; };

	void get(core::Size index, ObjexxFCL::FArray2D_double& result);
	std::string& get_tag(core::Size index);

	void pop_back();
	void clear();

	void set_reserve_size( int max_size ); //ek


	core::Size evaluate_and_add( ObjexxFCL::FArray2D_double& coords, std::string& best_decoy, core::Real& best_rmsd, core::Real transition_threshold  );
	core::Size evaluate( ObjexxFCL::FArray2D_double&, std::string& best_decoy, core::Real& best_rmsd ) const; //ek sorry had to move it from  private to protected
	core::Size evaluate( ObjexxFCL::FArray2D_double&, std::string& best_decoy, core::Real& best_rmsd, core::Size index) const;
	void add( ObjexxFCL::FArray2D_double const&, core::Size nres, std::string new_tag);

	//void add( ObjexxFCL::FArray2D_double const&, core::Size nres, core::Size num_added, utility::vector1<std::string> const& new_tag);
	void set_excluded_residues( utility::vector1< core::Size > const& excluded_residues );
protected:

	core::Size natom() const {return pool_.n_atoms();} //yuan
	ObjexxFCL::FArray3D_double const& coords() const {return pool_.coords();} //yuan
	std::string const &tag(core::Size index) const {return tags_[index];}


private:
	void fill_pool( std::string const& silentfile );

	toolbox::DecoySetEvaluation pool_;
	ObjexxFCL::FArray1D_double weights_;
	Tags tags_;
	core::Size reserve_size_;

	//residues excluded from RMSD calculation
	utility::vector1< core::Size > excluded_residues_;
};

class Pool_ConvergenceCheck : public moves::MonteCarloExceptionConverge {
public:
	Pool_ConvergenceCheck( Pool_RMSD_OP rmsd_pool_in, core::Real threshold )
	: threshold_( threshold ),
		rmsd_pool_( rmsd_pool_in ) {};

	virtual ~Pool_ConvergenceCheck();
	//throws EXCN_Pool_Converged if lowest_score pose is < threshold away from any pool structure
	virtual bool operator() ( const core::pose::Pose&, moves::MonteCarlo const& mc, bool /*reject*/ );

private:
	core::Real threshold_;
	Pool_RMSD_OP rmsd_pool_;
};

class Pool_Evaluator : public evaluation::PoseEvaluator {
public:
	virtual ~Pool_Evaluator();
	Pool_Evaluator( Pool_RMSD_OP rmsd_pool_in ) : rmsd_pool_( rmsd_pool_in ) {};
	/// from PoseEvaluator interface:

	/// @brief evaluate pose and store values in Silent_Struct
	virtual void apply( core::pose::Pose&, std::string tag, core::io::silent::SilentStruct &pss) const;
	virtual void apply( core::io::silent::SilentStruct &pss) const;

	// void apply( core::pose::Pose&, std::string tag, core::io::silent::SilentStruct &pss) const;
	virtual core::Size size() const { return 2; };
	virtual std::string name( core::Size i ) const {
		if ( i == 1 ) return "pool_converged_tag";
		if ( i == 2 ) return "pool_converged_rmsd";
		return "NO_TAG";
	}
private:
	Pool_RMSD_OP rmsd_pool_;
};

} // mc_convergence_check
} // moves
} // rosetta

#endif
