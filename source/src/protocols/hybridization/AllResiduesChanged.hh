// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_hybridization_AllResiduesChanged_hh
#define INCLUDED_protocols_hybridization_AllResiduesChanged_hh

#include <protocols/moves/WhileMover.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <utility/vector1.hh>
#include <set>

namespace protocols {
namespace hybridization {

/// @brief (helper) functor class which keeps track of initial phi/psi values.
/// @detail
/// calls of operator ( pose ) compare the initial phi/psi values
////to the current values of the given pose. Returns false once all phi/psi values
/// have been modified.
class AllResiduesChanged : public moves::PoseCondition {
public:
	AllResiduesChanged( core::pose::Pose const & pose,
		utility::vector1< core::Real > const residue_weights,
		utility::vector1< core::Size > const anchor_residues=utility::vector1< core::Size >(0)
	) :
		insert_pos_( pose.total_residue(), false )
	{
		set_initial_pose( pose );
		for ( core::Size i_anchor = 1; i_anchor <= anchor_residues.size(); ++i_anchor ) {
			anchor_residues_.insert(anchor_residues[i_anchor]);
		}
		compute_insert_pos( residue_weights );
	}

	AllResiduesChanged( core::pose::Pose const & pose ) :
		insert_pos_( pose.total_residue(), true )
	{
		set_initial_pose( pose );
	}

private:

	void compute_insert_pos(
		utility::vector1< core::Real > const residue_weights
	) {
		for ( core::Size ires=1; ires<= residue_weights.size(); ++ires ) {
			if ( ires > insert_pos_.size() ) break;
			if ( residue_weights[ires] > 0.0 && !anchor_residues_.count(ires) ) {
				insert_pos_[ires] = true;
			}
		}
	}

	void set_initial_pose( const core::pose::Pose & pose ) {
		initial_phis.resize(pose.total_residue());
		initial_psis.resize(pose.total_residue());
		for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
			if ( ! pose.residue(i).is_protein() ) continue;
			initial_phis[i] = pose.phi(i);
			initial_psis[i] = pose.psi(i);
		}

		original_sequence_ = pose.sequence();
	}

public:
	void show_unmoved( const core::pose::Pose & pose, std::ostream& out ) {
		runtime_assert( original_sequence_ == pose.sequence() );
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( ! pose.residue(i).is_protein() )  continue;
			if ( initial_phis[i] == pose.phi(i) && insert_pos_[ i ] ) {
				out << i << " ";
				continue;
			}
			if ( initial_psis[i] == pose.psi(i) && insert_pos_[ i ] ) {
				out << i << " ";
			}
		}
		out << std::endl;
	}

	virtual bool operator() ( const core::pose::Pose & pose ) {
		runtime_assert( original_sequence_ == pose.sequence() ); // imperfect attempt to check that Pose hasn't changed ...
		for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
			if ( ! pose.residue(i).is_protein() ) continue;
			if ( initial_phis[i] == pose.phi(i) && insert_pos_[ i ] ) {
				return false;
			}
			if ( initial_psis[i] == pose.psi(i) && insert_pos_[ i ] ) {
				return false;
			}
		}
		return true;
	}

private:
	utility::vector1< core::Real > initial_phis;
	utility::vector1< core::Real > initial_psis;

	std::set<core::Size> anchor_residues_;
	std::string original_sequence_;
	utility::vector1< bool > insert_pos_;
};

}
}

#endif
