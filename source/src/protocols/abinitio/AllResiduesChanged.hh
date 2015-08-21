// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_abinitio_AllResiduesChanged_hh
#define INCLUDED_protocols_abinitio_AllResiduesChanged_hh

#include <protocols/moves/WhileMover.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

//so far this class is only included by FragmentSampler.cc and ClassicAbinitio.cc
//change to a proper modul if you must


/// @brief (helper) functor class which keeps track of initial phi/psi values.
/// @detail
/// calls of operator ( pose ) compare the initial phi/psi values
////to the current values of the given pose. Returns false once all phi/psi values
/// have been modified.
class AllResiduesChanged : public moves::PoseCondition {
public:
	AllResiduesChanged( core::pose::Pose const & pose,
		core::fragment::InsertMap const& insert_map,
		core::kinematics::MoveMap const& mm
	) :
		insert_pos_( pose.total_residue(), false )
	{
		set_initial_pose( pose );
		compute_insert_pos( insert_map, mm );
	}

	AllResiduesChanged( core::pose::Pose const & pose ) :
		insert_pos_( pose.total_residue(), true )
	{
		set_initial_pose( pose );
	}

private:

	void compute_insert_pos( core::fragment::InsertMap const& insert_map,
		core::kinematics::MoveMap const& mm
	) {
		for ( core::fragment::InsertMap::const_iterator it = insert_map.begin(),
				eit = insert_map.end(); it != eit; ++it ) {
			Size const pos ( *it );
			if ( pos > insert_pos_.size() ) break;
			if ( mm.get_bb( pos ) ) {
				insert_pos_[ pos ] = true;
			}
		}
	}

	void set_initial_pose( const core::pose::Pose & pose ) {
		for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
			if ( ! pose.residue(i).is_protein() ) continue;
			initial_phis.push_back( pose.phi(i) );
			initial_psis.push_back( pose.psi(i) );
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

	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool initialized_;

	std::string original_sequence_;
	utility::vector1< bool > insert_pos_;
};

}
}

#endif
