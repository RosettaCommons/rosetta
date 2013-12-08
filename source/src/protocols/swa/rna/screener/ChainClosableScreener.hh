// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/screener/ChainClosableScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_screener_ChainClosableScreener_HH
#define INCLUDED_protocols_swa_rna_screener_ChainClosableScreener_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/rna/screener/ChainClosableScreener.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {
namespace screener {

	class ChainClosableScreener: public utility::pointer::ReferenceCount {

	public:

		//constructor
		ChainClosableScreener( Size const five_prime_chain_break_res, Size const gap_size );

		ChainClosableScreener( Size const five_prime_chain_break_res, Size const three_prime_chain_break_res, Size const gap_size );

		//destructor
		~ChainClosableScreener();

	public:

		bool
		check_screen( pose::Pose const & pose,
									bool const strict = false ) const;

		bool
		check_screen( pose::Pose const & moving_pose,
									pose::Pose const & reference_pose,
									bool const is_prepend,
									bool const strict = false ) const;

		bool
		check_screen(  pose::Pose const & pose,
									 utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
									 core::kinematics::Stub const & moving_res_base_stub,
									 Size const & reference_res ) const;

		bool
		check_screen(  utility::vector1< core::pose::PoseOP > const & pose_data_list,
									 utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
									 core::kinematics::Stub const & moving_res_base_stub,
									 Size const & reference_res ) const;

	private:

		bool
		check_chain_closable( numeric::xyzVector< core::Real > const & xyz_1, numeric::xyzVector< core::Real > const & xyz_2 ) const;

		bool
		check_chain_closable( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose ) const;

		bool
		check_chain_closable( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose, bool const strict ) const;

		bool
		check_chain_closable( core::Size const & reference_res,
													utility::vector1< core::pose::PoseOP > const & pose_data_list,
													utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
													core::kinematics::Stub const & moving_res_base_stub,
													bool const is_prepend ) const;
		bool
		check_chain_closable( core::Size const & reference_res,
													core::pose::Pose const & pose,
													utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, //this one correspond to the moving_base
													core::kinematics::Stub const & moving_res_base_stub,
													bool const is_prepend ) const;

		void
		get_specific_atom_coordinate( std::string const & atom_name,
																	numeric::xyzVector< core::Real > & atom_pos,
																	core::conformation::Residue const & rsd_at_origin,
																	core::kinematics::Stub const & moving_res_base_stub ) const;

		bool
		check_chain_closable_strict( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose ) const;

		void
		get_C4_C3_distance_range( conformation::Residue const & five_prime_rsd,
															conformation::Residue const & three_prime_rsd,
															Distance & C4_C3_dist_min,
															Distance & C4_C3_dist_max ) const;

	private:

		Size const five_prime_chain_break_res_;
		Size const three_prime_chain_break_res_;
		Size const gap_size_;

	};

} //screener
} //rna
} //swa
} //protocols

#endif
