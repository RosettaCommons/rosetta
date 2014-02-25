// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.hh>
#include <protocols/rotamer_sampler/rigid_body/FloatingBaseUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.rna.checker.ChainClosableGeometryChecker" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace checker {

	//Constructor
	ChainClosableGeometryChecker::ChainClosableGeometryChecker( Size const five_prime_chain_break_res, Size const gap_size ):
		five_prime_chain_break_res_( five_prime_chain_break_res ),
		three_prime_chain_break_res_( five_prime_chain_break_res + 1 ),
		gap_size_( gap_size )
	{
		initialize_distance_range();
	}

	ChainClosableGeometryChecker::ChainClosableGeometryChecker( Size const five_prime_chain_break_res,
																								Size const three_prime_chain_break_res,
																								Size const gap_size ):
		five_prime_chain_break_res_( five_prime_chain_break_res ),
		three_prime_chain_break_res_( three_prime_chain_break_res ),
		gap_size_( gap_size )
	{
		initialize_distance_range();
	}


	//Destructor
	ChainClosableGeometryChecker::~ChainClosableGeometryChecker()
	{}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	void
	ChainClosableGeometryChecker::initialize_distance_range(){
		Distance min_dist_( 0.0 ), max_dist_( 0.0 );
		get_possible_O3prime_C5prime_distance_range( gap_size_, min_dist_, max_dist_ );
		min_dist_squared_ = min_dist_ * min_dist_;
		max_dist_squared_ = max_dist_ * max_dist_;
		dist_squared_ = 0.0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_screen( pose::Pose const & pose,
																			 bool const strict /* = false */ ) const {

		return check_screen( pose, pose, true /* is_prepend -- does not matter */, strict );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_screen( pose::Pose const & moving_pose,
																			 pose::Pose const & reference_pose,
																			 bool const is_prepend,
																			 bool const strict /* = false */ ) const {
		if ( is_prepend ) return check_chain_closable_geometry( moving_pose, reference_pose, strict );
		return check_chain_closable_geometry( reference_pose, moving_pose, strict );
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_screen(  utility::vector1< core::pose::PoseOP > const & pose_data_list,
																				utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
																				core::kinematics::Stub const & moving_res_base_stub,
																				Size const & reference_res ) const {

		for ( Size n = 1; n <= pose_data_list.size(); n++ ){
			pose::Pose const & pose = ( *pose_data_list[n] );
			if ( check_screen( pose, rsd_at_origin_list, moving_res_base_stub, reference_res ) ) return true;
		}
		return false;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_screen(  pose::Pose const & pose,
																				utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
																				core::kinematics::Stub const & moving_res_base_stub,
																				Size const & reference_res ) const {
		//		std::cout << reference_res  << " " << five_prime_chain_break_res_ << std::endl;
		runtime_assert( reference_res >= five_prime_chain_break_res_ );
		bool const is_prepend = ( reference_res > five_prime_chain_break_res_ );
		return check_chain_closable_geometry( reference_res, pose, rsd_at_origin_list, moving_res_base_stub, is_prepend );

	}

	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry( numeric::xyzVector< core::Real > const & xyz_1,
																							numeric::xyzVector< core::Real > const & xyz_2 ) const {
		//Two possibilities
		//1. xyz_1-> five_prime_O3_xyz &&  xyz_2->three_prime_C5_xyz
		//2. xyz_2-> five_prime_O3_xyz &&  xyz_1->three_prime_C5_xyz
		Vector deviation = ( xyz_1 - xyz_2 );
		dist_squared_ = deviation.length_squared();
		return ( ( dist_squared_ > min_dist_squared_ ) &&
						 ( dist_squared_ < max_dist_squared_ ) );
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////
	// following are private functions
	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose, bool const strict ) const {

		if ( strict ) return check_chain_closable_geometry_strict( five_prime_pose, three_prime_pose );
		return check_chain_closable_geometry( five_prime_pose, three_prime_pose );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose ) const {
		return check_chain_closable_geometry( five_prime_pose.residue( five_prime_chain_break_res_  ),
																 three_prime_pose.residue( three_prime_chain_break_res_ ) );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry( core::conformation::Residue const & five_prime_residue,
																							core::conformation::Residue const & three_prime_residue ) const {
		return check_chain_closable_geometry( five_prime_residue.xyz( " O3'" ), three_prime_residue.xyz( " C5'" ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry( core::Size const & reference_res,
																							 utility::vector1< core::pose::PoseOP > const & pose_data_list,
																							 utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list,
																							 core::kinematics::Stub const & moving_res_base_stub,
																							 bool const is_prepend ) const {

		for ( Size n = 1; n <= pose_data_list.size(); n++ ){
			pose::Pose const & pose = ( *pose_data_list[n] );
			if ( check_chain_closable_geometry( reference_res, pose, rsd_at_origin_list, moving_res_base_stub, is_prepend ) ) return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry( core::Size const & reference_res,
																							 core::pose::Pose const & pose,
																							 utility::vector1 < core::conformation::ResidueOP > const & rsd_at_origin_list, //this one correspond to the moving_base
																							 core::kinematics::Stub const & moving_res_base_stub,
																							 bool const is_prepend ) const {
		using namespace core::conformation;

		for ( Size n = 1; n <= rsd_at_origin_list.size(); n++ ){

			Residue const & rsd_at_origin = ( *rsd_at_origin_list[n] );
			std::string const moving_atom_name    = ( is_prepend ) ? " O3'" : " C5'";
			std::string const reference_atom_name = ( is_prepend ) ? " C5'" : " O3'";

			numeric::xyzVector< core::Real > atom_coordinate;
			rotamer_sampler::rigid_body::get_specific_atom_coordinate( moving_atom_name, atom_coordinate, rsd_at_origin, moving_res_base_stub );

			if ( check_chain_closable_geometry( atom_coordinate, pose.residue( reference_res ).xyz( reference_atom_name ) ) ) {
				return true;
			}
		}

		return false;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Optimization for floating_base_chain_closure:
	//
	//   5'-residue                  3'-residue
	//   C3' - O3' - [ P - O5' - ] - C5' - C4' - C3'
	//   |     |                     |     |
	//   |      ---------------------      |
	//   |                                 |
	//    ---------------------------------
	//
	bool
	ChainClosableGeometryChecker::check_chain_closable_geometry_strict( pose::Pose const & five_prime_pose, pose::Pose const & three_prime_pose ) const {

		runtime_assert ( gap_size_ == 0 );

		/////////// C5'-O3' screen //////////////////////////////////////////////////
		Distance C5_O3_dist = ( three_prime_pose.residue( three_prime_chain_break_res_ ).xyz( "C5'" ) - five_prime_pose.residue( five_prime_chain_break_res_ ).xyz( "O3'" ) ).length();
		static Distance const C5_O3_min( 2.866000 ), C5_O3_max( 3.968000 ), leniency_dist( 0.0 );

		//		std::cout << "C5_O3_dist [new]  " << C5_O3_dist << " " << five_prime_chain_break_res_ << " " << three_prime_chain_break_res_ <<
		//			" " << three_prime_pose.residue( three_prime_chain_break_res_ ).xyz( "C5'" )[1] <<
		//			" " << five_prime_pose.residue( five_prime_chain_break_res_ ).xyz( "O3'" )[1] << std::endl;

		//basically cannot close chain if the C5_O3_distance is either too short or too long.
		if ( ( C5_O3_dist > ( C5_O3_max + leniency_dist ) ) || ( C5_O3_dist < ( C5_O3_min - leniency_dist ) ) ) return false;

		/////////// C3'-C4' screen //////////////////////////////////////////////////
		conformation::Residue const & five_prime_rsd  = five_prime_pose.residue( five_prime_chain_break_res_ );
		conformation::Residue const & three_prime_rsd = three_prime_pose.residue( three_prime_chain_break_res_ );
		Distance C4_C3_min( 0.0 ), C4_C3_max( 0.0 );
		get_C4_C3_distance_range( five_prime_rsd, three_prime_rsd, C4_C3_min, C4_C3_max );
		Distance C4_C3_dist = ( three_prime_pose.residue( three_prime_chain_break_res_ ).xyz( " C4'" ) - five_prime_pose.residue( five_prime_chain_break_res_ ).xyz( " C3'" ) ).length();
		if ( ( C4_C3_dist > ( C4_C3_max + leniency_dist ) ) || ( C4_C3_dist < ( C4_C3_min - leniency_dist ) ) ) return false;

		return true;
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Optimization for 'second screen' in floating_base_chain_closure:
	//
	//  Distance depends on vectors pointing into and out of chainbreak
	//
	//   5'-residue                  3'-residue
	//   C3' - O3' - [ P - O5' - ] - C5' - C4' - C3'
	//   | <---                       <---- |
	//    ----------------------------------
	//
	void
	ChainClosableGeometryChecker::get_C4_C3_distance_range( conformation::Residue const & five_prime_rsd,
																									 conformation::Residue const & three_prime_rsd,
																									 Distance & C4_C3_dist_min,
																									 Distance & C4_C3_dist_max ) const{

			numeric::xyzVector< Real > start_vector = five_prime_rsd.xyz( " O3'" ) - five_prime_rsd.xyz( " C3'" );
			numeric::xyzVector< Real > end_vector   = three_prime_rsd.xyz( " C4'" ) - three_prime_rsd.xyz( " C5'" );

			start_vector.normalize();
			end_vector.normalize();

			Real dot_product = dot( start_vector, end_vector );

			// awesome, parin-style.
			if ( dot_product > -1.00 && dot_product < -0.95  ) { C4_C3_dist_min = 2.428;  C4_C3_dist_max = 4.337; }
			else if ( dot_product > -0.95 && dot_product < -0.90  ) { C4_C3_dist_min = 2.238;  C4_C3_dist_max = 4.582; }
			else if ( dot_product > -0.90 && dot_product < -0.85  ) { C4_C3_dist_min = 2.064;  C4_C3_dist_max = 4.743; }
			else if ( dot_product > -0.85 && dot_product < -0.80  ) { C4_C3_dist_min = 1.979;  C4_C3_dist_max = 4.882; }
			else if ( dot_product > -0.80 && dot_product < -0.75  ) { C4_C3_dist_min = 1.833;  C4_C3_dist_max = 4.995; }
			else if ( dot_product > -0.75 && dot_product < -0.70  ) { C4_C3_dist_min = 1.735;  C4_C3_dist_max = 5.099; }
			else if ( dot_product > -0.70 && dot_product < -0.65  ) { C4_C3_dist_min = 1.659;  C4_C3_dist_max = 5.195; }
			else if ( dot_product > -0.65 && dot_product < -0.60  ) { C4_C3_dist_min = 1.590;  C4_C3_dist_max = 5.273; }
			else if ( dot_product > -0.60 && dot_product < -0.55  ) { C4_C3_dist_min = 1.500;  C4_C3_dist_max = 5.347; }
			else if ( dot_product > -0.55 && dot_product < -0.50  ) { C4_C3_dist_min = 1.418;  C4_C3_dist_max = 5.417; }
			else if ( dot_product > -0.50 && dot_product < -0.45  ) { C4_C3_dist_min = 1.337;  C4_C3_dist_max = 5.488; }
			else if ( dot_product > -0.45 && dot_product < -0.40  ) { C4_C3_dist_min = 1.282;  C4_C3_dist_max = 5.552; }
			else if ( dot_product > -0.40 && dot_product < -0.35  ) { C4_C3_dist_min = 1.223;  C4_C3_dist_max = 5.611; }
			else if ( dot_product > -0.35 && dot_product < -0.30  ) { C4_C3_dist_min = 1.145;  C4_C3_dist_max = 5.659; }
			else if ( dot_product > -0.30 && dot_product < -0.25  ) { C4_C3_dist_min = 1.075;  C4_C3_dist_max = 5.713; }
			else if ( dot_product > -0.25 && dot_product < -0.20  ) { C4_C3_dist_min = 1.022;  C4_C3_dist_max = 5.769; }
			else if ( dot_product > -0.20 && dot_product < -0.15  ) { C4_C3_dist_min = 0.963;  C4_C3_dist_max = 5.812; }
			else if ( dot_product > -0.15 && dot_product < -0.10  ) { C4_C3_dist_min = 1.019;  C4_C3_dist_max = 5.861; }
			else if ( dot_product > -0.10 && dot_product < -0.05  ) { C4_C3_dist_min = 1.331;  C4_C3_dist_max = 5.904; }
			else if ( dot_product > -0.05 && dot_product < 0.00  ) { C4_C3_dist_min = 1.532;  C4_C3_dist_max = 5.942; }
			else if ( dot_product > 0.00 && dot_product < 0.05  ) { C4_C3_dist_min = 1.768;  C4_C3_dist_max = 5.979; }
			else if ( dot_product > 0.05 && dot_product < 0.10  ) { C4_C3_dist_min = 1.953;  C4_C3_dist_max = 6.017; }
			else if ( dot_product > 0.10 && dot_product < 0.15  ) { C4_C3_dist_min = 2.121;  C4_C3_dist_max = 6.046; }
			else if ( dot_product > 0.15 && dot_product < 0.20  ) { C4_C3_dist_min = 2.292;  C4_C3_dist_max = 6.083; }
			else if ( dot_product > 0.20 && dot_product < 0.25  ) { C4_C3_dist_min = 2.424;  C4_C3_dist_max = 6.118; }
			else if ( dot_product > 0.25 && dot_product < 0.30  ) { C4_C3_dist_min = 2.563;  C4_C3_dist_max = 6.140; }
			else if ( dot_product > 0.30 && dot_product < 0.35  ) { C4_C3_dist_min = 2.726;  C4_C3_dist_max = 6.171; }
			else if ( dot_product > 0.35 && dot_product < 0.40  ) { C4_C3_dist_min = 2.849;  C4_C3_dist_max = 6.200; }
			else if ( dot_product > 0.40 && dot_product < 0.45  ) { C4_C3_dist_min = 2.998;  C4_C3_dist_max = 6.219; }
			else if ( dot_product > 0.45 && dot_product < 0.50  ) { C4_C3_dist_min = 3.128;  C4_C3_dist_max = 6.245; }
			else if ( dot_product > 0.50 && dot_product < 0.55  ) { C4_C3_dist_min = 3.261;  C4_C3_dist_max = 6.261; }
			else if ( dot_product > 0.55 && dot_product < 0.60  ) { C4_C3_dist_min = 3.380;  C4_C3_dist_max = 6.284; }
			else if ( dot_product > 0.60 && dot_product < 0.65  ) { C4_C3_dist_min = 3.523;  C4_C3_dist_max = 6.298; }
			else if ( dot_product > 0.65 && dot_product < 0.70  ) { C4_C3_dist_min = 3.658;  C4_C3_dist_max = 6.315; }
			else if ( dot_product > 0.70 && dot_product < 0.75  ) { C4_C3_dist_min = 3.785;  C4_C3_dist_max = 6.329; }
			else if ( dot_product > 0.75 && dot_product < 0.80  ) { C4_C3_dist_min = 3.914;  C4_C3_dist_max = 6.340; }
			else if ( dot_product > 0.80 && dot_product < 0.85  ) { C4_C3_dist_min = 4.065;  C4_C3_dist_max = 6.350; }
			else if ( dot_product > 0.85 && dot_product < 0.90  ) { C4_C3_dist_min = 4.209;  C4_C3_dist_max = 6.356; }
			else if ( dot_product > 0.90 && dot_product < 0.95  ) { C4_C3_dist_min = 4.374;  C4_C3_dist_max = 6.357; }
			else if ( dot_product > 0.95 && dot_product < 1.00  ) { C4_C3_dist_min = 4.570;  C4_C3_dist_max = 6.349; }
			else{
				TR << "dot_product = " << dot_product << std::endl;
				utility_exit_with_message( "Invalid dot_product!" );
			}
	}



} //checker
} //rna
} //sampling
} //stepwise
} //protocols
