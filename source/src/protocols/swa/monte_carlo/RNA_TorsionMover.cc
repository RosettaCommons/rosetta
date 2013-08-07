// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeleteMover
/// @brief Torsions an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/swa/monte_carlo/RNA_TorsionMover.hh>
#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloUtil.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/swa/monte_carlo/types.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/rna/RNA_Util.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <numeric/random/random.hh>

#include <basic/Tracer.hh>

#include <map>


static numeric::random::RandomGenerator RG(245099111);  // <- Magic number, do not change it!

using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose full_model_info object.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_torsion_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	//	RNA_TorsionMover::RNA_TorsionMover():
	//		Mover(),
	//		default_sample_range_( 10.0 )
	//	{}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_TorsionMover::~RNA_TorsionMover()
  {}

	//////////////////////////////////////////////////////////////////////
  void
  RNA_TorsionMover::apply( core::pose::Pose & pose )
	{
		std::string move_type = "";
		default_sample_range_ = 10; //not very elegant.
		apply( pose, move_type, default_sample_range_ );
	}

	//////////////////////////////////////////////////////////////////////
  void
  RNA_TorsionMover::apply( core::pose::Pose & pose, std::string & move_type, Real const & sample_range )
	{

		using namespace core::pose::full_model_info;

		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		utility::vector1< Size > const & moving_res_list = full_model_info.moving_res_list();
		random_torsion_move( pose, moving_res_list, move_type, sample_range );

	}

	//////////////////////////////////////////////////////////////////////
	void
	RNA_TorsionMover::random_torsion_move( pose::Pose & pose,
																				 utility::vector1< Size > const & moving_res_list,
																				 std::string & move_type,
																				 Real const & sample_range ){

		Size const random_idx = int( RG.uniform() * moving_res_list.size() ) + 1;
		Size const i = moving_res_list[ random_idx ];

		MovingResidueCase moving_residue_case = get_moving_residue_case( pose, i );

		if ( moving_residue_case == CHAIN_TERMINUS_3PRIME  || moving_residue_case == CHAIN_TERMINUS_5PRIME ){

			// an edge residue -- change both its nucleoside & suite -- can go crazy.
			Size const nucleoside_num = i;
			sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range);

			Size suite_num( 0 );
			if ( moving_residue_case == CHAIN_TERMINUS_3PRIME ) suite_num = i-1;
			else suite_num = i;

			if ( RG.uniform() < 0.5) {
				sample_near_suite_torsion( pose, suite_num, sample_range);
				move_type += "-nuc-suite";
			} else {
				crankshaft_alpha_gamma( pose, suite_num, sample_range);
				move_type += "-nuc-crank";
			}

		} else {
			runtime_assert( moving_residue_case == INTERNAL ); // cannot handle floating base yet.

			// don't do anything super-crazy -- either do previous suite, current nucleoside, or next suite.
			Real const random_number = RG.uniform();

			if ( random_number < 0.6 ){

				Size suite_num( 0 );
				if ( RG.uniform() < 0.5) {
					suite_num= i-1;
				} else {
					suite_num= i;
				}

				if ( RG.uniform() < 0.5) {
					sample_near_suite_torsion( pose, suite_num, sample_range);
					//move_type += "-suite" + string_of(suite_num);
					move_type += "-suite";
				} else {
					crankshaft_alpha_gamma( pose, suite_num, sample_range);
					//move_type += "-suite" + string_of(suite_num);
					move_type += "-crank";
				}
			} else {
				Size const nucleoside_num = i;
				sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range);
				//move_type += "-nuc" + string_of(nucleoside_num);
				move_type += "-nuc";
			}

		}

	}


//////////////////////////////////////////////////
	std::string
	RNA_TorsionMover::get_name() const {
		return "RNA_TorsionMover";
	}


//////////////////////////////////////////////////
void
RNA_TorsionMover::sample_near_suite_torsion(utility::vector1< Real > & torsion_list, Real const stddev) {
	//	static const Real delta_north = rna_fitted_torsion_info_.ideal_delta_north();
	//	static const Real delta_south = rna_fitted_torsion_info_.ideal_delta_south();

	torsion_list[1] += RG.gaussian() * stddev;
	torsion_list[2] += RG.gaussian() * stddev;
	torsion_list[3] += RG.gaussian() * stddev;
	torsion_list[4] += RG.gaussian() * stddev;
	torsion_list[5] += RG.gaussian() * stddev;

}


//////////////////////////////////////////////////
// copied from fang's turner_test_one_chain.cc
//  why not use principal_value  [-180 to 180]?
void
RNA_TorsionMover::sample_near_nucleoside_torsion(utility::vector1< Real > & torsion_list, Real const stddev) {
	static const Real delta_north = rna_fitted_torsion_info_.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info_.ideal_delta_south();

	if (RG.uniform() < 0.2) {
		torsion_list[1]  = (RG.uniform() < 0.5) ? delta_south : delta_north;
	}

	torsion_list[2] += RG.gaussian() * stddev;
	if (torsion_list[2] > 360) {
		torsion_list[2] -= 360;
	} else if (torsion_list[2] <=  0) {
		torsion_list[2] += 360;
	}

}

//////////////////////////////////
void
RNA_TorsionMover::apply_nucleoside_torsion( utility::vector1< Real > const & torsion_set,
													pose::Pose & pose,
													Size const moving_res){

	using namespace id;
	using namespace chemical::rna;

	Real delta, nu2, nu1;
	if (torsion_set[1] < 115) { //North pucker, [6] is delta angle (only pick one of the two states)
		delta = rna_fitted_torsion_info_.ideal_delta_north();
		nu2 = rna_fitted_torsion_info_.ideal_nu2_north();
		nu1 = rna_fitted_torsion_info_.ideal_nu1_north();
	} else { //South pucker
		delta = rna_fitted_torsion_info_.ideal_delta_south();
		nu2 = rna_fitted_torsion_info_.ideal_nu2_south();
		nu1 = rna_fitted_torsion_info_.ideal_nu1_south();
	}

	pose.set_torsion( TorsionID( moving_res, id::BB,  4 ), delta );
	pose.set_torsion( TorsionID( moving_res, id::CHI, 2 ), nu2 );
	pose.set_torsion( TorsionID( moving_res, id::CHI, 3 ), nu1 );
	pose.set_torsion( TorsionID( moving_res, id::CHI, 1 ), torsion_set[2] );
}

//////////////////////////////////
void
RNA_TorsionMover::apply_suite_torsion( utility::vector1< Real > const & torsion_set,
										 pose::Pose & pose,
										 Size const moving_suite ){

	using namespace id;
	using namespace chemical::rna;

	pose.set_torsion( TorsionID( moving_suite, id::BB, 5 ), torsion_set[1] );   //epsilon
	pose.set_torsion( TorsionID( moving_suite, id::BB, 6 ), torsion_set[2] );   //zeta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 1 ), torsion_set[3] ); //alpha
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 2 ), torsion_set[4] ); //beta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 3 ), torsion_set[5] ); //gamma

}


//////////////////////////////////
void
RNA_TorsionMover::apply_random_nucleoside_torsion( pose::Pose & pose,
																									 Size const moving_res ){

	utility::vector1< Real > torsion_set;

	bool north_pucker = (RG.uniform() < 0.5) ? true : false;

	Size chi_rotamer = 1;
	// could be syn if purine.
	if ( chemical::rna::is_purine( pose.residue( moving_res ) ) && RG.uniform() < 0.5 ) chi_rotamer = 2;

	if ( north_pucker ){
		torsion_set.push_back( rna_fitted_torsion_info_.ideal_delta_north() );
		torsion_set.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_chi_north()[chi_rotamer].center );
	} else {
		torsion_set.push_back( rna_fitted_torsion_info_.ideal_delta_south() );
		torsion_set.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_chi_south()[chi_rotamer].center );
	}

	apply_nucleoside_torsion( torsion_set, pose, moving_res );
}


//////////////////////////////////
void
RNA_TorsionMover::apply_random_suite_torsion( pose::Pose & pose,
														Size const moving_suite ){

	utility::vector1< Real > torsion_set;

	bool north_pucker = ( pose.delta( moving_suite ) < 115 );
	Real const epsilon =  ( north_pucker ) ? rna_fitted_torsion_info_.gaussian_parameter_set_epsilon_north()[1].center : rna_fitted_torsion_info_.gaussian_parameter_set_epsilon_south()[1].center;
	torsion_set.push_back( epsilon );

	Size const alpha_rotamer = int( 3 * RG.uniform() ) + 1;
	Real const alpha = rna_fitted_torsion_info_.gaussian_parameter_set_alpha()[ alpha_rotamer ].center;

	Size const zeta_rotamer = int( 2 * RG.uniform() ) + 1;
	Real zeta;
	if ( alpha_rotamer == 1 ){
		zeta = rna_fitted_torsion_info_.gaussian_parameter_set_zeta_alpha_sc_minus()[ zeta_rotamer ].center;
	} else if ( alpha_rotamer == 2 ){
		zeta = rna_fitted_torsion_info_.gaussian_parameter_set_zeta_alpha_sc_plus()[ zeta_rotamer ].center;
	} else {
		zeta = rna_fitted_torsion_info_.gaussian_parameter_set_zeta_alpha_ap()[ zeta_rotamer ].center;
	}

	torsion_set.push_back( zeta );
	torsion_set.push_back( alpha );
	torsion_set.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_beta()[1].center );

	Size const gamma_rotamer = int( 3 * RG.uniform() ) + 1;
	torsion_set.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_gamma()[gamma_rotamer].center );

	apply_suite_torsion( torsion_set, pose, moving_suite );

}


//////////////////////////////////
void
RNA_TorsionMover::apply_nucleoside_torsion_Aform(
																								 pose::Pose & pose,
																								 Size const moving_res ){

	utility::vector1< Real > ideal_A_form_torsions;
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.ideal_delta_north() );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_chi_north()[1].center );
	apply_nucleoside_torsion( ideal_A_form_torsions, pose, moving_res );
}


//////////////////////////////////
void
RNA_TorsionMover::apply_suite_torsion_Aform(
													pose::Pose & pose,
													Size const moving_suite ){

	utility::vector1< Real > ideal_A_form_torsions;
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_epsilon_north()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_alpha()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_beta()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info_.gaussian_parameter_set_gamma()[1].center );
	apply_suite_torsion( ideal_A_form_torsions, pose, moving_suite );
}


///////////////////////////////////////////////////
utility::vector1< Real>
RNA_TorsionMover::get_suite_torsion( pose::Pose const & pose, Size const moving_suite ){

	using namespace id;

	utility::vector1< Real > torsion_set;
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite, id::BB, 5 ) ) );   //epsilon
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite, id::BB, 6 ) ) );   //zeta
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite+1, id::BB, 1 ) ) ); //alpha
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite+1, id::BB, 2 ) ) ); //beta
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite+1, id::BB, 3 ) ) ); //gamma

	return torsion_set;
}

///////////////////////////////////////////////////
utility::vector1< Real>
RNA_TorsionMover::get_nucleoside_torsion( pose::Pose const & pose, Size const moving_nucleoside ){

	using namespace id;

	utility::vector1< Real > torsion_set;
	torsion_set.push_back(	pose.torsion( TorsionID( moving_nucleoside, id::BB, 4 ) ) );  //delta
	torsion_set.push_back(	pose.torsion( TorsionID( moving_nucleoside, id::CHI, 1 ) ) ); //chi

	return torsion_set;
}

///////////////////////////////////////////////////
void
RNA_TorsionMover::sample_near_suite_torsion( pose::Pose & pose, Size const moving_suite, Real const sample_range){
	utility::vector1< Real> torsion_set = get_suite_torsion( pose, moving_suite );
	sample_near_suite_torsion( torsion_set, sample_range );
	apply_suite_torsion( torsion_set, pose, moving_suite );
}

///////////////////////////////////////////////////
void
RNA_TorsionMover::sample_near_nucleoside_torsion( pose::Pose & pose, Size const moving_res, Real const sample_range){
	utility::vector1< Real> torsion_set = get_nucleoside_torsion( pose, moving_res );
	sample_near_nucleoside_torsion( torsion_set, sample_range);
	apply_nucleoside_torsion( torsion_set, pose, moving_res );
}
///////////////////////////////////////////////////
void
RNA_TorsionMover::crankshaft_alpha_gamma( pose::Pose & pose, Size const moving_suite, Real const sample_range){

	using namespace id;

	TorsionID alpha_torsion_id( moving_suite+1, id::BB, 1 );
	TorsionID gamma_torsion_id( moving_suite+1, id::BB, 3 );

	Real alpha = pose.torsion( alpha_torsion_id );
	Real gamma = pose.torsion( gamma_torsion_id );
	Real const perturb = RG.gaussian() * sample_range;

	alpha += perturb;
	gamma -= perturb;

	pose.set_torsion( alpha_torsion_id, alpha );
	pose.set_torsion( gamma_torsion_id, gamma );

}


}
}
}
