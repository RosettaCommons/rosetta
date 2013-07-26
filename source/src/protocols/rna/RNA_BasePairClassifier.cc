// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief Implementation of Leontis/Westhof nucleic acid base-pair classification
/// @detailed
/// @author Rhiju Das


// Unit headers
#include <protocols/rna/RNA_BasePairClassifier.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyzMatrix.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

// External library headers

//C++ headers
#include <vector>
#include <string>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


using namespace core;
using basic::T;

static basic::Tracer TR( "protocols.rna.rna_base_pair_classifier" ) ;

namespace protocols {
namespace rna {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Implementation of Leontis/Westhof  base pair classifier (RNA, 2001)
//
// Partly follows RNAVIEW algorithm, though no use of CH..O bonds, and
// all base pairs get a Watson-Crick,Hoogsteen,Sugar classification.
//
// This should probably be made into an object.
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////
void
update_edge_hbond_numbers(
   conformation::Residue const & rsd,
	 Size const & atm,
	 Size & N_W,
	 Size & N_H,
	 Size & N_S )
{
	using namespace core::chemical;

	//std::cout << atm << std::endl;
	std::string atom_name = rsd.atom_name( atm );
	//	std::cout << atom_name << std::endl;

	if ( rsd.aa() == na_rad ) {

		if ( atom_name == " N1 "  ||
				 atom_name == " C2 "  ||
				 atom_name == " N6 "  ) N_W++;

		if ( atom_name == " N6 "  ||
				 atom_name == " C5 "  ||
				 atom_name == " C8 "  ||
				 atom_name == " N7 "  ) N_H++;

		if ( atom_name == " N3 "  ||
				 atom_name == " C2 "  ||
				 atom_name == " C4 "  ||
				 atom_name == " C1'"  ||
				 atom_name == " C3'"  ||
				 atom_name == " O3'"  ||
				 atom_name == " O2'"   ) N_S++;

	} else if ( rsd.aa() == na_rcy ) {

		if ( atom_name == " O2 "  ||
				 atom_name == " N3 "  ||
				 atom_name == " N4 " ) N_W++;

		if ( atom_name == " N4 "  ||
				 atom_name == " C5 "  ||
				 atom_name == " C6 "  ) N_H++;

		if ( atom_name == " O2 "  ||
				 atom_name == " N1 "  ||
				 atom_name == " C1'"  ||
				 atom_name == " C3'"  ||
				 atom_name == " O3'"  ||
				 atom_name == " O2'"   ) N_S++;

	} else if ( rsd.aa() == na_rgu ) {

		if ( atom_name == " N1 "  ||
				 atom_name == " N2 "  ||
				 atom_name == " O6 "  ) N_W++;

		if ( atom_name == " O6 "  ||
				 atom_name == " C5 "  ||
				 atom_name == " C8 "  ||
				 atom_name == " N7 "  ) N_H++;

		if ( atom_name == " N3 "  ||
				 atom_name == " N2 "  ||
				 atom_name == " C4 "  ||
				 atom_name == " N9 "  ||
				 atom_name == " C1'"  ||
				 atom_name == " O2'"   ) N_S++;

	} else if ( rsd.aa() == na_ura ) {

		if ( atom_name == " O2 "  ||
				 atom_name == " N3 "  ||
				 atom_name == " O4 " ) N_W++;

		if ( atom_name == " O4 "  ||
				 atom_name == " C5 "  ||
				 atom_name == " C6 "  ) N_H++;

		if ( atom_name == " O2 "  ||
				 atom_name == " N1 "  ||
				 atom_name == " C1'"  ||
				 atom_name == " C3'"  ||
				 atom_name == " O3'"  ||
				 atom_name == " O2'"   ) N_S++;

	} else {
		std::cout << "PROBLEM !!!! " << rsd.aa() << std::endl;
		utility_exit_with_message( "Problem with base classification, residue " );
	}
}

///////////////////////////////////
void
update_edge_hbond_numbers_careful_hydrogen(
   conformation::Residue const & rsd,
	 Size const & atm,
   conformation::Residue const & other_rsd,
	 Size const & other_atm,
	 Size & N_W,
	 Size & N_H,
	 Size & N_S )
{
	using namespace core::chemical;

	std::string atom_name = rsd.atom_name( atm );

	if (rsd.aa() == na_rad && atom_name == " N6 ") {
		//std::cout << "CHECKING " << rsd.seqpos() << std::endl;
		if ( (rsd.xyz( rsd.atom_index(" H61") ) - other_rsd.xyz( other_atm )).length()  <
				 (rsd.xyz( rsd.atom_index(" H62") ) - other_rsd.xyz( other_atm )).length() ) {
			N_W++;
		} else {
			N_H++;
		}
	}

	if (rsd.aa() == na_rcy && atom_name == " N4 ") {
		//		TR << "cyt check " << rsd.seqpos() << " to " << other_rsd.seqpos() << " atom " << other_rsd.atom_name( other_atm ) <<
		//			"  dist to H42 " << (rsd.xyz( rsd.atom_index(" H42") ) - other_rsd.xyz( other_atm )).length() <<
		//			" dist to H41 " << (rsd.xyz( rsd.atom_index(" H41") ) - other_rsd.xyz( other_atm )).length()  << std::endl;
		if ( (rsd.xyz( rsd.atom_index(" H42") ) - other_rsd.xyz( other_atm )).length()  <
				 (rsd.xyz( rsd.atom_index(" H41") ) - other_rsd.xyz( other_atm )).length() ) {
			N_H++;
		} else {
			N_W++;
		}
	}

	if (rsd.aa() == na_rgu && atom_name == " N2 ") {
		if ( (rsd.xyz( rsd.atom_index(" H22") ) - other_rsd.xyz( other_atm )).length()  <
				 (rsd.xyz( rsd.atom_index(" H21") ) - other_rsd.xyz( other_atm )).length() ) {
			N_S++;
		} else {
			N_W++;
		}
	}

}

//////////////////////////////////
bool
atom_is_polar( core::conformation::Residue const & rsd, Size const & atm )
{
	for ( chemical::AtomIndices::const_iterator
					anum  = rsd.Hpos_polar().begin(),
					anume = rsd.Hpos_polar().end(); anum != anume; ++anum ) {
		Size const H_atm( *anum );
		if ( H_atm == atm ) return true;
	}
	return false;

}

//////////////////////////////////
bool
heavy_atom_is_polar( core::conformation::Residue const & rsd, Size const & atm )
{
	for ( chemical::AtomIndices::const_iterator
					anum  = rsd.Hpos_polar().begin(),
					anume = rsd.Hpos_polar().end(); anum != anume; ++anum ) {
		Size const H_atm( *anum );
		if ( rsd.atom_base( H_atm ) == atm ) return true;
	}
	return false;

}

//////////////////////////////////
bool
atom_is_acceptor( core::conformation::Residue const & rsd, Size const & atm )
{
	for ( chemical::AtomIndices::const_iterator
					anum  = rsd.accpt_pos().begin(),
					anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const H_atm( *anum );
		if ( H_atm == atm ) return true;
	}
	return false;

}


////////////////////////////////////////////////////////
// Look at contacts made by i to base on residue j
// Tabulate any contacts (not just H-bonds and possible CH-bonds)
////////////////////////////////////////////////////////
void
figure_out_number_base_contacts(
	conformation::Residue const & rsd_i,
	conformation::Residue const & rsd_j,
	Size & edge_classification )
{
	using namespace core::scoring::rna;

	static Real const DIST_CUTOFF(  4.2 );
	//	static Real const ANGLE_CUTOFF( numeric::conversions::radians( 100.0 ) );
	//Size n_hbonds( 0 );
	Size N_W( 0 ), N_H( 0 ), N_S( 0 );

	//	Size const i = rsd_i.seqpos();
	//	Size const j = rsd_j.seqpos();

	//std::cout << i << " <--> " << j << std::endl;
	// heavy atom on base j; heavy atom on i.
	for ( Size k = rsd_j.first_sidechain_atom()+1 ; k <= rsd_j.nheavyatoms(); k++ ) {

		//		if ( k <= rsd_j.first_sidechain_atom() ) continue;

		for ( Size m = 1; m <= rsd_i.nheavyatoms(); m++ ) {

			Real const dist_ij = ( rsd_i.xyz( m ) - rsd_j.xyz( k ) ).length();

			//			Real const angle_ij = angle_radians( rsd_i.xyz( acc_atm ),
			//																					 rsd_j.xyz( don_h_atm ),
			//																					 rsd_j.xyz( don_atm ) );
			if ( dist_ij < DIST_CUTOFF  /*&& angle_ij > ANGLE_CUTOFF*/ )				{
				update_edge_hbond_numbers( rsd_i, m, N_W, N_H, N_S );

				if ( atom_is_acceptor( rsd_j, k ) && heavy_atom_is_polar( rsd_i, m) ) {
					update_edge_hbond_numbers_careful_hydrogen( rsd_i, m, rsd_j, k, N_W, N_H, N_S ); //helps resolve confusion for A, C, and G bifurcated
				}

			}

		}

	}


	//std::cout << i << " <--> " << j << std::endl;

	// acceptor on base j; donor on i.
	//	for ( Size acc_atm = 1; acc_atm <= rsd_j.nheavyatoms(); acc_atm++ ) {
	//
	//		if ( acc_atm <= rsd_j.first_sidechain_atom() ) continue;
	//
	//		for ( Size k = rsd_i.nheavyatoms() + 1; k <= rsd_i.natoms(); k++ ) {
	//
	//			Size const & don_h_atm = k;
	//			Size const & don_atm = rsd_i.atom_base( k );
	//
	//			if ( ( rsd_j.xyz( acc_atm ) - rsd_i.xyz( don_h_atm ) ).length() < DIST_CUTOFF  /*&&
	//					 angle_radians( rsd_j.xyz( acc_atm ),
	//													rsd_i.xyz( don_h_atm ),
	//													rsd_i.xyz( don_atm ) ) > ANGLE_CUTOFF*/ )				{
	//				update_edge_hbond_numbers( rsd_i, don_atm, N_W, N_H, N_S );
	//				//if ( atom_is_polar( rsd_i, k ) && atom_is_acceptor( rsd_j, acc_atm) &&
	//				//						 don_atm > rsd_i.first_sidechain_atom() )  {
	//					//					if( (i == 4 && j ==9) ||  (i == 9 && j ==4)  )std::cout << "2. FOUND H-BOND " << i << " " << j << " " << rsd_i.atom_name( k ) << " " << rsd_j.atom_name( acc_atm ) << std::endl;
	//				//					n_hbonds++;
	//				//				}
	//			}
	//
	//		}
	//	}

	if ( N_W >= N_H && N_W >= N_S ) {
		edge_classification = WATSON_CRICK;
	} else if (N_H >= N_S ) {
		edge_classification = HOOGSTEEN;
	} else {
		edge_classification = SUGAR;
	}

	//	return n_hbonds;

}



/////////////////////////////////////////////////////////////////////////////
typedef  numeric::xyzMatrix< Real > Matrix;

/////////////////////////////////////////////////////////////////////////////
Size
figure_out_base_pair_orientation(
  core::pose::Pose & pose,
	Size const & i,
	Size const & j )
{
	using namespace core::scoring::rna;

	RNA_ScoringInfo  & rna_scoring_info( nonconst_rna_scoring_info_from_pose( pose ) );
	RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );

	//utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & z_i = M_i.col_z();

	kinematics::Stub const & stub_j( base_stubs[j] );
	Matrix const & M_j( stub_j.M );
	Vector const & z_j = M_j.col_z();
	Real const cos_theta = dot_product( z_i, z_j );

	Size const theta_bin = (cos_theta < 0) ?  1 : 2;
	return theta_bin;

}

////////////////////////////////////////////////////////////////
bool
residue_is_bulge( pose::Pose const & pose, Size const i )
{

	conformation::Residue const & rsd_i ( pose.residue( i ) ) ;
	static Real const DIST_CUTOFF( 4.0 );

	for ( Size k = rsd_i.first_sidechain_atom()+1; k <= rsd_i.nheavyatoms(); k++ ) {
		for ( Size j = 1 ; j <= pose.total_residue(); j++ ) {
			if ( i == j ) continue;
			for ( Size m = 1; m <= pose.residue( j ).nheavyatoms(); m++ ) {
				if ( ( rsd_i.xyz( k ) - pose.residue( j ).xyz( m ) ).length() < DIST_CUTOFF ) {
					//					std::cout << "Residue " << i << " --> " << j << " " << pose.residue(j).atom_name( m ) << std::endl;
					return false;
				}
			}
		}
	}

	return true;

}

/////////////////////////////////////////////////////////////////////////////
Size
bases_form_a_hydrogen_bond( core::scoring::hbonds::HBondSetOP const & hbond_set,
														core::pose::Pose & pose,
														Size const & i,
														Size const & j )
{

	Size num_hbonds = 0;

	static Real const HBOND_CUTOFF( -0.1 );
	for (Size n = 1; n <= hbond_set->nhbonds(); n++ ) {
		core::scoring::hbonds::HBond const & hbond( hbond_set->hbond( n ) );

		Size const & don_res_num = hbond.don_res();
		Size const & don_hatm = hbond.don_hatm();

		Size const & acc_res_num = hbond.acc_res();
		Size const & acc_atm = hbond.acc_atm();

		if ( don_res_num == i && acc_res_num == j ) {
			if ( pose.residue( i ).atom_base( don_hatm ) > pose.residue( i ).first_sidechain_atom() &&
					 acc_atm > pose.residue( j ).first_sidechain_atom() ) {

				//if ( ( i == 4  && j == 10 ) || ( i==10 && j==4) )		std::cout << i << "--" << j << "  -> " << pose.residue( i ).atom_name( don_hatm) << " " << pose.residue( j ).atom_name( acc_atm ) << hbond.energy() << std::endl;
				if ( hbond.energy() <= HBOND_CUTOFF )	num_hbonds++;

			}
		}

		if ( don_res_num == j && acc_res_num == i ) {
			if ( pose.residue( j ).atom_base( don_hatm ) > pose.residue( j ).first_sidechain_atom() &&
					 acc_atm > pose.residue( i ).first_sidechain_atom() ) {

				//if ( ( i == 4  && j == 10 ) || ( i==10 && j==4) )		 std::cout << j << "--" << i << "  -> " << pose.residue( j ).atom_name( don_hatm) << " " << pose.residue( i ).atom_name( acc_atm ) << hbond.energy() << std::endl;
				if ( hbond.energy() <= HBOND_CUTOFF )	num_hbonds++;

			}
		}

	}

	if (num_hbonds > 0 ) return num_hbonds;

	// Last check -- there may be unusually *close* contacts in lo-res structures that don't get good hbond scores.
	static Real const DIST_CUTOFF( 3.0 );
	conformation::Residue const & rsd_i( pose.residue( i ) ) ;
	conformation::Residue const & rsd_j( pose.residue( j ) ) ;

	for ( chemical::AtomIndices::const_iterator
					anum  = rsd_i.accpt_pos().begin(),
					anume = rsd_i.accpt_pos().end(); anum != anume; ++anum ) {

		Size const aatm( *anum ) ;
		if ( aatm <= rsd_i.first_sidechain_atom() ) continue;

		for ( chemical::AtomIndices::const_iterator
						hnum  = rsd_j.Hpos_polar().begin(),
						hnume = rsd_j.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const hatm( *hnum );

			if ( rsd_j.atom_base( hatm ) <= rsd_j.first_sidechain_atom() ) continue;

			if ( ( rsd_i.xyz( aatm ) - rsd_j.xyz( hatm ) ).length() < DIST_CUTOFF  ) {
				//				if ( ( i == 4  && j == 10 ) || ( i==10 && j==4) )		 std::cout << j << "--" << i << "  -> " << pose.residue( j ).atom_name( hatm) << " " << pose.residue( i ).atom_name( aatm ) << " blah" << std::endl;
				return 1;
			}
		}
	}


	for ( chemical::AtomIndices::const_iterator
					anum  = rsd_j.accpt_pos().begin(),
					anume = rsd_j.accpt_pos().end(); anum != anume; ++anum ) {

		Size const aatm( *anum ) ;
		if ( aatm <= rsd_j.first_sidechain_atom() ) continue;

		for ( chemical::AtomIndices::const_iterator
						hnum  = rsd_i.Hpos_polar().begin(),
						hnume = rsd_i.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const hatm( *hnum );

			if ( rsd_i.atom_base( hatm ) <= rsd_i.first_sidechain_atom() ) continue;

			if ( ( rsd_j.xyz( aatm ) - rsd_i.xyz( hatm ) ).length()  < DIST_CUTOFF  ) {
				//				if ( ( i == 4  && j == 10 ) || ( i==10 && j==4) )		 std::cout << i << "--" << j << "  -> " << pose.residue( i ).atom_name( hatm) << " " << pose.residue( j ).atom_name( aatm ) << " blah " << dist << std::endl;
				return 1;
			}
		}
	}


	return false;
}

//////////////////////////////////////////////////
bool
bases_are_coplanar(
										core::pose::Pose & pose,
										Size const & i,
										Size const & j )
{

	using namespace core::scoring::rna;

	RNA_ScoringInfo  & rna_scoring_info( nonconst_rna_scoring_info_from_pose( pose ) );
	RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.update( pose );

	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	//Vector const & x_i = M_i.col_x();
	//Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	Vector const & centroid_j( base_centroids[j] );
	kinematics::Stub const & stub_j( base_stubs[j] );

	Vector d_ij = centroid_j - centroid_i;
	//Real const dist_x = dot_product( d_ij, x_i );
	//	Real const dist_y = dot_product( d_ij, y_i );
	Real const dist_z = dot_product( d_ij, z_i );

	Matrix const & M_j( stub_j.M );
	Vector const & z_j = M_j.col_z();
	Real const cos_theta = dot_product( z_i, z_j );

	//	if ( i == 8  && j == 11 )		std::cout << "DIST_Z COS_THETA " << dist_z << " " << cos_theta << std::endl;
	//	if ( i == 5  && j == 8 )		std::cout << "DIST_Z COS_THETA " << dist_z << " " << cos_theta << std::endl;
	//	if ( j == 5  && i == 8 )		std::cout << "DIST_Z COS_THETA " << dist_z << " " << cos_theta << std::endl;


	static Real const rna_basepair_stagger_cutoff_( 2.8 );
	static Real const COS_THETA_CUTOFF( 0.6 );
	if ( std::abs(dist_z) < rna_basepair_stagger_cutoff_  && std::abs( cos_theta ) > COS_THETA_CUTOFF ) return true;

	//
	//	static Real const rna_basepair_stagger_cutoff_loose_( 3.5 );
	//	static Real const COS_THETA_CUTOFF_STRICT( 0.65 );
	//	if ( std::abs(dist_z) < rna_basepair_stagger_cutoff_loose_ && std::abs( cos_theta ) > COS_THETA_CUTOFF_STRICT ) return true;

	return false;
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void
classify_base_pairs(
    core::pose::Pose const & pose_input,
		utility::vector1< core::scoring::rna::Base_pair> & base_pair_list,
		utility::vector1< bool > & is_bulged
)
{
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::chemical;

	base_pair_list.clear();

	pose::Pose pose = pose_input;

	//////////////////////////////////////////////////////////////
	ObjexxFCL::FArray1D < bool > is_base_paired( pose.total_residue(), false );

	// Get hydrogen bond list.
	//ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS ) );
	ScoreFunctionOP score_fxn = new ScoreFunction;
	score_fxn->set_weight( hbond_sc, 1.0);
	(*score_fxn)(pose);

	hbonds::HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
	hbond_options->use_hb_env_dep( false );
	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( hbond_options ));

	hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

	//	std::cout << "---------" << std::endl;

	//////////////////////////////////////////////////////////////
	for (Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( ! pose.residue(i).is_RNA()  ) continue;
		for (Size j = i+1; j <= pose.total_residue(); j++ ) {
			if ( ! pose.residue(j).is_RNA()  ) continue;

			Size const num_hbonds = bases_form_a_hydrogen_bond( hbond_set, pose, i, j );
			if ( num_hbonds == 0  ) continue;
			if ( ! bases_are_coplanar( pose, i, j ) || ! bases_are_coplanar( pose, j, i ) ) continue;

			Size edge_classification_i( 0 ), edge_classification_j( 0 );

			figure_out_number_base_contacts(  pose.residue( i ), pose.residue( j ), edge_classification_i );
			figure_out_number_base_contacts(  pose.residue( j ), pose.residue( i ), edge_classification_j );

			Size const orientation = figure_out_base_pair_orientation( pose, i, j );

			//These pernicious bifurcated hydrogen bonds.
			//			if ( num_hbonds == 1 &&
			//					 ( pose.residue(i).aa() == na_rad  || pose.residue(i).aa() == na_rgu ) &&
			//					 ( pose.residue(j).aa() == na_rad  || pose.residue(j).aa() == na_rgu ) ) {
			//				edge_classification_i = WATSON_CRICK;
			//				edge_classification_j = WATSON_CRICK;
			//			}


			//			if ( n_i > 0 && n_j > 0 ) {
			base_pair_list.push_back( core::scoring::rna::Base_pair( i, j, edge_classification_i, edge_classification_j , orientation ) );
			if ( false ) std::cout << pose.residue( i ).name1() << i << " " << pose.residue(j).name1() << j << "  " << get_edge_from_num( edge_classification_i ) << " " << get_edge_from_num( edge_classification_j ) << " " << orientation << std::endl;

				//			}

		}
	}

	//////////////////////////////////////////////////////////////
	is_bulged.clear();
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		bool const check_bulge = residue_is_bulge( pose, i );
		is_bulged.push_back( check_bulge );
		if ( check_bulge && false ) {
			std::cout << " BULGE " << i << std::endl;
		}
	}

}

//////////////////////////////////////////////////////////////////////
Size
get_number_base_stacks(
	 core::pose::Pose const & pose_input
)
{
	using namespace core::scoring;
	using namespace core::scoring::rna;

	ScoreFunctionOP denovo_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_LORES_WTS ) );

	pose::Pose pose = pose_input;

	(*denovo_scorefxn)( pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_stack_list const & scored_base_stack_list( rna_filtered_base_base_info.scored_base_stack_list() );

	return scored_base_stack_list.size();

}

} // namespace rna
} // namespace protocols
