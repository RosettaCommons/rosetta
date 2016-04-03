// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffractiobn/FiberDiffractionEnergy.cc
/// @brief  FiberDiffraction scoring and derivatives (all-atom)
/// @author Wojciech Potrzebowski and Ingemar Andre


// Unit headers
#include <core/scoring/fiber_diffraction/FiberDiffractionEnergy.hh>
#include <core/scoring/fiber_diffraction/FiberDiffraction.hh>
#include <core/scoring/fiber_diffraction/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/fiber_diffraction/xray_scattering.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>

#include <core/id/NamedAtomID.hh>
#include <core/io/pdb/build_pose_as_is.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/prof.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/basic_sys_util.hh>

// C++
#include <iomanip>
#include <string>
//#include <sys/time.h>

#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

namespace core {
namespace scoring {
namespace fiber_diffraction {

// tracer
static basic::Tracer TR("core.scoring.fiber_diffraction.FiberDiffractionEnergy");

ScoreTypes FiberDiffractionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fiberdiffraction );
	return sts;
}

methods::EnergyMethodOP FiberDiffractionEnergyCreator::create_energy_method( methods::EnergyMethodOptions const &) const {
	return methods::EnergyMethodOP( new FiberDiffractionEnergy() );
}

FiberDiffractionEnergy::FiberDiffractionEnergy() :
	parent( methods::EnergyMethodCreatorOP( new FiberDiffractionEnergyCreator ) ) {
	chi2_=0;
	dchi2_d.clear();
	dchi2_d_cross_R.clear();
}

/// clone
methods::EnergyMethodOP FiberDiffractionEnergy::clone() const {
	return methods::EnergyMethodOP( new FiberDiffractionEnergy( *this ) );
}

void FiberDiffractionEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {

	if ( !pose.is_fullatom() ) return;

	if ( !core::pose::symmetry::is_symmetric(pose) ) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	// load fiber diffraction data and bessel functions orders.
	// Bessel orders satisfies helix selection rule
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_I;
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_R;
	utility::vector0< utility::vector0 < int > >::iterator nvals;
	core::Size lmax, Rmax;

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user() ) {
		a_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user() ) {
		b_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
	}

	if ( !(a_ > 0 ) ) {
		utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
	}

	if ( !(b_ > 0 ) ) {
		utility_exit_with_message("The number of turns per repeat, score::fiber_diffraction::b, must be set!");
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ].user() ) {
		p_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ]();
	} else {
		find_pitch( pose, p_ );
	}
	c_ = p_*a_;

	core::Real  res_cutoff_low_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_low ]();
	core::Real  res_cutoff_high_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_high ]();

	core::Real b_factor_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor ]();
	core::Real b_factor_solv_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor_solv ]();
	core::Real b_factor_solv_K_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor_solv_K ]();

	getFiberDiffractionData(c_, res_cutoff_high_, res_cutoff_low_).getAllFiberData(layer_lines_I, layer_lines_R, nvals, lmax, Rmax);
	TR << "Rmax : lmax " << "( " << Rmax << " : " << lmax << " )" << std::endl;

	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_(
		setup_form_factors( pose, lmax, layer_lines_R, c_, b_factor_, b_factor_solv_, b_factor_solv_K_ ));
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator form_factors(form_factors_.begin());

	core::Real max_r_value;
	find_max_r( pose, max_r_value );
	TR<<"Max r value "<<max_r_value<<std::endl;
	I.resize(lmax+1);
	for ( Size l=0; l<= lmax; ++l ) {
		I[l].resize(Rmax);
		for ( Size i=1; i<= Rmax; ++i ) I[l][i] = 0.0;
	}

	TR << "Preparing fiber model data..." << std::endl;
	Size natoms(0);
	utility::vector1< Size > atom_type_number;
	utility::vector1< Real > phi, z, r, bfactors;
	setup_cylindrical_coords( pose, natoms, atom_type_number, AtomID_to_atomnbr_, phi, z, r, bfactors);
	TR << "Model contains " << natoms << " atoms" << std::endl;

	TR << "Calculating Chi2..." << std::endl;
	Size n_iter(0);
	Real Rinv(0);

	for ( Size l=0; l <= lmax; ++l ) {
		Size max_b_order( nvals[l].size() );
		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {

			int n( nvals[l][b_order-1] );
			int abs_n( std::abs(n) );
			// Precalculate phases. Initialize vector
			utility::vector1< utility::vector1 < Real > > phases;
			phases.resize(natoms);
			for ( Size atom=1; atom <= natoms; ++atom ) phases[atom].resize(natoms);

			// store phases for each value of n
			for ( Size atom1=1; atom1 <= natoms; ++atom1 ) {
				for ( Size atom2=atom1+1; atom2 <= natoms; ++atom2 ) {
					Real phase ( cos( n*(phi[atom2]-phi[atom1] )+2*M_PI*l/c_*( z[atom1]-z[atom2] ) ) );
					phases[atom1][atom2] = phase;
				}
			}

			// loop over all R values
			Size max_R_values( layer_lines_R[l].size() );
			for ( Size R=1; R<= max_R_values; ++R ) {
				Rinv=layer_lines_R[l][R];
				Real x_factor( 2*M_PI*Rinv );
				// precalculate jn for each n,R
				utility::vector1< Real > jn_vec;
				jn_vec.resize(natoms);
				// Calculate scattering...

				for ( Size atom1=1; atom1 <= natoms; ++atom1 ) {
					// Calculate only the jn values once.
					if ( atom1 == 1 ) {
						for ( Size at=1; at <= natoms; ++at ) {
							jn_vec[at] = jn(n,x_factor*r[at]);
						}
					}
					Real X1 (x_factor*r[atom1]);
					// calculate R_min
					if ( abs_n > X1 +2 ) continue;
					Real jn1 (jn_vec[atom1]);
					Size atom_type_index1 ( atom_type_number[atom1] );
					Real dummy(  form_factors[l][atom_type_index1][R]*jn1 );
					Real dummy_sum( 0 );
					I[l][R] += dummy*form_factors[l][atom_type_index1][R]*jn1;
					for ( Size atom2=atom1+1; atom2 <= natoms; ++atom2 ) {
						Real X2 (x_factor*r[atom2]);
						if ( abs_n > X2 +2 ) continue;
						Real jn2 (jn_vec[atom2]);
						Size atom_type_index2 ( atom_type_number[ atom2 ] );
						dummy_sum += form_factors[l][atom_type_index2][R]*jn2*phases[atom1][atom2];
						++n_iter;
					} // atoms2
					I[l][R] +=2*dummy*dummy_sum;
				} // atom1
			} // b_order
		}
	}


	Real prod(0);
	square_obs_ = 0;
	sum_obs_ = 0;

	for ( Size l=0; l <= lmax; ++l ) {
		for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
			if ( I[l][R] < -10.0 ) { utility_exit_with_message("Intensity below -10. Something went wrong...");}
			if ( I[l][R] < 0 ) { I[l][R]=0.0;}
			Real I_obs_ ( layer_lines_I[l][R]*layer_lines_I[l][R] );
			prod += I[l][R]*I_obs_;
			square_obs_ += I_obs_*I_obs_;
			sum_obs_ +=  I_obs_;
		}
	}

	scale_factor_ = square_obs_/prod;

	bool output_fiber_spectra_(false);
	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::output_fiber_spectra ].user() ) {
		output_fiber_spectra_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::output_fiber_spectra ]();
	}
	std::ofstream out;
	if ( output_fiber_spectra_ ) {
		std::string outfile = "Intensity.txt";
		out.open(outfile.c_str(), std::ios::out);
	}

	/////////////////////////////Rfacor//////////////////////////////////
	bool rfactor_refinement_=false;
	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ].user() ) {
		rfactor_refinement_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ]();
	}

	if ( rfactor_refinement_ ) {
		Rsum_obs = 0;
		Real Rfactor(0);
		Real Rprod(0);
		Real Rsquare_obs(0);
		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				if ( I[l][R] > 0.0 ) {
					Rprod += sqrt(I[l][R])*fabs(layer_lines_I[l][R]);
					Rsquare_obs += layer_lines_I[l][R]*layer_lines_I[l][R];
					Rsum_obs += fabs(layer_lines_I[l][R]);
				}
			}
		}
		Rscale_factor =Rsquare_obs/Rprod;
		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				if ( I[l][R] > 0.0 ) {
					Rfactor +=  fabs(Rscale_factor*sqrt(I[l][R])-fabs(layer_lines_I[l][R]));
				}
			}
		}
		TR<<"Rfactor (unnormalized): "<<Rfactor << std::endl;
		Rfactor /=Rsum_obs;
		TR<<"Rfactor: "<<Rfactor << " scale factor " << Rscale_factor <<" sum obsreved " << Rsum_obs << std::endl;
		scale_factor_ = Rscale_factor;
		chi2_ = Rfactor;
	} else {
		///////////////////////////END Rfactor/////////////////////////////
		chi2_=0;
		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				chi2_ +=  (scale_factor_*I[l][R]-layer_lines_I[l][R]*layer_lines_I[l][R])*(scale_factor_*I[l][R]-layer_lines_I[l][R]*layer_lines_I[l][R]);
				if ( output_fiber_spectra_ ) {
					out << scale_factor_*I[l][R] << " " << layer_lines_R[l][R] <<" "<< l << std::endl;
				}
			}
		}
		chi2_ /=square_obs_;
	}

	TR << "chi2 " << chi2_ << " sum_obs_ " << sum_obs_ << " scale_factor " <<  scale_factor_ << std::endl;

	TR << " number of iterations " << n_iter  << std::endl;

	////////////////////////////////////////////

	core::Size chi_iterations_(0);
	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::chi_iterations ].user() ) {
		chi_iterations_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::chi_iterations ]();
	}

	if ( chi_iterations_ > 0 ) {
		core::Real max_r_value;
		find_max_r( pose, max_r_value );

		utility::vector0< utility::vector1< core::Real > >::iterator  sampling_points_lE;
		utility::vector0< core::Size > ::iterator  sampling_points_number_l;
		utility::vector0< core::Size > ::iterator  lowest_bessel_orders_l;
		utility::vector0< core::Size > ::iterator  lowest_resolution_l;
		utility::vector0< core::Size > ::iterator  highest_resolution_l;
		utility::vector0< utility::vector1< core::Real > >::iterator  selected_Rinv_l;

		core::Real chi_free(0);

		chi_free = calculate_chi2_free(pose, chi_iterations_,lmax,\
			layer_lines_I,layer_lines_R,natoms,c_,nvals,\
			atom_type_number,phi,z,r,b_factor_, b_factor_solv_, b_factor_solv_K_);

		TR<<"Chi free: "<<chi_free<<" after "<<chi_iterations_<<" iterations."<<std::endl;
		///////////////////////////////////////////
	}

	if ( output_fiber_spectra_ ) {
		out.close();
	}
}

void FiberDiffractionEnergy::finalize_total_energy(pose::Pose & /*pose*/, ScoreFunction const &, EnergyMap & emap) const {
	// all the work is done in setup for scoring
	// just return results here
	emap[ fiberdiffraction ] += chi2_;
}

void FiberDiffractionEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {

	if ( !pose.is_fullatom() ) return;

	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_I;
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_R;
	utility::vector0< utility::vector0 < int > >::iterator nvals;
	core::Size lmax, Rmax;

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user() ) {
		a_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user() ) {
		b_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
	}

	if ( !(a_ > 0 ) ) {
		utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
	}

	if ( !(b_ > 0 ) ) {
		utility_exit_with_message("The number of turns per repeat, score::fiber_diffraction::b, must be set!");
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ].user() ) {
		p_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ]();
	} else {
		find_pitch( pose, p_ );
	}
	c_ = p_*a_;

	core::Real res_cutoff_low_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_low ]();
	core::Real res_cutoff_high_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::resolution_cutoff_high ]();

	core::Real b_factor_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor ]();
	core::Real b_factor_solv_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor_solv ]();
	core::Real b_factor_solv_K_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b_factor_solv_K ]();

	getFiberDiffractionData(c_, res_cutoff_high_, res_cutoff_low_).getAllFiberData(layer_lines_I, layer_lines_R, nvals, lmax, Rmax);
	TR << "Rmax : lmax " << "( " << Rmax << " : " << lmax << " )" << std::endl;


	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_(
		setup_form_factors( pose, lmax, layer_lines_R, c_, b_factor_, b_factor_solv_, b_factor_solv_K_ ));
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator form_factors(form_factors_.begin());

	TR << "Preparing fiber model data for deriv calculation..." << std::endl;

	Size natoms(0);
	utility::vector1< Size > atom_type_number;
	utility::vector1< Real > phi, z, r,bfactors;
	setup_cylindrical_coords( pose, natoms, atom_type_number, AtomID_to_atomnbr_, phi, z, r, bfactors);

	bool rfactor_refinement_=false;
	if ( basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ].user() ) {
		rfactor_refinement_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ]();
	}

	TR << "Calculating deriv Chi2..." << std::endl;
	if ( rfactor_refinement_ ) {
		TR << "Scale factor deriv? " << Rscale_factor <<" square_obs " << Rsum_obs << std::endl;
	} else {
		TR << "Scale factor deriv? " << scale_factor_ <<" square_obs " <<square_obs_ << std::endl;
	}
	dchi2_d.resize(natoms);
	dchi2_d_cross_R.resize(natoms);
	for ( Size i=1; i <= natoms; ++i ) {
		dchi2_d[i] = 0.0;
		dchi2_d_cross_R[i] = 0.0;
	}

	///////////////////////////////////////////
	Size n_iter(0);
	Real Rinv(0);
	for ( Size l=0; l <= lmax; ++l ) {
		Size max_b_order( nvals[l].size() );

		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
			int n( nvals[l][b_order-1] );
			int abs_n( abs(n) );
			// Precalculate phases. Initialize vector
			utility::vector1< utility::vector1 < Real > > phases, phases_prime;
			phases.resize(natoms);
			phases_prime.resize(natoms);
			for ( Size atom=1; atom <= natoms; ++atom ) {
				phases[atom].resize(natoms);
				phases_prime[atom].resize(natoms);
			}

			// store phases for each value of n
			for ( Size atom1=1; atom1 <= natoms; ++atom1 ) {
				for ( Size atom2=atom1+1; atom2 <= natoms; ++atom2 ) {
					Real dummy( n*(phi[atom2]-phi[atom1] )+2*M_PI*l/c_*( z[atom1]-z[atom2] ) );
					Real phase ( cos( dummy) );
					Real phase_prime (  -sin( dummy ) );
					phases[atom1][atom2] = phase;
					phases[atom2][atom1] = phase;
					phases_prime[atom1][atom2] = phase_prime;
					phases_prime[atom2][atom1] = -phase_prime;
				}
			}

			// loop over all R values
			Size max_R_values( layer_lines_R[l].size() );
			for ( Size R=1; R<= max_R_values; ++R ) {
				Rinv=layer_lines_R[l][R];
				Real x_factor( 2*M_PI*Rinv );

				// precalculate jn for each n,R
				utility::vector1< Real > jn_vec, jn_vec_plus_1;
				jn_vec.resize(natoms);
				jn_vec_plus_1.resize(natoms);

				// Calculate scattering...
				for ( Size atom1=1; atom1 <= natoms; ++atom1 ) {
					// Calculate only the jn values once.
					if ( atom1 == 1 ) {
						for ( Size at=1; at <= natoms; ++at ) {
							jn_vec[at] = jn(n, x_factor*r[at]);
							jn_vec_plus_1[at] = jn(n+1, x_factor*r[at]);
						}
					}

					Real X1 (x_factor*r[atom1]);

					// calculate R_min
					if ( abs_n > X1 +2 ) continue;
					Real jn1 (jn_vec[atom1]);
					Real jn1_plus_1 (jn_vec_plus_1[atom1]);
					Size atom_type_index1 ( atom_type_number[atom1] );
					Real Jn_deriv_atom1( -x_factor*jn1_plus_1+n*jn1/r[atom1] );

					if ( fabs(r[atom1] ) < 1e-2 ) Jn_deriv_atom1=0.0;

					// cartesian coords and unit vectors...
					numeric::xyzVector< core::Real > cartesian_coord_atom1( r[atom1]*cos(phi[atom1]), r[atom1]*sin(phi[atom1]), z[atom1] );
					numeric::xyzVector< core::Real > unit_r(cos(phi[atom1]), sin(phi[atom1]), 0 );
					numeric::xyzVector< core::Real > unit_x(1, 0, 0 );
					numeric::xyzVector< core::Real > unit_y(0, 1, 0 );
					numeric::xyzVector< core::Real > unit_z(0, 0, 1 );
					numeric::xyzVector< core::Real > D(0,0,0);
					numeric::xyzVector< core::Real > D_cross_R(0,0,0);

					Real tmp( 2*form_factors[l][atom_type_index1][R]*form_factors[l][atom_type_index1][R]*jn1*Jn_deriv_atom1);
					D += tmp*unit_r;

					//cross of ATOM1 X ATOM1?
					D_cross_R += D.cross(cartesian_coord_atom1);
					for ( Size atom2=1; atom2 <= natoms; ++atom2 ) {
						if ( atom2 == atom1 ) continue;

						Real X2 (x_factor*r[atom2]);
						if ( abs_n > X2 +2 ) continue;

						Real jn2 (jn_vec[atom2]);
						Size atom_type_index2 ( atom_type_number[ atom2 ] );
						Real fact( form_factors[l][atom_type_index1][R]*form_factors[l][atom_type_index2][R]*jn2 );

						//Jn_deriv_atom1
						Real dr( Jn_deriv_atom1*fact*phases[atom1][atom2] );
						Real dphi( n*jn1*fact*phases_prime[atom1][atom2] );
						Real dz ( 2*M_PI*l/c_*jn1*fact*phases_prime[atom1][atom2] );

						numeric::xyzVector< core::Real > D_tmp(0,0,0);
						numeric::xyzVector< core::Real > dphi_vec(0,0,0);

						if ( r[atom1] >= 1e-2 ) {
							if ( fabs(sin(phi[atom1])) > 1e-3 ) {
								dphi_vec = (unit_x-cos(phi[atom1])*unit_r)/(sin(phi[atom1])*r[atom1]);
							} else {
								dphi_vec = numeric::xyzVector< core::Real > (0,-1/r[atom1],0);
							}
						}
						D_tmp += 2*(dr*unit_r + dz*unit_z + dphi*dphi_vec);

						numeric::xyzVector< core::Real > t1( dr*unit_r);
						numeric::xyzVector< core::Real > t2( dz*unit_z );
						numeric::xyzVector< core::Real > t3( dphi_vec  );
						numeric::xyzVector< core::Real > t4( (unit_x-cos(phi[atom1]+0e-6)*unit_r)  );

						D += D_tmp;
						D_cross_R += D_tmp.cross(cartesian_coord_atom1);
						++n_iter;
					} // atoms2
					if ( rfactor_refinement_ ) {
						if ( I[l][R]>0.0 ) {
							Real F_diff (Rscale_factor*sqrt(I[l][R]) - fabs(layer_lines_I[l][R]));
							D *= scale_factor_*F_diff/(2*Rsum_obs*fabs(F_diff)*sqrt(I[l][R]));
							D_cross_R *= scale_factor_*F_diff/(2*Rsum_obs*fabs(F_diff)*sqrt(I[l][R]));
						}
					} else {
						Real dummy3( layer_lines_I[l][R]*layer_lines_I[l][R] );
						Real I_diff (( scale_factor_*I[l][R] - dummy3 ));
						D *= 2*scale_factor_*I_diff/square_obs_;
						D_cross_R *= 2*scale_factor_*I_diff/square_obs_;
					}
					dchi2_d[atom1] += D;
					dchi2_d_cross_R[atom1] += D_cross_R;
				} // atom1
			} // b_order
		}
	}
	TR << "number of iterations " << n_iter << std::endl;
}

void FiberDiffractionEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {

	if ( !pose.is_fullatom() ) return;

	int resid = id.rsd();
	int atmid = id.atomno();
	core::conformation::Residue const &rsd_i = pose.residue(resid);
	numeric::xyzVector<core::Real> X = pose.xyz(id);

	if ( rsd_i.aa() != core::chemical::aa_vrt && !rsd_i.atom_type(atmid).is_heavyatom() ) return;

	std::map<  core::id::AtomID, core::Size >::iterator it( AtomID_to_atomnbr_.find( id ) );
	if ( it == AtomID_to_atomnbr_.end() ) return;
	Size atomnbr( it->second );

	numeric::xyzVector<core::Real> f1( dchi2_d_cross_R[atomnbr] );
	numeric::xyzVector<core::Real> f2( dchi2_d[atomnbr] );

	F1 += weights[ fiberdiffraction ] * f1;
	F2 += weights[ fiberdiffraction ] * f2;
}

core::Size
FiberDiffractionEnergy::version() const
{
	return 1; // Initial versioning
}

} // fiber_diffraction
} // scoring
} // core

