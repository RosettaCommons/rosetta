// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffractiobn/FiberDiffractionEnergyGpu.cc
/// @brief  FiberDiffraction scoring on GPU
/// @author Wojciech Potrzebowski and Ingemar Andre

// Unit headers
#ifdef USECUDA
#include <core/scoring/fiber_diffraction/FiberDiffractionEnergyGpu.hh>
#include <core/scoring/fiber_diffraction/FiberDiffractionKernelGpu.hh>
#include <core/scoring/fiber_diffraction/FiberDiffraction.hh>
#include <core/scoring/fiber_diffraction/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethod.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
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

namespace core {
namespace scoring {
namespace fiber_diffraction {

// tracer
static basic::Tracer TR("core.scoring.fiber_diffraction.FiberDiffractionEnergyGpu");

ScoreTypes FiberDiffractionEnergyGpuCreator::score_types_for_method() const {
	ScoreTypes sts;
	//Registering fiberdiffractiongpu weight
	sts.push_back( fiberdiffractiongpu );
	return sts;
}

methods::EnergyMethodOP FiberDiffractionEnergyGpuCreator::create_energy_method( methods::EnergyMethodOptions const &) const {
	return methods::EnergyMethodOP(new FiberDiffractionEnergyGpu() );
}

FiberDiffractionEnergyGpu::FiberDiffractionEnergyGpu() :
	parent( methods::EnergyMethodCreatorOP( new FiberDiffractionEnergyGpuCreator ) ) {
	chi2_=0;
	dchi2_d.clear();
	dchi2_d_cross_R.clear();
}

/// clone
methods::EnergyMethodOP FiberDiffractionEnergyGpu::clone() const {
	return methods::EnergyMethodOP( new FiberDiffractionEnergyGpu( *this ) );
}

void FiberDiffractionEnergyGpu::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {

	if (!pose.is_fullatom()) return;

	if (!core::pose::symmetry::is_symmetric(pose)) {
		utility_exit_with_message("Structure needs to be symmetric! Aborting...");
	}

	// load fiber diffraction data
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_I;
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_R;
	utility::vector0< utility::vector0 < int > >::iterator nvals;
	core::Size lmax, Rmax;

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user()) {
		a_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
	}

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user()) {
		b_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
	}

	if ( !(a_ > 0 ) ){
		utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
	}

	if ( !(b_ > 0 ) ){
		utility_exit_with_message("The number of turns in one repeat, score::fiber_diffraction::b, must be set!");
	}

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ].user()) {
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
	TR << "Trimmed Rmax : lmax " << "( " << Rmax << " : " << lmax << " )" << std::endl;

	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_(
	setup_form_factors( pose, lmax, layer_lines_R, c_, b_factor_, b_factor_solv_, b_factor_solv_K_ ));
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator form_factors(form_factors_.begin());

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

	int gpu_processor_=0;
	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::gpu_processor ].user()) {
		gpu_processor_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::gpu_processor ]();
	}

	//Caluculating all intensities on the gpu.
	calculate_intensity_gpu(lmax, natoms, nvals, layer_lines_R, I, form_factors, 
				 phi, z, r, atom_type_number, c_, res_cutoff_low_, res_cutoff_high_,
				gpu_processor_);

	Real prod(0);
	square_obs_ = 0;
	sum_obs_ = 0;

	for ( Size l=0; l <= lmax; ++l ) {
		for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
			if ( I[l][R] < -10.0 ) {
				utility_exit_with_message("Intensity lower than -10.0. Something went wrong...");
			}
			if ( I[l][R] < 0 ) {
				TR << "Due to rounding error, intenisty is slightly bellow zero..., so it is set to 0.0" << std::endl;
				I[l][R]=0.0;
			}
			Real I_obs_ ( layer_lines_I[l][R]*layer_lines_I[l][R] );
			prod += I[l][R]*I_obs_;
			square_obs_ += I_obs_*I_obs_;
			sum_obs_ +=  I_obs_;
		}
	}
	//Scale factor is crucial for numerical derivatives check
	//Even small differences cause significant dicsrapancies in derivatives
	//We found that by fixing scale factor you may correct numeric vs. analytical deriv check

	scale_factor_ = square_obs_/prod;

	bool output_fiber_spectra_(false);
	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::output_fiber_spectra ].user()) {
		output_fiber_spectra_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::output_fiber_spectra ]();
	}

	std::ofstream out;
	if ( output_fiber_spectra_ ) {
		std::string outfile = "IntensityGpu.txt";
		out.open(outfile.c_str(), std::ios::out);
	}

	bool rfactor_refinement_=false;
	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ].user()) {
		rfactor_refinement_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ]();
	}

	if (rfactor_refinement_) {
		Rsum_obs = 0;
		Real Rfactor(0);
		Real Rprod(0);
		Real Rsquare_obs(0);
		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				if (I[l][R] <= 0.0) continue;
				
				Rprod += sqrt(I[l][R])*fabs(layer_lines_I[l][R]);
				Rsquare_obs += layer_lines_I[l][R]*layer_lines_I[l][R];
				Rsum_obs += fabs(layer_lines_I[l][R]);
			}
		}
		Rscale_factor =Rsquare_obs/Rprod;

		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				if (I[l][R] > 0.0) {
					Rfactor +=  fabs(Rscale_factor*sqrt(I[l][R])-fabs(layer_lines_I[l][R]));
				}
			}
		}
		TR<<"Rfactor (unnormalized): "<<Rfactor << std::endl;
		Rfactor /=Rsum_obs;
		TR<<"Rfactor: "<<Rfactor << " scale factor " << Rscale_factor <<" sum obsreved " << Rsum_obs << std::endl;
		scale_factor_ = Rscale_factor;
		chi2_ = Rfactor;
		square_obs_ = Rsum_obs;
	}
	///////////////////////////END Rfactor/////////////////////////////     
	else {
		chi2_=0;
		for ( Size l=0; l <= lmax; ++l ) {
			for ( Size R=1; R<=layer_lines_R[l].size(); ++R ) {
				chi2_ +=  (scale_factor_*I[l][R]-layer_lines_I[l][R]*layer_lines_I[l][R])*(scale_factor_*I[l][R]-layer_lines_I[l][R]*layer_lines_I[l][R]);
				if (output_fiber_spectra_) {
					out << scale_factor_*I[l][R] << " " << layer_lines_R[l][R] <<" "<< l << std::endl;
				}
			}
		}
		chi2_ /=square_obs_;
	}

	TR << "Chi2: " << chi2_<< " sum_obs_ " << sum_obs_ <<" scale factor " << scale_factor_ << std::endl;
	
	if (output_fiber_spectra_) {
		out.close();
	}
}

void FiberDiffractionEnergyGpu::finalize_total_energy(pose::Pose & /*pose*/, ScoreFunction const &, EnergyMap & emap) const {
	emap[ fiberdiffractiongpu ] += chi2_;
}

void FiberDiffractionEnergyGpu::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const {

	if (!pose.is_fullatom()) return;

	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_I;
	utility::vector0< utility::vector1< core::Real > >::iterator layer_lines_R;
	utility::vector0 < utility::vector0 < int > >::iterator nvals;
	core::Size lmax, Rmax;

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ].user()) {
		a_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::a ]();
	}

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ].user()) {
 		b_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::b ]();
	}

	if ( !(a_ > 0 ) ){
		utility_exit_with_message("The number of subunits per repeat, score::fiber_diffraction::a, must be set!");
	}

	if ( !(b_ > 0 ) ){
		utility_exit_with_message("The number of turns per repeat, score::fiber_diffraction::b, must be set!");
	}

	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::p ].user()) {
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
	TR << "Trimmed Rmax : lmax " << "( " << Rmax << " : " << lmax << " )" << std::endl;

	utility::vector0< utility::vector1< utility::vector1< core::Real > > > form_factors_(
	setup_form_factors( pose, lmax, layer_lines_R, c_, b_factor_, b_factor_solv_, b_factor_solv_K_ ));
	utility::vector0< utility::vector1< utility::vector1< core::Real > > >::iterator form_factors(form_factors_.begin());

	bool rfactor_refinement_=false;
	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ].user()) {
		rfactor_refinement_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::rfactor_refinement ]();
	}

	int gpu_processor_=0;
	if (basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::gpu_processor ].user()) {
		gpu_processor_ = basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::gpu_processor ]();
	}

	TR << "Preparing fiber model data for deriv calculation..." << std::endl;

	Size natoms(0);
	utility::vector1< Size > atom_type_number;
	utility::vector1< Real > phi, z, r, bfactors;
	setup_cylindrical_coords( pose, natoms, atom_type_number, AtomID_to_atomnbr_, phi, z, r, bfactors);


	TR << "Calculating deriv Chi2..." << std::endl;
	TR << "Scale factor deriv? " << scale_factor_ << std::endl;
	dchi2_d.resize(natoms);
	dchi2_d_cross_R.resize(natoms);
	for ( Size i=1; i <= natoms; ++i) {
		dchi2_d[i] = 0.0;
		dchi2_d_cross_R[i] = 0.0;
	}

	/////////////////GPU action//////////////////////////
	calculate_derivatives_gpu(lmax, natoms, nvals, layer_lines_R, layer_lines_I, I, form_factors, 
				phi, z, r, atom_type_number, dchi2_d, dchi2_d_cross_R, c_, 
				res_cutoff_low_, res_cutoff_high_, scale_factor_, square_obs_,
				gpu_processor_, rfactor_refinement_);	

}

void FiberDiffractionEnergyGpu::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {

	if (!pose.is_fullatom()) return;

	int resid = id.rsd();
	int atmid = id.atomno();
	core::conformation::Residue const &rsd_i = pose.residue(resid);
	numeric::xyzVector<core::Real> X = pose.xyz(id);

	// if (hydrogen) return
	if ( rsd_i.aa() != core::chemical::aa_vrt && !rsd_i.atom_type(atmid).is_heavyatom() ) return;

	std::map<  core::id::AtomID, core::Size >::iterator it( AtomID_to_atomnbr_.find( id ) );
	if ( it == AtomID_to_atomnbr_.end() ) return;
	Size atomnbr( it->second );

	numeric::xyzVector<core::Real> f1( dchi2_d_cross_R[atomnbr] );
	numeric::xyzVector<core::Real> f2( dchi2_d[atomnbr] );

	F1 += weights[ fiberdiffractiongpu ] * f1;
	F2 += weights[ fiberdiffractiongpu ] * f2;
}

core::Size
FiberDiffractionEnergyGpu::version() const
{
	return 1; // Initial versioning
}

} // fiber_diffraction
} // scoring
} // core
#endif //end USECUDA
