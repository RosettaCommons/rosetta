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


#ifndef INCLUDED_core_pack_dunbrack_cenrot_SingleResidueCenrotLibrary_hh
#define INCLUDED_core_pack_dunbrack_cenrot_SingleResidueCenrotLibrary_hh

// Unit Headers
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.fwd.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>

#include <basic/basic.hh>

// Numeric Headers
#include <numeric/types.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

/// @brief Simple class storing all the data for one centroid-rotamer well
class CentroidRotamerSampleData
{
private:
	static Size const NUMBER_OF_PARAMS;

public:
	typedef utility::fixedsizearray1< Real, 3 > DOF3;

	//prob, dis, ang, dih, sd_dis, sd_ang, sd_dih
	Real *data_;
	Real *deriv_phi_;
	Real *deriv_psi_;

private:
	Real prob_;		/// probability for sidechain showing in this well
	Real distance_;	/// CEN-CB bond length
	Real angle_;	/// CEN-CB-CA bond angle
	Real dihedral_;	/// CEN-CB-CA-N dihedral angle chi1
	Real sd_dis_; 		/// standard deviation of dis
	Real sd_ang_; 		/// standard deviation of ang
	Real sd_dih_; 		/// standard deviation of dih

	Real energy_;	/// -ln(P)
	Real norm_factor_; /// normalization factor sigma=sqrt(2pi)

public:
	CentroidRotamerSampleData(){
		data_ = new Real [NUMBER_OF_PARAMS];
		deriv_phi_ = new Real [NUMBER_OF_PARAMS];
		deriv_psi_ = new Real [NUMBER_OF_PARAMS];
	}
	CentroidRotamerSampleData(Real p, Real d, Real a, Real w, Real vd, Real va, Real vw)
	:prob_(p),distance_(d),angle_(a),dihedral_(w), sd_dis_(vd), sd_ang_(va), sd_dih_(vw) {
		energy_ = -log(p);
		norm_factor_ = 1.0; ///sqrt(2.0*3.14159265) no need
		data_ = new Real [NUMBER_OF_PARAMS];
		deriv_phi_ = new Real [NUMBER_OF_PARAMS];
		deriv_psi_ = new Real [NUMBER_OF_PARAMS];
	}

	CentroidRotamerSampleData( const CentroidRotamerSampleData &cs )
	{
		prob_		=cs.prob();
		energy_ 	=cs.energy();
		distance_	=cs.distance();
		angle_		=cs.angle();
		dihedral_	=cs.dihedral();
		sd_dis_	=cs.sd_dis();
		sd_ang_	=cs.sd_ang();
		sd_dih_	=cs.sd_dih();
		norm_factor_ = cs.norm_factor();

		data_ = new Real [NUMBER_OF_PARAMS];
		deriv_phi_ = new Real [NUMBER_OF_PARAMS];
		deriv_psi_ = new Real [NUMBER_OF_PARAMS];

		for (Size i=0; i<NUMBER_OF_PARAMS; i++) {
			data_[i] = cs.data_[i];
			deriv_phi_[i] = cs.deriv_phi_[i];
			deriv_psi_[i] = cs.deriv_psi_[i];
		}
	}

	~CentroidRotamerSampleData(){
		delete [] data_;
		delete [] deriv_phi_;
		delete [] deriv_psi_;
	}

	Real distance() const { return distance_; }
	Real angle() const { return angle_; }
	Real dihedral() const { return dihedral_; }
	Real sd_dis() const { return sd_dis_; }
	Real sd_ang() const { return sd_ang_; }
	Real sd_dih() const { return sd_dih_; }
	Real prob() const { return prob_; }
	Real energy() const {return energy_; }
	Real norm_factor() const {
		return 1.0;
		//return norm_factor_;
		//return 15.74961; //sqrt(2pi)^3
	}

	void set_distance( Real d ){ distance_ = d; }
	void set_angle( Real a ){ angle_ = a; }
	void set_dihedral( Real w ){ dihedral_ = w; }
	void set_prob( Real p ){ prob_ = p; energy_=-log(prob_); }
	void set_sd_dis( Real s ){ sd_dis_ = s; }
	void set_sd_ang( Real s ){ sd_ang_ = s; }
	void set_sd_dih( Real s ){ sd_dih_ = s; }

	void private_data_to_public_array(); //init before interp
	void public_array_to_private_data(); //after interp

	/// DOF3 sample: (dis, ange, dih)
	/// calculate the distance between this rot and given CEN
	Real cal_distance( const DOF3 &sample, bool use_xyz=false ) const;
	Real cal_distance_squared( const DOF3 &sample, bool use_xyz=false ) const;
	Real cal_distance_squared( const conformation::Residue & rsd ) const;

	/// return the value of angle (in rad)
	Real cal_delta_internal_coordinates_squared(
		const conformation::Residue & rsd,
		Real & d_sq, Real & a_sq, Real & w_sq ) const;
	Real cal_delta_internal_coordinates(
		const conformation::Residue & rsd,
		Real & delta_d, Real & delta_a, Real & delta_w ) const;
	
	/// generate a random rot inside the well
	void assign_random_rotamer( DOF3 &sample, numeric::random::RandomGenerator &RG ) const;
	/// generate the best rot (mean of the well)
	void assign_best_rotamer( DOF3 &sample ) const;
};


//////////////////////////////////////////////////////////////////////////////
class SingleResidueCenrotLibrary : public core::pack::rotamers::SingleResidueRotamerLibrary {

public:
typedef chemical::AA AA;

public:
SingleResidueCenrotLibrary(AA const aa);
virtual ~SingleResidueCenrotLibrary();

public:
AA aa() const { return aa_; }
std::string read_from_file(
	utility::io::izstream & infile,
	bool first_line_three_letter_code_already_read );

const utility::vector1< CentroidRotamerSampleData >
get_rotamer_samples(conformation::Residue const & rsd) const;

/// Virtual functions required by the base classes
virtual
Real rotamer_energy(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const;

virtual
Real rotamer_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const;

//eval cart version
Real eval_rotameric_energy_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch,
	bool eval_deriv
) const;

//eval internal version (bb only)
Real eval_rotameric_energy_bb_dof_deriv(
	conformation::Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const;

/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
/// (based on e.g. its current phi and psi values).
/// If curr_rotamer_only is true, then consider only the idealized version of the
/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
virtual
Real best_rotamer_energy(
	conformation::Residue const & rsd,
	bool curr_rotamer_only,
	RotamerLibraryScratchSpace & scratch
) const;

virtual void
assign_random_rotamer_with_bias(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	RotamerLibraryScratchSpace & scratch,
	numeric::random::RandomGenerator & RG,
	ChiVector & new_chi_angles,
	bool perturb_from_rotamer_center
) const;

virtual void
fill_rotamer_vector(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue,
	utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
	bool buried,
	rotamers::RotamerVector & rotamers
) const;

virtual void write_to_file( utility::io::ozstream &out ) const;

CentroidRotamerSampleData const & get_closest_rotamer(
	conformation::Residue const & rsd, Size &nrot, Real &dis) const;

private:
void get_phipsi_bins(
	Real phi,
	Real psi,
	Size & phibin,
	Size & psibin,
	Size & phibin_next,
	Size & psibin_next,
	Real & phi_alpha,
	Real & psi_alpha
) const;

void get_phipsi_bins(
	Real phi,
	Real psi,
	Size & phibin,
	Size & psibin ) const;

Real get_phi_from_rsd(
	conformation::Residue const & rsd
) const;

Real get_psi_from_rsd(
	conformation::Residue const & rsd
) const;

void verify_phipsi_bins(
	Real phi,
	Real psi,
	Size const phibin,
	Size const psibin,
	Size const phibin_next,
	Size const psibin_next
) const;

/// @brief This is not the right place for this code, but the numeric interpolation library
/// uselessly indexes by 0 and the basic functions aren't inlined...
inline
void bin_angle(
	Real const angle_start,
	Real const angle_step,
	Real const ASSERT_ONLY( angle_range ),
	Size const nbins,
	Real const ang,
	Size & bin_lower,
	Size & bin_upper,
	Real & angle_alpha
) const {
	/// very, very rarely, periodic_range( angle, 360 ) will return 180 instead of -180.
	/// though it is supposed to return values in the range [-180, 180).
debug_assert( angle_start <= ang && ang <= angle_start + angle_range );
debug_assert( std::abs( nbins * angle_step - angle_range ) < 1e-15 );

	Real real_bin_lower = ( ang - angle_start ) / angle_step;
	Size bin_prev = static_cast< Size > ( real_bin_lower );
	bin_lower = 1 + numeric::mod( bin_prev, nbins );
	bin_upper = numeric::mod( bin_lower, nbins ) + 1;
	angle_alpha = ( (ang - angle_start ) - ( bin_prev * angle_step ) ) / angle_step;
}

protected:
	static Size const N_PHIPSI_BINS;
	static Real const PHIPSI_BINRANGE;
	static Size const RSD_PHI_INDEX;
	static Size const RSD_PSI_INDEX;

	//Real const SingleResidueDunbrackLibrary::NEUTRAL_PHI = -90;
	// R++ value.  Roland Dunbrack suggests -60.
	//Real const SingleResidueDunbrackLibrary::NEUTRAL_PSI = 130;
	// R++ value.  Roland Dunbrack suggests  60.

	static Real const NEUTRAL_PHI;
	static Real const NEUTRAL_PSI;

	static Real const MAX_ROT_ENERGY;
	static Real const MIN_ROT_PROB;

private:
	//utility::vector1< CentroidRotamerSampleData > all_rots_;
	ObjexxFCL::FArray2D< utility::vector1< CentroidRotamerSampleData > > all_rots_bb_;
	utility::vector1< CentroidRotamerSampleData > dummy_sample_;
	AA aa_;
	Size max_rot_num;
	Real ref_energy_; //ref energy for rot, -log(1/N)
	// ref_energy_ is too rude
	//pointed by hahnbeom, the ref energy should be treated as bb dependent
	//E_ref(phi, psi) = <log Pi(phi, psi))> = SUM Pi*logPi / SUM Pi
	ObjexxFCL::FArray2D< Real > entropy_; /// E_ref is the same as shannon entropy

	void setup_entropy_correction();
};

} // cenrot
} // dunbrack
} // pack
} // core

#endif // INCLUDED_core_pack_dunbrack_cenrot_SingleResidueCenrotLibrary_HH

