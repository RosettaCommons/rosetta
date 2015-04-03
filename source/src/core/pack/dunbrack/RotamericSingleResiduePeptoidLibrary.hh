// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/dunbrack/RotamericSingleResiduePeptoidLibrary.hh
/// @brief   Declaration of rotameric peptoid rotamer libraries
/// @brief   Similar to RotamericSingleResidueDunbrackLibrary
/// @author  P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResiduePeptoidLibrary_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResiduePeptoidLibrary_hh

// Unit Headers
#include <core/pack/dunbrack/RotamericSingleResiduePeptoidLibrary.fwd.hh>

// Package Headers
#include <core/pack/dunbrack/SingleResiduePeptoidLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray4D.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>

namespace core {
namespace pack {
namespace dunbrack {

template < Size T, Size N >
class RotamericSingleResiduePeptoidLibrary : public SingleResiduePeptoidLibrary
{
public:
	typedef SingleResiduePeptoidLibrary parent;

private:
	/// No default ctor
	///RotamericSingleResiduePeptoidLibrary();

public:
	RotamericSingleResiduePeptoidLibrary();

	virtual ~RotamericSingleResiduePeptoidLibrary();

	friend class SingleResiduePeptoidLibrary;

public:
	/// Virtual functions required by the base classes

	virtual
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	virtual
	Real
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
	/// (based on e.g. its current phi and psi values).
	/// If curr_rotamer_only is true, then consider only the idealized version of the
	/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
	virtual
	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		bool curr_rotamer_only,
		RotamerLibraryScratchSpace & scratch
	) const;

	virtual
	void
	assign_random_rotamer_with_bias(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const;

	virtual
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		RotamerVector & rotamers
	) const;

	/// DOUG DOUG DOUG The definitions of these need to change (and have been changed)
	/// @brief Return all of the rotamer sample data given a particular phi/psi.
	/// For N-terminus residues, hand in the phi value SingleResiduePeptoidLibrary::PHI_NEUTRAL and
	/// for C-terminus residues, hand in the psi value SingleResiduePeptoidLibrary::PSI_NEUTRAL.
	/// The returned samples should be in semi-decrasing order by probability; semi, because the
	/// rotamers are constructed in sorted order by their probability in the lower phi-psi bin that
	/// the input phi/psi perscribes.
	virtual
	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples(
		Real omg,
		Real phi,
		Real psi
	) const;

	virtual
	Real
	get_probability_for_rotamer(
		Real omg,
		Real phi,
		Real psi,
		Size rot_ind
	) const;

	virtual
	DunbrackRotamerSampleData
	get_rotamer(
		Real omg,
		Real phi,
		Real psi,
		Size rot_ind
	) const;

	virtual
	Size nchi() const;

	virtual
	Size n_rotamer_bins() const;

	virtual
	void
	write_to_file( utility::io::ozstream &out ) const;

	virtual void write_to_binary( utility::io::ozstream & out ) const;
	virtual void read_from_binary( utility::io::izstream & in );

	virtual
	void
	get_rotamer_from_chi(
		ChiVector const & chi,
		RotVector & rot
	) const;

protected:

	/// @brief When given a statically sized fixedsizearray, use this method.
	void
	get_rotamer_from_chi_static(
		ChiVector const & chi,
		Size4 & rot
	) const;

	/// @brief When given a statically sized fixedsizearray, use this method.
	void
	get_rotamer_from_chi_static(
		Real4 const & chi,
		Size4 & rot
	) const;


public:

	/// @brief Read the rot lib from the input stream
	void
	read_from_file(
		utility::io::izstream & in
	);

protected:
	/// DOUG DOUG DOUG This use to be N_PHIPSI_BINS and PHIPSI_BINRANGE

	static Size const N_OMG_BINS = 7;
	static Size const N_PHI_BINS = 36;
	static Size const N_PSI_BINS = 36;
	static Real const OMG_BINRANGE;
	static Real const PHI_BINRANGE;
	static Real const PSI_BINRANGE;
	static Real const CIS_OMG_LOWER_RANGE;
	static Real const CIS_OMG_UPPER_RANGE;
	static Real const TRANS_OMG_LOWER_RANGE;
	static Real const TRANS_OMG_UPPER_RANGE;
	static Real const CIS_UPPER_TRANS_LOWER_SPLIT;
	static Real const TRANS_UPPER_CIS_LOWER_SPLIT;

protected:
	/// Read and write access for derived classes

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T, N > > const &
	rotamers( Real omg_angle ) const;

	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T, N > > &
	rotamers( Real omg_angle );

	ObjexxFCL::FArray4D< Size > const &
	packed_rotno_2_sorted_rotno( Real omg_angle ) const;

	ObjexxFCL::FArray4D< Size > &
	packed_rotno_2_sorted_rotno( Real omg_angle );

protected:
	/// Worker functions

	virtual Size memory_usage_static() const;
	virtual Size memory_usage_dynamic() const;

	/// @brief  Evaluates the score and chi-deviation penalty for the rotameric
	/// chi (in this class, that means all the chi) and stores the answers in
	/// the scratch object.  If eval_deriv is true, then at the end of this
	/// function, scratch contains up-to-date dchidevpen_dbb, dchidevpen_dchi,
	/// chimean, chisd, chidev, chidevpen, dchimean_d(phi/psi), dchisd_d(phi/psi)
	/// rotwell and rotprob data.
	Real
	eval_rotameric_energy_deriv(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		bool eval_deriv
	) const;

	Size
	find_another_representative_for_unlikely_rotamer(
		conformation::Residue const & rsd,
		Size4 & rotwell
	) const;


	void
	interpolate_rotamers(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch,
		Size const packed_rotno,
		PackedDunbrackRotamer< T, N, Real > & interpolated_rotamer
	) const;

	void
	interpolate_rotamers(
		typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T, N > > const & rotamers,
		ObjexxFCL::FArray4D< Size > const & packed_rotno_2_sorted_rotno,
		RotamerLibraryScratchSpace & scratch,
		Size const packed_rotno,
		Size const omgbin,
		Size const phibin,
		Size const psibin,
		Size const omgbin_next,
		Size const phibin_next,
		Size const psibin_next,
		Real const omg_alpha,
		Real const phi_alpha,
		Real const psi_alpha,
		PackedDunbrackRotamer< T, N, Real > & interpolated_rotamer
	) const;

	/// @brief Assigns random chi angles and returns the packed_rotno for the chosen random rotamer.
	void
	assign_random_rotamer(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center,
		Size & packed_rotno
	) const;

	void
	assign_chi_for_interpolated_rotamer(
		PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer,
		conformation::Residue const & rsd,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const;

	void
	correct_termini_derivatives(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	/// DOUG DOUG DOUG Make public just for now
public:
	void
	get_omgphipsi_bins(
		Real omg,
		Real phi,
		Real psi,
		Size & omgbin,
		Size & phibin,
		Size & psibin,
		Size & omgbin_next,
		Size & phibin_next,
		Size & psibin_next,
		Real & omg_alpha,
		Real & phi_alpha,
		Real & psi_alpha
	) const;

	void
	get_omgphipsi_bins(
		Real omg,
		Real phi,
		Real psi,
		Size & omgbin,
		Size & phibin,
		Size & psibin
	) const;

	/// DOUG DOUG DOUG
protected:
	Real
	get_omg_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const;

	Real
	get_phi_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const;

	Real
	get_psi_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const;


private:
	/// @brief After an interpolated rotamer has been created, this method appends new
	/// rotamers to the vector of rotamers, taking extra samples at designated chi intervals
	/// if instructed to do so in the extra_chi_steps argument.
	void
	build_rotamers(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		RotamerVector & rotamers,
		PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer
	) const;


	void verify_omgphipsi_bins(
		Real omg,
		Real phi,
		Real psi,
		Size const omgbin,
		Size const phibin,
		Size const psibin,
		Size const omgbin_next,
		Size const phibin_next,
		Size const psibin_next
	) const;

protected:

	template< class P >
	DunbrackRotamer< T, N, P >
	packed_rotamer_2_regular_rotamer(
		PackedDunbrackRotamer< T, N, P > const & packedrot
	) const;

	/// @brief This member function constructs a list of all combinations of chi angles
	/// for a rotamer sample.  It relies on a virtual function chisamples_for_rotamer_chi
	/// that may be overridden by derived classes.  With a list of samples for each chi,
	/// this function then enumerates all combinations.
	void
	enumerate_chi_sets(
		chemical::ResidueType const & rsd_type,
		pack::task::PackerTask const & task,
		Size const seqpos,
		bool buried,
		RotamericData< T, N > const & rotamer_data,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		utility::vector1< ChiSetOP > & chi_set_vector
	) const;

	/// @brief Used in tandem with enumerate_chi_sets, this function pushes
	/// back chi sample data into the four input vectors: the 1. chi value sample,
	/// 2. the rotamer well #, 3. a discription of the kind of sample
	/// (how far off the ideal rotamer) and 4. the probability of observing this
	/// chi at that sample.  The rotamer_data that's input is the same data that
	/// was handed to the enumerate_chi_sets call (may be downcast as needed).
	virtual
	void
	chisamples_for_rotamer_and_chi(
		chemical::ResidueType const & rsd_type,
		pack::task::ResidueLevelTask const & rtask,
		bool buried,
		Size const chi_index,
		RotamericData< T, N > const & rotamer_data,
		utility::vector1< Real > const & extra_steps,
		utility::vector1< Real > & total_chi,
		utility::vector1< int  > & total_rot,
		utility::vector1< Real > & total_ex_steps,
		utility::vector1< Real > & chisample_prob
	) const;

	/// @brief Once all the chi have been enumerated, building the rotamers is a trivial task.
	/// This function is protected so that derived classes may simply enumerate their chi and then
	/// invoke this function.
	///  This arguably should be moved into the SingleResidueRotamerLibrary base class.
	void
	create_rotamers_from_chisets(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< ChiSetOP > const & chi_set_vector,
		RotamerVector & rotamers
	) const;

	// /// @brief returns if the omega angle is in the range defined by the cis omega part of the rotlib
	// bool
	// is_trans_omg(
	// 	Real omg_angle
	// ) const
	// {
	// 	Real angle( numeric::principal_angle_degrees( omg_angle ) ); // returns angle ( -180.0, 180.0 ]
	// 	if ( angle >= TRANS_OMG_LOWER_RANGE && angle <= TRANS_OMG_UPPER_RANGE ) return true;
	// 	return false;
	// }

	// /// @brief returns true if the omega angle is in the range defined by the trans omega part of the rotlib
	// bool
	// is_cis_omg(
	// 	Real omg_angle
	// ) const
	// {
	// 	Real angle( numeric::principal_angle_degrees( omg_angle ) ); // returns angle ( -180.0, 180.0 ]
	// 	if ( angle >= CIS_OMG_LOWER_RANGE && angle <= CIS_OMG_UPPER_RANGE) return true;
	// 	return false;
	// }

	/// @brief returns if the omega angle is in the range defined by the cis omega part of the rotlib
	bool
	is_trans_omg(
		Real omg_angle
	) const
	{
		Real angle( numeric::nonnegative_principal_angle_degrees( omg_angle ) );
		if ( angle >= 90.00 && angle < 270.00 ) return true;
		return false;
	}

	/// @brief returns true if the omega angle is in the range defined by the trans omega part of the rotlib
	bool
	is_cis_omg(
		Real omg_angle
	) const
	{
		Real angle( numeric::nonnegative_principal_angle_degrees( omg_angle ) );
		if ( angle >= 270 || angle < 90 ) return true;
		return false;
	}

private:

	/// The (chi_mean, chi_sd, packed_rotno, and prob) data for the chi dihedrals
	/// The FArray4D is indexed into by (omg, phi, psi, sorted_index ), where
	/// sorted index simply means the order for a particular packed_rotno in the
	/// list of rotamers sorted by probability.
	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T, N > > trans_rotamers_;
	typename ObjexxFCL::FArray4D< PackedDunbrackRotamer< T, N > > cis_rotamers_;

	/// Quick lookup that lists the sorted position for the packed rotamer number
	/// given a omg/phi/psi.  Indexed by (omg, phi, psi, packed_rotno ).
	ObjexxFCL::FArray4D< Size > trans_packed_rotno_2_sorted_rotno_;
	ObjexxFCL::FArray4D< Size > cis_packed_rotno_2_sorted_rotno_;
};


} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_RotamericSingleResiduePeptoidLibrary_HH
