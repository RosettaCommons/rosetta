// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/RotamericSingleResidueDunbrackLibrary.hh
/// @brief   Declaration of rotameric rotamer libraries from Jun08 (as opposed to semi-rotameric)
/// @author  Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibrary_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibrary_hh

// Unit Headers
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.fwd.hh>

// Package Headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
//#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

template < Size T, Size N >
class RotamericSingleResidueDunbrackLibrary : public SingleResidueDunbrackLibrary
{

	friend class RotamericSingleResidueDunbrackLibraryParser;

public:
	typedef SingleResidueDunbrackLibrary parent;

	// PyRosetta debug helper
	//std::string who_am_i() override { return "RotamericSingleResidueDunbrackLibrary<" + std::to_string(T) + ", " + std::to_string(N) + ">"; }


private:
	/// No default ctor
	///RotamericSingleResidueDunbrackLibrary();

public:

	/// @brief Creator.
	RotamericSingleResidueDunbrackLibrary(
		chemical::ResidueType const & rt, // mostly use for aa.
		bool dun02,
		bool use_bicubic,
		bool entropy_correction,
		core::Real prob_buried, // 0.98
		core::Real prob_nonburied, // 0.95
		bool const reduced_resolution_library=false
	);

	virtual ~RotamericSingleResidueDunbrackLibrary() throw();

	friend class SingleResidueDunbrackLibrary;

public:
	/// Virtual functions required by the base classes

	virtual
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch
	) const override;

	virtual
	Real
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch
	) const override;

	/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
	/// (based on e.g. its current phi and psi values).
	/// If curr_rotamer_only is true, then consider only the idealized version of the
	/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
	virtual
	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		bool curr_rotamer_only,
		RotamerLibraryScratchSpace & scratch
	) const override;

	virtual
	void
	assign_random_rotamer_with_bias(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const override;

	virtual
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		utility::graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		rotamers::RotamerVector & rotamers
	) const override;

	/// @brief Return all of the rotamer sample data given a particular backbone bin.
	/// For N-terminus residues, hand in the first bb value SingleResidueDunbrackLibrary::NEUTRAL_PHI and
	/// for C-terminus residues, hand in the final bb value SingleResidueDunbrackLibrary::NEUTRAL_PSI.
	/// The returned samples should be in semi-decrasing order by probability; semi, because the
	/// rotamers are constructed in sorted order by their probability in the lower torsion bin that
	/// the input backbone torsions prescribe.
	virtual
	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples(
		Real5 bbs
	) const override;

	virtual
	Real
	get_probability_for_rotamer(
		utility::fixedsizearray1< Real, N > bbs,
		Size rot_ind
	) const;

	virtual
	DunbrackRotamerSampleData
	get_rotamer(
		utility::fixedsizearray1< Real, N > bbs,
		Size rot_ind
	) const;

	virtual
	Real
	get_probability_for_rotamer(
		Real phi,
		Real psi,
		Size rot_ind
	) const override;

	virtual
	DunbrackRotamerSampleData
	get_rotamer(
		Real phi,
		Real psi,
		Size rot_ind
	) const override;

	virtual
	Size nchi() const override;

	virtual
	Size nbb() const override;

	virtual
	Size n_rotamer_bins() const override;

	virtual
	void
	write_to_file( utility::io::ozstream &out ) const override;

	virtual void write_to_binary( utility::io::ozstream & out ) const override;
	virtual void read_from_binary( utility::io::izstream & in ) override;

	/// @brief Comparison operator, mainly intended to use in ASCII/binary comparsion tests
	/// Values tested should parallel those used in the read_from_binary() function.
	virtual
	bool
	operator ==( SingleResidueRotamerLibrary const & ) const override;

	virtual
	void
	get_rotamer_from_chi(
		ChiVector const & chi,
		RotVector & rot
	) const override;

protected:

	void
	initialize_bicubic_splines();

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

	/// @brief Read from input stream; stream may contain data for other amino acids.
	/// Quit once another amino acid is specified in the input file, returning the
	/// name of the next amino acid specifed (since it's already been extracted from
	/// the input stream).  Return the empty string if no other amino acid is specified.
	/// @details Returns the three letter string of the next amino acid specified in the
	/// input library.
	/// @author Rewritten on 8 January 2018 by Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	read_from_file(
		utility::io::izstream & in,
		bool first_line_three_letter_code_already_read
	);

public:

	utility::fixedsizearray1< Size, N > N_BB_BINS;
	utility::fixedsizearray1< Real, N > BB_BINRANGE;

	utility::fixedsizearray1< utility::vector1< core::Size >, N> BIN_EQUIVALENTS;

	// Public by convention!
	// This array describes function objects that, given a residue and some context, return some
	// computable value on which this rotamer library depends.
	// Typical rotamer libraries depend on backbone dihedrals, but not all!
	// Note: would use lambdas -- minutely lighter-weight, better syntax for types -- but we need
	// C++14 for the full deduction.
	utility::fixedsizearray1< std::function< Real( conformation::Residue const & rsd, pose::Pose const & pose ) >, N > IVs;

protected:
	/// Read and write access for derived classes and parser class

	typename ObjexxFCL::FArray2D< PackedDunbrackRotamer< T, N > > const &
	rotamers() const {
		return rotamers_;
	}

	typename ObjexxFCL::FArray2D< PackedDunbrackRotamer< T, N > > &
	rotamers() {
		return rotamers_;
	}

	ObjexxFCL::FArray2D< Size > const &
	packed_rotno_2_sorted_rotno() const {
		return packed_rotno_2_sorted_rotno_;
	}

	ObjexxFCL::FArray2D< Size > &
	packed_rotno_2_sorted_rotno() {
		return packed_rotno_2_sorted_rotno_;
	}

protected:
	/// Worker functions

	virtual Size memory_usage_static() const override;
	virtual Size memory_usage_dynamic() const override;

	/// @brief  Evaluates the score and chi-deviation penalty for the rotameric
	/// chi (in this class, that means all the chi) and stores the answers in
	/// the scratch object.  If eval_deriv is true, then at the end of this
	/// function, scratch contains up-to-date dchidevpen_dbb, dchidevpen_dchi,
	/// chimean, chisd, chidev, chidevpen, dchimean_d(phi/psi), dchisd_d(phi/psi)
	/// rotwell and rotprob data.
	Real
	eval_rotameric_energy_deriv(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch,
		bool eval_deriv
	) const;

public:
	Size
	find_another_representative_for_unlikely_rotamer(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Size4 & rotwell
	) const;


	void
	interpolate_rotamers(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		RotamerLibraryScratchSpace & scratch,
		Size packed_rotno,
		PackedDunbrackRotamer< T, N, Real > & interpolated_rotamer
	) const;

protected:

	void
	interpolate_rotamers(
		RotamerLibraryScratchSpace & scratch,
		Size packed_rotno,
		utility::fixedsizearray1< Size, N > const & bb_bin,
		utility::fixedsizearray1< Size, N > const & bb_bin_next,
		utility::fixedsizearray1< Real, N > const & bb_alpha,
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

	Real
	get_phi_from_rsd(
		conformation::Residue const & rsd
	) const override;

	Real
	get_psi_from_rsd(
		conformation::Residue const & rsd
	) const override;


public:

	void
	get_bb_bins(
		utility::fixedsizearray1< Real, N > const & bbs,
		utility::fixedsizearray1< Size, N > & bb_bin,
		utility::fixedsizearray1< Size, N > & bb_bin_next,
		utility::fixedsizearray1< Real, N > & bb_alpha
	) const;

	void
	get_bb_bins(
		utility::fixedsizearray1< Real, N > const & bbs,
		utility::fixedsizearray1< Size, N > & bb_bin
	) const;

	/// @brief Version for external, non-template classes to access.
	/// @details Converts to fixedsizearrays internally.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void
	get_bb_bins(
		utility::vector1< core::Real > const & bbs,
		utility::vector1< core::Size > & bb_bin
	) const;

protected:

	utility::fixedsizearray1< Real, N >
	get_IVs_from_rsd(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const;

	Real
	get_IV_from_rsd(
		core::Size bbn,
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
		utility::graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		rotamers::RotamerVector & rotamers,
		PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer
	) const;


	void verify_phipsi_bins(
		Real phi,
		Real psi,
		Size const phibin,
		Size const psibin,
		Size const phibin_next,
		Size const psibin_next
	) const;

	void verify_bb_bins(
		utility::fixedsizearray1< Real, N > const & bbs,
		utility::fixedsizearray1< Size, N > const & bb_bin,
		utility::fixedsizearray1< Size, N > const & bb_bin_next
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
		utility::graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< ChiSetOP > const & chi_set_vector,
		rotamers::RotamerVector & rotamers
	) const;

	void setup_entropy_correction();
	void setup_entropy_correction() const;

private:

	/// The (chi_mean, chi_sd, packed_rotno, and prob) data for the chi dihedrals
	/// The FArray3D is indexed into by (bb_bin_index, sorted_index ), where
	/// sorted index simply means the order for a particular packed_rotno in the
	/// list of rotamers sorted by probability and the bb_bin_index is a composite
	/// of what you would get from essentially expressing the backbone torsions as
	/// a number in base N_BB_BINS (often 36).
	typename ObjexxFCL::FArray2D< PackedDunbrackRotamer< T, N > > rotamers_;
	/// Quick lookup that lists the sorted position for the packed rotamer number
	/// given a phi/psi.  Indexed by (bb_bin_index, packed_rotno ).
	ObjexxFCL::FArray2D< Size > packed_rotno_2_sorted_rotno_;

	// Entropy correction
	utility::fixedsizearray1< ObjexxFCL::FArray1D< Real >, ( 1 << N ) > ShannonEntropy_n_derivs_;

	// AMW: one peptoid flag needed
	bool peptoid_ = false;

};


} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibrary_HH
