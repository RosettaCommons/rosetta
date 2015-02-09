// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh
/// @brief   Declaration of semi-rotameric rotamer libraries from Jun08
/// @author  Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_hh
#define INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_hh

// Unit Headers
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>

// Package Headers
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
//#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

template < class P = DunbrackReal >
class BBDepNRChiSample;

template < class P = DunbrackReal >
class BBIndNRChiSample;

/// P for precision
template < class P >
class BBDepNRChiSample
{
public:
	BBDepNRChiSample() :
		packed_rotno_( 0 ),
		nrchi_bin_( 0 ),
		nrchi_mean_( 0.0 ),
		nrchi_sd_( 0.0 ),
		prob_( 0.0 )
	{}

	bool
	operator==( BBDepNRChiSample<P> const & other ) const;

	/// The packed_rotno_ information is sufficient to determine the
	/// mean and standard deviation of the rotameric chis; it is used
	/// as an index into another table containing this data.
	Size packed_rotno_;
	Size nrchi_bin_;
	P nrchi_mean_;
	P nrchi_sd_;
	P prob_;
};

template < class P >
class BBIndNRChiSample
{
public:
	BBIndNRChiSample() :
		left_(   P( 0.0 ) ),
		median_( P( 0.0 ) ),
		right_(  P( 0.0 ) ),
		prob_(   P( 0.0 ) )
	{}

	bool
	operator==( BBIndNRChiSample<P> const & other ) const;

	P left_;
	P median_;
	P right_;
	P prob_;
};

/// @brief A class to hold rotamer building data on the stack and yet have it
/// accessible to derived classes when invoking base class functions.  An
/// alternative would have been to store mutable member data in the Library
/// class itself. This option, however, is not thread safe.
/// This data is used by the SemiRotamericSRDL class for when building
/// backbone dependent rotamers.
template < Size T, Size N >
class BBDepSemiRotamericData : public RotamericData< T, N >
{
public:
	typedef RotamericData< T, N > parent;

public:
	BBDepSemiRotamericData(
		DunbrackRotamer< T, N, Real > const & rotamer_in,
		BBDepNRChiSample< Real > const & bbdep_nrchi_sample_in
	)
	:
		parent( rotamer_in ),
		bbdep_nrchi_sample_( bbdep_nrchi_sample_in )
	{}

	virtual ~BBDepSemiRotamericData() {}

	BBDepNRChiSample< Real > const &
	bbdep_nrchi_sample() const {
		return bbdep_nrchi_sample_;
	}

private:
	BBDepNRChiSample< Real > bbdep_nrchi_sample_;

};

template < Size N >
struct BBDepScoreInterpData
{
public:
    BBDepScoreInterpData()
    {}
	
	inline bool operator==( BBDepScoreInterpData< N > const & other ) const;
    
	utility::fixedsizearray1< Real, ( 1 << (N+1) ) > n_derivs_;

};


/// @brief A class to hold rotamer building data on the stack and yet have it
/// accessible to derived classes when invoking base class functions.  An
/// alternative would have been to store mutable member data in the Library
/// class itself. This option, however, is not thread safe.
/// This data is used by the SemiRotamericSRDL class for when building
/// backbone independent rotamers.
template < Size T, Size N >
class BBIndSemiRotamericData : public RotamericData< T, N >
{
public:
	typedef RotamericData< T, N > parent;
public:
	BBIndSemiRotamericData(
		DunbrackRotamer< T, N, Real > const & rotamer,
		BBIndNRChiSample<> const & bbind_nrchi_sample_in,
		Size const nrchi_bin_id_in
	)
	:
		parent( rotamer ),
		bbind_nrchi_sample_( bbind_nrchi_sample_in ),
		nrchi_bin_id_( nrchi_bin_id_in )
	{}

	virtual ~BBIndSemiRotamericData() {}

	BBIndNRChiSample<> const &
	bbind_nrchi_sample() const {
		return bbind_nrchi_sample_;
	}

	Size
	nrchi_bin_id() const {
		return nrchi_bin_id_;
	}

private:

	BBIndNRChiSample<> bbind_nrchi_sample_;
	Size nrchi_bin_id_;
};

class ProbSortClass {
public:
	Real probability_;
	Size index_;
};

bool
psc_compare( ProbSortClass left, ProbSortClass right );


/// @brief  This class is meant to represent the non-rotameric chi observed in
/// several amino acids (asn, asp, gln, glu, his, phe, trp, tyr ) which are rotameric
/// for the chi closest to the backbone and non rotameric for exactly one chi angle.
/// This non-rotameric chi (abv. nrchi) is the last chi for each of these 8 amino acids
/// except tyrosine, where this chi is the last heavy-atom chi.  The last chi on tyrosine
/// governs a hydroxyl.
/// Unlike in the fully rotameric residues, the last heavyatom chi in semi-rotameric residues
/// do not "communicate" to the rotameric chi.  That is, in the rotameric chi, the mean chi1 value
/// is sensitive to the chi3 value.  If the third diherdal switches from trans to g+, then chi1
/// would shift in response.  Changes to the non-rotameric chi do not effect the rotameric chi.
/// The data structure here is good for this model but no other.
template < Size T, Size N >
class SemiRotamericSingleResidueDunbrackLibrary : public RotamericSingleResidueDunbrackLibrary< T, N >
{
public:
	typedef RotamericSingleResidueDunbrackLibrary< T, N > parent;
	typedef SingleResidueDunbrackLibrary grandparent;

public:

	/// @brief The constructor determines the path the library takes: whether
	/// it uses a backbone dependent or independent score function for the
	/// non-rotameric chi (the rotameric chi score function is always backbone dependent)
	/// and whether it uses a backbone dependent or independent rotamer sampling
	/// scheme.  All four combinations are possible, though backbone independent scoring
	/// and backbone dependent rotamer sampling seems a poor combination.
	SemiRotamericSingleResidueDunbrackLibrary(
		chemical::AA const aa_in,
		bool const backbone_independent_scoring,         // true uses less memory
		bool const backbone_independent_rotamer_sampling // true uses less memory
	);

	virtual ~SemiRotamericSingleResidueDunbrackLibrary() throw();

	friend class SingleResidueDunbrackLibrary;
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

	/// @brief Return all of the rotamer sample data given a particular phi/psi.
	/// For N-terminus residues, hand in the phi value SingleResidueDunbrackLibrary::PHI_NEUTRAL and
	/// for C-terminus residues, hand in the psi value SingleResidueDunbrackLibrary::PSI_NEUTRAL.
	/// The returned samples should be in semi-decrasing order by probability; semi, because the
	/// rotamers are constructed in sorted order by their probability in the lower phi-psi bin that
	/// the input phi/psi perscribes.
	virtual
	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples(
		utility::fixedsizearray1< Real, N > bbs
	) const;

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
	Size nchi() const;
	
	virtual
	Size nbb() const;

	virtual
	Size n_rotamer_bins() const;

	//XRW_B_T1
	/*
	// stubbed out
	virtual
	SingleResidueRotamerLibraryOP
	coarsify(coarse::Translator const &map) const;
	*/
	//XRW_E_T1

	// stubbed out
	virtual
	void
	write_to_file( utility::io::ozstream & out ) const;

	virtual void write_to_binary( utility::io::ozstream & out ) const;

	/// @brief Initialize either a backbone-independent or a backbone-dependent SRSRDL
	/// from the set of four files which describe both (not all files are read).
	void read_from_files(
		utility::io::izstream & in_rotdef,
		utility::io::izstream & in_rotameric,
		utility::io::izstream & in_continmin_bbdep,
		utility::io::izstream & in_continmin_bbind
	);

	/// @brief Initialize a backbone-dependent SRSRDL from the set of three files
	/// which describe it.
	void
	read_from_files(
		utility::io::izstream & in_rotdef,
		utility::io::izstream & in_rotameric,
		utility::io::izstream & in_continmin_bbdep
	);

	virtual void read_from_binary( utility::io::izstream & in );

	/// @brief Comparison operator, mainly intended to use in ASCII/binary comparsion tests
	/// Values tested should parallel those used in the read_from_binary() function.
	virtual
	bool
	operator ==( SingleResidueRotamerLibrary const & ) const;

	virtual
	void
	get_rotamer_from_chi(
		ChiVector const & chi,
		RotVector & rot
	) const;

public:
	/// Initialization functions that must be called before reading the input libraries.


	/// @brief For the non rotameric chi, how many asymmetric degrees are there?  (e.g. 360 for asn, 180 for asp)
	void
	set_nrchi_periodicity( Real angle_in_degrees );

	/// @brief What angle do the interpolation data start from?
	void
	set_nonrotameric_chi_start_angle( Real angle_in_degrees );

	/// @brief What is the angular step size of the bbdep score?
	void
	set_nonrotameric_chi_bbdep_scoring_step_size( Real step_size_in_degrees );

	/// @brief What is the angular step size of the bbind score?
	void
	set_nonrotameric_chi_bbind_scoring_step_size( Real step_size_in_degrees );

protected:
    
	virtual Size memory_usage_static() const;
	virtual Size memory_usage_dynamic() const;

	Real
	rotamer_energy_deriv_bbdep(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		bool eval_deriv
	) const;

	Real
	rotamer_energy_deriv_bbind(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		bool eval_deriv
	) const;

	void
	assign_random_rotamer_with_bias_bbind(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const;

	void
	assign_random_rotamer_with_bias_bbdep(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const;

	Real
	bbind_nrchi_score(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		Real & dnrchi_score_dnrchi
	) const;

	Real
	bbdep_nrchi_score(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		Real & dnrchi_score_dnrchi,
        utility::fixedsizearray1< Real, N > & dnrchi_score_dbb
	) const;


	void
	fill_rotamer_vector_bbind(
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

	void
	fill_rotamer_vector_bbdep(
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

	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples_bbind(
		utility::fixedsizearray1<Real, N> bbs
	) const;

	utility::vector1< DunbrackRotamerSampleData >
	get_all_rotamer_samples_bbdep(
		utility::fixedsizearray1<Real, N> bbs
	) const;

	Real
	get_probability_for_rotamer_bbind(
		utility::fixedsizearray1<Real, N> bbs,
		Size rot_ind
	) const;

	Real
	get_probability_for_rotamer_bbdep(
		utility::fixedsizearray1<Real, N> bbs,
		Size rot_ind
	) const;

	DunbrackRotamerSampleData
	get_rotamer_bbind(
		utility::fixedsizearray1<Real, N> bbs,
		Size rot_ind
	) const;


	DunbrackRotamerSampleData
	get_rotamer_bbdep(
		utility::fixedsizearray1<Real, N> bbs,
		Size rot_ind
	) const;

	/// @brief overrides parent class for the non-rotameric chi, but falls back
	/// on the parent class functionality for all other chi.
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


private:
	
	void
	read_rotamer_definitions( utility::io::izstream & in_rotdef );

	void
	read_bbdep_continuous_minimization_data( utility::io::izstream & in_contmin );

	void
	read_bbind_continuous_minimization_data( utility::io::izstream & in_contmin );

	void
	read_rotameric_data( utility::io::izstream & in_rotameric );

	void
	build_bbdep_rotamers(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer,
		BBDepNRChiSample< Real > const interpolated_sample,
		RotamerVector & rotamers
	) const;

	void
	build_bbind_rotamers(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		PackedDunbrackRotamer< T, N, Real > const & interpolated_rotamer,
		Size const nrchi_rotno,
		BBIndNRChiSample<> const & interpolated_sample,
		RotamerVector & rotamers
	) const;

	void
	bbind_chisamples_for_rotamer_chi(
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

	void
	bbdep_chisamples_for_rotamer_chi(
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

	Real clip_to_nrchi_range( Real chi ) const;

	void get_bbdep_nrchi_bin(
		Real nrchi,
		Size & bin_lower,
		Size & bin_upper,
		Real & nrchi_alpha
	) const;

	void get_bbind_nrchi_bin(
		Real nrchi,
		Size & bin_lower,
		Size & bin_upper,
		Real & nrchi_alpha
	) const;
    
    BBDepNRChiSample< Real >
	interpolate_bbdep_nrchi_sample(
                                   Size const packed_rotno,
                                   Size const nrchi_bin,
                                   utility::fixedsizearray1< Size, N > const bb_bin,
                                   utility::fixedsizearray1< Size, N > const bb_bin_next,
                                   utility::fixedsizearray1< Real, N > const bb_alpha
                                   ) const;
    
	BBDepNRChiSample< Real >
	interpolate_bbdep_nrchi_sample(
                                   utility::vector1< BBDepNRChiSample<> > const & nrchi_sample,
                                   utility::fixedsizearray1< Real, N > const bb_alpha
                                   ) const;


private:
	/// @brief no default ctor
	///SemiRotamericSingleResidueDunbrackLibrary();

private:

	bool const bbind_nrchi_scoring_;  // Score non-rotameric chi independently of the backbone?
	bool const bbind_nrchi_sampling_; // Sample non-rotameric chi independently of the backbone?

	Real nrchi_periodicity_; // 360 w/o symmetry; 180 w/.
	Real nrchi_lower_angle_; // Starting angle for both bbdep and bbind nrchi data

	/// The non rotameric chi is n_rotameric_chi + 1;
	/// The FArray3D is indexed as (phi, psi, nrchi).  Probabilities are input; the
	/// parameters for tricubic interpolation are stored.  The supporting information is used to snap chi values stored in a
	/// residue to a meaningful periodic range and to make tricubic interpolation possible.
	/// This data is used only if bbind_nrchi_scoring is false.
    // amw now indexed as index, nrchi
	utility::vector1< ObjexxFCL::FArray2D< BBDepScoreInterpData<N> > > bbdep_nrc_interpdata_;
	Size bbdep_nrchi_nbins_;
	Real bbdep_nrchi_binsize_;

	/// If bbind_nrchi_scoring_ is true, then this array holds the continuous minimization data for
	/// each of the packed rotamers.
	ObjexxFCL::FArray2D< Real > bbind_non_rotameric_chi_scores_;
	Size bbind_nrchi_nbins_;
	Real bbind_nrchi_binsize_;


	/// For rotamer sampling, "pseudo rotamers" are defined for the non rotameric chi.  These
	/// rotamers are arbitrary bins in chi space, and each bin gives a left and right edge
	/// as well as a median (mode?) angle which is used for rotamer sampling iff bbind_nrchi_sampling_
	/// is true.  The bins do not have a uniform width.

	// Used in both bbind and bbdep rotamer building.  Read from rotamer definition file.
	Size n_nrchi_sample_bins_;


	/// This variable is used iff bbind_nrchi_sampling_ is true;
	/// Space is not allocated if it is false.
	/// The algorithm for creating the top 95% most probable rotamers relies on a
	/// sorted order for these rotamers, which is why we need bbind_rotamers_sorted_by_probability_.
	ObjexxFCL::FArray2D< BBIndNRChiSample<> > bbind_rotamers_to_sample_;
	ObjexxFCL::FArray2D< Size > bbind_rotamers_sorted_by_probability_;

	/// This variable is used iff bbind_nrchi_sampling_ is false;
	/// Space is not allocated if it is true.
    // amw reduced
	ObjexxFCL::FArray2D< BBDepNRChiSample<> > bbdep_rotamers_to_sample_;
	ObjexxFCL::FArray3D< Size > bbdep_rotsample_sorted_order_;

};


} // namespace dunbrack
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_dunbrack_SemiRotamericSingleResidueDunbrackLibrary_HH
