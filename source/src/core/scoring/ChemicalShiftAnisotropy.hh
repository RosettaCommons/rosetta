// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ChemicalShiftAnisotropy.hh
/// @brief  Uses NMR CSA for scoring
/// @author Lei Shi

#ifndef INCLUDED_core_scoring_ChemicalShiftAnisotropy_hh
#define INCLUDED_core_scoring_ChemicalShiftAnisotropy_hh

#include <core/scoring/ChemicalShiftAnisotropy.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/datacache/CacheableData.hh>
#include <numeric/numeric.functions.hh>
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <utility>
#include <utility/vector1_bool.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

void store_CSA_in_pose(ChemicalShiftAnisotropyOP, core::pose::Pose&);
ChemicalShiftAnisotropyOP retrieve_CSA_from_pose(core::pose::Pose&);
ChemicalShiftAnisotropyCOP retrieve_CSA_from_pose(core::pose::Pose const&);

/// @brief ChemicalShiftAnisotropys are mainly handled by this class
/// @detail related classed: CSA --- a single line in an CSA file - representing a single csa coupling
///                         ChemicalShiftAnisotropyEnergy -- an energy method which triggers computations handled by this class.
///
class ChemicalShiftAnisotropy: public basic::datacache::CacheableData {
	// friend class ChemicalShiftAnisotropyEnergy;
public:
	// typedefs
	typedef core::Size Size;
	typedef utility::vector1<core::scoring::CSA> CSA_lines;

public:
	//need to make sure it reads as the flags
	/// @brief standard c'stor -- will access option -in:file:csa to read CSA data
	ChemicalShiftAnisotropy( std::string const& filename ="" ) {
		if ( filename != "" ) {
			read_CSA_file( filename );
		} else {
			read_CSA_file( );
		}
	}

	//explicit copy c'stor to initialize buffers
	ChemicalShiftAnisotropy(ChemicalShiftAnisotropy const& other);

	//explicit assignment operator to initialize buffers
	ChemicalShiftAnisotropy& operator=(ChemicalShiftAnisotropy const & other);

	//explicit destructor because we use raw pointers for buffers
	~ChemicalShiftAnisotropy() override = default;

	//this class lives in the PoseCache.... need to provide clone()
	basic::datacache::CacheableDataOP clone() const override {
		return basic::datacache::CacheableDataOP( new ChemicalShiftAnisotropy(*this) );
	}

	/// @brief compute csa score for given pose (non-constant due to membrane)
	core::Real compute_csascore(core::pose::Pose & pose);

	void show(std::ostream&) const;

	/// @brief get the raw CSA data
	inline CSA_lines const& get_CSA_data() const {
		return All_CSA_lines_;
	}

private:
	/// @brief read CSA data from file
	void read_CSA_file( std::string const& filename );
	void read_CSA_file( );

private:
	/// some internal buffers in
	CSA_lines All_CSA_lines_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/////////////////////////////////////////////////
//@brief short class that stores the CSA data lines
/////////////////////////////////////////////////
class CSA {

public:
	CSA() :
		CSAval_computed_(),
		res1_(),
		res2_(),
		res3_(),
		sigma1_(),
		sigma2_(),
		sigma3_(),
		alpha_(),
		beta_(),
		gamma_(),
		CSAval_(),
		CSAerr_(),
		weight_()
	{
	}

	CSA(Size res1, std::string  atom1, Real sigma1, Real sigma2, Real sigma3, Real CSAval, Real CSAerr, Real weight) :
		CSAval_computed_(-999),f1ij_(0.0),f2ij_(0.0),f3ij_(0.0),
		res1_(res1), res2_(res1-1), res3_(res1),
		atom1_(std::move(atom1)), atom2_("C"), atom3_("CA"),
		sigma1_(sigma1), sigma2_(sigma2), sigma3_(sigma3),
		alpha_(0), beta_(105), gamma_(0),
		CSAval_(CSAval), CSAerr_(CSAerr),
		weight_(weight)
	{}

	inline Size res1() const {
		return res1_;
	}

	inline Size res2() const {
		return res2_;
	}

	inline Size res3() const {
		return res3_;
	}

	std::string const& atom1() const {
		return atom1_;
	}

	std::string const& atom2() const {
		return atom2_;
	}

	std::string const& atom3() const {
		return atom3_;
	}

	inline Real sigma1() const {
		return sigma1_;
	}

	inline Real sigma2() const {
		return sigma2_;
	}

	inline Real sigma3() const {
		return sigma3_;
	}

	inline Real alpha() const {
		return alpha_;
	}

	inline Real beta() const {
		return beta_;
	}

	inline Real gamma() const {
		return gamma_;
	}

	inline Real CSAval() const {
		return CSAval_;
	}

	inline Real CSAerr() const {
		return CSAerr_;
	}

	Real const& CSAcomputed() const {
		return CSAval_computed_;
	}

	Real& CSAcomputed() {
		return CSAval_computed_;
	}

	Vector f1ij() const {
		return f1ij_;
	}

	Vector f2ij() const {
		return f2ij_;
	}

	Vector f3ij() const {
		return f3ij_;
	}


	Real weight() const {
		return weight_;
	}

	void set_weight(Real w_in) {
		weight_ = w_in;
	}

	friend class ChemicalShiftAnisotropy;

	void show(std::ostream&) const;

public :
	Real CSAval_computed_;
	core::Vector f1ij_;
	core::Vector f2ij_;
	core::Vector f3ij_;

private:
	Size res1_, res2_, res3_;
	std::string atom1_,atom2_,atom3_;
	Real sigma1_, sigma2_, sigma3_;
	Real alpha_, beta_, gamma_;
	Real CSAval_, CSAerr_;
	Real weight_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

extern std::ostream& operator<<(std::ostream&, ChemicalShiftAnisotropy const&);
extern std::ostream& operator<<(std::ostream&, CSA const&);

} //scoring
} //core
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_ChemicalShiftAnisotropy )
#endif // SERIALIZATION


#endif
