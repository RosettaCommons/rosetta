// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/DipolarCoupling.hh
/// @brief  Uses NMR DC for scoring
/// @author Lei Shi

#ifndef INCLUDED_core_scoring_DipolarCoupling_hh
#define INCLUDED_core_scoring_DipolarCoupling_hh

#include <core/scoring/DipolarCoupling.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/datacache/CacheableData.hh>
#include <numeric/numeric.functions.hh>
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <utility/vector1_bool.hh>

namespace core {
namespace scoring {

void store_DC_in_pose(DipolarCouplingOP, core::pose::Pose&);
DipolarCouplingOP retrieve_DC_from_pose(core::pose::Pose&);
DipolarCouplingCOP retrieve_DC_from_pose(core::pose::Pose const&);

/// @brief DipolarCouplings are mainly handled by this class
/// @detail related classed: DC --- a single line in an DC file - representing a single dc coupling
///                         DipolarCouplingEnergy -- an energy method which triggers computations handled by this class.
///
class DipolarCoupling: public basic::datacache::CacheableData {
	//	friend class DipolarCouplingEnergy;
public:
	// typedefs
	typedef core::Size Size;
	typedef utility::vector1<core::scoring::DC> DC_lines;

public:
//need to make sure it reads as the flags
	/// @brief standard c'stor -- will access option -in:file:dc to read DC data
	DipolarCoupling( std::string const& filename ="" ) {
		if (filename != "" ) {
			read_DC_file( filename );
			}else {
			read_DC_file( );
		}
	}

	//explicit copy c'stor to initialize buffers
	DipolarCoupling(DipolarCoupling const& other);

	//explicit assignment operator to initialize buffers
	DipolarCoupling& operator=(DipolarCoupling const & other);

	//explicit destructor because we use raw pointers for buffers
	~DipolarCoupling() {}

	//this class lives in the PoseCache.... need to provide clone()
	basic::datacache::CacheableDataOP clone() const {
		return basic::datacache::CacheableDataOP( new DipolarCoupling(*this) );
	}

	/// @brief compute dc score for given pose (non-constant due to membrane)
	core::Real compute_dcscore(core::pose::Pose & pose);

	void show(std::ostream&) const;

  /// @brief get the raw DC data
  inline DC_lines const& get_DC_data() const {
    return All_DC_lines_;
  }

private:
	/// @brief read DC data from file
	void read_DC_file( std::string const& filename );
	void read_DC_file( );

private:
	/// some internal buffers in
	DC_lines All_DC_lines_;
};


/////////////////////////////////////////////////
//@brief short class that stores the DC data lines
/////////////////////////////////////////////////
class DC {

public:
	enum DC_TYPE {
		DC_TYPE_NH = 1, DC_TYPE_NC, DC_TYPE_CH, DC_TYPE_CC
	};

	DC() :
		DCval_computed_(),
		type_(),
		res1_(),
		res2_(),
		DCval_(),
		DCerr_(),
		weight_()
	{}

	DC(Size res1, std::string const& atom1, Size res2, std::string const& atom2, Real DCval, Real DCerr, Real weight) :
		DCval_computed_(-999), f1ij_(0.0), f2ij_(0.0),
		type_(get_DC_data_type(atom1, atom2)),
		res1_(res1), res2_(res2),
		atom1_(atom1), atom2_(atom2),
		DCval_(DCval), DCerr_(DCerr),
		weight_(weight)
	{}

  DC_TYPE type() const {
                return type_;
  }

  DC_TYPE get_DC_data_type(std::string const & atom1, std::string const & atom2);

  inline Size res1() const {
    return res1_;
  }

  inline Size res2() const {
    return res2_;
  }

  std::string const& atom1() const {
    return atom1_;
  }

  std::string const& atom2() const {
    return atom2_;
  }

  inline Real DCval() const {
    return DCval_;
  }

  inline Real DCerr() const {
    return DCerr_;
  }

  Real const& DCcomputed() const {
    return DCval_computed_;
  }

  Real& DCcomputed() {
    return DCval_computed_;
  }

  Vector f1ij() const {
    return f1ij_;
  }

  Vector f2ij() const {
    return f2ij_;
  }


  Real weight() const {
    return weight_;
  }

  void set_weight(Real w_in) {
    weight_ = w_in;
  }


//Dcconst()/3/r^3/2
  inline Real Dconst() const {
          core::Real Dcnst(0.0);
          if (type() == DC_TYPE_NH)
                  Dcnst = 10.735/2;
                  //Dcnst = 36.5089/3/1.04/1.04/1.04/2;
          if (type() == DC_TYPE_NC)
                  Dcnst = 9.179532/3/1.341/1.341/1.341/2;
          if (type() == DC_TYPE_CH)
                  Dcnst = 90.55324/3/1.08/1.08/1.08/2;
          if (type() == DC_TYPE_CC)
                  Dcnst = 22.76805/3/1.525/1.525/1.525/2;
          runtime_assert( Dcnst != 0 );
          return Dcnst;
  }

  friend class DipolarCoupling;

  void show(std::ostream&) const;

public :
  Real DCval_computed_;
  core::Vector f1ij_;
  core::Vector f2ij_;

private:
  DC_TYPE type_;
  Size res1_, res2_;
  std::string atom1_,atom2_;
  Real DCval_, DCerr_;
  Real weight_;

};

extern std::ostream& operator<<(std::ostream&, DipolarCoupling const&);
extern std::ostream& operator<<(std::ostream&, DC const&);

} //scoring
} //core
#endif
