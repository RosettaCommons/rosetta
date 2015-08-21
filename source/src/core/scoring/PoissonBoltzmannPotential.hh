// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/PoissonBoltzmannPotential.hh
/// @brief  Poisson Boltzmann potential class delcaration
/// @author Yifan Song (yfsong@uw.edu)

#ifndef INCLUDED_core_scoring_PoissonBoltzmannPotential_HH
#define INCLUDED_core_scoring_PoissonBoltzmannPotential_HH

// Unit Headers
#include <core/scoring/PoissonBoltzmannPotential.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/file/PathName.hh>

#include <map>
#include <vector>

namespace core {
namespace scoring {

class PoissonBoltzmannPotential : public utility::pointer::ReferenceCount
{
public:
	static const std::string APBS_CONFIG_EXT;
	static const std::string APBS_PQR_EXT;
	static const std::string APBS_DX_EXT;
	static const std::string DEFAULT_APBS_PATH;


	PoissonBoltzmannPotential();

	virtual ~PoissonBoltzmannPotential(); // auto-removing definition from header{}

	core::Real get_potential(ObjexxFCL::FArray3D< core::Real > const & potential,
		numeric::xyzVector<core::Real> const & cartX) const;
	void
	eval_PB_energy_residue(
		core::conformation::Residue const & rsd,
		Real & PB_energy_residue,
		Real & PB_energy_backbone,
		Real & PB_energy_sidechain,
		Real const & PB_burial_weight
	) const;

	//////////////////////////////////
	//////////////////////////////////
	// functions to convert between indices and cartesian coords
	inline void cart2idx( numeric::xyzVector<core::Real> const & cartX, numeric::xyzVector<core::Real> & idxX ) const {
		idxX = c2i_*(cartX-lower_bound_) + numeric::xyzVector<core::Real> (1,1,1);
	}

	template<class Q>
	inline void idx2cart( numeric::xyzVector<Q> const & idxX , numeric::xyzVector<core::Real> &cartX ) const {
		cartX = i2c_*numeric::xyzVector<core::Real>(idxX - numeric::xyzVector<Q> (1,1,1)) + lower_bound_;
	}
	numeric::xyzVector< core::Real > lower_bound() const {
		return lower_bound_;
	}
	numeric::xyzVector< core::Real > upper_bound() const {
		return upper_bound_;
	}

	bool out_of_bounds(numeric::xyzVector< core::Real > const cartX) const {
		for ( core::Size i=0; i<3; ++i ) {
			if ( cartX[i] < lower_bound()[i] ) return true;
			if ( cartX[i] > upper_bound()[i] ) return true;
		}
		return false;
	}

	/// Execute ABPS to freshly compute the electrotatic field.
	/// @param pose  The pose
	/// @param state_tag Arbitrary string for generating APBS files.  e.g. The current energy state.
	/// @param is_residue_charged_by_name Which residues are charged?  The key is the residue name.
	void solve_pb( core::pose::Pose const & pose,
		std::string const & state_tag,
		std::map<std::string, bool> const & is_residue_charged_by_name );
private:
	numeric::xyzMatrix< core::Real > i2c_, c2i_;
	numeric::xyzVector< core::Real > lower_bound_;
	numeric::xyzVector< core::Real > upper_bound_;
	numeric::xyzVector< core::Real > grid_spacing_;
	numeric::xyzVector< core::Size > n_grid_;
	ObjexxFCL::FArray3D< core::Real > potential_;

	std::string config_filename_;
	std::string pqr_filename_;
	std::string dx_filename_;
	std::string apbs_path_;  // full path name to the APBS executable: e.g. /usr/bin/apbs.exe
	bool calcenergy_;

	/// Prepare ABPS - generate .in and .pqr
	void write_config (
		core::pose::Pose const & pose) const;

	/// Read & load the APBS results
#ifdef LINK_APBS_LIB
	void load_potential( const double grid_meta[],
											 const double pot[] );
#else
	void load_APBS_potential();
#endif

	/// Write out .pqr
	void write_pqr( core::pose::Pose const & pose,
		std::map<std::string, bool> const & is_residue_charged_by_name) const;

};

}
}

#endif
