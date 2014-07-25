// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/membrane/MembraneData.hh
///
///	@brief		Membrane Scorefunction Statistics Class
///	@details	Access to centroid rotamer pair potential and membrane
///				environemnt potential statistics. Loads data in database tables
///				in construction. All immutable data - no setters please and const!.
///				Last Modified: 5/13/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_membrane_MembraneData_hh
#define INCLUDED_core_scoring_membrane_MembraneData_hh

// Unit Headers
#include <core/scoring/membrane/MembraneData.fwd.hh>
#include <core/scoring/EnvPairPotential.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace scoring {
namespace membrane {

/// @brief Mmebrane Data Class
/// @details Stores membrane potential statistics and cenlist access
class MembraneData : public core::scoring::EnvPairPotential {
	
public:
	
	/// @brief Default Constructor
	MembraneData();
	
	/// @brief Destructor
	~MembraneData();
	
	/// @brief Finalize Setup of MP Base Potential Class
	virtual void finalize( pose::Pose & pose ) const;
	
	/// @brief Access Cenlist from Pose
	/// @details Pose must already contain a cenlist object or this method will fail.
	CenListInfo const &
	get_cenlist_from_pose( pose::Pose const & pose ) const;
	
private: // database io
	
	/// @brief Database IO Helper Methods for Membrane
	void load_menv_info();
	
public: // accessors
	
	/// @brief Membrane Base Potential Statistics
	ObjexxFCL::FArray3D< Real > mem_env_log6() const;
	ObjexxFCL::FArray3D< Real > mem_env_log10() const;
	
	ObjexxFCL::FArray1D< Real > mem_cbeta_den6() const;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den12() const;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den6() const;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den12() const;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den6() const;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den12() const;
	
	ObjexxFCL::FArray4D< Real > mem_pair_log() const;
	
	/// @brief Env Pair Potential Statistics
	Real cen_dist5_pad() const;
	Real cen_dist6_pad() const;
	Real cen_dist7_pad() const;
	Real cen_dist10_pad() const;
	Real cen_dist12_pad() const;
	
	Real cen_dist5_pad_plus() const;
	Real cen_dist6_pad_plus() const;
	Real cen_dist7_pad_plus() const;
	Real cen_dist10_pad_plus() const;
	Real cen_dist12_pad_plus() const;
	
	Real cen_dist5_pad_minus() const;
	Real cen_dist7_pad_minus() const;
	Real cen_dist10_pad_minus() const;
	Real cen_dist12_pad_minus() const;
	
	Real cen_dist5_pad_hinv() const;
	Real cen_dist6_pad_hinv() const;
	Real cen_dist7_pad_hinv() const;
	Real cen_dist10_pad_hinv() const;
	Real cen_dist12_pad_hinv() const;
	

private: // data
	
	/// Membrane Environment Pair Potential Statistics
	ObjexxFCL::FArray3D< Real > mem_env_log6_;
	ObjexxFCL::FArray3D< Real > mem_env_log10_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_den12_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den12_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den6_;
	ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den12_;
	ObjexxFCL::FArray4D< Real > mem_pair_log_;
	
	/// Centroid Rotamer Pair Potential Statistics
	Real const cen_dist5_pad_;
	Real const cen_dist6_pad_;
	Real const cen_dist7_pad_;
	Real const cen_dist10_pad_;
	Real const cen_dist12_pad_;
	
	Real const cen_dist5_pad_plus_;
	Real const cen_dist6_pad_plus_;
	Real const cen_dist7_pad_plus_;
	Real const cen_dist10_pad_plus_;
	Real const cen_dist12_pad_plus_;
	
	Real const cen_dist5_pad_minus_;
	Real const cen_dist7_pad_minus_;
	Real const cen_dist10_pad_minus_;
	Real const cen_dist12_pad_minus_;
	
	Real const cen_dist5_pad_hinv_;
	Real const cen_dist6_pad_hinv_;
	Real const cen_dist7_pad_hinv_;
	Real const cen_dist10_pad_hinv_;
	Real const cen_dist12_pad_hinv_;
	
}; // MembraneData
	
} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_membrane_MembraneData_hh

