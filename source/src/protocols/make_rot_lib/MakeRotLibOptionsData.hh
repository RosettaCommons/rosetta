// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJobInputter.hh
/// @brief  Header file for MakeRotLibOptionsData class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )


#ifndef INCLUDED_protocols_make_rot_lib_MakeRotLibOptionsData_hh
#define INCLUDED_protocols_make_rot_lib_MakeRotLibOptionsData_hh

// unit headers
#include <protocols/make_rot_lib/MakeRotLibOptionsData.fwd.hh>

// core headers
#include <core/types.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// basic headers
#include <basic/options/option.hh>

// c++ headers
#include <string>

namespace protocols {
namespace make_rot_lib {

// Describe range of torsion angle values
struct TorsionRange {
  core::Real low;
  core::Real high;
  core::Real step;
};

typedef utility::vector1< TorsionRange > TorsionRangeVec;
typedef utility::vector1< TorsionRangeVec > TorsionRangeVecVec;

// Hold a centoids starting value and rotamer bin number
struct CentroidRotNum {
  core::Real angle;
  core::Size rot_num;
};

typedef utility::vector1< CentroidRotNum > CentroidRotNumVec;
typedef utility::vector1< CentroidRotNumVec > CentroidRotNumVecVec;

enum MakeRotLibPolymerType { PEPTIDE, PEPTOID };

///@details Stores all options contained in a MakeRotLib option file
class MakeRotLibOptionsData : public utility::pointer::ReferenceCount
{
public:

	/// @brief ctor
  MakeRotLibOptionsData( std::string filename );

	/// @breif dtor
	~MakeRotLibOptionsData() {}

	/// @brief acessors
	std::string get_name() const { return name_; }
	core::Size get_n_chi() const { return n_chi_; }
	core::Size get_n_centroids() const { return n_centroids_; }
	TorsionRange get_omg_range() const { return omg_range_; }
	TorsionRange get_phi_range() const { return phi_range_; }
	TorsionRange get_psi_range() const { return psi_range_; }
	TorsionRange get_eps_range() const { return eps_range_; }
	TorsionRangeVec get_chi_data() const { return chi_ranges_; }
	CentroidRotNumVecVec get_centroid_data() const { return centroid_data_; }
	MakeRotLibPolymerType get_polymer_type() const { return polymer_type_; }

private:

  // AA name
  std::string name_;

  // number of chi torsions
  core::Size n_chi_;

  // store sets of torsion ranges for each backbone torsion
  TorsionRange omg_range_;
  TorsionRange phi_range_;
  TorsionRange psi_range_;
  TorsionRange eps_range_;

  // store sets of torsion ranges for each chi torsion
  TorsionRangeVec chi_ranges_;

  // number of centroids
  core::Size n_centroids_;

  // centroid data
  CentroidRotNumVecVec centroid_data_;

	// polymer type
	MakeRotLibPolymerType polymer_type_;

}; // MakeRotLibOptionsData

} // namespace make_rot_lib
} // namespace protocols

#endif //INCLUDED_protocols_make_rot_lib_MakeRotLibJobOptionsData_hh
