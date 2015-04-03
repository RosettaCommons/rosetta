// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for generating statistics from orbitals
///
/// @details
/// This is an attempt to be transparent in how the orbital-hydrogen interactions were converted into
/// a KBP. In general terms, a set of proteins from PDB were taken and the shortest distance from an
/// orbital to a polar or aromatic hydrogen was recorded and binned. In addition to the shortest distance,
/// an angle was taken. The angle considered is the angle between the base atom that contains the orbital, the
/// orbital, and the hydrogen. For polar hydrogens, only sidechain polar hydrogens were considered.
///
/// For protein interactions, there are 7 classes (orbital types) of orbitals that statistics are generated.
/// These 7 types are mapped using a map to enum data structure. For each sidechain interaction, only the shortest
/// distance between any given orbital type to a hydrogen is calculated. That means, that for each sidechain interaction
/// only 1 distance, 1 angle, and 1 class is recorded.
///
/// Bin sizes were calculated by .1A for distance and cos of 1 for angles. See below:
///
///           angle
///       -1 -.9 -.8 -.7..........
/// d .1 | 0   0   0   0
/// i .2 | 500 0   0   0
/// s .3 | 25  0   0   0
/// t .4 |  0  30  5   0
///
/// This is not the original code that was used to generate the statistics. The original code was much more
/// convoluted than this because I had no idea how to program. I wrote this piece for clarity. I have tested
/// it and it produces the same results.
///
/// @author
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_core_scoring_orbitals_OrbitalsStatistics_hh
#define INCLUDED_core_scoring_orbitals_OrbitalsStatistics_hh

#include <core/pose/Pose.fwd.hh>
#include <numeric/histograms/TwoDHistogram.hh>
#include <utility/vector1.hh>
#include <map>
#include <core/types.hh>
#include <fstream>

//Auto Headers
namespace core{
namespace scoring{
namespace orbitals{

class OrbitalsStatistics {
public:

	enum orbital_type_name{
		C_pi_sp2=1,
		N_pi_sp2,
		N_p_sp2,
		O_pi_sp2,
		O_p_sp2,
		O_p_sp3,
		S_p_sp3
	};


	OrbitalsStatistics();

	//void orbital_orbital(core::pose::Pose & pose);

	void sc_H_orbital(core::pose::Pose & pose);
	/// Undefined, commenting out to fix PyRosetta build  void bb_stats(core::pose::Pose & pose);

	void increment_histogram_bin(
			core::Real & distance,
			core::Real & angle,
			numeric::histograms::TwoDHistogram<core::Size, core::SSize> & histogram
	);

	void bb_stats(
			core::pose::Pose & pdb
	);

	numeric::histograms::TwoDHistogram<core::Size, core::SSize> get_2D_histogram();
	utility::vector1< numeric::histograms::TwoDHistogram<core::Size, core::SSize> >  get_histogram_vector();
	core::Size get_number_of_histograms();


private:
	OrbitalsStatistics(OrbitalsStatistics const & );  /// Non-copyable due to std::ofstream member, we need so PyRosetta builder
	                                                  /// can figure out not to try to create copy constructor.

	utility::vector1< numeric::histograms::TwoDHistogram<core::Size, core::SSize> > histogram_vector_;
	numeric::histograms::TwoDHistogram<core::Size, core::SSize> twoD_histogram_;
	core::Size number_of_histograms_;
	std::map<std::string, orbital_type_name> orbital_type_2_enum_;
	std::ofstream statistics_output_;


};


}
}
}

#endif /* INCLUDED_core_scoring_orbitals_generate_statistics_OrbitalStatistics_hh */
