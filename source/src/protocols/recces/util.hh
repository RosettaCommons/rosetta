// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/recces/util.hh
/// @brief util functions for RECCES and thermal_sampler
/// @details
///
/// @author Andy Watkins
/// @author Fang-Chieh Chou
/// @author Rhiju Das


#ifndef INCLUDED_protocols_recces_util_HH
#define INCLUDED_protocols_recces_util_HH

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <utility/io/ozstream.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>
#include <map>
#include <utility/vector1.hh>


namespace protocols {
namespace recces {

//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::scoring::ScoreType> const & get_scoretypes();

core::Size data_dim( utility::vector1< core::scoring::ScoreType > const & score_types );

core::Size data_dim();

//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn
);

//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	utility::vector1< core::scoring::ScoreType > const & scoretypes
);

//////////////////////////////////////////////////////////////////////////////
void fill_data(
	utility::vector1<float> & data,
	core::Size const count,
	utility::vector1<float> const & scores
);

//////////////////////////////////////////////////////////////////////////////
template<typename T>
void vector2disk_in1d(
	std::string const & out_filename,
	utility::vector1<T> const & out_vector
) {
	utility::io::ozstream out( out_filename.c_str(), std::ios::out | std::ios::binary );
	if ( out_vector.size() != 0 )	out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
	out.close();
}
//////////////////////////////////////////////////////////////////////////////
template<typename T>
void vector2disk_in2d(
	std::string const & out_filename,
	core::Size const dim1,
	core::Size const dim2,
	utility::vector1<T> const & out_vector
) {
	//std::cout << "dim1 " << dim1 << " dim2 " << dim2 << " out " << out_vector.size() << std::endl;
	utility::io::ozstream out( out_filename.c_str(), std::ios::out | std::ios::binary );
	runtime_assert( dim1 * dim2 == out_vector.size() );
	out.write( (const char*) &dim1, sizeof(core::Size) );
	out.write( (const char*) &dim2, sizeof(core::Size) );
	if ( out_vector.size() == 0 ) {
		std::cout << "Warning: no data available for the requested condition. Output file may be malformed.\n";
		// noop
	} else {
		out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
	}
	out.close();
}

// @brief used to compute moments of inertia, phase space volume
// exact copy of function print_base_centroid_Atoms in rb_entropy -- delete one of these!!
void
print_base_centroid_atoms_for_rb_entropy( core::conformation::Residue const & rsd, std::string filename_xyz  );

/// @brief used to output torsions from pose -- useful for clustering states, etc.
utility::vector1<core::Real>
get_torsions(
		utility::vector1<core::id::TorsionID> const & torsion_ids,
		core::pose::Pose const & pose	);

} //recces
} //protocols

#endif
