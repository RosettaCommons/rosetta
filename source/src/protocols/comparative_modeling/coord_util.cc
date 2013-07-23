// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/comparative_modeling/coord_util.cc
/// @author James Thompson

// AUTO-REMOVED #include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/conformation/Residue.hh>

#include <ObjexxFCL/FArray2D.hh>

#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

void gather_coords(
	core::pose::Pose const & model,
	core::pose::Pose const & native,
	core::sequence::SequenceAlignment const & aln,
	int & natoms,
	ObjexxFCL::FArray2D< core::Real > & p1a,
	ObjexxFCL::FArray2D< core::Real > & p2a,
	std::string const & atom_name = "CA"
) {
	using core::Size;
	using namespace core::id;
	using namespace core::sequence;
	SequenceMapping mapping( aln.sequence_mapping(1,2) );

	natoms = 0;
	for ( Size ii = 1; ii <= model.total_residue(); ++ii ) {
		Size const native_ii( mapping[ii] );
		bool skip(
			native_ii == 0 ||
			native_ii > native.total_residue() ||
			(!model.residue(ii).is_protein())
		);
		if ( !skip ) ++natoms;
	}
	p1a.dimension(3,natoms);
	p2a.dimension(3,natoms);

	Size n_gap(0);
	for ( Size ii = 1; ii <= model.total_residue(); ++ii ) {
		Size const native_ii(mapping[ii]);
		bool skip(
			native_ii == 0 ||
			native_ii > native.total_residue() ||
			(!model.residue(ii).is_protein())
		);
		if ( skip ) {
			n_gap++;
		} else {
			using core::Real;
			numeric::xyzVector< Real > model_xyz ( model.residue(ii).xyz(atom_name) );
			numeric::xyzVector< Real > native_xyz( native.residue(native_ii).xyz(atom_name) );
			for ( Size jj = 1; jj <= 3; ++jj ) {
				p1a(jj,ii - n_gap) = native_xyz[jj-1];
				p2a(jj,ii - n_gap) = model_xyz [jj-1];
			}
		}
	}
} // gather_coords

}
}
