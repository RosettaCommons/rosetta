// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Align_RotamerEvaluator
/// @author James Thompson

#include <protocols/comparative_modeling/Align_RotamerEvaluator.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>

#include <ObjexxFCL/string.functions.hh>
#include <numeric/angle.functions.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.comparative_modeling.Align_RotamerEvaluator" );

namespace protocols {
namespace comparative_modeling {

Align_RotamerEvaluator::~Align_RotamerEvaluator() {}

Align_RotamerEvaluator::Align_RotamerEvaluator(
	core::pose::PoseCOP native_pose,
	std::string tag,
	core::Real const chi_dev,
	core::sequence::SequenceAlignmentOP aln
) :
	AlignEvaluator( native_pose, tag, true, aln ),
	chi_dev_(chi_dev)
{}

utility::vector1< core::pack::dunbrack::RotVector >
rots_from_pose(
	core::pose::Pose const & pose
) {
	using namespace core::pack::dunbrack;
	using utility::vector1;

	vector1< RotVector > rots;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		RotVector rot;
		rotamer_from_chi( pose.residue(ii), rot );
		rots.push_back( rot );
	}
	return rots;
}

utility::vector1< utility::vector1< core::Real > >
chis_from_pose(
	core::pose::Pose const & pose
) {
	using core::Size;
	using core::Real;
	using utility::vector1;

	vector1< vector1< Real > > chis;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		vector1< Real > chi_vec;
		core::Size const n_chi( pose.residue_type(ii).nchi() );
		for ( Size jj = 1; jj <= n_chi; ++jj ) {
			chi_vec.push_back( pose.chi(jj,ii) );
		}
		chis.push_back( chi_vec );
	}
	return chis;
}

void
Align_RotamerEvaluator::apply(
	core::pose::Pose & pose,
	std::string /*tag*/,
	core::io::silent::SilentStruct & ss
) const {
	using core::Real;
	using core::Size;
	using utility::vector1;
	using core::id::SequenceMapping;
	using namespace core::pack::dunbrack;

	SequenceMapping mapping( get_alignment(pose)->sequence_mapping(1,2) );

	tr.flush_all_channels();
	static Size const MAX_N_CHI (4); // hacky
	vector1< Size > n_rots_equiv( MAX_N_CHI, 0 );
	vector1< Size > n_chis_equiv( MAX_N_CHI, 0 );

	Size n_ali(0);
	Real const max_dev( chi_dev() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const native_ii( mapping[ii] );
		bool skip(
			native_ii == 0 ||
			native_ii > native_pose()->total_residue()
		);
		if ( !skip ) {
			++n_ali;
			RotVector model_rots, native_rots;
			rotamer_from_chi( pose.residue(ii), model_rots );
			rotamer_from_chi( native_pose()->residue(native_ii), native_rots );
			if ( model_rots.size() == native_rots.size() ) {
				for ( Size chi_idx = 1; chi_idx <= model_rots.size(); ++chi_idx ) {
					// rotamer equivalence
					if ( model_rots[chi_idx] == native_rots[chi_idx] ) {
						n_rots_equiv[chi_idx]++;
					}

					// chi deviation
					//Real const model_chi ( pose.chi(ii,chi_idx) );
					//Real const native_chi( native_pose()->chi(native_ii,chi_idx) );
					Real const model_chi ( pose.chi(chi_idx,ii) );
					Real const native_chi( native_pose()->chi(chi_idx,native_ii) );
					using numeric::nearest_angle_degrees;
					if ( nearest_angle_degrees(model_chi,native_chi) <= max_dev ) {
						n_chis_equiv[chi_idx]++;
					}
				}
			}
		} // residue is aligned
	} // for residues

	Real const coverage(
		(Real) (n_ali) / (Real) (native_pose()->total_residue())
	);
	Size const native_nres( native_pose()->total_residue() );
	for ( Size ii = 1; ii <= n_chis_equiv.size(); ++ii ) {
		using std::string;
		using ObjexxFCL::string_of;

		string const chi_tag( "chi" + string_of(ii) + "_equiv_" + tag() );
		Real const chi_equiv( coverage * (Real) n_chis_equiv[ii] / (Real) native_nres );
		ss.add_energy( chi_tag, chi_equiv );

		string const rot_tag( "rot" + string_of(ii) + "_equiv_" + tag() );
		Real const rot_equiv( coverage * (Real) n_rots_equiv[ii] / (Real) native_nres );
		ss.add_energy( rot_tag, rot_equiv );
	}
} // apply

core::Real Align_RotamerEvaluator::chi_dev() const {
	return chi_dev_;
}

} // comparative_modeling
} // protocols
