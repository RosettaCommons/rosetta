// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Align_RmsdEvaluator
/// @author James Thompson

// Unit Headers
#include <protocols/comparative_modeling/Align_RmsdEvaluator.hh>

// Package Headers
#include <protocols/comparative_modeling/coord_util.hh>

// Project Headers
#include <core/id/SequenceMapping.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <string>

static thread_local basic::Tracer tr( "protocols.comparative_modeling.Align_RmsdEvaluator" );

namespace protocols {
namespace comparative_modeling {

Align_RmsdEvaluator::~Align_RmsdEvaluator() {}

Align_RmsdEvaluator::Align_RmsdEvaluator(
	core::pose::PoseCOP native_pose,
	std::string tag,
	bool calc_gdt,
	core::sequence::SequenceAlignmentOP aln,
	bool gdt_by_TM
) :
	AlignEvaluator( native_pose, tag, true, aln ),
	calc_gdt_ (calc_gdt),
	gdt_by_TM_( gdt_by_TM ),
	report_gdt_components_(false)
{}

void
Align_RmsdEvaluator::apply(
	core::pose::Pose & pose,
	std::string /*tag*/,
	core::io::silent::SilentStruct & ss
) const {
	using core::Real;
	using ObjexxFCL::FArray2D;
	using core::scoring::xyz_gdtmm;
	using core::scoring::xyz_gdttm;

	int n_atoms;
	FArray2D< Real > p1a, p2a;
	tr.Debug << "gathering xyz coordinates ... ";
	protocols::comparative_modeling::gather_coords(
		pose, *native_pose(),
		*get_alignment(pose),
		n_atoms, p1a, p2a
	);
	tr.Debug << "gathered " << n_atoms << "coordinates" << std::endl;
	tr.flush_all_channels();

	Real rmsd( numeric::model_quality::rms_wrapper( n_atoms, p1a, p2a ) );
	ss.add_energy( "rms_" + tag(), rmsd );

  Real const coverage(
		(Real) (n_atoms) / (Real) (native_pose()->total_residue())
	);

	if ( report_aln_components() ) {
		ss.add_energy( "nres_ali_" + tag(), n_atoms );
		ss.add_energy( "coverage_" + tag(), coverage );
	}

	if ( calc_gdt() ) {

		// required for high-accuracy statistics
		core::id::SequenceMapping mapping = get_alignment(pose)->sequence_mapping(1, 2);
		std::map<Size, Size> residues;
		for (Size idx_mod = 1; idx_mod <= mapping.size1(); ++idx_mod) {
			Size idx_ref = mapping[idx_mod];
			if (idx_ref > 0) {
				residues[idx_ref] = idx_mod;
			}
		}

		if( gdt_by_TM() ){
			tr.Debug << "computing gdttm for " << tag() << std::endl;
			Real gdttm, gdtha;
			xyz_gdttm( p1a, p2a, gdttm, gdtha );

			ss.add_energy( "gdttm_" + tag(), coverage * gdttm );
			ss.add_energy( "gdtha_" + tag(), coverage * gdtha );

		} else {
			tr.Debug << "computing gdtmm for " << tag() << std::endl;
			Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
			Real gdtmm = xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
			ss.add_energy ( "gdtmm_" + tag(), coverage * gdtmm );
			if ( report_gdt_components() ) {
				ss.add_energy( "m11", m_1_1 );
				ss.add_energy( "m22", m_2_2 );
				ss.add_energy( "m33", m_3_3 );
				ss.add_energy( "m43", m_4_3 );
				ss.add_energy( "m74", m_7_4 );
			}

			tr.Debug << "computing high-accuracy statistics for " << tag() << std::endl;
			ss.add_energy("gdtha_" + tag(), core::scoring::gdtha(*native_pose(), pose, residues));
			ss.add_energy("gdtsc_" + tag(), core::scoring::gdtsc(*native_pose(), pose, residues));

			tr.Debug << "computing maxsub for " << tag() << std::endl;
			ss.add_energy( "maxsub_" + tag(), core::scoring::xyz_maxsub( p1a, p2a, n_atoms ) );
		}
	}
}

} // comparative_modeling
} // protocols
