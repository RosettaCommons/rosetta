// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CompositeFragmentSelector.cc
/// @brief provides a selector that picks fragments diferent enough from the fragments that have been already picked
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/DiversifyDihedralsSelector.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <cmath>

namespace protocols {
namespace frag_picker {

using namespace core;

static basic::Tracer trDiversifyDihedralsSelector(
	"protocols.frag_picker.DiversifyDihedralsSelector");


/// @brief  Selects desired number of fragments from a given set of candidates
void DiversifyDihedralsSelector::select_fragments(
	ScoredCandidatesVector1 const& in,
	ScoredCandidatesVector1& out )
{

	if ( in.size()==0 ) return;

	core::Size len = in[1].first->get_length();

	if ( phi_.size() < len ) {
		phi_.resize(len);
		psi_.resize(len);
	}

	out.push_back( in[1] );
	for ( core::Size i=2; i<=in.size(); i++ ) {
		if ( out.size() >= frags_per_pos() ) break;
		bool is_ok = true;
		for ( core::Size j=1; j<=out.size(); j++ ) {
			core::Real rms = dihedral_rmsd(in[i].first, out[j].first);
			if ( rms<cutoff_ ) {
				is_ok = false;
				trDiversifyDihedralsSelector.Trace<<"Phi-Psi rmsd is "<<rms<<" fragment "<< *in[i].first<<" denied"<<std::endl;
				break;
			}
		}
		if ( is_ok ) {
			out.push_back( in[i] );
			trDiversifyDihedralsSelector.Trace<<"Fragment "<< *in[i].first<<" passed"<<std::endl;
		}
	}
	trDiversifyDihedralsSelector<<out.size()<<" fragments passed through DiversifyDihedralsSelector at query position "
		<< in[1].first->get_first_index_in_query()<<std::endl;
}

core::Real DiversifyDihedralsSelector::dihedral_rmsd(FragmentCandidateOP f1,FragmentCandidateOP f2) {

	core::Real rms = 0.0;
	debug_assert ( f1->get_length() == f2->get_length() );
	for ( core::Size k=1; k<=f1->get_length(); k++ ) {
		core::Real d = f1->get_residue(k)->phi() - f2->get_residue(k)->phi();
		rms += d*d;
		d = f1->get_residue(k)->psi() - f2->get_residue(k)->psi();
		rms += d*d;
		d = f1->get_residue(k)->omega() - f2->get_residue(k)->omega();
		rms += d*d;
	}

	return sqrt(rms) / (3.0 * f1->get_length());
}

} // frag_picker
} // protocols
