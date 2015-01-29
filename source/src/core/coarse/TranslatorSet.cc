// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ChemicalManager.cc
/// @brief  Chemical manager class
/// @author Oliver Lange (olange@u.washington.edu)

#include <core/coarse/TranslatorSet.hh>

// Package headers
#include <core/coarse/Translator.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace core {
namespace coarse {

static thread_local basic::Tracer TR( "core.coarse.TranslatorSet" );

using namespace chemical;

TranslatorSet::TranslatorSet(
	const RuleSet &rules,
	ResidueTypeSetCAP residue_set,
	ResidueTypeSetAP coarse_set)
: residue_set_(residue_set),
	coarse_residue_set_(coarse_set)
{
	using namespace std;
	TR.Error << "constructing TranslatorSet..." << "\n";
	for ( ResidueTypeSet::const_residue_iterator it=residue_set_->all_residues_begin(),
				eit=residue_set_->all_residues_end();
				it!=eit;
				 ++it )
		{
			const ResidueType& res = *(it->second);
			T("coarse.setup") << "in loop " << res.name() << " " << res.name3() << "  aa: " << res.aa() <<  "\n";
			if ( coarse_residue_set_->has_name( res.name() ) ) {
	TranslatorOP cmap_ptr =
		new Translator( rules,ResidueTypeCOP( res ),
			ResidueTypeAP( const_cast<ResidueType&> (coarse_residue_set_->name_map( res.name() )) )
		); //the non-const object is only needed in the construction of Translator (fix_geometry)
	//, a CAP is stored for later...
	cmap_ptr->pretty_print(cout);
	coarse_maps_.insert(TranslatorMap::value_type(it->first,cmap_ptr));
			} else {
	Warning() << res.name() << " not found in coarse_residue set... ignored..." << "\n";
			}
		}
}

void
TranslatorSet::pretty_print( std::ostream &os ) const
{
	for ( const_iterator it=begin(), eit=end(); it!=eit; ++it ) {
		it->second->pretty_print(os);
		os << std::endl;
	}
}

TranslatorCOP const&
TranslatorSet::default_for_aa( chemical::AA aa ) const
{
	chemical::ResidueTypeCOPs types=residue_set_->aa_map(aa);
	if (types.size() > 1) {
		Warning() << "WARNING: coarse/Translator.cc: using the first residue in the set with AA "
				<< aa << "\n" << "do not know if that is the best choice. Might be good to add " <<
			"knowledge about the 'default' residue_type to the residueset " << "\n";
	};
	return (coarse_maps_[types.front()->name()]);
}

void
 TranslatorSet::coarsify( pose::Pose &pose_out, pose::Pose const &pose_in ) const
{
	pose_out.clear();
	for (uint seqpos = 1;seqpos<=pose_in.total_residue();seqpos++) {
		pose_out.append_residue_by_bond( *coarsify(pose_in.residue(seqpos)) );
	};
	pose_out.fold_tree( kinematics::FoldTree( pose_out.total_residue() ) );
}


conformation::ResidueOP TranslatorSet::coarsify(const conformation::Residue& fine_rsd) const {
debug_assert( has(fine_rsd.name()) ); // silly check
	return coarse_maps_[fine_rsd.name()]->coarsify(fine_rsd);
}

bool TranslatorSet::has(TranslatorSet::ResName name) const {
	return ( coarse_maps_.find(name) != coarse_maps_.end() );
}



}
}
