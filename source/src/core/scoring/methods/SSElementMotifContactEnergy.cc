// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/motif/MotifScore.cc
/// @brief  Will's motif score to determine how well packed the protein core is
/// @author TJ Brunette
///
// Unit headers
#include <core/scoring/methods/SSElementMotifContactEnergyCreator.hh>
#include <core/scoring/methods/SSElementMotifContactEnergy.hh>


// Project Headers

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/EnergyMap.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/motif/reference_frames.hh>
#include <numeric/xyzTransform.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <string>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace methods {

using utility::vector1;
using core::Size;
using core::Real;
using std::set;
using namespace core::scoring::motif;

static THREAD_LOCAL basic::Tracer TR( "core.scoring.SSElementMotifContactEnergy" );

using namespace core::scoring::motif;

methods::EnergyMethodOP
SSElementMotifContactEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new SSElementMotifContactEnergy );
}

ScoreTypes
SSElementMotifContactEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( ss_contact_worst );
	return sts;
}

SSElementMotifContactEnergy::SSElementMotifContactEnergy():
	parent(methods::EnergyMethodCreatorOP( new SSElementMotifContactEnergyCreator ) )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	mman_ = core::scoring::motif::MotifHashManager::get_instance();
}


/// @brief get ss_elements
vector1<std::pair<Size,Size> > SSElementMotifContactEnergy::get_ss_elements(const Pose & pose) const{
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	char lastSecStruct = dssp.get_dssp_secstruct( 1 );
	Size startSS = 0;
	Size endSS = 0;
	vector1<std::pair<Size,Size> > ss_elements;
	if ( dssp.get_dssp_secstruct(1)=='H' || dssp.get_dssp_secstruct(1)=='E' ) {
		startSS=  1;
	}
	for ( core::Size ii = 2; ii <= pose.size(); ++ii ) {
		if ( (dssp.get_dssp_secstruct(ii) == 'H' || dssp.get_dssp_secstruct(ii)) && lastSecStruct == 'L' ) {
			startSS = ii;
		}
		if ( dssp.get_dssp_secstruct(ii)!=lastSecStruct && (lastSecStruct ==  'H' || lastSecStruct ==  'E') ) {
			endSS = ii-1;
			if ( endSS-startSS >= 2 ) {
				ss_elements.push_back(std::make_pair(startSS,endSS));
			}
		}
		lastSecStruct = dssp.get_dssp_secstruct(ii);
	}
	return(ss_elements);
}

Size SSElementMotifContactEnergy::which_ssElement(Size res,vector1<std::pair<Size,Size> > ssElements) const{
	for ( Size ii=1; ii<=ssElements.size(); ++ii ) {
		if ( res>= ssElements[ii].first && res<=ssElements[ii].second ) {
			return(ii);
		}
	}
	return(0);
}

Size SSElementMotifContactEnergy::get_SSelements_in_contact(Size element,vector1<std::pair<Size,Size> > ssElements, const Pose & pose) const{
	using core::kinematics::FoldTree;
	set<Size> ssElements_in_contact;
	core::scoring::dssp::Dssp dssp( pose );
	const FoldTree& tree = pose.fold_tree();
	for ( size_t ir = ssElements[element].first; ir <= ssElements[element].second; ++ir ) {
		Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
		char ss1 = dssp.get_dssp_secstruct( ir );
		char aa1 = pose.residue(ir).name1();
		for ( size_t jr = 1; jr <= pose.size(); ++jr ) {
			if ( !tree.is_jump_point(ir) && !tree.is_jump_point(jr) ) {
				Real dist = pose.residue(ir).xyz("CA").distance(pose.residue(jr).xyz("CA"));
				if ( dist < 12 ) {
					char ss2 = dssp.get_dssp_secstruct( jr );
					char aa2 = pose.residue(jr).name1();
					Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
					Xform const Xbb = ibb_stub.inverse() * jbb_stub;
					core::scoring::motif::XformScoreCOP xs_bb_fxn1(mman_->get_xform_score_BB_BB(ss1,ss2,aa1,aa2));
					core::scoring::motif::XformScoreCOP xs_bb_fxn2(mman_->get_xform_score_BB_BB(ss2,ss1,aa2,aa1));
					Real tmpScore = 0;
					if ( xs_bb_fxn1 != NULL ) {
						tmpScore += xs_bb_fxn1->score_of_bin(Xbb);
						tmpScore += xs_bb_fxn2->score_of_bin(Xbb.inverse());
					}
					if ( tmpScore>0 ) {
						Size resHit = which_ssElement(jr,ssElements);
						if ( resHit != 0 && resHit !=element ) { //must be in SS element and not the current element
							ssElements_in_contact.insert(resHit);
						}
					}
				}
			}
		}
	}
	return(ssElements_in_contact.size());
}

//I have implemented this as a WholeStructureEnergy because the DSSP calls would waste time. But it may be useful in the future to develop this term over each residue
void SSElementMotifContactEnergy::finalize_total_energy( pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	vector1<std::pair<Size,Size> > ss_elements = get_ss_elements(pose);
	Size ignore_terminal_ss = option[OptionKeys::score::ignore_terminal_ss_elements]();
	Size startElement = 1;
	Size endElement = ss_elements.size();
	if ( ignore_terminal_ss>0 ) {
		startElement+=ignore_terminal_ss;
		endElement-=ignore_terminal_ss;
	}
	Real tmpScore = 9999;
	for ( Size ii=startElement; ii<=endElement; ++ii ) {
		Real tmp=get_SSelements_in_contact(ii,ss_elements,pose);
		if ( tmp<tmpScore ) {
			tmpScore=tmp;
		}
	}
	if ( tmpScore>option[OptionKeys::score::max_contacting_ss]() ) {
		tmpScore=option[OptionKeys::score::max_contacting_ss]();
	}
	totals[ss_contact_worst] =(-1*tmpScore);
}


core::Size
SSElementMotifContactEnergy::version() const
{
	return 1; // Initial versioning
}

}//motif
}//scoring
}//core
