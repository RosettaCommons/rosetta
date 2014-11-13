// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_scoring_motif_motif_hash_stuff_FWD_hh
#define INCLUDED_core_scoring_motif_motif_hash_stuff_FWD_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace motif {


class ResPairMotif;
class Xfres;
class Xfrag;

class XfragSet;
typedef utility::pointer::shared_ptr< XfragSet       > XfragSetOP;
typedef utility::pointer::shared_ptr< XfragSet const > XfragSetCOP;
typedef utility::pointer::weak_ptr< XfragSet const > XfragSetCAP;

class XformScore;
typedef utility::pointer::shared_ptr< XformScore       > XformScoreOP;
typedef utility::pointer::shared_ptr< XformScore const > XformScoreCOP;
typedef utility::pointer::weak_ptr< XformScore const > XformScoreCAP;

class MotifHash;
typedef utility::pointer::shared_ptr< MotifHash       > MotifHashOP;
typedef utility::pointer::shared_ptr< MotifHash const > MotifHashCOP;
typedef utility::pointer::weak_ptr< MotifHash const > MotifHashCAP;

class MotifHits;
class ResPairMotifQuery;

enum RM_Type {
	RM_Type_NONE,
	RM_SC,
	RM_BB,
	RM_PH,
	RM_PO,
	NUM_RM_TYPES=RM_PO
};
enum RPM_Type {
	RPM_Type_NONE,
	SC_SC,
	SC_BB,
	SC_PH,
	SC_PO,
	BB_BB,
	BB_PH,
	BB_PO,
	PH_PO,
	NUM_RPM_TYPES=PH_PO
};


}
}
}

#endif

