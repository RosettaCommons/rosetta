// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Oliver Lange



#include "core/coarse/Rules.hh"
#ifdef WIN32
#include "core/pose/Pose.hh"
#include "core/conformation/ResidueFactory.hh"
#endif
#include "utility/assert.hh"

//Auto Headers
#include <core/conformation/Residue.fwd.hh>




using namespace core;
using namespace coarse;
using namespace conformation;
using namespace chemical;
using namespace std;


const std::string Rule::FULL_ATOM = "FULL_ATOM";
const std::string Rule::REST_SIDECHAIN = "REST_SC";
const std::string Rule::REST_ALL = "REST";


void Rule::pretty_print(std::ostream &os) {
	os << "Rule for " << my_AA_ << " : "<< endl; // should use name3(my_AA) if it becomes available
	for (ConstBeadIterator it=begin(), eit=end(); it!=eit; ++it) {
		os << "coarse atom " << it->first <<" { ";
		for (ConstTokenIterator tit=it->second.begin(), etit=it->second.end(); tit!=etit; ++tit) {
			os << *tit << ", ";
		};
		os << " }" << endl;
	};
}

RuleSet::RuleSet(std::string ASSERT_ONLY(tag) ) {
	//depending on tag we could have different rulesets or load them from file...
	// now we have only the soup-of-the-day
debug_assert( tag == "coarse_two_bead" );
	create_rules();
}

void RuleSet::create_rules() {
	const int n_restype = aa_tyr;

	/* first of all add Cbeta to bead 1 for all aminoacids apart from aa_gly */
	for (int iaa=aa_ala; iaa<=n_restype;++iaa) {
		AA aa=static_cast<AA>(iaa);
		Rule * pcr = new Rule(aa);
		RuleOP cr(pcr);
		rules_[aa]=cr;
		if (aa!=aa_gly) {
			rules_[aa]->add_to_bead("B1") << "CB";
		}
	};
	/*hydrogens are added automatically...atoms which are not available will be ignored */
	/* add other atoms to bead1 manually -- there are not many */
	rules_[aa_ile]->add_to_bead("B1") << "CG2";
	rules_[aa_arg]->add_to_bead("B1") << "CG" << "CD";
	rules_[aa_glu]->add_to_bead("B1") << "CG";
	rules_[aa_gln]->add_to_bead("B1") << "CG";
	rules_[aa_lys]->add_to_bead("B1") << "CG";
	rules_[aa_pro]->add_to_bead("B1") << "CG" << "CD";
	rules_[aa_met]->add_to_bead("B1") << "CG";
	for (int iaa=aa_ala;iaa<=n_restype;iaa++) {
		AA aa=static_cast<AA>(iaa);
		//add everything that hasn't been bound to bead1
		if (aa!=aa_gly && aa!=aa_ala && aa!=aa_pro) {
			rules_[aa]->add_to_bead("B2") << Rule::REST_SIDECHAIN;
		}
		rules_[aa]->add_to_bead(Rule::FULL_ATOM) << Rule::REST_ALL;
	};

}

bool RuleSet::has(chemical::AA aa) const {
	return rules_.find(aa)!=rules_.end();
}

void RuleSet::pretty_print (ostream &os) {
	typedef RuleMap::iterator MIT;
	for (MIT it=rules_.begin(),eit=rules_.end(); it!=eit; ++it) {
		it->second->pretty_print(os);
		os << endl;
	};
}



