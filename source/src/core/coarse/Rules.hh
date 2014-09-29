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

#ifndef INCLUDED_core_coarse_Rules_hh
#define INCLUDED_core_coarse_Rules_hh

// Unit headers
#include <core/coarse/Rules.fwd.hh>
// you cannot #include yourself #include <core/coarse/Rules.hh>


// Project headers
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>


// std headers
#include <iostream>

//Auto Headers
#include <utility/vector1_bool.hh>
#include <map>


namespace core {
namespace coarse {
// class Rule;
// typedef utility::pointer::owning_ptr< Rule > RuleOP;
// typedef utility::pointer::owning_ptr< Rule const > RuleCOP;
// typedef utility::pointer::access_ptr< Rule > RuleAP;
// typedef utility::pointer::access_ptr< Rule const > RuleCAP;


// class Rule  : public utility::pointer::ReferenceCount {
// public:
// 	typedef std::string AtomToken;
// 	typedef std::vector<std::string> Tokens;
// 	typedef std::string BeadName;
// 	typedef std::map<BeadName, Tokens> BeadRules;
// 	typedef BeadRules::const_iterator ConstBeadIterator;
// 	typedef Tokens::const_iterator ConstTokenIterator;

// 	class BeadAdder {
// 	public:
// 		BeadAdder(Tokens& tokens) : tokens_(tokens) {};
// 		BeadAdder& operator << (std::string token) {tokens_.push_back(token); return *this;};
// 	private:
// 		Tokens& tokens_;
// 	};

// public:
// 	Rule() : my_AA_(chemical::aa_unk) {};
// 	Rule(chemical::AA _AA) : my_AA_(_AA) {};
// 	void pretty_print(std::ostream &os);
// 	static const BeadName FULL_ATOM;
// 	static const AtomToken REST_SIDECHAIN;
// 	static const AtomToken REST_ALL;
// 	//      void add_token(std::string bead,std::string token);
// 	BeadAdder add_to_bead(BeadName bead) {return BeadAdder(rules_[bead]); };
// 	ConstBeadIterator begin() const {return rules_.begin(); };
// 	ConstBeadIterator end() const {return rules_.end(); }
// private:
// 	chemical::AA my_AA_;
// 	BeadRules rules_;
// };

// class GenericRule : public Rule {
// 	// the generic rule: momentarily in doubt all atoms remain full-atom
// public:
// 	GenericRule(chemical::AA _AA) : Rule(_AA) {
// 		add_to_bead(Rule::FULL_ATOM) << Rule::REST_ALL;
// 	};
// };

// class RuleSet : public utility::pointer::ReferenceCount {
// 	typedef std::map<chemical::AA,RuleOP> RuleMap;
// 	// might become AA-> vector<Rule> if we have more than 1 possible rule per AA-type
// public:
// 	RuleSet() {};
// 	RuleSet(std::string tag);
// 	void create_rules();
// 	void pretty_print(std::ostream &os);
// 	//void read_rules_from_file(std::istream&) {};
// 	bool has(chemical::AA aa) const;
// 	RuleCOP operator[] (chemical::AA aa) const  {
// 		if (has(aa))
// 			return rules_[aa];
// 		else {
// 			std::cerr << "WARNING: no rule for AA type " << aa << std::endl << "returning generic rule" << std::endl;
// 			return RuleCOP(new GenericRule(aa) );
// 		}
// 	}

// private:
// 	mutable RuleMap rules_; //because is not a const function operator[]

// };

} //namespace coarse
} // namespace core
#endif
