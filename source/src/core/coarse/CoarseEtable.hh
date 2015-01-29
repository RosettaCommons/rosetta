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

#ifndef INCLUDED_core_coarse_CoarseEtable_hh
#define INCLUDED_core_coarse_CoarseEtable_hh


// Project headers
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// unit headers
#include <core/coarse/CoarseEtable.fwd.hh>

// Objexx Headers
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>

// std headers
// AUTO-REMOVED #include <ostream>


// AUTO-REMOVED #include <basic/Tracer.hh>

//Auto Headers
#include <core/chemical/AtomTypeSet.fwd.hh>





namespace core {
namespace coarse {
// Unit headers

class CoarseEtable : public utility::pointer::ReferenceCount {
public:
	CoarseEtable(chemical::AtomTypeSetCAP atom_set, std::string tag);
	void dump_oldstyle_type_table(std::ostream &os,const chemical::ResidueTypeSet&);
	void read_files(std::string resolve, std::string etable, std::string dtable);

	/// @brief setup before atom_pair functions can be called
	void set_residue_pair(conformation::Residue const &rsd1,conformation::Residue const &rsd2) const;

	void print_residue_info(conformation::Residue const &rsd1,conformation::Residue const &rsd2) const;

	bool handles(conformation::Atom const &atom1, conformation::Atom const &atom2) const
	{ return get_eID(atom1,atom2,seq_dist_)>0; }

	/// @brief atom_pair evaluations
	/// they return false if atom_pair is not evaluated by coarse table --> use normal etable instead
	/// WARNING, atom_pair evaluations require to call set_residue_pair first !
	/// Okay that is a bit dirty, but otherwise we need to change quite a lot of the mini-scoring methods
	/// --- do that when this stuff is established
	bool atom_pair_energy(
		int disbin,
		Real frac,
		conformation::Atom const &atom1,
		conformation::Atom const &atom2,
		Energy &bb
	) const;

	Real eval_dE_dR(
		int disbin,
		Real frac,
		conformation::Atom const &atom1,
		conformation::Atom const &atom2,
		const scoring::EnergyMap &weights
	) const;

	void check_atomset_compatibility(chemical::AtomTypeSetCAP normal) const;

	chemical::AtomTypeSetCAP
	atom_set() const
	{
		return atom_set_;
	}


private:

	/// @brief get entry index in table
	int get_eID(
		conformation::Atom const &atom1,
		conformation::Atom const &atom2,
		int seq_dist
	 ) const
	{
		if ((atom1.type()<=maxType) && (atom2.type()<=maxType)) {
			int const dist = seq_dist > maxDist ? 0 : seq_dist;
		debug_assert(dist>=0);

			/*						basic::T("coarse.scoring") << "coarse score requested for "
																<< (*atom_set_)[atom1.type()].name() << ' '
																<< (*atom_set_)[atom2.type()].name() << ' '
																			<< dist << " ( " << maxDist << " ) " << std::endl;
						basic::T("coarse.scoring") << resolve_(atom1.type(),atom2.type(),dist+1) <<	std::endl;
			*/
			return resolve_(atom1.type(),atom2.type(),dist+1);
		} else return 0;
	}


public:
	chemical::AtomTypeSetCAP atom_set_;
private:
	ObjexxFCL::FArray3D< uint > resolve_; // find entry in energy table
	ObjexxFCL::FArray2D< Real > etable_;  // energies
	ObjexxFCL::FArray2D< Real > dtable_;  // derivatives
	int maxID; int maxType; int maxDist;
	mutable int seq_dist_;

};



} // namespace coarse
} // namespace core

#endif
