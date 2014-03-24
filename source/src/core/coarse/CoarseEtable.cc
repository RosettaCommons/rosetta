// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/CoarseEtable.cc
/// @brief
/// @author Oliver Lange

// Unit headers
#include <core/coarse/CoarseEtable.hh>

/// Project Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <basic/database/open.hh>

#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
#include <fstream>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>



using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace coarse;
using namespace std;

CoarseEtable::CoarseEtable( chemical::AtomTypeSetCAP atom_set, std::string tag )
	: atom_set_(atom_set) {
	if ( tag == chemical::COARSE_TWO_BEAD )
		{
			read_files(
		 basic::database::full_name( "scoring/score_functions/etable/resolve_etable.twobead.txt" ),
		 basic::database::full_name( "scoring/score_functions/etable/etable.twobead.lj.dat" ),
		 basic::database::full_name( "etable.twobead.dlj.dat" )
			);
		} else {
		std::string msg = "do not know how to load etable "+tag;
		utility_exit_with_message( msg );
	}
}


void CoarseEtable::check_atomset_compatibility(chemical::AtomTypeSetCAP normal) const {
	using namespace std;
	chemical::AtomTypeSetCAP coarse = atom_set_; //our own copy of AtomSetCOP

	for (Size i=1;i<=normal->n_atomtypes();i++) {
		string name = (*normal)[i].name();
		Size ii=coarse->atom_type_index(name);
		if (i!=ii) {
			cerr << " --------------- WARNING --------------------- " << endl;
			cerr << " incompatible atomtype indexes -- reorder coarse atomset or change program " << endl;
			cerr << " if you go on, the wrong entries from the normal full-atom etable might be used " << endl;
			cerr << " ------WARNING-----------------WARNING---------" << endl;
		}
		cerr << i << ' ' << ii << ' ' << name << endl;
	}
}

/// @brief setup before atom_pair functions can be called
void CoarseEtable::set_residue_pair(conformation::Residue const &rsd1,conformation::Residue const &rsd2) const {
	seq_dist_=rsd1.polymeric_sequence_distance(rsd2);
}


void CoarseEtable::print_residue_info(conformation::Residue const &rsd1,conformation::Residue const &rsd2) const {
		T("coarse.scoring") << "prepare coarse scoring for " << rsd1.type().name() << " " << rsd2.type().name() << "\n";


}

void CoarseEtable::dump_oldstyle_type_table(std::ostream &os,const chemical::ResidueTypeSet& rsd_set) {
	for (Size i=1;i<=atom_set_->n_atomtypes();++i) {
		os << i << ' ' << (*atom_set_)[i].name() << endl;
	}
	const int nBeads = 2;
	for (int aa = 1; aa<=20; aa++) {
		for (int ibead = 1;ibead<=nBeads;ibead++) {
			os  << -((aa-1)*nBeads+ibead) << ' '
		<< rsd_set.aa_map(static_cast<chemical::AA>(aa)).front()->name3() << "_B" << ibead+1 << endl;
		}
	}
}

struct Entry {
	int type[2];
	int seq_dist_;
	int ID;
};
typedef  std::vector< Entry > Entries;

void CoarseEtable::read_files(std::string fn_resolve,std::string fn_etable, std::string fn_dtable) {

 //the tables have this length
	const int TABLE_LENGTH ( 721 );

	std::ifstream resolve( fn_resolve.c_str() ); // change to izstream!
	// parse the header line
	Entries entries;
	maxID = 0;
	maxType = 0;
	maxDist = 0;
	{ // scope
		std::string line, tag, tag2;
		while (getline( resolve, line )) {
			std::istringstream l( line );
			Entry entry;

			//read entry from line
			for (int i=0;i<2;i++) {
				l >> tag;
				if ( l.fail() ) utility_exit_with_message("CoarseEtable::read_file: bad line: "+  line );
				entry.type[i] = atom_set_->atom_type_index( tag ); //fails if type not available
				maxType = std::max(maxType,entry.type[i]);
			};
			l >> entry.seq_dist_;
			if ( l.fail() ) utility_exit_with_message("CoarseEtable::read_file: bad line: "+  line );
			l >> entry.ID;
			if ( l.fail() ) utility_exit_with_message("CoarseEtable::read_file: bad line: "+  line );
			entries.push_back(entry);
			//      cerr << "entry " << entry.type[0] << ' ' << entry.type[1] << ' ' << entry.seq_dist_ << ' ' << entry.ID << endl;
			maxID = std::max(maxID,entry.ID);
			maxDist= std::max(maxDist,entry.seq_dist_);
		};
		cerr << "creating resolve_table with " << maxType << ' ' << maxDist+1 << endl;
		resolve_.dimension(maxType,maxType,maxDist+1);
		resolve_ = 0;
		cerr << "fill resolve_table " << endl;
		for (Entries::const_iterator it=entries.begin(), eit=entries.end(); it!=eit; ++it) {
			resolve_(it->type[0],it->type[1],it->seq_dist_+1)=it->ID;
			resolve_(it->type[1],it->type[0],it->seq_dist_+1)=it->ID;

			//      std::cout << "resolve_( " << it->type[ 0 ] << ", " << it->type[ 1 ] << ", " << it->seq_dist_ + 1 << ") = " << it->ID << std::endl;
		}

		cerr << "fill missing entries in resolve table ... " << endl;
		//fix non-specific seq_dist_ances  -- replace by default behaviour
		for (int ti=1;ti<=maxType;ti++) {
			for (int tj=1;tj<=maxType;tj++) {
	for (int d=1;d<=maxDist;d++) {
		//  cerr << ti << ' ' << tj << ' ' << d << endl;
		if (resolve_(ti,tj,d+1)==0) resolve_(ti,tj,d+1)=resolve_(ti,tj,1);
	};
			}
		}

		etable_.dimension(TABLE_LENGTH,maxID);
		dtable_.dimension(TABLE_LENGTH,maxID);
		std::ifstream eet( fn_etable.c_str() ); // change to izstream!
		std::ifstream det( fn_dtable.c_str() ); // change to izstream!

		cerr << "read etable data " << endl;
		for (int i=1;i<=maxID;i++) {
			for (int k=1;k<=TABLE_LENGTH;k++) {
	//	cerr << i << ' ' << k << ' '  << endl;
	eet >> etable_(k,i);
	det >> dtable_(k,i);
			}
		}
#if 0
		cerr << "check table " << maxID << endl;
		for (int k=1;k<TABLE_LENGTH;k++) {
			cerr << etable_(5,k) << ' ';
		}
		cerr << endl;
#endif

	}
}



bool
CoarseEtable::atom_pair_energy(
						 int disbin,
						 Real frac,
						 conformation::Atom const &atom1,
						 conformation::Atom const &atom2,
						 core::Energy &bb
) const
{
	int const eID = get_eID(atom1,atom2,seq_dist_);
	bb = 0.0;
	if (eID>0) {
		//    std::cerr << __FILE__<< ' ' << __LINE__ << ' ' << disbin << ' ' << eID << std::endl;
		int const l1 = etable_.index( disbin, eID ), l2 = l1 + 1;
		Real e1 = etable_[ l1 ];
		bb = ( e1 + frac * ( etable_[ l2 ] - e1 ) );
		return true;
	}
	return false;
}

Real CoarseEtable::eval_dE_dR(
						int disbin,
						Real frac,
						conformation::Atom const &atom1,
						conformation::Atom const &atom2,
						scoring::EnergyMap const &weights
) const
{
	Real deriv = 0.0;
	int const eID = get_eID(atom1,atom2,seq_dist_);
	//std::cout << "coarse dE_dR requested for " << (*atom_set_)[atom1.type()].name() << ' '<< (*atom_set_)[atom2.type()].name()  << ' ' << eID << std::endl;

	if (eID>0) {
		int const l1 = dtable_.index( disbin, eID), l2 = l1 + 1;
		Real e1 = dtable_[ l1 ];
		deriv = weights[ scoring::coarse_beadlj ] * ( e1 + frac * ( dtable_[ l2 ] - e1 ) );
	}
	return deriv;
}




