// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/chemical/sdf/ctab_typer.hh
///
/// @brief determine the type of an atom
/// @author Robert Carroll


#ifndef INCLUDED_core_chemical_sdf_ctab_typer_HH
#define INCLUDED_core_chemical_sdf_ctab_typer_HH

#include <string>
#include <map>
#include <vector>
#include <core/chemical/sdf/ctab_typer.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueType.fwd.hh>



namespace core {
namespace chemical {
namespace sdf {
enum typerBond { Csingle=0,
				Cdouble=1,
				Ctriple=2,
				Caro=3,
				Nsingle=4,
				Ndouble=5,
				Ntriple=6,
				Naro=7,
				Osingle=8,
				Odouble=9,
				Oaro=11,
				Hsingle=12,
				AroBonds=13};

class atomTyper {
public:
	atomTyper(core::Size atomno, core::chemical::ResidueTypeOP &molecule_container_);
	std::string getType();
	//core::Size getNumBonds(typerBond i);
	std::string get_element();
	core::Size getNumBonds();

	core::Size get_bond_count(std::string const element, core::chemical::BondName const bond_type) const;

	core::Size get_bondtype_count(core::chemical::BondName const bond_type) const;

	core::Size get_bondelement_count(std::string const element) const;

private:
	core::Size atomno_;
	std::string atomname_;
	core::chemical::ResidueTypeOP molecule_container_;
	std::string type_;
	std::map< std::pair<std::string,core::chemical::BondName>,core::Size> numBonds_;
	//core::Size numBonds_[14];
	std::string element_;

	void set_bond_count(std::string const element, core::chemical::BondName const bond_type,core::Size const bond_count);



	bool hasbbN(); //For CObb (bb N means 2x C bonds)
	bool hasCarbonylC(); //for Cabb and NH2O, Nbb, Npro
	bool hasNinRing(); //for Nhis
	bool hasGuanidiniumC(); //for Narg
	bool hasAmideN(); //For ONH2
	bool hasCarboxylC(); //For OOC
	bool hasCwithbbN(); //for OCbb
};


}
}
}


#endif /* CTAB_TYPER_HH_ */
