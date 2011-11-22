// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/chemical/sdf/ctab_typer.cc
///
/// @brief
/// @author Robert Carroll

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <core/chemical/sdf/ctab_typer.hh>
#include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
//#include <core/chemical/ResidueType.hh>
//#include <protocols/ligand_docking/ColoredGraph.hh>
#include <utility/exit.hh>
#include <numeric/xyzVector.hh>
#include <vector>
#include <algorithm>
#include <map>
// AUTO-REMOVED #include <sstream>
// AUTO-REMOVED #include <cstring>
// AUTO-REMOVED #include <stdlib.h>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


namespace core {
namespace chemical {
namespace sdf {



std::string atomTyper::getType()
{
	if(type_=="")
	{
		//Get the default atom type
		type_=molecule_container_->atom_type(atomno_).name();
		core::Size nonO_aro_bonds = get_bond_count("C",core::chemical::AromaticBond)+
						get_bond_count("N",core::chemical::AromaticBond)+
						get_bond_count("S",core::chemical::AromaticBond);

		core::Size H_single_bonds = get_bond_count("H",core::chemical::SingleBond);

		core::Size O_single_bonds = get_bond_count("O",core::chemical::SingleBond);
		core::Size O_double_bonds = get_bond_count("O",core::chemical::DoubleBond);

		core::Size N_single_bonds = get_bond_count("N",core::chemical::SingleBond);
		core::Size N_double_bonds = get_bond_count("N",core::chemical::DoubleBond);

		core::Size C_single_bonds = get_bond_count("C",core::chemical::SingleBond);
		core::Size C_aromatic_bonds = get_bond_count("C",core::chemical::AromaticBond);

		if(element_=="C")
		{

			if(nonO_aro_bonds>0) {
				type_="aroC";
			}else if(O_single_bonds > 0 && O_double_bonds > 0) {
				type_="COO";
			}else if(N_single_bonds == 1) {
				if(O_double_bonds == 1 && C_single_bonds == 1 && hasbbN()) {
					type_="CObb";
				} else if(O_double_bonds == 1 && hasCarbonylC()) {
					type_="CAbb";
				} else {
					type_="CNH2";
				}
			}else if(N_single_bonds == 2 && N_double_bonds == 1) {
				type_="aroC"; //Guanidinium
			}else {
				switch(H_single_bonds){
				case 1: type_="CH1";break;
				case 2: type_="CH2";break;
				case 3: type_="CH3";break;
				}
			}
		} else if (element_ == "N")
		{
			if(nonO_aro_bonds==2) {
				if(hasNinRing()) {
					type_="Nhis";
				} else {
					type_="Ntrp";
				}
			} else if(hasGuanidiniumC()) {
				type_="Narg";
			} else if(molecule_container_->atomic_charge(atomno_) == 1) {
				type_="Nlys";
			} else if(hasCarbonylC()){
				if(H_single_bonds == 2) {
					type_="NH2O";
				} else if (H_single_bonds == 1 && C_single_bonds == 1) {
					type_="Nbb";
				} else if (C_single_bonds == 2) {
					type_="Npro";
				}
			}
		} else if (element_ == "O")
		{
			if(nonO_aro_bonds > 0) {
				type_="Oaro";
			} else if(H_single_bonds == 2) {
				type_="HOH";
			} else if(H_single_bonds == 1 && C_single_bonds == 1) {
				type_="OH";
			}else if(hasCarboxylC()) {
				type_="OOC";
			} else if(hasCwithbbN()) {
				type_="OCbb";
			} else if(hasAmideN()) {
				type_="ONH2";
			}
		} else if (element_ == "H")
		{

			core::Size NOS_bonds = get_bondelement_count("N")+get_bondelement_count("O")+get_bondelement_count("S");

			if(NOS_bonds >= 1)
			{
				type_="Hpol";
			}else if(C_aromatic_bonds >=1)
			{
	 			type_="Haro";
			}else
			{
				type_="Hapo";
			}

			/*
			if(C_single_bonds == 1) {
				type_="Hapo";
			} else {
				type_="Hpol";
			}
			*/
		} else if (element_ == "Fe")
		{
			//Read the charge, and change to +2/+3 if applicable.
			if(molecule_container_->atomic_charge(atomno_)==2)
			{
				molecule_container_->set_atom_type(atomname_,"Fe2p");
			}
			else if(molecule_container_->atomic_charge(atomno_)==3)
			{
				molecule_container_->set_atom_type(atomname_,"Fe3p");
			}
		} else
		{
			//It is properly set already
		}
	}
	return type_;
}

atomTyper::atomTyper(core::Size atomno, core::chemical::ResidueTypeOP &molecule_container_)
{
	this->atomno_=atomno;
	this->molecule_container_=molecule_container_;
	this->atomname_=molecule_container_->atom_name(atomno);
	element_ = molecule_container_->atom_type(atomno).element();


	type_="";
	/*
	for(core::Size i=0; i<14; i++) {
		numBonds_[i]=0;
	}
	*/
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size(); i++) {
		//typerBond bond;
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		core::chemical::BondName  bondtype = bondtypes[i];
		core::Size bond_count = get_bond_count(element,bondtype);
		set_bond_count(element,bondtype,bond_count+1);
		//std::cout <<element << " " << bondtype << " "<<get_bond_count(element,bondtype) <<std::endl;
	}
}

bool atomTyper::hasbbN()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="N") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("C",core::chemical::SingleBond)==2) {
				found=true;
			}
		}
	}
	return found;
}

std::string atomTyper::get_element()
{
	return element_;
}

bool atomTyper::hasCarbonylC()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="C") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("O",core::chemical::DoubleBond )==1) {
				found=true;
			}
		}
	}
	return found;
}

bool atomTyper::hasNinRing()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="C") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("N",core::chemical::SingleBond)+
				neighbor.get_bond_count("N",core::chemical::DoubleBond)+
				neighbor.get_bond_count("N",core::chemical::AromaticBond)>=2) {
				found=true;
			}
		}
	}
	return found;
}

bool atomTyper::hasGuanidiniumC()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="C") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("N",core::chemical::SingleBond)+
				neighbor.get_bond_count("N",core::chemical::DoubleBond)+
				neighbor.get_bond_count("N",core::chemical::AromaticBond)>=3) {
				found=true;
			}
		}
	}
	return found;
}

bool atomTyper::hasAmideN()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="C") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("N",core::chemical::SingleBond) >= 1) {
				found=true;
			}
		}
	}
	return found;
}

bool atomTyper::hasCarboxylC()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="C") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("O",core::chemical::DoubleBond) == 1 &&
				neighbor.get_bond_count("O",core::chemical::SingleBond) == 1) {
				found=true;
			}
		}
	}
	return found;
}

bool atomTyper::hasCwithbbN()
{
	bool found=false;
	//Bonded to this atom
	core::chemical::AtomIndices bonded = molecule_container_->bonded_neighbor(atomno_);
	utility::vector1<core::chemical::BondName> bondtypes = molecule_container_->bonded_neighbor_types(atomno_);

	//Save the C with an N
	core::Size iC=0;

	//Both are vector1
	for(core::Size i=1; i<=bonded.size()&&!found; i++) {
		std::string element = molecule_container_->atom_type(bonded[i]).element();
		if(element=="C") {
			atomTyper neighbor(bonded[i],molecule_container_);
			if(neighbor.get_bond_count("N",core::chemical::SingleBond)>=1) {
				found=true;
				iC=bonded[i];
			}
		}
	}
	if(found) {
		found=false;
		bonded = molecule_container_->bonded_neighbor(iC);
		bondtypes = molecule_container_->bonded_neighbor_types(iC);
		//Both are vector1
			for(core::Size i=1; i<=bonded.size()&&!found; i++) {
				std::string element = molecule_container_->atom_type(bonded[i]).element();
				if(element=="N") {
					atomTyper neighbor(bonded[i],molecule_container_);
					if(neighbor.get_bond_count("C",core::chemical::SingleBond)>=2) {
						found=true;
					}
				}
			}
	}
	return found;
}

/*
core::Size atomTyper::getNumBonds(typerBond i)
{
	return numBonds_[i];
}
*/

core::Size atomTyper::getNumBonds()
{
	core::Real total=0;
	total += 1*(get_bond_count("C",core::chemical::SingleBond)+get_bond_count("N",core::chemical::SingleBond)+
			get_bond_count("H",core::chemical::SingleBond));
	//total+=1*(numBonds_[Csingle]+numBonds_[Nsingle]+numBonds_[Osingle]+numBonds_[Hsingle]);
	total += 2*(get_bond_count("C",core::chemical::DoubleBond)+get_bond_count("N",core::chemical::DoubleBond)+
			get_bond_count("O",core::chemical::DoubleBond));
	//total+=2*(numBonds_[Cdouble]+numBonds_[Ndouble]+numBonds_[Odouble]);
	total += 3*(get_bond_count("C",core::chemical::TripleBond)+get_bond_count("N",core::chemical::TripleBond));
	//total+=3*(numBonds_[Ctriple]+numBonds_[Ntriple]);
	total+=1.5*(get_bondtype_count(core::chemical::AromaticBond));

	return static_cast<core::Size>(total);
}

void atomTyper::set_bond_count(std::string const element, core::chemical::BondName const bond_type,core::Size const bond_count)
{
	std::pair<std::string,core::chemical::BondName> key(element,bond_type);
	numBonds_[key] = bond_count;
}

core::Size atomTyper::get_bond_count(std::string const element, core::chemical::BondName const bond_type) const
{
	std::pair<std::string,core::chemical::BondName> key(element,bond_type);
	std::map< std::pair<std::string,core::chemical::BondName>,core::Size>::const_iterator bond_iter;
	bond_iter = numBonds_.find(key);
	if(bond_iter == numBonds_.end())
	{
		return 0;
	}else
	{
		return bond_iter->second;
	}
}

core::Size atomTyper::get_bondtype_count(core::chemical::BondName const bond_type) const
{
	core::Size total = 0;
	std::map< std::pair<std::string,core::chemical::BondName>,core::Size >::const_iterator bond_iter;
	for(bond_iter = numBonds_.begin();bond_iter != numBonds_.end();++bond_iter)
	{
		std::pair<std::string,core::chemical::BondName> key = bond_iter->first;
		if(key.second == bond_type)
		{
			total += bond_iter->second;
		}
	}
	return total;
}

core::Size atomTyper::get_bondelement_count(std::string const element) const
{
	core::Size total = 0;
	std::map< std::pair<std::string,core::chemical::BondName>,core::Size >::const_iterator bond_iter;
	for(bond_iter = numBonds_.begin();bond_iter != numBonds_.end();++bond_iter)
	{
		std::pair<std::string, core::chemical::BondName> key = bond_iter->first;
		if(key.first == element)
		{
			total += bond_iter->second;
		}
	}
	return total;
}

} // sdf
} // io
} // core
