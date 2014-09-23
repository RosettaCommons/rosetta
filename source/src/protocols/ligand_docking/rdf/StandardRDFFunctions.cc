// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/rdf/StandardRDFFunctions.hh
///
/// @brief headers for standard RDF Functions used by RosettaHTS
/// @author Sam DeLuca


#include <core/pose/Pose.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/orbitals/OrbitalsScore.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <protocols/ligand_docking/rdf/StandardRDFFunctions.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>


namespace protocols {
namespace ligand_docking {
namespace rdf {


RDFBaseOP RDFEtableCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFEtableFunction );
}

std::string RDFEtableCreator::type_name() const
{
	return "RDFEtableFunction";
}

RDFEtableFunction::RDFEtableFunction() : RDFBase("RDFEtableFunction"), etable_evaluator_(/* NULL */)
{
	this->add_function_name("atr");
	this->add_function_name("rep");
	this->add_function_name("solv");
}

RDFEtableFunction::~RDFEtableFunction()
{

}

void RDFEtableFunction::parse_my_tag(utility::tag::TagCOP tag, basic::datacache::DataMap & data_map)
{
	if(!tag->hasOption("scorefxn"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'RDFEtableFunction' requires 'scorefxn' tag");
	}

	std::string scorefxn_name = tag->getOption<std::string>("scorefxn");
	core::scoring::ScoreFunctionOP scorefxn(data_map.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_name));
	core::scoring::methods::EnergyMethodOptions options(scorefxn->energy_method_options());
	core::scoring::etable::EtableCOP etable(core::scoring::ScoringManager::get_instance()->etable( options.etable_type()));
	etable_evaluator_ = core::scoring::etable::AnalyticEtableEvaluatorOP( new core::scoring::etable::AnalyticEtableEvaluator(*etable) );

}

RDFResultList RDFEtableFunction::operator()(AtomPairData const & atom_data )
{
	core::Real weight = 1.0;
	core::Real atr_energy = 0.0;
	core::Real rep_energy = 0.0;
	core::Real solv_energy = 0.0;
	core::Real d2 = 0.0;

	etable_evaluator_->atom_pair_energy(
		atom_data.protein_atom,
		atom_data.ligand_atom,
		weight,
		atr_energy,
		rep_energy,
		solv_energy,
		d2);

	RDFResultList results;
	results.push_back(std::make_pair("atr",atr_energy));
	results.push_back(std::make_pair("rep", rep_energy));
	results.push_back(std::make_pair("solv", solv_energy));
	return results;
}

RDFBaseOP RDFElecCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFElecFunction );
}

std::string RDFElecCreator::type_name() const
{
	return "RDFElecFunction";
}

RDFElecFunction::RDFElecFunction() : RDFBase("RDFElecFunction"), coloumb_(/* NULL */)
{
	this->add_function_name("elec");
}

RDFElecFunction::~RDFElecFunction()
{

}
void RDFElecFunction::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map)
{
	if(!tag->hasOption("scorefxn"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'RDFEtableFunction' requires 'scorefxn' tag");
	}

	std::string scorefxn_name = tag->getOption<std::string>("scorefxn");
	core::scoring::ScoreFunctionOP scorefxn(data_map.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_name));
	core::scoring::methods::EnergyMethodOptions options(scorefxn->energy_method_options());
	coloumb_ = core::scoring::etable::coulomb::CoulombOP( new core::scoring::etable::coulomb::Coulomb(options) );
	coloumb_->initialize();

}

RDFResultList RDFElecFunction::operator()(AtomPairData const & atom_data )
{
	RDFResultList results;
	core::Real fa_elec = coloumb_->eval_atom_atom_fa_elecE(
		atom_data.protein_atom_coords,
		atom_data.protein_atom_charge,
		atom_data.ligand_atom_coords,
		atom_data.ligand_atom_charge);

	results.push_back(std::make_pair("elec", fa_elec));
	return results;
}

RDFBaseOP RDFChargeCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFChargeFunction );
}

std::string RDFChargeCreator::type_name() const
{
	return "RDFChargeFunction";
}

RDFChargeFunction::RDFChargeFunction() : RDFBase("RDFChargeFunction")
{

}

RDFChargeFunction::~RDFChargeFunction()
{

}
void RDFChargeFunction::parse_my_tag(
								   utility::tag::TagCOP tag,
								   basic::datacache::DataMap &)
{
	//Sign mode can be "ligand_plus" or "ligand_minus" or "same_sign
	std::string sign_mode = tag->getOption<std::string>("sign_mode","same_sign");
	if(sign_mode == "ligand_plus")
	{
		function_sign_ = LigandPlusProteinMinus;
		function_name_ = "charge_plus";

	}else if(sign_mode == "ligand_minus")
	{
		function_sign_ = LigandMinusProteinPlus;
		function_name_ = "charge_minus";
	}else if(sign_mode == "same_sign")
	{
		function_sign_ = SameSign;
		function_name_ = "charge_unsigned";
	}else
	{
		utility::excn::EXCN_RosettaScriptsOption("RDFHbondFunction sign_mode can only be ligand_acceptor or ligand_donor");
	}
	this->add_function_name(function_name_);
}

RDFResultList RDFChargeFunction::operator()(AtomPairData const & atom_data )
{
	RDFResultList results;

	core::Real charge_product = atom_data.protein_atom_charge*atom_data.ligand_atom_charge;
	if(function_sign_ == LigandMinusProteinPlus &&
		atom_data.ligand_atom_charge < 0 &&
		atom_data.protein_atom_charge >= 0
	   )
	{
		results.push_back(std::make_pair(function_name_, charge_product));
	}else if(function_sign_ == LigandPlusProteinMinus &&
		atom_data.ligand_atom_charge >= 0 &&
		atom_data.protein_atom_charge < 0
			 )
	{
		results.push_back(std::make_pair(function_name_, charge_product));
	}else if(function_sign_ == SameSign &&
			 ((atom_data.ligand_atom_charge < 0 &&
			 atom_data.protein_atom_charge < 0) ||
			 (atom_data.ligand_atom_charge >= 0 &&
			 atom_data.protein_atom_charge >= 0)) )
	{
		results.push_back(std::make_pair(function_name_, charge_product));
	}

	return results;
}

RDFBaseOP RDFHbondCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFHbondFunction );
}

std::string RDFHbondCreator::type_name() const
{
	return "RDFHbondFunction";
}

RDFHbondFunction::RDFHbondFunction() : RDFBase("RDFHbondFunction"),hbond_set_(/* NULL */)
{

}

RDFHbondFunction::~RDFHbondFunction()
{

}

void RDFHbondFunction::parse_my_tag(
	utility::tag::TagCOP tag ,
	basic::datacache::DataMap & )
{
	//Sign mode can be "ligand_acceptor" or "ligand_donor"
	std::string sign_mode = tag->getOption<std::string>("sign_mode","unsigned");
	if(sign_mode == "ligand_acceptor")
	{
		function_sign_ = LigandPlusProteinMinus;
		function_name_ = "hbond_acceptor";

	}else if(sign_mode == "ligand_donor")
	{
		function_sign_ = LigandMinusProteinPlus;
		function_name_ = "hbond_donor";
	}else if(sign_mode == "unsigned")
	{
		function_sign_ = SameSign;
		function_name_ = "hbond_unsigned";
	}else
	{
		utility::excn::EXCN_RosettaScriptsOption("RDFHbondFunction sign_mode can only be ligand_acceptor or ligand_donor");
	}
	this->add_function_name(function_name_);
}

void RDFHbondFunction::preamble(core::pose::Pose & pose)
{
	hbond_set_ = core::scoring::hbonds::HBondSetOP( new core::scoring::hbonds::HBondSet );
	pose.update_residue_neighbors();
	hbond_set_->setup_for_residue_pair_energies(pose,false,false);
}

RDFResultList RDFHbondFunction::operator()(AtomPairData const & atom_data )
{
	RDFResultList results;
	//get hbond energy if availible
	utility::vector1< core::scoring::hbonds::HBondCOP > hbond_list(hbond_set_->atom_hbonds(atom_data.protein_atom_id));	//Most atoms will have zero or some small number of hbonds, this is a short loop
	for(core::Size hbond_index = 1; hbond_index <= hbond_list.size();++hbond_index)
	{
		core::scoring::hbonds::HBondCOP current_bond(hbond_list[hbond_index]);
		core::id::AtomID acc_id(current_bond->acc_atm(),current_bond->acc_res());
		core::id::AtomID don_id(current_bond->don_hatm(),current_bond->don_res());

		if((function_sign_ == SameSign) && ((acc_id == atom_data.ligand_atom_id) || (don_id == atom_data.ligand_atom_id)))
		{
			//Unsigned function
			core::Real hbond = current_bond->energy();
			results.push_back(std::make_pair(function_name_, hbond));
			return results;
		}else if( (function_sign_ == LigandPlusProteinMinus) && (acc_id == atom_data.ligand_atom_id) )
		{
			//ligand_acceptor function
			core::Real hbond = current_bond->energy();
			results.push_back(std::make_pair(function_name_, hbond));
			return results;
		}else if( (function_sign_ == LigandMinusProteinPlus) && (don_id == atom_data.ligand_atom_id) )
		{
			//ligand_donor function
			core::Real hbond = current_bond->energy();
			results.push_back(std::make_pair(function_name_, hbond));
			return results;
		}
	}

	results.push_back(std::make_pair(function_name_, 0.0));
	return results;

}

RDFBaseOP RDFBinaryHbondCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFBinaryHbondFunction );
}

std::string RDFBinaryHbondCreator::type_name() const
{
	return "RDFBinaryHbondFunction";
}

RDFBinaryHbondFunction::RDFBinaryHbondFunction() : RDFBase("RDFBinaryHbindFunction")
{

}

RDFBinaryHbondFunction::~RDFBinaryHbondFunction()
{

}

void RDFBinaryHbondFunction::parse_my_tag(utility::tag::TagCOP tag,basic::datacache::DataMap &)
{
	//Sign mode can be "ligand_acceptor" or "ligand_donor", or matching pair
	std::string sign_mode = tag->getOption<std::string>("sign_mode","matching_pair");
	if(sign_mode == "ligand_acceptor")
	{
		function_sign_ = LigandPlusProteinMinus;
		function_name_ = "hbond_binary_acceptor";
	}else if(sign_mode == "ligand_donor")
	{
		function_sign_ = LigandMinusProteinPlus;
		function_name_ = "hbond_binary_donor";

	}else if(sign_mode == "matching_pair")
	{
		function_sign_ = SameSign;
		function_name_ = "hbond_matching_pair";
	}else
	{
		utility::excn::EXCN_RosettaScriptsOption("RDFHbondFunction sign_mode can only be ligand_acceptor or ligand_donor or matching pair");
	}
	this->add_function_name(function_name_);
}

RDFResultList RDFBinaryHbondFunction::operator()(AtomPairData const & atom_data)
{
	RDFResultList results;
	//is the current pair of atoms a donor and acceptor?
	core::Real hbond_binary = 0.0;
	if((function_sign_ == LigandPlusProteinMinus) &&
		( atom_data.protein_atom_type.is_donor() && atom_data.ligand_atom_type.is_acceptor()) )
	{

		hbond_binary = 1.0;
	}else if((function_sign_ == LigandMinusProteinPlus) &&
		(atom_data.protein_atom_type.is_acceptor() && atom_data.ligand_atom_type.is_donor()))
	{
		hbond_binary = 1.0;
	}else if( (function_sign_ == SameSign) &&
		( (atom_data.protein_atom_type.is_donor() && atom_data.ligand_atom_type.is_donor() ) ||
		  (atom_data.protein_atom_type.is_acceptor() && atom_data.ligand_atom_type.is_acceptor()))
		)
	{
		hbond_binary = 1.0;
	}
	else{
		hbond_binary = 0.0;
	}

	results.push_back(std::make_pair(function_name_, hbond_binary));
	return results;
}

RDFBaseOP RDFOrbitalFunctionCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFOrbitalFunction );
}

std::string RDFOrbitalFunctionCreator::type_name() const
{
	return "RDFOrbitalFunction";
}

RDFOrbitalFunction::RDFOrbitalFunction() : RDFBase("RDFOrbitalFunction"),pose_(/* NULL */),orbital_score_(NULL)
{
	this->add_function_name("pci_cation_pi");
	this->add_function_name("pci_pi_pi");
	this->add_function_name("pci_hbond");
	this->add_function_name("orbitals_hpol_bb");
	this->add_function_name("pci_salt_bridge");
}

RDFOrbitalFunction::~RDFOrbitalFunction()
{

}

void RDFOrbitalFunction::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap & )
{

}

RDFResultList RDFOrbitalFunction::operator()(AtomPairData const & atom_data )
{
	core::pose::PoseOP pose( pose_ );

	core::scoring::EnergyMap emap;
	emap.clear();
	orbital_score_->get_E_haro_one_way(*pose,atom_data.protein_atom_id,atom_data.ligand_atom_id,emap);
	orbital_score_->get_E_haro_one_way(*pose,atom_data.ligand_atom_id,atom_data.protein_atom_id,emap);

	orbital_score_->get_E_hpol_one_way(*pose,atom_data.protein_atom_id,atom_data.ligand_atom_id,emap);
	orbital_score_->get_E_hpol_one_way(*pose,atom_data.ligand_atom_id,atom_data.protein_atom_id,emap);

	orbital_score_->get_orb_orb_E(*pose,atom_data.protein_atom_id,atom_data.ligand_atom_id,emap);
	orbital_score_->get_orb_orb_E(*pose,atom_data.ligand_atom_id,atom_data.protein_atom_id,emap);

	RDFResultList results;

	results.push_back(std::make_pair("pci_cation_pi", emap[core::scoring::pci_cation_pi]));
	results.push_back(std::make_pair("pci_pi_pi", emap[core::scoring::pci_pi_pi]));
	results.push_back(std::make_pair("pci_hbond", emap[core::scoring::pci_hbond]));
	results.push_back(std::make_pair("orbitals_hpol_bb", emap[core::scoring::orbitals_hpol_bb]));
	results.push_back(std::make_pair("pci_salt_bridge", emap[core::scoring::pci_salt_bridge]));

	return results;

}

void RDFOrbitalFunction::preamble(core::pose::Pose & pose)
{
	pose_ = pose.get_self_ptr();
	orbital_score_ = core::scoring::orbitals::OrbitalsScoreOP( new core::scoring::orbitals::OrbitalsScore() );
	pose.update_residue_neighbors();
}

RDFBaseOP RDFBinaryOrbitalFunctionCreator::create_rdf_function() const
{
	return RDFBaseOP( new RDFBinaryOrbitalFunction );
}

std::string RDFBinaryOrbitalFunctionCreator::type_name() const
{
	return "RDFBinaryOrbitalFunction";
}


RDFBinaryOrbitalFunction::RDFBinaryOrbitalFunction() : RDFBase("RDFBinaryOrbitalFunction"),pose_(/* NULL */)
{

	this->add_function_name("pi_pi_counts");
	this->add_function_name("salt_bridge_counts");
	this->add_function_name("hbond_orbital_counts");
	this->add_function_name("cation_pi_counts");

}

RDFBinaryOrbitalFunction::~RDFBinaryOrbitalFunction()
{

}

void RDFBinaryOrbitalFunction::parse_my_tag(utility::tag::TagCOP, basic::datacache::DataMap &) {}

RDFResultList RDFBinaryOrbitalFunction::operator()(AtomPairData const & atom_data )
{
	core::pose::PoseOP pose( pose_ );
	core::conformation::Residue const & protein_residue(pose->conformation().residue(atom_data.protein_atom_id.rsd()));
	core::conformation::Residue const & ligand_residue(pose->conformation().residue(atom_data.ligand_atom_id.rsd()));

	utility::vector1< core::Size > const & protein_orbitals(protein_residue.bonded_orbitals(atom_data.protein_atom_id.atomno()));
	utility::vector1< core::Size > const & ligand_orbitals(ligand_residue.bonded_orbitals(atom_data.ligand_atom_id.atomno()));
	core::Real pi_pi_counts = 0.0;
	core::Real salt_bridge_counts = 0.0;
	core::Real hbond_orbital_counts = 0.0;
	core::Real cation_pi_counts = 0.0;
	for(
		utility::vector1< core::Size >::const_iterator
		protein_orb_it = protein_orbitals.begin(),
		protein_orb_end = protein_orbitals.end();
		protein_orb_it != protein_orb_end; ++protein_orb_it)
	{
		std::string protein_orb = protein_residue.orbital_type(*protein_orb_it).name();
		for(
			utility::vector1< core::Size >::const_iterator
			ligand_orb_it = ligand_orbitals.begin(),
			ligand_orb_end = ligand_orbitals.end();
			ligand_orb_it != ligand_orb_end; ++ligand_orb_it)
		{
			std::string ligand_orb = ligand_residue.orbital_type(*ligand_orb_it).name();
			if(
			   (ligand_orb == "C.pi.sp2" && protein_orb == "C.pi.sp2") ||
			   (ligand_orb == "C.pi.sp2" && protein_orb == "N.pi.sp2") ||
			   (ligand_orb == "N.pi.sp2" && protein_orb == "C.pi.sp2"))
			{
				pi_pi_counts += 1.0;
			}
		}
	}

	std::string protein_atom_name = protein_residue.type().atom_type((atom_data.protein_atom_id.atomno())).name();
	std::string ligand_atom_name = ligand_residue.type().atom_type((atom_data.ligand_atom_id.atomno())).name();
	if(protein_atom_name == "Hpol" || protein_atom_name == "Haro")
	{
		for(
			utility::vector1< core::Size >::const_iterator
			ligand_orb_it = ligand_orbitals.begin(),
			ligand_orb_end = ligand_orbitals.end();
			ligand_orb_it != ligand_orb_end; ++ligand_orb_it)
		{
			std::string ligand_orb = ligand_residue.orbital_type(*ligand_orb_it).name();
			if (protein_atom_name == "Hpol")
			{
				if(ligand_orb == "O.p.sp2")
				{
					salt_bridge_counts += 1.0;
				}else if(ligand_orb == "O.p.sp3")
				{
					hbond_orbital_counts += 1.0;
				}else if(ligand_orb == "C.pi.sp2")
				{
					cation_pi_counts += 1.0;
				}
			}else if(protein_atom_name == "Haro")
			{
				if(ligand_orb == "C.pi.sp2")
				{
					pi_pi_counts += 1.0;
				}
			}
		}
	}else if(ligand_atom_name == "Hpol" || ligand_atom_name == "Haro")
	{
		for(
			utility::vector1< core::Size >::const_iterator
			protein_orb_it = protein_orbitals.begin(),
			protein_orb_end = protein_orbitals.end();
			protein_orb_it != protein_orb_end; ++protein_orb_it)
		{
			std::string protein_orb = protein_residue.orbital_type(*protein_orb_it).name();
			if (ligand_atom_name == "Hpol")
			{
				if(protein_orb == "O.p.sp2")
				{
					salt_bridge_counts += 1.0;
				}else if(protein_orb == "O.p.sp3")
				{
					hbond_orbital_counts += 1.0;
				}else if(protein_orb == "C.pi.sp2")
				{
					cation_pi_counts += 1.0;
				}
			}else if(ligand_atom_name == "Haro")
			{
				if(protein_orb == "C.pi.sp2")
				{
					pi_pi_counts += 1.0;
				}
			}
		}
	}

	RDFResultList results;
	results.push_back(std::make_pair("pi_pi_counts", pi_pi_counts));
	results.push_back(std::make_pair("salt_bridge_counts",salt_bridge_counts));
	results.push_back(std::make_pair("hbond_orbital_counts",hbond_orbital_counts));
	results.push_back(std::make_pair("cation_pi_counts",cation_pi_counts));
	return results;

}

void RDFBinaryOrbitalFunction::preamble(core::pose::Pose & pose)
{
	pose_ = pose.get_self_ptr();
}

}
}
}
