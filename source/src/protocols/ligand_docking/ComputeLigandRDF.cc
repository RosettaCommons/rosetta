// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ComputeLigandRDF.cc
///
/// @brief A mover for computing RDF functions
/// @author Sam DeLuca (sam@decarboxy.com)

#include <protocols/ligand_docking/ComputeLigandRDF.hh>
#include <protocols/ligand_docking/ComputeLigandRDFCreator.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/orbitals/OrbitalsScore.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>


#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/orbitals/OrbitalType.hh>

#include <core/id/AtomID.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <basic/Tracer.hh>

#include <map>
#include <string>

namespace protocols {
namespace ligand_docking {

static basic::Tracer compute_rdf_tracer("protocols.ligand_docking.ComputeLigandRDF");

std::string
ComputeLigandRDFCreator::keyname() const
{
	return ComputeLigandRDFCreator::mover_name();
}

protocols::moves::MoverOP
ComputeLigandRDFCreator::create_mover() const {
	return new ComputeLigandRDF;
}

std::string
ComputeLigandRDFCreator::mover_name()
{
	return "ComputeLigandRDF";
}


ComputeLigandRDF::ComputeLigandRDF() : mode_(""),ligand_chain_(""),bin_count_(100), smoothing_factor_(100.0)
{
	
}

ComputeLigandRDF::~ComputeLigandRDF()
{
	
}

ComputeLigandRDF::ComputeLigandRDF(ComputeLigandRDF const & that) :
	Mover(that),
	mode_(that.mode_),
	ligand_chain_(that.ligand_chain_),
	bin_count_(that.bin_count_),
	smoothing_factor_(that.smoothing_factor_),
	score_fxn_(that.score_fxn_)
{
	
}

protocols::moves::MoverOP ComputeLigandRDF::clone() const {
	return new ComputeLigandRDF( *this );
}

protocols::moves::MoverOP ComputeLigandRDF::fresh_instance() const {
	return new ComputeLigandRDF;
}

void ComputeLigandRDF::apply( core::pose::Pose & pose )
{

	std::map<std::string, utility::vector1<core::Real> > rdf_map;

	if(mode_ == "pocket")
	{
		compute_rdf_tracer << "computing ligand pocket RDFs" <<std::endl;
		rdf_map = protein_protein_rdf(pose);
	}else if(mode_ == "interface")
	{
		compute_rdf_tracer << "computing ligand_interface RDFs" <<std::endl;
		rdf_map = ligand_protein_rdf(pose);
	}

	compute_rdf_tracer << "Done computing RDF, appending data to job" <<std::endl;
	std::string rep_rdf_string(utility::join(rdf_map["fa_rep"]," "));
	std::string atr_rdf_string(utility::join(rdf_map["fa_atr"]," "));
	std::string sol_rdf_string(utility::join(rdf_map["fa_sol"]," "));
	std::string A_A1_rdf_string(utility::join(rdf_map["A_A1"]," "));
	std::string fa_elec_rdf_string(utility::join(rdf_map["fa_elec"]," "));
	std::string hbond_rdf_string(utility::join(rdf_map["hbond"]," "));
	std::string hbond_binary_rdf_string(utility::join(rdf_map["hbond_binary"]," "));
	std::string pci_cation_pi_rdf_string(utility::join(rdf_map["pci_cation_pi"]," " ));
	std::string pci_pi_pi_rdf_string(utility::join(rdf_map["pci_pi_pi"]," " ));
	std::string pci_hbond_rdf_string(utility::join(rdf_map["pci_hbond"]," " ));
	std::string orbitals_hpool_bb_rdf_string(utility::join(rdf_map["orbitals_hpol_bb"]," " ));
	std::string pci_salt_bridge_rdf_string(utility::join(rdf_map["pci_salt_bridge"]," " ));
	std::string pi_pi_binary_rdf_string(utility::join(rdf_map["pi_pi_binary"]," "));
	std::string salt_bridge_binary_rdf_string(utility::join(rdf_map["salt_bridge_binary"]," "));
	std::string hbond_orbital_binary_rdf_string(utility::join(rdf_map["hbond_orbital_binary"]," "));
	std::string cation_pi_binary_rdf_string(utility::join(rdf_map["cation_pi_binary"]," "));

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());

	job->add_string_string_pair("fa_rep_"+mode_+"_rdf",rep_rdf_string);
	job->add_string_string_pair("fa_atr_"+mode_+"_rdf",atr_rdf_string);
	job->add_string_string_pair("fa_sol_"+mode_+"_rdf",sol_rdf_string);
	job->add_string_string_pair("A_A1_"+mode_+"_rdf",A_A1_rdf_string);
	job->add_string_string_pair("fa_elec_"+mode_+"_rdf",fa_elec_rdf_string);
	job->add_string_string_pair("hbond_"+mode_+"_rdf",hbond_rdf_string);
	job->add_string_string_pair("hbond_binary_"+mode_+"_rdf",hbond_binary_rdf_string);
	job->add_string_string_pair("pci_cation_pi_"+mode_+"_rdf",pci_cation_pi_rdf_string);
	job->add_string_string_pair("pci_pi_pi_"+mode_+"_rdf",pci_pi_pi_rdf_string);
	job->add_string_string_pair("pci_hbond_"+mode_+"_rdf",pci_hbond_rdf_string);
	job->add_string_string_pair("orbitals_hpol_bb_"+mode_+"_rdf",orbitals_hpool_bb_rdf_string);
	job->add_string_string_pair("pci_salt_bridge_"+mode_+"_rdf",pci_salt_bridge_rdf_string);
	job->add_string_string_pair("pi_pi_binary_"+mode_+"_rdf",pi_pi_binary_rdf_string);
	job->add_string_string_pair("salt_bridge_binary_"+mode_+"_rdf",salt_bridge_binary_rdf_string);
	job->add_string_string_pair("hbond_orbital_binary_"+mode_+"_rdf",hbond_orbital_binary_rdf_string);
	job->add_string_string_pair("cation_pi_binary_"+mode_+"_rdf",cation_pi_binary_rdf_string);
}

std::string ComputeLigandRDF::get_name() const
{
	return "ComputeLigandRDF";
}
	
void ComputeLigandRDF::parse_my_tag
(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if(! tag->hasOption("mode"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("Mover 'ComputeLigandRDF' needs option 'mode'");
	}
	
	//right now, we assume that all things that aren't the ligand chain are the protein chain
	if(!tag->hasOption("ligand_chain"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("Mover 'ComputeLigandRDF' needs option 'ligand_chain'");
	}
	
	if ( ! tag->hasOption("scorefxn") ) throw utility::excn::EXCN_RosettaScriptsOption("'ComputeLigandRDF' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= data_map.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name);
	
	mode_ = tag->getOption<std::string>("mode");
	ligand_chain_ = tag->getOption<std::string>("ligand_chain");

	if(mode_ != "pocket" && mode_ != "interface")
	{
		throw utility::excn::EXCN_RosettaScriptsOption("Mover 'ComputeLigandRDF' must have a mode of either 'pocket' or 'interface'");
	}
}

//@brief do the work of computing an RDF.  look at Equation 18 of the ADRIANA.Code Manual (http://mol-net.com/files/docs/adrianacode/adrianacode_manual.pdf)
std::map<std::string, utility::vector1<core::Real> > ComputeLigandRDF::ligand_protein_rdf(core::pose::Pose & pose)
{
	core::scoring::TenANeighborGraph const & graph(pose.energies().tenA_neighbor_graph());
	int chain_id = core::pose::get_chain_id_from_chain(ligand_chain_, pose);
	core::Size ligand_start = pose.conformation().chain_begin(chain_id);
	core::Size ligand_end = pose.conformation().chain_end(chain_id);

	compute_rdf_tracer << "Identifying ligand interface atom pairs"<<std::endl;
	utility::vector1<std::pair<core::id::AtomID, core::id::AtomID> > atom_pair_list;
	//outer loop over ligand residue(s)
	for(core::Size ligand_resnum = ligand_start; ligand_resnum <= ligand_end; ++ligand_resnum)
	{
		core::conformation::Residue ligand_residue(pose.conformation().residue(ligand_resnum));

		for(core::Size ligand_atom_index = 1; ligand_atom_index <= ligand_residue.natoms();++ligand_atom_index)
		{
			core::id::AtomID ligand_atom_id(ligand_atom_index,ligand_resnum);
			//every residue within 12 A of ligand
			core::Vector ligand_atom_coords(ligand_residue.xyz(ligand_atom_index));
			for ( core::graph::Graph::EdgeListConstIter
				ir  = graph.get_node( ligand_resnum )->const_edge_list_begin(),
				ire = graph.get_node( ligand_resnum )->const_edge_list_end();
				ir != ire; ++ir )
			{
				core::Size protein_resnum = (*ir)->get_other_ind(ligand_resnum);
				if(pose.chain(protein_resnum) == chain_id)
				{
					continue;
				}
				core::conformation::Residue protein_residue(pose.conformation().residue((*ir)->get_other_ind(ligand_resnum)));
				for(core::Size protein_atom_index = 1; protein_atom_index <= protein_residue.natoms();++protein_atom_index)
				{
					core::id::AtomID protein_atom_id(protein_atom_index,protein_resnum);
					atom_pair_list.push_back(std::make_pair(ligand_atom_id,protein_atom_id));
				}
			}
		}
	}
	return compute_rdf(pose,atom_pair_list);
}

std::map<std::string, utility::vector1<core::Real> > ComputeLigandRDF::protein_protein_rdf(core::pose::Pose & pose)
{
	core::scoring::TenANeighborGraph const & graph(pose.energies().tenA_neighbor_graph());
	int chain_id = core::pose::get_chain_id_from_chain(ligand_chain_, pose);
	core::Size ligand_start = pose.conformation().chain_begin(chain_id);
	core::Size ligand_end = pose.conformation().chain_end(chain_id);

	utility::vector1<std::pair<core::id::AtomID, core::id::AtomID> > atom_pair_list;

	compute_rdf_tracer << "Identifying ligand pocket atom pairs"<<std::endl;
	utility::vector1<core::id::AtomID> binding_pocket_atoms;
	for(core::Size ligand_resnum = ligand_start; ligand_resnum <= ligand_end; ++ligand_resnum)
	{
		for ( core::graph::Graph::EdgeListConstIter
			ir  = graph.get_node( ligand_resnum )->const_edge_list_begin(),
			ire = graph.get_node( ligand_resnum )->const_edge_list_end();
			ir != ire; ++ir )
		{
			core::Size protein_resnum = (*ir)->get_other_ind(ligand_resnum);
			if(pose.chain(protein_resnum) == chain_id)
			{
				continue;
			}
			core::conformation::Residue protein_residue(pose.conformation().residue((*ir)->get_other_ind(ligand_resnum)));
			for(core::Size protein_atom_index = 1; protein_atom_index <= protein_residue.natoms();++protein_atom_index)
			{
				core::id::AtomID protein_atom_id(protein_atom_index,protein_resnum);
				binding_pocket_atoms.push_back(protein_atom_id);
			}
		}
	}

	for(core::Size outer_atom_index = 1; outer_atom_index <= binding_pocket_atoms.size();++outer_atom_index)
	{
		for(core::Size inner_atom_index = outer_atom_index; inner_atom_index <= binding_pocket_atoms.size();++inner_atom_index)
		{
			if(outer_atom_index == inner_atom_index)
			{
				continue;
			}
			atom_pair_list.push_back(std::make_pair(binding_pocket_atoms[outer_atom_index],binding_pocket_atoms[inner_atom_index]));
		}
	}
	return compute_rdf(pose,atom_pair_list);

}

std::map<std::string, utility::vector1<core::Real> > ComputeLigandRDF::compute_rdf(
	core::pose::Pose & pose,
	utility::vector1<std::pair<core::id::AtomID, core::id::AtomID> > const & atom_pairs )
{
	compute_rdf_tracer << "Computing " << bin_count_ << " bin RDF using " << atom_pairs.size() << " atom pairs" <<std::endl;
	core::Real step_size = 10.0/static_cast<core::Real>(bin_count_);

	std::map<std::string, utility::vector1<core::Real> > rdf_map;
	utility::vector1<core::Real> rdf_vector(bin_count_,0);
	rdf_map["fa_atr"] = rdf_vector;
	rdf_map["fa_rep"] = rdf_vector;
	rdf_map["fa_sol"] = rdf_vector;
	rdf_map["A_A1"] = rdf_vector;
	rdf_map["fa_elec"] = rdf_vector;
	rdf_map["hbond"] = rdf_vector;
	rdf_map["hbond_binary"] = rdf_vector;
	rdf_map["pci_cation_pi"] = rdf_vector;
	rdf_map["pci_pi_pi"] = rdf_vector;
	rdf_map["pci_hbond"] = rdf_vector;
	rdf_map["orbitals_hpol_bb"] = rdf_vector;
	rdf_map["pci_salt_bridge"] = rdf_vector;
	rdf_map["pi_pi_binary"] = rdf_vector;
	rdf_map["salt_bridge_binary"] = rdf_vector;
	rdf_map["hbond_orbital_binary"] = rdf_vector;
	rdf_map["cation_pi_binary"] = rdf_vector;

	core::scoring::methods::EnergyMethodOptions options(score_fxn_->energy_method_options());
	core::scoring::etable::EtableCAP etable(core::scoring::ScoringManager::get_instance()->etable( options.etable_type()));
	core::scoring::etable::AnalyticEtableEvaluator etable_evaluator(*etable);
	core::scoring::etable::coulomb::Coulomb coloumb(options);
	coloumb.initialize();
	
	core::scoring::hbonds::HBondSet hbond_set;
	pose.update_residue_neighbors();
	hbond_set.setup_for_residue_pair_energies(pose,false,false);
	
	
	for(core::Size pair_index = 1; pair_index <= atom_pairs.size();++pair_index)
	{
		core::id::AtomID protein_atom_id(atom_pairs[pair_index].first);
		core::id::AtomID ligand_atom_id(atom_pairs[pair_index].second);

		core::Vector protein_atom_coords(pose.conformation().xyz(protein_atom_id));
		core::Vector ligand_atom_coords(pose.conformation().xyz(ligand_atom_id));

		core::Real atom_atom_distance =  ligand_atom_coords.distance(protein_atom_coords);

		if(atom_atom_distance > 10.0)
		{
			continue;
		}

		std::map<std::string,core::Real> rdf_energies;

		core::conformation::Atom protein_atom(pose.conformation().residue(protein_atom_id.rsd()).atom(protein_atom_id.atomno()));
		core::conformation::Atom ligand_atom(pose.conformation().residue(ligand_atom_id.rsd()).atom(ligand_atom_id.atomno()));

		core::Real protein_atom_charge(pose.conformation().residue(protein_atom_id.rsd()).type().atom(protein_atom_id.atomno()).charge());
		core::Real ligand_atom_charge(pose.conformation().residue(ligand_atom_id.rsd()).type().atom(ligand_atom_id.atomno()).charge());

		core::chemical::AtomType protein_atom_type(pose.conformation().residue_type(protein_atom_id.rsd()).atom_type(protein_atom_id.atomno()));
		core::chemical::AtomType ligand_atom_type(pose.conformation().residue_type(ligand_atom_id.rsd()).atom_type(ligand_atom_id.atomno()));
		
		get_etable_energies(etable_evaluator,protein_atom,ligand_atom,rdf_energies);

		get_fa_elec_energy(coloumb,protein_atom_coords,protein_atom_charge,ligand_atom_coords,ligand_atom_charge,rdf_energies);
		
		get_hbond_energies(hbond_set,ligand_atom_id,protein_atom_id,rdf_energies);
		
		get_binary_hbond_energies(protein_atom_type,ligand_atom_type,rdf_energies);
		
		get_orbital_energies(pose,protein_atom_id,ligand_atom_id,rdf_energies);

		get_orbital_pair_counts(pose,protein_atom_id,ligand_atom_id,rdf_energies);

		//compute the RDF component for each bin
		for(core::Size bin_index = 1; bin_index <= rdf_vector.size(); bin_index++)
		{
			core::Real shell_radius = bin_index*step_size;
			core::Real rdf_bin_value = exp(-smoothing_factor_*(pow((shell_radius-atom_atom_distance),2)));

			
			rdf_map["fa_atr"][bin_index] += rdf_energies["atr_energy"]*rdf_bin_value;
			rdf_map["fa_rep"][bin_index] += rdf_energies["rep_energy"]*rdf_bin_value;
			rdf_map["fa_sol"][bin_index] += rdf_energies["solv_energy"]*rdf_bin_value;
			rdf_map["fa_elec"][bin_index] += rdf_energies["fa_elec_energy"]*rdf_bin_value;
			rdf_map["hbond"][bin_index]	 += rdf_energies["hbond_energy"]*rdf_bin_value;
			rdf_map["A_A1"][bin_index] += rdf_bin_value;
			rdf_map["hbond_binary"][bin_index] += rdf_energies["hbond_binary_energy"]*rdf_bin_value;
			rdf_map["pci_cation_pi"][bin_index] += rdf_energies["pci_cation_pi_energy"]*rdf_bin_value;
			rdf_map["pci_pi_pi"][bin_index] += rdf_energies["pci_pi_pi_energy"]*rdf_bin_value;
			rdf_map["pci_hbond"][bin_index] += rdf_energies["pci_hbond_energy"]*rdf_bin_value;
			rdf_map["orbitals_hpol_bb"][bin_index] += rdf_energies["orbitals_hpol_bb"]*rdf_bin_value;
			rdf_map["pci_salt_bridge"][bin_index] += rdf_energies["pci_salt_bridge"]*rdf_bin_value;
			rdf_map["pi_pi_binary"][bin_index] += rdf_energies["pi_pi_counts"]*rdf_bin_value;
			rdf_map["salt_bridge_binary"][bin_index] += rdf_energies["salt_bridge_counts"]*rdf_bin_value;
			rdf_map["hbond_orbital_binary"][bin_index] += rdf_energies["hbond_orbital_counts"]*rdf_bin_value;
			rdf_map["cation_pi_binary"][bin_index] += rdf_energies["cation_pi_counts"]*rdf_bin_value;
		}
	}
	//Normalize the rdf bin
	std::map<std::string, utility::vector1<core::Real> >::iterator rdf_map_it;
	for(rdf_map_it = rdf_map.begin();rdf_map_it != rdf_map.end();++rdf_map_it)
	{
		for(core::Size bin_index = 1; bin_index <= rdf_map_it->second.size();++bin_index)
		{
			rdf_map_it->second[bin_index] = rdf_map_it->second[bin_index]/2;
		}
	}

	return rdf_map;
}

void ComputeLigandRDF::get_etable_energies(
		core::scoring::etable::AnalyticEtableEvaluator const & etable_evaluator,
		core::conformation::Atom const & protein_atom,
		core::conformation::Atom const & ligand_atom,
		std::map<std::string,core::Real> & rdf_energies)
{
	core::Real weight = 1.0;
	core::Real atr_energy = 0.0;
	core::Real rep_energy = 0.0;
	core::Real solv_energy = 0.0;
	core::Real d2 = 0.0;
	etable_evaluator.atom_pair_energy(protein_atom,ligand_atom,weight,atr_energy,rep_energy,solv_energy,d2);

	rdf_energies["atr_energy"] = atr_energy;
	rdf_energies["rep_energy"] = rep_energy;
	rdf_energies["solv_energy"] = solv_energy;
}

void ComputeLigandRDF::get_fa_elec_energy(
		core::scoring::etable::coulomb::Coulomb const & coloumb,
		core::Vector const & protein_atom_coords,
		core::Real const & protein_atom_charge,
		core::Vector const & ligand_atom_coords,
		core::Real const & ligand_atom_charge,
		std::map<std::string,core::Real> & rdf_energies)
{
	rdf_energies["fa_elec_energy"] = coloumb.eval_atom_atom_fa_elecE(protein_atom_coords,protein_atom_charge,ligand_atom_coords,ligand_atom_charge);

}

void ComputeLigandRDF::get_hbond_energies(
	core::scoring::hbonds::HBondSet const & hbond_set,
	core::id::AtomID const & ligand_atom_id,
	core::id::AtomID const & protein_atom_id,
	std::map<std::string,core::Real> & rdf_energies)
{
	//get hbond energy if availible
	utility::vector1< core::scoring::hbonds::HBondCOP > hbond_list(hbond_set.atom_hbonds(protein_atom_id));
	core::Real hbond_energy = 0.0;
	//Most atoms will have zero or some small number of hbonds, this is a short loop
	for(core::Size hbond_index = 1; hbond_index <= hbond_list.size();++hbond_index)
	{
		core::scoring::hbonds::HBondCOP current_bond(hbond_list[hbond_index]);
		core::id::AtomID acc_id(current_bond->acc_atm(),current_bond->acc_res());
		core::id::AtomID don_id(current_bond->don_hatm(),current_bond->don_res());
		if((acc_id == ligand_atom_id) || (don_id == ligand_atom_id))
		{
			rdf_energies["hbond_energy"] = current_bond->energy();
			return;
		}
	}
	rdf_energies["hbond_energy"] = 0.0;

}

void ComputeLigandRDF::get_binary_hbond_energies(
	core::chemical::AtomType const & protein_atom_type,
	core::chemical::AtomType const & ligand_atom_type,
	std::map<std::string,core::Real> & rdf_energies)
{
	//is the current pair of atoms a donor and acceptor?
	if(protein_atom_type.is_donor() && ligand_atom_type.is_acceptor())
	{
		rdf_energies["hbond_binary_energy"] = 1.0;
	}else if(protein_atom_type.is_acceptor() && ligand_atom_type.is_donor())
	{
		rdf_energies["hbond_binary_energy"] = 1.0;
	}else{
		rdf_energies["hbond_binary_energy"] = 0.0;
	}


}

void ComputeLigandRDF::get_orbital_energies(
	core::pose::Pose & pose,
	core::id::AtomID const & protein_atom_id,
	core::id::AtomID const & ligand_atom_id,
	std::map<std::string,core::Real> & rdf_energies)
{
	core::scoring::EnergyMap emap;
	emap.clear();
	core::scoring::orbitals::OrbitalsScore orbital_score;

	pose.update_residue_neighbors();
	orbital_score.get_E_haro_one_way(pose,protein_atom_id,ligand_atom_id,emap);
	orbital_score.get_E_haro_one_way(pose,ligand_atom_id,protein_atom_id,emap);

	orbital_score.get_E_hpol_one_way(pose,protein_atom_id,ligand_atom_id,emap);
	orbital_score.get_E_hpol_one_way(pose,ligand_atom_id,protein_atom_id,emap);

	orbital_score.get_orb_orb_E(pose,protein_atom_id,ligand_atom_id,emap);
	orbital_score.get_orb_orb_E(pose,ligand_atom_id,protein_atom_id,emap);

	rdf_energies["pci_cation_pi_energy"] = emap[core::scoring::pci_cation_pi];
	rdf_energies["pci_pi_pi_energy"] = emap[core::scoring::pci_pi_pi];
	rdf_energies["pci_hbond_energy"] = emap[core::scoring::pci_hbond];
	rdf_energies["orbitals_hpol_bb"] = emap[core::scoring::orbitals_hpol_bb];
	rdf_energies["pci_salt_bridge"] = emap[core::scoring::pci_salt_bridge];


}

//this isn't the most efficient thing in the world right now but im not sure it matters
void ComputeLigandRDF::get_orbital_pair_counts(
	core::pose::Pose & pose,
	core::id::AtomID const & protein_atom_id,
	core::id::AtomID const & ligand_atom_id,
	std::map<std::string,core::Real> & rdf_energies
)
{
	core::conformation::Residue const & protein_residue(pose.conformation().residue(protein_atom_id.rsd()));
	core::conformation::Residue const & ligand_residue(pose.conformation().residue(ligand_atom_id.rsd()));

	utility::vector1< core::Size > const & protein_orbitals(protein_residue.bonded_orbitals(protein_atom_id.atomno()));
	utility::vector1< core::Size > const & ligand_orbitals(ligand_residue.bonded_orbitals(ligand_atom_id.atomno()));
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

	std::string protein_atom_name = protein_residue.type().atom_type((protein_atom_id.atomno())).name();
	std::string ligand_atom_name = ligand_residue.type().atom_type((ligand_atom_id.atomno())).name();
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


	rdf_energies["pi_pi_counts"] = pi_pi_counts;
	rdf_energies["salt_bridge_counts"] = salt_bridge_counts;
	rdf_energies["hbond_orbital_counts"] = hbond_orbital_counts;
	rdf_energies["cation_pi_counts"] = cation_pi_counts;
}

}
}
