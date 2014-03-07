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
#include <protocols/ligand_docking/rdf/RDFFunctionFactory.hh>
#include <protocols/ligand_docking/rdf/RDFBase.hh>

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


ComputeLigandRDF::ComputeLigandRDF() : mode_(""),ligand_chain_(""),bin_count_(100), smoothing_factor_(100.0),range_squared_(100.0)
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
	range_squared_(that.range_squared_),
	rdf_functions_(that.rdf_functions_)
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
	
	std::map<std::string, utility::vector1<core::Real> >::iterator it;
	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	for(it = rdf_map.begin(); it != rdf_map.end();++it)
	{
		std::string rdf_string(utility::join(it->second," "));
		std::string rdf_tag(it->first+"_"+mode_+"_rdf");
		job->add_string_string_pair(rdf_tag, rdf_string);
	}
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
	
	mode_ = tag->getOption<std::string>("mode");
	ligand_chain_ = tag->getOption<std::string>("ligand_chain");

	if(mode_ != "pocket" && mode_ != "interface")
	{
		throw utility::excn::EXCN_RosettaScriptsOption("Mover 'ComputeLigandRDF' must have a mode of either 'pocket' or 'interface'");
	}
	
	bin_count_ = tag->getOption<core::Size>("bin_count",100);
	smoothing_factor_ = tag->getOption<core::Real>("smoothing_factor",100);
	range_squared_ = pow(tag->getOption<core::Real>("range",10.0),2.0);
	
	utility::vector0< TagCOP >::const_iterator begin=tag->getTags().begin();
	utility::vector0< TagCOP >::const_iterator end=tag->getTags().end();
	
	for(; begin != end; ++begin){
		TagCOP function_tag= *begin;
		//	BOOST_FOREACH(TagCOP const & feature_tag, tag->getTags()){
		
		if(function_tag->getName() != "RDF"){
			compute_rdf_tracer.Error << "Please include only tags with name 'RDF' as subtags of ReportToDB" << std::endl;
			compute_rdf_tracer.Error << "Tag with name '" << function_tag->getName() << "' is invalid" << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("");
		}
		
		rdf::RDFBaseOP rdf_function(rdf::RDFFunctionFactory::get_instance()->get_rdf_function(function_tag, data_map));
		
		rdf_functions_.push_back(rdf_function);

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
					core::Vector protein_atom_coords(protein_residue.xyz(protein_atom_index));
					if(protein_atom_coords.distance_squared(ligand_atom_coords) <= range_squared_)
					{
						core::id::AtomID protein_atom_id(protein_atom_index,protein_resnum);
						atom_pair_list.push_back(std::make_pair(ligand_atom_id,protein_atom_id));
					}
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
		core::Vector outer_atom_coords(pose.xyz(binding_pocket_atoms[outer_atom_index]));
		for(core::Size inner_atom_index = outer_atom_index; inner_atom_index <= binding_pocket_atoms.size();++inner_atom_index)
		{
			if(outer_atom_index == inner_atom_index)
			{
				continue;
			}
			core::Vector inner_atom_coords(pose.xyz(binding_pocket_atoms[inner_atom_index]));
			
			if(outer_atom_coords.distance_squared(inner_atom_coords) <= range_squared_)
			{
				atom_pair_list.push_back(std::make_pair(binding_pocket_atoms[outer_atom_index],binding_pocket_atoms[inner_atom_index]));
			}
		}
	}
	return compute_rdf(pose,atom_pair_list);

}

std::map<std::string, utility::vector1<core::Real> > ComputeLigandRDF::compute_rdf(
	core::pose::Pose & pose,
	utility::vector1<std::pair<core::id::AtomID, core::id::AtomID> > const & atom_pairs )
{
	compute_rdf_tracer << "Computing " << bin_count_ << " bin RDF using " << atom_pairs.size() << " atom pairs" <<std::endl;
	core::Real step_size = sqrt(range_squared_)/static_cast<core::Real>(bin_count_);

	std::map<std::string, utility::vector1<core::Real> > rdf_map;
	utility::vector1<core::Real> rdf_vector(bin_count_,0);
	
	//each RDF function object stores a list of the terms it generates
	//we fill rdf map with an vector of 0.0 RDF values for each generated term
	//the empty vectors will be populated further down.
	
	for(core::Size rdf_count = 1; rdf_count <= rdf_functions_.size();++rdf_count)
	{
		//While we're looping through all our RDFs we can also run the "preamble" function for each
		//Which will set up various data that only need to be computed once per pose
		rdf_functions_[rdf_count]->preamble(pose);
		utility::vector1<std::string> function_names(rdf_functions_[rdf_count]->get_function_names());
		
		for(core::Size func_count = 1; func_count <= function_names.size();++func_count)
		{
			rdf_map[function_names[func_count]] = rdf_vector;
		}
	}
	
	pose.update_residue_neighbors();
	
	for(core::Size pair_index = 1; pair_index <= atom_pairs.size();++pair_index)
	{
		// The AtomPairData object is a struct which is used to pass around derived
		// data that is needed by the various RDF functions.
		
		core::id::AtomID protein_atom_id(atom_pairs[pair_index].first);
		core::id::AtomID ligand_atom_id(atom_pairs[pair_index].second);
		
		core::chemical::AtomType protein_atom_type =pose.conformation().residue_type(protein_atom_id.rsd()).atom_type(protein_atom_id.atomno());
		core::chemical::AtomType ligand_atom_type =pose.conformation().residue_type(ligand_atom_id.rsd()).atom_type(ligand_atom_id.atomno());
		
		rdf::AtomPairData atom_data(protein_atom_type,ligand_atom_type);
		atom_data.protein_atom_id = protein_atom_id;
		atom_data.ligand_atom_id = ligand_atom_id;

		atom_data.protein_atom_coords = pose.conformation().xyz(atom_data.protein_atom_id);
		atom_data.ligand_atom_coords = pose.conformation().xyz(atom_data.ligand_atom_id);

		atom_data.atom_atom_distance =  atom_data.ligand_atom_coords.distance(atom_data.protein_atom_coords);

		// We use a 10 A cutoff for everything
		if(atom_data.atom_atom_distance > 10.0)
		{
			continue;
		}

		std::map<std::string,core::Real> rdf_energies;

		atom_data.protein_atom = pose.conformation().residue(atom_data.protein_atom_id.rsd()).atom(atom_data.protein_atom_id.atomno());
		
		atom_data.ligand_atom = pose.conformation().residue(atom_data.ligand_atom_id.rsd()).atom(atom_data.ligand_atom_id.atomno());

		atom_data.protein_atom_charge = pose.conformation().residue(atom_data.protein_atom_id.rsd()).type().atom(atom_data.protein_atom_id.atomno()).charge();
		atom_data.ligand_atom_charge = pose.conformation().residue(atom_data.ligand_atom_id.rsd()).type().atom(atom_data.ligand_atom_id.atomno()).charge();

		
		// Compute all the RDF data for this pair of atoms
		rdf::RDFResultList total_results;
		for(core::Size rdf_count = 1; rdf_count <= rdf_functions_.size();++rdf_count)
		{
			rdf::RDFBaseOP current_rdf = rdf_functions_[rdf_count];
			rdf::RDFResultList current_results(current_rdf->operator()(atom_data));
			total_results.insert(total_results.end(), current_results.begin(), current_results.end());
			
		}
		
		
		//compute the RDF component for each bin
		for(core::Size bin_index = 1; bin_index <= rdf_vector.size(); bin_index++)
		{
			core::Real shell_radius = bin_index*step_size;
			core::Real rdf_bin_value = exp(-smoothing_factor_*(pow((shell_radius-atom_data.atom_atom_distance),2)));

			for(rdf::RDFResultList::iterator it = total_results.begin(); it != total_results.end();++it)
			{
				std::string term_name(it->first);
				core::Real term_value(it->second);
				rdf_map[term_name][bin_index] += term_value*rdf_bin_value;
			}

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

}
}
