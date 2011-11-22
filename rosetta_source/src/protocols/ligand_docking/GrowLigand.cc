// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/GrowLigand.cc
///
/// @brief
/// @Gordon Lemmon

#include <protocols/ligand_docking/GrowLigand.hh>
#include <protocols/ligand_docking/GrowLigandCreator.hh>

#include <protocols/ligand_docking/LigandDesign.hh> // For helper functions.  Refactor this
#include <core/pose/util.hh>


#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/docking.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <utility/tag/Tag.hh>

//Auto Headers
#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/random/random_permutation.hh>
#include <boost/foreach.hpp>

//Auto Headers
#include <core/conformation/Conformation.hh>
#define foreach BOOST_FOREACH

namespace protocols {
namespace ligand_docking {

static basic::Tracer grow_ligand_tracer("protocols.ligand_docking.GrowLigand", basic::t_debug);

std::string
GrowLigandCreator::keyname() const
{
	return GrowLigandCreator::mover_name();
}

protocols::moves::MoverOP
GrowLigandCreator::create_mover() const {
	return new GrowLigand;
}

std::string
GrowLigandCreator::mover_name()
{
	return "GrowLigand";
}

GrowLigand::GrowLigand():
		Mover("GrowLigand"),
		chain_("")
{
	set_fragments();
}

GrowLigand::GrowLigand(std::string chain):
		Mover("GrowLigand"),
		chain_(chain)
{
	set_fragments();
}

GrowLigand::GrowLigand(GrowLigand const & that):
	    //utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_),
		fragments_(that.fragments_)
{}

GrowLigand::~GrowLigand() {}

void
GrowLigand::set_fragments(){
	core::chemical::ResidueSelector rs;
	rs.set_property("FRAGMENT");
	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::ResidueTypeSetCAP rsd_set= cm->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical::AtomTypeSetCAP atom_set= cm->atom_type_set( core::chemical::FA_STANDARD );
	core::Size atom_set_size= atom_set->n_atomtypes();
	fragments_.assign(atom_set_size, utility::vector1< std::pair<core::conformation::ResidueCOP, core::Size> >() );
	core::chemical::ResidueTypeCAPs fragment_types= rs.select( *rsd_set );
	grow_ligand_tracer<< fragment_types.size()<< " fragment_types"<< std::endl;

	foreach(core::chemical::ResidueTypeCAP fragment_type, fragment_types){
		core::conformation::Residue* temp= new core::conformation::Residue(*fragment_type, true);
		for(core::Size i=1; i<= temp->n_residue_connections(); ++i){
			core::chemical::ResidueConnection const & res_conn= temp->residue_connection(i);
			int atom_index_number= res_conn.atomno();
			int atom_type_index= temp->atom_type_index(atom_index_number);
			std::pair<core::conformation::ResidueCOP, core::Size> pair(temp, i);
			fragments_[atom_type_index].push_back(pair);
		}

		grow_ligand_tracer<< "frag_name: "<< temp->name()<< std::endl;
	}
}

protocols::moves::MoverOP GrowLigand::clone() const {
	return new GrowLigand( *this );
}

protocols::moves::MoverOP GrowLigand::fresh_instance() const {
	return new GrowLigand;
}

std::string GrowLigand::get_name() const{
	return "GrowLigand";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
GrowLigand::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "GrowLigand" ) {
		grow_ligand_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( tag->hasOption("chain") ) {
		chain_ = tag->getOption<std::string>("chain");
	}else{
		utility_exit_with_message("HeavyAtom filter needs a 'chain' option");
	}
}

void
GrowLigand::apply( core::pose::Pose & pose )
{
	assert(!fragments_.empty());
	assert(chain_.size() == 1);

	utility::vector1<core::Size> unconnected_residues;
	{
		core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);;
		core::Size const start = pose.conformation().chain_begin(chain_id);
		core::Size const end = pose.conformation().chain_end(chain_id);
		unconnected_residues=find_unconnected_residues(pose, start, end);
	}

	numeric::random::random_permutation(unconnected_residues, numeric::random::RG);// shuffle vector

	foreach(core::Size grow_from, unconnected_residues){
		core::Size grow_from_connection= random_connection(&pose.residue(grow_from));
		core::Size atom_type_index= pose.residue(grow_from).residue_connection(grow_from_connection).atom_type_index();
		utility::vector1< std::pair<core::conformation::ResidueCOP, core::Size > >
				residue_connection_pairs= fragments_[atom_type_index];
		if( residue_connection_pairs.size() > 0 ){
			std::pair<core::conformation::ResidueCOP, core::Size > const growth= numeric::random::random_element(residue_connection_pairs);
			bool build_ideal_geometry= true;
			pose.append_residue_by_bond(*(growth.first), build_ideal_geometry, growth.second, grow_from, grow_from_connection);
			return;
		}
	}

}

void GrowLigand::fragments_to_string() const{
	for(core::Size i=1; i <= fragments_.size(); ++i){
		utility::vector1< std::pair<core::conformation::ResidueCOP, core::Size> >::const_iterator  begin= fragments_[i].begin();
		for(; begin != fragments_[i].end(); ++begin){
			core::conformation::Residue const & res= *(begin->first);
			core::Size connect_id= begin->second;
			std::string name= res.name();
			grow_ligand_tracer<< "atom_type, res_name, connection"<< i << " "<< name << " "<< connect_id<< std::endl;
		}
	}
}

void GrowLigand::add_scores_to_job(
	core::pose::Pose & /*pose*/
)
{
}

} // namespace ligand_docking
} // namespace protocols
