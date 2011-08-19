// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/AddHydrogens.cc
///
/// @brief
/// @Gordon Lemmon

#include <protocols/ligand_docking/AddHydrogens.hh>
#include <protocols/ligand_docking/AddHydrogensCreator.hh>
#include <protocols/ligand_docking/AddHydrogen.hh>
#include <protocols/ligand_docking/LigandDesign.hh>

#include <core/pose/util.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// option key includes
#include <utility/tag/Tag.hh>

//Auto Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
//#include <numeric/random/random.hh>

namespace protocols {
namespace ligand_docking {

static basic::Tracer add_hydrogens_tracer("protocols.ligand_docking.LigandDesign", basic::t_debug);

std::string
AddHydrogensCreator::keyname() const
{
	return AddHydrogensCreator::mover_name();
}

protocols::moves::MoverOP
AddHydrogensCreator::create_mover() const {
	return new AddHydrogens;
}

std::string
AddHydrogensCreator::mover_name()
{
	return "AddHydrogens";
}

AddHydrogens::AddHydrogens():
		//utility::pointer::ReferenceCount(),
		Mover("AddHydrogens"),
		chain_("X")
{
	Mover::type( "AddHydrogens" );
}

AddHydrogens::AddHydrogens(AddHydrogens const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_)
{}

AddHydrogens::~AddHydrogens() {}

protocols::moves::MoverOP AddHydrogens::clone() const {
	return new AddHydrogens( *this );
}

protocols::moves::MoverOP AddHydrogens::fresh_instance() const {
	return new AddHydrogens;
}

std::string AddHydrogens::get_name() const{
	return "AddHydrogens";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
AddHydrogens::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	assert( tag->getName() != "AddHydrogens");
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'AddHydrogens' requires 'chain' tag");
	chain_ = tag->getOption<std::string>("chain");
}

void
AddHydrogens::apply( core::pose::Pose & pose )
{
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const chain_begin= pose.conformation().chain_begin(chain_id);
	core::Size const chain_end= pose.conformation().chain_end(chain_id);
	utility::vector1<core::Size> unconnected= find_unconnected_residues(pose, chain_begin, chain_end);

	utility::vector1<core::Size>::iterator unconnected_begin= unconnected.begin();
	utility::vector1<core::Size>::iterator const unconnected_end= unconnected.end();
	for(; unconnected_begin != unconnected_end ; ++unconnected_begin){
		core::conformation::Residue const & res_to_fix= pose.residue(*unconnected_begin);

		utility::vector1<core::Size> connect_ids= get_incomplete_connections(&res_to_fix);

		utility::vector1<core::Size>::iterator connect_ids_begin= connect_ids.begin();
		utility::vector1<core::Size>::iterator const connect_ids_end= connect_ids.end();
		for(; connect_ids_begin != connect_ids_end; ++connect_ids_begin){
			core::Size connect_id= *connect_ids_begin;

			AddHydrogen AH(*unconnected_begin, connect_id);
			AH.apply(pose);
		}
	}
}

} // namespace ligand_docking
} // namespace protocols
