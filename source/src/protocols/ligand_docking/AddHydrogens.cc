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

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer add_hydrogens_tracer( "protocols.ligand_docking.LigandDesign", basic::t_debug );

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
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*datamap*/,
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
	utility::vector1<core::Size> unconnected_ids= find_unconnected_residues(pose, chain_begin, chain_end);

	BOOST_FOREACH(core::Size unconnected_id, unconnected_ids){
		core::conformation::Residue const & res_to_fix= pose.residue(unconnected_id);
		utility::vector1<core::Size> connect_ids= get_incomplete_connections(&res_to_fix);

		BOOST_FOREACH(core::Size connect_id, connect_ids){
			AddHydrogen AH(unconnected_id, connect_id);
			AH.apply(pose);
		}
	}
}

} // namespace ligand_docking
} // namespace protocols
