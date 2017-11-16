// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/qsar/qsarMover.cc
/// @author Sam DeLuca

#include <protocols/qsar/qsarMover.hh>
#include <protocols/qsar/qsarMoverCreator.hh>
#include <protocols/qsar/qsarMap.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/jd2/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace qsar {
static basic::Tracer TR( "protocols.ligand_docking.qsar.qsarMover" );


std::string
qsarCreator::keyname() const
{
	return qsarCreator::mover_name();
}

protocols::moves::MoverOP
qsarCreator::create_mover() const {
	return protocols::moves::MoverOP( new qsarMover );
}

std::string
qsarCreator::mover_name()
{
	return "qsar";
}

qsarMover::qsarMover():
	qsar_map_(/* 0 */),
	chain_()
{}

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
qsarMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	grid_set_prototype_ = scoring_grid::parse_grid_set_from_tag( tag, datamap );

	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'qsar' mover requires chain tag");
	chain_= tag->getOption<std::string>("chain");

	if ( ! tag->hasOption("grids") ) throw utility::excn::EXCN_RosettaScriptsOption("'qsar' mover requires grids tag");

	std::string grids_string= tag->getOption<std::string>("grids");
	utility::vector1<std::string> grid_strings(utility::string_split(grids_string, ','));
	grids_to_use_ = grid_strings;
}

void qsarMover::apply(core::pose::Pose & pose)
{
	debug_assert( grid_set_prototype_ != nullptr );

	if ( grids_to_use_.size()==0 ) {
		TR.Warning << "no grids specified, QSAR scoring function will be empty!!" <<std::endl;
		return;
	}

	debug_assert( chain_.size() == 1 );
	utility::vector1< core::Size > chain_residues( core::pose::get_resnums_for_chain( pose, chain_[0] ) );
	if ( chain_residues.size() != 1 ) {
		utility_exit_with_message("The qsarMover can only operate on single residue chains.");
	}

	core::Size resnum( chain_residues[1] );
	// need to make a copy, as qsarMap takes an OP to it.
	core::conformation::ResidueOP residue( new core::conformation::Residue(pose.residue(resnum)) );

	if ( qsar_map_ == nullptr ) {

		qsar_map_ = qsarMapOP( new qsarMap("default", residue) );

		qsar_map_->fill_with_value(1,grids_to_use_);

		scoring_grid::GridSetOP mod_prototype( grid_set_prototype_->clone() );
		mod_prototype->set_qsar_map(qsar_map_);
		grid_set_prototype_ = mod_prototype;
	}

	core::Vector center( core::pose::all_atom_center(pose, chain_residues) );

	scoring_grid::GridSetCOP grid_set = scoring_grid::GridManager::get_instance()->get_grids( *grid_set_prototype_, pose, center, chain_ );

	scoring_grid::GridSet::ScoreMap grid_scores( grid_set->grid_scores(*residue) );
	core::Real total_score(grid_set->total_score(*residue));
	TR.Debug << "total score is " << total_score <<std::endl;

	for ( scoring_grid::GridSet::ScoreMap::value_type const & pair: grid_scores ) {
		protocols::jd2::add_string_real_pair_to_current_job( "grid_"+pair.first, pair.second );
	}
	protocols::jd2::add_string_real_pair_to_current_job( "grid_total", total_score );

	//grid_manager->write_grids("test_");

}

std::string qsarMover::get_name() const
{
	return "qsarMover";
}

}
}
