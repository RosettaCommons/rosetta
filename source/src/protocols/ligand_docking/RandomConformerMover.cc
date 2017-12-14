// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/RandomConformerMover.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/RandomConformerMover.hh>

#include <utility/graph/Graph.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <numeric/random/random.hh>

// option key includes

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {


RandomConformerMover::RandomConformerMover(core::Size resid):
	Mover(),
	resid_(resid)
{
	Mover::type( "RandomConformerMover" );
}


RandomConformerMover::~RandomConformerMover() = default;


void RandomConformerMover::apply( core::pose::Pose & pose )
{
	using core::conformation::ResidueOP;
	using namespace core::pack::task;
	utility::vector1< ResidueOP > conformers;
	// Dummy parameters that the ligand rotamer library doesn't use:
	core::scoring::ScoreFunction dummy_scorefxn;
	PackerTaskOP dummy_pack_task = TaskFactory::create_packer_task(pose);
	dummy_pack_task->initialize_from_command_line(); // -ex1 -ex2  etc.
	utility::vector1< utility::vector1< core::Real > > dummy_extra_chi_steps;
	utility::graph::GraphCOP dummy_graph( utility::graph::GraphOP( new utility::graph::Graph() ) );

	// Retrieve conformers
	core::pack::rotamers::SingleResidueRotamerLibraryCOP reslib = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue_type(resid_) );
	if ( ! reslib ) return;
	reslib->fill_rotamer_vector(
		pose,
		dummy_scorefxn,
		*dummy_pack_task,
		dummy_graph,
		pose.residue_type_ptr(resid_), //ResidueTypeCOP
		pose.residue(resid_),
		dummy_extra_chi_steps,
		true /* sure, let's pretend it's buried */,
		conformers // output appended here
	);
	// If -include_current push back current conformer
	// This is rarely necessary...
	//if( option[ OptionKeys::packing::use_input_sc ] ) {
	// ResidueOP curr_copy = new Residue( pose.residue(resid_) );
	// conformers.push_back(curr_copy);
	//}
	// Choose one at random
	ResidueOP selected_res = conformers[ numeric::random::rg().random_range(1, conformers.size()) ];
	// Residue library has already superimpose residues appropriately, so don't orient again
	pose.replace_residue(resid_, *selected_res, false /*orient backbone*/);
}

std::string
RandomConformerMover::get_name() const {
	return "RandomConformerMover";
}


} // namespace ligand_docking
} // namespace protocols
