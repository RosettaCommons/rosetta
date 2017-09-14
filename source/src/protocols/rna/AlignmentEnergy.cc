// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/AlignmentEnergy.cc
/// @brief An RMSD-type energy to a reference pose, complete with derivatives hacked out of coordinate constraints
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/rna/AlignmentEnergy.hh>
#include <protocols/rna/AlignmentEnergyCreator.hh>

// Unit headers
#include <protocols/rna/AlignmentEnergy.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

// Project Headers
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh> // should NOT be in here!
#include <numeric/xyzVector.hh>

#include <core/scoring/EnergyMap.hh>

// Basic/Utility headers
#include <utility/vector1.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.AlignmentEnergy" );

namespace protocols {
namespace rna {

using namespace core;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::scoring::func;

/*AlignmentEnergy::AlignmentEnergy():
core::scoring::methods::WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new AlignmentEnergyCreator ) ),
align_pose_( nullptr ),
pose_aligner_( nullptr )
{}*/

AlignmentEnergy::AlignmentEnergy(
	core::scoring::methods::EnergyMethodOptions const & //options //core::pose::PoseOP const & align_pose,
	//utility::vector1< Size > const & moving_res_list
):
	core::scoring::methods::WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new AlignmentEnergyCreator ) )//,
	//align_pose_( align_pose ),
	//moving_res_list_( moving_res_list ),
{
	// Do the worst possible thing:
	using namespace basic::options;
	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	// Keep at nullptr and hope it gets initialized later if there's no align_pdb
	if (  option[ OptionKeys::stepwise::new_align_pdb ].user() ) {
		align_pose_ = core::import_pose::get_pdb_with_full_model_info(  option[ OptionKeys::stepwise::new_align_pdb ](), rsd_set );
		pose_aligner_ = stepwise::modeler::align::StepWisePoseAlignerOP( new stepwise::modeler::align::StepWisePoseAligner( *align_pose_ ) );
	}

	if ( option[ OptionKeys::stepwise::rmsd_screen ].user() ) {
		rmsd_screen_ = option[ OptionKeys::stepwise::rmsd_screen ].value();
		func_ = core::scoring::func::FuncOP( new core::scoring::func::FlatHarmonicFunc( 0, 1, rmsd_screen_ ) );
	}
}

void AlignmentEnergy::rmsd_screen( core::Real const setting ) {
	rmsd_screen_ = setting;
	// Update func too...
	func_ = core::scoring::func::FuncOP( new core::scoring::func::FlatHarmonicFunc( 0, 1, rmsd_screen_ ) );
}

// copy constructor (not needed unless you need deep copies)
//AlignmentEnergy::AlignmentEnergy( AlignmentEnergy const & src );

// destructor (important for properly forward-declaring smart-pointer members)
//AlignmentEnergy::~AlignmentEnergy();

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP
AlignmentEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new AlignmentEnergy( *this ) );
}

/// @brief Indicate required setup steps for scoring
void
AlignmentEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	pose.update_residue_neighbors();
	//moving_res_list_ = core::pose::full_model_info::get_moving_res_from_full_model_info( pose );
}

/// @brief Is the score context dependent or context independent?
void
AlignmentEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const {

}

void AlignmentEnergy::align_pose( core::pose::PoseOP const & align_pose ) {
	align_pose_ = align_pose;
	pose_aligner_ = stepwise::modeler::align::StepWisePoseAlignerOP( new stepwise::modeler::align::StepWisePoseAligner( *align_pose_ ) );
}


void
AlignmentEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const & ,//domain_map,
	ScoreFunction const & ,//sfxn,
	EnergyMap const & ,//emap,
	Vector & F1,
	Vector & F2
) const {
	Vector our_f1( 0 );
	Vector our_f2( 0 );

	utility::vector1< Size > root_partition_res = stepwise::modeler::figure_out_root_partition_res( pose, core::pose::full_model_info::get_moving_res_from_full_model_info_const( pose ) );
	pose_aligner_->set_root_partition_res( root_partition_res );
	pose_aligner_->initialize( pose );
	if ( pose_aligner_->superimpose_atom_id_map().size() == 0 ) return;

	Real const rmsd = core::scoring::rms_at_corresponding_atoms_no_super( pose, *align_pose_, pose_aligner_->superimpose_atom_id_map() );
	if ( rmsd == 0 ) return;
	Size const n = pose_aligner_->superimpose_atom_id_map().size();
	Real const dfunc = func_->dfunc( rmsd );
	if ( dfunc == 0 ) return;
	Real const factor =  dfunc / ( 2 * n * rmsd );
	if ( factor == 0 ) return;

	auto func = FuncOP( new HarmonicFunc( 0, 1 ) );
	for ( auto const & elem : pose_aligner_->superimpose_atom_id_map() ) {

		if ( id.rsd() != elem.first.rsd() ) continue;
		if ( id.atomno() != elem.first.atomno() ) continue;

		Vector const & xyz1( pose.residue( elem.first.rsd() ).xyz( elem.first.atomno() ) ), xyz2( align_pose_->residue( elem.second.rsd() ).xyz( elem.second.atomno() ) );

		Vector const f2( xyz1 - xyz2 );
		if ( f2.length() != 0.0 ) {
			Vector const f1( xyz1.cross( xyz2 ) );
			our_f1 += 2 * f1;
			our_f2 += 2 * f2;
		}
	}

	our_f1 *= factor;
	our_f2 *= factor;

	F1 += our_f1;
	F2 += our_f2;
}


/// @brief Indicates the current version of this score term
core::Size
AlignmentEnergy::version() const {
	return 1;
}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
void
AlignmentEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	// OK. The energy ought to be a FlatHarmonicFunc of RMSD to the reference pose.
	// No penalty from 0 to rmsd_screen_.
	// wtf might
	/*moving_res_list_ = */

	utility::vector1< Size > root_partition_res = stepwise::modeler::figure_out_root_partition_res( pose, core::pose::full_model_info::get_moving_res_from_full_model_info( pose ) );
	//TR << "root_partition_res is " << root_partition_res << std::endl;
	pose_aligner_->set_root_partition_res( root_partition_res );

	//Pose pose_save = pose;
	// don't apply!!!
	pose_aligner_->initialize( pose );
	if ( pose_aligner_->superimpose_atom_id_map().size() == 0 ) return;
	Real const rmsd = core::scoring::rms_at_corresponding_atoms_no_super( pose, *align_pose_, pose_aligner_->superimpose_atom_id_map() );

	totals[ alignment ] = func_->func( rmsd );
}

core::scoring::methods::EnergyMethodOP
AlignmentEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & options ) const {
	// return bunk one maybe for now...
	return core::scoring::methods::EnergyMethodOP( new AlignmentEnergy( options ) );
}

core::scoring::ScoreTypes
AlignmentEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes st;
	st.push_back( alignment );
	return st;
}


} //protocols
} //rna

