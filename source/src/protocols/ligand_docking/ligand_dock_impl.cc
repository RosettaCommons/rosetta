// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ligand_dock_impl.cc
///
/// @brief
/// @author Rocco Moretti (rmoretti@u.washington.edu); Ian Davis (ian.w.davis@gmail.com)

// Unit Headers

#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/moves/Mover.hh>

// Package Headers

#include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/jd2/util.hh>

#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>

// Project Headers

#include <core/types.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //for addding constraints if demanded by user
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableString.hh>

// Utility Headers

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer TR( "protocols.ligand_docking.main" );

class LigandDockMain : public protocols::ligand_docking::LigandDockProtocol {
public:
	LigandDockMain():
		LigandDockProtocol()
	{
		if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			native_ = core::pose::PoseOP( new core::pose::Pose() );
			core::import_pose::pose_from_pdb( *native_, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() );
		}

		if( basic::options::option[ basic::options::OptionKeys::enzdes::cstfile].user() ){
			//we need the residue type set, assuming FA standard is used
			core::chemical::ResidueTypeSetCAP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
			basic::options::option[ basic::options::OptionKeys::run::preserve_header ].value(true);
			constraints_ = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO( restype_set ) );
			constraints_->read_enzyme_cstfile( basic::options::option[ basic::options::OptionKeys::enzdes::cstfile ] );

			core::scoring::ScoreFunctionOP scorefunc = this->get_scorefxn();
			if (scorefunc->has_zero_weight( core::scoring::coordinate_constraint ) ) scorefunc->set_weight(core::scoring::coordinate_constraint, 1.0 ) ;
			if (scorefunc->has_zero_weight( core::scoring::atom_pair_constraint ) ) scorefunc->set_weight(core::scoring::atom_pair_constraint, 1.0 ) ;
			if (scorefunc->has_zero_weight( core::scoring::angle_constraint ) ) scorefunc->set_weight(core::scoring::angle_constraint, 1.0 ) ;
			if (scorefunc->has_zero_weight( core::scoring::dihedral_constraint ) ) scorefunc->set_weight(core::scoring::dihedral_constraint, 1.0 ) ;
		}

	}

	virtual
	~LigandDockMain(){};

	virtual
	void
	apply( core::pose::Pose & pose ) {

		if( constraints_ ){
			constraints_->add_constraints_to_pose( pose, this->get_scorefxn(), false );
		}

		LigandDockProtocol::apply( pose );

		// Score new structure and add to silent file
		std::map< std::string, core::Real > scores;
		core::import_pose::atom_tree_diffs::map_of_weighted_scores(pose, *(this->get_scorefxn()), scores);
		if( native_ ) {
			this->append_ligand_docking_scores(*native_, pose, this->get_scorefxn(), scores, constraints_);
		} else { // input serves as native
			core::pose::PoseCOP input_pose(protocols::jd2::JobDistributor::get_instance()->current_job()->get_pose());
			this->append_ligand_docking_scores(*input_pose, pose,	this->get_scorefxn(), scores, constraints_);
		}

		// Attach relevant scores to current job (Assumes use of JD2)
		{
		protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );
		std::map< std::string, core::Real >::const_iterator curr( scores.begin() );
		std::map< std::string, core::Real >::const_iterator end( scores.end() );
		for (; curr != end; ++curr) {
			curr_job->add_string_real_pair(curr->first, curr->second);
		}
		}

		/* TODO: From old jd1 code:
		// Want to recycle the reference poses from the input silent file, or else every entry becomes a new reference pose!
		if( use_silent_in ) {
			//std::cout << "Ref addr " << (atdiff->ref_pose_for( curr_job->input_tag() ))() << std::endl;
			jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *(atdiff->ref_pose_for( curr_job->input_tag() )), *the_pose );
		}
		else jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *input_pose, *the_pose );
		*/


	}

	virtual
	std::string
	get_name() const {
		return "LigandDockMain";
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new LigandDockMain );
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

private:

	core::pose::PoseOP native_;
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP constraints_;
};

typedef utility::pointer::shared_ptr< LigandDockMain > LigandDockMainOP;

} //namespace ligand_docking
} //namespace protocols

//////////////////////////////////////////////////////////////////////////////
/// @details Assumes option system has already been initialized!
int
ligand_dock_main()
{
	protocols::ligand_docking::LigandDockMainOP ligand_dock( new protocols::ligand_docking::LigandDockMain );

	// Change precision of AtomTreeDiffJobOutputter, for backward compatability reasons
	// makes silent file much smaller, ~3x vs. default 6,4,2
	// but has too much atom positioning error (~0.1 Ang) to be reliable for rescoring
	if ( ! basic::options::option[ basic::options::OptionKeys::out::file::atom_tree_diff_bb].user() &&
			! basic::options::option[ basic::options::OptionKeys::out::file::atom_tree_diff_sc].user() &&
			! basic::options::option[ basic::options::OptionKeys::out::file::atom_tree_diff_bl].user()) {
		protocols::ligand_docking::TR << "Changing atom_tree_diff output precision to 6,3,1." << std::endl;
		basic::options::option[ basic::options::OptionKeys::out::file::atom_tree_diff_bb].value(6);
		basic::options::option[ basic::options::OptionKeys::out::file::atom_tree_diff_sc].value(3);
		basic::options::option[ basic::options::OptionKeys::out::file::atom_tree_diff_bl].value(1);
	}

 	protocols::jd2::JobOutputterOP 	job_outputter(protocols::jd2::JobDistributor::get_instance()->job_outputter());
	protocols::jd2::AtomTreeDiffJobOutputterOP atd_outputter( utility::pointer::dynamic_pointer_cast< protocols::jd2::AtomTreeDiffJobOutputter > ( job_outputter ) );
	if ( atd_outputter ) {
		atd_outputter->use_input_for_ref(true); // match original ligand_dock application
	}

	protocols::jd2::JobDistributor::get_instance()->go(ligand_dock);

	return 0;
}

