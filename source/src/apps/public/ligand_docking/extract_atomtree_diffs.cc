// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/public/ligand_dockiing/extract_atomtree_diffs.cc
///
/// @brief  Convert atomtree_diff silent files into normal PDBs.
/// @author Rocco Moretti (rmoretti@u.washington.edu); Ian Davis (ian.w.davis@gmail.com)

// Note: the original (JD1) extract_atomtree_diffs application can be found (at least for the time being)
// as apps/pilot/rmoretti/extract_atomtree_diffs_jd1.cc

// Unit Headers

// Package Headers

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

// Project Headers

#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

//for addding constraints if demanded by user
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/jd2/AtomTreeDiffJobInputter.hh>
#include <core/chemical/ChemicalManager.hh>

// Utility Headers

#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

// Numeric Headers

// C++ headers

static basic::Tracer TR( "extract_atomtree_diffs.main" );

class ExtractATD : public protocols::moves::Mover {
public:
	ExtractATD() {

		// Create a protocol object just to steal its score function.
		scorefxn_ = protocols::ligand_docking::LigandBaseProtocol().get_scorefxn();
	}

	virtual
	~ExtractATD(){};

	virtual
	void
	apply( core::pose::Pose & pose ) {

		core::pack::dunbrack::load_unboundrot( pose ); // adds scoring bonuses for the "unbound" rotamers, if any
		(*scorefxn_)( pose );
	}

	virtual
	std::string
	get_name() const {
		return "ExtractATD";
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new ExtractATD );
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	core::scoring::ScoreFunctionOP
	scorefxn() { return scorefxn_; }

private:

	core::scoring::ScoreFunctionOP scorefxn_;
};

typedef utility::pointer::shared_ptr< ExtractATD > ExtractATDOP;


int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);
		protocols::jd2::register_options();

		// Munge the input options to be compatible with the previous versions of extract_atomtree_diffs
		if ( basic::options::option[ basic::options::OptionKeys::in::file::l ].user() || basic::options::option[ basic::options::OptionKeys::in::file::list ].user() ) {
			TR.Error << "Sorry, cannot use -in:file:l or -in::file::list with extract_atomtree_diffs" << std::endl;
			utility_exit_with_message("Invalid option -in:file:l or -in::file::list passed to extract_atomtree_diffs.");
		}

		if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
			if ( basic::options::option[ basic::options::OptionKeys::in::file::atom_tree_diff ].user() ) {
				TR.Error << "Sorry, cannot use both -in:file:s and -in::file::atom_tree_diff with extract_atomtree_diffs" << std::endl;
				utility_exit_with_message("Invalid option combination -in:file:s and -in::file::atom_tree_diff passed to extract_atomtree_diffs.");
			} else {
				TR << "Converting -in:file:s to -in::file::atom_tree_diff." << std::endl;
				basic::options::option[ basic::options::OptionKeys::in::file::atom_tree_diff ].value(basic::options::option[ basic::options::OptionKeys::in::file::s ]);
				basic::options::option[ basic::options::OptionKeys::in::file::s ].deactivate();
			}
		}

		// Match output tagging status of previous version of extract_atomtree_diffs, unless overwritten
		if ( ! basic::options::option[ basic::options::OptionKeys::out::no_nstruct_label ].user() ) {
			basic::options::option[ basic::options::OptionKeys::out::no_nstruct_label ].value(true);
			TR << "Setting no_nstruct_label" << std::endl;
		}

		if ( ! basic::options::option[ basic::options::OptionKeys::run::preserve_header].user() ) {
			basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);
		}

		ExtractATDOP extract_mover( new ExtractATD );

		// Setup constraints on reference poses
		// Reading the enzdes constraints creates new residue types we may need, in the case of covalent constraints.
		protocols::jd2::JobInputterOP inputter(protocols::jd2::JobDistributor::get_instance()->job_inputter());
		protocols::jd2::AtomTreeDiffJobInputterOP atd_inputter( utility::pointer::dynamic_pointer_cast< protocols::jd2::AtomTreeDiffJobInputter > ( inputter ) );
		if ( atd_inputter && basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ) {
			protocols::toolbox::match_enzdes_util::EnzConstraintIOOP constraints;
			//we need the residue type set, assuming FA standard is used
			core::chemical::ResidueTypeSetCOP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
			constraints = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO( restype_set ) );
			constraints->read_enzyme_cstfile(basic::options::option[basic::options::OptionKeys::enzdes::cstfile]);

			utility::vector1< core::pose::PoseOP > const & ref_poses = atd_inputter->all_ref_poses();
			for ( core::Size i = 1; i <= ref_poses.size(); ++i ) {
				constraints->add_constraints_to_pose( *ref_poses[i], extract_mover->scorefxn(), false );
			}
		}


		protocols::jd2::JobDistributor::get_instance()->go(extract_mover);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

