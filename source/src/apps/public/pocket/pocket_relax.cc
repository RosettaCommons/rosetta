// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/import_pose/import_pose.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/util.hh>

#include <core/types.hh>

#include <devel/init.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/pockets/PocketGrid.hh>

#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace core::scoring;

using utility::vector1;

static basic::Tracer TR( "pocket_relax" );
OPT_KEY( String, exemplar_target_pdb_num )

//This is copy/pasted from src/apps/public/minimize.cc, but has had the constraints removed
class NCMinimize : public moves::Mover {

public:
	NCMinimize();

	~NCMinimize() override;

	MoverOP clone() const override;
	MoverOP fresh_instance() const override;

	void apply( Pose & pose ) override;
	std::string get_name() const override;
	void test_move( Pose & pose ) override
	{
		apply(pose);
	}

private:
	ScoreFunctionOP score_function_;
};

NCMinimize::NCMinimize() :
	Mover( "benchmark" ),
	score_function_( get_score_function() )
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

}

NCMinimize::~NCMinimize() = default;

MoverOP NCMinimize::clone() const {
	return MoverOP( new NCMinimize( *this ) );
}

MoverOP NCMinimize::fresh_instance() const {
	return MoverOP( new NCMinimize );
}

void NCMinimize::apply( Pose & pose ) {
	using namespace pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	score_function_->set_weight ( pocket_constraint, 0);
	std::string min_type =  option[ OptionKeys::run::min_type ]();
	core::Real min_tol =  option[ OptionKeys::run::min_tolerance ]();
	core::optimization::MinimizerOptions options( min_type, min_tol, true, false );
	core::kinematics::MoveMap final_mm;
	final_mm.set_chi( true );
	final_mm.set_bb( true );

	/*core::Real start_score =*/ (*score_function_)(pose);
	core::Size repeats = 1;
	for ( core::Size i = 0; i < repeats; i++ ) {
		core::optimization::AtomTreeMinimizer().run( pose, final_mm, *score_function_, options );
		TR << "Score: " << i << "  " <<  (*score_function_)(pose) << std::endl;
	}

	core::Real final_score = (*score_function_)(pose);
	TR << "FinalScore: " << final_score << std::endl;

}

std::string NCMinimize::get_name() const {
	return "NCMinimize";
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace core;
		using namespace protocols::moves;
		using namespace scoring;
		using namespace basic::options;
		using namespace protocols::jobdist;
		using namespace protocols::relax;
		using namespace basic::options::OptionKeys;
		using protocols::moves::MoverOP;

		NEW_OPT( exemplar_target_pdb_num, "Target residue(s) for exemplar generation", "");






		protocols::relax::ClassicRelax::register_options();
		jd2::register_options();
		option.add_relevant( OptionKeys::in::file::fullatom );
		option.add_relevant( OptionKeys::in::file::movemap );
		option.add_relevant( OptionKeys::relax::fast );
		devel::init(argc, argv);

		std::string const resid_c = option[exemplar_target_pdb_num];

		//validate that the exemplar_target_pdb_num is properly formatted
		if ( resid_c.length() ) {
			const std::string & delimiters = ";";
			std::string::size_type lastPos = resid_c.find_first_not_of(delimiters, 0);

			// Find first "non-delimiter".
			std::string::size_type pos = resid_c.find_first_of(delimiters, lastPos);
			while ( std::string::npos != pos || std::string::npos != lastPos ) {
				std::string const & resid_list = resid_c.substr( lastPos, pos - lastPos );
				lastPos = resid_c.find_first_not_of(delimiters, pos);
				pos = resid_c.find_first_of(delimiters, lastPos);

				const std::string & subdelimiters = ",";
				std::string::size_type lastSubPos = resid_list.find_first_not_of(subdelimiters, 0);
				// Find first "non-delimiter".
				std::string::size_type subpos = resid_list.find_first_of(subdelimiters, lastSubPos);
				while ( std::string::npos != subpos || std::string::npos != lastSubPos ) {
					std::string const & resid = resid_list.substr( lastSubPos, subpos - lastSubPos );
					lastSubPos = resid_list.find_first_not_of(subdelimiters, subpos);
					subpos = resid_list.find_first_of(subdelimiters, lastSubPos);
					std::size_t fpos( resid.find(':') );
					if ( fpos != std::string::npos ) {
						if ( fpos != resid.length() -2 || fpos == 0 ) {
							throw ( CREATE_EXCEPTION(utility::excn::BadInput, " invalid exemplar_target_pdb_num " + resid_c));
						}
						char chain = resid[resid.length() -1];
						if ( !( (chain >='A' && chain <='Z') || (chain >='a' && chain <='z') ) ) {
							throw ( CREATE_EXCEPTION(utility::excn::BadInput, " invalid exemplar_target_pdb_num " + resid_c));
						}
						for ( int i = 0; i< (int) fpos; ++i ) {
							if ( !(resid[i] >='0' && resid[i] <='9') ) {
								throw ( CREATE_EXCEPTION(utility::excn::BadInput, " invalid exemplar_target_pdb_num " + resid_c));
							}
						}
					} else {
						for ( char i : resid ) {
							if ( !(i >='0' && i <='9') ) {
								throw ( CREATE_EXCEPTION(utility::excn::BadInput, " invalid exemplar_target_pdb_num " + resid_c));
							}
						}
					}
				}
			}
		}


		//return relax::Relax_main( false );

		protocols::moves::MoverOP protocol = generate_relax_from_cmd();
		protocols::jd2::set_native_in_mover( *protocol );

		// add constraints from cmd line
		if ( option[ OptionKeys::constraints::cst_fa_file ].user() || option[ OptionKeys::constraints::cst_file ].user() ) {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
			if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
			} else {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_file_option() );
			}
			seqmov->add_mover( loadCsts );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		// set pose for density scoring if a map was input
		// (potentially) dock map into density
		if ( option[ OptionKeys::edensity::mapfile ].user() ) {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			seqmov->add_mover( MoverOP( new protocols::electron_density::SetupForDensityScoringMover ) );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		// setup symmetry mover ... this should happen _before_ SetupForDensityScoringMover
		//   to avoid adding extra VRTs
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			seqmov->add_mover( MoverOP( new protocols::simple_moves::symmetry::SetupForSymmetryMover ) );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		// superimpose input model to the native structure or to a supplied PDB file
		if ( option[ OptionKeys::relax::superimpose_to_file ].user() ||
				option[ OptionKeys::relax::superimpose_to_native ].user()
				) {
			core::pose::Pose ref_pose;
			std::string ref_filename;
			if (  option[ OptionKeys::relax::superimpose_to_file ].user() ) ref_filename = option[ basic::options::OptionKeys::relax::superimpose_to_file ]();
			if (  option[ OptionKeys::relax::superimpose_to_native ].user() ) ref_filename =  option[ basic::options::OptionKeys::in::file::native ]();
			core::import_pose::pose_from_file( ref_pose, ref_filename , core::import_pose::PDB_file);
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			protocols::simple_moves::SuperimposeMoverOP sm( new protocols::simple_moves::SuperimposeMover );
			sm->set_reference_pose( ref_pose );
			seqmov->add_mover( sm );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
		seqmov->add_mover( protocol );


		MoverOP minprotocol( new NCMinimize() );
		seqmov->add_mover( minprotocol );

		protocol = seqmov;

		protocols::jd2::JobDistributor::get_instance()->go( protocol );


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
