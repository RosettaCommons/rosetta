// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/Align_RmsdEvaluator.hh>
#include <protocols/evaluation/StructuralSimilarityEvaluator.hh>
#include <protocols/evaluation/Align_RotamerEvaluator.hh>
#include <protocols/evaluation/RDC_Evaluator.hh>
#include <protocols/evaluation/ScoreEvaluator.hh>
#include <protocols/evaluation/JScoreEvaluator.hh>
#include <protocols/evaluation/ContactMapEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
#include <protocols/evaluation/ConstraintEvaluator.hh>
#include <protocols/evaluation/CombinedConstraintEvaluator.hh>
#include <protocols/evaluation/ChemicalShiftEvaluator.hh>
#include <protocols/evaluation/RPF_ScoreEvaluator.hh>
#include <protocols/evaluation/CamShiftEvaluator.hh>
#include <protocols/evaluation/PalesEvaluator.hh>
#include <protocols/evaluation/PredictedBurialEvaluator.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>

// Project Headers
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentStructFactory.hh>
//// C++ headers

// due to template function
#include <core/io/silent/SilentStruct.hh>


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/fragment/SecondaryStructure.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/evaluation/util.hh>
#include <utility/vector0.hh>

static basic::Tracer tr("protocols.evalution");

namespace protocols {
namespace evaluation {
using namespace core;

static bool options_registered_=false;

void register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;
	OPT( evaluation::rmsd );
	OPT( evaluation::gdtmm );
	OPT( evaluation::rdc );
	OPT( evaluation::pool );
	OPT( evaluation::constraints );
	OPT( in::file::native );
}

void invert_include_residues( Size nres, core::scoring::ResidueSelectionVector const& include_list, core::scoring::ResidueSelectionVector& exclude_list ) {

	exclude_list.clear();

	for ( Size ir = 1; ir <= nres; ++ir ) {
		bool include_residue = false;
		for ( Size ex = 1; ex <= include_list.size(); ex ++ ) {
			if ( include_list[ex] == ir ) {
				include_residue = true;
				break;
			}
		}

		if ( !include_residue ) {
			exclude_list.push_back( ir );
		}
	} // for ( Size ir = 1; ir <= native_pose.total_residue(); ++ir )
}


void read_common_evaluator_options( MetaPoseEvaluator& eval ) {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  if ( option[ OptionKeys::evaluation::rmsd_target ].user() ) {
    /*
      this creates Evaluators to get RMSD and GDTMM for as many structures as you specify in this FileVector option.
      pls: provide also as many column-names to match your rmsd-targets.
      ---
      missing density in the input file will be ignored, e.g., this is your way to signal on which residues the rmsd should be computed
    */
    utility::vector1< std::string > const & rmsd_target  ( option[ OptionKeys::evaluation::rmsd_target ]() );
		utility::vector1< std::string >         rmsd_col_name;

		if ( option[ OptionKeys::evaluation::rmsd_column ].user() ){
			rmsd_col_name = option[ OptionKeys::evaluation::rmsd_column ]();
		}else{
			// make up default names: firone is just empty, leading to the columns being correctly named "rms", "gdtmm" etc.. the
			// subsequent columns then become "_2", "_3" etc..
			rmsd_col_name.push_back("");
			for( core::Size j=2; j<=rmsd_target.size(); ++j){
				rmsd_col_name.push_back("_" + ObjexxFCL::string_of(j));
			}
		}
    for ( Size ct = 1; ct <= rmsd_target.size(); ct ++ ) {
      pose::PoseOP rmsd_pose = new pose::Pose;
      core::import_pose::pose_from_pdb( *rmsd_pose, rmsd_target[ ct ] );
      std::string tag( ObjexxFCL::string_of( ct ) );
      if ( rmsd_col_name.size() >= ct ) tag = rmsd_col_name[ ct ];
      eval.add_evaluation( new evaluation::SelectRmsdEvaluator( rmsd_pose, tag ) );
      if ( option[ OptionKeys::evaluation::gdtmm ]() ) eval.add_evaluation( new evaluation::SelectGdtEvaluator( rmsd_pose, tag ) );
			if ( option[ OptionKeys::evaluation::score_with_rmsd ]() ){
				core::scoring::ResidueSelection selection;
				find_existing_residues( rmsd_pose, tag, selection );
				core::scoring::ResidueSelectionVector vector;
				std::copy( selection.begin(), selection.end(), std::back_inserter( vector ) );
				eval.add_evaluation( new evaluation::TruncatedScoreEvaluator( tag, vector ) );
			}
			if ( option[ OptionKeys::evaluation::symmetric_rmsd ]() ) {
				eval.add_evaluation( new evaluation::SymmetricRmsdEvaluator( rmsd_pose, tag ) );
			}
    }
  }

	if ( option[ OptionKeys::evaluation::predicted_burial_fn ].user() ) {
		std::string const fn( option[ OptionKeys::evaluation::predicted_burial_fn ]() );
		eval.add_evaluation( new PredictedBurialEvaluator(fn) );
	}

	if ( option[ OptionKeys::evaluation::align_rmsd_target ].user() ) {
		using std::string;
		using utility::vector1;
    vector1< string > const & align_rmsd_target(
			option[ OptionKeys::evaluation::align_rmsd_target ]()
		);
    vector1< string > const & align_rmsd_col_names(
			option[ OptionKeys::evaluation::align_rmsd_column ]()
		);

    vector1< string > align_rmsd_fns;
		if ( option[ OptionKeys::evaluation::align_rmsd_fns ].user() ) {
			align_rmsd_fns = option[ OptionKeys::evaluation::align_rmsd_fns ]();
		}
		runtime_assert( align_rmsd_target.size() == align_rmsd_col_names.size() );
		for ( Size ii = 1; ii <= align_rmsd_target.size(); ++ii ) {
      pose::PoseOP rmsd_pose = new pose::Pose;
      core::import_pose::pose_from_pdb( *rmsd_pose, align_rmsd_target[ii] );
			//string const tag( align_rmsd_target[ii] );
			string const tag( align_rmsd_col_names[ii] );
			core::sequence::SequenceAlignmentOP aln(0);
			if ( align_rmsd_fns.size() >= ii ) {
				*aln = core::sequence::read_aln(
					align_rmsd_fns[ii], option[ OptionKeys::evaluation::align_rmsd_format ]()
				).front();
			}
			eval.add_evaluation( new Align_RmsdEvaluator(rmsd_pose,tag,true,aln) );
			eval.add_evaluation( new Align_RotamerEvaluator(rmsd_pose,tag,true,aln) );
		}
	}

	if ( option[ OptionKeys::evaluation::structural_similarity ].user() ) {
		using std::string;
		using core::import_pose::pose_stream::SilentFilePoseInputStream;

		SilentFilePoseInputStream silent_input(
			option[ OptionKeys::evaluation::structural_similarity ]()
		);
		core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::rsd_set_from_cmd_line();
		utility::vector1< core::pose::Pose > poses;
		while( silent_input.has_another_pose() ) {
			core::pose::Pose pose;
			silent_input.fill_pose( pose, *rsd_set );
			poses.push_back( pose );
		}
		eval.add_evaluation( new StructuralSimilarityEvaluator(poses) );
	}

	if ( option[ OptionKeys::evaluation::jscore_evaluator ].user() ) {
		using std::string;
		using utility::vector1;
		using core::scoring::ScoreFunctionOP;
		vector1< string > const & tags( option[ OptionKeys::evaluation::jscore_evaluator ]() );

		for ( Size ii = 1; ii <= tags.size() - 1; ii += 2 ) {
			string scorefxn_name( tags[ii]   );
			string rsd_set_name ( tags[ii+1] );

			eval.add_evaluation( new JScoreEvaluator( scorefxn_name, rsd_set_name ) );
		}
	}

	core::pose::PoseOP native_pose = NULL;
	if ( option[ in::file::native ].user() ) {
		native_pose = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
	}

	if ( option[ OptionKeys::evaluation::contact_map ] ) {
		if ( !native_pose ) {
			tr.Error << "Error: -evaluation::contact_map must be specified with a native!\n";
		} else {
			core::Real max_dist(12);
			core::Size min_seqsep(5);
			eval.add_evaluation(
				new ContactMapEvaluator( *native_pose, max_dist, min_seqsep )
			);
		}
	}

	// set rmsd native
	if ( native_pose && option[ in::file::native ].user() ) {
		if ( option[ in::file::native_exclude_res ].user() ) {
			eval.add_evaluation( new SelectRmsdEvaluator(
						native_pose,
						core::scoring::invert_exclude_residues( native_pose->total_residue(), option[ in::file::native_exclude_res ]()),
						"" )
			);
			if ( option[ OptionKeys::evaluation::gdtmm ]() ) {
				eval.add_evaluation( new SelectGdtEvaluator(
						native_pose,
						core::scoring::invert_exclude_residues( native_pose->total_residue(), option[ in::file::native_exclude_res ]()),
						"" )
				);
			}
		} else if ( option[ OptionKeys::abinitio::rmsd_residues ].user() ){
			core::Size start = option[ OptionKeys::abinitio::rmsd_residues ]()[ 1 ];
			Size end = option[ OptionKeys::abinitio::rmsd_residues ]()[ 2 ];
			eval.add_evaluation( new RmsdEvaluator( native_pose, start, end,  "", option[ OptionKeys::abinitio::bGDT ]() ) );
		} else {
			eval.add_evaluation( new SelectRmsdEvaluator( native_pose, "" ) );
			if ( option[ OptionKeys::evaluation::gdtmm ]() ) eval.add_evaluation( new SelectGdtEvaluator( native_pose, "" ) );
      eval.add_evaluation( new SelectMaxsubEvaluator( native_pose, "" ) );
		}
	} // if ( native_pose_ )

	if ( option[ OptionKeys::evaluation::pool ].user() ) {
		using namespace protocols::moves::mc_convergence_checks;
		Pool_RMSD_OP pool_ptr = new Pool_RMSD( option[ OptionKeys::evaluation::pool ]() );
		eval.add_evaluation( new Pool_Evaluator( pool_ptr ) );
	}

	if ( option[ OptionKeys::evaluation::rmsd ].user() ) {
    typedef utility::vector1< std::string > RmsdVector;
		RmsdVector const& rmsd( option[ OptionKeys::evaluation::rmsd ]() );

		//the ultimate in RMSD technology:
		//you can have any number of tripletts:
		// each triplett contains  the
		//    target: file or special tags: NATIVE (use in:file:native ) IRMS (use structure in Job object)
		///   columnname: XXX creates: rms_XXX  gdtmm_XXX
		///   selection: file or special tag: FULL
		///   modifier: EXCLUDE ( take the inverse of the selection )
		///                INLINE r1 r2 r3 ... rm END_INLINE   -- select residues r1, r2, r3, ..., rm
		for ( RmsdVector::const_iterator it=rmsd.begin(); it!=rmsd.end(); ++it ) {
			core::pose::PoseOP target_pose = NULL;
			std::string fname( *it );
			std::string column;
			std::string selection_file;
			bool CA( true );
			bool invert( false );
			bool loop_rms( false );
			bool superimpose_for_looprms( false );
			++it;
			loops::Loops core;
			while ( it != rmsd.end() && (*it).find('=')!=std::string::npos ) {
				Size pos( it->find('=') );
				std::string key=it->substr(0,pos);
				std::string value=it->substr(pos+1);
				if ( key=="heavy" && value=="yes" ) {
					tr.Info << "switching to full-atom RMSD " << std::endl;
					CA=false;
				} else if ( key == "loop") {
					loop_rms = true;
				} else if ( key == "superimpose") {
					superimpose_for_looprms = (value=="yes");
				} else if ( key == "core" ) {
					core.read_loop_file( value, false, "RIGID" );
			} else {
					utility_exit_with_message( "key not recognized: "+key+" possible keys: { heavy }" );
				}
				tr.Info << "detected option: " << key << " : " << value << std::endl;
				++it;
			}
			if ( it != rmsd.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
							 "need to specify tripletts <target> [key/value pairs] <column> <selection/FULL> with option -evaluation:rmsd   last read: "+fname );
			}
			++it;
			if ( it != rmsd.end() ) {
				selection_file = *it;
			} else {
				utility_exit_with_message(
							 "need to specify tripletts <target> <column> <selection/FULL> with option -evaluation:rmsd   last read: "+column );
			}
			if ( fname == "NATIVE" ) target_pose = native_pose;
			else if ( fname != "IRMS" ) {
				target_pose = new pose::Pose;
				core::import_pose::pose_from_pdb( *target_pose, fname );
			}

			if ( selection_file == "EXCLUDE" ) {
				invert = true;
				++it;
				if ( it != rmsd.end() ) {
					selection_file = *it;  //read next tag
				} else {
					utility_exit_with_message(
						"need to specify a <selection/FULL> after 'EXCLUDE' with option -evaluation:rmsd    last read: "+column );
				}
			}
			loops::Loops loops;
			utility::vector1< Size> selection; //figure out selection
			if ( selection_file == "INLINE" ) {
				std::string next_tag( "" );
				++it;
				if ( it != rmsd.end() ) {
					next_tag = *it;
				} else {
					utility_exit_with_message(
						"need to find END_INLINE after INLINE in option -evaluation:rmsd    last read: "+column );
				} //error condition
				while( next_tag != "END_INLINE" ) {
					selection.push_back( utility::string2int( next_tag ) );
					++it;
					if ( it != rmsd.end() ) {
						next_tag = *it;
					} else {
						utility_exit_with_message(
							"need to find END_INLINE after INLINE in option -evaluation:rmsd    last read: "+column );
					} //error condition
				}
			} else if ( selection_file != "FULL" ) {
				if ( loop_rms ) {
					loops.read_loop_file( selection_file, false, "LOOP" );
				} else {
					loops::Loops core;
					core.read_loop_file( selection_file, false, "RIGID" );
					core.get_residues( selection );
				}
			}
			if ( invert ) {
				utility::vector1< Size > inverted_selection;
				for ( Size i = 1; i<=target_pose->total_residue(); ++i) {
					bool found( false );
					for ( Size j = 1; j<=selection.size() && !found; ++j ) {
						if ( selection[ j ]==i ) {
							found = true;
						}
					}
					if ( !found) inverted_selection.push_back( i );
				}
				selection = inverted_selection;
			}
			if ( selection_file != "FULL" ) {
				if ( !loop_rms ) {
					eval.add_evaluation( new evaluation::SelectRmsdEvaluator( target_pose, selection, column, CA ) );
				} else {
					eval.add_evaluation( new evaluation::LoopRmsdEvaluator( target_pose, loops, core, column, CA, superimpose_for_looprms ) );
				}
				if ( option[ OptionKeys::evaluation::gdtmm ]() ) {
					eval.add_evaluation( new evaluation::SelectGdtEvaluator( target_pose, selection, column) );
				}
			} else {
				eval.add_evaluation( new evaluation::SelectRmsdEvaluator( target_pose, column, CA ) );
				if ( option[ OptionKeys::evaluation::gdtmm ]() ) {
					eval.add_evaluation( new evaluation::SelectGdtEvaluator( target_pose, column ) );
				}
			} // no selection
		} // iterate over tripletts in option -rmsd
	} // option -rmsd

	if ( option[ OptionKeys::evaluation::rmsd_select ].user() ) {
    utility::vector1< utility::file::FileName > const& rmsd_core( option[ OptionKeys::evaluation::rmsd_select ]() );
		if ( !option[ in::file::native ].user() ) utility_exit_with_message( "need to specify in:file:native together with rmsd_select " );

	  for ( Size ct = 1; ct <= rmsd_core.size(); ct ++ ) {
			loops::Loops core;
			core.read_loop_file( rmsd_core[ ct ], false, "RIGID" );
			utility::vector1< Size> selection;
			core.get_residues( selection );
			if ( native_pose ) eval.add_evaluation( new evaluation::SelectRmsdEvaluator( native_pose, selection, rmsd_core[ ct ].base() ) );
			if ( option[ OptionKeys::evaluation::score_with_rmsd ]() ) {
				eval.add_evaluation( new evaluation::TruncatedScoreEvaluator( rmsd_core[ ct ].base(), selection ) );
			}
		}
	}

  if ( option[ OptionKeys::evaluation::rdc_target ].user() ) {
    /*
      this creates Evaluators to get RMSD and GDTMM for as many structures as you specify in this FileVector option.
      pls: provide also as many column-names to match your rmsd-targets.
      ---
      missing density in the input file will be ignored, e.g., this is your way to signal on which residues the rmsd should be computed
    */
    utility::vector1< std::string > const& rdc_target( option[ OptionKeys::evaluation::rdc_target ]() );
    utility::vector1< std::string > const& rdc_col_name( option[ OptionKeys::evaluation::rdc_column ]() );
    for ( Size ct = 1; ct <= rdc_target.size(); ct ++ ) {
      pose::PoseOP rdc_pose = new pose::Pose;
      core::import_pose::pose_from_pdb( *rdc_pose, rdc_target[ ct ] );
      std::string tag( ObjexxFCL::string_of( ct ) );
      if ( rdc_col_name.size() >= ct ) tag = rdc_col_name[ ct ];
      eval.add_evaluation( new evaluation::SelectRDC_Evaluator( rdc_pose, tag ) );
    }
  }

	if ( option[ OptionKeys::evaluation::rdc_select ].user() ) {
    utility::vector1< utility::file::FileName > const& rdc_core( option[ OptionKeys::evaluation::rdc_select ]() );

	  for ( Size ct = 1; ct <= rdc_core.size(); ct ++ ) {
			loops::Loops core;
			core.read_loop_file( rdc_core[ ct ], false, "RIGID" );
			utility::vector1< Size> selection;
			core.get_residues( selection );
			eval.add_evaluation( new evaluation::SelectRDC_Evaluator( selection, rdc_core[ ct ].base() ) );
		}
	}

	if ( option[ OptionKeys::evaluation::rdc ].user() ) {
    typedef utility::vector1< std::string > RdcVector;
		RdcVector const& rdc( option[ OptionKeys::evaluation::rdc ]() );
		utility::vector1< core::Size> empty_selection;
		for ( RdcVector::const_iterator it=rdc.begin(); it!=rdc.end(); ++it ) {
			std::string fname( *it );
			std::string column;
			++it;
			if ( it != rdc.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
							 "need to specify dupletts <rdcs> <column> with option -evaluation:rdc   last read: "+fname );
			}
			eval.add_evaluation( new evaluation::SelectRDC_Evaluator( empty_selection, column, fname ) );
		} // iterate over tripletts in option -rmsd
	}


  if ( option[ OptionKeys::evaluation::constraints ].user() ) {
    /*
      this creates Evaluators to evaluate different constraint sets against your decoys
      pls: provide also as many column-names to match your constraint sets
      ---
    */
    utility::vector1< std::string > const& cst_target( option[ OptionKeys::evaluation::constraints ]() );
    utility::vector1< std::string > const& cst_col_name( option[ OptionKeys::evaluation::constraints_column ]() );
    for ( Size ct = 1; ct <= cst_target.size(); ct ++ ) {
      std::string tag( ObjexxFCL::string_of( ct ) );
      if ( cst_col_name.size() >= ct ) tag = cst_col_name[ ct ];
      eval.add_evaluation( new evaluation::ConstraintEvaluator( tag, cst_target[ ct ] ) );
    }
  }

  if ( option[ OptionKeys::evaluation::combined_constraints ].user() ) {
    utility::vector1< std::string > const& cst_target( option[ OptionKeys::evaluation::combined_constraints ]() );
    utility::vector1< std::string > const& cst_col_name( option[ OptionKeys::evaluation::combined_constraints_column ]() );
    for ( Size ct = 1; ct <= cst_target.size(); ct ++ ) {
      std::string tag( ObjexxFCL::string_of( ct ) );
      if ( cst_col_name.size() >= ct ) tag = cst_col_name[ ct ];
      eval.add_evaluation( new evaluation::CombinedConstraintEvaluator( tag, cst_target[ ct ], 2, option[ OptionKeys::evaluation::combine_statistics ] ) );
    }
  }

// 	if ( option[ OptionKeys::evaluation::rpf ].user() ) {
// 		eval.add_evaluation( new evaluation::RPF_ScoreEvaluator( "rpf_score", 5 ) );
// 	}

	if ( option[ OptionKeys::evaluation::chemical_shifts ].user() ) {
    typedef utility::vector1< std::string > CSVector;
		CSVector const& cs_shifts( option[ OptionKeys::evaluation::chemical_shifts ]() );

		for ( CSVector::const_iterator it=cs_shifts.begin(); it!=cs_shifts.end(); ++it ) {
			std::string fname( *it );
			std::string column;
			++it;
			if ( it != cs_shifts.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
							 "need to specify dupletss <cs_shifts> <column> with option -evaluation:chemical_shifts   last read: "+fname );
			}
			eval.add_evaluation( new ChemicalShiftEvaluator( column, fname ) );
		}
	}


  if ( option[ OptionKeys::evaluation::cam_shifts ].user() ) {
    typedef utility::vector1< std::string > CSVector;
    CSVector const& cs_shifts( option[ OptionKeys::evaluation::cam_shifts ]() );

    for ( CSVector::const_iterator it=cs_shifts.begin(); it!=cs_shifts.end(); ++it ) {
      std::string fname( *it );
      std::string column;
      ++it;
      if ( it != cs_shifts.end() ) {
        column = *it;
      } else {
        utility_exit_with_message(
               "need to specify dupletss <cs_shifts> <column> with option -evaluation:cam_shifts   last read: "+fname );
      }
      eval.add_evaluation( new CamShiftEvaluator( column, fname ) );
    }
  }

	if ( option[ OptionKeys::evaluation::pales ].user() ) {
    typedef utility::vector1< std::string > CSVector;
    CSVector const& pales( option[ OptionKeys::evaluation::pales ]() );

    for ( CSVector::const_iterator it=pales.begin(); it!=pales.end(); ++it ) {
      std::string fname( *it );
      std::string column;
      ++it;
      if ( it != pales.end() ) {
        column = *it;
      } else {
        utility_exit_with_message(
               "need to specify dupletss <pales_rdcs> <column>  with option -evaluation:pales   last read: "+fname );
      }
      eval.add_evaluation( new PalesEvaluator( column, fname ) );
    }
  }


/*  if ( option[ OptionKeys::evaluation::cam_shifts ].user() ) {
    eval.add_evaluation( new CamShiftEvaluator( "cs_score", option[ OptionKeys::evaluation::cam_shifts ]() ) );
  }
*/


	if ( option[ OptionKeys::evaluation::extra_score ].user() ) {
		using namespace core::scoring;
		utility::vector1< std::string > const& extra_scores( option[ OptionKeys::evaluation::extra_score ]() );
    utility::vector1< std::string > const& extra_score_names( option[ OptionKeys::evaluation::extra_score_column]() );
		if ( extra_scores.size() != extra_score_names.size() ) {
			utility_exit_with_message("-extra_score: you need to provide as much extra_score_names as extra_scores! ");
		}
		for ( Size ct = 1; ct <= extra_scores.size(); ct ++ ) {
			std::string const& tag = extra_score_names[ ct ];
			std::string patch( "NOPATCH" );
			if ( option[ OptionKeys::evaluation::extra_score_patch ].user() ) {
				if ( option[ OptionKeys::evaluation::extra_score_patch ]().size() != extra_scores.size() ) {
					utility_exit_with_message("-extra_score: you need to provide as much extra_score_patch(es) as \
                    extra_scores! use NOPATCH as placeholder");
				}
				patch = option[ OptionKeys::evaluation::extra_score_patch ]()[ ct ];
			}
			ScoreFunctionOP scfxn( NULL );

			if ( patch != "NOPATCH" ) {
				scfxn = ScoreFunctionFactory::create_score_function( extra_scores[ ct ], patch );
			} else {
				scfxn = ScoreFunctionFactory::create_score_function( extra_scores[ ct ] );
			}

			std::string name( extra_scores[ ct ] );
			if ( (name == "score0") ||
				(name == "score2") ||
				(name == "score3") ||
				(name == "score5") ) {
				core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scfxn );
			} else {
				core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scfxn );
			}

			std::string select_string( "SELECT_ALL" );
			if ( option[ OptionKeys::evaluation::extra_score_select ].user() ) {
				if ( option[ OptionKeys::evaluation::extra_score_select ]().size() != extra_scores.size() ) {
					utility_exit_with_message("-extra_score: you need to provide as much extra_score_patch(es) as \
                    extra_scores! use SELECT_ALL as placeholder");
				}
				select_string = option[ OptionKeys::evaluation::extra_score_select ]()[ ct ];
			}
			if ( select_string != "SELECT_ALL" ) {
				loops::Loops core;
				core.read_loop_file( select_string, false, "RIGID" );
				utility::vector1< Size> selection;
				core.get_residues( selection );
				eval.add_evaluation( new evaluation::TruncatedScoreEvaluator( tag, selection, scfxn ) );
			} else {
				eval.add_evaluation( new ScoreEvaluator( tag, scfxn ) );
			}
		}
	}


	if ( option[ OptionKeys::evaluation::jump_nr ]() ) {
		eval.add_evaluation( new JumpNrEvaluator );
	}
  //return eval;
} //read_common_options


//@detail find residues that don't have missing density
void find_existing_residues(  core::pose::PoseCOP pose, std::string tag, core::scoring::ResidueSelection& selection ) {
	for ( Size pos = 1; pos <= pose->total_residue(); pos++ ) {
		if ( pose->residue_type( pos ).is_protein() && pose->residue_type( pos ).has("CA") ) {
			numeric::xyzVector< core::Real> ca_pos = pose->residue( pos ).atom("CA").xyz();
			bool good ( true );
			for ( Size j=1; j<= pose->residue( pos ).natoms(); ++j ) {
				if ( ( ca_pos - pose->residue( pos ).atom(j).xyz() ).length() > 20 ) {
					good = false;
				}
			}
			if ( good ) selection.push_back( pos );
		}
	}
	if ( tr.Trace.visible() ) {
		tr.Trace << "selection of residues for rmsd of " << tag << std::endl;
		for ( std::list< core::Size >::const_iterator it = selection.begin(), eit = selection.end();
					it != eit; ++it ) {
			tr.Trace << " " << *it;
		}
		tr.Trace << std::endl;
	}
}

void evaluate_pose( core::pose::Pose& pose, PoseEvaluator& eval, std::ostream& os ) {
		//		ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		pss->fill_struct( pose, "eval" );
		eval.apply( pose, "eval", *pss );
		os << "\n";
		pss->print_score_header( os );
		os << "\n";
		pss->print_scores( os );
		os << std::endl;
}


void define_scorable_core_from_secondary_structure(
   core::fragment::SecondaryStructure const& ss_def,
	 protocols::loops::Loops& scored_core )
{
	using namespace core;
	using namespace basic::options;
	//	Size const max_loop_size( 3 );
	//	Size const max_short_helix( 5 );
	Size const max_loop_size( option[ OptionKeys::evaluation::score_sscore_maxloop ]() );
	Size const max_short_helix( option[ OptionKeys::evaluation::score_sscore_short_helix ]() );

	//find residues that are part of a short helix -- less than or equal to 5 residues
	utility::vector1< bool > short_helix( ss_def.total_residue(), false );

	//selection of loop definitions...
	//these loops define regions that are scored. Add all loops that are 4 residues or longer.
	//subsequently add also helices that have fewer than 6 residues if they terminated a long loop (>=4)
	loops::Loops unscored_loops;

	for ( Size pos=1; pos <= ss_def.total_residue(); pos++ ) {

		//detect loops
		if ( ss_def.loop_fraction( pos ) > 0.1 ) {
			//go to end of loop
			Size lpos = 1;
			for ( ; ( lpos+pos <= ss_def.total_residue() ) && ( ss_def.loop_fraction( pos+lpos ) > 0.1); ++lpos ) {}
			if ( lpos > max_loop_size ) { //this loop has 4 or more residues
				unscored_loops.add_loop( pos, pos+lpos-1 );
			}
			pos+=lpos-1;
		} // have found a loop

		// look for short helices and store in short_helix
		if ( ss_def.helix_fraction( pos ) > 0.1 ) {
			Size hpos = 1;
			for ( ; ( hpos+pos <= ss_def.total_residue() ) && ( ss_def.helix_fraction( pos+hpos ) > 0.1); ++hpos ) {}
			if ( hpos <= max_short_helix   ) { //this helix has 5 or fewer residues
				for ( Size ipos = 0; ipos < hpos; ++ipos ) {
					short_helix[ pos+ipos] = true;
				}
			}
		}

		//finished parsing secondary structure definition
	}

	//elongate loops if they are terminated by a short helix
	loops::Loops removed_short_helices( unscored_loops );
	for ( loops::Loops::const_iterator it=unscored_loops.begin(); it != unscored_loops.end(); ++it ) {
		Size npos( it->stop() + 1 );
		while ( short_helix[ npos ] ) {
			removed_short_helices.add_loop( npos-1, npos );
			npos++;
		}
	}

	scored_core =	removed_short_helices.invert( ss_def.total_residue() );
}


}
}
