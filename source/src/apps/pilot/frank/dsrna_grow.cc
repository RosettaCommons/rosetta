/// @file
/// @brief


#include <devel/init.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
//#include <protocols/rna/movers/RNA_HelixAssembler.hh>
#include <core/import_pose/RNA_HelixAssembler.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>

#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

using namespace basic;
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::chemical;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;



OPT_KEY( Boolean,  init )
OPT_KEY( Boolean,  premin )
OPT_KEY( Integer,  blocksize )
OPT_KEY( Integer,  melt )
OPT_KEY( Integer,  growDS )
OPT_KEY( Integer,  growUS )
OPT_KEY( Integer,  minevery )

static basic::Tracer TR("dsRNA_grow");

class dsRNA_grow : public protocols::moves::Mover {
private:
	core::chemical::rna::RNA_FittedTorsionInfo const torsion_info_;
	core::chemical::ResidueTypeSetCOP rsd_set_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::optimization::AtomTreeMinimizerOP minimizer_;
	core::optimization::MinimizerOptionsOP minimizer_options_;

public:
	dsRNA_grow() {
		using namespace core::optimization;
		minimizer_ = core::optimization::AtomTreeMinimizerOP( new AtomTreeMinimizer );
		minimizer_options_ = core::optimization::MinimizerOptionsOP( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.00001, false, false, false ) );
		minimizer_options_->nblist_auto_update( true );
	}

	virtual std::string get_name() const { return "dsRNA_grow"; }

	// stolen from Rhiju's RNA helix assembler
	void
	constrain_and_minimize( pose::Pose & pose, utility::vector1< std::pair< Size, Size > > const & pairings ) const {
		using namespace core::scoring;
		using namespace core::optimization;

		scoring::constraints::ConstraintSetOP new_cst_set;
		pose.constraint_set( new_cst_set ); //blank out cst set.
		core::pose::rna::setup_base_pair_constraints( pose, pairings );

		kinematics::MoveMap mm;
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( true );
		for ( Size i = 1; i <= pairings.size(); i++ ) {
			mm.set_bb( pairings[i].first, true );
			mm.set_chi( pairings[i].first, true );
			mm.set_bb( pairings[i].second, true );
			mm.set_chi( pairings[i].second, true );
		}

		Real const cst_weight = scorefxn_->get_weight( base_pair_constraint );
		//Real const dens_weight = scorefxn_->get_weight( elec_dens_fast );

		//scorefxn_->set_weight( elec_dens_fast, 0.0 );
		minimizer_->run( pose, mm, *scorefxn_, *minimizer_options_ );

		(*scorefxn_)( pose );
		//scorefxn_->set_weight( elec_dens_fast, dens_weight );
		scorefxn_->set_weight( base_pair_constraint, 0.0 );
		minimizer_->run( pose, mm, *scorefxn_, *minimizer_options_ );

		scorefxn_->set_weight( base_pair_constraint, cst_weight );
		(*scorefxn_)( pose );
	}

	void
	minimize( pose::Pose & pose ) const {
		using namespace core::scoring;
		using namespace core::optimization;

		scoring::constraints::ConstraintSetOP new_cst_set;
		pose.constraint_set( new_cst_set ); //blank out cst set.

		kinematics::MoveMap mm;
		mm.set_bb( true );
		mm.set_chi( true );
		mm.set_jump( true );
		minimizer_->run( pose, mm, *scorefxn_, *minimizer_options_ );
		(*scorefxn_)( pose );
	}

	// stolen from Rhiju's RNA helix assembler
	void
	set_Aform_torsions( pose::Pose & pose, Size const & n ) const
	{
		using namespace core::id;
		using namespace core::pose::rna;
		using namespace core::chemical::rna;
		pose.set_torsion( TorsionID( n, BB, 1), torsion_info_.alpha_aform() );
		pose.set_torsion( TorsionID( n, BB, 2), torsion_info_.beta_aform() );
		pose.set_torsion( TorsionID( n, BB, 3), torsion_info_.gamma_aform() );
		pose.set_torsion( TorsionID( n, BB, 4), torsion_info_.delta_north() );
		pose.set_torsion( TorsionID( n, BB, 5), torsion_info_.epsilon_aform() );
		pose.set_torsion( TorsionID( n, BB, 6), torsion_info_.zeta_aform() );

		apply_pucker(pose, n, NORTH, false /*skip_same_state*/, false);
	}

	// stolen from Rhiju's RNA helix assembler
	void
	append_Aform_residue( pose::Pose & pose, Size const & n, std::string const & nt ) const {
		using namespace core::conformation;
		using namespace core::chemical;
		using namespace core::id;

		runtime_assert( pose.fold_tree().is_cutpoint(n) || n == pose.size() );

		ResidueOP rsd1 = core::conformation::ResidueFactory::create_residue(
			*core::pose::residue_types_from_sequence( nt, *rsd_set_, false /*auto_termini*/ )[1] );
		runtime_assert( rsd1->is_NA() );
		pose.conformation().safely_append_polymer_residue_after_seqpos( *rsd1, n, true /*build_ideal_geometry*/ );
		pose.set_torsion( TorsionID( n, BB, 5), torsion_info_.epsilon_aform());
		pose.set_torsion( TorsionID( n, BB, 6), torsion_info_.zeta_aform());
		set_Aform_torsions( pose, n+1 );
	}


	// stolen from Rhiju's RNA helix assembler
	void
	prepend_Aform_residue( pose::Pose & pose, Size const & n, std::string const & nt ) const {
		using namespace core::conformation;
		using namespace core::chemical;
		using namespace core::id;

		runtime_assert( n == 1 || pose.fold_tree().is_cutpoint(n-1) );

		ResidueOP rsd2 = core::conformation::ResidueFactory::create_residue(
			*core::pose::residue_types_from_sequence( nt, *rsd_set_, false /*auto_termini*/ )[1] );
		runtime_assert( rsd2->is_NA() );
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *rsd2, n, true /*build_ideal_geometry*/ );

		set_Aform_torsions( pose, n );
		pose.set_torsion( TorsionID( n+1, BB, 1), torsion_info_.alpha_aform());
		pose.set_torsion( TorsionID( n+1, BB, 2), torsion_info_.beta_aform());
		pose.set_torsion( TorsionID( n+1, BB, 3), torsion_info_.gamma_aform());
	}


	void apply( core::pose::Pose & pose) {
		rsd_set_ = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		scorefxn_ = core::scoring::get_score_function();
		if ( !scorefxn_->has_nonzero_weight( core::scoring::base_pair_constraint ) ) {
			TR << TR.Magenta << "Scorefunction does not have base_pair_constraint on -- setting weight to 5.0." << TR.Reset << std::endl;
			scorefxn_->set_weight( core::scoring::base_pair_constraint, 5.0 );
		}

		core::Size blocksize_ = option[ blocksize ];
		core::Size melt_ = option[ melt ];
		core::Size growDS_ = option[ growDS ]() / blocksize_;
		core::Size growUS_ = option[ growUS ]() / blocksize_;
		core::Size minevery_ = option[ minevery ];

		if ( option[init] ) {
			std::string sequence_to_build;
			for ( core::Size i=1; i<blocksize_/2; ++i ) {
				sequence_to_build += "au";
			}
			pose::make_pose_from_sequence( pose, sequence_to_build, *rsd_set_ );

			core::import_pose::RNA_HelixAssembler rna_helix_assembler;
			rna_helix_assembler.set_minimize_all( false );
			rna_helix_assembler.set_dump( false );
			rna_helix_assembler.use_phenix_geo( true );

			rna_helix_assembler.set_scorefxn ( scorefxn_ );
			rna_helix_assembler.set_model_and_remove_capping_residues( true );

			rna_helix_assembler.apply( pose, sequence_to_build );
		} else {
			// now add density scores from cmd line
			core::Size sub=0;
			core::Size N = pose.total_residue();

			if ( option[ edensity::mapfile ].user() ) {
				core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_ );
				protocols::electron_density::SetupForDensityScoringMoverOP edens( new protocols::electron_density::SetupForDensityScoringMover );
				edens->apply( pose );
				sub=1;
			}

			core::kinematics::FoldTree f = pose.fold_tree();
			if ( option[ edensity::mapfile ].user() ) {
				ObjexxFCL::FArray2D<core::Size> Fjumps(2, 2);
				ObjexxFCL::FArray1D<core::Size> Fcuts(2);
				Fjumps(1, 1) = N/4;
				Fjumps(2, 1) = N+1;
				Fjumps(1, 2) = N/4;
				Fjumps(2, 2) = 3*N/4;
				Fcuts(1) = N;
				Fcuts(2) = N/2;
				f.tree_from_jumps_and_cuts(N+1, 2, Fjumps, Fcuts, N+1, true);
			} else {
				f.slide_jump( 1, N/4, 3*N/4);
			}
			pose.fold_tree(f);

			if ( option[premin] ) {
				minimize( pose );
			}

			for ( core::Size i=1; i<=growDS_; ++i ) {
				for ( core::Size j=1; j<=blocksize_; ++j ) {
					N = pose.total_residue() - sub;
					append_Aform_residue( pose, N, j%2==1?"a":"u" );
					prepend_Aform_residue( pose, 1, j%2==1?"u":"a" );
				}
				N = pose.total_residue() - sub;
				utility::vector1< std::pair< Size, Size > > pairings;
				for ( core::Size j=1; j<=blocksize_+melt_; ++j ) {
					pairings.push_back( std::make_pair( j, N-j+1) );
				}
				constrain_and_minimize( pose, pairings );
				if ( minevery_>0 && i%minevery_ == 0 ) {
					minimize( pose );
				}
				if ( i%10 == 0 ) {
					pose.dump_pdb( "termgrow"+utility::to_string(i)+".pdb");
				}
			}
			for ( core::Size i=1; i<=growUS_; ++i ) {
				for ( core::Size j=1; j<=blocksize_; ++j ) {
					N = (pose.total_residue() - sub) / 2;
					append_Aform_residue( pose, N, j%2==1?"a":"u" );
					prepend_Aform_residue( pose, N+2, j%2==1?"u":"a" );
				}

				N = (pose.total_residue() - sub) / 2;
				utility::vector1< std::pair< Size, Size > > pairings;
				for ( core::Size j=1; j<=blocksize_+melt_; ++j ) {
					pairings.push_back( std::make_pair( N-j+1, N+j) );
				}
				constrain_and_minimize( pose, pairings );

				if ( minevery_>0 && i%minevery_ == 0 ) {
					minimize( pose );
				}
				if ( i%10 == 0 ) {
					pose.dump_pdb( "midgrow"+utility::to_string(i)+".pdb");
				}
			}

			minimize( pose );
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	using namespace protocols::moves;
	using namespace protocols;
	using namespace protocols::jd2;

	try {
		NEW_OPT( init, "Initialize?", false );
		NEW_OPT( premin, "Preminimize?", false );
		NEW_OPT( blocksize , "Bps to add in each grow step", 8 );
		NEW_OPT( melt, "Number of prior bps to minimize", 4 );
		NEW_OPT( minevery, "Full min every n cycles", 0 );
		NEW_OPT( growDS, "Number of bps to grow downstream", 60 );
		NEW_OPT( growUS, "Number of bps to grow upstream", 60 );

		devel::init( argc, argv );
		SequenceMoverOP seq( new SequenceMover() );
		seq->add_mover( MoverOP(new dsRNA_grow()) );
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}

