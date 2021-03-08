// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    glycomutagenesis.cc
/// @brief   This application performs a simulated shotgun glycomutagenesis experiment of an input structure.
/// @author  Labonte <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <protocols/enzymatic_movers/GlycosyltransferaseMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/enzymes/EnzymeManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>


using namespace std;
using namespace utility;
using namespace core;
using namespace protocols;


// Class Definitions //////////////////////////////////////////////////////////

/// @brief  This protocol exhaustively mutates and glycosylates every possible position of a given protein.
class GlycomutagenesisProtocol : public moves::Mover {
public:  // Standard methods
	/// @brief  Default constructor.
	GlycomutagenesisProtocol() : moves::Mover()
	{
		init();
	}

	/// @brief  Copy constructor.
	GlycomutagenesisProtocol( GlycomutagenesisProtocol const & object_to_copy ) : Mover( object_to_copy )
	{
		copy_data( *this, object_to_copy );
	}

	// Assignment operator
	GlycomutagenesisProtocol &
	operator=( GlycomutagenesisProtocol const & object_to_copy )
	{
		// Abort self-assignment.
		if ( this != &object_to_copy ) {
			moves::Mover::operator=( object_to_copy );
			copy_data( *this, object_to_copy );
		}
		return *this;
	}

	// Destructor
	virtual ~GlycomutagenesisProtocol() {}


public:  // Standard Rosetta methods
	// General methods
	/// @brief  Register options with the option system.
	static void
	register_options()
	{
		//using namespace basic::options;
		using namespace protocols::minimization_packing;

		//option.add_relevant( OptionKeys::rings::idealize_rings );

		// Call register_options() on all other Movers used by this class.
		simple_moves::MutateResidue::register_options();
		enzymatic_movers::GlycosyltransferaseMover::register_options();
		PackRotamersMover::register_options();
		MinMover::register_options();
		simple_moves::SmallMover::register_options();
	}

	/// @brief  Generate string representation of GlycomutagenesisProtocol for debugging purposes.
	virtual
	void
	show( std::ostream & output=std::cout ) const
	{
		moves::Mover::show( output );  // name, type, tag
	}


	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual
	std::string
	get_name() const
	{
		return type();
	}

	virtual
	moves::MoverOP
	clone() const
	{
		return protocols::moves::MoverOP( utility::pointer::make_shared< GlycomutagenesisProtocol >( *this ) );
	}

	virtual
	moves::MoverOP
	fresh_instance() const
	{
		return protocols::moves::MoverOP( utility::pointer::make_shared< GlycomutagenesisProtocol >() );
	}

	/// @brief    Apply the corresponding protocol to <pose>.
	/// @details  On each apply, the Mover glycosylates the pose at the residue position corresponding to the current
	/// model as designated by the job distributor, using the consensus sequence from c. jejuni PglB.  This is
	/// currently hard-coded into the protocol, which is not yet intended for public use.
	virtual
	void
	apply( core::pose::Pose & pose )
	{
		using namespace enzymes;

		show( cout );

		Size const n_res_pre_glycosylation( pose.size() );

		//sf_->show( cout, pose );


		cout << "Mutating..." << endl;

		jd2::JobCOP job( jd2::JobDistributor::get_instance()->current_job() );
		core::uint const seqpos( job->nstruct_index() );

		cout << " Installing Glycan at Position " << seqpos << "..." << endl;

		string const & enzyme_family( glycosyltransferase_->get_enzyme_family() );

		cout << " Sequon: " <<
			EnzymeManager::get_consensus_sequence( enzyme_family, species_name_, enzyme_name_ ) << endl;

		vector1< vector1< string > > const & consensus_residues(
			EnzymeManager::get_consensus_residues( enzyme_family, species_name_, enzyme_name_ ) );
		core::uint const site_residue_position( EnzymeManager::get_reactive_residue_consensus_sequence_position(
			enzyme_family, species_name_, enzyme_name_ ) );
		Size const n_consensus_residues( consensus_residues.size() );
		Size const n_residues_left_of_site( site_residue_position - 1 );
		Size const n_residues_right_of_site( n_consensus_residues - site_residue_position );

		if ( ( seqpos < site_residue_position ) || ( seqpos + n_residues_right_of_site > n_res_pre_glycosylation ) ) {
			cout << "  Sequon will not fit here." << endl;
			minimizer_->score_function( sf_ );
			minimizer_->max_iter( 1 );
			minimizer_->apply ( pose );
			return;
		}

		disulfide_treatment( pose, seqpos, n_residues_left_of_site, n_residues_right_of_site );

		for ( core::uint i( seqpos - n_residues_left_of_site ); i <= seqpos + n_residues_right_of_site; ++ i ) {
			cout << "  Mutating position " << i;
			mutator_->set_target( i );

			string const & new_res( consensus_residues[ i - seqpos + site_residue_position ][ 1 ] );
			cout << " to " << new_res << endl;
			mutator_->set_res_name( chemical::aa_from_name( new_res ) );
			mutator_->apply( pose );
		}

		cout << "  Packing & minimizing..." << endl;

		packer_->apply( pose );
		minimizer_->apply ( pose );


		cout << "Glycosylating..." << endl;

		glycosyltransferase_->apply( pose );

		//sf_->show( cout, pose );


		cout << " Minimizing..." << endl;

		Size const n_res_post_glycosylation( pose.size() );

		mm_->set_bb_true_range( n_res_pre_glycosylation + 1, n_res_post_glycosylation );
		mm_->set_chi( true );
		mm_->set_branches( true );

		minimizer_->movemap( mm_ );
		minimizer_->score_function( sf_sugar_only_ );
		minimizer_->apply ( pose );
		minimizer_->score_function( sf_ );
		minimizer_->apply ( pose );

		//sf_->show( cout, pose );


		cout << " Packing..." << endl;

		packer_->apply( pose );

		minimizer_->apply ( pose );

		//sf_->show( cout, pose );


		cout << "Refining..." << endl;

		mc_ = moves::MonteCarloOP( utility::pointer::make_shared< moves::MonteCarlo >( pose, *sf_, kt_ ) );

		for ( core::uint cycle( 1 ); cycle <= n_cycles_; ++cycle ) {
			small_mover_->apply( pose );
			packer_->apply( pose );
			minimizer_->apply ( pose );

			// Metropolis criterion.
			mc_->boltzmann( pose );
			if ( cycle % 5 == 0 ) {
				cout << " Cycle " << cycle << "\tCurrent Score:" << ( *sf_ )( pose ) << endl;
			}
		}

		cout << " Final score for decoy:" << endl;
		sf_->show( cout, pose );
	}


private:  // Private methods
	// Set command-line options.  (Called by init())
	void
	set_commandline_options()
	{
		//using namespace basic::options;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( option[ OptionKeys::run::n_cycles ] > 1 ) {
			n_cycles_ = option[ OptionKeys::run::n_cycles ];
		}
	}


	// Initialize data members from arguments.
	void
	init()
	{
		using namespace scoring;
		using namespace pack::task;
		using namespace simple_moves;
		using namespace minimization_packing;

		type( "GlycomutagenesisProtocol" );

		sf_ = get_score_function();
		sf_->set_weight( fa_intra_rep_nonprotein, 0.550 );
		sf_sugar_only_ = ScoreFunctionOP( utility::pointer::make_shared< ScoreFunction >() );
		sf_sugar_only_->set_weight( sugar_bb, 1.0 );

		mm_ = kinematics::MoveMapOP( utility::pointer::make_shared< kinematics::MoveMap >() );

		// Instantiate the Movers.
		mutator_ = MutateResidueOP( utility::pointer::make_shared< MutateResidue >() );

		glycosyltransferase_ =
			enzymatic_movers::GlycosyltransferaseMoverOP( utility::pointer::make_shared< enzymatic_movers::GlycosyltransferaseMover >() );
		glycosyltransferase_->set_species( species_name_ );
		glycosyltransferase_->set_enzyme( enzyme_name_ );
		glycosyltransferase_->set_efficiency( 1.0 );

		small_mover_ = SmallMoverOP( utility::pointer::make_shared< SmallMover >( mm_, kt_, 1 ) );

		minimizer_ = MinMoverOP( utility::pointer::make_shared< MinMover >() );

		TaskFactoryOP tf( utility::pointer::make_shared< TaskFactory >() );
		tf->push_back( operation::RestrictToRepackingOP( utility::pointer::make_shared< operation::RestrictToRepacking >() ) );
		tf->push_back( operation::IncludeCurrentOP( utility::pointer::make_shared< operation::IncludeCurrent >() ) );
		packer_ = PackRotamersMoverOP( utility::pointer::make_shared< PackRotamersMover >( sf_ ) );
		packer_->task_factory( tf );

		set_commandline_options();
	}

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void
	copy_data( GlycomutagenesisProtocol & object_to_copy_to, GlycomutagenesisProtocol const & object_to_copy_from )
	{
		//object_to_copy_to.species_name_ = object_to_copy_from.species_name_;
		//object_to_copy_to.enzyme_name_ = object_to_copy_from.enzyme_name_;

		object_to_copy_to.sf_ = object_to_copy_from.sf_;
		object_to_copy_to.sf_sugar_only_ = object_to_copy_from.sf_sugar_only_;
		object_to_copy_to.mm_ = object_to_copy_from.mm_;

		object_to_copy_to.mutator_ = object_to_copy_from.mutator_;
		object_to_copy_to.glycosyltransferase_ = object_to_copy_from.glycosyltransferase_;
		object_to_copy_to.small_mover_ = object_to_copy_from.small_mover_;
		object_to_copy_to.minimizer_ = object_to_copy_from.minimizer_;
		object_to_copy_to.packer_ = object_to_copy_from.packer_;

		object_to_copy_to.mc_ = object_to_copy_from.mc_;
	}


	// Breaking disulfide to avoid runtime errors. 
	// name is "disulfide_treatment: because other alternate option "ignore" can also be used. 
	void 
	disulfide_treatment( core::pose::Pose & pose, core::uint const seqpos, Size const n_residues_left_of_site, Size const n_residues_right_of_site )
	{
		for ( core::uint i( seqpos - n_residues_left_of_site ); i <= seqpos + n_residues_right_of_site; ++i ) {

			if ( pose.residue_type( i ).is_disulfide_bonded() ) {
				Size res1_disulf_atom_num = pose.residue_type( i ).atom_index( pose.residue_type( i ).get_disulfide_atom_name() );
				Size res1_disulf_partner_id = pose.residue_type( i ).residue_connection_id_for_atom( res1_disulf_atom_num );
				Size res2_disulf_partner_num = pose.residue( i ).residue_connection_partner( res1_disulf_partner_id );
				cout << "  Disulfide found! Breaking cysteine dimer: CYS[" << i << "]-CYS[" << res2_disulf_partner_num << "]" << endl;
				change_cys_state( i, "CYS", pose.conformation() );
				change_cys_state( res2_disulf_partner_num, "CYS", pose.conformation() );
			}
		}
	}


private:  // Private data
	string const species_name_ = "c_jejuni";
	string const enzyme_name_ = "PglB";

	scoring::ScoreFunctionOP sf_;
	scoring::ScoreFunctionOP sf_sugar_only_;
	kinematics::MoveMapOP mm_;

	Real const kt_ = 0.8;
	Size n_cycles_ = 100;

	// Movers
	simple_moves::MutateResidueOP mutator_;
	enzymatic_movers::GlycosyltransferaseMoverOP glycosyltransferase_;
	minimization_packing::MinMoverOP minimizer_;
	minimization_packing::PackRotamersMoverOP packer_;
	simple_moves::SmallMoverOP small_mover_;

	protocols::moves::MonteCarloOP mc_;
};  // class


// Constants & Type Definitions ///////////////////////////////////////////////
int const SUCCESS( 0 );
int const FAILURE( -1 );

typedef utility::pointer::shared_ptr< GlycomutagenesisProtocol > GlycomutagenesisProtocolOP;


// Main ///////////////////////////////////////////////////////////////////////
int
main( int argc, char *argv[] )
{
	try {
		// Initialize Rosetta.
		cout << "Initializing Rosetta..." << endl;
		devel::init( argc, argv );

		// Construct the protocol.
		GlycomutagenesisProtocolOP protocol( utility::pointer::make_shared< GlycomutagenesisProtocol >() );

		// Distribute the mover.
		jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return FAILURE;
	}
	return SUCCESS;
}
