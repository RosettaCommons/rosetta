// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseScreener
/// @brief Makes a list of (phi, psi, omega) at moving_residues that
///              could be useful for full-atom packing
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseScreener.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/stepwise/MainChainTorsionSet.hh>

//////////////////////////////////
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/id/AtomID.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/pose/Pose.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/Jump.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/exit.hh>

#include <string>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>

//Auto Headers



using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.StepWiseScreener" ) ;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {



  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseScreener::StepWiseScreener(
																		 utility::vector1< Size > const & moving_residues
																		 ):
		moving_residues_( moving_residues ),
		n_sample_( 18 /* Corresponds to 20 degree bins */ ),
		rmsd_cutoff_( -1.0 ),
		centroid_scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "score3.wts"  )/* score3 */ ),
		silent_file_( "" ),
		filter_native_big_bins_( false ),
		centroid_screen_( false ),
		centroid_score_ref_( 999999999999999.99),
		centroid_score_diff_cut_( 20.0 ),
		apply_vdw_cut_( false ),
		centroid_vdw_ref_( 9999999999999999.999 ),
		nstruct_centroid_( 0 ),
		ghost_loops_( false )
  {
		for ( Size n = 1; n <= moving_residues.size(); n++ ) 	{
			main_chain_torsion_set_for_moving_residues_.push_back( MainChainTorsionSet( 0.0, 0.0, 0.0 ) );
		}

  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseScreener::~StepWiseScreener()
  {}

  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseScreener::apply( core::pose::Pose & pose )
	{

		Size which_res( 1 );
		Size count( 1 );

		clock_t const time_start( clock() );

		main_chain_torsion_sets_for_moving_residues_.clear();
		centroid_scores_.clear();

		// convert to centroid
		pose::Pose pose_save = pose;
		if ( centroid_screen_ ) convert_to_centroid( pose );
		if ( ghost_loops_ ) prepare_ghost_pose( pose ); // pose with loops excised

		sample_residues_recursively( which_res, count, pose );

		std::cout << "Total time in StepWiseScreener: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

		pose = pose_save;

		if ( nstruct_centroid_ > 0 && centroid_screen_ ) filter_main_chain_torsion_sets();

	}

	std::string
	StepWiseScreener::get_name() const {
		return "StepWiseScreener";
	}

	///////////////////////////////////////////////////////////////////////////
	// pose with loops excised. centroid level.
	void
	StepWiseScreener::prepare_ghost_pose( core::pose::Pose const & pose ){

		using namespace core::kinematics;
		using namespace core::scoring;

		// Figure out sequence of a pose without the loops (i.e., the "moving residues" );
		std::string const full_sequence = pose.sequence();
		ObjexxFCL::FArray1D<bool> moving_array( full_sequence.size(), false );
		for ( Size n = 1; n <= moving_residues_.size(); n++ ) moving_array( moving_residues_[n] ) = true;

		ghost_map_.clear();
		Size count( 0 );
		std::string desired_sequence  = "";
		for ( Size n = 1; n <= full_sequence.size(); n++ ){
			if ( moving_array(n) ) continue;
			count++;
			desired_sequence += full_sequence[n-1];
			ghost_map_[ count ] = n;
		}

		FoldTree f = figure_out_fold_tree( ghost_map_ ); //chainbreaks, etc.

		std::cout << "FULL SEQUENCE " << full_sequence << std::endl;
		std::cout << "DESIRED SEQUENCE " << desired_sequence << std::endl;

		initialize_ghost_pose( ghost_pose_, desired_sequence, pose, ghost_map_, f);
		ghost_pose_->dump_pdb( "GHOST_START.pdb" );

		if ( get_native_pose() ) {
			Pose native_centroid_pose =  *get_native_pose();
			convert_to_centroid( native_centroid_pose );
			initialize_ghost_pose( ghost_native_pose_, desired_sequence, native_centroid_pose, ghost_map_, f);
			ghost_native_pose_->dump_pdb( "GHOST_NATIVE_START.pdb" );
		}

		// Also need to figure out reference energy.
		Pose ghost_pose_blowup = *ghost_pose_;
		for ( Size n = 1; n <= ghost_pose_blowup.num_jump(); n++ ) {
			Jump jump( ghost_pose_blowup.jump( n ) );
			jump.set_translation( Vector( 100.0, 0.0, 0.0 ) ); //This is a little dangerous.
			ghost_pose_blowup.set_jump( n, jump );
		}
		ghost_pose_blowup.dump_pdb( "GHOST_BLOWUP.pdb" );
		centroid_score_ref_ = (*centroid_scorefxn_)( ghost_pose_blowup );
		centroid_vdw_ref_ = ghost_pose_blowup.energies().total_energies()[ vdw ];

		centroid_scorefxn_->show( std::cout, ghost_pose_blowup );
		std::cout << " REFERENCE SCORES " << centroid_score_ref_ << " " << centroid_vdw_ref_ << std::endl;

	}

	///////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::initialize_ghost_pose(
		  core::pose::PoseOP & ghost_pose,
			std::string const & desired_sequence,
			core::pose::Pose const & template_pose,
			ResMap const & ghost_map,
			core::kinematics::FoldTree f )
	{
		using namespace core::chemical;
		using namespace core::pose;
		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( CENTROID );

		ghost_pose = new Pose;
		core::pose::make_pose_from_sequence( *ghost_pose, desired_sequence, *rsd_set );
		copy_coords( *ghost_pose, template_pose, ghost_map );
		ghost_pose->fold_tree( f );
	}

	///////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::copy_coords( pose::Pose & pose, pose::Pose const & template_pose, ResMap const & ghost_map ) const
	{
		using namespace core::id;

		template_pose.residue( 1 ); // force a refold.

		for ( ResMap::const_iterator it=ghost_map.begin(); it != ghost_map.end(); it++ ) {
			Size const & res( it->first );
			Size const & template_res( it->second );

			for( Size j = 1; j <= pose.residue_type( res ).natoms(); j++ ) {

				if ( pose.residue_type( res ).atom_name( j ) !=
						 template_pose.residue_type( template_res ).atom_name( j ) ) {
					std::cout << "PROBLEM! " << res << " " << pose.residue( res ).atom_name( j ) <<  "  !=  " <<
						template_res << " " << template_pose.residue( template_res ).atom_name( j ) <<  std::endl;
					utility_exit_with_message( "mismatch in ghost pose" );
				}

				pose.set_xyz( AtomID( j, res ), template_pose.xyz( AtomID( j, template_res ) ) );

			}

		}

		pose.residue( 1 ); // force a refold.
	}

	///////////////////////////////////////////////////////////////////////////
	core::kinematics::FoldTree
	StepWiseScreener::figure_out_fold_tree( ResMap const & ghost_map ) const {

		Size prev_res( 0 ), total_res( 0 );
		utility::vector1< Size > cutpoints;

		for ( ResMap::const_iterator it=ghost_map.begin(); it != ghost_map.end(); it++ ) {
			Size const & res( it->first );
			Size const & template_res( it->second );

			std::cout << "MAPPING " << res << " --> " << template_res << std::endl;

			if ( (template_res-1) != prev_res  ) cutpoints.push_back( res-1 );
			prev_res = template_res;

			total_res = res;
		}

		core::kinematics::FoldTree f( total_res );

		for ( Size n = 1; n <= cutpoints.size(); n++ ) {
			std::cout << "Adding jump across: " << cutpoints[ n ] << std::endl;
			f.new_jump( cutpoints[ n ], cutpoints[ n ] +1, cutpoints[ n ] );
		}

		return f;

	}


	///////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::convert_to_centroid( core::pose::Pose & pose ) {

		core::pose::remove_variant_type_from_pose_residue( pose, "N_ACETYLATION", 1 );
		core::pose::add_variant_type_to_pose_residue( pose, "LOWER_TERMINUS", 1 );

		core::pose::remove_variant_type_from_pose_residue( pose, "C_METHYLAMIDATION", pose.total_residue() );
		core::pose::add_variant_type_to_pose_residue( pose, "UPPER_TERMINUS", pose.total_residue() );

		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

		// get DSSP, assign secondary structure
		core::scoring::dssp::Dssp dssp_obj( pose );
		dssp_obj.insert_ss_into_pose( pose );
		//pose.dump_pdb( "CENTROID.pdb" );

	}

	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::sample_residues_recursively(
								Size const which_res,
								Size & count,
								core::pose::Pose & pose )
	{

		 using namespace core::chemical;
		 using namespace core::scoring;

		 Size const n = moving_residues_[ which_res ];

		 MainChainTorsionSetList main_chain_torsion_set_list;
		 get_main_chain_torsion_set_list( n, pose, main_chain_torsion_set_list );

		 for ( Size k = 1; k <= main_chain_torsion_set_list.size(); k++ ) {

			 MainChainTorsionSet const & main_chain_torsion_set( main_chain_torsion_set_list[ k ] );

			 pose.set_phi( n, main_chain_torsion_set.phi() );
			 pose.set_psi( n, main_chain_torsion_set.psi() );
			 pose.set_omega( n, main_chain_torsion_set.omega() );

			 main_chain_torsion_set_for_moving_residues_[ which_res ] = main_chain_torsion_set;

			 if ( which_res == moving_residues_.size() ) {

				 count++;
				 std::string const tag = "S_"+ ObjexxFCL::lead_zero_string_of( count, 5 );
				 filter_and_output( pose, tag );

			 } else {
				 sample_residues_recursively( which_res+1, count, pose );
			 }

		 }

	 }

	 ///////////////////////////////////////////////////////////////////////////////
	 void
	 StepWiseScreener::filter_and_output( core::pose::Pose & pose,
																				std::string const & tag	 )
	 {
		 using namespace core::scoring;
		 using namespace core::chemical;
		 using namespace core::pose;

		 if ( get_native_pose() && rmsd_cutoff_ > 0.0 ) {
			 Real const rmsd = rmsd_with_super( pose, *get_native_pose(), is_protein_backbone_including_O );
			 TR << "CHECKING: decoy " << tag << ": " <<  rmsd << " vs. filter " << rmsd_cutoff_ << std::endl;
			 if ( rmsd > rmsd_cutoff_ ) return;
		 }

		 Real centroid_score( 0.0 );
		 if ( centroid_screen_ ) {

			 if ( ghost_loops_ ) {
				 //this may be a little inefficient, since all internal dofs needs to be recalculated. Its safe, though, I think.
				 copy_coords( *ghost_pose_, pose, ghost_map_ );

				 centroid_score = (*centroid_scorefxn_)( *ghost_pose_ ); // pose with loops excised

				 Real const centroid_vdw = ghost_pose_->energies().total_energies()[ vdw ];
				 if ( apply_vdw_cut_ && centroid_vdw > centroid_vdw_ref_ ) return;

			 } else {

				 centroid_score = (*centroid_scorefxn_)( pose );
				 // Keep running tabs on best score seen so far.
				 if ( centroid_score < centroid_score_ref_ ) centroid_score_ref_ = centroid_score;
			 }

			 Real const centroid_score_diff = centroid_score - centroid_score_ref_;

			 if (  centroid_score_diff >  centroid_score_diff_cut_ )  return;

			 //std::cout << "COMPARING " << centroid_score << " " << centroid_score_ref_ << std::endl;

		 }

		 main_chain_torsion_sets_for_moving_residues_.push_back( main_chain_torsion_set_for_moving_residues_ );
		 centroid_scores_.push_back( centroid_score );

		 if ( silent_file_.size() > 0 ) {

			 if ( ghost_loops_ ) {
				 output_silent_struct( *ghost_pose_, ghost_native_pose_, silent_file_, tag );
			 } else {
				 output_silent_struct( pose, get_native_pose(), silent_file_, tag );
			 }

		 }


 }


 //////////////////////////////////////////////////////////////////////
 // This could even be its own class...
	 void
	 StepWiseScreener::get_main_chain_torsion_set_list(
					 Size const & n,
					 pose::Pose const & pose,
					 MainChainTorsionSetList & main_chain_torsion_set_list )
	 {
		 using namespace core::chemical;

		 main_chain_torsion_set_list.clear();

		 // following is totally hacky... would be much better to
		 // figure out which bins would sum to give probability > 99%, or something like that.
		 Real best_energy_cutoff = 0.8;
		 // Not really necessary...
		 //		 if ( pose.aa( n ) == aa_gly ) best_energy_cutoff = 0.8;
		 //		if ( pose.aa( n ) == aa_pro ) best_energy_cutoff = 0.0;

		 // A few special cases first.
		 if ( n == 1 && !pose.residue( 1 ).has_variant_type( "N_ACETYLATION" ) ){

			 // If we're at the N-terminal, there's no point in sampling phi -- just a dinky hydrogen.
			 // basically slave the phi to the psi.
			 get_main_chain_torsion_set_list_n_terminus( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

		 } else if ( n == pose.total_residue() && !pose.residue( n ).has_variant_type( "C_METHYLAMIDATION" ) ){

			 // If we're at the C-terminal, there's no point in sampling psi -- just an oxygen.
			 // basically slave the psi to the phi.
			 get_main_chain_torsion_set_list_c_terminus( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

		 } else if ( n > moving_residues_[ 1 ] &&
								 n < moving_residues_[ moving_residues_.size() ] ) {

			 // Trying coarse sample for internal residues -- otherwise number of conformations really blows up.

			 get_main_chain_torsion_set_list_coarse( n, pose, main_chain_torsion_set_list );

			 //get_main_chain_torsion_set_list_full( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

		 } else {

			 /////////////////
			 // General case.
			 /////////////////
			 get_main_chain_torsion_set_list_full( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

		 }

		 // CHEAT -- later can use secondary structure.
		 if ( filter_native_big_bins_ )	 filter_native_BIG_BINS( n, main_chain_torsion_set_list );

	 }

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::get_main_chain_torsion_set_list_coarse(
																						Size const & n,
																						pose::Pose const & pose,
																						MainChainTorsionSetList & main_chain_torsion_set_list )
	 {
		 using namespace core::chemical;

		 // Later, can individually dial in cluster centers based on careful
		 // inspection of residue-specific rama plots. For now, this is basically a quick hack.
		 if ( pose.aa( n ) == aa_pro ) {

			 main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, -30.0 ) );
			 main_chain_torsion_set_list.push_back( MainChainTorsionSet(  60.0, 140.0 ) );

		 } else if ( pose.aa( n ) == aa_gly ) {

			 main_chain_torsion_set_list.push_back( MainChainTorsionSet( -80.0, -10.0 ) );
			 main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, 160.0 ) );
			 main_chain_torsion_set_list.push_back( MainChainTorsionSet(  90.0,  10.0 ) );
			 main_chain_torsion_set_list.push_back( MainChainTorsionSet( 100.0, 180.0 ) );

		 } else {

			 main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, -30.0 ) );
			 main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, 140.0 ) );
			 main_chain_torsion_set_list.push_back( MainChainTorsionSet(  50.0,  50.0 ) );

		 }

	 }


	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::get_main_chain_torsion_set_list_full( core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
																											MainChainTorsionSetList & main_chain_torsion_set_list )

	{
		//generic -- a residue in the middle of the loop
		for (Size i = 1; i <= n_sample_; i++ ) {
			for (Size j = 1; j <= n_sample_; j++ ) {
				Real const phi = get_rotamer_angle( i, n_sample_ );
				Real const psi = get_rotamer_angle( j, n_sample_ );

				if ( ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  )  > best_energy_cutoff ) {
					continue;
				}

				main_chain_torsion_set_list.push_back( MainChainTorsionSet( phi, psi ) );

			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::get_main_chain_torsion_set_list_n_terminus( core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
																														MainChainTorsionSetList & main_chain_torsion_set_list )
	{
		assert( n == 1);
		for (Size j = 1; j <= n_sample_; j++ ) {
			Real const psi = get_rotamer_angle( j, n_sample_ );
			Real best_phi = 0.0;
			Real best_energy = 999999.9999;
			for (Size i = 1; i <= n_sample_; i++ ) {
				Real const phi = get_rotamer_angle( i, n_sample_ );
				Real temp_rama = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
				if ( temp_rama  < best_energy ) {
					best_energy = temp_rama;
					best_phi = phi;
				}
			}
			if ( best_energy < best_energy_cutoff ){
				main_chain_torsion_set_list.push_back( MainChainTorsionSet( best_phi, psi ) );
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::get_main_chain_torsion_set_list_c_terminus( core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
																														MainChainTorsionSetList & main_chain_torsion_set_list )
	{
		for (Size i = 1; i <= n_sample_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_ );

			Real best_psi = 0.0;
			Real best_energy = 999999.9999;
			for (Size j = 1; j <= n_sample_; j++ ) {
				Real const psi = get_rotamer_angle( j, n_sample_ );
				Real temp_rama = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
				if ( temp_rama  < best_energy ) {
					best_energy = temp_rama;
					best_psi = psi;
				}
			}
			if ( best_energy < best_energy_cutoff ){
				main_chain_torsion_set_list.push_back( MainChainTorsionSet( phi, best_psi ) );
			}
		}
	}


	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::get_main_chain_torsion_set_list_sample_phi_only( core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
																																 MainChainTorsionSetList & main_chain_torsion_set_list )
	{
		TR << "JUNCTION RESIDUE --> PREPEND " << n << std::endl;

		// we are prepending and this is the junction residue. sample phi only!
		Real const psi = pose.psi( n );
		for (Size i = 1; i <= n_sample_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_ );
			if ( ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  )  > best_energy_cutoff ) {
				continue;
			}
			main_chain_torsion_set_list.push_back( MainChainTorsionSet( phi, psi ) );
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::get_main_chain_torsion_set_list_sample_psi_only( core::Size const & n, core::pose::Pose const & pose, core::Real const & best_energy_cutoff,
																																 MainChainTorsionSetList & main_chain_torsion_set_list )
	{
		TR << "JUNCTION RESIDUE --> APPEND " << n << std::endl;

		// we are appending and this is the junction residue. sample psi only!
		Real const phi = pose.phi( n );
		for (Size j = 1; j <= n_sample_; j++ ) {
			Real const psi = get_rotamer_angle( j, n_sample_ );
			if ( ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  )  > best_energy_cutoff ) {
				continue;
			}
			main_chain_torsion_set_list.push_back( MainChainTorsionSet( phi, psi ) );
		}

	}

	//////////////////////////////////////////////////////////////////////////////
	// Coarse -- see set_ss_from_phipsi in core/pose/util.cc.
	//////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseScreener::get_big_bin( Real const phi, Real const psi ) const
	{
		if ( phi < -20.0 && psi > -90.0 && psi < -10.0 ) {
			return 1;
		} else if ( phi < -20.0 && (psi > 20.0 || psi < -170.0) ) {
			return 2;
		}
		return 3;
	}


	///////////////////////////////////////////////////////////////
	// Check set_ss_from_phipsi -- pretty coarse. Could
	// later use some guess for strand/helix propensity as an
	// energy term or something.
	///////////////////////////////////////////////////////////////
	 void
	 StepWiseScreener::filter_native_BIG_BINS(
					 Size const & n,
					 MainChainTorsionSetList & main_chain_torsion_set_list )
	 {

		 if ( get_native_pose() == 0 ) {
			 utility_exit_with_message(  "Filter based on native big bins, but not native specified?" );
		 }

		 //big bin not really specified for end residues...
		 pose::Pose const & native_pose( *get_native_pose() );

		 if ( n == 1 ) return;
		 if ( n == native_pose.total_residue() ) return;

		 Size const big_bin =	 get_big_bin( native_pose.phi(n), native_pose.psi(n) );
		 if ( big_bin == 3 ) return; // If loop, no constraints.

		 MainChainTorsionSetList main_chain_torsion_set_list_new;

		 for ( Size n = 1; n <= main_chain_torsion_set_list.size(); n++ ) {
			 Real const phi = main_chain_torsion_set_list[ n ].phi();
			 Real const psi = main_chain_torsion_set_list[ n ].psi();
			 if ( ( big_bin == 1 && get_big_bin( phi, psi ) == 1 ) ||
						( big_bin == 2 && get_big_bin( phi, psi ) == 2 ) ) {
				 main_chain_torsion_set_list_new.push_back( main_chain_torsion_set_list[ n ] );
			 }
		 }

		 main_chain_torsion_set_list = main_chain_torsion_set_list_new;

	 }

	///////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::filter_main_chain_torsion_sets(){

		std::list< std::pair< Real, Size > > energy_index_list;
		for (Size n = 1; n <= centroid_scores_.size(); n++ ) {
			energy_index_list.push_back( std::make_pair( centroid_scores_[n], n ) ) ;
		}
		energy_index_list.sort();

		utility::vector1< MainChainTorsionSetList > main_chain_torsion_sets_for_moving_residues_new;

		for (Size n = 1; n <= centroid_scores_.size(); n++ ) {
			if ( n > nstruct_centroid_ ) break;
			main_chain_torsion_sets_for_moving_residues_new.push_back(
					 main_chain_torsion_sets_for_moving_residues_[ n ] );
		}

		main_chain_torsion_sets_for_moving_residues_ =  main_chain_torsion_sets_for_moving_residues_new;
	}


  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseScreener::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_rmsd_cutoff( core::Real const & setting ){
		rmsd_cutoff_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_n_sample( core::Size const & setting ){
		n_sample_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_nstruct_centroid( core::Size const & setting ){
		nstruct_centroid_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_centroid_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		centroid_scorefxn_ = scorefxn;
	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_filter_native_big_bins( bool const & setting ){
		filter_native_big_bins_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_centroid_screen( bool const & setting ){
		centroid_screen_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_ghost_loops( bool const & setting ){
		ghost_loops_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_apply_vdw_cut( bool const & setting ){
		apply_vdw_cut_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseScreener::set_centroid_score_diff_cut( core::Real const & setting ){
		centroid_score_diff_cut_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	utility::vector1< MainChainTorsionSetList > const &
	StepWiseScreener::main_chain_torsion_set_lists() const
	{
		return main_chain_torsion_sets_for_moving_residues_;
	}



} //protein
} //enumerate
} //stepwise
} //protocols
