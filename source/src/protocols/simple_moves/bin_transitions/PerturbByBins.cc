// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/bin_transitions/PerturbByBins.cc
/// @brief  This mover takes a stretch of backbone and perturbs its mainchain torsions based on the probabilities of transitions from
/// one torsion bin to another.
/// @details Bin transitions are read from database files.  The algorithm is: set the first residue based on the probability of a residue
/// being in a bin.  Set subsequent residues based on the probability of a residue being in a bin given that the previous residue is in
/// a particular bin.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit headers
#include <protocols/simple_moves/bin_transitions/PerturbByBins.hh>
#include <protocols/simple_moves/bin_transitions/PerturbByBinsCreator.hh>

// Bin transition calculator headers:
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>
#include <core/scoring/bin_transitions/BinTransitionData.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
//parsing
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <numeric/random/random.hh>

namespace protocols {
	namespace simple_moves {
		namespace bin_transitions {

			static thread_local basic::Tracer TR( "protocols.simple_moves.bin_transitions.PerturbByBins" );

			std::string
			PerturbByBinsCreator::keyname() const
			{
				return PerturbByBinsCreator::mover_name();
			}

			protocols::moves::MoverOP
			PerturbByBinsCreator::create_mover() const {
				return protocols::moves::MoverOP( new PerturbByBins );
			}

			std::string
			PerturbByBinsCreator::mover_name()
			{
				return "PerturbByBins";
			}

			///@brief Default constructor
			///
			PerturbByBins::PerturbByBins() : //TODO: initialize variables here!
				protocols::moves::Mover("PerturbByBins"),
				start_res_(0),
				end_res_(0),
				binfile_("ABBA"),
				binfile_loaded_(false),
				bin_transition_calculator_( new core::scoring::bin_transitions::BinTransitionCalculator ),
				repeats_(1),
				must_switch_bins_(false)
			{}
			
			/// @brief Copy constructor.
			///
			PerturbByBins::PerturbByBins( PerturbByBins const &src ) :
				protocols::moves::Mover( src ),
				start_res_(src.start_res_),
				end_res_(src.end_res_),
				binfile_(src.binfile_),
				binfile_loaded_(src.binfile_loaded_),
				bin_transition_calculator_( src.bin_transition_calculator_->clone() ), //CLONE this object.
				repeats_(src.repeats_),
				must_switch_bins_(src.must_switch_bins_)
			{}
			
			/// @brief Destructor.
			///
			PerturbByBins::~PerturbByBins() {}
			
			std::string
			PerturbByBins::get_name() const {
				return PerturbByBinsCreator::mover_name();
			}

			/// @brief Apply the mover to a pose.
			///
			void PerturbByBins::apply( core::pose::Pose & pose ) {			
				//Check whether we've loaded bin transitions.
				if( !bin_transition_calculator_->bin_params_loaded() ) utility_exit_with_message(
					"In protocols::simple_moves::bin_transitions::PerturbByBins::apply(): Bin transition probability parameters must be loaded before calling the apply() function!");

				//if(TR.visible()) TR << bin_transition_calculator_->summarize_stored_data(false); //DELETE ME!
			
				//Number of residues in the pose
				core::Size const nres( pose.n_residue() );
				if( nres<1 ) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::PerturbByBins::apply(): The pose has no residues!" );
				
				core::Size const startres( start_res_==0 ? 2 : start_res_ );
				core::Size const endres ( end_res_==0 ? nres-1 : end_res_ );
				if( endres < startres ) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::PerturbByBins::apply(): The last residue cannot come before the first!" );
				if( startres > nres) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::PerturbByBins::apply(): The first residue is out of range!  It must lie within [1, n_residues]." );
				if( endres > nres) utility_exit_with_message( "In protocols::simple_moves::bin_transitions::PerturbByBins::apply(): The last residue is out of range!  It must lie within [1, n_residues]." );

				// ITERATIONS of application:
				for(core::Size irepeat=1, nrepeats=repeats(); irepeat<=nrepeats; ++irepeat) {				
					utility::vector1 < core::Real > mainchain_torsions; //Vector of mainchain torsion values for a single residue that will be populated by the BinTransitionCalculator object.
					core::Size res_index( numeric::random::rg().random_range(startres, endres) ); //Index of the residue that will be set my the mover this iteration.  This is randomly drawn from the range [startres, endres] inclusive.
				
					//Actually generate the mainchain torsion values:
					bin_transition_calculator_->random_mainchain_torsions_using_adjacent_bins( pose.conformation(), res_index, must_switch_bins(), mainchain_torsions) ; //pose.conformation() and res_index are const inputs; mainchain_torsions is the output.
								
					{ //Scope to set the mainchain torsions:
						assert(mainchain_torsions.size()==pose.residue(res_index).mainchain_torsions().size()); //Should be true at this point.
						for(core::Size j=1, jmax=mainchain_torsions.size(); j<=jmax; ++j) { //Set mainchain torsions
							pose.set_torsion( core::id::TorsionID( res_index, core::id::BB, j ) , mainchain_torsions[j] );
						}				
					}
				} //End ITERATIONS of application
				
				//Final housekeeping:
				pose.update_residue_neighbors();
				if(TR.visible()) TR.flush();
				
				return;
			} //apply()

			/// @brief Parse XML for RosettaScripts.
			///
			void PerturbByBins::parse_my_tag( utility::tag::TagCOP tag,
					basic::datacache::DataMap &,
					protocols::filters::Filters_map const &,
					protocols::moves::Movers_map const &,
					Pose const &//pose
			)
			{
				if ( tag->getName() != "PerturbByBins" ){
					throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
				}
				
				//Get the bin params file:
				set_binfile_and_load( tag->getOption< std::string >( "bin_params_file", "ABBA" ) );
				
				//Get the residue ranges:
				set_residue_range( tag->getOption< core::Size >("start", 0), tag->getOption< core::Size >("end", 0) );
				
				//Get the number of repeats:
				set_repeats( tag->getOption< core::Size >("repeats", 1) );
				
				//Set whether the perturber always changes the bin of the residue it's acting on or not.
				set_must_switch_bins( tag->getOption< bool >("must_switch_bins", false) );
				
				if(TR.visible()) TR.flush();

				return;
			} //parse_my_tag()

			/// @brief Set the bin transition probability file.
			///
			void PerturbByBins::set_binfile_and_load( std::string const &name ) {
				if(binfile_loaded_==true) utility_exit_with_message(
					"In protocols::simple_moves::bin_transitions::PerturbByBins::set_binfile_and_load(): The bin params file was already loaded!  This operation cannot be repeated.");
				binfile_=name; 
				if(TR.visible()) TR << "Set bin params file name to \"" << binfile_ << "\"." << std::endl;
				bin_transition_calculator_->load_bin_params(binfile_);
				binfile_loaded_=true;
				return;
			} //set_binfile_and_load()
			
			/// @brief Set the residue ranges.  If set to (0,0), the
			/// start and end of the pose are used as the range bounds.
			void PerturbByBins::set_residue_range( core::Size const start, core::Size const end )
			{
				if( (end != 0) && (end < start) ) utility_exit_with_message(
					"In protocols::simple_moves::bin_transitions::PerturbByBins::set_residue_range(): The end residue cannot be before the start.");
				start_res_=start;
				end_res_=end;
				if(TR.visible()) {
					TR << "Set bin start and end ranges to ";
					if(start_res_!=0) TR << start_res_; else TR << "[start of pose+1]";
					TR << ", ";
					if(end_res_!=0) TR << end_res_; else TR << "[end of pose-1]";
					TR << std::endl;
				}
				return;
			} //set_residue_range
			
			/// @brief Set the number of repeats.
			/// @details A value of 1 means that a single residue in the range will randomly be selected and
			/// flipped to another bin (with probabilities based on its neighbours and the bin transition
			/// probabilities.  Higher values mean that this operation will be repeated.  If set to 0, defaults
			/// to 1.
			void PerturbByBins::set_repeats( core::Size const repeats_in ) {
				if(repeats_in == 0) {
					repeats_=1;
					if(TR.Warning.visible()) {
						TR.Warning << "Warning!  The number of repeats was set to 0.  Setting to 1." << std::endl; 
						TR.Warning.flush();
					}
				} else {
					repeats_=repeats_in;
					if(TR.visible()) {
						TR << "Setting number of repeats to " << repeats_ << "." << std::endl;
					}
				}
				return;
			} //set_repeats()
			
			/// @brief Sets whether the residue that is being perturbed can stay within its own bin (in which case new mainchain
			/// torsions are drawn from within the bin), or whether it must jump to a different bin.  True means it must jump.
			void PerturbByBins::set_must_switch_bins( bool const val ) {
				must_switch_bins_=val;
				if(TR.visible()) {
					if(must_switch_bins_) TR << "Perturbations set to ALWAYS change the bin of the residue being switched." << std::endl;
					else TR << "Perturbations set so that the residue being switched has some probability of staying within the bin in which it started (in which case its mainchain torsions are chosen randomly within the bin)." << std::endl;
					TR.flush();
				}
				return;
			} //set_must_switch_bins()
			
		} // bin_transitions
	} // simple_moves
} // protocols
