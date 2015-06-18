// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/bin_transitions/BinTransitionData.hh
/// @brief  Headers for BinTransitionData class.
/// @details This class stores data associated with transitions from one mainchain torsion bin to another (e.g. ABEGO bins, OO-ABBA bins, etc.)
/// for ONE specific type of transition (ith residue has certain properties, i+1st residue has certain other properties).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_scoring_bin_transitions_BinTransitionData_hh
#define INCLUDED_core_scoring_bin_transitions_BinTransitionData_hh

//BinTransitionData owning pointers header:
#include <core/scoring/bin_transitions/BinTransitionData.fwd.hh>

//Other headers:
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

//Ramachandran -- needed for sampling
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/AA.hh>

///////////////////////////////////////////////////////////////////////

namespace core {
	namespace scoring {
		namespace bin_transitions {
		
			/// @brief Enum for the properties that can be required or prohibited at the i or i+1 position.
			/// @details These reflect properties that residues can have.  Whenever a property is added, add
			/// it (1) here, (2) to the get_property_effect_name function, and (3) to the has_property function.
			enum BT_PROPERTIES {
				BT_PROTEIN=1,
				BT_L_AA,
				BT_D_AA,
				BT_ALPHA_AA,
				BT_BETA_AA,
				BT_POLAR,
				BT_METALBINDING,
				BT_CHARGED,
				BT_AROMATIC,
				BT_DISULFIDE_BONDED,
				BT_FORMS_DISULFIDE_BOND,
				BT_CYCLIC,
				
				BT_UNKNOWN_PROPERTY, //Keep this second-to-last
				BT_END_OF_LIST //Keep this last
			};
			
			/// @brief Enum for the way in which bins will be divided into subbins.
			/// @details  Default is BTSB_NONE, which means that the bin has no subbins.  As subbin types are
			/// added, add them (1) here, and (2) to the get_subbin_type_name function, and (3) to the setup_subbin_type
			/// function, then add a function called by the setup_subbin_type function.
			enum BTSB_SUBBIN_TYPE {
				BTSB_NONE=1,
				BTSB_L_AA,
				BTSB_D_AA,
				BTSB_L_PRO,
				BTSB_D_PRO,
				BTSB_GLY,
				
				BTSB_UNKNOWN, //Keep second-to-last
				BTSB_END_OF_LIST //Keep last
			};

			class BinTransitionData : public utility::pointer::ReferenceCount {
				public: //Constructors and destructors:
					/// @brief Default constructor for BinTransitionData
					///
					BinTransitionData();

					/// @brief Copy constructor for BinTransitionData
					///
					BinTransitionData( BinTransitionData const &src );

					/// @brief Default destructor for BinTransitionData
					///
					~BinTransitionData();

					/// @brief Clone operation for BinTransitionData.
					/// @details Returns an owning pointer to a copy of this object.
					BinTransitionDataOP clone() const;
					
				private: //Private functions:

					/// @brief Am I within the bounds of a 1D bin?
					/// @details If min > max, I want to be OUTSIDE of the range [min, max].
					bool in_bin( core::Real const &min, core::Real const &max, core::Real const &val ) const {
						if(min < max) {
							return ((val >= min) && (val < max));
						} else {
							return ((val < max) || (val >= min));
						}
						return false; //Should never reach here.
					}

					/// @brief Adjust an angle to lie in the range (-180,180).
					///				
					core::Real set_in_range( core::Real const &angle) const {
						core::Real returnval = angle;
						while(returnval > 180.0) returnval-=360.0;
						while(returnval < -180.0) returnval+=360.0;
						return returnval;
					}
					
					/// @brief Check whether a value is in a list
					///
					bool is_in_list( BT_PROPERTIES const val, utility::vector1 < BT_PROPERTIES > const &list ) const  {
						for (core::Size i=1, imax=list.size(); i<=imax; ++i) {
							if(list[i]==val) return true;
						}
						return false;
					}
					
					/// @brief Check whether a string is in a list
					///
					bool is_in_list( std::string const &str, utility::vector1 < std::string > const &list ) const {
						for (core::Size i=1, imax=list.size(); i<=imax; ++i) {
							if(list[i]==str) return true;
						}
						return false;
					}

					/// @brief Set up the sub-bins and their cumulative probability distributions (if appropriate).
					///
					void set_up_subbins();
					
					/// @brief Set up the cumulative probability distributions for each bin at position i+1 given a bin at position i, and vice versa.
					///
					void set_up_bin_cdfs();

				public: //Public functions -- calculators and evaluators:
				
					/// @brief Ensure that a sub-bin doesn't exceed the bounds of a bin.  If it does, scale its back, and reduce its contribution
					/// to the cumulative distribution function proportionately.
					void trim_subbin_edges_and_rescale_subbin(
						core::Real &curbin_cdf_val,
						std::pair <core::Real,core::Real> &tors_range,
						core::Real const &phimin,
						core::Real const &phimax
					) const;

					/// @brief Given a sub-bin type enum, return its name as it would appear in
					/// a residue params file.				
					std::string get_subbin_type_name (BTSB_SUBBIN_TYPE const type) const;
					
					/// @brief Given a sub-bin type name, return the sub-bin type enum value.
					///	@details Returns BTSB_UNKNOWN if not identifiable.
					BTSB_SUBBIN_TYPE get_subbin_type_from_name( std::string const &name ) const;
				
					/// @brief Given a property enum value, return the property name as it would appear in a
					/// residue params file.
					std::string get_property_effect_name( BT_PROPERTIES const property) const;

					/// @brief Given a property name, return the property enum value.
					///	@details Returns BT_UNKNOWN_PROPERTY if not identifiable.
					BT_PROPERTIES get_property_from_name( std::string const &name ) const;

					/// @brief Given a property, check whether a residue has that property.
					///	@details This function provides the link to the suitable property lookup function in ResidueType.
					bool has_property( BT_PROPERTIES const property, core::conformation::Residue const &rsd ) const;
					
					/// @brief Given a residue, check whether its properties match the required and prohibited properties lists for residue i,
					/// and whether its name matches the required and prohibited names list for residue i.
					bool criteria_match_i( core::conformation::Residue const &rsd ) const;
					
					/// @brief Given a residue, check whether its properties match the required and prohibited properties lists for residue i+1,
					/// and whether its name matches the required and prohibited names list for residue i+1.
					bool criteria_match_iplus1( core::conformation::Residue const &rsd ) const;
					
					/// @brief Select a random bin from the bins allowed for residue i, biased by the relative counts for that bin.
					///
					core::Size random_bin_i() const;
					
					/// @brief Select a random bin from the bins allowed for residue i+1, biased by the relative counts for that bin.
					///
					core::Size random_bin_iplus1() const;
					
					/// @brief Select a random bin from the bins allowed for residue i+1, biased by the relative counts for that bin given that residue i is in a particular bin.
					///
					core::Size random_bin_iplus1_given_i( core::Size const bin_i ) const;
					
					/// @brief Given a vector of mainchain torsions for a particular residue (at position i), figure out which bin the torsions lie in.
					///
					core::Size which_bin_i( utility::vector1 < core::Real > const &mainchain_torsions ) const;
					
					/// @brief Given a vector of mainchain torsions for a particular residue (at position i+1), figure out which bin the torsions lie in.
					///
					core::Size which_bin_iplus1( utility::vector1 < core::Real > const &mainchain_torsions ) const;
					
					/// @brief Is a given residue within the bounds of a given bin?
					/// @details Uses bin definitions for the ith residue.
					bool in_bin_i( core::Size const bin_index, core::conformation::Residue const &rsd  ) const;

					/// @brief Is a given residue within the bounds of a given bin?
					/// @details Uses bin definitions for the i+1st residue.
					bool in_bin_iplus1( core::Size const bin_index, core::conformation::Residue const &rsd  ) const;

				public: //Public functions -- cleanups and checks:
				
				/// @brief Do final post-load calculations (e.g. precomputing the sum of the transition probability matrix entries,
				/// the sum of the columns, the sum of the rows, etc.)
				void finalize();
				
				/// @brief Writes a report summarizing the data stored in this BinTransitionData object.
				/// @details If verbose is true, the full set of sub-bins is written out, too.
				std::string summarize_data( bool const verbose ) const;

				/// @brief Checks through the required and prohibited properties lists, and asserts that there is no overlap.
				/// @details Throws an error if there is overlap.				
				void check_property_overlap() const {
					if(properties_i_.size()>0 && prohibited_properties_i_.size()>0) {
						std::string const errmsg( "In core::scoring::bin_transitions::BinTransitionData::check_property_overlap(): The same property has been added to the required and prohibited properties lists for residue i." );
						for(core::Size i=1, imax=properties_i_.size(); i<=imax; ++i) {
							for(core::Size j=1,jmax=prohibited_properties_i_.size(); j<=jmax; ++j) {
								runtime_assert_string_msg(properties_i_[i]!=prohibited_properties_i_[j], errmsg);
							}
						}
					}
					if(properties_iplus1_.size()>0 && prohibited_properties_iplus1_.size()>0) {
						std::string const errmsg( "In core::scoring::bin_transitions::BinTransitionData::check_property_overlap(): The same property has been added to the required and prohibited properties lists for residue i+1." );
						for(core::Size i=1, imax=properties_iplus1_.size(); i<=imax; ++i) {
							for(core::Size j=1,jmax=prohibited_properties_iplus1_.size(); j<=jmax; ++j) {
								runtime_assert_string_msg(properties_iplus1_[i]!=prohibited_properties_iplus1_[j], errmsg);
							}
						}
					}
				} //check_property_overlap
				
				/// @brief Checks through the required and prohibited residue identities lists, and asserts that there is no overlap.
				/// @details Throws an error if there is overlap.
				void check_residentity_overlap() const {
					if(res_identities_i_.size()>0 && prohibited_res_identities_i_.size()>0) {
						std::string const errmsg( "In core::scoring::bin_transitions::BinTransitionData::check_residentity_overlap(): The same property has been added to the required and prohibited res identity lists for residue i." );
						for(core::Size i=1, imax=res_identities_i_.size(); i<=imax; ++i) {
							for(core::Size j=1,jmax=prohibited_res_identities_i_.size(); j<=jmax; ++j) {
								runtime_assert_string_msg(res_identities_i_[i]!=prohibited_res_identities_i_[j], errmsg);
							}
						}
					}
					if(res_identities_iplus1_.size()>0 && prohibited_res_identities_iplus1_.size()>0) {
						std::string const errmsg( "In core::scoring::bin_transitions::BinTransitionData::check_residentity_overlap(): The same property has been added to the required and prohibited res identity lists for residue i+1." );
						for(core::Size i=1, imax=res_identities_iplus1_.size(); i<=imax; ++i) {
							for(core::Size j=1,jmax=prohibited_res_identities_iplus1_.size(); j<=jmax; ++j) {
								runtime_assert_string_msg(res_identities_iplus1_[i]!=prohibited_res_identities_iplus1_[j], errmsg);
							}
						}
					}
				} //check_residentities_overlap
				
				public: //Public functions -- getters:

					/// @brief Have the n_bins_i_ and n_bins_iplus1_ variables been set, so that the
					/// transition probability matrix could be initialized?
					bool matrix_initialized() const { return matrix_initialized_; }
					
					/// @brief Have final calculations been done after initializing the probability matrix?
					///
					bool matrix_finalized() const { return matrix_finalized_; }
					
					/// @brief Given a bin name, get the bin index.
					/// @details Returns 0 if the name is not found.  Checks the list of bin names
					/// for the ith residue.
					core::Size binname_index_from_string_i( std::string const &bin_name ) const {
						for(core::Size i=1, imax=binnames_i_.size(); i<=imax; ++i) {
							if(binnames_i_[i]==bin_name) return i;
						}
						return 0;
					}

					/// @brief Given a bin name, get the bin index.
					/// @details Returns 0 if the name is not found.  Checks the list of bin names
					/// for the i+1st residue.
					core::Size binname_index_from_string_iplus1( std::string const &bin_name ) const {
						for(core::Size i=1, imax=binnames_iplus1_.size(); i<=imax; ++i) {
							if(binnames_iplus1_[i]==bin_name) return i;
						}
						return 0;
					}

					/// @brief Get the number of bins for the ith residue.
					///
					core::Size n_bins_i() const { return n_bins_i_; }
					
					/// @brief Get the number of bins for the i+1st residue.
					///
					core::Size n_bins_iplus1() const { return n_bins_iplus1_; }

					/// @brief Get the number of mainchain torsions for the ith residue.
					///
					core::Size n_mainchain_torsions_i() const { return n_mainchain_torsions_i_; }

					/// @brief Get the number of mainchain torsions for the i+1st residue.
					///
					core::Size n_mainchain_torsions_iplus1() const { return n_mainchain_torsions_iplus1_; }
					
					/// @brief Get the bin start and end boundaries, in degrees, for the nth bin and the mth mainchain torsion for residue i.
					/// @details If the end is less than the start, it means that the bin runs through 180 degrees and wraps around.
					std::pair <core::Real, core::Real> bin_boundaries_i( core::Size const n, core::Size const m ) const {
						if(n<1 || n>binranges_i_.size()) utility_exit_with_message( "In BinTransitionData::bin_boundaries_i(): The index given for the bin is out of the range [1,number_of_bins].\n" );
						if(m<1 || m>binranges_i_[n].size()) utility_exit_with_message( "In BinTransitionData::bin_boundaries_i(): The index given for the mainchain torsion is out of the range [1,number_of_mainchain_torsions].\n" );
						return binranges_i_[n][m];
					}
					
					/// @brief Get the bin start and end boundaries, in degrees, for the nth bin and the mth mainchain torsion for residue i+1.
					/// @details If the end is less than the start, it means that the bin runs through 180 degrees and wraps around.
					std::pair <core::Real, core::Real> bin_boundaries_iplus1( core::Size const n, core::Size const m ) const {
						if(n<1 || n>binranges_iplus1_.size()) utility_exit_with_message( "In BinTransitionData::bin_boundaries_iplus1(): The index given for the bin is out of the range [1,number_of_bins].\n" );
						if(m<1 || m>binranges_iplus1_[n].size()) utility_exit_with_message( "In BinTransitionData::bin_boundaries_iplus1(): The index given for the mainchain torsion is out of the range [1,number_of_mainchain_torsions].\n" );
						return binranges_iplus1_[n][m];
					}

					
					/// @brief Get an entry in the transition probability matrix.
					/// @details [n][m] = bin of ith residue, bin of i+1st residue.
					core::Real probability_matrix( core::Size const n, core::Size const m ) const {
						runtime_assert_string_msg( matrix_initialized(), "In core::scoring::bin_transitions::BinTransitionData::probability_matrix(): The matrix must be initialized before calling this function.");
						runtime_assert_string_msg( n>0 && m>0, "In core::scoring::bin_transitions::BinTransitionData::probability_matrix(): The matrix indices must be greater than zero.");
						runtime_assert_string_msg( probability_matrix_.size() >= n,
							"In core::scoring::bin_transitions::BinTransitionData::probability_matrix(): The matrix has fewer lines than the line specified." );
						runtime_assert_string_msg( probability_matrix_[n].size() >= m,
							"In core::scoring::bin_transitions::BinTransitionData::probability_matrix(): The matrix has fewer columns than the column specified." );
						return probability_matrix_[n][m];
					}
					
					/// @brief Get the sub-bin type enum for residue i.
					///
					BTSB_SUBBIN_TYPE subbin_type_i() const { return subbin_type_i_; }
					
					/// @brief Get the sub-bin type enum for residue i+1.
					///
					BTSB_SUBBIN_TYPE subbin_type_iplus1() const { return subbin_type_iplus1_; }
					
					/// @brief Get the sub-bin cumulative distribution function for bin n for the ith residue.
					/// @details Const-access only.
					utility::vector1 < core::Real > const & subbin_cdf_i( core::Size const n ) const {
						if(n < 1 || n > subbin_cdf_i_.size()) utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_cdf_i(): The bin index is out of range.\n" );
						return subbin_cdf_i_[n];
					}
					
					/// @brief Get the sub-bin cumulative distribution function for bin n for the i+1st residue.
					/// @details Const-access only.
					utility::vector1 < core::Real > const & subbin_cdf_iplus1( core::Size const n ) const {
						if(n < 1 || n > subbin_cdf_iplus1_.size()) utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_cdf_iplus1(): The bin index is out of range.\n" );
						return subbin_cdf_iplus1_[n];
					}
					
					/// @brief Get the start of the sub-bin torsion range for the ith bin, the jth sub-bin, and the kth mainchain torsion.
					///
					core::Real subbin_start_i( core::Size const bin_i, core::Size const subbin_j, core::Size const mainchain_torsion_k) const
					{
						if( bin_i < 1 || bin_i > subbin_ranges_i_.size() || bin_i > n_bins_i())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_start_i(): The bin index is out of range.\n" );
						if( subbin_j < 1 || subbin_j > subbin_ranges_i_[bin_i].size())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_start_i(): The sub-bin index is out of range.\n" );
						if( mainchain_torsion_k < 1 || mainchain_torsion_k > subbin_ranges_i_[bin_i][subbin_j].size() || mainchain_torsion_k > n_mainchain_torsions_i())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_start_i(): The mainchain torsion index is out of range.\n" );
						return subbin_ranges_i_[bin_i][subbin_j][mainchain_torsion_k].first;
					}
					
					/// @brief Get the end of the sub-bin torsion range for the ith bin, the jth sub-bin, and the kth mainchain torsion.
					///
					core::Real subbin_end_i( core::Size const bin_i, core::Size const subbin_j, core::Size const mainchain_torsion_k) const
					{
						if( bin_i < 1 || bin_i > subbin_ranges_i_.size() || bin_i > n_bins_i())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_end_i(): The bin index is out of range.\n" );
						if( subbin_j < 1 || subbin_j > subbin_ranges_i_[bin_i].size())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_end_i(): The sub-bin index is out of range.\n" );
						if( mainchain_torsion_k < 1 || mainchain_torsion_k > subbin_ranges_i_[bin_i][subbin_j].size() || mainchain_torsion_k > n_mainchain_torsions_i())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_end_i(): The mainchain torsion index is out of range.\n" );
						return subbin_ranges_i_[bin_i][subbin_j][mainchain_torsion_k].second;
					}
					
					/// @brief Get the start of the sub-bin torsion range for the ith bin, the jth sub-bin, and the kth mainchain torsion.
					///
					core::Real subbin_start_iplus1( core::Size const bin_i, core::Size const subbin_j, core::Size const mainchain_torsion_k) const
					{
						if( bin_i < 1 || bin_i > subbin_ranges_iplus1_.size() || bin_i > n_bins_iplus1())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_start_iplus1(): The bin index is out of range.\n" );
						if( subbin_j < 1 || subbin_j > subbin_ranges_iplus1_[bin_i].size())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_start_iplus1(): The sub-bin index is out of range.\n" );
						if( mainchain_torsion_k < 1 || mainchain_torsion_k > subbin_ranges_iplus1_[bin_i][subbin_j].size() || mainchain_torsion_k > n_mainchain_torsions_iplus1())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_start_iplus1(): The mainchain torsion index is out of range.\n" );
						return subbin_ranges_iplus1_[bin_i][subbin_j][mainchain_torsion_k].first;
					}
					
					/// @brief Get the end of the sub-bin torsion range for the ith bin, the jth sub-bin, and the kth mainchain torsion.
					///
					core::Real subbin_end_iplus1( core::Size const bin_i, core::Size const subbin_j, core::Size const mainchain_torsion_k) const
					{
						if( bin_i < 1 || bin_i > subbin_ranges_iplus1_.size() || bin_i > n_bins_iplus1())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_end_iplus1(): The bin index is out of range.\n" );
						if( subbin_j < 1 || subbin_j > subbin_ranges_iplus1_[bin_i].size())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_end_iplus1(): The sub-bin index is out of range.\n" );
						if( mainchain_torsion_k < 1 || mainchain_torsion_k > subbin_ranges_iplus1_[bin_i][subbin_j].size() || mainchain_torsion_k > n_mainchain_torsions_iplus1())
							utility_exit_with_message( "In core::scoring::bin_transitions::BinTransitionData::subbin_end_iplus1(): The mainchain torsion index is out of range.\n" );
						return subbin_ranges_iplus1_[bin_i][subbin_j][mainchain_torsion_k].second;
					}
					
					/// @brief Get the sum of a row of the probability matrix.
					/// @details This is the sum of all counts(bin_i, bin_iplus1), summed across all bin_iplus1.  This is pre-calculated and stored,
					/// so this is a fast lookup.
					core::Real binsums_i(core::Size const bin_index_i) const {
						if(bin_index_i<1 || bin_index_i>binsums_i_.size()) utility_exit_with_message(
							"In core::scoring::bin_transitions::BinTransitionData::binsums_i(): The index is out of range.");
						return binsums_i_[bin_index_i];
					}
					
					/// @brief Get the sum of a column of the probability matrix.
					/// @details This is the sum of all counts(bin_i, bin_iplus1), summed across all bin_i.  This is pre-calculated and stored,
					/// so this is a fast lookup.
					core::Real binsums_iplus1(core::Size const bin_index_iplus1) const {
						if(bin_index_iplus1<1 || bin_index_iplus1>binsums_iplus1_.size()) utility_exit_with_message(
							"In core::scoring::bin_transitions::BinTransitionData::binsums_iplus1(): The index is out of range.");
						return binsums_iplus1_[bin_index_iplus1];
					}


				public: //Public functions -- setters:

					/// @brief Set the number of mainchain torsions for residue i.
					///
					void set_n_mainchain_torsions_i( core::Size const val ) {
						runtime_assert_string_msg( val > 0, "In core::scoring::bin_transitions::BinTransitionData::set_n_mainchain_torsions_i(): The input value must be greater than 0." );
						n_mainchain_torsions_i_=val;
						return;
					};

					/// @brief Set the number of mainchain torsions for residue i+1.
					///
					void set_n_mainchain_torsions_iplus1( core::Size const val ) {
						runtime_assert_string_msg( val > 0, "In core::scoring::bin_transitions::BinTransitionData::set_n_mainchain_torsions_iplus1(): The input value must be greater than 0." );
						n_mainchain_torsions_iplus1_=val;
						return;
					};

					/// @brief Set the number of bins for the ith and i+1st residues.
					/// @details This also initializes the probability matrix (to all zeros), the binnames vectors
					/// (to vectors of empty strings), and the binranges vectors (to all zeros).
					void set_n_bins( core::Size const n_bins_i, core::Size const n_bins_iplus1 );
					
					/// @brief Set the name of bin n for residue i.
					/// @details This requires that the number of bins has already been set.
					void set_binname_i( core::Size const bin, std::string const name ) {
						runtime_assert_string_msg( matrix_initialized(), "In core::scoring::bin_transitions::BinTransitionData::set_binname_i(): This function requires that the set_n_bins() function be called first." );
						runtime_assert_string_msg( bin<=n_bins_i(), "In core::scoring::bin_transitions::BinTransitionData::set_binname_i(): The specified bin is greater than the number of bins for residue i." );
						binnames_i_[bin] = name;
						return;
					}
					
					/// @brief Set the name of bin n for residue i+1.
					/// @details This requires that the number of bins has already been set.
					void set_binname_iplus1( core::Size const bin, std::string const name ) {
						runtime_assert_string_msg( matrix_initialized(), "In core::scoring::bin_transitions::BinTransitionData::set_binname_iplus1(): This function requires that the set_n_bins() function be called first." );
						runtime_assert_string_msg( bin<=n_bins_iplus1(), "In core::scoring::bin_transitions::BinTransitionData::set_binname_iplus1(): The specified bin is greater than the number of bins for residue i+1." );
						binnames_iplus1_[bin] = name;
						return;
					}
					
					/// @brief Set the torsion ranges for bin n, mainchain torsion m of residue i.
					/// @details This requires that the number of bins and the number of mainchain torsions have already been set.
					void set_binrange_i( core::Size const bin, core::Size const torsion, core::Real const &start, core::Real const &end ) {
						runtime_assert_string_msg( matrix_initialized(), "In core::scoring::bin_transitions::BinTransitionData::set_binrange_i(): This function requires that the set_n_bins() function be called first." );
						runtime_assert_string_msg( bin<=n_bins_i(), "In core::scoring::bin_transitions::BinTransitionData::set_binrange_i(): The specified bin is greater than the number of bins for residue i." );
						runtime_assert_string_msg( torsion<=n_mainchain_torsions_i(), "In core::scoring::bin_transitions::BinTransitionData::set_binrange_i(): The specified mainchain torsion is greater than the number of mainchain torsions for residue i." );
						
						binranges_i_[bin][torsion].first=set_in_range(start); //Adjust the start angle to lie between (-180, 180)
						binranges_i_[bin][torsion].second=set_in_range(end); //Adjust the end angle to lie between (-180, 180)

						return;
					}
					
					/// @brief Set the torsion ranges for bin n, mainchain torsion m of residue i+1.
					/// @details This requires that the number of bins and the number of mainchain torsions have already been set.
					void set_binrange_iplus1( core::Size const bin, core::Size const torsion, core::Real const &start, core::Real const &end ) {
						runtime_assert_string_msg( matrix_initialized(), "In core::scoring::bin_transitions::BinTransitionData::set_binrange_iplus1(): This function requires that the set_n_bins() function be called first." );
						runtime_assert_string_msg( bin<=n_bins_iplus1(), "In core::scoring::bin_transitions::BinTransitionData::set_binrange_iplus1(): The specified bin is greater than the number of bins for residue i+1." );
						runtime_assert_string_msg( torsion<=n_mainchain_torsions_iplus1(), "In core::scoring::bin_transitions::BinTransitionData::set_binrange_iplus1(): The specified mainchain torsion is greater than the number of mainchain torsions for residue i+1." );
						
						binranges_iplus1_[bin][torsion].first=set_in_range(start); //Adjust the start angle to lie between (-180, 180)
						binranges_iplus1_[bin][torsion].second=set_in_range(end); //Adjust the end angle to lie between (-180, 180)

						return;
					}
					
					/// @brief Copy the bin names and torsion ranges for residue i to those for residue i+1.
					///
					void copy_i_bins_to_iplus1() {
						runtime_assert_string_msg( n_bins_iplus1()==n_bins_i(), "In core::scoring::bin_transitions::BinTransitionData::copy_i_bins_to_iplus1(): The \"IPLUS1_BINS_COPY_I\" option cannot be used if the i and i+1 residues have different numbers of bins!" );
						runtime_assert_string_msg( n_mainchain_torsions_iplus1()==n_mainchain_torsions_i(), "In core::scoring::bin_transitions::BinTransitionData::copy_i_bins_to_iplus1(): The \"IPLUS1_BINS_COPY_I\" option cannot be used if the i and i+1 residues have different numbers of mainchain torsions!" );
						binranges_iplus1_ = binranges_i_;
						binnames_iplus1_=binnames_i_;
						subbin_ranges_iplus1_=subbin_ranges_i_;
						subbin_cdf_iplus1_=subbin_cdf_i_;
						subbin_type_iplus1_=subbin_type_i_;
						return;
					}

					/// @brief Set a value in the transition probability matrix.
					///
					void set_matrix_entry(core::Size const n, core::Size const m, core::Real const &val)
					{
						runtime_assert_string_msg( matrix_initialized(), "In core::scoring::bin_transitions::BinTransitionData::set_matrix_entry(): The matrix must be initialized before calling this function.");
						runtime_assert_string_msg( n>0 && m>0, "In core::scoring::bin_transitions::BinTransitionData::set_matrix_entry(): The matrix indices must be greater than zero.");
						runtime_assert_string_msg( probability_matrix_.size() >= n,
							"In core::scoring::bin_transitions::BinTransitionData::set_matrix_entry(): The matrix has fewer lines than the line specified." );
						runtime_assert_string_msg( probability_matrix_[n].size() >= m,
							"In core::scoring::bin_transitions::BinTransitionData::set_matrix_entry(): The matrix has fewer columns than the column specified." );
						runtime_assert_string_msg( val >= 0.0, "In core::scoring::bin_transitions::BinTransitionData::set_matrix_entry(): All values must be positive." );
						probability_matrix_[n][m]=val;
						return;
					}
					
					/// @brief Add a property to the list of properties that the ith residue MUST have.
					///
					void add_property_i( std::string const &property ) {
						BT_PROPERTIES prop = get_property_from_name(property);
						runtime_assert_string_msg(prop!=BT_UNKNOWN_PROPERTY,
							"In core::scoring::bin_transitions::BinTransitionData::add_property_i(): The property to be added could not be identified from the string \"" + property + "\".");
						if(!is_in_list(prop, properties_i_)) properties_i_.push_back(prop);
						return;
					};

					/// @brief Add a property to the list of properties that the ith residue MUST NOT have.
					///
					void prohibit_property_i( std::string const &property ) {
						BT_PROPERTIES prop = get_property_from_name(property);
						runtime_assert_string_msg(prop!=BT_UNKNOWN_PROPERTY,
							"In core::scoring::bin_transitions::BinTransitionData::prohibit_property_i(): The property to be added could not be identified from the string \"" + property + "\".");
						if(!is_in_list(prop, prohibited_properties_i_)) prohibited_properties_i_.push_back(prop);
						return;
					};
					
					/// @brief Add a property to the list of properties that the i+1st residue MUST have.
					///
					void add_property_iplus1( std::string const &property ) {
						BT_PROPERTIES prop = get_property_from_name(property);
						runtime_assert_string_msg(prop!=BT_UNKNOWN_PROPERTY,
							"In core::scoring::bin_transitions::BinTransitionData::add_property_iplus1(): The property to be added could not be identified from the string \"" + property + "\".");
						if(!is_in_list(prop, properties_iplus1_)) properties_iplus1_.push_back(prop);
						return;
					};

					/// @brief Add a property to the list of properties that the i+1st residue MUST NOT have.
					///
					void prohibit_property_iplus1( std::string const &property ) {
						BT_PROPERTIES prop = get_property_from_name(property);
						runtime_assert_string_msg(prop!=BT_UNKNOWN_PROPERTY,
							"In core::scoring::bin_transitions::BinTransitionData::prohibit_property_iplus1(): The property to be added could not be identified from the string \"" + property + "\".");
						if(!is_in_list(prop, prohibited_properties_iplus1_)) prohibited_properties_iplus1_.push_back(prop);
						return;
					};
					
					/// @brief Add a residue identity to the list of residue identities that the ith residue MUST have.
					/// @details If the list is left empty, the residue identities can be anything.  Residue names are
					/// provided as three-letter codes.
					void add_res_identity_i( std::string const &resname ) {
						runtime_assert_string_msg( resname.length()<=3,
							"In core::scoring::bin_transitions::BinTransitionData::add_res_identity_i(): The residue names must be provided as three-letter codes.");
						if(!is_in_list(resname, res_identities_i_)) res_identities_i_.push_back(resname);
						return;
					}
					
					/// @brief Add a residue identity to the list of residue identities that the ith residue MUST NOT have.
					/// @details If the list is left empty, the residue identities can be anything.  Residue names are
					/// provided as three-letter codes.
					void prohibit_res_identity_i( std::string const &resname ) {
						runtime_assert_string_msg( resname.length()<=3,
							"In core::scoring::bin_transitions::BinTransitionData::prohibit_res_identity_i(): The residue names must be provided as three-letter codes.");
						if(!is_in_list(resname, prohibited_res_identities_i_)) prohibited_res_identities_i_.push_back(resname);
						return;
					}

					/// @brief Add a residue identity to the list of residue identities that the i+1st residue MUST have.
					/// @details If the list is left empty, the residue identities can be anything.  Residue names are
					/// provided as three-letter codes.
					void add_res_identity_iplus1( std::string const &resname ) {
						runtime_assert_string_msg( resname.length()<=3,
							"In core::scoring::bin_transitions::BinTransitionData::add_res_identity_iplus1(): The residue names must be provided as three-letter codes.");
						if(!is_in_list(resname, res_identities_iplus1_)) res_identities_iplus1_.push_back(resname);
						return;
					}
					
					/// @brief Add a residue identity to the list of residue identities that the i+1st residue MUST NOT have.
					/// @details If the list is left empty, the residue identities can be anything.  Residue names are
					/// provided as three-letter codes.
					void prohibit_res_identity_iplus1( std::string const &resname ) {
						runtime_assert_string_msg( resname.length()<=3,
							"In core::scoring::bin_transitions::BinTransitionData::prohibit_res_identity_iplus1(): The residue names must be provided as three-letter codes.");
						if(!is_in_list(resname, prohibited_res_identities_iplus1_)) prohibited_res_identities_iplus1_.push_back(resname);
						return;
					}
					
					/// @brief Set the sub-bin type for residue i.
					/// @details Each bin can be divided into sub-bins in one of several ways.  For example,
					/// phi/psi/omega bins (like the ABEGO bins) can be divided into 5x5x90 bins corresponding
					/// to the bins used by the Rama scorefunction.  This is useful for drawing random phi/psi
					/// angles from within the bin.
					void set_subbin_type_i( std::string const &type ) {
						BTSB_SUBBIN_TYPE sbtype = get_subbin_type_from_name( type );
						runtime_assert_string_msg( sbtype != BTSB_UNKNOWN, "In core::scoring::bin_transitions::BinTransitionData::set_subbin_type_i(): Could not get type from name \"" + type + "\"." );
						subbin_type_i_=sbtype;
						return;
					} //set_subbin_type_i
					
					/// @brief Set the sub-bin type for residue i+1.
					/// @details Each bin can be divided into sub-bins in one of several ways.  For example,
					/// phi/psi/omega bins (like the ABEGO bins) can be divided into 5x5x90 bins corresponding
					/// to the bins used by the Rama scorefunction.  This is useful for drawing random phi/psi
					/// angles from within the bin.
					void set_subbin_type_iplus1( std::string const &type ) {
						BTSB_SUBBIN_TYPE sbtype = get_subbin_type_from_name( type );
						runtime_assert_string_msg( sbtype != BTSB_UNKNOWN, "In core::scoring::bin_transitions::BinTransitionData::set_subbin_type_iplus1(): Could not get type from name \"" + type + "\"." );
						subbin_type_iplus1_=sbtype;
						return;
					} //set_subbin_type_iplus1

				private: //Private variables:

					/// @brief Number of mainchain torsions in the ith residue.
					///
					core::Size n_mainchain_torsions_i_;

					/// @brief Number of mainchain torsions in the i+1st residue.
					///
					core::Size n_mainchain_torsions_iplus1_;

					/// @brief Number of bins for the ith residue.
					///
					core::Size n_bins_i_;

					/// @brief Torsion ranges for each bin for the ith residue.
					/// @details The outer vector is the bins.  The inner vector is the mainchain torsions.
					/// For each mainchain torsion, a pair is defined of (start of range, end of range).  If
					/// the end of the range is smaller than the start of the range, it is understood that the
					/// bin wraps around (e.g. 90, -90 means a bin running from 90, through 180, back to -180,
					/// and to -90).  Each bin INCLUDES the start of its range and EXCLUDES the end of its range
					/// (i.e. 25 35 means 25 <= angle < 35).
					utility::vector1 < utility::vector1 < std::pair< core::Real, core::Real > > > binranges_i_;

					/// @brief Bin names for the ith residue.
					///
					utility::vector1 <std::string> binnames_i_;

					/// @brief Number of bins for the i+1st residue.
					///
					core::Size n_bins_iplus1_;

					/// @brief Torsion ranges for each bin for the i+1st residue.
					/// @details The outer vector is the bins.  The inner vector is the mainchain torsions.
					/// For each mainchain torsion, a pair is defined of (start of range, end of range).  If
					/// the end of the range is smaller than the start of the range, it is understood that the
					/// bin wraps around (e.g. 90, -90 means a bin running from 90, through 180, back to -180,
					/// and to -90).  Each bin INCLUDES the start of its range and EXCLUDES the end of its range
					/// (i.e. 25 35 means 25 <= angle < 35).
					utility::vector1 < utility::vector1 < std::pair< core::Real, core::Real > > > binranges_iplus1_;

					/// @brief Bin names for the i + 1st residue.
					///
					utility::vector1 <std::string> binnames_iplus1_;

					/// @brief Have the n_bins_i_ and n_bins_iplus1_ variables been set, so that the
					/// transition probability matrix could be initialized?
					bool matrix_initialized_;

					/// @brief Have final calculations (e.g. normalization factor, column and row sums) been
					/// done after loading the transition probability matrix?
					bool matrix_finalized_;

					/// @brief The matrix of probability values.
					/// @brief probability_matrix_[x][y] represents the probability of the ith residue
					/// being in torsion bin x and the i+1st residue being in torsion bin y.  These are
					/// unnormalized probabilities.
					utility::vector1 < utility::vector1 < core::Real > > probability_matrix_;
					
					/// @brief List of properties that the ith residue MUST have.
					///
					utility::vector1 < BT_PROPERTIES > properties_i_;
					
					/// @brief List of properties that the ith residue MUST NOT have.
					///
					utility::vector1 < BT_PROPERTIES > prohibited_properties_i_;
					
					/// @brief List of properties that the i+1st residue MUST have.
					///
					utility::vector1 < BT_PROPERTIES > properties_iplus1_;
					
					/// @brief List of properties that the i+1st residue MUST NOT have.
					///
					utility::vector1 < BT_PROPERTIES > prohibited_properties_iplus1_;
					
					/// @brief List of residue identities that the ith residue MUST have.
					/// @details Vector of strings of three-letter codes.
					utility::vector1 <std::string> res_identities_i_;
					
					/// @brief List of residue identities that the ith residue MUST NOT have.
					/// @details Vector of strings of three-letter codes.
					utility::vector1 <std::string> prohibited_res_identities_i_;

					/// @brief List of residue identities that hte i+1st residue MUST have.
					/// @details Vector of strings of three-letter codes.
					utility::vector1 <std::string> res_identities_iplus1_;
					
					/// @brief List of residue identities that hte i+1st residue MUST NOT have.
					/// @details Vector of strings of three-letter codes.
					utility::vector1 <std::string> prohibited_res_identities_iplus1_;
					
					/// @brief Sum of counts for the ith residue bins (across all i+1st residue bins).
					///
					utility::vector1 <core::Real> binsums_i_;
					
					/// @brief Sum of counts for the i+1st residue bins (across all ith residue bins).
					///
					utility::vector1 <core::Real> binsums_iplus1_;

					/// @brief Sum of all entries in the probability matrix.
					///
					core::Real total_counts_;
					
					/// @brief Cumulative probability distribution for the ith residue bins.
					/// @details This ranges from 0 to 1.
					utility::vector1 <core::Real> binsums_i_cdf_;
					
					/// @brief Cumulative probability distribution for the i+1st residue bins.
					/// @details This ranges from 0 to 1.
					utility::vector1 <core::Real> binsums_iplus1_cdf_;
					
					/// @brief The type of sub-bins for residue i.
					/// @details Each bin can be divided into sub-bins in one of several ways.  For example,
					/// phi/psi/omega bins (like the ABEGO bins) can be divided into 5x5x90 bins corresponding
					/// to the bins used by the Rama scorefunction.  This is useful for drawing random phi/psi
					/// angles from within the bin.
					BTSB_SUBBIN_TYPE subbin_type_i_;
					
					/// @brief The type of sub-bins for residue i+1.
					/// @details Each bin can be divided into sub-bins in one of several ways.  For example,
					/// phi/psi/omega bins (like the ABEGO bins) can be divided into 5x5x90 bins corresponding
					/// to the bins used by the Rama scorefunction.  This is useful for drawing random phi/psi
					/// angles from within the bin.
					BTSB_SUBBIN_TYPE subbin_type_iplus1_;
										
					/// @brief Cumulative probability distribution for the ith residue's bins,
					/// WITHIN each bin.  (For example, based on ramachandran probabilities
					/// within each ABEGO bin).
					/// @details The outer vector represents the bin, and the inner vector
					/// represents the cumulative probability distribution for each sub-bin
					/// within that bin.
					utility::vector1 /*bins*/ < utility::vector1 /*sub-bins*/ < core::Real >  > subbin_cdf_i_;
					
					/// @brief The sub-bin identities for the ith residue's bins.
					/// @details For example, we might divide ABEGO bin A into sub-bins
					/// set up as a five-degree grid.  The outer vector represents the
					/// bin, the middle vector represents sub-bins, the inner vector represents
					/// mainchain torsions, and the innermost pairs represent torsion ranges for
					/// that mainchain torsion of that sub-bin.
					utility::vector1 /*bins*/ < utility::vector1 /*sub-bins*/ < utility::vector1 /*mainchain torsions*/ < std::pair < core::Real, core::Real> /*Start and end ranges*/ > > > subbin_ranges_i_;
					
					/// @brief Cumulative probability distribution for the i+1st residue's bins,
					/// WITHIN each bin.  (For example, based on ramachandran probabilities
					/// within each ABEGO bin).
					/// @details The outer vector represents the bin, and the inner vector
					/// represents the cumulative probability distribution for each sub-bin
					/// within that bin.
					utility::vector1 /*bins*/ < utility::vector1 /*sub-bins*/ < core::Real >  > subbin_cdf_iplus1_;
					
					/// @brief The sub-bin identities for the i+1st residue's bins.
					/// @details For example, we might divide ABEGO bin A into sub-bins
					/// set up as a five-degree grid.  The outer vector represents the
					/// bin, the middle vector represents sub-bins, the inner vector represents
					/// mainchain torsions, and the innermost pairs represent torsion ranges for
					/// that mainchain torsion of that sub-bin.
					utility::vector1 /*bins*/ < utility::vector1 /*sub-bins*/ < utility::vector1 /*mainchain torsions*/ < std::pair < core::Real, core::Real> /*Start and end ranges*/ > > > subbin_ranges_iplus1_;
					
					/// @brief Pre-calculated cumulative distribution functions for the probability of a bin at position i+1 (inner vector) given particular bins at position i (outer vector).
					///
					utility::vector1 /*bin_i*/ < utility::vector1 /*bin_i+1_cdf*/ < core::Real > > bin_iplus1_cdf_given_i_;

					/// @brief Pre-calculated cumulative distribution functions for the probability of a bin at position i (inner vector) given particular bins at position i+1 (outer vector).
					/// @details NOTE that indices are backwards to what's usual, here!
					utility::vector1 /*bin_i+1*/ < utility::vector1 /*bin_i_cdf*/ < core::Real > > bin_i_cdf_given_iplus1_;

			}; //BinTransitionData class

		} //namespace bin_transitions
	}//namespace scoring
}//namespace core

#endif //INCLUDED_core_scoring_bin_transitions_BinTransitionData_hh
