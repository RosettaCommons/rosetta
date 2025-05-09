// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/ReactionBasedAnalogSampler.hh
/// @brief Reaction-based and similarity guided analoging in a given chemical library
/// @author Yidan Tang (yidan.tang@vanderbilt.edu)

#ifndef INCLUDED_protocols_drug_design_ReactionBasedAnalogSampler_hh
#define INCLUDED_protocols_drug_design_ReactionBasedAnalogSampler_hh

#include <protocols/drug_design/ReactionBasedAnalogSampler.fwd.hh>
#include <protocols/chemistries/Chemistry.hh>
#include <protocols/drug_design/ChemicalReaction.hh>

#include <rdkit/DataStructs/ExplicitBitVect.h>
#include <rdkit/GraphMol/Fingerprints/MorganFingerprints.h>
#include <rdkit/GraphMol/Fingerprints/Fingerprints.h>

#include <core/chemical/rdkit/RDKit.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <numeric/random/WeightedSampler.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <list>

namespace protocols {
namespace drug_design {

class ReactionBasedAnalogSampler : public protocols::chemistries::Chemistry {
	struct Reagent
	{
		::RDKit::ROMolOP rdmol_;
		std::shared_ptr<ExplicitBitVect> fp_;
		std::string rxn_;
		core::Size no_;		// component number

		bool operator<( const Reagent& r ) const {	// for sorting components before running the reaction
			return no_ < r.no_;
		}
	};

	struct Product
	{
		::RDKit::RWMolOP rdmol_;
		std::string smiles_;
		std::shared_ptr<ExplicitBitVect> fp_;
		utility::vector1< core::Size > frags_;

		bool operator==( const Product& p ) const {
			return smiles_ == p.smiles_;
		}

		struct HashFunction
		{
			size_t operator()( const Product& p ) const {
				return std::hash<std::string>()( p.smiles_ );
			}
		};
	};

public:
	ReactionBasedAnalogSampler();
	virtual ~ReactionBasedAnalogSampler();

	void apply( core::chemical::MutableResidueType & ) override;

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & datacache
	) override;

	core::chemical::VDVDMapping
	get_mapping() const override;

	void reset_spl_ratio() { geo_spl_ratio_ = dynamic_spl_ratios_[1]; }

	/// @brief load reaction files from a designated path (deprecated)
	void load_reactions( std::string const & reaction_dir, std::string const & filename );

	/// @brief load reactions from a single file
	void load_reactions( std::string const & filename );

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	/// @brief load reagents through corresponding reactions (deprecated)
	void load_all_reagents();

	/// @brief load reagents from a single file
	void load_all_reagents( std::string const & filename );

	/// @brief check if the reagent contains bad structures that would later fail RDKit
	bool check_reagent_validity( ::RDKit::ROMOL_SPTR reag ) const;

	/// @brief Draw samples given a list of fragment similarity to reference input
	utility::vector1< Product > sample( std::list< std::pair< core::Real, core::Size > > & score_idx_set ) const;
	/// @brief Draw samples given lists of fragment similarity to first round product reagents
	utility::vector1< Product > sample_fragment( utility::vector1< std::list< std::pair< core::Real, core::Size > > > & score_idx_set ) const;
	/// @brief Draw a product from the given candidate set
	Product sample_candidate( const utility::vector1< Product > & candidate_set, const ::RDKit::RWMOL_SPTR rdmol );

	/// @brief A helper function to easily switch between using ECFP and FCFP
	std::shared_ptr<::RDKit::SparseIntVect<unsigned int>> getMorganFingerprint( ::RDKit::ROMol const & mol, bool useFeatures = false ) const;

	/// @brief Find a common substructure mapping using RDKit Open3DAlign method
	core::chemical::IndexIndexMapping find_O3A_mapping( ::RDKit::ROMOL_SPTR from, ::RDKit::ROMOL_SPTR to ) const;

	/// @brief Setup weights for a sampler. Weighted sampling with geometric distribution; common ratio dynamically decided by population size
	void geometric_sampling( numeric::random::WeightedSampler & sampler, core::Size N );

	/// @brief Return analogs of the chosen product, can be itself
	utility::vector1< Product > analog_search( const ReactionBasedAnalogSampler::Product & prod ) const;

	/// @brief Find all pair of reagents that can undergo the specific reaction from the current sets.
	void pair( core::Size r_no, utility::vector1< utility::vector1< core::Size > > const & sets, utility::vector1< Product > & candidates, int & revisits ) const;
	/// @brief Helper function for pairing. Use a DFS to generate all pairs.
	void pair( core::Size curr_no, core::Size r_no, utility::vector1< core::Size > single_pair, utility::vector1< utility::vector1< core::Size > > & all_pairs, utility::vector1< utility::vector1< core::Size > > const & sets ) const;

	/// @brief Run the reaction given the list of reagent indices
	Product run_reaction( utility::vector1< core::Size > const & reagents ) const;

	/// @brief helper function for comparing score_index paires
	static bool sortbySim( const std::pair< core::Real, core::Size >& a, const std::pair< core::Real, core::Size >& b );

private:

	// Products that have been visited in previous MC loops
	std::unordered_set< Product, Product::HashFunction > visited_;

	// history parameter which records the input structure from last MC cycle
	std::string last_smiles_;

	// The reactions objects
	std::unordered_map< std::string, ChemicalReactionOP > rxns_;
	utility::vector1< Reagent > reagents_;
	numeric::random::WeightedSampler reagent_sampler_;

	// parameter for geometric sampling that mimics the Pareto principle. Ex. ratio = 0.2 means first 20% of the fragments accounts for 80% of the weights.
	// Should not exceeds 0.25 at any point so that the weight still sums to 1
	core::Real geo_spl_ratio_;
	// parameters for dynamic sampling, in the order: min, max, step, OFF_after_n_step, base
	bool dynamic_sampling_;
	utility::vector1< core::Real > dynamic_spl_ratios_;

	// The minimum number of candidates required in the sampling set before the output is drawn
	// This should be at least 4 times of the number of top samples in geometric sampling to ensure validity of the assumption that weights sum to 1.
	core::Size minCandidates_;

	core::chemical::VDVDMapping mapping_;

};

} // namespace drug_design
} // namespace protocols

#endif
