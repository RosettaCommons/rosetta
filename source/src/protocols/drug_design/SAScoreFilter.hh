// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/SAScoreFilter.hh
/// @brief A filter which computes the Novartis SA_Score metric
/// @author Rocco Moretti (rmorettiase@gmail.com)

/// @details
///  calculation of synthetic accessibility score as described in:
///
/// Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
/// Peter Ertl and Ansgar Schuffenhauer
/// Journal of Cheminformatics 1:8 (2009)
/// http://www.jcheminf.com/content/1/1/8

#ifndef INCLUDED_protocols_drug_design_SAScoreFilter_hh
#define INCLUDED_protocols_drug_design_SAScoreFilter_hh

//unit headers
#include <protocols/drug_design/SAScoreFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <string>

#include <utility/SingletonBase.hh>

#include <boost/cstdint.hpp>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/DataStructs/SparseIntVect.h>

namespace protocols {
namespace drug_design {

class SAScoreFilter : public filters::Filter
{
public:
	SAScoreFilter():
		Filter( class_name() )
	{}

	SAScoreFilter( std::string const & residue ):
		Filter( class_name() ),
		residue_( residue )
	{}

	core::Real threshold() const { return threshold_; }

	void threshold(core::Real setting) { threshold_ = setting; }

	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new SAScoreFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new SAScoreFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const &pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	core::Size
	num_spiro(RDKit::ROMol const & mol) const;

	core::Size
	num_bridgeheads(RDKit::ROMol const & mol) const;

	/// @details nonconst mol as it recalculates annotations
	core::Size
	num_chiral(RDKit::ROMol & mol, bool includeUnassigned=false ) const;

	RDKit::SparseIntVect<boost::uint32_t> *
	get_morgan_fingerprint(RDKit::ROMol const & mol, int radius) const;

public: // public for direct access in tests and other movers

	// @details nonconst mol because it calls num_chiral()
	core::Real
	calculate_rdkit(RDKit::ROMol & mol) const;

private:
	/// @brief Which residue to calculate the score for
	std::string residue_;

	/// @brief Threshold (maximum) for truth value context
	core::Real threshold_;

};

///////////////////////////////////////////////////////////////////////////////
class SAScoreData : public utility::SingletonBase< SAScoreData >
{

public:
	friend class utility::SingletonBase< SAScoreData >;

	float
	operator[] ( boost::uint32_t index ) const;

private:
	SAScoreData();
	SAScoreData( SAScoreData const & ) = delete; // unimplemented
	SAScoreData const & operator = ( SAScoreData const & ) = delete; // unimplemented

	/// @brief function that actually loads the data
	void
	load_data_from_file( std::string const & filename );

private: // Data
	// float to save space
	typedef std::map< boost::uint32_t, float > fscores_t;
	/// @brief The data for computation
	fscores_t fscores_;
	// @brief The default value to use.
	float const default_;

};


}
}

#endif
