// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/electron_density/DensityZscores.hh
/// @brief protocol to score local density-fit
/// @author Gabriella Reggiano (reggiano@uw.edu)

#ifndef INCLUDED_protocols_electron_density_DensityZscores_HH
#define INCLUDED_protocols_electron_density_DensityZscores_HH

// Unit headers
#include <protocols/electron_density/DensityZscores.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>
#include <utility/vector1.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

#include <map> // AUTO IWYU For map, allocator

namespace protocols {
namespace electron_density {

///@brief protocol to score local density-fit
class DensityZscores : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	DensityZscores();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~DensityZscores() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	core::Real
	calc_avg_residue_bfac( core::Size const resi,
		core::pose::Pose const & pose ) const;

	core::Real
	calc_avg_nbrhood_bfac( core::pose::Pose const & pose,
		core::Vector const x) const;

	utility::vector1< core::Real >
	get_win1_denscc() const { return win1_denscc_; };

	utility::vector1< core::Real >
	get_win3_denscc()const { return win3_denscc_; };

	utility::vector1< core::Real >
	get_win1_dens_zscore() const { return win1_dens_zscore_; };

	utility::vector1< core::Real >
	get_win3_dens_zscore() const { return win3_dens_zscore_; };

	utility::vector1< core::Real >
	get_nbrhood_bfacs() const { return nbrhood_bfacs_; };

	utility::vector1< core::Real >
	get_res_bfacs() const { return res_bfacs_; };

	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//DensityZscores & operator=( DensityZscores const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Function overrides needed for the citation manager:

	/// @brief This mover is unpublished.  It returns Gabriella Reggiano as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;


private: // data

	core::Real const NBRHD_DIST_ = 8.0;
	core::Size const LG_WIN_SIZE_ = 3;
	core::Size const SM_WIN_SIZE_ = 1;

	// containers for scores
	utility::vector1< core::Real > win1_denscc_;
	utility::vector1< core::Real > win3_denscc_;
	utility::vector1< core::Real > win1_dens_zscore_;
	utility::vector1< core::Real > win3_dens_zscore_;
	utility::vector1< core::Real > nbrhood_bfacs_;
	utility::vector1< core::Real > res_bfacs_;

};

std::ostream &
operator<<( std::ostream & os, DensityZscores const & mover );

} //electron_density
} //protocols

#endif //protocols_electron_density_DensityZscores_HH
