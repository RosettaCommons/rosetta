// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/HBNetScore.hh
/// @brief header file for HBNet filter
/// @details Calculates HBNet's score used for evaulation and scoring. The score is equal to ( sum of all hbond energies in network / num residues in network ) + penalty for unsatisfaction. The best networks get scores < 0 but there are still some good networks with slightly positive scores ( 0 - 5.0 ) so a threshold of 0 is on the conservative side. More benchmarking needs to be done.
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_protocols_hbnet_HBNetScore_HH
#define INCLUDED_protocols_hbnet_HBNetScore_HH

//#include <utility/pointer/ReferenceCount.hh>
#include <protocols/hbnet/HBNetScore.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace hbnet {

class HBNetScore : public protocols::filters::Filter {

public:

	//constructor
	HBNetScore();
	HBNetScore( protocols::filters::Filter const & );
	HBNetScore( HBNetScore const & );

	//destructor
	~HBNetScore();

	protocols::filters::FilterOP clone() const override{
		return protocols::filters::FilterOP( new HBNetScore( *this ) );
	}

	protocols::filters::FilterOP fresh_instance() const override{
		return protocols::filters::FilterOP( new HBNetScore );
	}

public:
	void report( std::ostream &, core::pose::Pose const & ) const override;

	///@brief just returns get_score()
	inline core::Real report_sm( core::pose::Pose const & pose) const override{
		core::pose::PoseOP clone = pose.clone();
		return get_score( *clone );
	}

	bool apply( core::pose::Pose const & pose ) const override;

	///@brief just returns get_score()
	inline core::Real score( core::pose::Pose & pose) override{
		return get_score( pose );
	}

	core::Real get_score( core::pose::Pose & pose ) const;

	//virtual
	inline std::string name() const override{ return "HBNetScore"; }

	inline void clear() override{
		threshold_ = 0;
		hbond_threshold_ = -0.25;
	}

public:

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

private:
	core::Real threshold_;
	core::Real hbond_threshold_;
};

} //hbnet
} //protocols

#endif
