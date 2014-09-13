// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/DisulfideEntropyFilter.cc
/// @brief Filter on the entropic effect of disulfide linkage
/// @author Gabriel Rocklin (grocklin@gmail.com)

/// Computes the change in entropy of the denatured state caused by forcing certain disulfide
/// bonds, according to random flight configurational statistics (Gaussian approximation).
///
/// See "Analysis and Classification of Disulphide Connectivity in Proteins: The Entropic
/// effect of Cross-Linkage", PM Harrison & MJE Sternberg, J Mol Biol 1994 244, 448-463

//Unit Headers
#include <protocols/simple_filters/DisulfideEntropyFilter.hh>
#include <protocols/simple_filters/DisulfideEntropyFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>
//Project Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
//Utility Headers
#include <utility/string_util.hh>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
 
namespace bnu = boost::numeric::ublas;


namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_filters.DisulfideEntropyFilter" );

protocols::filters::FilterOP
DisulfideEntropyFilterCreator::create_filter() const { return new DisulfideEntropyFilter; }

std::string
DisulfideEntropyFilterCreator::keyname() const { return "DisulfideEntropy"; }


//default ctor
DisulfideEntropyFilter::DisulfideEntropyFilter() :
	protocols::filters::Filter( "DisulfideEntropy" ),
	tightness_(0.0),
	lower_bound_(0.0)
{}

DisulfideEntropyFilter::DisulfideEntropyFilter(core::Real tightness, core::Real lower_bound) :
	protocols::filters::Filter( "DisulfideEntropy" ),
	tightness_(tightness),
	lower_bound_(lower_bound)
{}

DisulfideEntropyFilter::DisulfideEntropyFilter(
	DisulfideEntropyFilter const & src
) :
	protocols::filters::Filter( "DisulfideEntropy" ),
	tightness_(src.tightness_),
	lower_bound_(src.lower_bound_)

{}

DisulfideEntropyFilter::~DisulfideEntropyFilter() {}

filters::FilterOP
DisulfideEntropyFilter::clone() const {
	return new DisulfideEntropyFilter( *this );
}

filters::FilterOP
DisulfideEntropyFilter::fresh_instance() const {
	return new DisulfideEntropyFilter();
}

core::Real
DisulfideEntropyFilter::lower_bound() const {
	return lower_bound_;
}

void
DisulfideEntropyFilter::lower_bound(
	core::Real value
) {
	lower_bound_ = value;
}

core::Real
DisulfideEntropyFilter::tightness() const {
	return tightness_;
}

void
DisulfideEntropyFilter::tightness(
	core::Real value
) {
	tightness_ = value;
}




void
DisulfideEntropyFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {

	if(tag->hasOption("tightness")){
		tightness(tag->getOption< core::Real >("tightness",0));
		TR << "setting tightness" << std::endl;
	}


	if(tag->hasOption("lower_bound")){
		lower_bound(tag->getOption< core::Real >("lower_bound",0));
		TR << "setting lower_bound" << std::endl;
	}

}


bool
DisulfideEntropyFilter::apply(
	core::pose::Pose const & pose
) const {
	
	//calculate the path length
	core::Real entropy = compute(pose);

	//must be below the maximum path length
	if (entropy < lower_bound_) {
		TR << "Failed disulfide entropy filter (current: " << entropy << " RT/mol | user inputted maximum " << lower_bound_ << ")" << std::endl;
		return false;
	}

	//must also be below the path tightness.
	//calculate number of disulfides
	int n_disulfides = 0;
	for ( core::Size i=1; i != pose.total_residue(); ++i) {
		for ( core::Size j=i + 2; j < pose.total_residue() + 1; ++j) {
			if (pose.residue(i).is_bonded(pose.residue(j)) ){
				n_disulfides++;
			}
		}
	}

	core::Real expected_entropy = ((0.1604 * pose.total_residue()) + (1.7245 * n_disulfides) + 5.1477);


	if (entropy < (expected_entropy - tightness_) ) {
		TR << "Failed disulfide entropy filter (current: " << entropy << " RT/mol | expected base ";
		TR << expected_entropy << " - tightness " << tightness_ << ")" << std::endl;
		return false;
	}

	TR << "Disulfide entropy filter success (current: " << entropy << " RT/mol | expected base ";
	TR << expected_entropy << " - tightness " << tightness_ << ")" << std::endl;

	return true;
}

void
DisulfideEntropyFilter::report(
	std::ostream & out,
	core::pose::Pose const & pose
) const {
	out << "Entropic effect of disulfide bonds: " << compute( pose ) << " RT/mol" << std::endl;
}

core::Real
DisulfideEntropyFilter::report_sm(
	core::pose::Pose const & pose
) const {
	return compute( pose );
}

int determinant_sign(const bnu::permutation_matrix<std::size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}
 
core::Real determinant( bnu::matrix<core::Real> const & m ) {
    bnu::permutation_matrix<std::size_t> pm(m.size1());
    core::Real det = 1.0;
    bnu::matrix<core::Real> m_copy = m;
    if( bnu::lu_factorize(m_copy,pm) ) {
        det = 0.0;
    } else {
        for(core::Size i = 0; i < m_copy.size1(); i++) 
            det *= m_copy(i,i); // multiply by elements on diagonal
        det = det * determinant_sign( pm );
    }
    return det;
}

core::Real
DisulfideEntropyFilter::compute_residual(
    core::pose::Pose const & pose
) const {
    core::Real entropy = compute(pose);

    //must also be below the path tightness.
    //calculate number of disulfides
    int n_disulfides = 0;
    for ( core::Size i=1; i != pose.total_residue(); ++i) {
        for ( core::Size j=i + 2; j < pose.total_residue() + 1; ++j) {
            if (pose.residue(i).is_bonded(pose.residue(j)) ){
                n_disulfides++;
            }
        }
    }

    core::Real expected_entropy = ((0.1604 * pose.total_residue()) + (1.7245 * n_disulfides) + 5.1477);

    return entropy - expected_entropy;
}

core::Real
DisulfideEntropyFilter::compute(
	core::pose::Pose const & pose
) const {
	
	//determine the current disulfide configuration
	utility::vector1< std::pair<Size,Size> > disulf_config;

	for ( core::Size i=1; i != pose.total_residue(); ++i) {
		for ( core::Size j=i + 2; j < pose.total_residue() + 1; ++j) {
			if (pose.residue(i).is_bonded(pose.residue(j)) ){
				std::pair< Size, Size > temp_pair;
				temp_pair = std::make_pair( i, j );
				disulf_config.push_back( temp_pair );
                //TR << i << "  " << j << std::endl;
			}
		}
	}

    //TR << "disulfide_config size " << disulf_config.size() << std::endl;
	if (disulf_config.size() == 0) {
		return 0;
	}

	bnu::matrix<double> m(disulf_config.size(), disulf_config.size());
    for (Size i = 0; i < m.size1() ; ++i) {
        for (Size j = i; j < m.size2() ; ++j) {
            if (disulf_config[j+1].first > disulf_config[i+1].second) {
            	m(i,j) = 0;
                m(j,i) = 0;
                //TR << "0   ";
            } else {
            	if (disulf_config[i+1].second < disulf_config[j+1].second) {
                    m(i, j) = disulf_config[i+1].second - disulf_config[j+1].first; // fill matrix                    
                } else {
                    m(i, j) = disulf_config[j+1].second - disulf_config[j+1].first;
                }

                m(j, i) = m(i, j);
                //TR << m(i,j) << "  ";
            }
        //TR << std::endl;
        }
    }

    // Calculating Eqs. 6 and 7 from 
    // "Analysis and Classification of Disulphide Connectivity in Proteins: The Entropic
    // effect of Cross-Linkage", PM Harrison & MJE Sternberg, J Mol Biol 1994 244, 448-463
    //
    // Using parameters:
    // deltaV = 29.65 angstrom^3
    // b = 3.8 angstrom
    // n = number of disulfides 
    //
    // Create a prefactor "k" in eq. 6 which takes all the terms raised to the n power
    //
    // k = deltaV ( (3/(2*pi*b^2))^(3/2) ) = 0.17827
    //
    // P_n = k^n det(m)^(-3/2)
    //
    // Return dimensionless entropy ln (P_n) rather than multiplying by R
    //

    //TR << "t1 " << static_cast<float>(0.17827) << std::endl;
    //TR << "t2 " << static_cast<int>(disulf_config.size()) << std::endl;
    //TR << "pow1 " << std::pow(static_cast<float>(0.17827), static_cast<int>(disulf_config.size())) << std::endl;
    //TR << "m " << m << std::endl;
    //TR << "t3 " << static_cast<double>(determinant(m)) << std::endl;
    //TR << "m " << m << std::endl;
    //TR << "t4 " << static_cast<double>(-1.5) << std::endl;
    //TR << "pow2 " << std::pow(static_cast<double>(determinant(m)), static_cast<double>(-1.5)) << std::endl;
    //TR << "-log " << -log(std::pow(static_cast<float>(0.17827), static_cast<int>(disulf_config.size())) * std::pow(static_cast<double>(determinant(m)), static_cast<double>(-1.5))) << std::endl;

    return -log(std::pow(static_cast<float>(0.17827), static_cast<int>(disulf_config.size())) * std::pow(static_cast<double>(determinant(m)), static_cast<double>(-1.5)));


    //
}


} // namespace
} // namespace
