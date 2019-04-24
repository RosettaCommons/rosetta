// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SecretionPredictionFilter
/// @brief Filter for greasy helices  https://www.nature.com/articles/nature06387
/// @author Designed by John Wang(jyjwang@uw.edu) converted to a filter by TJ Brunette(tjbrunette@gmail.com)

// Unit Headers
#include <protocols/simple_filters/SecretionPredictionFilter.hh>
#include <protocols/simple_filters/SecretionPredictionFilterCreator.hh>

// Project Headers
#include <core/pose/chains_util.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/constants.hh>
#include <utility/tag/Tag.hh>

// Parser headers
#include <protocols/filters/Filter.hh>


// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
#include <map>
#include <algorithm>
#include <unordered_map>

static basic::Tracer tr("protocols.filters.SecretionPredictionFilter");

namespace protocols {
namespace simple_filters {

using utility::vector1;
using core::Size;
using core::Real;
using std::set;
using namespace core;



SecretionPredictionFilter::SecretionPredictionFilter():
	Filter( "SecretionPredictionFilter" )
{
	window_size_=19;
	dG_ins_threshold_=2.7;
	keep_ = "DEKRQN";
}


// @brief copy constructor
SecretionPredictionFilter::SecretionPredictionFilter( SecretionPredictionFilter const & rval ):
	Super( rval ),
	threshold_(rval.threshold_),
	window_size_(rval.window_size_),
	dG_ins_threshold_(rval.dG_ins_threshold_),
	keep_(rval.keep_)
{
}

// @brief destructor
SecretionPredictionFilter::~SecretionPredictionFilter() = default;


/// @brief
SecretionPredictionFilter::Real
SecretionPredictionFilter::report_sm( const Pose & pose ) const
{
	return  compute( pose );
}

/// @brief
void
SecretionPredictionFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "SectionScore" <<  compute( pose ) << std::endl;
}





//@brief Build a table, the hard way
double SecretionPredictionFilter::get_coeff(const std::string & str1_const, const std::string & str2_const) const
//had trouble making making unordered_map of unordered_maps const. So I refactored.
{
	std::unordered_map<std::string, std::unordered_map<std::string, double>> coeffs;
	coeffs["A"]["a0"] = 1.27e-01; coeffs["A"]["a1"] = 2.15e-02; coeffs["A"]["a2"] = 0; coeffs["A"]["a3"] = 0; coeffs["A"]["a4"] =0;
	coeffs["C"]["a0"] = -7.65e-02; coeffs["C"]["a1"] = 9.94e-02; coeffs["C"]["a2"] = 0; coeffs["C"]["a3"] = 0; coeffs["C"]["a4"] =0;
	coeffs["D"]["a0"] = 1.79e00; coeffs["D"]["a1"] = 1.73e-02; coeffs["D"]["a2"] = 0; coeffs["D"]["a3"] = 0; coeffs["D"]["a4"] =0;
	coeffs["E"]["a0"] = 1.42e00; coeffs["E"]["a1"] = 8.94e-03; coeffs["E"]["a2"] = 0; coeffs["E"]["a3"] = 0; coeffs["E"]["a4"] =0;
	coeffs["F"]["a0"] = -2.77e-01; coeffs["F"]["a1"] = 1.03e-03; coeffs["F"]["a2"] = 0; coeffs["F"]["a3"] = 0; coeffs["F"]["a4"] =0;
	coeffs["G"]["a0"] = 4.81e-01; coeffs["G"]["a1"] = 4.72e-03; coeffs["G"]["a2"] = 0; coeffs["G"]["a3"] = 0; coeffs["G"]["a4"] =0;
	coeffs["H"]["a0"] = 1.20e00; coeffs["H"]["a1"] = 8.01e-03; coeffs["H"]["a2"] = 0; coeffs["H"]["a3"] = 0; coeffs["H"]["a4"] =0;
	coeffs["I"]["a0"] = -4.60e-01; coeffs["I"]["a1"] = 1.81e-02; coeffs["I"]["a2"] = 0; coeffs["I"]["a3"] = 0; coeffs["I"]["a4"] =0;
	coeffs["K"]["a0"] = 1.85e00; coeffs["K"]["a1"] = 2.18e-02; coeffs["K"]["a2"] = 0; coeffs["K"]["a3"] = 0; coeffs["K"]["a4"] =0;
	coeffs["L"]["a0"] = -4.28e-01; coeffs["L"]["a1"] = 2.38e-03; coeffs["L"]["a2"] = 0; coeffs["L"]["a3"] = 0; coeffs["L"]["a4"] =0;
	coeffs["M"]["a0"] = -7.75e-02; coeffs["M"]["a1"] = 9.84e-02; coeffs["M"]["a2"] = 0; coeffs["M"]["a3"] = 0; coeffs["M"]["a4"] =0;
	coeffs["N"]["a0"] = 1.33e00; coeffs["N"]["a1"] = 9.24e-03; coeffs["N"]["a2"] = 0; coeffs["N"]["a3"] = 0; coeffs["N"]["a4"] =0;
	coeffs["P"]["a0"] = 1.09e00; coeffs["P"]["a1"] = 1.01e-02; coeffs["P"]["a2"] = 0; coeffs["P"]["a3"] = 0; coeffs["P"]["a4"] =0;
	coeffs["Q"]["a0"] = 1.33e00; coeffs["Q"]["a1"] = 1.12e-02; coeffs["Q"]["a2"] = 0; coeffs["Q"]["a3"] = 0; coeffs["Q"]["a4"] =0;
	coeffs["R"]["a0"] = 1.65e00; coeffs["R"]["a1"] = 5.12e-02; coeffs["R"]["a2"] = 0; coeffs["R"]["a3"] = 0; coeffs["R"]["a4"] =0;
	coeffs["S"]["a0"] = 7.02e-01; coeffs["S"]["a1"] = 7.77e-03; coeffs["S"]["a2"] = 0; coeffs["S"]["a3"] = 0; coeffs["S"]["a4"] =0;
	coeffs["T"]["a0"] = 5.27e-01; coeffs["T"]["a1"] = 3.12e-02; coeffs["T"]["a2"] = 0; coeffs["T"]["a3"] = 0; coeffs["T"]["a4"] =0;
	coeffs["V"]["a0"] = -2.45e-01; coeffs["V"]["a1"] = 9.79e-02; coeffs["V"]["a2"] = 0; coeffs["V"]["a3"] = 0; coeffs["V"]["a4"] =0;
	coeffs["W"]["a0"] = 2.91e-01; coeffs["W"]["a1"] = 1.89e-02; coeffs["W"]["a2"] = -5.48e-01; coeffs["W"]["a3"] = 9.30e-02; coeffs["W"]["a4"] =6.47e00;
	coeffs["Y"]["a0"] = 6.28e-01; coeffs["Y"]["a1"] = 1.04e-02; coeffs["Y"]["a2"] = -5.74e-01; coeffs["Y"]["a3"] = 9.48e-02; coeffs["Y"]["a4"] =6.92e00;
	std::string str1=str1_const;
	std::string str2=str2_const;
	return(coeffs[str1][str2]);
}


//@brief Calculate the dG_insertion for a given window
core::Real SecretionPredictionFilter::dG_ins_for_window( const std::string & window ) const
{
	using numeric::constants::d::pi;
	core::Real dG_tot = 0;
	double dG_aa = 0;
	double mu_tot = 0;
	double length_term = 0;
	double k = window.length();
	double mu0 = 0;
	double mu1 = 0;

	for ( unsigned int j = 0; j < window.length(); j++ ) {
		double i = 9.0*((j/(k-1.0))*2.0 - 1.0);
		std::string AA = window.substr(j,1);
		double a0 = get_coeff(AA,"a0");
		double a1 = get_coeff(AA,"a1");
		double a2 = get_coeff(AA,"a2");
		double a3 = get_coeff(AA,"a3");
		double a4 = get_coeff(AA,"a4");
		double dG_aa_i = a0*exp(-a1*i*i) + a2*exp(-(a3*(i-a4)*(i-a4)))+a2*exp(-a3*(i+a4)*(i+a4));
		mu0 = mu0 + dG_aa_i*sin(100.0*(j+1)*pi/180);
		mu1 = mu1 + dG_aa_i*cos(100.0*(j+1)*pi/180);

		dG_aa = dG_aa + dG_aa_i;
	}
	mu_tot = 0.27*sqrt(mu0*mu0 + mu1*mu1);
	length_term = 9.29 - 0.645*k + 0.00822*k*k;
	dG_tot = dG_aa + mu_tot + length_term;

	return dG_tot;
}


std::string
SecretionPredictionFilter::get_window_at_index(const Pose & pose , const core::Size & wlen, const core::Size & start ) const
{
	std::string seq = pose.sequence();
	std::string out = seq.substr(start,wlen);
	return out;
}


//@brief This identifies TM regions in the pose
vector1<SecretionPredictionFilter::tm_region>
SecretionPredictionFilter::find_tm_regions( std::ostream & out,  const Pose & pose  ) const
{
	//it'll return output, which is a subset of all of the tm energies in the pose
	vector1<tm_region> output;
	//energies will have ALL of the tm regions in the sequence
	std::vector<tm_region> energies;
	//handle symmetry by just counting the tm regions in the asu (first x residues)
	//core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	core::Size indep = pose.chain_sequence( core::pose::get_chain_id_from_chain("A",pose) ).length();
	core::Size total_scan = indep - window_size_;
	//i'm not good at c++ so i instantiate a temporary tm_region struct just to make life a little easier
	tm_region temp;
	for ( core::Size q = 0; q < total_scan; q++ ) {
		temp.index = q + 1; //dammit, should have used proper indexing
		temp.sequence = get_window_at_index(pose,window_size_,q);
		temp.dG_ins = dG_ins_for_window(temp.sequence);
		energies.push_back({temp.index, temp.sequence, temp.dG_ins, 0, temp.possible_mutants});
	}
	std::vector<tm_region> energies_sorted = energies;
	std::sort(energies_sorted.begin(), energies_sorted.end());

	int tot = (int) total_scan;
	for ( core::Size q = 0; q < (total_scan)/2; q++ ) {//i just scan the top half of the scoring regions for local minima. there might be a better way to do this.
		int i = energies_sorted[q].index - 1;
		int iplus = i + 1;
		int iminus = i - 1;
		int iplus2 = i + 2;
		int iminus2 = i - 2;
		if ( iplus >= tot ) { iplus = i; }
		if ( iplus2 >= tot ) { iplus2 = i; }
		if ( iminus <= 0 ) { iminus = i; }
		if ( iminus2 <= 0 ) { iminus2 = i; }
		if ( energies[i].dG_ins <= energies[iplus].dG_ins && energies[i].dG_ins <= energies[iminus].dG_ins && energies[i].dG_ins <= dG_ins_threshold_ && energies[i].dG_ins <= energies[iplus2].dG_ins && energies[i].dG_ins <= energies[iminus2].dG_ins && energies[i].dG_ins <= dG_ins_threshold_ ) {
			out << "Found a local minimum at index " << energies[i].index << " with sequence "<< energies[i].sequence <<  " with dG_ins_pred of: " << energies[i].dG_ins << '\n';
			for ( core::Size j = energies[i].index; j <= energies[i].index + window_size_ - 1 ; ++j ) {
				energies[i].possible_mutants.push_back({j, keep_, 0, 0});
			}
			output.push_back({energies[i].index,energies[i].sequence,energies[i].dG_ins, 0, energies[i].possible_mutants});
		}//this is a little complicated. energies and output are vectors of tm_region. tm_region also contains vectors of mutts as their possible_mutants. for now, just fill in keep_ (designable residues) for the possible mutants of each output tm_region.
	}
	return output;
}




/// @brief
Real
SecretionPredictionFilter::compute( const Pose & pose ) const
{
	vector1<SecretionPredictionFilter::tm_region> found_regions = SecretionPredictionFilter::find_tm_regions(tr, pose);
	Real min_score = 99;
	for ( Size ii=1; ii<=found_regions.size(); ++ii ) {
		if ( min_score>found_regions[ii].dG_ins ) {
			min_score=found_regions[ii].dG_ins;
		}
	}
	return(min_score);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SecretionPredictionFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	tr << "value" << value << "filtered_value_" << threshold_ << std::endl;
	if ( value <= threshold_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << threshold_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
SecretionPredictionFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	// set threshold
	threshold_ = tag->getOption<Real>("threshold",-999999);
	window_size_=19;
}

// XRW TEMP filters::FilterOP
// XRW TEMP SecretionPredictionFilterCreator::create_filter() const { return protocols::filters::FilterOP(new SecretionPredictionFilter); }

// XRW TEMP std::string
// XRW TEMP SecretionPredictionFilterCreator::keyname() const { return "SSBisectddGFilter"; }

std::string SecretionPredictionFilter::name() const {
	return class_name();
}

std::string SecretionPredictionFilter::class_name() {
	return "SecretionPredictionFilter";
}

void SecretionPredictionFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "XRW TO DO", "-999999");
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SecretionPredictionFilterCreator::keyname() const {
	return SecretionPredictionFilter::class_name();
}

protocols::filters::FilterOP
SecretionPredictionFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( utility::pointer::make_shared<SecretionPredictionFilter>() );
}

void SecretionPredictionFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SecretionPredictionFilter::provide_xml_schema( xsd );
}

} // filters
} // protocols
