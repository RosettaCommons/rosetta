// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/NMerPSSMEnergyFilter.cc
/// @brief
/// @author Chris King (chrisk1@uw.edu)


//Unit Headers
#include <protocols/simple_filters/NMerPSSMEnergyFilter.hh>
#include <protocols/simple_filters/NMerPSSMEnergyFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/NMerPSSMEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/exit.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/format.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/excn/Exceptions.hh>
#include <core/pose/selection.hh>

namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

static thread_local basic::Tracer TR( "protocols.simple_filters.NMerPSSMEnergyFilter" );

protocols::filters::FilterOP
NMerPSSMEnergyFilterCreator::create_filter() const { return protocols::filters::FilterOP( new NMerPSSMEnergyFilter ); }

std::string
NMerPSSMEnergyFilterCreator::keyname() const { return "NMerPSSMEnergy"; }

//default ctor
NMerPSSMEnergyFilter::NMerPSSMEnergyFilter() :
	protocols::filters::Filter( "NMerPSSMEnergy" )
{}

//full ctor, default ctor defined in header file
//TODO: allow setting energy_method_ funxns w/ input params (like in parse ctor)
NMerPSSMEnergyFilter::NMerPSSMEnergyFilter(
	core::Real const score_type_threshold,
	std::string const string_resnums
) :
	protocols::filters::Filter( "NMerPSSMEnergy" )
{
	score_type_threshold_ = score_type_threshold;
	string_resnums_ = string_resnums;
}

NMerPSSMEnergyFilter::~NMerPSSMEnergyFilter() {}

void
NMerPSSMEnergyFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & /*data*/, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	if ( ! tag->hasOption( "threshold" ) ) throw utility::excn::EXCN_RosettaScriptsOption("Must specify 'threshold' for NMerPSSMEnergyFilter.");
	score_type_threshold_ = tag->getOption< core::Real >( "threshold" );

	if ( tag->hasOption( "pssm_fname" ) ) energy_method_.read_nmer_pssm( tag->getOption< std::string >( "pssm_fname" ) );
	if ( tag->hasOption( "pssm_list_fname" ) ) energy_method_.read_nmer_pssm_list( tag->getOption< std::string >( "pssm_list_fname" ) );
	//default blank string --> empty res_set_vec --> incl all residues
	string_resnums_ = tag->getOption< std::string >( "resnums", "" );// these are kept in memory until the pose is available (at apply time)
	energy_method_.nmer_length( tag->getOption< core::Size >( "nmer_length", 9 ) );
	energy_method_.nmer_pssm_scorecut( tag->getOption< core::Real >( "nmer_pssm_scorecut", 0.0 ) );
	energy_method_.gate_pssm_scores( tag->getOption< bool >( "gate_pssm_scores", false ) );
}

bool
NMerPSSMEnergyFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	TR << "NMerPSSM score is " << score << ". ";
	if ( score <= score_type_threshold_ ) {
		TR<<"passing." << std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
NMerPSSMEnergyFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"NMerPSSM score of " << compute( pose )<<'\n';
}

core::Real
NMerPSSMEnergyFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
NMerPSSMEnergyFilter::compute_residue(
	core::pose::Pose const & pose,
	core::Size const seqpos
) const {
	assert( seqpos <= pose.total_residue() );
	//TR<< "Calculating nmer_pssm energies at seqpos: " << seqpos << std::endl;
	using namespace core::scoring;
	EnergyMap emap; //we need to zero out the emap each time!
	energy_method_.residue_energy( pose.residue( seqpos ), pose, emap );
	return( emap[ nmer_pssm ] );
}

core::Real
NMerPSSMEnergyFilter::compute(
	core::pose::Pose const & pose
) const {

	utility::vector1< core::Size > const res_set_vec( core::pose::get_resnum_list_ordered( string_resnums_, pose ) );

	// sum over all res pos
	core::Real score( 0. );
	if ( res_set_vec.size() > 0 ) {
		for ( core::Size i_res_vec = 1; i_res_vec <= res_set_vec.size(); ++i_res_vec ) {
			Size const seqpos( res_set_vec[ i_res_vec ] );
			score += NMerPSSMEnergyFilter::compute_residue( pose, seqpos );
		}
	} else {
		for ( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ) {
			score += NMerPSSMEnergyFilter::compute_residue( pose, seqpos );
		}
	}
	return( score );
}

}
}
