// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file   core/energy_methods/DEPC_MS_Energy.cc
/// @brief  Score term that uses DEPC MS labeling data for Ser, Thr, and Tyr residues to reward based on SASA value, labeling status, and number of hydrophobic neighboring residues and for labeled His and Lys residues to reward based on SASA value. Corroborates with 2021 manuscript by Biehn, Limpikirati, Vachet, and Lindert.
/// @author Sarah Biehn (biehn.4@osu.edu) with help from Sten Heinze (heinze.18@osu.edu)

// Package headers
#include <core/energy_methods/DEPC_MS_Energy.hh>
#include <core/energy_methods/DEPC_MS_EnergyCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/sasa/util.hh>

//Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <basic/prof.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/Tracer.hh>

//Utility headers
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>

static basic::Tracer TR( "core.energy_methods.DEPC_MS_Energy" );

namespace core {
namespace energy_methods {

// In order to speed up code - save residue SASA vector in pose data
struct CacheableSasaVector : public basic::datacache::CacheableData {
	CacheableSasaVector( utility::vector1< core::Real > const & sasa_to_cache ):
		sasa( sasa_to_cache )
	{}

	basic::datacache::CacheableDataOP
	clone() const override {
		return utility::pointer::make_shared< CacheableSasaVector >(*this);
	}

	utility::vector1< core::Real > sasa;

};

typedef utility::pointer::shared_ptr< CacheableSasaVector const > CacheableSasaVectorCOP;

core::scoring::methods::EnergyMethodOP
DEPC_MS_EnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< DEPC_MS_Energy >( options );
}

core::scoring::ScoreTypes
DEPC_MS_EnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back(depc_ms);
	return sts;
}

void
DEPC_MS_Energy::setup_for_scoring(
	pose::Pose & pose, core::scoring::ScoreFunction const &
) const {
	prepare_for_scoring( pose );
}

void
DEPC_MS_Energy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const {
	prepare_for_scoring( pose );
}

void
DEPC_MS_Energy::setup_for_minimizing (
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	kinematics::MinimizerMapBase const &
) const {
	prepare_for_scoring( pose );
}

//constructor
DEPC_MS_Energy::DEPC_MS_Energy( core::scoring::methods::EnergyMethodOptions const & options ):
	parent(utility::pointer::make_shared< DEPC_MS_EnergyCreator >() ),
	depc_ms_input_file_ ( options.depc_ms_input() )
{
	init_from_file();
}

//clone
core::scoring::methods::EnergyMethodOP
DEPC_MS_Energy::clone() const {
	return utility::pointer::make_shared< DEPC_MS_Energy > ( *this );
}

//scoring
void
DEPC_MS_Energy::residue_energy(
	core::conformation::Residue const & residue,
	core::pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const {
	core::Size const number_residues (pose.total_residue());
	core::Real hydrophobic_neighbor_count (0.0);
	static const core::Real labeled_average_hnc (4.42); // Value was optimized in Biehn 2021 DEPC manuscript
	static const core::Real unlabeled_average_hnc (3.39); // Value was optimized in Biehn 2021 DEPC manuscript
	static const core::Real labeled_hk_sasa_midpt (0.650);
	static const core::Real sty_sasa_min (0.045); // Value was optimized in Biehn 2021 DEPC manuscript
	static const core::Real sty_sasa_max (0.355); //Value was optimized in Biehn 2021 DEPC manuscript
	core::Real dist (0.0);

	basic::datacache::CacheableDataCOP dat( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::DEPC_MS_SASA_POSE_INFO ) );
	CacheableSasaVectorCOP cached_sasa = utility::pointer::dynamic_pointer_cast<CacheableSasaVector const>( dat );
	utility::vector1< core::Real > const & sasa = cached_sasa->sasa;

	for ( core::Size j=1; j <= input_res_.size(); j++ ) {
		char res_id = residue.type().name1();
		if ( res_id != 'S' && res_id != 'T' && res_id != 'Y' && res_id != 'H' && res_id !='K' ) { continue; }
		if ( input_res_[j].first == residue.seqpos() ) {
			char label_status = input_res_[j].second;
			switch ( res_id ) {
			case 'S': case 'T': case 'Y' :
				if ( sasa[residue.seqpos()] > sty_sasa_min && sasa[residue.seqpos()] < sty_sasa_max ) { //start looking at STY residues with 5-35% relative SASA
					for ( core::Size res_count_neighbor = 1; res_count_neighbor <= number_residues; res_count_neighbor++ ) {
						char neighbor_res_id = pose.residue(res_count_neighbor).type().name1();
						if ( residue.seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
							std::string neighbor_atom = "CB";
							std::string target_atom;
							switch ( neighbor_res_id ) {
							case 'F': case 'I': case 'W': case 'L': case 'V': case 'M': case 'Y': case 'A': case 'P' :
								if ( res_id == 'S' ) { target_atom = "OG"; } //assign correct atom name of serine hydroxyl oxygen
								else if ( res_id == 'T' ) { target_atom = "OG1"; } //assign correct atom name of threonine hydroxyl oxygen
								else if ( res_id == 'Y' ) { target_atom = "OH"; } //assign correct atom name of tyrosine hydroxyl oxygen
								dist = residue.xyz(target_atom).distance(pose.residue(res_count_neighbor).xyz(neighbor_atom));
								hydrophobic_neighbor_count += 1.0/(1.0 + std::exp(2.0*(dist-8.0)));
								break;
							default : break;
							}
						}
					}
					if ( label_status == 'L' ) { emap[ core::scoring::depc_ms ] += ((-1.0) + 1.0/(1.0 + std::exp(8.0*( hydrophobic_neighbor_count - labeled_average_hnc ))) );
					}
					if ( label_status == 'U' ) { // consider unlabeled residues for unlabeled portion of term
						emap[ core::scoring::depc_ms ] += -1.0/(1.0 + std::exp(8.0*( hydrophobic_neighbor_count - unlabeled_average_hnc )));
					}
				}
				break;
			case 'H': case 'K' :
				if ( label_status == 'L' ) {
					emap[ core::scoring::depc_ms ] += -1.0 + 1.0/(1.0 + std::exp(2.0*( sasa[residue.seqpos()] - labeled_hk_sasa_midpt )));
				}
				break;
			default : break;
			}
		}
	}
}

void
DEPC_MS_Energy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ core::scoring::twelve_A_neighbor_graph ] = false;
}

core::Size
DEPC_MS_Energy::version() const
{
	return 1;
}

void DEPC_MS_Energy::init_from_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::io::izstream input(depc_ms_input_file_);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + depc_ms_input_file_);
		utility_exit_with_message( msg );
	}
	std::string line;

	while ( getline( input, line ) ) {
		if ( line.empty() ) continue;
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream ss(line);
		core::Size resi;
		char label;

		ss >> resi >> label;

		runtime_assert_string_msg( !(ss.fail() || ss.bad() ), "Error in DEPC_MS_Energy::init_from_file(): Could not parse line \"" + line + "\"." );

		if ( label != 'L' && label != 'U' ) { // deliver error message if residue is not specified as labeled/unlabeled
			runtime_assert_string_msg( !(ss.fail() || ss.bad() ), "Error in DEPC_MS_Energy::init_from_file(): Must specify labeled ('L') or unlabeled ('U') - please indicate L or U for each residue in input file.");
		}

		input_res_.push_back( std::pair< core::Size, char >(resi, label) );

	}
}

void DEPC_MS_Energy::prepare_for_scoring( pose::Pose & pose ) const {
	pose.update_residue_neighbors();
	utility::pointer::shared_ptr<CacheableSasaVector> cacheable_sasa_op( utility::pointer::make_shared<CacheableSasaVector >( core::scoring::sasa::rel_per_res_sc_sasa( pose ) ) ); // get relative per residue SASA - SASA vector goes over all residues
	pose.data().set( core::pose::datacache::CacheableDataType::DEPC_MS_SASA_POSE_INFO, cacheable_sasa_op );
}

} //energy_methods
} //core
