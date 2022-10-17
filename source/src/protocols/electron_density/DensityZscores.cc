// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/electron_density/DensityZscores.cc
/// @brief protocol to score local density-fit
/// @author Gabriella Reggiano (reggiano@uw.edu)

// Unit headers
#include <protocols/electron_density/DensityZscores.hh>
#include <protocols/electron_density/DensityZscoresCreator.hh>

//Package headers
#include <core/scoring/DensityZscoresStatsSetup.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// C++ headers
#include <string>
#include <cmath>

static basic::Tracer TR( "protocols.electron_density.DensityZscores" );

namespace protocols {
namespace electron_density {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
DensityZscores::DensityZscores():
	protocols::moves::Mover( DensityZscores::mover_name() ){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
DensityZscores::~DensityZscores(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief get the average bfactor for heavy atoms
core::Real
DensityZscores::calc_avg_residue_bfac(
	core::Size const resi,
	core::pose::Pose const & pose ) const
{
	core::Real avg_bfac = 0.0;
	core::Size count = 0;
	for ( core::Size atmi=1; atmi <= pose.residue(resi).nheavyatoms(); atmi ++ ) {
		if ( pose.residue(resi).atom_type(atmi).is_virtual() ) continue;
		// we have problems where the OXT artificially lowers the average bfactor
		std::string atm_name_str = pose.residue(resi).atom_name(atmi);
		if ( atm_name_str == " OXT" ) continue;
		avg_bfac += pose.pdb_info()->temperature(resi, atmi);
		count++;
	}
	avg_bfac /= count;
	return avg_bfac;
}

/// @brief get the average bfactor for a neighborhood
core::Real
DensityZscores::calc_avg_nbrhood_bfac(
	core::pose::Pose const & pose,
	core::Vector const x) const
{
	core::Real avg_bfac = 0.0;
	core::Size count = 0;
	for ( core::Size resi=1; resi <= pose.total_residue(); resi++ ) {
		if ( pose.residue(resi).is_virtual_residue() ) continue;
		core::Vector nbr_x = pose.residue(resi).nbr_atom_xyz();
		if ( (x - nbr_x).length() < NBRHD_DIST_ ) {
			for ( core::Size atmi=1; atmi<=pose.residue(resi).nheavyatoms(); atmi++ ) {
				if ( pose.residue(resi).atom_type(atmi).is_virtual() ) continue;
				// we have problems where the OXT artificially lowers the average bfactor
				std::string atm_name_str = pose.residue(resi).atom_name(atmi);
				if ( atm_name_str == " OXT" ) continue;
				avg_bfac += pose.pdb_info()->temperature( resi, atmi );
				count++;
			}
		}
	}
	avg_bfac /= count;
	return avg_bfac;
}

/// @brief Apply the mover
void
DensityZscores::apply( core::pose::Pose& pose ){

	core::Size nres = pose.size();

	res_bfacs_.resize( nres, 0.0 );
	nbrhood_bfacs_.resize( nres, 0.0 );
	win1_denscc_.resize( nres, 0.0 );
	win1_dens_zscore_.resize( nres, 0.0 );
	win3_denscc_.resize( nres, 0.0 );
	win3_dens_zscore_.resize( nres, 0.0 );

	core::scoring::electron_density::ElectronDensity & edm = core::scoring::electron_density::getDensityMap();
	edm.setScoreWindowContext( true );
	edm.setWindow( LG_WIN_SIZE_ );

	// score the pose
	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( core::scoring::elec_dens_window, 1.0 );

	if ( pose.is_fullatom() ) {
		myscore->set_weight( core::scoring::fa_rep, 10e-30 );
	} else {
		myscore->set_weight( core::scoring::vdw, 10e-30 );
	}

	(*myscore)(pose);

	// set up the database reader
	core::scoring::ScoringManager* manager( core::scoring::ScoringManager::get_instance() );
	core::scoring::DensityZscoresStatsSetup const & zscore_stats_setup = manager->get_edens_stats_table();

	core::conformation::symmetry::SymmetryInfoCOP symminfo = nullptr;

	for ( core::Size resi=1; resi <= nres; resi++ ) {
		if ( !pose.residue(resi).is_protein() ) continue;

		std::string resn = pose.residue(resi).name3();

		nbrhood_bfacs_[resi] = calc_avg_nbrhood_bfac(pose, pose.residue(resi).nbr_atom_xyz());
		res_bfacs_[resi] = calc_avg_residue_bfac(resi, pose);

		// collect raw density cross correlation then calc zscore
		edm.setWindow( SM_WIN_SIZE_ );
		win1_denscc_[resi] = edm.matchRes( resi, pose.residue(resi), pose, symminfo, false, false);
		win1_dens_zscore_[resi] = zscore_stats_setup.eval_zscore_by_resn_bfactor( win1_denscc_[resi], resn, nbrhood_bfacs_[resi], SM_WIN_SIZE_ );

		edm.setWindow( LG_WIN_SIZE_ );
		win3_denscc_[resi] = edm.matchRes(resi, pose.residue(resi), pose, symminfo, false, false);
		win3_dens_zscore_[resi] = zscore_stats_setup.eval_zscore_by_resn_bfactor( win3_denscc_[resi], resn, nbrhood_bfacs_[resi], LG_WIN_SIZE_ );

	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
DensityZscores::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
DensityZscores::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void DensityZscores::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "protocol to score local density-fit", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
DensityZscores::fresh_instance() const
{
	return utility::pointer::make_shared< DensityZscores >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
DensityZscores::clone() const
{
	return utility::pointer::make_shared< DensityZscores >( *this );
}

std::string DensityZscores::get_name() const {
	return mover_name();
}

std::string DensityZscores::mover_name() {
	return "DensityZscores";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
DensityZscoresCreator::create_mover() const
{
	return utility::pointer::make_shared< DensityZscores >();
}

std::string
DensityZscoresCreator::keyname() const
{
	return DensityZscores::mover_name();
}

void DensityZscoresCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DensityZscores::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Gabriella Reggiano as its author.
void
DensityZscores::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"DensityZscores", basic::citation_manager::CitedModuleType::Mover,
		"Gabriella Reggiano",
		"University of Washington",
		"reggiano@uw.edu",
		"Wrote the DensityZscores."
		)
	);
}


std::ostream &
operator<<( std::ostream & os, DensityZscores const & mover )
{
	mover.show(os);
	return os;
}


} //electron_density
} //protocols
