// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @detailed
///
///
/// @author Oliver Lange



// Unit Headers
#include <protocols/simple_filters/ChiWellRmsdEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>
#include <protocols/jd2/util.hh>
#include <protocols/loops/loops_main.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <basic/MetricValue.hh>
// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/evaluation/util.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

//Auto Headers
static basic::Tracer tr("protocols.evaluation.RMSD");
namespace protocols {
namespace simple_filters {
using namespace core;


ChiWellRmsdEvaluator::ChiWellRmsdEvaluator( core::pose::PoseCOP pose, core::Size nchi_max,  core::Real sasa_threshold, utility::vector1< Size> const& selection, std::string tag )
  : evaluation::SingleValuePoseEvaluator< Real >( tag ),
		rmsd_pose_( pose ),
		nchi_max_( nchi_max ),
		sasa_threshold_( sasa_threshold ),
		tag_( tag )
{
	copy( selection.begin(), selection.end(), std::back_inserter( selection_ ) );
}

ChiWellRmsdEvaluator::ChiWellRmsdEvaluator( core::pose::PoseCOP pose, core::Size nchi_max, core::Real sasa_threshold, std::string tag )
	: evaluation::SingleValuePoseEvaluator< Real >( tag ),
		rmsd_pose_( pose ),
		nchi_max_( nchi_max ),
		sasa_threshold_( sasa_threshold ),
		tag_( tag )
{
	if ( pose ) evaluation::find_existing_residues( pose, tag, selection_ );
}

Real
ChiWellRmsdEvaluator::apply( core::pose::Pose& pose ) const {
	core::pose::PoseCOP target_pose = rmsd_pose_;
	if ( !target_pose ) {
		runtime_assert( jd2::jd2_used() );
		target_pose = jd2::get_current_jobs_starting_pose();
	}
	if ( !target_pose ) utility_exit_with_message(" no target pose for rmsd simple_filters "+tag_ );

  core::Size correct( 0 );
	core::Size  total( 0 );

	core::pose::metrics::simple_calculators::SasaCalculatorLegacy sasa_calc;
	basic::MetricValue<utility::vector1< Real > > residue_sasa;

	sasa_calc.get( "residue_sasa", residue_sasa, pose );


	for (	core::scoring::ResidueSelection::const_iterator it = selection_.begin(); it != selection_.end(); ++it ) {
		Size seqpos( *it );
		if ( 	residue_sasa.value()[ seqpos ] > sasa_threshold_ ) {
			tr.Debug << "residue " << seqpos << " is solvent exposed (SASA: " << residue_sasa.value()[ seqpos ] << " ) ignored... " << std::endl;
			continue;
		}
		conformation::Residue const & rsd( pose.residue(seqpos) );
		conformation::Residue const & target_rsd( target_pose->residue(seqpos) );
		if ( rsd.type().name() != target_rsd.type().name() ) continue; // residue types must match

		pack::dunbrack::RotVector rot;
		pack::dunbrack::rotamer_from_chi( rsd, rot );

		pack::dunbrack::RotVector target_rot;
		pack::dunbrack::rotamer_from_chi( target_rsd, target_rot );

		Size nchi = rsd.chi().size();
		bool good( true );
		for ( Size i=1; i<= std::min( nchi, nchi_max_ ) && good; ++i ) {
			tr.Debug << "Chi " << i << " " << rot[ i ]  << " target: " << target_rot[ i ] << std::endl;
			good = rot[ i ] == target_rot[ i ];
		}
		if ( good ) correct+=1;
		total += 1;
		tr.Debug << "pos: " << seqpos << " correct: " << correct << " total: " << total << std::endl;
	}
	return 1.0/total*correct;
}

}
}
