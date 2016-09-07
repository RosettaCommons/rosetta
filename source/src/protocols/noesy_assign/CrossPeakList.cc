// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/CrossPeakList.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
#include <protocols/noesy_assign/PeakCalibrator.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/PeakFileFormat.hh>

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/Pose.hh>

#include <core/id/Exceptions.hh>

// for switching residue type set to centroid
#include <core/chemical/ChemicalManager.fwd.hh>
// for error output
#include <core/chemical/ResidueType.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

//// C++ headers
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

using core::Real;
using namespace core;
using namespace basic;

namespace protocols {
namespace noesy_assign {

CrossPeakList::CrossPeakList() :
	peaks_()
{}

CrossPeakList::~CrossPeakList() = default;

Size CrossPeakList::count_assignments() const {
	Size total_size( 0 );
#ifndef WIN32
	for (const auto & it : *this) {
		total_size+=it->assignments().size();
	}
#endif
	return total_size;
}

void CrossPeakList::delete_diagonal_peaks() {
	tr.Info << "remove diagonal peaks..." << std::endl;
#ifndef WIN32

	for ( auto it = begin(); it != end(); ) {
		bool delete_peak( false );
		for ( auto ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			if ( (*ait)->resonance_id( 1 )==(*ait)->resonance_id( 2 ) ) {
				delete_peak = true;
				break;
			}
		}
		if ( delete_peak ) {
			tr.Debug << "remove peak because it might be a diagonal peak " << (*it)->peak_id() << std::endl;
			it = peaks_.erase( it );
			continue;
		}
		++it;
	}
	assignments_ = nullptr;
#endif
}

void CrossPeakList::read_from_stream( std::istream& is, PeakFileFormat& input_adaptor, ResonanceListOP resonances  ) {
	std::string new_peak_line ="";
	while ( is.good() ) {
		input_adaptor.read_header( is, new_peak_line );
		input_adaptor.output_diagnosis( tr.Debug );
		CrossPeakOP cp;
		if ( !is.good() ) {
			utility_exit_with_message( "nothing to read after peak-file header" );
		}
		while ( is.good() && ( cp = input_adaptor.read_peak( is, new_peak_line )) ) {
			cp->set_resonances( resonances );
			input_adaptor.write_peak( tr.Debug, cp->peak_id(), *cp );
			tr.Debug << std::endl;
			if ( std::abs( cp->volume() ) < input_adaptor.minimum_peak_intensity() ) {
				cp->set_volume( input_adaptor.minimum_peak_intensity() );
			}
			if ( std::abs( cp->volume() ) < 0.01 ) {
				tr.Warning << "ignored peak: zero intensity for peak " << cp->peak_id() << std::endl;
				continue;
			}
			if ( cp->volume() < 0 && input_adaptor.ignore_negative_intensity() ) {
				tr.Warning << "ignored peak [#IGNORE_NEGATIVE_INTENSITY]: negative intensity for peak " << cp->peak_id() << std::endl;
				continue;
			}
			peaks_.push_back( cp );
		}
	}

	if ( is.fail() && is.eof() && is.bad() ) {
		tr.Error << "[ERROR WHILE READING]" << std::endl;
	}
}

void CrossPeakList::write_to_stream( std::ostream& os, PeakFileFormat& output_adaptor  ) const {
	output_adaptor.set_format_from_peak( **peaks_.begin() );
	output_adaptor.write_header( os );
	Size last_peak_id( 0 );
	for (const auto & peak : peaks_) {
		if ( last_peak_id > peak->peak_id() ) {
			output_adaptor.set_format_from_peak( *peak );
			output_adaptor.write_header( os );
		}
		last_peak_id=peak->peak_id();
		output_adaptor.write_peak( os, peak->peak_id(), *peak );
		os << std::endl;
	}
}

void open_stream( utility::io::ozstream& os, std::string const& prefix, std::string const& peak_filename )  {
	utility::file::FileName pre_fn( prefix );
	std::string filename( pre_fn.base() + "_" +peak_filename +".peaks" );
	os.open(filename);
}

void CrossPeakList::write_peak_files( std::string const& prefix, PeakFileFormat& output_adaptor ) const {
	output_adaptor.set_format_from_peak( **peaks_.begin() );
	std::string filename( (**peaks_.begin()).filename() );

	utility::io::ozstream os;
	open_stream( os, prefix, filename );

	output_adaptor.write_header( os );
	Size last_peak_id( 0 );
	for (const auto & peak : peaks_) {
		if ( last_peak_id > peak->peak_id() ) {
			output_adaptor.set_format_from_peak( *peak );
			std::string filename( (*peak).filename() );
			os.close();
			open_stream( os, prefix, filename );
			output_adaptor.write_header( os );
		}
		last_peak_id=peak->peak_id();
		output_adaptor.write_peak( os, peak->peak_id(), *peak );
		os << std::endl;
	}
}

void CrossPeakList::find_assignments() {
	tr.Info << "determine initial assignments..." << std::endl;

	for ( CrossPeaks::const_iterator it = peaks_.begin(); it != peaks_.end(); ++it ) {
		(*it)->find_assignments();
	}
	assignments_ = nullptr;
}

void CrossPeakList::update_assignment_list() {
	assignments_ = PeakAssignmentResidueMapOP( new PeakAssignmentResidueMap() );
	assignments_->add( *this );

	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( params.network_use_all_covalent_atoms_ ) { //mainly for backwards compatibility in network analysis.  put all atoms in list for cov. analysis
		if ( peaks().begin() != peaks().end() ) {
			assignments_->add_all_atoms( (*peaks().begin())->resonances() );
		}
	}
}

void CrossPeakList::update_peak_volumina() {
#ifndef WIN32
	PROF_START( NOESY_ASSIGN_UPDATE_PEAK_VOL );
	tr.Info << "update peak volumina..."  << std::endl;
	for ( auto it = begin(); it != end(); ++it ) {
		Real sum_volume( 0.0 );
		for ( auto ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			sum_volume+=(*ait)->peak_volume();
		}
		(*it)->set_cumulative_peak_volume( sum_volume );
	}
	PROF_STOP( NOESY_ASSIGN_UPDATE_PEAK_VOL );
#endif
}

void CrossPeakList::update_chemshiftscore() {
#ifndef WIN32
	tr.Info << "compute chemical shift score..." << std::endl;
	for ( auto it = begin(); it != end(); ++it ) {
		for ( auto ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->update_chemshiftscore_from_peak();
		}
	}
#endif
}

void CrossPeakList::update_upperdistance_score() {
#ifndef WIN32
	tr.Info << "compute uppder distance (covalent) score " << std::endl;
	for ( auto it = begin(); it != end(); ++it ) {
		for ( auto ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->update_upperdistance_score();
		}
	}
#endif
}


void CrossPeakList::set_trivial_decoy_compatibility_score() {
#ifndef WIN32
	tr.Info << "set trivial 1/n decoy compatibility" << std::endl;
	for ( auto it = begin(); it != end(); ++it ) {
		Real invn = 1.0/(*it)->n_assigned();
		for ( auto ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->set_decoy_compatibility( invn );
		}
	}
#endif
}


void CrossPeakList::update_symmetry_score() {
	tr.Info << " symmetry score " << std::endl;
	if ( !assignments_ ) update_assignment_list();
	runtime_assert( assignments_ != nullptr );
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	Real const min_sym_cont( params.min_contribution_symmetric_peaks_ );
	assignments().check_for_symmetric_peaks( *this, min_sym_cont < 0.99 );
}

void CrossPeakList::network_analysis() { //ResonanceList const& resonances ) {
	tr.Info << " network analysis ... " << std::endl;
	if ( !assignments_ ) update_assignment_list();
	runtime_assert( assignments_ != nullptr );
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( params.network_mode_ == "orig" ) {
		Size n_assignments( count_assignments() );
		assignments_->network_analysis( n_assignments );
	} else if ( params.network_mode_ == "clean" ) {
		assignments_->network_analysis2();
	} else {
		utility_exit_with_message(" network mode " + params.network_mode_ + " is unknown " );
	}
}

void CrossPeakList::eliminate_spurious_peaks() {
	for ( auto it = begin(); it != end(); ++it ) {
		(*it)->eliminated( true /*recompute*/ );
	}
}

Real CrossPeakList::calibrate( PeakCalibrator const& calibrator ) {
	Real average_dist( 0.0 );
	Size ct( 0 );
	for ( auto it = begin(); it != end(); ++it ) {
		PeakCalibrator::TypeCumulator calibration_types;
		(*it)->calibrate( calibrator, calibration_types );
		average_dist += (*it)->distance_bound();
		ct += ( (*it)->distance_bound() > 0.0 );
	}
	if ( ct == 0 ) return 0.0;
	return average_dist / ct;
}

void CrossPeakList::generate_fa_and_cen_constraints(
	core::scoring::constraints::ConstraintSetOP fa_set,
	core::scoring::constraints::ConstraintSetOP cen_set,
	core::pose::Pose const& pose,
	core::pose::Pose const& centroid_pose,
	core::Size min_seq_separation,
	core::Size min_quali,
	core::Size max_quali,
	core::Real padding,
	bool ignore_elimination_candidates, /*default = true */
	bool elimination_candidates /* default = false */
) const {
	//count for normalization:
	core::Size ct( 0 );
	for (const auto & it : *this) {
		if ( it->eliminated() ) continue;
		if ( it->min_seq_separation_residue_assignment( 0.1 ) < min_seq_separation ) continue; //ignore peaks that have confident intra-residue assignment
		//  if ( !(*it)->has_inter_residue_assignment( resonances(), 0.1 ) ) continue; //ignore peaks that have confident intra-residue assignment
		if ( max_quali+1 < CrossPeak::MAX_CLASS ) {
			Size const quality( it->quality_class() );
			if ( quality < min_quali || quality > max_quali ) continue;
		}
		++ct;
	}
	for (const auto & it : *this) {
		if ( it->eliminated() ) continue;
		if ( it->min_seq_separation_residue_assignment( 0.1 ) < min_seq_separation ) continue; //ignore peaks that have confident intra-residue assignment
		if ( max_quali+1 < CrossPeak::MAX_CLASS ) {
			Size const quality( it->quality_class() );
			if ( quality < min_quali || quality > max_quali ) continue;
		}
		if ( !ignore_elimination_candidates ) {
			if ( it->is_elimination_candidate() != elimination_candidates ) continue;
		}
		try {
			core::scoring::constraints::ConstraintOP fa_cst;
			core::scoring::constraints::ConstraintOP cen_cst;
			it->create_fa_and_cen_constraint( fa_cst, cen_cst, pose, centroid_pose, ct, padding );
			runtime_assert( fa_cst && cen_cst );
			fa_set->add_constraint( fa_cst );
			cen_set->add_constraint( cen_cst );
		} catch ( core::id::EXCN_AtomNotFound& excn ) {
			tr.Error << "failed to generate " << "constraint for peak: " << (*it) << std::endl;
			tr.Error  << excn << std::endl;
			tr.Info << " residue-type in pose: " << pose.residue_type( excn.atom().rsd() ).name3() << " " << excn.atom().rsd() << std::endl;
			tr.Debug << " with these atoms ";
			pose.residue_type( excn.atom().rsd() ).show_all_atom_names( tr.Debug );

			tr.Info << " residue-type in centroid_pose: " << centroid_pose.residue_type( excn.atom().rsd() ).name3() << " " << excn.atom().rsd() << std::endl;
			tr.Debug << " with these atoms ";
			centroid_pose.residue_type( excn.atom().rsd() ).show_all_atom_names( tr.Debug );
		}
	}
}

}
}
