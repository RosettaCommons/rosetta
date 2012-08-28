// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
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
// AUTO-REMOVED #include <protocols/noesy_assign/DistanceScoreMover.hh>

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
#include <core/id/Exceptions.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh> // REQUIRED FOR WINDOWS

// for switching residue type set to centroid
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
// for error output
#include <core/chemical/ResidueType.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

//// C++ headers
// AUTO-REMOVED #include <math.h> //for isnan

// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>



static basic::Tracer tr("protocols.noesy_assign.crosspeaks");

using core::Real;
using namespace core;
using namespace basic;

namespace protocols {
namespace noesy_assign {

CrossPeakList::CrossPeakList() :
	assignments_()
{}

CrossPeakList::~CrossPeakList() {}

Size CrossPeakList::count_assignments() const {
	Size total_size( 0 );
	for ( const_iterator it = begin(); it != end(); ++it ) {
		total_size+=(*it)->assignments().size();
	}
	return total_size;
}

void CrossPeakList::delete_diagonal_peaks() {
	tr.Info << "remove diagonal peaks..." << std::endl;

	for ( iterator it = begin(); it != end(); ) {
		bool delete_peak( false );
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
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
	assignments_ = NULL;
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
			if ( cp->volume() <= 0 ) {
				tr.Warning << "ignored peak: zero or negative peak intensity for peak " << cp->peak_id() << std::endl;
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
  for ( CrossPeaks::const_iterator it = peaks_.begin(); it != peaks_.end(); ++it ) {
		if ( last_peak_id > (*it)->peak_id() ) {
			output_adaptor.set_format_from_peak( **it );
			output_adaptor.write_header( os );
		}
		last_peak_id=(*it)->peak_id();
    output_adaptor.write_peak( os, (*it)->peak_id(), **it );
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
  for ( CrossPeaks::const_iterator it = peaks_.begin(); it != peaks_.end(); ++it ) {
		if ( last_peak_id > (*it)->peak_id() ) {
			output_adaptor.set_format_from_peak( **it );
			std::string filename( (**it).filename() );
			os.close();
			open_stream( os, prefix, filename );
			output_adaptor.write_header( os );
		}
		last_peak_id=(*it)->peak_id();
    output_adaptor.write_peak( os, (*it)->peak_id(), **it );
    os << std::endl;
  }
}

void CrossPeakList::find_assignments() {
	tr.Info << "determine initial assignments..." << std::endl;

  for ( CrossPeaks::const_iterator it = peaks_.begin(); it != peaks_.end(); ++it ) {
    (*it)->find_assignments();
  }
	assignments_ = NULL;
}

void CrossPeakList::update_assignment_list() {
	assignments_ = new PeakAssignmentResidueMap();
	assignments_->add( *this );

	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( params.network_use_all_covalent_atoms_ ) { //mainly for backwards compatibility in network analysis.  put all atoms in list for cov. analysis
		if ( peaks().begin() != peaks().end() ) {
			assignments_->add_all_atoms( (*peaks().begin())->resonances() );
		}
	}
}

void CrossPeakList::update_peak_volumina() {
	PROF_START( NOESY_ASSIGN_UPDATE_PEAK_VOL );
	tr.Info << "update peak volumina..."  << std::endl;
	for ( CrossPeakList::iterator it = begin(); it != end(); ++it ) {
    Real sum_volume( 0.0 );
    for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			sum_volume+=(*ait)->peak_volume();
    }
		(*it)->set_cumulative_peak_volume( sum_volume );
	}
	PROF_STOP( NOESY_ASSIGN_UPDATE_PEAK_VOL );
}

void CrossPeakList::update_chemshiftscore() {
	tr.Info << "compute chemical shift score..." << std::endl;
	for ( CrossPeakList::iterator it = begin(); it != end(); ++it ) {
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->update_chemshiftscore_from_peak();
		}
	}
}

void CrossPeakList::update_upperdistance_score() {
	tr.Info << "compute uppder distance (covalent) score " << std::endl;
	for ( CrossPeakList::iterator it = begin(); it != end(); ++it ) {
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->update_upperdistance_score();
		}
	}
}



void CrossPeakList::set_trivial_decoy_compatibility_score() {
	tr.Info << "set trivial 1/n decoy compatibility" << std::endl;
	for ( CrossPeakList::iterator it = begin(); it != end(); ++it ) {
		Real invn = 1.0/(*it)->n_assigned();
		for ( CrossPeak::iterator ait = (*it)->begin(); ait != (*it)->end(); ++ait ) {
			(*ait)->set_decoy_compatibility( invn );
		}
	}
}


void CrossPeakList::update_symmetry_score() {
	tr.Info << " symmetry score " << std::endl;
	if ( !assignments_ ) update_assignment_list();
	runtime_assert( assignments_ );
	assignments().check_for_symmetric_peaks( *this );
}

void CrossPeakList::network_analysis() {
	tr.Info << " network analysis ... " << std::endl;
	if ( !assignments_ ) update_assignment_list();
	runtime_assert( assignments_ );
	Size n_assignments( count_assignments() );
	update_peak_volumina();
	assignments_->network_analysis( n_assignments );
	update_peak_volumina();

 	assignments_->network_analysis( n_assignments );
 	update_peak_volumina();
 	assignments_->network_analysis( n_assignments );
 	update_peak_volumina();
}

void CrossPeakList::eliminate_spurious_peaks() {
	for ( CrossPeakList::iterator it = begin(); it != end(); ++it ) {
		(*it)->eliminated( true /*recompute*/ );
	}
}

Real CrossPeakList::calibrate( PeakCalibrator const& calibrator ) {
	Real average_dist( 0.0 );
	Size ct( 0 );
	for ( CrossPeakList::iterator it = begin(); it != end(); ++it ) {
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
		 core::Size min_seq_separation
) const {
	//count for normalization:
	core::Size ct( 0 );
	for ( CrossPeakList::const_iterator it = begin(); it != end(); ++it ) {
		if ( (*it)->eliminated() ) continue;
		if ( (*it)->min_seq_separation_residue_assignment( 0.1 ) < min_seq_separation ) continue; //ignore peaks that have confident intra-residue assignment
		//		if ( !(*it)->has_inter_residue_assignment( resonances(), 0.1 ) ) continue; //ignore peaks that have confident intra-residue assignment
		++ct;
	}

	for ( CrossPeakList::const_iterator it = begin(); it != end(); ++it ) {
		if ( (*it)->eliminated() ) continue;
		if ( (*it)->min_seq_separation_residue_assignment( 0.1 ) < min_seq_separation ) continue; //ignore peaks that have confident intra-residue assignment
		try {
			core::scoring::constraints::ConstraintOP fa_cst;
			core::scoring::constraints::ConstraintOP cen_cst;
			(*it)->create_fa_and_cen_constraint( fa_cst, cen_cst, pose, centroid_pose, ct );
			runtime_assert( fa_cst && cen_cst );
			fa_set->add_constraint( fa_cst );
			cen_set->add_constraint( cen_cst );
		} catch ( core::id::EXCN_AtomNotFound& excn ) {
			tr.Error << "failed to generate " << "constraint for peak: " << (**it) << std::endl;
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


core::scoring::constraints::ConstraintSetOP CrossPeakList::generate_constraints( core::pose::Pose const& pose, bool centroid, core::Size min_seq_separation ) const {
	using namespace core::scoring::constraints;
	ConstraintSetOP cstset = new ConstraintSet;
	core::pose::Pose centroid_pose;

	if ( tr.Debug.visible() ) pose.dump_pdb("pose_in_generate_constraints.pdb");

	if ( centroid ) {
		centroid_pose = pose;
		core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID );
		if ( tr.Debug.visible() ) {
			centroid_pose.dump_pdb( "centroid_pose.pdb" );
		}
	}
	//count for normalization:
	core::Size ct( 0 );
	for ( CrossPeakList::const_iterator it = begin(); it != end(); ++it ) {
		if ( (*it)->eliminated() ) continue;
		if ( (*it)->min_seq_separation_residue_assignment( 0.1 ) < min_seq_separation ) continue; //ignore peaks that have confident intra-residue assignment
		//		if ( !(*it)->has_inter_residue_assignment( resonances(), 0.1 ) ) continue; //ignore peaks that have confident intra-residue assignment
		++ct;
	}


	for ( CrossPeakList::const_iterator it = begin(); it != end(); ++it ) {
		if ( (*it)->eliminated() ) continue;
		if ( (*it)->min_seq_separation_residue_assignment( 0.1 ) < min_seq_separation ) continue; //ignore peaks that have confident intra-residue assignment
		try {
			if ( !centroid ) cstset->add_constraint( (*it)->create_constraint( pose, ct ) );
			else cstset->add_constraint( (*it)->create_centroid_constraint( pose, centroid_pose, ct ) );
		} catch ( core::id::EXCN_AtomNotFound& excn ) {
			tr.Error << "failed to generate " << ( centroid ?  "centroid " : "full-atom " ) << "constraint for peak: " << (**it) << std::endl;
			tr.Error  << excn << std::endl;
			core::pose::Pose const* current_pose = ( centroid ? &centroid_pose : &pose );
			tr.Info << " residue-type in pose: " << current_pose->residue_type( excn.atom().rsd() ).name3() << " " << excn.atom().rsd() << std::endl;
			tr.Debug << " with these atoms ";
			current_pose->residue_type( excn.atom().rsd() ).show_all_atom_names( tr.Debug );
		}
	}
	return cstset;
}

}
}
