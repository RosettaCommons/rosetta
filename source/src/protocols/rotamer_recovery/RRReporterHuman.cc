// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRReporterHuman.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// James Thompson's code.

/// and apps::pilot::doug::rotamer_prediction_benchmark()


// Unit Headers
#include <protocols/rotamer_recovery/RRReporter.hh>
#include <protocols/rotamer_recovery/RRReporterHuman.hh>


// Project Headers
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/statistics/functions.hh>

// ObjecxxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//C++ Headers
#include <ostream>
#include <string>
#include <map>

using std::endl;
using std::ostream;
using std::string;
using std::map;
using basic::Tracer;
using basic::datacache::CacheableString;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pack::dunbrack::RotVector;
using core::pack::dunbrack::rotamer_from_chi;
using core::pose::Pose;
using core::pose::tag_from_pose;
using utility::vector1;
using numeric::statistics::mean;
using numeric::statistics::std_dev_with_provided_mean;

using ObjexxFCL::string_of;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;

namespace protocols {
namespace rotamer_recovery {

static THREAD_LOCAL Tracer TR("protocol.rotamer_recovery.RRReporterHuman");


PerNativeRRReporterHuman::PerNativeRRReporterHuman() :
	initialized_( false ),
	native_pose_(),
	nat_bb_bins_(),
	nat_rots_(),
	nat_chis_(),
	res_scores_(),
	res_recovered_()
{}

PerNativeRRReporterHuman::PerNativeRRReporterHuman(
	Pose const & native_pose
) :
	initialized_( false ),
	native_pose_(),
	nat_bb_bins_(),
	nat_rots_(),
	nat_chis_(),
	res_scores_(),
	res_recovered_()
{
	set_native( native_pose );
}

PerNativeRRReporterHuman::~PerNativeRRReporterHuman() = default;

PerNativeRRReporterHuman::PerNativeRRReporterHuman(
	PerNativeRRReporterHuman const & ) = default;


/// @detail returns
/// 'O' <- cis-omega
/// 'G' <- alpha-L
/// 'E' <- E
/// 'A' <- helical
/// 'B' <- extended
/// 'X' <- otherwise
char
PerNativeRRReporterHuman::torsion2big_bin(
	Real const phi,
	Real const psi,
	Real const omega
) {
	if ( std::abs( omega ) < 90 ) {
		return 'O'; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100 < psi && psi <= 100 ) {
			return 'G'; // alpha-L
		} else {
			return 'E'; // E
		}
	} else {
		if ( -125 < psi && psi <= 50 ) {
			return 'A'; // helical
		} else {
			return 'B'; // extended
		}
	}
	return 'X';
}


void
PerNativeRRReporterHuman::set_native(
	Pose const & native_pose) {
	if ( !native_pose_  && !initialized_ ) {
		native_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new Pose( native_pose ) ) ); // deep copy please
	} else {
		TR.Fatal << "Attempting to set the native pose in PerNativeRRReporterHuman after it has already been initialized." << endl;
		utility_exit();
	}

	for ( Size ii = 1; ii <= native_pose.size(); ++ii ) {
		Residue const & res( native_pose.residue(ii) );
		nat_bb_bins_.push_back(
			torsion2big_bin(
			res.mainchain_torsion(1),
			res.mainchain_torsion(2),
			res.mainchain_torsion(3)
			)
		);

		RotVector rot;
		rotamer_from_chi( res, rot );
		nat_rots_.push_back( rot );

		//vector1< Real > chi_vec;
		core::pack::dunbrack::RotVector chi_vec;

		for ( Size jj = 1; jj <= res.nchi(); ++jj ) {
			chi_vec.push_back( core::Size( native_pose.chi(jj,ii) ) );
		}
		nat_chis_.push_back( chi_vec );

	}

	res_scores_.resize( native_pose.size() );
	res_recovered_.resize( native_pose.size() );

	initialized_ = true;
}

void
PerNativeRRReporterHuman::report_rotamer_recovery(
	Pose const &,
	Residue const & res,
	Real const score,
	bool const recovered)
{
	if ( !initialized_ ) {
		utility_exit_with_message("Attempting to report rotamer recovery when the native has not be set");
	}

	res_scores_[ res.seqpos() ].push_back( score );
	res_recovered_[ res.seqpos() ].push_back( recovered );
}

bool
PerNativeRRReporterHuman::initialized(
) const {
	return initialized_;
}

void
PerNativeRRReporterHuman::show(
	ostream & out,
	Size const column_width,
	Size const precision
) const {
	out << "#Structure: " << tag_from_pose( *native_pose_ ) << endl;

	for ( Size ii=1; ii <= native_pose_->size(); ++ii ) {
		if ( res_scores_[ii].size() == 0 ) continue;

		Real mean_score = mean(
			res_scores_[ii].begin(), res_scores_[ii].end(), Real(0) );

		Real std_dev_score = std_dev_with_provided_mean(
			res_scores_[ii].begin(), res_scores_[ii].end(), mean_score );

		out << A(column_width,string_of(ii));
		out << A(column_width,nat_bb_bins_[ii]);
		for ( Size i=1; i <= 4; ++i ) {
			if ( i <= nat_rots_[ii].size() ) {
				out << F(column_width,
					precision,
					static_cast< long double >(nat_rots_[ii][i]));
			} else {
				out << A(column_width, "");
			}
		}
		out << F(column_width,precision,mean_score);
		out << F(column_width,precision,std_dev_score);
		out << endl;
	}
}


RRReporterHuman::RRReporterHuman() :
	protocol_name_(),
	protocol_params_(),
	comparer_name_(),
	comparer_params_(),
	column_width_(12),
	precision_(4),
	per_native_recovery_(),
	residues_considered_(0),
	rotamers_recovered_(0),
	recovery_score_mean_(0),
	recovery_score_m2_(0)
{}

RRReporterHuman::RRReporterHuman( RRReporterHuman const & src ) :
	RRReporter(),
	protocol_name_(src.protocol_name_),
	protocol_params_(src.protocol_params_),
	comparer_name_(src.comparer_name_),
	comparer_params_(src.comparer_params_),
	column_width_(src.column_width_),
	precision_(src.precision_),
	per_native_recovery_(src.per_native_recovery_),
	residues_considered_(src.residues_considered_),
	rotamers_recovered_(src.rotamers_recovered_),
	recovery_score_mean_(src.recovery_score_mean_),
	recovery_score_m2_(src.recovery_score_m2_)
{}

RRReporterHuman::~RRReporterHuman() = default;

void
RRReporterHuman::set_protocol_info(
	string const & protocol_name,
	string const & protocol_params
){
	protocol_name_ = protocol_name;
	protocol_params_ = protocol_params;
}

void
RRReporterHuman::set_comparer_info(
	string const & comparer_name,
	string const & comparer_params
){
	comparer_name_ = comparer_name;
	comparer_params_ = comparer_params;
}

void
RRReporterHuman::reset_recovery(){

	per_native_recovery_.clear();
	residues_considered_=0;
	rotamers_recovered_=0;
	recovery_score_mean_=0;
	recovery_score_m2_=0;
}

void
RRReporterHuman::report_rotamer_recovery(
	Pose const & pose1,
	Pose const & pose2,
	Residue const & /*res1*/,
	Residue const & res2,
	Real const score,
	bool const recovered
) {
	residues_considered_++;
	rotamers_recovered_ += recovered;

	// Update mean and M2 with online algorithm
	// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
	Real delta = score - recovery_score_mean_;
	recovery_score_mean_ += delta/residues_considered_;
	recovery_score_m2_ += delta*(score - recovery_score_mean_); // uses new mean

	string pname = tag_from_pose( pose1 );
	PerNativeRRReporterHuman & native_recovery( per_native_recovery_[ pname ] );
	if ( ! native_recovery.initialized() ) {
		native_recovery.set_native( pose1 );
	}
	native_recovery.report_rotamer_recovery( pose2, res2, score, recovered );

}

void
RRReporterHuman::write_header( ostream & out ) const {

	out << "#Rotamer Recovery Benchmark" << endl;
	out << "#Protocol: " << protocol_name_ << endl;
	out << "#Protocol Paremeters: " << protocol_params_ << endl;
	out << "#Comparer: " << comparer_name_ << endl;
	out << "#Comparer Paremeters: " << comparer_params_ << endl;
	out << "#Number of native structures: " << per_native_recovery_.size() << endl;
	out << "#Number of residues in decoys: " << residues_considered_ << endl;

	Real recovery_rate =   static_cast< Real >(rotamers_recovered_) / static_cast< Real >(residues_considered_);

	out << "#Recovery rate: " << recovery_rate << endl;
	out << "#" << comparer_name_ << " Recovery score mean: " << recovery_score_mean_ << endl;

	Real recovery_score_sample_variance = recovery_score_m2_ / static_cast< Real >(residues_considered_);

	out << "#" << comparer_name_ << " Recovery score sample variance: " << recovery_score_sample_variance << endl;

	out << endl << endl;

	out << A( column_width_, "resi_idx" );
	out << A( column_width_, "nat_bb_bin" );
	out << A( column_width_, "nat_rot1" );
	out << A( column_width_, "nat_rot2" );
	out << A( column_width_, "nat_rot3" );
	out << A( column_width_, "nat_rot4" );
	out << A( column_width_, "E[score]" );
	out << A( column_width_, "SD[score]" );
	out << endl;

}

void
RRReporterHuman::show( ostream & out ) const {

	write_header( out );
	out << endl;

	for ( auto const & nat_it : per_native_recovery_ ) {
		nat_it.second.show( out );
		out << endl;
	}

}

void
RRReporterHuman::show( ) const {
	TR << "Recovered " << rotamers_recovered_ << " rotamers"
		<< " at " << residues_considered_ << " residues"
		<< " for a recovery rate of " << recovery_rate() << "." << endl;
}

Real
RRReporterHuman::recovery_rate() const {
	return Real(rotamers_recovered_) / Real(residues_considered_);
}


} // rotamer_recovery
} // protocols

