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
#include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/ResonanceList.hh>

// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/PeakCalibrator.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
// Utility headers

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/exit.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
// #include <core/options/keys/abinitio.OptionKeys.gen.hh>
// #include <core/options/keys/run.OptionKeys.gen.hh>
//#include <core/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <iostream>
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>

//Auto Headers
#include <cmath>

static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

using core::Real;
using namespace core;
using namespace basic;


namespace protocols {
namespace noesy_assign {

CrossPeak::Spin::Spin( Real freq ) : freq_ ( freq ) {}
CrossPeak::Spin::Spin() { }
CrossPeak::Spin::~Spin() {}
core::Size CrossPeak::Spin::assignment_index( core::Size assignment ) const {
	core::Size ct( 1 );
	for ( SpinAssignments::const_iterator it = assignments_.begin(); it != assignments_.end(); ++it ) {
		if ( assignment == *it ) return ct;
		++ct;
	}
	return 0;
}

CrossPeak::CrossPeak( Spin const& sp1, Spin const& sp2, Real strength ) :
	proton1_( sp1 ),
	proton2_( sp2 ),
	volume_( strength ),
	cumulative_peak_volume_( 1.0 ),
	distance_bound_( 100 ), //what would be a good initial value? -1 ?
	eliminated_( NOT_ELIMINATED ),
	eliminated_due_to_dist_violations_( false ),
	elimination_candidate_( false )
{}

CrossPeak::CrossPeak() :
	cumulative_peak_volume_( 1.0 ),
	distance_bound_( 100 ),
	eliminated_( NOT_ELIMINATED ),
	eliminated_due_to_dist_violations_( false ),
	elimination_candidate_( false )
{}

CrossPeak::~CrossPeak() {}

/// @detail use for reading of assignments from file
/// pass as spin1, spin2, label1, label2 (indices 1..4)
void CrossPeak::add_full_assignment( Size res_ids[] ) {
	//  std::cerr << "CrossPeak::add_full_assignment stubbed out " << res_ids[ 1 ] << std::endl;
	//   std::cerr << "run with 'ignore_assignments'" << std::endl;
	//  utility_exit_with_message( "stubbed function" );
	Size ind1 = assign_spin( 1, res_ids );
	Size ind2 = assign_spin( 2, res_ids );
#ifndef WIN32
	assignments_.push_back( protocols::noesy_assign::PeakAssignmentOP( new PeakAssignment( this, ind1, ind2 ) ) );
#endif
	// for ( Size i = 1;
	//need to find resonances in Spins and add them if they are still missing.
	//get id for spin1 and spin2
	//then make PeakAssignment() and add to assignments_
}

bool CrossPeak::has_proton( core::Size select ) const {
	return info( select ).proton_tolerance() < 99;
}

/// @brief find all possible assignments based on chemical shifts and tolerances
void CrossPeak::find_assignments( ) {

	if ( proton1_.n_assigned() && proton2_.n_assigned() ) return; //if assignments are already present do nothing

	runtime_assert( proton1_.n_assigned() == proton2_.n_assigned() );
	assign_spin( 1 );
	assign_spin( 2 );

	Size const n_assigned_1( proton( 1 ).n_assigned() );
	Size const n_assigned_2( proton( 2 ).n_assigned() );

	for ( Size ct1 = 1; ct1 <= n_assigned_1; ++ct1 ) {
		for ( Size ct2 = 1; ct2 <= n_assigned_2; ++ct2 ) {
#ifndef WIN32
			assignments_.push_back( protocols::noesy_assign::PeakAssignmentOP( new PeakAssignment( this, ct1, ct2 ) ) );
#endif
		}
	}
	//write_to_stream( tr.Debug );
	//  tr.Debug << std::endl;
}

void CrossPeak::print_peak_info( std::ostream& os ) const {
	os << "peak: " << peak_id() << " ";
	for ( Size i = 1; i<=2 ; i++ ) {
		os << info( i ).main_atom();
		if ( info( i ).has_label() ) os << "-" << info( i ).label_atom_type();
		if ( i == 1 ) os << " <-> ";
	}
}

/// @brief assign protons based on chemical shifts and tolerances
void CrossPeak::assign_spin( Size iproton ) {
	//base-class: disregard label
	Real const my_freq( proton( iproton ).freq() );
	Real const my_tolerance( info( iproton ).proton_tolerance() );
	for ( ResonanceList::const_iterator it = resonances().begin(); it != resonances().end(); ++it )  {
		//    if ( std::abs( fold_resonance( it->second->freq(), iproton ) - my_freq ) < std::max( my_tolerance, it->second->error() ) ) {
		if ( it->second->match( my_freq, my_tolerance, folder( iproton ) ) ) {
			proton( iproton ).add_assignment( it->first );
		}
	}
}

/// @brief assign protons ass pre-determined
core::Size CrossPeak::assign_spin( Size iproton, Size res_id[] ) {
	Size ind = proton( iproton ).assignment_index( res_id[ iproton ] );
	if ( ind ) return ind;
	proton( iproton ).add_assignment( res_id[ iproton ] );
	return proton( iproton ).n_assigned();
}

core::Real round( core::Real d, core::Size digits ) {
	for ( Size i=1; i<=digits; ++i ) {
		d*=10;
	}
	d = floor( d + 0.5 );
	for ( Size i=1; i<=digits; ++i ) {
		d/=10;
	}
	return d;
}

void
CrossPeak::create_fa_and_cen_constraint(
	core::scoring::constraints::ConstraintOP& fa_cst,
	core::scoring::constraints::ConstraintOP& cen_cst,
	pose::Pose const& pose,
	pose::Pose const& centroid_pose,
	core::Size normalization,
	core::Real padding,
	bool fa_only
) const {
#ifndef WIN32
	core::Size const round_digits( 2 );
	using namespace core::scoring::constraints;
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	core::Real inv_weight( sqrt( 1.0*normalization )/params.cst_strength_ );

	core::scoring::func::FuncOP func( new BoundFunc(1.5,
		round(distance_bound()+padding,round_digits),
		round(inv_weight,round_digits),
		"automatic NOE Peak "+ObjexxFCL::string_of( peak_id() )+" "+filename()+" Volume: "+ObjexxFCL::string_of( volume() )
		) );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool const map_to_CEN( params.map_to_cen_atom_ );
	// tr.Debug << "MAPPING TO " << ( map_to_CEN ? "CEN-ATOM" : "C-BETA" ) << std::endl;

	QualityClass const qc( quality_class() );
	Real min_vc = params.min_volume_;
	if ( qc == HI_NEAR_UNAMBIG ) {
		min_vc = std::max( min_vc, 0.1 );
	}

	Size ct_ambiguous( 0 );
	PeakAssignments::const_iterator first_valid = end();
	for ( PeakAssignments::const_iterator it = begin(); it != end(); ++it ) {
		if ( (*it)->normalized_peak_volume() <= min_vc ) continue; //not enough contribution
		if ( !ct_ambiguous ) first_valid = it;
		ct_ambiguous += (*it)->float_ambiguity();
	}

	if ( ct_ambiguous > 1 ) {
		/// create Ambiguous constraint with BoundFunc describing the whole thing...
		AmbiguousNMRConstraintOP
			my_fa_cst( new AmbiguousNMRConstraint( func ) );

		AmbiguousNMRConstraintOP
			my_cen_cst( new AmbiguousNMRConstraint( func ) );

		//  core::Size n( n_assigned() );
		//mjo commenting out 'sd' because it is unused and causes a warning
		//core::Real sd( 1.0/params.cst_strength_ );
		core::Real eff_single_dist( pow( pow( distance_bound(), -6 ) / ct_ambiguous, -1.0/6 ) ); //assuming equal contribution of all which is of course wrong
		core::Real cum_new_distance( 0 );
		core::Size max_maps=0;
		// add individual constraints --- their potential does not matter ... it is ignored in evaluation...
		for ( PeakAssignments::const_iterator it = first_valid; it != end(); ++it ) {
			if ( (*it)->normalized_peak_volume() <= min_vc ) continue; //not enough contribution
			core::Size const float_ambiguity( (*it)->float_ambiguity() );
			for ( core::Size ifloat = 1; ifloat <=  float_ambiguity; ++ifloat ) {
				AmbiguousNMRDistanceConstraintOP new_cst( (*it)->create_constraint( pose, ifloat ) );
				my_fa_cst->add_individual_constraint( new_cst );
				fa_cst = my_fa_cst;
				if ( !fa_only ) {
					Size number_of_maps; //are 0, 1 or 2 sidechains replaced by CB
					if ( map_to_CEN ) {
						my_cen_cst->add_individual_constraint( new_cst->map_to_CEN( pose, centroid_pose, number_of_maps, "CEN" ) );
						max_maps = std::max( number_of_maps, max_maps );
					} else {
						my_cen_cst->add_individual_constraint( new_cst->map_to_CEN( pose, centroid_pose, number_of_maps, "CB" ) );
						cum_new_distance += pow( (eff_single_dist + params.centroid_mapping_distance_padding_*number_of_maps )*pow( new_cst->multiplicity(), 1.0/6 ), -6 );
					}
				}
			}
		}
		if ( my_fa_cst->member_constraints().size()<=1 ) {
			fa_cst = my_fa_cst->member_constraints()[ 1 ]->clone( my_fa_cst->get_func().clone() );
		}
		if ( !fa_only ) {
			Real mapped_upl;
			Real mapped_inv_weight;
			if ( map_to_CEN ) {
				mapped_upl = distance_bound() + max_maps; //add 0, 1 or 2 A
				mapped_inv_weight = inv_weight * pow( 2.0, double(max_maps) );
			} else {
				mapped_upl = pow( cum_new_distance, -1.0/6 );
				mapped_inv_weight = inv_weight;
			}
			core::scoring::func::FuncOP centroid_bound_func( new BoundFunc( 1.5,
				round(mapped_upl,round_digits),
				round(mapped_inv_weight,round_digits),
				"CEN mapped automatic NOE: Peak "+ ObjexxFCL::string_of( peak_id() )
				) );
			if ( my_cen_cst->member_constraints().size()<=1 ) {
				cen_cst = my_cen_cst->member_constraints()[ 1 ]->clone( centroid_bound_func );
			} else {
				cen_cst = my_cen_cst->clone( centroid_bound_func );
			}
		}

	} else { // not ambiguous
		if ( first_valid == end() ) return;
		AmbiguousNMRDistanceConstraintOP my_fa_cst = (*first_valid)->create_constraint( pose, 1 /*ifloat*/, func ); //first one should be only one,
		fa_cst = my_fa_cst;
		if ( !fa_only ) {
			Size number_of_maps; //are 0, 1 or 2 sidechains replaced by CB
			core::Size n( n_assigned() );
			core::Real eff_single_dist( pow( pow( distance_bound(), -6 ) / n, -1.0/6 ) ); //assuming equal contribution of all which is of course wrong
			core::Real cum_new_distance( 0 );
			ConstraintOP my_cen_cst;
			core::Real mapped_upl;
			core::Real mapped_inv_weight;
			if ( map_to_CEN ) {
				my_cen_cst= my_fa_cst->map_to_CEN( pose, centroid_pose, number_of_maps, "CEN" );
				mapped_upl = distance_bound() + number_of_maps;
				mapped_inv_weight = inv_weight*pow( 2.0, double(number_of_maps) );
			} else {
				my_cen_cst= my_fa_cst->map_to_CEN( pose, centroid_pose, number_of_maps, "CB" );
				cum_new_distance += pow( (eff_single_dist + params.centroid_mapping_distance_padding_*number_of_maps )*pow( my_fa_cst->multiplicity(), 1.0/6 ), -6 );
				mapped_upl = pow( cum_new_distance, -1.0/6 );
				mapped_inv_weight = inv_weight;
			}
			std::string const comment( "CEN mapped automatic NOE: Peak "+ ObjexxFCL::string_of( peak_id() ) );
			cen_cst = my_cen_cst->clone( core::scoring::func::FuncOP( new BoundFunc( 1.5, round(mapped_upl,round_digits), round(mapped_inv_weight,round_digits), comment ) ) );
		}
		// tr.Trace << "constraint for " << peak_id() << " finished " << std::endl;
	}

#endif //WIN32
}

Real sigmoid( Real x, Real tau, Real m, int sign = 1 ) {
	if ( sign > 0 ) {
		return 1.0/(1+exp(-1.0/tau*5.0*(x-m)));
	}
	return 1.0-1./(1+exp(-1.0/tau*5.0*(x-m)));
}

core::Real CrossPeak::probability() const {
#ifndef WIN32
	Real max_vc_cs( 0.0 );
	Real max_vc_sym( 0.0 );
	Real max_vc( 0.0 );
	Real const overall_vc( cumulative_peak_volume() );
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	for ( const_iterator ait = begin(); ait != end(); ++ait ) {
		if ( max_vc < (*ait)->normalized_peak_volume() ) {
			max_vc = (*ait)->normalized_peak_volume();
			max_vc_sym = (*ait)->symmetry_compliance();
			max_vc_cs = (*ait)->chemshift_compliance();
		}
	}
	utility::vector1< Real > const& w( params.prob_sigmoid_w_ );
	utility::vector1< Real > const& tau( params.prob_sigmoid_tau_ );
	utility::vector1< Real > const& m( params.prob_sigmoid_m_ );
	Real s[6];
	s[1] = sigmoid( max_vc_cs, tau[1], m[1], 1);
	s[2] = sigmoid( max_vc_sym, tau[2], m[2], 1);
	s[3] = -sigmoid( max_vc_sym, tau[3], m[3], -1);
	s[4] = sigmoid( max_vc, tau[4], m[4], 1 );
	s[5] = -sigmoid( overall_vc, tau[5], m[5], -1 );
	return (w[1]*s[1]+w[2]*s[2]+w[3]*s[3]+w[4]*s[4]+w[5]*s[5]+0.4)/1.6;
#else
	return 0.;
#endif
}

Real CrossPeak::smallest_native_violation() const {
#ifndef WIN32
	Real viol( 100000 );
	QualityClass pclass( quality_class() );
	for ( const_iterator ait = begin(); ait != end(); ++ait ) {
		Real const vc( (*ait)->normalized_peak_volume() );
		if ( pclass==HI_UNAMBIG && vc < 0.1 ) continue;
		if ( pclass<=MED_AMBIG && vc < 0.01 ) continue;
		if ( viol > (*ait)->native_distance_viol() ) viol=(*ait)->native_distance_viol();
	}
	return viol;
#else
	return 0.;
#endif
}

std::string CrossPeak::quality_class_str() const {
	char const strings[][50] = {"HI_UNAMBIG", "HI_NEAR_UNAMBIG", "HI_AMBIG", "MED_AMBIG", "MED_UNAMBIG", "LOW_AMBIG" };
	return strings[ quality_class() ];
}

CrossPeak::QualityClass CrossPeak::quality_class() const {
#ifndef WIN32
	Size count_vc_0p1( 0 );
	Size count_vc_0p01( 0 );
	for ( const_iterator ait = begin(); ait != end(); ++ait ) {
		Real const vc( (*ait)->normalized_peak_volume() );
		if ( vc > 0.1 ) ++count_vc_0p1;
		if ( vc > 0.01 ) ++count_vc_0p01;
	}
	enum ProbClass {
		HI = 0,
		MED,
		LOW
	};

	Real const prob( probability() );
	ProbClass pclass;
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	if ( prob > params.prob_level_[1] ) {
		pclass = HI;
	} else if ( prob > params.prob_level_[2] ) {
		pclass = MED;
	} else {
		pclass = LOW;
	}

	if ( pclass == HI ) {
		if ( count_vc_0p01 == 1 ) {
			return HI_UNAMBIG;
		} else if ( count_vc_0p1 == 1 ) {
			return HI_NEAR_UNAMBIG;
		} else {
			return HI_AMBIG;
		}
	} else if ( pclass == MED ) {
		if ( count_vc_0p1 == 1 || count_vc_0p01 == 1 ) {
			return UNAMBIG_MED_PROB;
		}
		return MED_AMBIG;
	} else {
		return BAD_LOW_PROB;
	}
#else
	return BAD_LOW_PROB;
#endif
}

core::Real CrossPeak::max_volume_contribution() const {
#ifndef WIN32
	Real max_volume( 0.0 );
	for ( const_iterator ait = begin(); ait != end(); ++ait ) {
		max_volume= std::max( max_volume, (*ait)->normalized_peak_volume() );
	}
	return max_volume;
#else
	return 0.;
#endif
}

/// @brief do we have a inter residue assignment with at least volume_threshold contribution ?
Size CrossPeak::min_seq_separation_residue_assignment( Real volume_threshold ) const {
#ifndef WIN32
	Size min_seq( 99999 );
	for ( PeakAssignments::const_iterator it = begin(); it != end(); ++it ) {
		if ( (*it)->normalized_peak_volume() <= volume_threshold ) continue; //not enough contribution
		Size res1( (*it)->resid( 1 ) );
		Size res2( (*it)->resid( 2 ) );
		Size diff( res1 < res2 ? res2-res1 : res1-res2 );
		min_seq = min_seq < diff ? min_seq : diff;
	}
	return min_seq;
#else
	return 0;
#endif
}

std::string CrossPeak::elimination_reason() const {
	static std::string const str_distviol( "DistViol" );
	static std::string const str_network( "Network" );
	static std::string const str_minvol( "MinPeakVol");
	static std::string const str_maxassign( "MaxAssign");
	static std::string const empty_str( "" );
	if ( eliminated_ == EL_DISTVIOL ) return str_distviol + " " + elimination_comment_;
	if ( eliminated_ == EL_NETWORK ) return str_network;
	if ( eliminated_ == EL_MINPEAKVOL ) return str_minvol;
	if ( eliminated_ == EL_MAXASSIGN ) return str_maxassign;
	return empty_str;
}

bool CrossPeak::eliminated( bool recompute, bool do_not_compute ) const {
#ifndef WIN32
	if ( recompute && !do_not_compute ) eliminated_ = eliminated_due_to_dist_violations_ ? EL_DISTVIOL : NOT_ELIMINATED;
	//if dist_cut problem eliminated_ should already be set to false.
	// tr.Trace << "elimination check for peak " << peak_id() << "...";
	// if ( eliminated_ ) tr.Trace << "eliminated from cached value" << std::endl;
	if ( eliminated_ ) return true;
	if ( do_not_compute ) return false; //not eliminated as of now...

	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	bool min_peak_volume_reached( false );
	for ( PeakAssignments::const_iterator it = begin(); it != end() && !min_peak_volume_reached ; ++it ) {
		min_peak_volume_reached = (*it)->normalized_peak_volume() > params.min_volume_;
	}
	if ( !min_peak_volume_reached ) {
		eliminated_ = EL_MINPEAKVOL;
		//tr.Trace << "eliminated (min_peak_volume)" << std::endl;
		return true;
	}

	if ( assignments().size() > params.nmax_ ) {
		eliminated_ = EL_MAXASSIGN;
		// tr.Trace << "eliminated (too many assignments)" << std::endl;
		return true;
	}

	//network anchoring
	Real N_atom_sum( 0.0 );
	Real N_res_sum( 0.0 );
	for ( PeakAssignments::const_iterator it = begin(); it != end(); ++it ) {
		Real vol( (*it)->normalized_peak_volume() );
		N_atom_sum +=  vol * (*it)->network_anchoring();
		N_res_sum += vol * (*it)->network_anchoring_per_residue();
	}

	if ( N_res_sum <= params.network_reswise_high_  &&
			( N_res_sum < params.network_reswise_min_ ||
			N_atom_sum < params.network_atom_min_ ) ) {
		eliminated_ = EL_NETWORK;
		// tr.Trace << "eliminated (network)" << std::endl;
		return true;
	}

#endif
	// tr.Trace << "passed" << std::endl;
	return false;
}

core::Size CrossPeak::n_Vmin_assignments() {
#ifndef WIN32
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	Size ct( 0 );
	for ( PeakAssignments::const_iterator it = begin(); it != end(); ++it ) {
		Real vol( (*it)->normalized_peak_volume() );
		ct += vol > params.min_volume_;
	}
	return ct;
#else
	return 0;
#endif
}

void CrossPeak::nudge_distance_bound( core::Real offset ) {
	distance_bound_ += offset;
}

void CrossPeak::calibrate( PeakCalibrator const& calibrator, PeakCalibrator::TypeCumulator& calibration_types ) {
#ifndef WIN32
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	Real sum( 0.0 );
	Size ct( 0 );
	// if ( volume_ <= 0.0 ) throw utility::excn::EXCN_BadInput("Peak intensity negative or zero for "+ObjexxFCL::string_of( peak_id_ ) );
	for ( PeakAssignments::const_iterator it = begin(); it != end(); ++it ) {
		CALIBRATION_ATOM_TYPE type1, type2;
		type1 = (*it)->calibration_atom_type( 1 );
		type2 = (*it)->calibration_atom_type( 2 );
		calibration_types.set( type1 );
		calibration_types.set( type2 );
		Real const cal( sqrt( calibrator( type1 ) * calibrator( type2 ) ) );
		Real vol( (*it)->normalized_peak_volume() );
		Real int_factor(1.0);
		int_factor*=resonances()[ (*it)->resonance_id( 1 ) ].intensity();
		int_factor*=resonances()[ (*it)->resonance_id( 2 ) ].intensity();
		sum += (vol > params.min_volume_ ? vol : 0 ) / cal / int_factor;
		ct += vol > params.min_volume_;
	}
	if ( ct > 0 ) distance_bound_ = pow( sum*std::abs( volume_ ), -1.0/6.0 );
	else distance_bound_ = 0.0;

	core::Real max_dist( info( 1 ).max_noe_distance() );
	if ( max_dist > 0.01 ) {
		distance_bound_ = std::min( distance_bound_, max_dist );
	}
#endif
}

/// @brief assign protons ass pre-determined
Size CrossPeak3D::assign_spin( Size iproton, Size res_id[] ) {
	Size ind = CrossPeak::assign_spin( iproton, res_id );
	if ( iproton == 1  && ind > label( iproton ).n_assigned() ) label( iproton ).add_assignment( res_id[ iproton+2 ] );
	return ind;
}

void CrossPeak3D::assign_spin( Size iproton ) {
	tr.Trace << "assign_spin for Peak " << peak_id() << std::endl;
	if ( iproton == 2 ) CrossPeak::assign_spin( iproton ); //base-class: no label
	else assign_labelled_spin( iproton );
}

void CrossPeak3D::assign_labelled_spin( Size iproton ) {
	runtime_assert( has_label( iproton ));
	Real const proton_freq( proton( iproton ).freq() );
	Real const label_freq( label( iproton ).freq() );
	Real const proton_tolerance( info( iproton ).proton_tolerance() );
	Real const label_tolerance( info( iproton ).label_tolerance() );

	for ( ResonanceList::const_iterator it = resonances().begin(); it != resonances().end(); ++it ) {
		//  tr.Trace << "match2D for peak " << peak_id() << " with resonance " << it->second->atom() << std::endl;
		Resonance::ResonancePairs matches;
		if ( it->second->match2D( proton_freq, proton_tolerance, folder( iproton ), label_freq, label_tolerance, folder( iproton+2 ), matches ) ) {
			tr.Debug << "successful 2D match for peak " << peak_id() << " with resonance " << it->second->atom() << std::endl;
			for ( Resonance::ResonancePairs::const_iterator mit=matches.begin(); mit != matches.end(); ++mit ) {
				proton( iproton ).add_assignment( mit->first );
				label( iproton ).add_assignment( mit->second );
			}
		}
	}
	tr.Trace << "no more matches" << std::endl;
	/*
	/// if we have pseudo 4D spectrum we speed things up a bit by filtering out non-protons here
	if ( my_tolerance > 99 ) {
	for ( ResonanceList::const_iterator it = resonances().begin(); it != resonances().end(); ++it )  {
	//   if ( std::abs( fold_resonance( it->second->freq(), iproton + 2 ) - my_label_freq ) < std::max( my_label_tolerance, it->second->error() ) ) {
	if ( it->second->match( my_label_freq, my_label_tolerance, folder( iproton+2 ) ) ) {
	//now find all proton-resonances that are bound to this label atom
	core::Size resid( it->second->resid() );
	std::string const& label_name( it->second->name() );
	ResonanceList::Resonances const& residue_list( resonances().resonances_at_residue( resid ) );
	for ( ResonanceList::Resonances::const_iterator rit = residue_list.begin(); rit != residue_list.end(); ++rit ) {
	if ( (*rit)->name()[ 0 ]=='Q' || (*rit)->name().find("H") != std::string::npos ) {
	//      tr.Debug << "resid: " << resid<< " test label: " << label_name << " with proton " << (*rit)->name() << std::endl;
	try {
	std::string possible_label( info( iproton ).
	label_atom_name( (*rit)->name(), resonances().aa_from_resid( resid ) ) );
	//tr.Debug << "found possible label " << possible_label << std::endl;
	if ( possible_label == label_name ) {
	label( iproton ).add_assignment( it->first );
	proton( iproton ).add_assignment( (*rit)->label() );
	//we have found a proton that can be attached to our label
	}
	} catch ( EXCN_UnknownAtomname& exception ) {
	continue;
	}
	} // if rit is proton
	} // for rit
	} // if matched resonance
	}// all resonances
	} else {
	for ( ResonanceList::const_iterator it = resonances().begin(); it != resonances().end(); ++it )  {
	if ( it->second->match( my_freq, my_tolerance, folder( iproton ) ) ) {
	Size resid( it->second->atom().rsd() );
	//maybe also map resonance by resid?
	try {
	id::NamedAtomID atomID( info( iproton ).label_atom_name( it->second->atom().atom(), resonances().aa_from_resid( resid ) ), resid );
	Resonance const& label_reso ( resonances()[ atomID ] );
	//     if ( tr_labels.Trace.visible() ) {
	//      tr_labels.Trace << "trying to match " << atomID << " as label to " << it->second->atom() << std::endl;
	//     }
	if ( label_reso.match( my_label_freq, my_label_tolerance, folder( iproton+2 ) ) ) {
	proton( iproton ).add_assignment( it->first );
	label( iproton ).add_assignment( label_reso.label() );
	}
	} catch ( EXCN_UnknownResonance& exception ) {
	if ( !unknown_resonances_.count( exception.atom() ) ) {
	unknown_resonances_.insert( exception.atom() );
	exception.show( tr.Warning );
	tr.Warning << " as label for atom " << it->second->atom().atom() << " " <<  resonances().aa_from_resid( resid ) << std::endl;
	//    if ( tr.Debug.visible() ) exception.show( tr.Debug );
	}
	continue; //if no label is known we don't assign this proton
	} catch ( EXCN_UnknownAtomname& exception ) { //this happens if we try to assign a proton that can't have a label: i.e., a H in a HCH spectrum
	if ( tr_labels.Trace.visible() ) {
	tr_labels.Trace << "cannot find label atom for resid: " + it->second->atom().atom() + " " + ObjexxFCL::string_of( resid ) + " --- ignore proton assignment" << std::endl;
	}
	continue;
	}
	}
	}
	}
	*/
}

CrossPeak3D::CrossPeak3D( Spin const& sp1, Spin const& sp2, Spin const& label1, Real strength ) :
	CrossPeak( sp1, sp2, strength ),
	label1_( label1 )
{}

CrossPeak3D::CrossPeak3D() {}
CrossPeak3D::~CrossPeak3D() {}

CrossPeak4D::CrossPeak4D( Spin const& sp1, Spin const& sp2, Spin const& label1, Spin const& label2, Real strength ) :
	CrossPeak3D( sp1, sp2, label1, strength ),
	label2_( label2 )
{}

CrossPeak4D::CrossPeak4D() {}
CrossPeak4D::~CrossPeak4D() {}

void CrossPeak4D::assign_spin( Size iproton ) {
	assign_labelled_spin( iproton );
}

Size CrossPeak4D::assign_spin( Size iproton, Size res_id[] ) {
	Size ind = CrossPeak::assign_spin( iproton, res_id );
	if ( ind > label( iproton ).n_assigned() ) label( iproton ).add_assignment( res_id[ iproton+2 ] );
	return ind;
}

// void CrossPeak::read_from_stream( std::istream& is ) {
//   Real freq;
//   is >> freq;
//   proton1_ = Spin( freq );
//   if ( has_label( 1 ) ) {
//     is >> freq;
//     label( 1 ) = Spin( freq );
//   }

//   is >> freq;
//   proton2_ = Spin( freq );
//   if ( has_label( 2 ) ) {
//     is >> freq;
//     label( 2 ) = Spin( freq );
//   }

//   Size number;
//   char letter;
//   is >> number >> letter; //no idea what these mean ...it reads now 3 U or 4 U
//   is >> volume_;
//   Real error_of_strength; //... ???
//   is >> error_of_strength;

//   char letter_e;
//   Size number_0;
//   is >> letter_e;
//   is >> number_0;
//   add_assignment_from_stream( is );
// }

// void CrossPeak::add_assignment_from_stream( std::istream& is ) {
//   Size id;
//   is >> id;
//   if ( id )  proton1_.add_assignment( id );
//   if ( has_label( 1 ) ) {
//     is >> id;
//     if ( id ) label( 1 ).add_assignment( id );
//   }

//   is >> id;
//   if ( id ) proton2_.add_assignment( id );
//   if ( has_label( 2 ) ) {
//     is >> id;
//     if ( id ) label( 2 ).add_assignment( id );
//   }
// }


// void CrossPeak::write_to_stream( std::ostream& os ) const {
//   os << ObjexxFCL::format::F( 8, 3, proton1_.freq() ) << " ";
//   if ( has_label( 1 ) ) {
//     os << ObjexxFCL::format::F( 8, 3, label( 1 ).freq() ) << " ";
//   }

//   os << ObjexxFCL::format::F( 8, 3, proton2_.freq() ) << " ";
//   if ( has_label( 2 ) ) {
//     os << ObjexxFCL::format::F( 8, 3, label( 2 ).freq() ) << " ";
//   }

//   os << ObjexxFCL::format::E( 10, 3, strength_ ) << " " << ObjexxFCL::format::E( 10, 3, 0.0 ) << " ";

//   Size assignments_written( 0 );
//   //  while ( assignments_written < proton1_.n_assigned() ) {
//   for ( PeakAssignments::const_iterator it = assignments_.begin(); it != assignments_.end(); ++it ) {
//     ++assignments_written;
//     if ( assignments_written > 1 ) os << std::endl << "                                            ";
//     for ( Size i=1; i<=2; i++ ) {
//       os << ObjexxFCL::format::RJ( 6, (*it)->resonance_id( i ) ) << " ";
//       if ( has_label( i ) ) {
//  os << ObjexxFCL::format::RJ( 6, (*it)->label_resonance_id( i ) << " ";
//       }
//     }
//   }
// }


} //noesy_assign
} //devel
