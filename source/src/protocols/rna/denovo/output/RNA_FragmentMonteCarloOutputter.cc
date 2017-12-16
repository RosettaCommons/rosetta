// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/output/RNA_FragmentMonteCarloOutputter.cc
/// @brief running output of jumps/scores during rna_denovo/farfar
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/denovo/output/RNA_FragmentMonteCarloOutputter.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/EulerAngles.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/MathNTensor_io.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>

static basic::Tracer TR( "protocols.rna.denovo.output.RNA_FragmentMonteCarloOutputter" );

using namespace core;
using utility::vector1;
using utility::tools::make_vector1;
using utility::tools::make_vector;
using namespace protocols::rna::denovo::options;

//////////////////////////////////////////////////////////////////////
/// @detailed
///
///   Stick everything related to running scores, pose info,
///    6D jump information, etc. here.
///
///                -- rhiju, 2016-2017.
///
//////////////////////////////////////////////////////////////////////

namespace protocols {
namespace rna {
namespace denovo {
namespace output {

//Constructor
RNA_FragmentMonteCarloOutputter::RNA_FragmentMonteCarloOutputter( RNA_FragmentMonteCarloOptionsCOP options,
	core::pose::PoseCOP align_pose ):
	options_(std::move( options ))
{
	initialize( align_pose );
}

//Destructor
RNA_FragmentMonteCarloOutputter::~RNA_FragmentMonteCarloOutputter() = default;

void
RNA_FragmentMonteCarloOutputter::apply( core::pose::Pose & )
{
	// no op. Maybe this should not be a mover.
}

void
RNA_FragmentMonteCarloOutputter::initialize( core::pose::PoseCOP align_pose ) {

	using namespace numeric;
	if ( options_->output_score_frequency() > 0 ) {
		if ( options_->output_score_file() != "none" )  {
			TR << "Opening file for output of running scores every " << options_->output_score_frequency() << " cycles: " << options_->output_score_file() << std::endl;
			running_score_output_.open_append( options_->output_score_file() );
		}
		if ( options_->output_jump_res().size() > 0 ) {
			using namespace core::pose::rna;
			runtime_assert( options_->output_jump_res().size() == 2 );
			output_stub_stub_type_ = BASE_CENTROID;
			if ( options_->output_jump_o3p_to_o5p() ) output_stub_stub_type_ = O3P_TO_O5P;
			if ( options_->output_jump_chainbreak() ) output_stub_stub_type_ = CHAINBREAK;
			if ( options_->output_jump_reference_RT_string().size() > 0 ) {
				reference_RT_ = kinematics::RTOP( new kinematics::RT );
				std::stringstream rt_stream( options_->output_jump_reference_RT_string() );
				kinematics::RT const align_RT = get_output_jump_RT( align_pose );
				TR << TR.Green << "Trying to read in RT from -output_jump_reference_RT. If this fails, you may want to use the following string, based on -align_pdb or -native: \"" << align_RT << "\" " << std::endl;
				rt_stream >> *reference_RT_;
			}
		}
		if ( options_->save_jump_histogram() ) {
			runtime_assert( options_->output_rotation_vector() );
			core::Real const & bxs( options_->jump_histogram_boxsize() );
			core::Real const & bw ( options_->jump_histogram_binwidth() );
			core::Real const & bwr( options_->jump_histogram_binwidth_rotvector() );
			jump_histogram_min_        = make_vector1( -bxs, -bxs, -bxs, -180.0, -180.0, -180.0 );
			jump_histogram_max_        = make_vector1( +bxs, +bxs, +bxs, +180.0, +180.0, +180.0 );
			jump_histogram_bin_width_  = make_vector1(   bw,   bw,   bw,   bwr,    bwr,   bwr );
			vector1< Size > jump_n_bins;
			for ( Size n = 1; n <= 6; n++ ) jump_n_bins.push_back( static_cast<Size>( (jump_histogram_max_[n] - jump_histogram_min_[n]) / jump_histogram_bin_width_[n] + 1.0 ) );
			if ( jump_histogram_ == nullptr ) {
				jump_histogram_ = MathNTensorOP< Size, 6>( new MathNTensor< Size, 6>( jump_n_bins, 0 ) );
			} else {
				for ( Size n = 1; n <= 6; n++ ) runtime_assert( jump_histogram_->n_bins( n ) == jump_n_bins[ n ] );
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOutputter::output_running_info(
	Size const & r,
	Size const & i,
	pose::Pose & pose,
	core::scoring::ScoreFunctionCOP working_denovo_scorefxn )
{
	if ( options_->output_score_frequency() != 0 &&
			i % options_->output_score_frequency() == 0 ) {
		running_score_output_ << r << ' ' << i << " " << ( *working_denovo_scorefxn )( pose );
		output_jump_information( pose );
		running_score_output_ << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////////
core::kinematics::RT
RNA_FragmentMonteCarloOutputter::get_output_jump_RT( pose::PoseCOP pose ) const
{
	using namespace core::kinematics;
	if ( pose == nullptr ) return RT();
	Stub stub1, stub2;
	get_output_jump_stub_stub( *pose, stub1, stub2 );
	Jump const j( stub1, stub2 );
	return j.rt();
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOutputter::get_output_jump_stub_stub(
	pose::Pose const & pose,
	kinematics::Stub & stub1,
	kinematics::Stub & stub2 ) const
{
	using namespace core::conformation;
	Residue const & rsd1 = pose.residue( options_->output_jump_res()[ 1 ] );
	Residue const & rsd2 = pose.residue( options_->output_jump_res()[ 2 ] );
	core::pose::rna::get_stub_stub( rsd1, rsd2, stub1, stub2, output_stub_stub_type_ );
}

//////////////////////////////////////////////////////////////////////////////////
// @details
// allows visualization of base-base rigid body transform distributions -- used for
// defining 'loop_close' energy for RNA, which in turn is important for free energy
// estimation of RNA motifs with unstructured nucleotides (e.g., extrahelical bulges).
void
RNA_FragmentMonteCarloOutputter::output_jump_information( pose::Pose const & pose)
{
	using namespace core::kinematics;
	using namespace core::chemical::rna;
	using namespace core::conformation;

	if ( options_->output_jump_res().size() == 0 ) return;
	Stub stub1, stub2;
	get_output_jump_stub_stub( pose, stub1, stub2 );
	if ( reference_RT_ != nullptr ) reference_RT_->make_jump( stub1 /*start*/, stub1 /*end*/ );

	Jump const j( stub1, stub2 );
	Vector const & t( j.get_translation() );
	vector1< Real > outvals( make_vector1( t.x(), t.y(), t.z() ) );

	if ( options_->output_rotation_vector() ) {
		Vector const rotation_vector( numeric::rotation_axis_angle( j.get_rotation() ) * (180.0 / numeric::constants::r::pi) );
		// magnitude will be angle in degrees, maximum of 180.
		outvals.append( make_vector1( rotation_vector.x(), rotation_vector.y(), rotation_vector.z() ) );
	} else {
		numeric::EulerAngles< Real > euler( j.get_rotation() ); // assuming ZXZ convention!
		outvals.append( make_vector1( euler.phi_degrees(), euler.theta_degrees(), euler.psi_degrees() ) );
	}
	if ( running_score_output_.good() ) {
		for ( auto const & outval : outvals ) running_score_output_ << ' ' << outval;
	}

	if ( options_->save_jump_histogram() ) {
		vector1< Size > outbins;
		for ( Size n = 1; n <= 6; n++ ) {
			// round to *closest* bin  by adding 0.5 before conversion to int.
			auto outbin = static_cast<int>( 0.5 + ( outvals[ n ] - jump_histogram_min_[ n ] ) / jump_histogram_bin_width_[ n ] );
			outbin = std::min( std::max( 0, outbin ), int(jump_histogram_->n_bins( n )) - 1 ); // zero-indexed
			outbins.push_back( Size( outbin ) );
		}
		(*jump_histogram_)( outbins )++;
	}
}


//////////////////////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOutputter::finalize( core::scoring::ScoreFunctionCOP denovo_scorefxn )
{

	if ( options_->output_score_frequency() == 0 ) return;

	if ( running_score_output_.good() ) {
		running_score_output_.close();
		TR << "Created running score file at: " << options_->output_score_file() << std::endl;
	}
	if ( jump_histogram_ != nullptr  ) {
		using namespace utility::json_spirit;
		std::vector< Value > n_bins;
		for ( auto const & v : jump_histogram_->n_bins() ) n_bins.emplace_back(boost::uint64_t(v) );
		std::vector< Value > minval, maxval, binwidth;
		for ( auto const & v : jump_histogram_min_ ) minval.emplace_back(v );
		for ( auto const & v : jump_histogram_max_ ) maxval.emplace_back(v );
		for ( auto const & v : jump_histogram_bin_width_ ) binwidth.emplace_back(v );
		Object json( make_vector( Pair( "n_bins", n_bins ),
			Pair( "type", "uint64" ),
			Pair( "minval",  minval ),
			Pair( "maxval",  maxval ),
			Pair( "binwidth",binwidth ) ) );
		if ( denovo_scorefxn->has_nonzero_weight( core::scoring::rna_stub_coord_hack ) ) { // biasing information
			json.push_back( Pair( "xyz_bias_weight", denovo_scorefxn->get_weight( core::scoring::rna_stub_coord_hack ) ) );
			std::vector< Value > target_xyz( make_vector( Value(options_->output_jump_target_xyz().x()),
				Value(options_->output_jump_target_xyz().y()),
				Value(options_->output_jump_target_xyz().z()) ) );
			json.push_back( Pair( "target_xyz", target_xyz ) );
		}
		write_tensor_to_file( options_->output_histogram_file(), *jump_histogram_, json);
	}
}

} //output
} //denovo
} //rna
} //protocols
