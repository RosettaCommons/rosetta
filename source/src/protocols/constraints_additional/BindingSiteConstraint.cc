// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// Unit headers
#include <protocols/constraints_additional/BindingSiteConstraint.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <numeric/model_quality/rms.hh>

#include <utility/string_util.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>



static basic::Tracer tr("core.io.constraints");

namespace protocols {
namespace constraints_additional {

// database of binding sites
std::map< core::id::AtomID , numeric::xyzVector< core::Real > > protocols::constraints_additional::BindingSiteConstraint::rot_db;

///c-tor
BindingSiteConstraint::BindingSiteConstraint(
	utility::vector1< AtomID > const & atms,
	core::pose::Pose const &start_pose,
	core::scoring::ScoreType scoretype   /// ? TO DO -- give own scoretype
):
Constraint( scoretype ),
atms_(atms) {
	init ( start_pose );
}

void BindingSiteConstraint::init( core::pose::Pose const &start_pose ) {
	// grab binding-site atom positions from starting pose
	// error-checking is handled in parser -- here we assume #constr >= 2 and all atmids exist
	numeric::xyzVector< core::Real > com(0,0,0);
	int natoms = atms_.size();

	tgt_pos_.dimension( 3, natoms );
	for ( int i = 1; i <= natoms; ++i ) {
		numeric::xyzVector< core::Real > x_i = start_pose.xyz( atms_[i] );
		com += x_i;
		for ( int k = 0; k < 3; ++k )
			tgt_pos_(k+1,i) = x_i[k];
	}

	// center at the origin
	com /= atms_.size();
	for (int i=1; i<=(int)atms_.size(); ++i)
		for ( int k = 0; k < 3; ++k )
			tgt_pos_(k+1,i) -= com[k];

	// also create centroid constraints
	core::pose::Pose start_pose_copy = start_pose;
	if ( start_pose_copy.residue(1).residue_type_set().name() != core::chemical::CENTROID ){
		core::util::switch_to_residue_type_set( start_pose_copy, core::chemical::CENTROID );
	}
	com = numeric::xyzVector< core::Real >(0,0,0);

	tgt_pos_centroid_.dimension( 3, natoms );
	for ( int i = 1; i <= natoms; ++i ) {
		core::id::AtomID atm_i;
		if ( atms_[i].atomno() > 5 )
			atm_i = core::id::AtomID( start_pose_copy.residue( atms_[i].rsd()  ).natoms() , atms_[i].rsd() );
		else
			atm_i = atms_[i];
		numeric::xyzVector< core::Real > x_i = start_pose_copy.xyz( atm_i );
		com += x_i;
		for ( int k = 0; k < 3; ++k )
			tgt_pos_centroid_(k+1,i) = x_i[k];
	}

	// center at the origin
	com /= atms_.size();
	for (int i=1; i<=(int)atms_.size(); ++i)
		for ( int k = 0; k < 3; ++k )
			tgt_pos_centroid_(k+1,i) -= com[k];
}

/// ctor from a vector of atom positions (in lieu of a pose)
BindingSiteConstraint::BindingSiteConstraint(
	utility::vector1< AtomID > const & atms,
	ObjexxFCL::FArray2D< core::Real >  tgt_pos,
	ObjexxFCL::FArray2D< core::Real >  tgt_pos_centroid,
	core::scoring::ScoreType scoretype   /// ? TO DO -- give own scoretype
):
	Constraint( scoretype ),
	atms_(atms),
	tgt_pos_(tgt_pos),
	tgt_pos_centroid_(tgt_pos_centroid)
{}

///
void
BindingSiteConstraint::score( core::scoring::constraints::XYZ_Func const&, core::scoring::EnergyMap const &, core::scoring::EnergyMap & emap ) const {
	// filler
	//std::cerr << "BindingSiteConstraint::score( core::scoring::constraints::XYZ_Func const & xyz, core::scoring::EnergyMap const &, core::scoring::EnergyMap & emap ) " << std::endl;

	// get optimally-aligned RMS
	core::Real bs_score = 0.0;
	for (int i=1; i<=(int)atms_.size(); ++i) {
		core::Real this_dist = rot_db[ atms_[i] ].length();
		bs_score += this_dist*this_dist;
		//std::cerr << "   " << bs_score << " <= +" << this_dist << "^2\n";
	}

	emap[ this->score_type() ] += bs_score;
}

// do some pre-scoring calculations
void
BindingSiteConstraint::setup_for_scoring( core::scoring::constraints::XYZ_Func const & xyz, core::scoring::ScoreFunction const & ) const {
	// filler
	//std::cerr << "BindingSiteConstraint::setup_for_scoring() " << std::endl;

	// align the >target atoms< to the >template atoms<
	utility::vector1< numeric::xyzVector< core::Real > > templ_atms;
	utility::vector1< bool > align_centroid;

	for (int i=1; i<=(int)atms_.size(); ++i) {
		core::chemical::ResidueType const & rsd_type = xyz.residue( atms_[i].rsd()  ).type();
		core::id::AtomID atm_i;
		bool aln_cent_i = false;

		// is pose centroid?? then we need to remap atmids
		if ( rsd_type.residue_type_set().name() == core::chemical::CENTROID && atms_[i].atomno() > 5 ) {
			//std::cerr << "ATOM " << atms_[i].rsd() << " , " << atms_[i].atomno() << "  (" << pose.residue( atms_[i].rsd()  ).natoms() << ")" << std::endl;
			//std::cerr << "Remapping ATOM " << atms_[i].rsd() << " , " << atms_[i].atomno() << " to CENTROID" << std::endl;
			atm_i = core::id::AtomID( xyz.residue( atms_[i].rsd()  ).natoms() , atms_[i].rsd() );
			aln_cent_i = true;
		} else {
			atm_i = atms_[i];
		}

		templ_atms.push_back( xyz( atm_i ) );
		align_centroid.push_back( aln_cent_i );
	}
	pre_align (templ_atms, align_centroid);
}

// align the atoms
//   ... placing a vector  -- from each atom to the the rotated >target< atoms -- in the database
void
BindingSiteConstraint::pre_align(
                utility::vector1< numeric::xyzVector< core::Real > > const & templ_atms,
                utility::vector1< bool > const & align_centroid
                                ) const {
	int natoms = templ_atms.size();

	if ( natoms != (int)atms_.size() ) {
		utility_exit_with_message( "BindingSiteConstraint::align() bad argument" );
	}

	// copy coords into Real arrays
	ObjexxFCL::FArray2D< core::Real > tmpl_pos( 3, natoms );
	numeric::xyzVector< core::Real > com(0,0,0);
	for ( int i = 1; i <= natoms; ++i ) {
		numeric::xyzVector< core::Real > const &x_i = templ_atms[i];
		com += x_i;
		for ( int k = 0; k < 3; ++k )
			tmpl_pos(k+1,i) = x_i[k];
	}

	// center at the origin
	com /= natoms;
	for (int i=1; i<=natoms; ++i)
		for ( int k = 0; k < 3; ++k )
			tmpl_pos(k+1,i) -= com[k];

	// get optimal superposition
	// rotate TARGET to the TEMPLATE
	ObjexxFCL::FArray1D< numeric::Real > ww( natoms, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;
	ObjexxFCL::FArray2D< core::Real >  tgt_pos_mixed( 3 , natoms );

	for (int i=1; i<=natoms; ++i) {
		if (align_centroid[i]) {
			for ( int k = 1; k <= 3; ++k ) {
				tgt_pos_mixed(k,i) = tgt_pos_centroid_(k,i);
			}
		} else {
			for ( int k = 1; k <= 3; ++k ) {
				tgt_pos_mixed(k,i) = tgt_pos_(k,i);
			}
		}
	}
	ObjexxFCL::FArray2D< core::Real >  tgt_pos_copy( 3 , natoms );
	for (int i=1; i<=natoms; ++i)
		for ( int k = 1; k <= 3; ++k )
			tgt_pos_copy(k,i) = tgt_pos_mixed(k,i);

	numeric::model_quality::findUU( tgt_pos_mixed, tmpl_pos, ww, natoms, uu, ctx );

	numeric::xyzMatrix< core::Real > R;
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	// apply rotation
	// if one atom appears in multiple templates this will overwrite previous atom's position ...
	//   perhaps this should be checked for in the constraint parser
	for ( int i = 1; i <= natoms; ++i ) {
		numeric::xyzVector< core::Real > x_i( tgt_pos_copy(1,i) , tgt_pos_copy(2,i) , tgt_pos_copy(3,i) );
		numeric::xyzVector< core::Real > rx_i = R*x_i;
		numeric::xyzVector< core::Real > y_i( tmpl_pos(1,i) , tmpl_pos(2,i) , tmpl_pos(3,i) );
		//rot_db[ atms_[i] ] = rx_i;
		rot_db[ atms_[i] ] = (rx_i - y_i);  // store the direction template to target

		///////////////
		///////////////
		//std::cerr << "   tgt: ATOM " << i << ": " << tgt_pos_copy(1,i) << " , "
		//																					<< tgt_pos_copy(2,i) << " , "
		//																					<< tgt_pos_copy(3,i) << std::endl;
		//std::cerr << "   tmpl:ATOM " << i << ": " << tmpl_pos(1,i) << " , "
		//																					<< tmpl_pos(2,i) << " , "
		//																					<< tmpl_pos(3,i) << std::endl;
		//std::cerr << "   del: ATOM " << i << ": " << rot_db[ atms_[i] ][0] << " , "
		//                                          << rot_db[ atms_[i] ][1] << " , "
		//                                          << rot_db[ atms_[i] ][2] << "\n\n";
		///////////////
		///////////////
	}
}

// call the setup_for_derivatives for each constraint
void
BindingSiteConstraint::setup_for_derivatives( core::scoring::constraints::XYZ_Func const & xyz, core::scoring::ScoreFunction const &scfxn ) const {
	// filler
	//std::cerr << "BindingSiteConstraint::setup_for_derivatives() " << std::endl;
	setup_for_scoring( xyz, scfxn );
}

// atom deriv
void
BindingSiteConstraint::fill_f1_f2(
	AtomID const & atom,
	core::scoring::constraints::XYZ_Func const & xyz,
	core::Vector & F1,
	core::Vector & F2,
	core::scoring::EnergyMap const & weights
) const	{
	// filler
	//std::cerr << "BindingSiteConstraint::fill_f1_f2() " << std::endl;
	if ( std::find( atms_.begin() , atms_.end() , atom ) == atms_.end()) return;

	numeric::xyzVector< core::Real > atom_x = -xyz(atom);
	numeric::xyzVector< core::Real > const f2( rot_db[ atom ] );
	numeric::xyzVector< core::Real > atom_y = f2 + atom_x;   // a "fake" atom in the direcion of the gradient
	numeric::xyzVector< core::Real > const f1( atom_x.cross( atom_y ) );

	//std::cerr << "   f1 = " << f1[0] << " , " << f1[1] << " , " << f1[2] << std::endl;
	//std::cerr << "   f2 = " << f2[0] << " , " << f2[1] << " , " << f2[2] << std::endl;

	F1 -= f1*weights[ this->score_type() ]*2;
	F2 -= f2*weights[ this->score_type() ]*2;
}

std::string
BindingSiteConstraint::type() const {
	return "BindingSite";
}

///
core::Size
BindingSiteConstraint::natoms() const
{
	return atms_.size();
}

core::scoring::constraints::ConstraintOP
BindingSiteConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	utility::vector1< AtomID > atms_remap;

	for (int i=1; i<=(int) atms_.size(); ++i) {
		if ( seqmap[atms_[i].rsd()] != 0 ) {
			atms_remap.push_back( AtomID( atms_[i].atomno(), seqmap[atms_[i].rsd()] ) );
		}
	}

	if (atms_remap.size() > 2)
		return core::scoring::constraints::ConstraintOP( new BindingSiteConstraint( atms_remap , tgt_pos_, tgt_pos_centroid_ ) );
	else
		return NULL;
}

///
core::id::AtomID const &
BindingSiteConstraint::atom( core::Size const n ) const
{
	if ( n >= 1 && n <= atms_.size() )
		return atms_[n];
	else
		utility_exit_with_message( "BindingSiteConstraint::atom() bad argument" );

	return atms_[1]; // stop the compiler from complaining
}

void
BindingSiteConstraint::show( std::ostream& out ) const {
	out << "BindingSiteConstraint::show() " << std::endl;
}

void
BindingSiteConstraint::show_def( std::ostream&  out , core::pose::Pose const& /* pose */ ) const {
	out << "BindingSiteConstraint::show_def() " << std::endl;
}

core::Size
BindingSiteConstraint::show_violations( std::ostream& /*out*/, core::pose::Pose const& /*pose*/, core::Size /*verbose_level*/, core::Real /*threshold*/ ) const {
	/*
	if (verbose_level > 80) {
		out << "AtomPairConstraint ("
				<< core::pose.residue_type(atom1_.rsd() ).atom_name( atom1_.atomno() ) << ":" << atom1_.atomno() << "," << atom1_.rsd() << "-"
				<< core::pose.residue_type(atom2_.rsd() ).atom_name( atom2_.atomno() ) << ":" << atom2_.atomno() << "," << atom2_.rsd() << ") " ;
	};
	return func_->show_violations( out, dist( pose ), verbose_level );
	*/
	return 0;
}

/*
///@details one line definition "BindingSite atom1 res1 ... atomN resN"
void
BindingSiteConstraint::read_def(
	std::istream& data,
	core::pose::Pose const& pose,
	core::scoring::constraints::FuncFactory const& func_factory
) {
	core::Size res;
	std::string name;

	std::string line;
	// = data.getline();
	getline( data, line );
	std::vector< std::string > tokens = utility::split( line );

	// to do? atm name -> centroid conversion
	tr.Debug << "read: ";
	int npairs = tokens.size() / 2;
	for (int i=0; i<npairs; ++i) {
		name = tokens[2*i];
		res = (core::Size) atoi(tokens[2*i + 1].c_str());
		tr.Debug << "    " << name << " " << res;

		atms_.push_back( core::id::AtomID( pose.residue_type( res ).atom_index( name ), res ) );
		if ( res > pose.total_residue() ) {
			tr.Debug << "** ignored **";
			continue;
		}
	}

	// get positions from start pose
	init ( pose );
}
*/

/// @brief Format should look like:
/// Dunbrack seqpos_ rot_vec_pos_ rot_bin_ bonus_
void
BindingSiteConstraint::read_def(
	std::istream & line_stream,
	core::pose::Pose const & pose,
	core::scoring::constraints::FuncFactory const & /* func_factory */
) {

	core::Size res; // ?
	std::string name; // ?

	// to do? atm name -> centroid conversion
	utility::vector1< core::id::AtomID > atms;
	tr.Debug << "read: ";
	while ( line_stream >> name >> res ) {
		tr.Debug << "   " << name << " " << res ;
		atms.push_back( core::id::AtomID( pose.residue_type( res ).atom_index( name ), res ) );
		if ( res > pose.total_residue() ) {
			tr.Debug << "** ignored **";
			continue;
		}
	}
	tr.Debug << std::endl;
	atms_ = atms;
	init( pose );
}



}
}
