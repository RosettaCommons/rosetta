// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/BackboneStubConstraint.cc
///
/// @brief
/// @author John Karanicolas, Sarel Fleishman


#include <core/scoring/constraints/BackboneStubConstraint.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/PeriodicFunc.hh>
#include <basic/Tracer.hh>

// used to make temporary alanines for gly cst's
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

//Auto Headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>

namespace core {
namespace scoring {
namespace constraints {


static thread_local basic::Tracer tr( "core.scoring.constraints.BackboneStubConstraint" );

utility::pointer::owning_ptr< AngleConstraint > BackboneStubConstraint::ang_cst_(0);


std::map< std::string, core::id::AtomID >
auto_detect_atoms(pose::Pose const & pose, Size const seqpos)
{
	//kdrew: stores four atoms ids
	std::map< std::string, core::id::AtomID > atom_ids;

	// constraint depends on CB, CA, C coordinates
	conformation::Residue const & rsd( pose.residue(seqpos) );

	//sidechain_atom_id_ = AtomID( rsd.atom_index(sidechain_atom_name), seqpos_ );
	//atom_ids_.push_back(sidechain_atom_id_);

	for ( Size j=1; j<= rsd.natoms(); ++j )
	{
			if( !rsd.atom_is_backbone(j) && !rsd.atom_is_hydrogen(j) )
			{
				core::chemical::AtomIndices sidechain_atom_bonded_neighbors = rsd.bonded_neighbor(j);
				for ( Size k=1; k <= sidechain_atom_bonded_neighbors.size(); ++k )
				{
					if( rsd.atom_is_backbone(sidechain_atom_bonded_neighbors[k] ) )
					{
						tr.Debug << "ATOM: " << j << ' ' << rsd.name3() << ' ' << rsd.atom_name(j) << " atom_is_backbone: " << rsd.atom_is_backbone(j) << " atom_is_hydrogen: "<< rsd.atom_is_hydrogen(j) << std::endl;
						//kdrew: store sidechain atom  (CB for peptides)
						atom_ids["sidechain_atom_id"] = core::id::AtomID( j, seqpos );

						tr.Debug << "\t\tbonded_neighbors: " << sidechain_atom_bonded_neighbors[k] << " name: " << rsd.atom_name(sidechain_atom_bonded_neighbors[k]) << std::endl;
						//kdrew: store primary backbone atom (CA for peptides)
						atom_ids["primary_backbone_atom_id"] = core::id::AtomID( sidechain_atom_bonded_neighbors[k], seqpos );

						if( rsd.has_lower_connect() && rsd.lower_connect_atom() == atom_ids["primary_backbone_atom_id"].atomno() )
						{
							//kdrew: primary_backbone_atom is lower connect, add upper connect of seqpos - 1 as backbone1
							conformation::Residue const & rsd_lower( pose.residue(seqpos-1));
							if( rsd_lower.has_upper_connect() ) {
								atom_ids["backbone_atom1_id"] = core::id::AtomID( rsd_lower.upper_connect_atom(), seqpos-1);
							}
						}
						if( rsd.has_upper_connect() && rsd.upper_connect_atom() == atom_ids["primary_backbone_atom_id"].atomno() )
						{
							//kdrew: primary_backbone_atom is upper connect, add lower connect of seqpos + 1 as backbone2
							conformation::Residue const & rsd_upper( pose.residue(seqpos+1));
							if( rsd_upper.has_lower_connect() ) {
								atom_ids["backbone_atom2_id"] = core::id::AtomID( rsd_upper.lower_connect_atom(), seqpos+1);
							}
						}

						core::chemical::AtomIndices primary_backbone_atom_bonded_neighbors = rsd.bonded_neighbor(sidechain_atom_bonded_neighbors[k]);
						for ( Size m=1; m <= primary_backbone_atom_bonded_neighbors.size(); ++m )
						{
							if( rsd.atom_is_backbone(primary_backbone_atom_bonded_neighbors[m] ) && !rsd.atom_is_hydrogen( primary_backbone_atom_bonded_neighbors[m]) )
							{
								if ( atom_ids.find("backbone_atom1_id") == atom_ids.end() )
								{
									atom_ids["backbone_atom1_id"] = core::id::AtomID( primary_backbone_atom_bonded_neighbors[m], seqpos );
								}
								else if ( atom_ids.find("backbone_atom2_id") == atom_ids.end() )
								{
									atom_ids["backbone_atom2_id"] = core::id::AtomID( primary_backbone_atom_bonded_neighbors[m], seqpos );

								}
								tr.Debug << "\t\t\tbonded_neighbors: " << primary_backbone_atom_bonded_neighbors[m] << " name: " << rsd.atom_name(primary_backbone_atom_bonded_neighbors[m]) << std::endl;

								/*
								if( rsd.lower_connect_atom() == primary_backbone_atom_bonded_neighbors[m] )
								{
                                	TR << "\t\t\tLOWER connect atom" << std::endl;
									//chemical::ResidueConnection lower_conn = rsd.lower_connect();
									conformation::Residue const & rsd_lower( pose.residue(seqpos_-1));
									TR << "\t\t\t\tconnected atom: " << rsd_lower.atom_name(rsd_lower.upper_connect_atom()) << std::endl;

								}
								if( rsd.upper_connect_atom() == primary_backbone_atom_bonded_neighbors[m] )
								{
                                	TR << "\t\t\tUPPER connect atom" << std::endl;
									conformation::Residue const & rsd_upper( pose.residue(seqpos_+1));
									TR << "\t\t\t\tconnected atom: " << rsd_upper.atom_name(rsd_upper.lower_connect_atom()) << std::endl;
								}
								core::chemical::AtomIndices secondary_backbone_atom_bonded_neighbors = rsd.bonded_neighbor(primary_backbone_atom_bonded_neighbors[m]);
								for ( Size n=1; n <= secondary_backbone_atom_bonded_neighbors.size(); ++n )
								{
									TR << "\t\t\t\tbonded_neighbors: " << secondary_backbone_atom_bonded_neighbors[n] << " name: " << rsd.atom_name(secondary_backbone_atom_bonded_neighbors[n]) << std::endl;

								}
								*/
							}

						}
					}
				}
				//break;
			}
	}

	return atom_ids;
}

//kdrew: this is the non-peptidic implementation (peptoids, beta-peptides, etc),
//kdrew: takes in a two residue target pose, auto detects atoms
BackboneStubConstraint::BackboneStubConstraint(
	pose::Pose const & pose,
	Size const seqpos,
	AtomID const & fixed_atom_id,
	pose::Pose const & target_pose,
	Size const target_seqpos,
	core::Real const & superposition_bonus,
	core::Real const & force_constant
):
	Constraint( core::scoring::backbone_stub_constraint ),
	superposition_bonus_( superposition_bonus ),
	force_constant_( force_constant ),
	seqpos_( seqpos ),
	fixed_atom_id_( fixed_atom_id )
{

	// store info about the target residue
	assert( target_pose.residue(target_seqpos).is_protein() );
	if ( (target_pose.residue(target_seqpos).aa() == chemical::aa_gly ) ) {
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
	}

	std::map<std::string, AtomID> target_atom_ids = auto_detect_atoms(target_pose, target_seqpos);

	sidechain_target_ = target_pose.residue(target_atom_ids["sidechain_atom_id"].rsd()).xyz(target_atom_ids["sidechain_atom_id"].atomno());
	primary_backbone_target_ = target_pose.residue(target_atom_ids["primary_backbone_atom_id"].rsd()).xyz(target_atom_ids["primary_backbone_atom_id"].atomno());
   	backbone1_target_ = target_pose.residue(target_atom_ids["backbone_atom1_id"].rsd()).xyz(target_atom_ids["backbone_atom1_id"].atomno());
   	backbone2_target_ = target_pose.residue(target_atom_ids["backbone_atom2_id"].rsd()).xyz(target_atom_ids["backbone_atom2_id"].atomno());

	//primary_backbone_target_ = target_pose.residue(target_seqpos).xyz(primary_backbone_atom_name);
	//backbone1_target_ = target_pose.residue(target_seqpos).xyz(backbone_atom1_name);
	//backbone2_target_ = target_pose.residue(target_seqpos).xyz(backbone_atom2_name);

	primary_sidechain_target_ = sidechain_target_ - primary_backbone_target_;
	backbone1_2_target_ = backbone2_target_ - backbone1_target_;


	conformation::Residue const & rsd( pose.residue(seqpos) );

	assert( rsd.is_protein() );
	if ( (rsd.aa() == chemical::aa_gly) ) {
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
	}

	std::map<std::string, AtomID> atom_ids = auto_detect_atoms(pose, seqpos);
	atom_ids_.push_back(atom_ids["sidechain_atom_id"]);
	atom_ids_.push_back(atom_ids["primary_backbone_atom_id"]);
	atom_ids_.push_back(atom_ids["backbone_atom2_id"]);
	atom_ids_.push_back(atom_ids["backbone_atom1_id"]);

	//sidechain_atom_id_ = AtomID( rsd.atom_index(sidechain_atom_name), seqpos_ );
	//atom_ids_.push_back(sidechain_atom_id_);
	//primary_atom_id_ = AtomID( rsd.atom_index(primary_backbone_atom_name), seqpos_ );
	//atom_ids_.push_back(primary_atom_id_);
	//backbone_atom2_id_ = AtomID( rsd.atom_index(backbone_atom2_name), seqpos_ );
	//atom_ids_.push_back(backbone_atom2_id_);
	//backbone_atom1_id_ = AtomID( rsd.atom_index(backbone_atom1_name), seqpos_ );
	//atom_ids_.push_back( backbone_atom1_id_ );



	// need a fixed reference point on the pose
	atom_ids_.push_back(fixed_atom_id_);
	// we don't actually need to save the coors of the reference point,
	// but will do so to ensure that it doesn't move
	fixed_reference_point_ = pose.xyz(fixed_atom_id_);

	// to get access to AngleConstraint derivatives
	if ( ang_cst_ == 0 ) {
		// note: PeriodicFunc has functional form y = ( k * cos(n * (x - x0) ) ) + C
		func::FuncOP cos_func = new func::PeriodicFunc(0., 1., 1., 0.);
		ang_cst_ = AngleConstraintOP( new AngleConstraint( cos_func ) );
	}
}

//kdrew: this is the non-peptidic implementation (peptoids, beta-peptides, etc),
//kdrew: specifies atom names
BackboneStubConstraint::BackboneStubConstraint(
	pose::Pose const & pose,
	Size const seqpos,
	AtomID const & fixed_atom_id,
	conformation::Residue const & target_rsd,
	core::Real const & superposition_bonus,
	core::Real const & force_constant,
	std::string const & sidechain_atom_name, // CB
	std::string const & primary_backbone_atom_name, //CA
	std::string const & backbone_atom1_name, //N
	std::string const & backbone_atom2_name //C
):
	Constraint( core::scoring::backbone_stub_constraint ),
	superposition_bonus_( superposition_bonus ),
	force_constant_( force_constant ),
	seqpos_( seqpos ),
	fixed_atom_id_( fixed_atom_id )
{

	// store info about the target residue
	assert( target_rsd.is_protein() );
	if ( (target_rsd.aa() == chemical::aa_gly ) ) {
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
	}

	sidechain_target_ = target_rsd.xyz(sidechain_atom_name);
	primary_backbone_target_ = target_rsd.xyz(primary_backbone_atom_name);
	backbone1_target_ = target_rsd.xyz(backbone_atom1_name);
	backbone2_target_ = target_rsd.xyz(backbone_atom2_name);
	primary_sidechain_target_ = sidechain_target_ - primary_backbone_target_;
	backbone1_2_target_ = backbone2_target_ - backbone1_target_;

	// constraint depends on CB, CA, C coordinates
	conformation::Residue const & rsd( pose.residue(seqpos_) );

	assert( rsd.is_protein() );
	if ( (rsd.aa() == chemical::aa_gly) ) {
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
	}
	sidechain_atom_id_ = AtomID( rsd.atom_index(sidechain_atom_name), seqpos_ );
	atom_ids_.push_back(sidechain_atom_id_);
	primary_atom_id_ = AtomID( rsd.atom_index(primary_backbone_atom_name), seqpos_ );
	atom_ids_.push_back(primary_atom_id_);
	backbone_atom2_id_ = AtomID( rsd.atom_index(backbone_atom2_name), seqpos_ );
	atom_ids_.push_back(backbone_atom2_id_);
	backbone_atom1_id_ = AtomID( rsd.atom_index(backbone_atom1_name), seqpos_ );
	atom_ids_.push_back( backbone_atom1_id_ );

	// need a fixed reference point on the pose
	atom_ids_.push_back(fixed_atom_id_);
	// we don't actually need to save the coors of the reference point,
	// but will do so to ensure that it doesn't move
	fixed_reference_point_ = pose.xyz(fixed_atom_id_);

	// to get access to AngleConstraint derivatives
	if ( ang_cst_ == 0 ) {
		// note: PeriodicFunc has functional form y = ( k * cos(n * (x - x0) ) ) + C
		func::FuncOP cos_func = new func::PeriodicFunc(0., 1., 1., 0.);
		ang_cst_ = AngleConstraintOP( new AngleConstraint( cos_func ) );
	}
}

//kdrew: this is the peptide centric implementation
BackboneStubConstraint::BackboneStubConstraint(
	pose::Pose const & pose,
	Size const seqpos,
	AtomID const & fixed_atom_id,
	conformation::Residue const & target_rsd,
	core::Real const & superposition_bonus,
	core::Real const & CB_force_constant
):
	Constraint( core::scoring::backbone_stub_constraint ),
	superposition_bonus_( superposition_bonus ),
	force_constant_( CB_force_constant ),
	seqpos_( seqpos ),
	fixed_atom_id_( fixed_atom_id )
{

	// store info about the target residue
	assert( target_rsd.is_protein() );
	if ( (target_rsd.aa() == chemical::aa_gly ) ) {
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
	}


	sidechain_target_ = target_rsd.xyz("CB");
	primary_backbone_target_ = target_rsd.xyz("CA");
	backbone2_target_ = target_rsd.xyz("C");
	backbone1_target_ = target_rsd.xyz("N");
	primary_sidechain_target_ = sidechain_target_ - primary_backbone_target_;
	backbone1_2_target_ = backbone2_target_ - backbone1_target_;

	// constraint depends on CB, CA, C coordinates
	conformation::Residue const & rsd( pose.residue(seqpos_) );

	assert( rsd.is_protein() );
	if ( (rsd.aa() == chemical::aa_gly) ) {
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
	}
	sidechain_atom_id_ = AtomID( rsd.atom_index("CB"), seqpos_ );
	atom_ids_.push_back(sidechain_atom_id_);
	primary_atom_id_ = AtomID( rsd.atom_index("CA"), seqpos_ );
	atom_ids_.push_back(primary_atom_id_);
	backbone_atom2_id_ = AtomID( rsd.atom_index("C"), seqpos_ );
	atom_ids_.push_back(backbone_atom2_id_);
	backbone_atom1_id_ = AtomID( rsd.atom_index("N"), seqpos_ );
	atom_ids_.push_back( backbone_atom1_id_ );

	// need a fixed reference point on the pose
	atom_ids_.push_back(fixed_atom_id_);
	// we don't actually need to save the coors of the reference point,
	// but will do so to ensure that it doesn't move
	fixed_reference_point_ = pose.xyz(fixed_atom_id_);

	// to get access to AngleConstraint derivatives
	if ( ang_cst_ == 0 ) {
		// note: PeriodicFunc has functional form y = ( k * cos(n * (x - x0) ) ) + C
		func::FuncOP cos_func = new func::PeriodicFunc(0., 1., 1., 0.);
		ang_cst_ = AngleConstraintOP( new AngleConstraint( cos_func ) );
	}
}

core::Size
BackboneStubConstraint::seqpos() const
{ return seqpos_; }

void BackboneStubConstraint::show( std::ostream& out ) const
{
	out << "BackboneStubCst Seqpos: " << seqpos_ << "    bonus: " << superposition_bonus_ << "    CB force constant: " << force_constant_ << std::endl;
}


bool
BackboneStubConstraint::operator == ( Constraint const & other_cst ) const
{

	if( !dynamic_cast< BackboneStubConstraint const * > ( &other_cst ) ) return false;

	BackboneStubConstraint const & other( static_cast< BackboneStubConstraint const & > (other_cst) );

	if( superposition_bonus_ != other.superposition_bonus_ ) return false;
	if( force_constant_ != other.force_constant_ ) return false;
	if( seqpos_ != other.seqpos_ ) return false;

	if( sidechain_atom_id_ != other.sidechain_atom_id_ ) return false;
	if( primary_atom_id_ != other.primary_atom_id_ ) return false;
	if( backbone_atom2_id_ != other.backbone_atom2_id_ ) return false;
	if( backbone_atom1_id_ != other.backbone_atom1_id_ ) return false;

	if( sidechain_target_ != other.sidechain_target_ ) return false;
	if( primary_backbone_target_ != other.primary_backbone_target_ ) return false;
	if( backbone2_target_ != other.backbone2_target_ ) return false;
	if( backbone1_target_ != other.backbone1_target_ ) return false;
	if( primary_sidechain_target_ != other.primary_sidechain_target_ ) return false;
	if( backbone1_2_target_ != other.backbone1_2_target_ ) return false;

	if( fixed_atom_id_ != other.fixed_atom_id_ ) return false;
	if( fixed_reference_point_ != other.fixed_reference_point_ ) return false;

	return true;
}



// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
// the constraint score, C, is given by C = (b + kx^2) * z, where b is the score bonus for the stub, k is the user-defined CB force constant
// z is the dot product between the Ca-Cb vectors of the stub and the current
// position on the scaffold times the C-N vectors of the stub and scaffold position.
// The intuition in this is that a stub that lands
// perfectly on the scaffold's bacbkone position AND for which the two vectors
// match the scaffold's perfectly would give maximal bonus. This choice of vectors is
// orthogonal and so provides unbiased resolution for conformation space.
// Note that we don't know in advance which phi/psi angles a scaffold position would have, and
// that may change the C_N vector by up to 25o from the stub's. Since the cosine of small angular
// differences is close to 1, this uncertainty would not cause large problems.
// Since the constraint does not produce positive values, the effect of a stub
// on the energy is only if the distance between the Cb's of stub and scaffold is
// sqrt( -b / k ). As a rule of thumb, b=-4, so k can be decided in a way that will
// determine the radius of the effect of a stub on pulling a scaffold towards it: k = 4/dist^2, where dist is the preferred radius
void
BackboneStubConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const {

	if ( weights[ this->score_type() ] == 0 ) return;

	core::conformation::Residue const & curr_rsd = xyz_func.residue(seqpos_);
	assert( curr_rsd.is_protein() );

	// verify that the fixed reference point is in the same place
	core::Vector curr_ref_location = xyz_func(fixed_atom_id_);
	core::Real ref_dist = curr_ref_location.distance_squared( fixed_reference_point_ );
	if ( ref_dist > 1E-8 ) {
		tr.Error << "ERROR - BackboneStubConstraint requires a fixed reference atom, but this atom has moved!!" << std::endl;
		tr.Error << "Reference location was " << fixed_reference_point_.x() << ", " << fixed_reference_point_.y() << ", " << fixed_reference_point_.z() << ", now it's " << curr_ref_location.x() << ", " << curr_ref_location.y() << ", " << curr_ref_location.z() << std::endl;
		utility_exit_with_message( "ERROR - BackboneStubConstraint requires a fixed reference atom, but this atom has moved!!" );
	}

	// return a value between superposition_bonus_ (-ve) and zero
	core::Real cst_val(superposition_bonus_);
	assert( cst_val < 0. );

	// apply a harmonic constraint on the CB's
	core::Vector CB_curr;
	if ( (curr_rsd.aa() == chemical::aa_gly) ) {
/// SJF I'm disabling the option for glycine backbone stub constraints b/c it hardly seems useful and is something of a pain
/// to deal with since changing curr_rsd to a Residue const &
/*		if( basic::options::option[basic::options::OptionKeys::hotspot::allow_gly]() ) {
			// make an alanine
			core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
			core::conformation::ResidueOP ala_copy = core::conformation::ResidueFactory::create_residue( alatype );
			// move ala on top of gly
			ala_copy->orient_onto_residue( curr_rsd );
			curr_rsd = *(ala_copy); // gly => ala
		}
		else { */
			tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
			return;
//		}
	}
	CB_curr = curr_rsd.xyz("CB");

	core::Real CB_d2 = CB_curr.distance_squared( sidechain_target_ );
	cst_val += force_constant_ * CB_d2;
 	if ( cst_val > 0. ) return;

///// SJF the following are the dot product results (or cosines) of the various
///// vectors in an amino-acid:
///// ca-cb x ca_c = 0.34
///// ca-cb x cb_c = 0.82
///// ca-cb x n-c = 0.02
///// So, by far the most orthogonal choice is to use ca-b and n-c

	// multiply by the cos of the angle between the CB-CA and CB'-CA' vectors
	core::Vector const CA_curr = curr_rsd.xyz("CA");
//	core::Vector const CB_CA_curr = CB_curr - CA_curr;
	cst_val *= ang_cst_->score( CB_curr, CA_curr, CA_curr + primary_sidechain_target_);
	if ( cst_val > -1E-10 ) return;

	// multiply by the cos of the angle between the C-N and C'-N' vectors
	core::Vector const C_curr = curr_rsd.xyz("C");
	core::Vector const N_curr = curr_rsd.xyz("N");
//	core::Vector const C_N_curr = C_curr - N_curr;

	cst_val *= ang_cst_->score( C_curr, N_curr, N_curr + backbone1_2_target_ );
	if ( cst_val > -1E-10 ) return;

	emap[ this->score_type() ] += cst_val;
}

void
BackboneStubConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	if ( weights[ this->score_type() ] == 0 ) return;

	core::conformation::Residue const & curr_rsd = xyz.residue( seqpos_ );

	assert( curr_rsd.is_protein() );

	if ( ( atom != sidechain_atom_id_ ) && ( atom != primary_atom_id_ ) && ( atom != backbone_atom2_id_ ) && ( atom != backbone_atom1_id_ ) ) {
		return;
	}

	// return a value between superposition_bonus_ (-ve) and zero
	core::Real cst_val(superposition_bonus_);
	assert( cst_val < 0. );

	// start by computing the cst value normally, collecting score components along the way. if cst_val is non-negative no derivative is computed

	if ( (curr_rsd.aa() == chemical::aa_gly) ) {
/*		if( basic::options::option[basic::options::OptionKeys::hotspot::allow_gly] ) {
			// make an alanine
			core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
			core::conformation::ResidueOP ala_copy = core::conformation::ResidueFactory::create_residue( alatype );
			// move ala on top of gly
			ala_copy->orient_onto_residue( *curr_rsd );
			curr_rsd = *ala_copy; // gly => ala
		}
		else { */
		tr.Warning << "WARNING - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
		return;
//		}
	}

	// apply a harmonic constraint on the CB's
	core::Vector const CB_curr( curr_rsd.xyz("CB") );

	core::Real const CB_d2 = CB_curr.distance_squared( sidechain_target_ );
	core::Real const CB_pos_term = force_constant_ * CB_d2;
	cst_val += CB_pos_term;
	if ( cst_val > -1E-10 ) return;

	// multiply by the cos of the angle between the CB-CA and CB'-CA' vectors
	core::Vector const CA_curr = curr_rsd.xyz("CA");
	core::Vector const CB_CA_curr = CB_curr - CA_curr;
	core::Real const CB_CA_angle_term( ang_cst_->score( CB_curr, CA_curr, CA_curr + primary_sidechain_target_ ) );
	if( CB_CA_angle_term <= 1E-10 ) return;

	// multiply by the cos of the angle between the C-N and C'-N' vectors
	core::Vector const C_curr = curr_rsd.xyz("C");
	core::Vector const N_curr = curr_rsd.xyz("N");
	core::Vector const C_N_curr = C_curr - N_curr;
	core::Real const C_N_angle_term( ang_cst_->score( C_curr, N_curr, N_curr + backbone1_2_target_ ) );
	if ( C_N_angle_term <= 1E-10 ) return;

	core::Real const constant_dist_term( weights[ this->score_type() ] );
	core::Real const constant_ang_term( constant_dist_term  * ( superposition_bonus_ + CB_pos_term ) );

	// contribution from differentiating CB_dist ( * CA_angle_term * C_angle_term from product rule)
	// and then adding the effect of the angular constraint on cb using the chain rule
	if ( atom == sidechain_atom_id_ ) {
// the effects of the coordinate constraint
		Vector const CB_f2( CB_curr - sidechain_target_ );
		core::Real const CB_dist( CB_f2.length() );
		core::Real const CB_deriv = 2. * force_constant_ * CB_dist;
		if ( CB_deriv != 0.0 && CB_dist != 0.0 ) {
			Vector const CB_f1( CB_curr.cross( sidechain_target_ ) );
			F1 += ( ( CB_deriv / CB_dist ) * CB_f1 ) * CB_CA_angle_term * C_N_angle_term * constant_dist_term;
			F2 += ( ( CB_deriv / CB_dist ) * CB_f2 ) * CB_CA_angle_term * C_N_angle_term * constant_dist_term;
		}

// the angular constraint on cb
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( CB_curr/*p1*/, CA_curr/*p2*/, CA_curr + primary_sidechain_target_/*p3*/, partial_F1, partial_F2 );
		F1 += partial_F1 * C_N_angle_term * constant_ang_term;
		F2 += partial_F2 * C_N_angle_term * constant_ang_term;
		return;
	}

// the angular constraint on ca
// See N atom for explanation on the derivs
	if ( atom == primary_atom_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( CA_curr, CA_curr - CB_CA_curr, CA_curr + primary_sidechain_target_ - CB_CA_curr, partial_F1, partial_F2 );
		F1 += -partial_F1 * C_N_angle_term * constant_ang_term;
		F2 += -partial_F2 * C_N_angle_term * constant_ang_term;
		return;
	}

// the angular constraint on c
	if ( atom == backbone_atom2_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( C_curr/*p1*/, N_curr/*p2*/, N_curr + backbone1_2_target_/*p3*/, partial_F1, partial_F2 );
		F1 += partial_F1 * CB_CA_angle_term * constant_ang_term;
		F2 += partial_F2 * CB_CA_angle_term * constant_ang_term;
		return;
	}

// the angular constraint on n
// This is tricky, and thanks to Frank for coming up with this!
// Ang_cst uses the actual coordinates of p1 in computing the derivatives, requiring it to be constant.
// But in our case, p1 changes with respect to p2->p3. Solution: translate the vectors by -C_N_curr.
// Now, p1 is at the site of the nonchanging vector, and so the coordinates are safe for ang_cst.
// However, the derivative for the angular constraint is reversed, so we multiply by -1
	if( atom == backbone_atom1_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( N_curr, N_curr - C_N_curr, N_curr + backbone1_2_target_ - C_N_curr, partial_F1, partial_F2 );
		F1 += -partial_F1 * CB_CA_angle_term * constant_ang_term;
		F2 += -partial_F2 * CB_CA_angle_term * constant_ang_term;
	}
	return;
}

ConstraintOP BackboneStubConstraint::clone() const
{
	return ConstraintOP( new BackboneStubConstraint( *this ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP BackboneStubConstraint::remapped_clone( pose::Pose const& /*src*/, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {

	core::Size new_seqpos = seqpos_;
	AtomID new_fixed_atom_id = fixed_atom_id_;

	if ( smap ) {
		new_seqpos = (*smap)[ seqpos_ ];
		new_fixed_atom_id.rsd() = (*smap)[ fixed_atom_id_.rsd() ];
		if( new_seqpos == 0 ) return NULL;
	}

	// make an alanine
	core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( alatype );
	ala->set_xyz("CB",sidechain_target_);
	ala->set_xyz("CA",primary_backbone_target_);
	ala->set_xyz("C",backbone2_target_);
	ala->set_xyz("N",backbone1_target_);

	return new BackboneStubConstraint(dest, new_seqpos, new_fixed_atom_id, *ala, superposition_bonus_, force_constant_ );
}


/*ConstraintOP
BackboneStubConstraint::remap_resid(
	core::id::SequenceMapping const & seqmap
) const {
	for(
  if ( seqmap[atom1_.rsd()] != 0 && seqmap[atom2_.rsd()] != 0 && seqmap[atom3_.rsd()] != 0 && seqmap[atom4_.rsd()] != 0 ) {
    AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] ),
      remap_a2( atom2_.atomno(), seqmap[atom2_.rsd()] ),
			remap_a3( atom3_.atomno(), seqmap[atom3_.rsd()] ),
			remap_a4( atom4_.atomno(), seqmap[atom4_.rsd()] );
    return ConstraintOP( new DihedralConstraint( remap_a1, remap_a2, remap_a3, remap_a4, this->func_ ) );
  } else {
    return NULL;
  }
}
*/


} // namespace constraints
} // namespace scoring
} // namespace core
