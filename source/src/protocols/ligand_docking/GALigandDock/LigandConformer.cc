// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <numeric/Quaternion.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.LigandConformer" );

LigandConformer::LigandConformer() {
	init_params();
}

LigandConformer::~LigandConformer() {}

LigandConformer::LigandConformer(
	core::pose::PoseCOP pose,
	core::Size ligid,
	utility::vector1< core::Size > movingscs
) {
	init_params();
	initialize( pose, ligid, movingscs );
}

void
LigandConformer::init_params() {
	torsmutationRate_ = 0.1;
	rtmutationRate_ = 0.5;
	transmutWidth_ = 2.0;
	rotmutWidth_ = 60.0;
	ligchimutWidth_ = 120.0;
	protchimutWidth_ = 120.0;
	generation_tag_ = "";
	rg_ = 0.0;
}

void
LigandConformer::initialize(
	core::pose::PoseCOP pose,
	core::Size ligid,
	utility::vector1< core::Size > movingscs
) {
	ligid_ = ligid;
	ref_pose_ = pose;
	movingscs_ = movingscs;
	update_conf( pose );

	// get ligand chi dependence... this is a stupid way, should have a better way ...
	core::conformation::Residue ligand( pose->residue( ligid ) );
	core::Size nligchi( ligand.nchi());
	ligandchi_downstream_.resize( nligchi );

	utility::vector1< core::Size > atmindex_defining_chi( nligchi, 0 );

	for ( core::Size ichi = 1; ichi <= nligchi; ++ichi ) {
		utility::vector1< core::Size > const &chiatms = ligand.chi_atoms( ichi );
		atmindex_defining_chi[ichi] = (ligand.atom_base(chiatms[2]) == chiatms[3])? chiatms[1] : chiatms[4];
	}

	for ( core::Size ichi = 1; ichi <= nligchi; ++ichi ) {
		core::Size iatm( atmindex_defining_chi[ichi] );
		core::Size ibase( ligand.atom_base(iatm) );
		while ( ibase != ligand.nbr_atom() ) { // recurrsive until reaches to nbr atom
			if ( atmindex_defining_chi.contains(iatm) ) {
				core::Size chi_parent_of_ichi = atmindex_defining_chi.index_of(iatm);
				if ( chi_parent_of_ichi != ichi ) ligandchi_downstream_[chi_parent_of_ichi].push_back( ichi );
			}
			iatm = ibase;
			ibase = ligand.atom_base(iatm);
		}
	}

	for ( core::Size ichi = 1; ichi <= nligchi; ++ichi ) {
		TR.Debug << "chi " << ichi << " influences downstream chis: ";
		for ( core::Size j = 1; j <= ligandchi_downstream_[ichi].size(); ++j ) {
			TR.Debug << " " << ligandchi_downstream_[ichi][j];
		}
		TR.Debug << std::endl;
	}

}


// update internal representation based on this conformation
void
LigandConformer::update_conf( core::pose::PoseCOP pose ) {

	rb_.resize(7);
	proteinchis_.resize(movingscs_.size());
	proteinrestypes_.resize(movingscs_.size());

	// skip if ligand-only case
	if ( pose->size() > 1 ) {
		core::Size jumpid = pose->fold_tree().get_jump_that_builds_residue( ligid_ );
		core::kinematics::Jump const &ligjump = pose->jump( jumpid );
		numeric::xyzVector< core::Real > T = ligjump.get_translation();
		numeric::Quaternion< core::Real > Q;
		numeric::R2quat( ligjump.get_rotation(), Q );
		rb_[1] = Q.w(); rb_[2] = Q.x(); rb_[3] = Q.y(); rb_[4] = Q.z();
		rb_[5] = T.x(); rb_[6] = T.y(); rb_[7] = T.z();

		// movable sidechains
		for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
			core::Size nchi = pose->residue(movingscs_[i]).nchi();
			utility::vector1< core::Real > chis_i(nchi,0.0);
			for ( core::Size j=1; j<=nchi; ++j ) {
				chis_i[j] = pose->chi( j, movingscs_[i] );
			}
			proteinchis_[i] = chis_i;
			proteinrestypes_[i] = pose->residue(movingscs_[i]).type_ptr();
		}
	}

	// internal torsions
	core::Size nchi = pose->residue(ligid_).nchi();
	ligandchis_.resize(nchi );
	for ( core::Size j=1; j<=nchi; ++j ) {
		ligandchis_[j] = pose->chi( j, ligid_ );
	}

	// copy ligandxyz
	core::conformation::Residue const &lig( pose->residue(ligid_));
	ligandxyz_.resize( lig.nheavyatoms() );

	core::Vector com( 0.0 );
	for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
		ligandxyz_[iatm] = lig.xyz( iatm );
		com += lig.xyz( iatm );
	}
	com /= lig.nheavyatoms();

	// for radius of gyration
	rg_ = 0.0;
	for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
		rg_ += com.distance_squared( ligandxyz_[iatm] );
	}
	rg_ /= lig.nheavyatoms();
	rg_ = std::sqrt(rg_);

	ligandxyz_synced_ = true;
}


// generate a pose based on this conformation
//    basically the reverse of update_conf
void
LigandConformer::to_pose( core::pose::PoseOP pose ) const {
	if ( pose->total_residue() == 0 ) *pose = *ref_pose_;

	// jump
	numeric::Quaternion< core::Real > Q( rb_[1], rb_[2], rb_[3], rb_[4] );
	numeric::xyzVector< core::Real >  T( rb_[5], rb_[6], rb_[7] );
	numeric::xyzMatrix< core::Real >  R;
	numeric::quat2R( Q, R );

	if ( pose->size() > 1 ) {
		core::Size jumpid = pose->fold_tree().get_jump_that_builds_residue( ligid_ );
		core::kinematics::Jump ligjump = pose->jump( jumpid );
		ligjump.set_translation( T );
		ligjump.set_rotation( R );
		pose->set_jump( jumpid, ligjump );

		// protein torsions
		for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
			core::Size resid_i = movingscs_[i];

			if ( pose->residue_type(resid_i).name() != proteinrestypes_[i]->name() ) {   // is there a better way to compare restypes?
				core::conformation::Residue newres(proteinrestypes_[i], true);
				pose->replace_residue( resid_i, newres, true );
			}
			for ( core::Size j=1; j<=proteinchis_[i].size(); ++j ) {
				pose->set_chi( j, resid_i, proteinchis_[i][j] );
			}
		}
	}

	// internal torsions
	core::Size nchi = pose->residue_type(ligid_).nchi();
	for ( core::Size j=1; j<=nchi; ++j ) {
		pose->set_chi( j, ligid_, ligandchis_[j] );
	}
}

///
/// create a reduced pose representation (to be used in minimization)
void
LigandConformer::to_minipose( core::pose::PoseOP pose, LigandConformer &minilig ) const {
	core::pose::PoseOP fullpose( new core::pose::Pose );
	to_pose( fullpose );

	minilig = *this;
	pose->clear();

	// sidechains
	for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
		if ( i == 1 ||  (movingscs_[i-1] == movingscs_[i]-1 && !fullpose->residue( movingscs_[i-1] ).is_upper_terminus() ) ) {
			pose->append_residue_by_bond( fullpose->residue( movingscs_[i] ) );
		} else {
			pose->append_residue_by_jump( fullpose->residue( movingscs_[i] ), pose->total_residue() );
		}
		minilig.movingscs_[i] = pose->total_residue();
	}

	// ligand
	pose->append_residue_by_jump( fullpose->residue( ligid_ ), pose->total_residue() );
	minilig.ligid_ = pose->total_residue();

	// root ligand on virtual if no residues to anchor jump
	// (so we have a jump to minimize)
	if ( pose->total_residue() == 1 ) {
		core::pose::addVirtualResAsRoot(*pose);
	}
}

///
/// update internal information from the reduced pose representation
void
LigandConformer::update_conf_from_minipose( core::pose::PoseCOP pose ) {
	core::pose::PoseOP fullpose( new core::pose::Pose );
	to_pose( fullpose );

	int resctr = 1;
	for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
		fullpose->replace_residue( movingscs_[i], pose->residue(resctr++), false );
	}
	fullpose->replace_residue( ligid_, pose->residue(resctr), false );

	update_conf( fullpose );
}

std::string
LigandConformer::to_string() const {
	using namespace ObjexxFCL::format;
	std::ostringstream oss;
	oss << "[[ rt: ";
	for ( core::Size i=1; i<=4; ++i ) {
		oss << F( 6, 3, rb_[i]) << " ";
	}
	for ( core::Size i=5; i<=7; ++i ) {
		oss << F( 6, 1, rb_[i]) << " ";
	}
	oss << "c: ";
	for ( core::Size i=1; i<=ligandchis_.size(); ++i ) {
		oss << F( 6, 1, ligandchis_[i]) << " ";
	}
	oss << " ]]";
	return oss.str();
}



void
LigandConformer::dump_pose( std::string pdbname ) const {
	core::pose::PoseOP pose( new core::pose::Pose() );
	to_pose( pose );
	pose->dump_pdb( pdbname );
}


// generate only the ligand residue based on this conformation
core::conformation::Residue
LigandConformer::ligand_residue( ) const {
	core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
	to_pose( pose );
	return (pose->residue(ligid_));
}

// generate only the ligand residue based on this conformation
core::conformation::Residue
LigandConformer::protein_residue( core::Size ires ) const {
	core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
	to_pose( pose );
	return (pose->residue( ires ));
}


// const access for distance calculation;
// ref_pose_ should be synced with rb,ligchi,and so on through update_conf
utility::vector1< core::Vector > const &
LigandConformer::ligand_xyz( ) {
	if ( !ligandxyz_synced_ ) {
		TR.Debug << "WARNING! Expensive ligand resynch" << std::endl; // it's not that expensive
		core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
		to_pose( pose );
		core::conformation::Residue const &lig( pose->residue(ligid_));
		ligandxyz_.resize( lig.nheavyatoms() );
		core::Vector com( 0.0 );
		for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
			ligandxyz_[iatm] = lig.xyz( iatm );
			com += lig.xyz( iatm );
		}
		com /= lig.nheavyatoms();

		// for radius of gyration
		rg_ = 0.0;
		for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
			rg_ += com.distance_squared( ligandxyz_[iatm] );
		}
		rg_ /= lig.nheavyatoms();
		rg_ = std::sqrt(rg_);

		ligandxyz_synced_ = true;
	}
	return ligandxyz_;
}

// the initial perturbation randomizes ligand conf
void
LigandConformer::randomize( core::Real transmax ) {
	ligand_xyz( ); // force synch

	numeric::Quaternion< core::Real > Q;
	numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
	numeric::R2quat( R, Q );
	rb_[1] = Q.w();
	rb_[2] = Q.x();
	rb_[3] = Q.y();
	rb_[4] = Q.z();

	core::Vector Taxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	core::Real Tlen = numeric::random::rg().uniform()*transmax;
	for ( core::Size k = 5; k <= 7; ++k ) {
		rb_[k] += Tlen*Taxis[k-5];
	}

	// ligand chis
	core::Size nligchi = ligandchis_.size();
	for ( core::Size j=1; j<=nligchi; ++j ) {
		core::Real angle_j = 360.0 * numeric::random::rg().uniform();
		ligandchis_[j] = angle_j;
	}

	ligandxyz_synced_ = false;
}

void
LigandConformer::assign_ligand_trans( core::Vector transv ) {
	for ( core::Size k = 5; k <= 7; ++k ) {
		rb_[k] = transv[k-5];
	}
}

// mutation
LigandConformer
mutate(LigandConformer const &l ) {
	using namespace ObjexxFCL::format;

	LigandConformer retval( l );
	std::string tag;

	// rigid body
	bool randomize_rt = (numeric::random::rg().uniform() < l.rtmutationRate_ );

	tag += " ";
	if ( randomize_rt ) {
		numeric::Quaternion< core::Real > Q;
		if ( l.rotmutWidth_ > 180.0 ) { // complete randomization
			numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
			numeric::R2quat( R, Q );
		} else {  // perturbation
			core::Vector Raxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
			core::Real angle = numeric::NumericTraits<core::Real>::deg2rad() * l.rotmutWidth_ * numeric::random::rg().gaussian();
			core::Real ca = cos(angle), sa=sin(angle);
			numeric::Quaternion< core::Real > Q1( ca, sa*Raxis[0], sa*Raxis[1], sa*Raxis[2]);
			Q = Q1*numeric::Quaternion< core::Real >(l.rb_[1],l.rb_[2],l.rb_[3],l.rb_[4]);
		}
		retval.rb_[1] = Q.w(); retval.rb_[2] = Q.x(); retval.rb_[3] = Q.y(); retval.rb_[4] = Q.z();

		core::Vector Taxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
		core::Real len = numeric::random::rg().uniform()*l.transmutWidth_;
		for ( core::Size k = 5; k <= 7; ++k ) {
			retval.rb_[k] += len*Taxis[k-5];
		}
		tag += "rt ";
	} else {
		// ligand chis
		core::Size nligchi = l.ligandchis_.size();
		core::Size ntors_changed(0);
		for ( core::Size j=1; j<=nligchi; ++j ) {
			if ( numeric::random::rg().uniform() < l.torsmutationRate_ ) {
				core::Real angle_j;
				if ( l.ligchimutWidth_ > 180.0 ) { // complete randomization
					angle_j = 360.0 * numeric::random::rg().uniform();
				} else {
					core::Real angleDel = l.ligchimutWidth_*(1.0 - 2.0*numeric::random::rg().uniform());
					angle_j = l.ligandchis_[j] + angleDel;
				}
				retval.ligandchis_[j] = angle_j;
				tag += "c"+utility::to_string(j)+" ";
				ntors_changed++;
			}
		}

		// make sure at least one torsion is mutated
		if ( ntors_changed == 0 && !randomize_rt && nligchi > 0 ) {
			core::Size j( numeric::random::rg().random_range(1,nligchi) );
			core::Real angle_j;
			if ( l.ligchimutWidth_ > 180.0 ) { // complete randomization
				angle_j = 360.0 * numeric::random::rg().uniform();
			} else {
				core::Real angleDel = l.ligchimutWidth_*(1.0 - 2.0*numeric::random::rg().uniform());
				angle_j = l.ligandchis_[j] + angleDel;
			}
			retval.ligandchis_[j] = angle_j;
			tag += "c"+utility::to_string(j)+" ";
		}
	}

	retval.set_generation_tag( tag );
	retval.ligandxyz_synced_ = false;

	return retval;
}

// crossover
LigandConformer
crossover(LigandConformer const &l1, LigandConformer const &l2) {
	using namespace ObjexxFCL::format;

	LigandConformer retval;
	std::string tag;

	// rigid body
	//   a) take rb as a single unit
	//   b) take sidechains from same gene we inherit rb (will be repacked though)
	bool takeRBfromIn1 = (numeric::random::rg().uniform() < 0.5 );
	if ( takeRBfromIn1 ) {
		retval = l1;
		tag += " rt:1 ";
	} else {
		retval = l2;
		tag += " rt:2 ";
	}

	// torsional
	core::Size nligchi = l1.ligandchis_.size();
	core::Size ntors_changed(0);
	tag += "c:";
	for ( core::Size j=1; j<=nligchi; ++j ) {
		bool takeChiJfromIn1 = (numeric::random::rg().uniform() < 0.5 );
		if ( takeChiJfromIn1 ) {
			retval.ligandchis_[j] = l1.ligandchis_[j];
			tag += "1";
		} else {
			ntors_changed++;
			retval.ligandchis_[j] = l2.ligandchis_[j];
			tag += "2";
		}
	}
	tag += " ";

	// make sure at least one torsion is mutated
	if ( ntors_changed == 0 && nligchi > 0 ) {
		core::Size j = core::Size(nligchi*numeric::random::rg().uniform())+1;
		retval.ligandchis_[j] = l2.ligandchis_[j];
		tag += "x"+utility::to_string(j)+" ";
	}

	// report distances
	//std::pair< core::Real, core::Real > internal_d = distance_internal( l1, l2 );
	//tag += " D"+utility::to_string(internal_d.first)+" "+utility::to_string(internal_d.second)+" ";

	retval.set_generation_tag( tag );

	retval.ligandxyz_synced_ = false;

	return retval;
}

LigandConformer
crossover_ft(LigandConformer const &l1, LigandConformer const &l2 ){
	// pick (randomly) one to be "master"
	//fd let's give every pose a shot to be master
	bool l1_is_master = true; // (numeric::random::rg().uniform() <= 0.5 );
	LigandConformer const &l_master = l1_is_master? l1 : l2;
	LigandConformer const &l_slave = l1_is_master? l2 : l1;

	LigandConformer retval(l_master);
	std::string tag;
	if ( l1_is_master ) tag = "master1 ";
	else tag = "master2 ";

	// pick ONE chi and set all downstream to l_slave
	core::Size nligchi = l1.ligandchis_.size();
	core::Size ichi_to_cross( numeric::random::rg().random_range(1,nligchi) );
	utility::vector1< core::Size > const &chis_down = l1.ligandchi_downstream_[ichi_to_cross];
	retval.ligandchis_[ichi_to_cross] = l_slave.ligandchis_[ichi_to_cross];
	tag += "flip:c"+utility::to_string(ichi_to_cross)+" ";
	for ( core::Size j = 1; j <= chis_down.size(); ++j ) {
		core::Size jchi = chis_down[j];
		tag += "down:c"+utility::to_string(jchi)+" ";
		retval.ligandchis_[jchi] = l_slave.ligandchis_[jchi];
	}

	retval.set_generation_tag( tag );
	retval.ligandxyz_synced_ = false;

	return retval;
}


core::Real
distance_slow( LigandConformer const &gene1, LigandConformer const &gene2 ){
	return core::scoring::automorphic_rmsd(
		gene1.ligand_residue(),
		gene2.ligand_residue(),
		false
	);
}

// this is MSD not RMSD!
core::Real
distance_fast( LigandConformer &gene1, LigandConformer &gene2 ) {
	core::Real d( 0.0 );
	core::Size n( 0 );
	utility::vector1< core::Vector >const &xyz1 = gene1.ligand_xyz();
	utility::vector1< core::Vector >const &xyz2 = gene2.ligand_xyz();
	for ( core::Size iatm = 1; iatm <= xyz1.size(); ++iatm ) {
		d += xyz1[iatm].distance_squared( xyz2[iatm] );
		n++;
	}

	return std::sqrt(d/n);
}

// return distance based on chi difference
// FD: NOTE THIS WILL NEED TO BE MODIFIED FOR DESIGN! (e.g. if nchi(gene1) != nchi(gene2) at some position)
std::pair< core::Real, core::Real >
distance_internal( LigandConformer const &gene1, LigandConformer const &gene2 ){
	std::pair< core::Real, core::Real > d;

	numeric::Quaternion< core::Real > q1 = gene1.quat();
	numeric::Quaternion< core::Real > q2 = gene2.quat();
	numeric::xyzVector< core::Real > t1 = gene1.trans();
	numeric::xyzVector< core::Real > t2 = gene2.trans();

	// more weight on rotation as ligand gets bigger
	// perhaps better estimation using radius of generation?
	//core::Real const r = gene1.ligand_residue().nbr_radius();
	core::Real const r = 0.5*(gene1.ligand_rg() + gene2.ligand_rg());

	numeric::Quaternion< core::Real > qdiff = q1.invert()*q2;
	d.first = t1.distance(t2) + r*qdiff.angle(); // angle is radian

	core::Size nligchi = gene1.get_ligandchis().size();
	core::Size ndof( nligchi );
	for ( core::Size ires = 1; ires <= gene1.moving_scs().size(); ++ires ) {
		ndof += gene1.get_protein_chis(ires).size();
	}
	//chidiff.resize(ndof);
	utility::vector1< core::Real > chidiff( ndof, 0.0 );

	for ( core::Size ichi = 1; ichi <= gene1.get_ligandchis().size(); ++ichi ) {
		core::Real diff( gene1.get_ligandchi(ichi) - gene2.get_ligandchi(ichi) );
		chidiff[ichi] = diff;
	}

	for ( core::Size ires = 1; ires <= gene1.moving_scs().size(); ++ires ) {
		utility::vector1< core::Real > const &reschi1 = gene1.get_protein_chis(ires);
		utility::vector1< core::Real > const &reschi2 = gene2.get_protein_chis(ires);
		for ( core::Size ichi = 1; ichi <= reschi1.size(); ++ichi ) {
			core::Real diff( reschi1[ichi] - reschi2[ichi] );
			chidiff[nligchi+ichi] = diff;
		}
	}

	for ( core::Size idof = 1; idof <= ndof; ++idof ) {
		d.second += std::abs(chidiff[idof]); // use Hamming instead of MSD or RMSD
	}

	return d;
}


}
}
}
