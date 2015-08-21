// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Mover.cc
///
/// @brief
/// @author
#include <devel/inv_kin_lig_loop_design/Mover.hh>
#include <devel/inv_kin_lig_loop_design/Fragment.hh>
//#include <devel/LoopDesign/Normal.hh>

// #include <Util/Util.hh>
// #include <Util/Fragment.hh>
// #include <Util/Normal.hh>
// #include <Util/KinPacker.hh>
// #include <Util/KinRotamer.hh>
// #include <Util/ErgRamaInfo.hh>
// #include <Util/Abng.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

/// Numeric headers
#include <numeric/random/random.hh>

/// C++ header
#include <cmath>

#include <utility/vector1.hh>


#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)


namespace devel {

namespace inv_kin_lig_loop_design {


Mover::Mover(core::pose::Pose* pose) :
	pose(pose),
	scratch_( core::pack::dunbrack::RotamerLibraryScratchSpaceOP( new core::pack::dunbrack::RotamerLibraryScratchSpace() ) )
{} // Mover::Mover

/// @details DANGER DANGER DANGER Pose is not deallocated
Mover::~Mover()
{}

void Mover::setPhiPsi(const int seqpos, const double phi, const double psi ) const {
	pose->set_phi( seqpos, phi );
	pose->set_psi( seqpos, psi );

	//cout << "set pose " << seqpos << " phipsi to " << phi << "," << psi << endl;


} // Mover::setPhiPsi

void Mover::setFragment(const int seqpos, const int frag_len, const Fragment::Entry& e) const { // XXX stupid interface really

	for ( int i = 0; i < frag_len; ++i ) {

		const double& phi = RAD2DEG*e.vResEntries[i].phi;
		const double& psi = RAD2DEG*e.vResEntries[i].psi;
		//const double& ohm = e.vResEntries[i].ohm;

		setPhiPsi(seqpos+i,phi,psi);

	} // i

} // Mover::setFragment

//void Mover::insertFragment(const char chain_id, const int res_num_a, const int res_num_b, const int frag_len ) const {

void Mover::insertFragment(const int seqpos_a, const int seqpos_b, const int frag_len ) const {

	// Get random res_num,res
	//const int seqpos = seqpos_a + random() % (seqpos_b - seqpos_a + 1 - frag_len + 1);// * ( res_num_b - frag_len - res_num_a);
	const int seqpos = numeric::random::random_range(seqpos_a,seqpos_b);
	assert( seqpos_a <= seqpos && seqpos <= seqpos_b );

	// Get random entry
	const Fragment::File* frag_file = getFragmentFile(seqpos,frag_len);
	const Fragment::Entry& e = frag_file->getEntry(frag_len);

	// Set the fragment
	setFragment(seqpos,frag_len,e);

} // Mover::insertFragment

void Mover::randomFragments(const int seqpos_a, const int seqpos_b, const int frag_len ) const {

	// pick offset

	int offset = numeric::random::random_range(0,frag_len-1);

	for ( int seqpos = seqpos_a + offset; seqpos <= seqpos_b - frag_len + 1; seqpos += frag_len ) {

		const Fragment::File* frag_file = getFragmentFile(seqpos,frag_len);
		const Fragment::Entry& e = frag_file->getEntry(frag_len);

		setFragment(seqpos,frag_len,e);
	} // res_num

} // Mover::randomFragments

void Mover::smallMoves(const int seqpos_a, const int seqpos_b, const double sigma) const {

	//Util::Normal Z(0,sigma);
	for ( int i = seqpos_a; i <= seqpos_b; ++i ) {
		//Peptide* p = const_cast<Peptide*>(getTube()->getResidue(chain_id,i)->get<Peptide>());

		const double old_phi = pose->phi(i);
		const double old_psi = pose->psi(i);

		const double z = sigma * numeric::random::gaussian();

		const double new_phi = old_phi + z;
		const double new_psi = old_psi - z;

		//cout << "Mover::smallMoves - " << REPORT4(RAD2DEG*old_phi,RAD2DEG*new_phi,RAD2DEG*old_psi,RAD2DEG*new_psi) << endl;
		//cout << "Mover::smallMoves - " << REPORT1(z) << endl;

		setPhiPsi(i,new_phi,new_psi);

	} // i
} // small_moves

//   const vector<Kin::Rotamer>& Mover::getRotamers(const char chain_id, const int res_num ) const {
//     //const vector<Rotamer> getRotamers(const Pdb::ResType& res_type, const double phi, const double psi, const int nb) const;

//     const Residue* r = getTube()->getResidue(Pdb::ResID(chain_id,res_num));
//     assert( r != 0 );
//     const Peptide* p = r->get<Peptide>();
//     assert( p != 0 );

//     if( mvRot.find(r) == mvRot.end() ) {
//       Kin::Packer packer;
//       //       packer.setEx1(true); // XXX don't know whether this is really useful or not
//       //       packer.setEx2(true);
//       //       packer.setEx1Aro(true);
//       mvRot[r] = packer.getRotamers(r); //->getResType(),p->getPhi(),p->getPsi(),40);
//       cout << "Mover::getRotamers - packer returned " << REPORT1(mvRot[r].size()) << endl;
//     }

//     return mvRot[r];

//   } // init_rotamers

//   void Mover::randomRotamer(const Residue* r) const {
//     randomRotamer(r->getChainID(),r->getResNum());
//   } // Mover::randomRotamer

void Mover::randomRotamer(const int seqpos) const {

	core::chemical::ResidueType const& res_type = pose->residue_type(seqpos);
	core::pack::dunbrack::ChiVector chivector( res_type.nchi() );

	cout << "inserting a random rotamer into " << res_type.name() << ":" << seqpos << endl;

	core::pack::rotamers::SingleResidueRotamerLibraryCOP rotamer_library = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( res_type );
	if ( ! rotamer_library ) return;

	core::pack::dunbrack::SingleResidueDunbrackLibraryCOP dun_rotlib =
		utility::pointer::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >( rotamer_library );
	//assert( dun_rotlib ); // WHOA! Where's the guarantee here?
	if ( ! dun_rotlib ) return; // Is this right? maybe we would want random rotamers for ligands?

	dun_rotlib->assign_random_rotamer_with_bias(
		pose->residue( seqpos ), *pose, *scratch_, numeric::random::rg(),
		chivector, true );

	for ( size_t k = 1; k <= res_type.nchi(); ++k ) {
		//res.set_chi( k, rotamer.chi_mean(k) + Z(0,1) * rotamer.chi_sd(k) );
		//pose->set_chi( k, seqpos, rotamer.chi_mean(k) + numeric::random::gaussian() * rotamer.chi_sd(k) );
		pose->set_chi( k, seqpos, chivector[ k ] );
	}
} // Mover::randomRotamer

// ____________________ redoing fragment insertions ____________________

//   void Mover::setFragmentFile(const Residue* r, const int len, const Fragment::File* frag_file) const {
//     mFragmentFiles[make_pair(r,len)] = frag_file;
//   } // Mover::setFragmentFile

namespace { /// What namespace is this in?  Is this in the global namespace?

const string get_ss(core::pose::Pose const* pose, int const seqpos, int const len ) {
	string rval;
	for ( int i = 0; i < len; ++i ) {
		rval.push_back(pose->secstruct(seqpos+i));
	}
	return rval;
} // get_ss

} // namespace

const Fragment::File* Mover::getFragmentFile(int seqpos, const int len) const {

	map<pair<int,int>,const Fragment::File*>::const_iterator i = mFragmentFiles.find(make_pair(seqpos,len));

	if ( i == mFragmentFiles.end() ) {

		// let's use two different ways of doing things

		const Fragment::File* frag_file = 0;

		if ( len == 1 ) {
			const char ss = pose->secstruct(seqpos);
			cout << "Mover::getFragmentFile - getting frag file '" << ss << "' for " << seqpos << endl;
			frag_file = Librarian::getFragmentFile(ss);
		} else if ( len == 3 ) {
			const string ss = get_ss(pose,seqpos,len);
			cout << "Mover::getFragmentFile - getting frag file '" << ss << "' for " << seqpos << endl;
			frag_file = Librarian::getFragmentFile(ss);
		} else {
			cout << "Mover::getFragmentFile - support for fragments of length " << len << " is not yet coded" << endl;
			assert( false );
		}

		assert( frag_file != 0 );
		mFragmentFiles[make_pair(seqpos,len)] = frag_file;
		return frag_file;
	} else { // i
		return i->second;
	}

} // Mover::getFragmentFile

//     void Mover::clearFragmentFiles() const {
//       mFragmentFiles.clear();
//     } // Mover::clearFragmentFiles

//   void Mover::setExtended(const Residue* a, const Residue* b) const {
//     static const double DEFAULT_PHI_E = -135*DEG2RAD;
//     static const double DEFAULT_PSI_E =  135*DEG2RAD;
//     //static const double DEFAULT_OHM_E =  180*DEG2RAD;

//     for( const Residue* r = a; r != b->getNext(); r = r->getNext() ) {
//       if( Peptide* p = const_cast<Peptide*>(r->get<Peptide>()) ) {
// setPhiPsi(p,DEFAULT_PHI_E,DEFAULT_PSI_E);
//       } // p
//     } // r

//   } // Mover::setExtended

//   void Mover::setHelix(const Residue* a, const Residue* b) const {
//     static const double DEFAULT_PHI_H = -47*DEG2RAD;
//     static const double DEFAULT_PSI_H =  57*DEG2RAD;
//     //static const double DEFAULT_OHM_H =  180*DEG2RAD;

//     for( const Residue* r = a; r != b->getNext(); r = r->getNext() ) {
//       if( Peptide* p = const_cast<Peptide*>(r->get<Peptide>()) ) {
// setPhiPsi(p,DEFAULT_PHI_H,DEFAULT_PSI_H);
//       } // p
//     } // r

//   } // Mover::setHelix

} // namespace LoopDesign

} // namespace devel
