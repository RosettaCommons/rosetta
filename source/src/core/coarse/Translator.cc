// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ChemicalManager.cc
/// @brief  Chemical manager class
/// @author Oliver Lange (olange@u.washington.edu)

#include <core/coarse/Translator.hh>
// AUTO-REMOVED #include <core/coarse/TranslatorSet.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/CoarseRotamer.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

#include <ObjexxFCL/string.functions.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/DOF_ID.hh>

#include <numeric/NumericTraits.hh>

using namespace core;
using namespace coarse;
using namespace std;
using namespace conformation;
using namespace chemical;
using namespace ObjexxFCL;

using namespace pack::dunbrack;

/// must be a better place for this,  already exists in ResidueType.cc!
inline
std::string
strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
}

void Translator::pretty_print(std::ostream& os) const {
	int bead_nr=0;
	os << "==============================" << endl;
	os << "Coarse-Graining for residue " << fine_res_type_->name() << endl;
	os << "-------------------------------" << endl;
	for (BeadList::const_iterator it=beads_.begin(), eit=beads_.end(); it!=eit;++it,++bead_nr) {
		os << bead_names_[bead_nr] << endl;
		for (AtomList::const_iterator ait=it->begin(), eait=it->end(); ait!=eait; ++ait) {
			os << ait->name_ << " ";
		};
		os << endl;
	};
	os << "==============================" << endl;
}
int
Translator::coarse_nchi() const {
	return coarse_res_type_->nchi();
}

std::string const &
Translator::name() const {
	return coarse_res_type_->name();
}


void
Translator::add_atom(AtomList &list, const ResidueType &res,  int pos) {
	add_atom(list,res,res.atom_name(pos));
}

void
Translator::add_atom(AtomList &list, const ResidueType &res,  const Rule::AtomToken &atom) {
	/* careful: although the ResidueType understands atomnames like "CB" its function atom_name yields " CB ",
		 so we store only the "stripped" version of atom_names in our list */
	int pos = -1;
	if (res.has(atom)) {
		pos=res.atom_index(atom);
		if (map_atom_to_bead(atom)<0) {
			list.push_back(BeadAtom(strip_whitespace(atom),1.0));
			for (uint na=res.attached_H_begin(pos);na<=res.attached_H_end(pos);na++) {
	if (map_atom_to_bead(res.atom_name(na))<0) {
		list.push_back(BeadAtom(strip_whitespace(res.atom_name(na)),0.0));
	}
			};
		};
	}
}


int Translator::map_atom_to_bead(std::string atom) const {
	// returns bead-position of atom at pos in full_atom
	// returns -1 if atom not found
	bool bfound=false;
	AtomList::const_iterator ait;
	int bead_nr=0;
	//WEIGHT: Find will not work anymore
	for (BeadList::const_iterator it=beads_.begin(), eit=beads_.end(); it!=eit && !bfound;++it,bead_nr++) {
		bfound=((ait=find(it->begin(),it->end(),strip_whitespace(atom)))!=it->end());
	};

	return bfound ? bead_nr : -1;
}


void
Translator::add_remaining_sidechain(AtomList &list, const ResidueType &res) {
	for (uint pos=res.first_sidechain_atom(); pos<=res.nheavyatoms(); ++pos) {
		add_atom(list,res,pos);
	};
}

void
Translator::add_all_remaining(AtomList &list, const ResidueType &res) {
	for (uint pos=1; pos<=res.nheavyatoms(); ++pos) {
		add_atom(list,res,pos);
	};
}

Translator::Translator(const RuleSet &rules, ResidueTypeCOP fine_res_ptr,  ResidueTypeAP coarse_res_ptr) {
	coarse_res_type_=coarse_res_ptr;
	coarse_res_type_->set_translator(TranslatorCAP(this));
	fine_res_type_=fine_res_ptr;
	ResidueType const &fine_res=(*fine_res_ptr);
	ResidueType const &coarse_res=(*coarse_res_ptr);
debug_assert( coarse_res.name() == fine_res.name() );

	chemical::AA aa = fine_res.aa();
	beads_.clear();
	beads_.resize(1); //make space for "full_atom" at pos 0;
	bead_names_.clear(); bead_names_.push_back(Rule::FULL_ATOM);

	RuleCOP rule = rules[aa];
	for (Rule::ConstBeadIterator it=rule->begin(), eit=rule->end(); it!=eit; ++it) {
		AtomList *pal;
		if (it->first==Rule::FULL_ATOM) {
			pal=&beads_.front();
		} else {
			// start "ingredient list" for the current coarse-bead
			beads_.push_back(AtomList());

			// put bead name into bead_names_   Full-Atom is already at position 0
			bead_names_.push_back(it->first);
			pal=&beads_.back();
		}


		cerr << "add bead " << it->first << endl;

		AtomList &al=*pal;
		for (Rule::ConstTokenIterator ait=it->second.begin(), eait=it->second.end(); ait!=eait; ++ait) {
			Rule::AtomToken token=*ait;
			if (token==Rule::REST_SIDECHAIN) {
	add_remaining_sidechain(al,fine_res);
			} else if (token==Rule::REST_ALL) {
	add_all_remaining(al,fine_res);
			} else
	add_atom(al,fine_res,token);
		}
	}

	// check for FULL_ATOM that both residues contain the atom name
	AtomList &full_atom=beads_.front();
	for (AtomList::const_iterator ait=full_atom.begin(),eait=full_atom.end(); ait!=eait; ++ait) {
		cerr << ait->name_ << endl;
	debug_assert( fine_res.has(ait->name_) && coarse_res.has(ait->name_) );
	};

	//fix bond lengths, angles, etc.
	//  ResidueType const * ptr = &(*coarse_res_ptr);
	//  ResidueType *non_cons = const_cast< ResidueType* > (ptr);
	//fix_coarsetype_geometry(ResidueTypeOP( non_cons ));
	fix_coarsetype_geometry( coarse_res_ptr );
}

void Translator::fix_coarsetype_geometry(ResidueTypeAP c_type) {
	//creates fine-residue instance, coarse grains it and reads out the coarse-geometry

	//create fine residue instance
	ResidueOP frsd = ResidueFactory::create_residue(*fine_res_type_);

	//coarse grain
	ResidueOP crsd = coarsify(*frsd);

	for (Size ai=1;ai<=crsd->natoms();ai++) {
		//    c_type->atom(ai).xyz( crsd->atom(ai).xyz() );
		// change in interface of residuetype
		c_type->set_xyz( ai, crsd->atom(ai).xyz() );
	}
}


ResidueOP Translator::coarsify(const Residue &fine) const {
	ResidueOP new_rsd( ResidueFactory::create_residue( *coarse_res_type_ ) );
	int bead_nr=0;
	for (BeadList::const_iterator it=beads_.begin(), eit=beads_.end(); it!=eit;++it,++bead_nr) {
		core::Vector cen(0.0,0.0,0.0);
		Real sum_weight=0;
		for (AtomList::const_iterator ait=it->begin(), eait=it->end(); ait!=eait; ++ait) {
			if (bead_nr==0) {
	//full_atom copy them
	//	cerr << "copy atom " << *ait << fine.atom(*ait).xyz()[1] << endl;
	new_rsd->atom(ait->name_).xyz(fine.atom(ait->name_).xyz());
	//	cout << "atom type of bead " << bead_names_[bead_nr] << " " << ait->name_ << " is " << new_rsd->atom(ait->name_).type() << endl;
			} else {
	//coarse_grain compute mean
	cen+=fine.atom(ait->name_).xyz()*ait->weight_;
	sum_weight+=ait->weight_;
	//	Vector tmp=cen/natom;
	//	cerr << "add bead atom " << *ait << fine.atom(*ait).xyz()[1] <<' ';
	//	cerr << tmp[1] << endl;

			};
		};
		if (bead_nr>0) {
			// coarse grained set bead coordinates
			if (sum_weight>0.0) {
	new_rsd->atom(bead_names_[bead_nr]).xyz(cen/sum_weight);
	//	cout << "atom type of bead " << bead_names_[bead_nr] << " is " << new_rsd->atom(bead_names_[bead_nr]).type() << endl;
			} else {
	cerr << "WARNING: empty bead " << bead_names_[bead_nr] << " in residue " << coarse_res_type_->name() << endl;
			};
		}
	}
	return new_rsd;
}

// some helper functions for the coarsification of SingleResidueDunbrackLibraries

bool match_mask(RotVector const& mask, int nchi, DunbrackRotamer< FOUR, Real > const & rotamer );

bool update_mask(RotVector& mask,int nchi,RotVector const &max_bins);

pose::PoseOP create_rotamer(
	Translator const& map,
	ResidueTypeCOP fine_res_type,
	DunbrackRotamer< FOUR, Real > const& rotamer
);

int find_most_frequent_rotamer(
	utility::vector1< DunbrackRotamer< FOUR > > const & fine_rotamers,
	RotVector const& mask,
	Size nchi,
	Real &pnew
);

void coarse_rotamer(
	Translator const& map,
	ResidueTypeCOP fine_res_type,
	ResidueTypeCOP coarse_res_type,
	DunbrackRotamer< FOUR, Real > const& rotamer,
	Size nchi,
	ChiVector &chi,
	AngleVector &angle
);

void average_rotamers(
	Translator const& map,
	ResidueTypeCOP fine_res_type,
	ResidueTypeCOP coarse_res_type,
	utility::vector1< DunbrackRotamer< FOUR > > const & fine_rotamers,
	RotVector const& mask,
	Size nchi,
	ChiVector &chi_mean,
	ChiVector &chi_std,
	AngleVector &angle_mean,
	AngleVector &angle_std
);


CoarseRotamerSetOP
Translator::coarsify(
	utility::vector1< DunbrackRotamer< FOUR > > const & fine_rotamers
) const
{

	using namespace scoring;
	using namespace pack::dunbrack;
	using namespace conformation;
	using namespace std;
	bool bAverage( false ); //use most frequent chi/angle for coarse rotamer

	CoarseRotamerSetOP coarse_rotset = new CoarseRotamerSet;
debug_assert(coarse_rotset);
	//coarse_rotset->tag="coarse";
	Size nchi=coarse_res_type_->nchi();
debug_assert(nchi==1); //this is only correct in the current version... a debuggin' assert

	RotVector mask(4);
	RotVector max_bins(nchi);

	{  //set mask to first rotamer_id
		for (Size chi=1;chi<=nchi; chi++) {
			mask[chi]=1;
			max_bins[chi]=3; // is that always so ?
		};
		for (Size chi=nchi+1;chi<=4; chi++) {
			mask[chi]=0;
		};
	}

	do { //while( update_mask(...) )
		Real pnew;
		int rot_freq = find_most_frequent_rotamer( fine_rotamers, mask, nchi, pnew );
	//debug_assert ( nchi == fine_rotamers[ rot_freq ].nchi_aa() ); //case of inequality might occur and then we need to think about it...
		ChiVector new_chi( 4, 0.0 ); //initialize with 0
		AngleVector new_angle( 4, 0.0 );
		DunbrackRotamer< FOUR, Real > fine_high_res = increase_rotamer_precision( fine_rotamers[ rot_freq ] );
		coarse_rotamer( *this, fine_res_type_, coarse_res_type_, fine_high_res, nchi, new_chi, new_angle );

		ChiVector chi_mean( 4, 0.0 ); //initialize with 0
		ChiVector chi_std( 4, 0.0 );
		AngleVector angle_mean( 4, 0.0 );
		AngleVector angle_std( 4, 0.0 );
		average_rotamers( *this, fine_res_type_, coarse_res_type_, fine_rotamers, mask, nchi, chi_mean, chi_std, angle_mean, angle_std );

		CoarseRotamerOP coarse_rot;
		if ( bAverage ) {
			coarse_rot = new CoarseRotamer( pnew, nchi, mask, chi_mean, chi_std, angle_mean, angle_std );
		} else {
			coarse_rot = new CoarseRotamer( pnew, nchi, mask, new_chi, chi_std, new_angle, angle_std );
		}
		coarse_rotset->push_back(coarse_rot);

	} while ( update_mask( mask, nchi, max_bins ) );
	// You could use a lexicographical iterator here...

	return coarse_rotset;
}

SingleResidueRotamerLibraryCOP
Translator::get_RotamerLibrary() const {
	std::cerr << "Translator::get_RotamerLibrary() for " << name() << std::endl;
	//SingleResidueRotamerLibraryCAP rotlib (fine_res_type_->get_RotamerLibrary());
	//if (rotlib) return rotlib->coarsify(*this);
	//else
	return NULL;
}


bool update_mask(RotVector& mask,int nchi,RotVector const &max_bins) {
	int chi=1;
	while (++mask[chi]>max_bins[chi]) {
		mask[chi++]=1;
		if (chi>nchi) return false;
	}
	return true;
}

// does a coarse rotamer 13xx match to the fine grain 1332 ?
// mask: 13xx
// nchi: length of mask (e.g, 2)
// rotid: 1321
bool match_mask(RotVector const& mask, int nchi, DunbrackRotamer< FOUR, Real > const & rotamer ) {
	int chi=1;
	while ( mask[chi] == rotamer.rotwell(chi) ) {
		chi++;
		if ( chi>nchi ) return true;
	};
	return false;
}

inline Real sqr ( Real x ) {
	return x*x;
}
inline Real sqr3 ( Real x ) {
	return x*x*x;
}

pose::PoseOP
create_rotamer(
	Translator const& map,
	ResidueTypeCOP fine_res_type,
	DunbrackRotamer< FOUR, Real > const & rotamer
)
{
	ResidueOP fine_res = ResidueFactory::create_residue( *fine_res_type );
	//real-world representation of rotamer
	for ( Size jj = 1; jj <= std::min( Size( 4 ), fine_res_type->nchi()); ++jj ) { // apl needs info from DunLib
		fine_res->set_chi( jj, rotamer.chi_mean( jj ) );
	}

	//coarse representation of rotamer
	ResidueOP coarse_res = map.coarsify(*fine_res);
	pose::PoseOP pose = new pose::Pose;
	pose->clear();
	pose->append_residue_by_bond(*coarse_res);
	pose->fold_tree( kinematics::FoldTree( pose->total_residue() ) );
	// @phil: might be nice to have a method update_internals() in Residue
					//      other thing: I think it is counter-intuitive that pose clones the residue
					//                that is passed in.
					//           might have good reasons -- might be a remnant of the old design, where you
					//            had to make clones from the "ideal" residues.
	return pose;
}


void coarse_rotamer(
	Translator const& map,
	ResidueTypeCOP fine_res_type,
	ResidueTypeCOP coarse_res_type,
	DunbrackRotamer< FOUR, Real > const & rotamer,
	Size nchi,
	ChiVector &chi,
	AngleVector &angle
) {
	const Real Pi = numeric::NumericTraits<Real>::pi(); 

	pose::PoseOP pose ( create_rotamer(map, fine_res_type, rotamer) );
	for ( Size jj = 1; jj<=nchi; ++jj) {
		chemical::AtomIndices chi_atoms = coarse_res_type->chi_atoms( jj );
		const Size seqpos ( 1 ); //only one residue in artifical pose
		chi[jj] =  pose->chi( jj, seqpos);
		angle[jj] = 180.0/Pi*pose->dof( id::DOF_ID( id::AtomID( chi_atoms[4], seqpos ),id::THETA) );
	}
}


void average_rotamers(
	Translator const& map,
	ResidueTypeCOP fine_res_type,
	ResidueTypeCOP coarse_res_type,
	utility::vector1< DunbrackRotamer< FOUR > > const & fine_rotamers,
	RotVector const& mask,
	Size nchi,
	ChiVector &chi_mean,
	ChiVector &chi_std,
	AngleVector &angle_mean,
	AngleVector &angle_std
)
{
	const Real Pi = numeric::NumericTraits<Real>::pi(); 
	//run thru all rotamers to get sdev's

	Real pnew( 0.0 );
	FArray2D< Real > xy_mean_chi( nchi, 2, 0.0 );
	FArray2D< Real > xy_mean_angle( nchi, 2, 0.0 );

	for ( Size i=1; i<=fine_rotamers.size(); i++) {
		DunbrackRotamer< FOUR, Real > const rotamer = increase_rotamer_precision( fine_rotamers[i] );
		if ( match_mask( mask, nchi, rotamer ) ) {

			Real p = rotamer.rotamer_probability();
			ChiVector newChi( 4, 0.0 );
			ChiVector newAngle( 4, 0.0 );
			coarse_rotamer( map, fine_res_type, coarse_res_type, rotamer, nchi, newChi, newAngle );
			for ( Size jj = 1; jj <= nchi; ++jj) {
				Real r_chi = 1-sqr3(rotamer.chi_sd(jj)/180.0*Pi)/2.0;
				xy_mean_chi( jj, 1)+= p * std::sin( newChi[jj] ) * r_chi;
				xy_mean_chi( jj, 2)+= p * std::cos( newChi[jj] ) * r_chi;

				xy_mean_angle( jj, 1)+= p * std::sin( newAngle[jj] );
				xy_mean_angle( jj, 2)+= p * std::cos( newAngle[jj] );
			};
			pnew+=p;
		} //match mask
	} // for loop

	for ( Size jj = 1; jj<=nchi; ++jj) {
		chi_mean[jj]=std::atan2( xy_mean_chi(jj,1), xy_mean_chi(jj,2) ) * 180/Pi;
		angle_mean[jj]=std::atan2( xy_mean_angle(jj,1), xy_mean_angle(jj,2) ) * 180/Pi;
		if (pnew>1e-7) {
			Real r_chi = std::sqrt( sqr( xy_mean_chi( jj, 1 ) )+sqr( xy_mean_chi( jj, 2) ) )/pnew;
			chi_std[jj] = pow( 2.0 * ( 1.0-r_chi ), 1/3.0 ) * 180 / Pi;
			Real r_ang = std::sqrt( sqr( xy_mean_angle( jj, 1 ) )+sqr( xy_mean_angle( jj, 2) ) )/pnew;
			angle_std[jj] = pow( 2.0 * (1.0 - r_ang), 1/3.0 ) * 180 / Pi;
		} else {
			chi_std[jj]=10; // it doesn't matter just avoid a nan...
			angle_std[jj]=10;
		}
	}
}

int find_most_frequent_rotamer(
	utility::vector1< DunbrackRotamer< FOUR > > const & fine_rotamers,
	RotVector const& mask,
	Size nchi,
	Real & pnew
)
{	// find most frequent rotamer
	using namespace pack::dunbrack;

	int rot_freq ( -1 );
	Real pmax ( -1.0 );
	pnew =  0.0;
	for (Size i=1; i<=fine_rotamers.size(); i++) {
		DunbrackRotamer< FOUR, Real > const rotamer = increase_rotamer_precision( fine_rotamers[i] );
		if ( match_mask( mask, nchi, rotamer ) ) {
			Real p = rotamer.rotamer_probability();
			if ( p > pmax ) {
				pmax = p;
				rot_freq=i;
				pnew+=p;
			}
		}
	}
	return rot_freq;
}

