// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GoapEnergy.cc
/// @brief  C++ implementaion of GOAP(Generalized Orientation-dependent, All-atom statistical Potential)
///         by Zhou H & Skolnick J, Biophys J 2011, 101(8):2043-52.
/// @author Hahnbeom Park

// Unit headers
#include <core/scoring/methods/GoapEnergy.hh>
#include <core/scoring/methods/GoapEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <basic/database/open.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

namespace core{
namespace scoring{
namespace methods{

static basic::Tracer TR("core.scoring.GoapEnergy");

//////////////////////
/// EnergyMethod Creator
methods::EnergyMethodOP
GoapEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new GoapEnergy( options );
}

ScoreTypes
GoapEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( goap );
	sts.push_back( goap_dist );
	sts.push_back( goap_angle );
	return sts;
}

//////////////////////////////////////////////////////
///////////// Start GoapRsdType

GoapRsdType::GoapRsdType(){}
GoapRsdType::~GoapRsdType(){}

void
GoapRsdType::setup_rsdtype( chemical::ResidueTypeCOP rsd ){

  Size const natom( rsd->natoms() );

	natom_ = natom;
  is_using_.resize( natom, false );
	atmid_.resize( natom, 0 );
  atmname_using_.resize( 0 );
  connected_by_twobonds_.resize( natom, false );
  i2_.resize( natom, 0 );
  i3_.resize( natom, 0 );
	name_ = rsd->name();
}

void 
GoapRsdType::setup_connectivity( chemical::ResidueType const &rsd )
{
  // Root and branch should be defined prior to this function

  for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
		Size i2, i3;
    bool connected_by_twobonds( false );
    std::string const atmname( rsd.atom_name( iatm ) );

    if( atmname.compare(" N  ") == 0 ){
      i2 = 999; // this means to search prvC
      i3 = rsd.atom_index(" CA ");
      connected_by_twobonds = true;

    } else if( atmname.compare(" CA ") == 0 ){
      i2 = rsd.atom_index(" N  ");
      i3 = rsd.atom_index(" C  ");
      connected_by_twobonds = true;

    } else if( atmname.compare(" C  ") == 0 ){
      i2 = rsd.atom_index(" CA ");
      i3 = rsd.atom_index(" O  ");
      connected_by_twobonds = true;

    } else if( atmname.compare(" O  ") == 0 ){
      i2 = rsd.atom_index(" C  ");
      i3 = rsd.atom_index(" CA "); 

    } else {
      i2 = root_atom_[iatm]; 
      i3 = branch_atom_[iatm];
      if( i3 > 0 ){ // If there is branch of this atom
				connected_by_twobonds = true;

			} else { // otherwise take root of the root connected by angle as i3
				i3 = angle_atom_[iatm];
				connected_by_twobonds = false;
			}
    }

    i2_[iatm] = i2;
    i3_[iatm] = i3;
    connected_by_twobonds_[iatm] = connected_by_twobonds;
  }
}

////////////// End GoapRsdType
//////////////////////////////////////////////////////


//////////////////////////////////////////////////////
///////////// Start GoapEnergy

GoapEnergy::GoapEnergy( GoapEnergy const & src ):
  parent( src ),
  max_dis_ ( src.max_dis_ ),
	continuous_ ( false )
{
  set_default();
  read_Goap_parameters();
}

GoapEnergy::GoapEnergy( EnergyMethodOptions const & ):
		parent( new GoapEnergyCreator ),
		continuous_ ( false )
{
	set_default();
	read_Goap_parameters();
}

GoapEnergy::~GoapEnergy(){}

//////////////////////////////////////////////////////
///////////// Setup functions

void 
GoapEnergy::set_default(){

  N_ATMTYPES = 167;
  N_DISTANCE_BINS = 20;
  N_ANGLE1_BINS = 5;
  N_ANGLE2_BINS = 12;

	// adjustable options
  //max_dis_ = option[ scoring::goap_max_dis ]();
	//MIN_SEQ_SEPARATIONS = option[ scoring::goap_min_seq_separation ]();
	max_dis_ = 15.0;
	MIN_SEQ_SEPARATIONS = 7;
}

void
GoapEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & )  const
{
	//std::cout << "Goap setup_for_scoring..." << std::endl;

  // Initialize dipole vectors given pose
  // These arrays are set as mutable to pass "const" declaration
  xn_.resize( pose.total_residue() );
  xd_.resize( pose.total_residue() );
	eval_res_.resize( pose.total_residue(), true );
  Vector const I( 0.0 );

  // Assign dipole vectors for current pose
  for( Size ires = 1; ires <= pose.total_residue(); ++ires ){

    conformation::Residue const &rsd = pose.residue( ires );
		//chemical::AA const &aa = pose.residue(ires).aa();
    chemical::ResidueType const &rsdtype = pose.residue(ires).type();

    // Energy won't work for terminus residues! this might change in the future,
    // but will require additional care for connectivity :)
		//std::cout << "ires " << ires << " " << rsdtype.aa() << " " << rsdtype.is_protein() << std::endl;
    if( rsdtype.is_terminus() || rsdtype.aa() > chemical::num_canonical_aas ){
			eval_res_[ires] = false;
			continue;
		}

    xn_[ires].resize( rsd.natoms(), I );
    xd_[ires].resize( rsd.natoms(), I );

    GoapRsdTypeMap::const_iterator it = rsdtypemap_.find( rsdtype.name3() );
		if( it == rsdtypemap_.end() ) continue;

    GoapRsdTypeOP goaptype = it->second;

    for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
      if( !goaptype->is_using( iatm ) ) continue;

      Vector xn1, xd1;
      bool result = calculate_dipoles( pose, rsd, goaptype, iatm, xn1, xd1 );
			if( !result ){
				eval_res_[ires] = false;
				TR.Debug << "Warning: Skipping " << rsd.seqpos() << " due to dipole setup failure for atom " << iatm << std::endl;
				break;
			}

			//std::cout << "xn1, xd1: " << std::setw(4) << rsd.name() << std::setw(5) << rsd.atom_name(iatm);
			//printf("%4d %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
			//			 ires, iatm, xn1[0], xn1[1], xn1[2], xd1[0], xd1[1], xd1[2]);
      xn_[ires][iatm] = xn1;
      xd_[ires][iatm] = xd1;
    }
  }

	//std::cout << "end setup!" << std::endl;
}

void 
GoapEnergy::read_Goap_parameters()
{
  // First, read atom definition and connectivity
  std::string const connection_file = 
    +"/scoring/score_functions/goap/angle_definitions.dat";
  
  read_angle_definitions( connection_file );

  std::string const distance_file = 
    +"/scoring/score_functions/goap/dist_table.dat.gz";
  std::string const angle_file = 
    +"/scoring/score_functions/goap/angle_table.dat.gz";

  read_potential_values( distance_file, angle_file );
}

void
GoapEnergy::read_angle_definitions( std::string const connection_file )
{
	Size i1;

	chemical::ResidueTypeSetCAP rsdtypeset =
	chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	utility::io::izstream instream;
	std::string line;

	TR << "Reading angle definition." << std::endl;
	
	basic::database::open( instream, connection_file );

	Size read_type( 0 );

	Size i_id( 0 );

	while( instream ){
		getline( instream, line );
			std::string s1("");
		std::istringstream linestream( line );
		linestream >> s1;

		if( s1.compare( "#ATOM" ) == 0 ){
			read_type = 1;
			continue;
		} else if( s1.compare( "#BOND" ) == 0 ){
			read_type = 2;
			continue;
		} else if( s1.compare( "#ANGLE" ) == 0 ){
			read_type = 3;
			continue;
			} else if( s1.compare( "#END" ) == 0 ){
				break;
		}

		chemical::ResidueType const &rsdtype = rsdtypeset->name_map( s1 );

			GoapRsdTypeMap::const_iterator it = rsdtypemap_.find( s1 );

			if( it == rsdtypemap_.end() ){
				GoapRsdTypeOP rsdtypeinfo = new GoapRsdType;
				//std::cout << "Adding new residue info for " << rsdtype.name() << std::endl;
				rsdtypeinfo->setup_rsdtype( &rsdtype );
				rsdtypemap_[ s1 ] = rsdtypeinfo;
				it = rsdtypemap_.find( s1 );
			}

			GoapRsdTypeOP rsdtypeinfo = it->second;

			// i1: Natoms
		linestream >> i1;
		for( Size iatm = 1; iatm <= i1; ++iatm ){
			if( read_type == 1 ) { // Atom order definition
				std::string s2("");
				linestream >> s2;
				assert( rsdtype.has( s2 ) );

				Size const atmno( rsdtype.atom_index( s2 ) );

				rsdtypeinfo->add_atmname_using( s2 );
				rsdtypeinfo->set_using( atmno, true );

				i_id++;
				// match rsdtype atmno with full Goap index (1~167)
				rsdtypeinfo->set_atmid( atmno, i_id );
				//std::cout << "idmap: " << i_id << " " << s1 << " " << s2 << std::endl;
			} else if( read_type == 2 ){ // Bond definition
				Size i2;
				linestream >> i2;

				// index will start from sidechain, so add 4 to skip backbone atoms
				std::string const atmname = rsdtypeinfo->atmname_using( iatm+4 );
				std::string const rootname = rsdtypeinfo->atmname_using( i2 );

				Size const atmno( rsdtype.atom_index( atmname ) ) ;
				Size const rootno( rsdtype.atom_index( rootname ) );

				// Note that root_atom_[1:4] are assigned as 0
				rsdtypeinfo->set_root_atom( atmno, rootno );
				rsdtypeinfo->set_branch_atom( rootno, atmno );
			} else if( read_type == 3 ){ // Angle definition
				Size i2;
				linestream >> i2;

				// index will start from sidechain, so add 4 to skip backbone atoms
				std::string const atmname = rsdtypeinfo->atmname_using( iatm+4 );
				std::string const anglename = rsdtypeinfo->atmname_using( i2 );

				Size const atmno( rsdtype.atom_index( atmname ) ) ;
				Size const angleno( rsdtype.atom_index( anglename ) );

				rsdtypeinfo->set_angle_atom( atmno, angleno );
			}
		}
	} //while instream

	// Setup connectivity as angles definitions are read
	GoapRsdTypeMap::const_iterator it;

	for( it = rsdtypemap_.begin(); it != rsdtypemap_.end(); ++it ){
		GoapRsdTypeOP rsdtypeinfo = it->second;
		chemical::ResidueType const &rsdtype = rsdtypeset->name_map( it->first );

		rsdtypeinfo->setup_connectivity( rsdtype );
	}
}

void
GoapEnergy::read_potential_values( std::string const distance_file,
				  std::string const angle_file )
{

	TR << "Reading energy table." << std::endl;

  utility::io::izstream instream;
  std::string line;

  distance_table_.dimension( N_ATMTYPES, N_ATMTYPES, N_DISTANCE_BINS );
  angle_table_.dimension( N_ATMTYPES, N_ATMTYPES, N_DISTANCE_BINS, N_ANGLE2_BINS, N_ANGLE1_BINS );

  // Convention: 
  // m: distance bin   mm, lx: angle bin
  // pp: score value   atm1, atm2: atom type index

  // Read database for distance part
  Size read_type( 0 );
  basic::database::open( instream, distance_file );

  while( instream ){
    getline( instream, line );
    std::istringstream linestream( line );
		std::string s1, s2;
		Size m, atm1, atm2, i1;
		Real pp, dis;

    linestream >> s1;

    if( s1.compare( "#MAP" ) == 0 ){
      read_type = 1;
      continue;
    } else if( s1.compare( "#DATA" ) == 0 ){
      read_type = 2;
      continue;
		} else if( s1.compare( "#END" ) == 0 ){
			break;
    }

    if( read_type == 1 ){
      linestream >> i1 >> dis; // dis = min distance
      Size const idis = (Size)( dis*2 );
      distbin_map_[ idis ] = i1; // Store which bin it is for Size(distance*2)

    } else if( read_type == 2 ){
      // originaly each of atm1, atm2 are two integers (like i,j for atm1 / k,l for atm2)
      linestream >> s2 >> m >> pp >> atm1 >> atm2;
      distance_table_(atm1,atm2,m) = pp;
      distance_table_(atm2,atm1,m) = pp;
    }
  }

  // Read database for angle part
  basic::database::open( instream, angle_file );
  read_type = 0;

  while( instream ){
    getline( instream, line );
    std::istringstream linestream( line );
		std::string s1;
    linestream >> s1;

    if( s1.compare( "#DATA" ) == 0){
      read_type = 1;
      continue;
    }
    
    if( read_type == 1 ){
			std::string s2, s3;
			Size m, atm1, atm2;
      linestream >> s2 >> m >> s3 >> atm1 >> atm2;

      // Angle1 
      for( Size mm = 1; mm <= N_ANGLE1_BINS; ++mm ){ //5
				getline( instream, line );
				std::istringstream linestream2( line );
				Size mx, m2;
				linestream2 >> m2 >> mx;

				// Angle2
				for ( Size lx = 1; lx <= N_ANGLE2_BINS; ++lx ){ //12
					Real pp2;
					linestream2 >> pp2;

					angle_table_(atm1,atm2,m,lx,mx) = (int)(pp2*1e4);
				}
      }

      for( Size lx = 1; lx <= N_ANGLE2_BINS; ++lx ){
				angle_table_(atm2,atm1,m,lx,1) = angle_table_(atm1,atm2,m,lx,3);
				angle_table_(atm2,atm1,m,lx,2) = angle_table_(atm1,atm2,m,lx,4);
				angle_table_(atm2,atm1,m,lx,3) = angle_table_(atm1,atm2,m,lx,1);
				angle_table_(atm2,atm1,m,lx,4) = angle_table_(atm1,atm2,m,lx,2);
				angle_table_(atm2,atm1,m,lx,5) = angle_table_(atm1,atm2,m,lx,5);
      }
    } // if read_type == 1
  } // while instream

} // read_database

//////////////////////////////////////////////////////
///////////// Energy evaluation functions
void
GoapEnergy::residue_pair_energy( conformation::Residue const &rsd1,
				conformation::Residue const &rsd2,
 			  pose::Pose const & , //pose,
				ScoreFunction const &,
				EnergyMap & emap
				) const
{

	//std::cout << "Respair start: " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

  // Passing terminus is original Goap convention
  if( !eval_res( rsd1.seqpos() ) || !eval_res( rsd2.seqpos() ) ) return;

	// Evaluate angle term only if sequence separation is long enough
	bool const evaluate_angle_term( std::abs(int(rsd1.seqpos() - rsd2.seqpos())) >= int(MIN_SEQ_SEPARATIONS) );

	// pass if sequence separation is under MIN_SEQ_SEPARATIONS
  GoapRsdTypeOP rsd1type = rsdtypemap_[ rsd1.name3() ];
  GoapRsdTypeOP rsd2type = rsdtypemap_[ rsd2.name3() ];

	Real const MAX_D2_CUT( max_dis_*max_dis_ );
	Size const res1( rsd1.seqpos() );
	Size const res2( rsd2.seqpos() );

  Real Edist_sum( 0.0 );
  Real Eang_sum( 0.0 );

  for( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ){
    if( atm1 > rsd1type->natom() ) continue;
		if( !rsd1type->is_using( atm1 ) ) continue;

		Size const atype1( rsd1type->atmid(atm1) );
		// Just make sure
    if( atype1 == 0 ) continue;

    Vector const &xn1 = xn( res1, atm1 );
    Vector const &xd1 = xd( res1, atm1 );
    Vector const &xyz1 = rsd1.xyz( atm1 );

    for( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ){
      if( atm2 > rsd2type->natom() ) continue;
			if( !rsd2type->is_using( atm2 ) ) continue;

			Size const atype2( rsd2type->atmid(atm2) );
			if( atype2 == 0 ) continue;

      Vector const &xn2 = xn( res2, atm2 );
      Vector const &xd2 = xd( res2, atm2 );
      Vector const &xyz2 = rsd2.xyz( atm2 );

			Vector const dxyz = xyz2 - xyz1;

      Real const d2 = dxyz.dot( dxyz );

      if ( !(d2 <= MAX_D2_CUT ) )continue;

      Real const dis( std::sqrt( d2 ) );

      Real const Edist = get_distance_score( dis, atype1, atype2 );
			Edist_sum += Edist;

			/*
			printf("atm1/atm2/type1/type2/d/Edist/Eang: %3d %3d %3d %3d", atm1, atm2, atype1, atype2);
			std::cout << std::setw(5) << rsd1.atom_name( atm1 ) << std::setw(5) << rsd2.atom_name( atm2 );
			printf(" %8.3f %8.3f", dis, Edist );
			std::cout << std::endl;
			*/

			if( evaluate_angle_term ){
				Real const Eang = get_angle_score( dis, atype1, atype2, xn1, xd1, xn2, xd2, xyz1, xyz2 );
				Eang_sum += Eang;
				/*
				printf("atm1/atm2/type1/type2/d/Edist/Eang: %3d %3d %3d %3d", atm1, atm2, atype1, atype2);
				std::cout << std::setw(5) << rsd1.atom_name( atm1 ) << std::setw(5) << rsd2.atom_name( atm2 );
				printf(" %8.3f %8.3f %8.3f\n", dis, Edist, Eang );
				*/
			}

    }
  }

	//std::cout << "res pair scoring b/w ";
	//printf("%4d %4d %8.3f %8.3f\n", rsd1.seqpos(), rsd2.seqpos(), Edist_sum, Eang_sum);

  emap[ scoring::goap       ] += Edist_sum+Eang_sum;
  emap[ scoring::goap_dist  ] += Edist_sum;
  emap[ scoring::goap_angle ] += Eang_sum;

}

Real
GoapEnergy::get_distance_score( Real const dist,
			       Size const atype1,
			       Size const atype2 ) const
{
	Size const kbin = (Size)(dist*2);
  Size const binno = distbin_map( kbin );
	if( binno == 0 ) return 0.0;

  if( continuous() ){
    // cubic spline
    return 0.0;

  } else {
		//std::cout << "DIST: " << atype1 << " " << atype2 << " " << dist;
		//std::cout << " " << binno << " " << distance_table_( atype1, atype2, binno ) << std::endl;
    return distance_table_( atype1, atype2, binno );
  }
}
  

Real
GoapEnergy::get_angle_score( Real const dist,
			     Size const atype1,
			     Size const atype2,
			     Vector const xn1,
			     Vector const xd1,
			     Vector const xn2,
			     Vector const xd2,
			     Vector const xyz1,
			     Vector const xyz2
			     ) const
{
  Real score( 0.0 );

	Size const kbin = (Size)(dist*2);
  Size const binno = distbin_map( kbin );
	if( binno == 0 ) return 0.0;

	Vector const dxyz( xyz2 - xyz1 );
  Vector const vt = dxyz/dist; // 1->2
  Vector const vt2 = -vt;

  // use cosine(ang), which is faster than angle and no worry about periodicity :)
  Real const cs1 = calc_cosineang( xn1, vt );
  Size const mm1 = std::min( int((cs1+1.001)*6.0)+1, 12);

  Real const phi1 = calc_phi( xn1, xd1, vt );
  Size const mm2 = std::min( int((phi1+180.001)/30.0)+1, 12 );

  Real const cs2 = calc_cosineang( xn2, vt2 );
  Size const mm3 = std::min( int((cs2+1.001)*6.0)+1, 12 );

  Real const phi2 = calc_phi( xn2, xd2, vt2 );
	Size const mm4 = std::min( int((phi2+180.001)/30.0)+1, 12 );

  Vector const vec_proj1 = xyz1 + xn2;
  Vector const vec_proj2 = xyz2 + xn1;

	// Convention is inverse
  Real const dih = -numeric::dihedral_degrees( vec_proj1, xyz1, xyz2, vec_proj2 );
  Size const mm5 = std::min( int((dih+180.001)/30.0) + 1, 12);

	/*
	printf("xn1,xn2: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
				 xn1[0],xn1[1],xn1[2],
				 xn2[0],xn2[1],xn2[2]);
	printf("xd1,xd2: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
				 xd1[0],xd1[1],xd1[2],
				 xd2[0],xd2[1],xd2[2]);
	printf("vt: %8.3f %8.3f %8.3f\n",
				 vt[0],vt[1],vt[2]);
	printf("angle: %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
				 phi1,phi2,cs1,cs2,dih);
	printf("mm: %3d %3d %3d %3d %3d\n",mm1,mm2,mm3,mm4,mm5);
	*/

  if( continuous() ){
    score += 0.0;

  } else {
    score += angle_table_(atype1, atype2, binno, mm1, 1);
    score += angle_table_(atype1, atype2, binno, mm2, 2);
    score += angle_table_(atype1, atype2, binno, mm3, 3);
    score += angle_table_(atype1, atype2, binno, mm4, 4);
    score += angle_table_(atype1, atype2, binno, mm5, 5);
  }

  return score*1.0e-4; // revert by scaling factor
}
bool
GoapEnergy::calculate_dipoles( pose::Pose const &pose, 
			      conformation::Residue const &rsd1,
			      GoapRsdTypeCOP rsdtype,
			      Size const atm1,
			      Vector &xn, 
			      Vector &xd
			      ) const
{

  Size const i2( rsdtype->i2(atm1) );
  Size const i3( rsdtype->i3(atm1) );
  bool connected_by_twobonds = rsdtype->connected_by_twobonds(atm1);

  Vector xyz1 = rsd1.xyz( atm1 );
  Vector xyz3 = rsd1.xyz( i3 );

  // careful; xyz2 can be from prv rsd
  Vector xyz2;
  if( i2 == 999 ){
    Size const i_prv( rsd1.seqpos() - 1 );

    if( i_prv == 0 ){
			return false;

		} else {
			conformation::Residue const &rsd_prv( pose.residue( i_prv ) );
			assert( rsd_prv.has(" C  ") );
			Size const iatm = rsd_prv.atom_index(" C  ");
			xyz2 = rsd_prv.xyz( iatm );
		}

  } else {
    xyz2 = rsd1.xyz( i2 );
  }

  if( !connected_by_twobonds ){
    xn = xyz2 - xyz1;
  } else {
    xn = xyz2 - xyz1 + xyz3 - xyz1;
  }

  Vector const v2 = xyz3 - xyz1;
  Vector const v1 = xn.cross( v2 );
  xd = v1.cross( xn );

  //normalize
  Real const xx1 = std::sqrt( xd.dot( xd ));
  Real const xx2 = std::sqrt( xn.dot( xn ));

  xd /= xx1;
  xn /= xx2;

	//Make sure xn & xd does not blow up - by making it an arbitrary normalized vector
	// This will invoke dummy outputs for unlikely poses, but will be safe from SegFault
	Vector I( 0.0 ); I[0] = 1.0;
	if( xx1 != xx1 || xx1 < 1.0e-6 ) xd = I;
	if( xx2 != xx2 || xx2 < 1.0e-6 ) xn = I;

	// xd, xn are returned Vectors
	return true;
}

//////////////////////////////////////////////////////
///////////// Derivatives

void
GoapEnergy::setup_for_derivatives( pose::Pose & ,
																	 ScoreFunction const & )  const
{}


void
GoapEnergy::eval_residue_pair_derivatives(
		conformation::Residue const & ,//rsd1,
		conformation::Residue const & ,//rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & , //min_data,
		pose::Pose const &,
		EnergyMap const & ,//weights,
		utility::vector1< DerivVectorPair > & ,//r1_atom_derivs,
		utility::vector1< DerivVectorPair > & //r2_atom_derivs
	     ) const
{}


//////////////////////////////////////////////////////
///////////// Packing - not implemented

void
GoapEnergy::setup_for_packing(
		 pose::Pose & ,//pose,
     utility::vector1< bool > const & , //residues_repacking,
     utility::vector1< bool > const &
     ) const
{}

} //namespace methods
} //namespace scoring
} //namespace core
