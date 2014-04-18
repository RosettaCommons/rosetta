#include <core/chemical/Elements.hh>
#include <ostream>
#include <map>
#include <string>
#include <utility/vector1.hh>
#include <utility/exit.hh>
namespace core {
namespace chemical {
namespace element {


utility::vector1< std::string > & element2name() {
	// static initialization only happens once
	static utility::vector1< std::string > * element2name_ = new utility::vector1< std::string >( setup_element2name() );
	return *element2name_;
}

std::string
name_from_elements(Elements element){
	if( element > total_number_elements ) return "ElementOutOfRange";
	return element2name()[ element ];
}

Elements
elements_from_name(std::string name) {
	std::map< std::string, Elements >::const_iterator iter = name2element().find( name );
	if ( iter == name2element().end() ) {
		utility_exit_with_message( "unrecognized element  " + name );
	}
	return iter->second;
}


inline
std::map< std::string, Elements > & name2element() {
	// static initialization only happens once
	static std::map< std::string, Elements > * name2element_ = new std::map< std::string, Elements >( setup_name2element() );
	return *name2element_;
}

/// @brief setup the map that converts string name to AA enum
std::map< std::string, Elements > setup_name2element() {

	std::map< std::string, Elements > n2element;

    // L-amino acid types
	n2element["H" ] = H;
	n2element["He" ] = He;
	n2element["Li" ] = Li;
	n2element["Be" ] = Be;
	n2element["B" ] = B;
	n2element["C" ] = C;
	n2element["N" ] = N;
	n2element["O" ] = O;
	n2element["F" ] = F;
	n2element["Ne" ] = Ne;
	n2element["Na" ] = Na;
	n2element["Mg" ] = Mg;
	n2element["Al" ] = Al;
	n2element["Si" ] = Si;
	n2element["P" ] = P;
	n2element["S" ] = S;
	n2element["Cl" ] = Cl;
	n2element["Ar" ] = Ar;
	n2element["K" ] = K;
	n2element["Ca" ] = Ca;
	n2element["Sc" ] = Sc;
	n2element["Ti" ] = Ti;
	n2element["V" ] = V;
	n2element["Cr" ] = Cr;
	n2element["Mn" ] = Mn;
	n2element["Fe" ] = Fe;
	n2element["Co" ] = Co;
	n2element["Ni" ] = Ni;
	n2element["Cu" ] = Cu;
	n2element["Zn" ] = Zn;
	n2element["Ga" ] = Ga;
	n2element["Ge" ] = Ge;
	n2element["As" ] = As;
	n2element["Se" ] = Se;
	n2element["Br" ] = Br;
	n2element["Kr" ] = Kr;
	n2element["Rb" ] = Rb;
	n2element["Sr" ] = Sr;
	n2element["Y" ] = Y;
	n2element["Zr" ] = Zr;
	n2element["Nb" ] = Nb;
	n2element["Mo" ] = Mo;
	n2element["Tc" ] = Tc;
	n2element["Ru" ] = Ru;
	n2element["Rh" ] = Rh;
	n2element["Pd" ] = Pd;
	n2element["Ag" ] = Ag;
	n2element["Cd" ] = Cd;
	n2element["In" ] = In;
	n2element["Sn" ] = Sn;
	n2element["Sb" ] = Sb;
	n2element["Te" ] = Te;
	n2element["I" ] = I;
	n2element["Xe" ] = Xe;
	n2element["Cs" ] = Cs;
	n2element["Ba" ] = Ba;
	n2element["La" ] = La;
	n2element["Ce" ] = Ce;
	n2element["Pr" ] = Pr;
	n2element["Nd" ] = Nd;
	n2element["Pm" ] = Pm;
	n2element["Sm" ] = Sm;
	n2element["Eu" ] = Eu;
	n2element["Gd" ] = Gd;
	n2element["Tb" ] = Tb;
	n2element["Dy" ] = Dy;
	n2element["Ho" ] = Ho;
	n2element["Er" ] = Er;
	n2element["Tm" ] = Tm;
	n2element["Yb" ] = Yb;
	n2element["Lu" ] = Lu;
	n2element["Hf" ] = Hf;
	n2element["Ta" ] = Ta;
	n2element["W" ] = W;
	n2element["Re" ] = Re;
	n2element["Os" ] = Os;
	n2element["Ir" ] = Ir;
	n2element["Pt" ] = Pt;
	n2element["Au" ] = Au;
	n2element["Hg" ] = Hg;
	n2element["Tl" ] = Tl;
	n2element["Pb" ] = Pb;
	n2element["Bi" ] = Bi;
	n2element["Po" ] = Po;
	n2element["At" ] = At;
	n2element["Rn" ] = Rn;
	n2element["Fr" ] = Fr;
	n2element["Ra" ] = Ra;
	n2element["Ac" ] = Ac;
	n2element["Th" ] = Th;
	n2element["Pa" ] = Pa;
	n2element["U" ] = U;
	n2element["Np" ] = Np;
	n2element["Pu" ] = Pu;
	n2element["Am" ] = Am;
	n2element["Cm" ] = Cm;
	n2element["Bk" ] = Bk;
	n2element["Cf" ] = Cf;
	n2element["Es" ] = Es;
	n2element["Fm" ] = Fm;
	n2element["Md" ] = Md;
	n2element["No" ] = No;
	n2element["Lr" ] = Lr;
	n2element["Rf" ] = Rf;
	n2element["Db" ] = Db;
	n2element["Sg" ] = Sg;
	n2element["Bh" ] = Bh;
	n2element["Hs" ] = Hs;
	n2element["Mt" ] = Mt;
	n2element["Ds" ] = Ds;
	n2element["Rg" ] = Rg;
	n2element["Cn" ] = Cn;
	n2element["Uut" ] = Uut;
	n2element["Fl" ] = Fl;
	n2element["Uup" ] = Uup;
	n2element["Lv" ] = Lv;
	n2element["Uus" ] = Uus;
	n2element["Uuo" ] = Uuo;

	return n2element;
}


/// @brief setup the vector that maps AA enum to string name
utility::vector1< std::string > setup_element2name() {

	utility::vector1< std::string > element2n( total_number_elements );

	for ( std::map< std::string, Elements >::const_iterator iter = name2element().begin(),
		iter_end = name2element().end(); iter != iter_end; ++iter ) {
		element2n[ iter->second ] = iter->first;
	}

	return element2n;
}

}
}
}


