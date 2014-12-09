// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

template< typename T>
class NonSymmetric;

template< typename T>
class Symmetric;

template< typename T>
class ToBeDetermined;

template< typename T>
class real_traits;

template< typename T, class MatrixType = NonSymmetric<T>, class Type = real_traits<T> >
class LapackFunctions;

template< typename T, class MatrixType = NonSymmetric<T>, class Type = real_traits<T> >
bool Seigsolv_(ObjexxFCL::FArray2D<T>& tmp,ObjexxFCL::FArray2D<T>& eigvec, ObjexxFCL::FArray1D<T>& eigval,bool ComputeVecs=true);


namespace ObjexxFCL {

 template< typename T>
 class FArray {
	 friend class LapackFunctions<T, Symmetric <T> >;
	 friend class LapackFunctions<T, NonSymmetric<T> >;
	 friend class LapackFunctions<T, ToBeDetermined<T> >;
	 friend bool Seigsolv_(ObjexxFCL::FArray2D<T>& tmp,ObjexxFCL::FArray2D<T>& eigvec, ObjexxFCL::FArray1D<T>& eigval,bool ComputeVecs=true)<T,Symmetric<T> >;
 protected:
	 T* data(){};
 };

template< typename T>
class FArray1D :public FArray<T> {

};

template< typename T>
class FArray2D :public FArray1D<T> {

};
};


template< typename T>
class Unknown {
 public:
	static bool isSymmetric(const ObjexxFCL::FArray2D<T>&) {return false;};
};

template< typename T>
class Symmetric {
 public:
		static bool isSymmetric(const ObjexxFCL::FArray2D<T>&) {return true;};
};

template< typename T>
class ToBeDetermined {
 public:
	static bool isSymmetric(const ObjexxFCL::FArray2D<T> &M) {
		std::cerr <<"testing matrix symmetry..." << std::endl;
		return M.symmetric();
	};
};

template < typename T>
class real_traits {
 public:
	static bool isfloat();
};




template< typename T, class Symmetrie , class Type >
	class LapackFunctions {
 public:
	/* get eigenvectors and eigenvalues of M --- input Array stays intact*/
	static    bool eigsolv(const ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray2D<T>& eigvec,ObjexxFCL::FArray1D<T>& eigval);

	/* get eigenvectors and eigenvalues of 1st Arg --- works in input Array -- destroys data*/
	static    bool eigsolv_destroy_input(ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray2D<T>& eigvec,ObjexxFCL::FArray1D<T>& eigval);

	/* compute eigenvalues only -- input Array stays intact */
	static    bool eigval(const ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray1D<T>& eigval);

	/* compute eigenvalues --- invalidates Input Array during call */
	static    bool eigval_destroy_input(ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray1D<T>& eigval);

	/* compute determinant */
	//    static    T det(const ObjexxFCL::FArray2D<T>&);
	//    static    T det(ObjexxFCL::FArray2D<T>&);
 protected:
	static bool eigsolv_(ObjexxFCL::FArray2D<T>& tmp,ObjexxFCL::FArray2D<T>& eigvec, ObjexxFCL::FArray1D<T>& eigval,bool ComputeVecs=true);
};


template<typename T , class MatrixType, class Type>
	bool LapackFunctions<T,MatrixType,Type>::eigsolv(
							 const ObjexxFCL::FArray2D<T>& M,
							 ObjexxFCL::FArray2D<T>&    eigvec,
							 ObjexxFCL::FArray1D<T>& eigval) {
	ObjexxFCL::FArray2D<T> tmp(M);
	return eigsolv_(tmp,eigvec,eigval);
}

template<typename T , class MatrixType, class Type>
	bool LapackFunctions<T,MatrixType,Type>::eigval(const ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray1D<T>& eigval) {
	ObjexxFCL::FArray2D<T> tmp(M);
	ObjexxFCL::FArray2D<T> eigvec;
	return eigsolv_(tmp,eigvec,eigval,false);
}

template<typename T , class MatrixType, class Type>
	bool LapackFunctions<T,MatrixType,Type>::eigval_destroy_input(ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray1D<T>& eigval) {
	ObjexxFCL::FArray2D<T> eigvec;
	return eigsolv_(M,eigvec,eigval,false);
}
template<typename T , class MatrixType, class Type>
	bool LapackFunctions<T,MatrixType,Type>::eigsolv_destroy_input(ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray2D<T>&    eigvec,ObjexxFCL::FArray1D<T>& eigval) {
	return eigsolv_(M,eigvec,eigval);
}


 template<>
bool real_traits<float>::isfloat() {return true;}

template<>
bool real_traits<double>::isfloat() {return false;}

template<typename T , class MatrixType, class Type>
	bool LapackFunctions<T, MatrixType ,Type>::eigsolv_(ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray2D<T>& eigvec,ObjexxFCL::FArray1D<T>& eigval,bool bComputeVectors) {
	T* pM=M.data();
	return false;};


template<typename T , class MatrixType, class Type>
	bool Seigsolv_(ObjexxFCL::FArray2D<T>& M,ObjexxFCL::FArray2D<T>& eigvec,ObjexxFCL::FArray1D<T>& eigval,bool bComputeVectors) {
	T* pM=M.data();
	return false;};
