// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Frank DiMaio


//  implementation loosely based on Kiss-FFT's fast fourier transform
//     for licensing info see external/kiss_fft_v1_2_8/COPYING

#include <numeric/fourier/kiss_fft.hh>
#include <cstdlib> //g++ 4.3.2 requires for exit()
//#include <string.h> //g++ 4.3.2 requires for memcpy()
// tracer
#include <iostream>


namespace numeric {
namespace fourier {


// replace the global variables with classes that protect access to buffer data
kiss_fft_cpx* get_scratch_buff( size_t nbuf ) {
	static kiss_fft_cpx *scratchbuf=nullptr;
	static size_t nscratchbuf=0;
	if ( nscratchbuf < nbuf ) {
		delete [] scratchbuf;
		scratchbuf = new kiss_fft_cpx[ nbuf ];
		nscratchbuf = nbuf;
	}
	// free memory
	// amw: do not check <= 0 as it is a size_t, eq is sufficient
	if ( nbuf == 0 ) {
		delete [] scratchbuf;
		scratchbuf = nullptr;
		nscratchbuf = 0;
	}
	return (scratchbuf);
}
kiss_fft_cpx* get_tmp_buff( size_t nbuf ) {
	static kiss_fft_cpx *tmpbuf=nullptr;
	static size_t ntmpbuf=0;
	if ( ntmpbuf < nbuf ) {
		delete [] tmpbuf;
		tmpbuf = new kiss_fft_cpx[ nbuf ];
		ntmpbuf = nbuf;
	}
	// free memory
	if ( nbuf == 0 ) {
		delete [] tmpbuf;
		tmpbuf = nullptr;
		ntmpbuf = 0;
	}
	return (tmpbuf);
}

///////////////////////////////////////
// 1d c->c fft
///////////////////////////////////////
void kf_bfly2(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	int m ) {
	kiss_fft_cpx * Fout2;
	kiss_fft_cpx * tw1 = st->twiddles();
	kiss_fft_cpx t;
	Fout2 = Fout + m;
	do{
		t = (*Fout2) * (*tw1);
		tw1 += fstride;
		*Fout2 = (*Fout) - (t);
		*Fout +=  t;
		++Fout2;
		++Fout;
	}while (--m);
}

void kf_bfly4(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	const size_t m ) {
	kiss_fft_cpx *tw1,*tw2,*tw3;
	kiss_fft_cpx scratch[6];
	size_t k=m;
	const size_t m2=2*m;
	const size_t m3=3*m;

	tw3 = tw2 = tw1 = st->twiddles();

	do {
		scratch[0] = Fout[m] * (*tw1);
		scratch[1] = Fout[m2] * (*tw2);
		scratch[2] = Fout[m3] * (*tw3);

		scratch[5] = (*Fout) - scratch[1];
		*Fout += scratch[1];
		scratch[3] = scratch[0] + scratch[2];
		scratch[4] = scratch[0] - scratch[2];
		Fout[m2] = (*Fout) - scratch[3];
		tw1 += fstride;
		tw2 += fstride*2;
		tw3 += fstride*3;
		*Fout += scratch[3];

		if ( st->inverse() ) {
			Fout[m]  = kiss_fft_cpx( scratch[5].real() - scratch[4].imag() , scratch[5].imag() + scratch[4].real());
			Fout[m3] = kiss_fft_cpx( scratch[5].real() + scratch[4].imag() , scratch[5].imag() - scratch[4].real());
		} else {
			Fout[m]  = kiss_fft_cpx( scratch[5].real() + scratch[4].imag() , scratch[5].imag() - scratch[4].real());
			Fout[m3] = kiss_fft_cpx( scratch[5].real() - scratch[4].imag() , scratch[5].imag() + scratch[4].real());
		}
		++Fout;
	}while(--k);
}

void kf_bfly3(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	size_t m
) {
	size_t k=m;
	const size_t m2 = 2*m;
	kiss_fft_cpx *tw1,*tw2;
	kiss_fft_cpx scratch[5];
	kiss_fft_cpx epi3;
	epi3 = st->twiddles()[fstride*m];

	tw1=tw2=st->twiddles();

	do{
		scratch[1] = Fout[m] * (*tw1);
		scratch[2] = Fout[m2] * (*tw2);

		scratch[3] = scratch[1]+scratch[2];
		scratch[0] = scratch[1]-scratch[2];
		tw1 += fstride;
		tw2 += fstride*2;

		Fout[m] = kiss_fft_cpx( Fout->real() - (scratch[3].real()/2) , Fout->imag() - (scratch[3].imag()/2) );
		scratch[0] *= epi3.imag();
		*Fout += scratch[3];

		Fout[m2] = kiss_fft_cpx( Fout[m].real() + scratch[0].imag() , Fout[m].imag() - scratch[0].real());
		Fout[m] = kiss_fft_cpx(  Fout[m].real() - scratch[0].imag() , Fout[m].imag() + scratch[0].real());

		++Fout;
	}while(--k);
}

void kf_bfly5(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	int m ) {
	kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
	int u;
	kiss_fft_cpx scratch[13];
	kiss_fft_cpx * twiddles = st->twiddles();
	kiss_fft_cpx *tw;
	kiss_fft_cpx ya,yb;
	ya = twiddles[fstride*m];
	yb = twiddles[fstride*2*m];

	Fout0=Fout;
	Fout1=Fout0+m;
	Fout2=Fout0+2*m;
	Fout3=Fout0+3*m;
	Fout4=Fout0+4*m;

	tw=st->twiddles();
	for ( u=0; u<m; ++u ) {
		scratch[0] = *Fout0;

		scratch[1] = *Fout1 * tw[u*fstride];
		scratch[2] = *Fout2 * tw[2*u*fstride];
		scratch[3] = *Fout3 * tw[3*u*fstride];
		scratch[4] = *Fout4 * tw[4*u*fstride];

		scratch[7] = scratch[1]+scratch[4];
		scratch[10] = scratch[1]-scratch[4];
		scratch[8] = scratch[2]+scratch[3];
		scratch[9] = scratch[2]-scratch[3];

		*Fout0 = kiss_fft_cpx( Fout0->real() + scratch[7].real() + scratch[8].real(),
			Fout0->imag() + scratch[7].imag() + scratch[8].imag() );

		scratch[5] = kiss_fft_cpx(  scratch[0].real() + scratch[7].real()*ya.real() + scratch[8].real()*yb.real() ,
			scratch[0].imag() + scratch[7].imag()*ya.real() + scratch[8].imag()*yb.real() );

		scratch[6] = kiss_fft_cpx(   (scratch[10].imag()*ya.imag()) + (scratch[9].imag()*yb.imag()) ,
			-(scratch[10].real()*ya.imag()) - (scratch[9].real()*yb.imag()));

		*Fout1 = scratch[5]-scratch[6];
		*Fout4 = scratch[5]+scratch[6];

		scratch[11] = kiss_fft_cpx( scratch[0].real() + (scratch[7].real()*yb.real()) + (scratch[8].real()*ya.real()),
			scratch[0].imag() + (scratch[7].imag()*yb.real()) + (scratch[8].imag()*ya.real()));
		scratch[12] = kiss_fft_cpx( -(scratch[10].imag()*yb.imag()) + (scratch[9].imag()*ya.imag()),
			(scratch[10].real()*yb.imag()) - (scratch[9].real()*ya.imag()));

		*Fout2 = scratch[11]+scratch[12];
		*Fout3 = scratch[11]-scratch[12];

		++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
	}
}

// perform the butterfly for one stage of a mixed radix FFT
void kf_bfly_generic(
	kiss_fft_cpx * Fout,
	const size_t fstride,
	const kiss_fft_cfg st,
	int m,
	int p
) {
	int u, q1, q;
	kiss_fft_cpx * twiddles = st->twiddles();
	kiss_fft_cpx t;
	int Norig = st->nfft();

	//CHECKBUF(scratchbuf,nscratchbuf,p);
	kiss_fft_cpx *scratchbuf = get_scratch_buff(p);

	for ( u = 0; u < m; ++u ) {
		int k = u;
		for ( q1 = 0 ; q1 < p; ++q1 ) {
			scratchbuf[ q1 ] = Fout[ k ];
			k += m;
		}

		k = u;
		for ( q1 = 0 ; q1 < p ; ++q1 ) {
			int twidx = 0;
			Fout[ k ] = scratchbuf[ 0 ];
			for ( q = 1; q < p; ++q ) {
				twidx += fstride * k;
				if ( twidx >= Norig ) twidx -= Norig;
				t = scratchbuf[ q ] * twiddles[ twidx ];
				Fout[ k ] += t;
			}
			k += m;
		}
	}
}

void kf_work(
	kiss_fft_cpx * Fout,
	const kiss_fft_cpx * f,
	const size_t fstride,
	int in_stride,
	int * factors,
	const kiss_fft_cfg st
) {
	kiss_fft_cpx * Fout_beg=Fout;
	const int p = *factors++; /* the radix  */
	const int m = *factors++; /* stage's fft length/p */
	const kiss_fft_cpx * Fout_end = Fout + p*m;

	if ( m==1 ) {
		do{
			*Fout = *f;
			f += fstride*in_stride;
		}while(++Fout != Fout_end );
	} else {
		do{
			// recursive call:
			// DFT of size m*p performed by doing
			// p instances of smaller DFTs of size m,
			// each one takes a decimated version of the input
			kf_work( Fout , f, fstride*p, in_stride, factors,st);
			f += fstride*in_stride;
		}while( (Fout += m) != Fout_end );
	}

	Fout=Fout_beg;

	// recombine the p smaller DFTs
	switch (p) {
	case 2 : kf_bfly2(Fout,fstride,st,m); break;
	case 3 : kf_bfly3(Fout,fstride,st,m); break;
	case 4 : kf_bfly4(Fout,fstride,st,m); break;
	case 5 : kf_bfly5(Fout,fstride,st,m); break;
	default : kf_bfly_generic(Fout,fstride,st,m,p); break;
	}
}

void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride) {
	if ( fin == fout ) {
		//CHECKBUF(tmpbuf,ntmpbuf,st->nfft);
		kiss_fft_cpx *tmpbuf = get_tmp_buff(st->nfft());
		kf_work(tmpbuf,fin,1,in_stride, st->factors(),st);
		//memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft());
		for ( int i=0; i<st->nfft(); ++i ) fout[i]=tmpbuf[i];
	} else {
		kf_work( fout, fin, 1,in_stride, st->factors(),st );
	}
}

void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout) {
	kiss_fft_stride(cfg,fin,fout,1);
}

void kiss_fft_split(
	kiss_fftsplit_cfg st,
	const kiss_fft_scalar *rin,
	const kiss_fft_scalar *iin,
	kiss_fft_scalar *rout,
	kiss_fft_scalar *iout,
	int fin_stride,
	int fout_stride) {
	// input buffer timedata is stored row-wise
	int k, ncfft;
	kiss_fft_cpx *tmpbuf = get_tmp_buff(st->nfft());

	ncfft = st->substate()->nfft();
	for ( k = 0; k < ncfft; ++k ) {
		st->tmpbuf()[k] = kiss_fft_cpx( rin[k*fin_stride],iin[k*fin_stride] );
	}
	kiss_fft (st->substate(), st->tmpbuf(), tmpbuf);
	for ( k = 0; k < ncfft; ++k ) {
		rout[k*fout_stride] = tmpbuf[k].real();
		iout[k*fout_stride] = tmpbuf[k].imag();
	}
}


// not really necessary to call, but if someone is doing in-place ffts, they may want to free the
//   buffers from CHECKBUF
void kiss_fft_cleanup(void) {
	get_scratch_buff(0);
	get_tmp_buff(0);
}

int kiss_fft_next_fast_size(int n) {
	while ( 1 ) {
		int m=n;
		while ( (m%2) == 0 ) m/=2;
		while ( (m%3) == 0 ) m/=3;
		while ( (m%5) == 0 ) m/=5;
		if ( m<=1 ) {
			break; // n is completely factorable by twos, threes, and fives
		}
		n++;
	}
	return n;
}


///////////////////////////////////////
/// real fft
///////////////////////////////////////
void kiss_fftr(kiss_fftr_cfg st, const kiss_fft_scalar *timedata, kiss_fft_cpx *freqdata) {
	// input buffer timedata is stored row-wise
	int k,ncfft;
	kiss_fft_cpx fpnk,fpk,f1k,f2k,tw,tdc;

	if ( st->substate()->inverse() ) {
		std::cerr << "kiss fft usage error: improper alloc\n";
		exit(1);
	}

	ncfft = st->substate()->nfft();

	// perform the parallel fft of two real signals packed in real,imag
	kiss_fft( st->substate() , (const kiss_fft_cpx*)timedata, st->tmpbuf() );
	// The real part of the DC element of the frequency spectrum in st->tmpbuf
	// contains the sum of the even-numbered elements of the input time sequence
	// The imag part is the sum of the odd-numbered elements
	//
	// The sum of tdc.r and tdc.i is the sum of the input time sequence.
	//   yielding DC of input time sequence
	// The difference of tdc.r - tdc.i is the sum of the input (dot product) [1,-1,1,-1...
	//   yielding Nyquist bin of input time sequence
	tdc = kiss_fft_cpx( st->tmpbuf()[0].real() , st->tmpbuf()[0].imag() );
	freqdata[0] = kiss_fft_cpx( tdc.real() + tdc.imag() , 0 );
	freqdata[ncfft] = kiss_fft_cpx( tdc.real() - tdc.imag() , 0 );

	//std::cout << ((kiss_fft_cpx*)timedata)[0] << "  " << ((kiss_fft_cpx*)timedata)[1] << std::endl;
	//std::cout << st->tmpbuf()[0] << "  " << st->tmpbuf()[1] << std::endl;

	for ( k=1; k <= ncfft/2 ; ++k ) {
		fpk = st->tmpbuf()[k];
		fpnk = kiss_fft_cpx( st->tmpbuf()[ncfft-k].real(), -st->tmpbuf()[ncfft-k].imag());

		f1k = fpk + fpnk;
		f2k = fpk - fpnk;
		tw  = f2k * st->super_twiddles()[k-1];

		freqdata[k] = kiss_fft_cpx( (f1k.real() + tw.real())/2.0 , (f1k.imag() + tw.imag())/2.0 );
		freqdata[ncfft-k] = kiss_fft_cpx( (f1k.real() - tw.real())/2.0 , (tw.imag() - f1k.imag())/2.0);
	}
}

void kiss_fftri(kiss_fftr_cfg st, const kiss_fft_cpx *freqdata, kiss_fft_scalar *timedata) {
	// input buffer timedata is stored row-wise
	int k, ncfft;

	if ( st->substate()->inverse() == 0 ) {
		std::cerr << "kiss fft usage error: improper alloc\n";
		exit (1);
	}

	ncfft = st->substate()->nfft();

	st->tmpbuf()[0] = kiss_fft_cpx( freqdata[0].real() + freqdata[ncfft].real() ,
		freqdata[0].real() - freqdata[ncfft].real() );

	for ( k = 1; k <= ncfft / 2; ++k ) {
		kiss_fft_cpx fk, fnkc, fek, fok, tmp;
		fk = freqdata[k];
		fnkc = kiss_fft_cpx( freqdata[ncfft - k].real() , -freqdata[ncfft - k].imag() );

		fek = fk + fnkc;
		tmp = fk - fnkc;
		fok = tmp * st->super_twiddles()[k-1];
		st->tmpbuf()[k] =  fek + fok;
		st->tmpbuf()[ncfft - k] = fek - fok;

		st->tmpbuf()[ncfft - k] = kiss_fft_cpx(  st->tmpbuf()[ncfft - k].real() ,
			-st->tmpbuf()[ncfft - k].imag() );
	}
	kiss_fft (st->substate(), st->tmpbuf(), (kiss_fft_cpx *) timedata);
}

///////////////////////////////////////
/// dct-II
///////////////////////////////////////
void kiss_dct(kiss_dct_cfg st, const kiss_fft_scalar *timedata, kiss_fft_scalar *freqdata) {
	// input buffer timedata is stored row-wise
	int k,ncfft;
	kiss_fft_cpx f1k,f2k;

	if ( st->substate()->inverse() ) {
		std::cerr << "kiss fft usage error: improper alloc\n";
		exit(1);
	}

	ncfft = st->substate()->nfft();

	// pack data
	//   ind = [(1:2:n) (n:-2:2)];
	//   a = a(ind);
	kiss_fft_cpx *tmpbuf = get_tmp_buff(ncfft);
	for ( k = 0; k < ncfft/2; ++k ) {
		tmpbuf[k] = kiss_fft_cpx(timedata[2*k],0.0);
	}
	for ( k = ncfft/2; k<ncfft; ++k ) {
		tmpbuf[k] = kiss_fft_cpx(timedata[2*(ncfft-k)-1],0.0);
	}

	// fft
	kiss_fft( st->substate() , tmpbuf, st->tmpbuf() );

	// twiddle
	for ( k=0; k<ncfft ; ++k ) {
		f1k = st->super_twiddles()[k];
		f2k = st->tmpbuf()[k];
		freqdata[k] = 2*(f1k.real()*f2k.real() - f1k.imag()*f2k.imag());
	}
}

///////////////////////////////////////
/// idct-II (dct-iii)
///////////////////////////////////////
void kiss_idct(kiss_dct_cfg st, const kiss_fft_scalar *freqdata, kiss_fft_scalar *timedata) {
	// input buffer timedata is stored row-wise
	int k,ncfft;
	kiss_fft_cpx f1k,f2k;
	kiss_fft_scalar x0 = freqdata[0];

	if ( !st->substate()->inverse() ) {
		std::cerr << "kiss fft usage error: improper alloc\n";
		exit(1);
	}

	ncfft = st->substate()->nfft();

	// twiddle
	kiss_fft_cpx *tmpbuf = get_tmp_buff(ncfft);
	for ( k=0; k<ncfft ; ++k ) {
		f1k = st->super_twiddles()[k];
		f2k = kiss_fft_cpx(freqdata[k],0.0);
		tmpbuf[k] = f1k*f2k;
	}

	// fft
	kiss_fft( st->substate() , tmpbuf, st->tmpbuf() );

	// unpack data
	//  tmp(1:2:n)=(1:n/2);
	//  tmp(2:2:n)=(n:-1:n/2+1);
	//  ind=tmp
	for ( k = 0; k < ncfft/2; ++k ) {
		timedata[2*k] = 2*st->tmpbuf()[k].real() - x0;
	}
	for ( k = ncfft/2; k<ncfft; ++k ) {
		timedata[2*(ncfft-k)-1] = 2*st->tmpbuf()[k].real() - x0;
	}
}


///////////////////////////////////////
/// multidim fft
///////////////////////////////////////
//  This works by tackling one dimension at a time.
//
//  In effect,
//  Each stage starts out by reshaping the matrix into a DixSi 2d matrix.
//  A Di-sized fft is taken of each column, transposing the matrix as it goes.
//
// Here's a 3-d example:
// Take a 2x3x4 matrix, laid out in memory as a contiguous buffer
//  [ [ [ a b c d ] [ e f g h ] [ i j k l ] ]
//    [ [ m n o p ] [ q r s t ] [ u v w x ] ] ]
//
// Stage 0 ( D=2): treat the buffer as a 2x12 matrix
//    [ [a b ... k l]
//      [m n ... w x] ]
//
//    FFT each column with size 2.
//    Transpose the matrix at the same time using kiss_fft_stride.
//
//    [ [ a+m a-m ]
//      [ b+n b-n]
//      ...
//      [ k+w k-w ]
//      [ l+x l-x ] ]
//
//    Note fft([x y]) == [x+y x-y]
//
// Stage 1 ( D=3) treats the buffer (the output of stage D=2) as an 3x8 matrix,
//    [ [ a+m a-m b+n b-n c+o c-o d+p d-p ]
//      [ e+q e-q f+r f-r g+s g-s h+t h-t ]
//      [ i+u i-u j+v j-v k+w k-w l+x l-x ] ]
//
//    And perform FFTs (size=3) on each of the columns as above, transposing
//    the matrix as it goes.  The output of stage 1 is
//        (Legend: ap = [ a+m e+q i+u ]
//                 am = [ a-m e-q i-u ] )
//
//    [ [ sum(ap) fft(ap)[0] fft(ap)[1] ]
//      [ sum(am) fft(am)[0] fft(am)[1] ]
//      [ sum(bp) fft(bp)[0] fft(bp)[1] ]
//      [ sum(bm) fft(bm)[0] fft(bm)[1] ]
//      [ sum(cp) fft(cp)[0] fft(cp)[1] ]
//      [ sum(cm) fft(cm)[0] fft(cm)[1] ]
//      [ sum(dp) fft(dp)[0] fft(dp)[1] ]
//      [ sum(dm) fft(dm)[0] fft(dm)[1] ]  ]
//
// Stage 2 ( D=4) treats this buffer as a 4*6 matrix,
//    [ [ sum(ap) fft(ap)[0] fft(ap)[1] sum(am) fft(am)[0] fft(am)[1] ]
//      [ sum(bp) fft(bp)[0] fft(bp)[1] sum(bm) fft(bm)[0] fft(bm)[1] ]
//      [ sum(cp) fft(cp)[0] fft(cp)[1] sum(cm) fft(cm)[0] fft(cm)[1] ]
//      [ sum(dp) fft(dp)[0] fft(dp)[1] sum(dm) fft(dm)[0] fft(dm)[1] ]  ]
//
//    Then FFTs each column, transposing as it goes.
//
//    The resulting matrix is the 3d FFT of the 2x3x4 input matrix.
//
//    Note as a sanity check that the first element of the final
//    stage's output (DC term) is
//    sum( [ sum(ap) sum(bp) sum(cp) sum(dp) ] )
//     , i.e. the summation of all 24 input elements.
void kiss_fftnd(kiss_fftnd_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout) {
	int i,k;
	const kiss_fft_cpx * bufin=fin;
	kiss_fft_cpx * bufout;

	// arrange it so the last bufout == fout
	if ( st->ndims() & 1 ) {
		bufout = fout;
		if ( fin==fout ) {
			//memcpy( st->tmpbuf(), fin, sizeof(kiss_fft_cpx) * st->dimprod() );
			for ( int i=0; i<st->dimprod(); ++i ) st->tmpbuf()[i]=fin[i];
			bufin = st->tmpbuf();
		}
	} else {
		bufout = st->tmpbuf();
	}

	for ( k=0; k < st->ndims(); ++k ) {
		int curdim = st->dims()[k];
		int stride = st->dimprod() / curdim;

		for ( i=0 ; i<stride ; ++i ) {
			kiss_fft_stride( st->states(k), bufin+i , bufout+i*curdim, stride );
		}

		// toggle back and forth between the two buffers
		if ( bufout == st->tmpbuf() ) {
			bufout = fout;
			bufin = st->tmpbuf();
		} else {
			bufout = st->tmpbuf();
			bufin = fout;
		}
	}
}


///////////////////////////////////////
// Multidim real FFT
///////////////////////////////////////
void kiss_fftndr(kiss_fftndr_cfg st, const kiss_fft_scalar *timedata, kiss_fft_cpx *freqdata) {
	int k1,k2;
	int dimReal = st->dimReal();
	int dimOther = st->dimOther();
	int nrbins = dimReal/2+1;

	kiss_fft_cpx * tmp1 = (kiss_fft_cpx*)st->tmpbuf();
	kiss_fft_cpx * tmp2 = tmp1 + std::max(nrbins,dimOther);

	// timedata is N0 x N1 x ... x Nk real
	// take a real chunk of data, fft it and place the output at correct intervals
	for ( k1=0; k1<dimOther; ++k1 ) {
		kiss_fftr( st->cfg_r(), timedata + k1*dimReal , tmp1 ); // tmp1 now holds nrbins complex points
		for ( k2=0; k2<nrbins; ++k2 ) {
			tmp2[ k2*dimOther+k1 ] = tmp1[k2];
		}
	}

	for ( k2=0; k2<nrbins; ++k2 ) {
		kiss_fftnd(st->cfg_nd(), tmp2+k2*dimOther, tmp1);  // tmp1 now holds dimOther complex points
		for ( k1=0; k1<dimOther; ++k1 ) {
			freqdata[ k1*(nrbins) + k2] = tmp1[k1];
		}
	}
}

void kiss_fftndri(kiss_fftndr_cfg st, const kiss_fft_cpx *freqdata, kiss_fft_scalar *timedata) {
	int k1,k2;
	int dimReal = st->dimReal();
	int dimOther = st->dimOther();
	int nrbins = dimReal/2+1;
	kiss_fft_cpx * tmp1 = (kiss_fft_cpx*)st->tmpbuf();
	kiss_fft_cpx * tmp2 = tmp1 + std::max(nrbins,dimOther);

	for ( k2=0; k2<nrbins; ++k2 ) {
		for ( k1=0; k1<dimOther; ++k1 ) {
			tmp1[k1] = freqdata[ k1*(nrbins) + k2 ];
		}
		kiss_fftnd(st->cfg_nd(), tmp1, tmp2+k2*dimOther);
	}

	for ( k1=0; k1<dimOther; ++k1 ) {
		for ( k2=0; k2<nrbins; ++k2 ) {
			tmp1[k2] = tmp2[ k2*dimOther+k1 ];
		}
		kiss_fftri( st->cfg_r(),tmp1,timedata + k1*dimReal);
	}
}

}
}
