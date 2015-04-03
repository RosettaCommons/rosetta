
use strict;
use constant PI    => 4 * atan2(1, 1);

return 1;


sub deep_copy {
	my $this = shift;
	if (not ref $this) {
		$this;
	} elsif (ref $this eq "ARRAY") {
		[map deep_copy($_), @$this];
	} elsif (ref $this eq "HASH") {
		+{map { $_ => deep_copy($this->{$_}) } keys %$this};
	} else { die "what type is $_?" }
}


# rotation from euler angles
sub euler {
	my ($aa, $bb, $gg) = @_;
	my $MM;

	$MM->[0][0] = (-sin($aa)*cos($bb)*sin($gg) + cos($aa)*cos($gg));
	$MM->[0][1] = ( cos($aa)*cos($bb)*sin($gg) + sin($aa)*cos($gg));
	$MM->[0][2] = ( sin($bb)*sin($gg));
	$MM->[1][0] = (-sin($aa)*cos($bb)*cos($gg) - cos($aa)*sin($gg));
	$MM->[1][1] = ( cos($aa)*cos($bb)*cos($gg) - sin($aa)*sin($gg));
	$MM->[1][2] = ( sin($bb)*cos($gg));
	$MM->[2][0] = ( sin($aa)*sin($bb));
	$MM->[2][1] = (-cos($aa)*sin($bb));
	$MM->[2][2] = ( cos($bb));

	return $MM;
}

# my ($X,$Y,$Z,$W)=R2quat($M)
sub R2quat {
	my $R = shift;
	my ($S,$X,$Y,$Z,$W);
	if ( $R->[0][0] > $R->[1][1] && $R->[0][0] > $R->[2][2] )  {
		$S  = sqrt( 1.0 + $R->[0][0] - $R->[1][1] - $R->[2][2] ) * 2;
		$X = 0.25 * $S;
		$Y = ($R->[1][0] + $R->[0][1] ) / $S;
		$Z = ($R->[2][0] + $R->[0][2] ) / $S;
		$W = ($R->[2][1] - $R->[1][2] ) / $S;
	} elsif ( $R->[1][1] > $R->[2][2] ) {
		$S  = sqrt( 1.0 + $R->[1][1] - $R->[0][0] - $R->[2][2] ) * 2;
		$X = ($R->[1][0] + $R->[0][1] ) / $S;
		$Y = 0.25 * $S;
		$Z = ($R->[2][1] + $R->[1][2] ) / $S;
		$W = ($R->[0][2] - $R->[2][0] ) / $S;
	} else {
		$S  = sqrt( 1.0 + $R->[2][2] - $R->[0][0] - $R->[1][1] ) * 2;
		$X = ($R->[0][2] + $R->[2][0] ) / $S;
		$Y = ($R->[2][1] + $R->[1][2] ) / $S;
		$Z = 0.25 * $S;
		$W = ($R->[1][0] - $R->[0][1]) / $S;
	}
	return ($X,$Y,$Z,$W);
}


sub R2angle {
	my $R = shift;
	my $x = $R->[2][1]-$R->[1][2];
	my $y = $R->[0][2]-$R->[2][0];
	my $z = $R->[1][0]-$R->[0][1];
	my $r = sqrt($x*$x + $y*$y + $z*$z);
	my $t = $R->[0][0]+$R->[1][1]+$R->[2][2];
	my $th=atan2($r,$t-1);
	return $th;
}


# my ($R)=R2quat($X,$Y,$Z,$W)
sub quat2R {
	my ($X,$Y,$Z,$W) = @_;
	my $xx = $X * $X; my $xy = $X * $Y; my $xz = $X * $Z;
	my $xw = $X * $W; my $yy = $Y * $Y; my $yz = $Y * $Z;
	my $yw = $Y * $W; my $zz = $Z * $Z; my $zw = $Z * $W;
	my $R = [ [ 1 - 2 * ( $yy+$zz ) ,     2 * ( $xy-$zw ) ,     2 * ( $xz+$yw ) ] ,
	          [     2 * ( $xy+$zw ) , 1 - 2 * ( $xx+$zz ) ,     2 * ( $yz-$xw ) ] ,
	          [     2 * ( $xz-$yw ) ,     2 * ( $yz+$xw ) , 1 - 2 * ( $xx+$yy ) ] ];
	return $R;
}

# my ($R)=R2quat($X,$Y,$Z,$W)
sub quatnorm {
	my ($X,$Y,$Z,$W) = @_;
	my $S = sqrt( $X*$X+$Y*$Y+$Z*$Z+$W*$W );
	return [ $X/$S , $Y/$S , $Z/$S , $W/$S ];
}

#####################################
#####################################

# vector addition
sub vadd {
	my ($x, $y) = @_;
	return [ $x->[0]+$y->[0], $x->[1]+$y->[1], $x->[2]+$y->[2] ];
}

# vector subtraction
sub vsub {
	my ($x, $y) = @_;
	return [ $x->[0]-$y->[0], $x->[1]-$y->[1], $x->[2]-$y->[2] ];
}

# mult vector by scalar
sub vscale {
	my ($x, $y) = @_;
	return [ $x*$y->[0], $x*$y->[1], $x*$y->[2] ];
}

# "min mod"
sub minmod {
	my ($x,$y) = @_;
	my $r = fmod($x,$y);
	if ($r < -fabs( $y/2.0 ) ) { $r += fabs( $y ); }
	elsif ($r >  fabs( $y/2.0 ) ) { $r -= fabs( $y ); }
	return $r;
}

# vector min-modulus
sub vminmod {
	my ($x,$y) = @_;
	return [ minmod($x->[0],$y->[0]), minmod($x->[1],$y->[1]), minmod($x->[2],$y->[2]) ];
}


#####################################
#####################################

# raise a matrix to a power
# dumb way of doing it
sub mpow {
	my ($mat, $pow) = @_;
	my $matpow = $mat;
	foreach my $i (2..$pow) {
		$matpow = mmult( $mat, $matpow );
	}
	return $matpow;
}

# matrix x vector mult
sub mapply {
	my ($rotmat, $cart) = @_;
	my $out = [0, 0, 0];
	my ($i, $j);
	for ($i=0; $i < 3; ++$i) {
		for ($j=0; $j < 3; ++$j) {
			$out->[$i] += $rotmat->[$i][$j] * $cart->[$j];
		}
	}
	return $out;
}

# matrix x matrix mult
sub mmult {
	my ($m1, $m2) = @_;
	my $out = [ [0,0,0], [0,0,0], [0,0,0] ];
	my ($i, $j, $k);
	for ($i=0; $i<3; ++$i) {
		for ($j=0; $j<3; ++$j) {
			for ($k=0; $k<3; ++$k) {
				$out->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];
			}
		}
	}
	return $out;
}


# matrix inversion
sub minv {
	my $M = shift;
	my $Minv = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $D = $M->[0][0] * ( $M->[1][1]*$M->[2][2] - $M->[2][1]*$M->[1][2] ) -
		    $M->[0][1] * ( $M->[1][0]*$M->[2][2] - $M->[1][2]*$M->[2][0] ) +
		    $M->[0][2] * ( $M->[1][0]*$M->[2][1] - $M->[1][1]*$M->[2][0] );
	if ($D == 0)  {
		print STDERR "ERROR ... Inversion of singular matrix!\n";
		exit -1;
	}

	$Minv->[0][0] =  ($M->[1][1]*$M->[2][2]-$M->[1][2]*$M->[2][1])/$D;
	$Minv->[0][1] = -($M->[0][1]*$M->[2][2]-$M->[0][2]*$M->[2][1])/$D;
	$Minv->[0][2] =  ($M->[0][1]*$M->[1][2]-$M->[0][2]*$M->[1][1])/$D;
	$Minv->[1][0] = -($M->[1][0]*$M->[2][2]-$M->[2][0]*$M->[1][2])/$D;
	$Minv->[1][1] =  ($M->[0][0]*$M->[2][2]-$M->[0][2]*$M->[2][0])/$D;
	$Minv->[1][2] = -($M->[0][0]*$M->[1][2]-$M->[0][2]*$M->[1][0])/$D;
	$Minv->[2][0] =  ($M->[1][0]*$M->[2][1]-$M->[2][0]*$M->[1][1])/$D;
	$Minv->[2][1] = -($M->[0][0]*$M->[2][1]-$M->[0][1]*$M->[2][0])/$D;
	$Minv->[2][2] =  ($M->[0][0]*$M->[1][1]-$M->[0][1]*$M->[1][0])/$D;

	return $Minv;
}

# vector norm
sub vnorm {
	my $x = shift;
	return sqrt( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# vector norm^2
sub vnorm2 {
	my $x = shift;
	return ( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# cart distance
sub vdist {
	my ($x1, $x2) = @_;
	return sqrt (
	           ($x1->[0]-$x2->[0])*($x1->[0]-$x2->[0]) +
	           ($x1->[1]-$x2->[1])*($x1->[1]-$x2->[1]) +
	           ($x1->[2]-$x2->[2])*($x1->[2]-$x2->[2])
	            );
}

# are two transformations the inverso of one another
sub is_inverse {
	my $tol = 1e-8;
	my ( $R_i,$T_i, $R_j,$T_j ) = @_;

	my $testR = mmult( $R_i , $R_j );
	my $testT = vadd( mapply( $R_j,$T_i ) , $T_j );

	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
	           square($testR->[0][1])   + square($testR->[0][2])   + square($testR->[1][2]) +
	           square($testR->[1][0])   + square($testR->[2][0])   + square($testR->[2][1]);
	my $errT = square($testT->[0])+square($testT->[1])+square($testT->[2]);

#print " (0) $errR   $errT\n";
	if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}

#
sub is_identity {
	my $testR = shift;
	my $tol = 1e-3;
	if (scalar( @_ ) >= 1) {
		$tol = shift;
	}
	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
	           square($testR->[0][1])   + square($testR->[0][2])   + square($testR->[1][2]) +
	           square($testR->[1][0])   + square($testR->[2][0])   + square($testR->[2][1]);

	if ($errR < $tol) { return 1; }
	return 0;
}

# is the transform (Rn,Tn) equivalent to the transform (Ri,Ti)->(Rj,Tj)
sub is_equivalent {
	my $tol = 1e-8;
	my ( $R_n,$T_n, $R_i,$T_i, $R_j,$T_j ) = @_;

	my $R_i_inv = minv( $R_i );
	my $T_i_inv = [ -$T_i->[0], -$T_i->[1], -$T_i->[2] ];

	my $R_test = mmult( $R_i_inv, $R_j );
	my $T_test = mapply( $R_i_inv, vsub( $T_j, $T_i ) );

	my $errR = square($R_test->[0][0]-$R_n->[0][0])+square($R_test->[0][1]-$R_n->[0][1])+square($R_test->[0][2]-$R_n->[0][2]) 
	         + square($R_test->[1][0]-$R_n->[1][0])+square($R_test->[1][1]-$R_n->[1][1])+square($R_test->[1][2]-$R_n->[1][2]) 
	         + square($R_test->[2][0]-$R_n->[2][0])+square($R_test->[2][1]-$R_n->[2][1])+square($R_test->[2][2]-$R_n->[2][2]);
	my $errT = square($T_test->[0]-$T_n->[0])+square($T_test->[1]-$T_n->[1])+square($T_test->[2]-$T_n->[2]);

	##
	## don't think we need to check T ...
	if ($errR < $tol) { return 1; }
	#if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}

###########
# f2c/c2f #
###########

sub d2r { return ($_[0]*PI/180); }
sub square { return $_[0]*$_[0]; }

# my ($f2c,$c2f) = crystparams($a,$b,$c,$alpha,$beta,$gamma)
sub crystparams {
	my ($a,$b,$c,$alpha,$beta,$gamma) = @_;

	if ($a*$b*$c == 0) {
		print STDERR "Must provide valid crystal parameters!\n";
		exit -1;
	}

	my $f2c = [ [0,0,0] , [0,0,0] , [0,0,0] ];

	my $ca = cos(d2r($alpha)); my $cb = cos(d2r($beta)); my $cg = cos(d2r($gamma));
	my $sa = sin(d2r($alpha)); my $sb = sin(d2r($beta)); my $sg = sin(d2r($gamma));
	$f2c->[0][0] = $a;
	$f2c->[0][1] = $b * $cg;
	$f2c->[0][2] = $c * $cb;
	$f2c->[1][1] = $b * $sg;
	$f2c->[1][2] = $c * ($ca - $cb*$cg) / $sg;
	$f2c->[2][2] = $c * $sb * sqrt(1.0 - square(($cb*$cg - $ca)/($sb*$sg)));

	my $c2f = minv($f2c);

	return ($f2c,$c2f);
}
