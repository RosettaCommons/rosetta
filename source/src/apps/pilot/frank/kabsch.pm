#################
# symmgp lookup #
#################

use strict;
use POSIX qw(ceil floor fmod fabs);
use List::Util qw(max min);

return 1;

####################################
####################################

# my ($R,$rmsd, $Ycom, $Ycom_to_Xcom) = rms_align( $x,$y );
sub rms_align {
	my ($X,$Y) = @_;

	my ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $X, $Y );
	my ($U,$residual) = calculate_rotation_matrix($R,$E0);
 
	my $rmsd = sqrt( fabs($residual*2.0/$nlist) );
	return ($U,$rmsd,$mov_com, $mov_to_ref);
}

# my $com = recenter( $x );
sub recenter {
	my ($X) = @_;
	my $Natms = scalar(@{ $X });
	my $com = [0,0,0];

	foreach my $n (0..$Natms-1) {
		foreach my $i (0..2) {
			$com->[$i] += $X->[$n][$i];
		}
	}

	foreach my $i (0..2) {
		$com->[$i] /= $Natms;
	}

	foreach my $n (0..$Natms-1) {
		foreach my $i (0..2) {
			$X->[$n][$i] -= $com->[$i];
		}
	}
	return $com;
}


# normalize($a)
sub normalize {
	my $a = shift;
	my $b = sqrt($a->[0]*$a->[0] + $a->[1]*$a->[1] + $a->[2]*$a->[2]);
	if ($b > 1e-6) {
		$a->[0] /= $b; $a->[1] /= $b; $a->[2] /= $b;
	}
}


# my $a_dot_b = dot($a,$b)
sub dot {
	my ($a,$b) = @_;
	return ($a->[0]*$b->[0] + $a->[1]*$b->[1] + $a->[2]*$b->[2]);
}


# my $a = cross ( b , c )
sub cross {
	my ($b,$c) = @_;
	my $a = [ $b->[1]*$c->[2] - $b->[2]*$c->[1] ,
	          $b->[2]*$c->[0] - $b->[0]*$c->[2] ,
	          $b->[0]*$c->[1] - $b->[1]*$c->[0] ];
	return $a;
}


# ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $ref_xlist, $mov_xlist )
sub setup_rotation {
	my ( $ref_xlist, $mov_xlist ) = @_;

	my $nlist = min( scalar(@{ $ref_xlist }) , scalar(@{ $mov_xlist }) );
	my $ref_com = [0,0,0];
	my $mov_com = [0,0,0];
	my $mov_to_ref = [0,0,0];

	foreach my $n (0..$nlist-1) {
		foreach my $i (0..2) {
			$mov_com->[$i] += $mov_xlist->[$n][$i];
			$ref_com->[$i] += $ref_xlist->[$n][$i];
		}
    }
	foreach my $i (0..2) {
		$mov_com->[$i] /= $nlist;
		$ref_com->[$i] /= $nlist;
		$mov_to_ref->[$i] = $ref_com->[$i] - $mov_com->[$i];
	}


	# shift mov_xlist and ref_xlist to centre of mass */
	foreach my $n (0..$nlist-1) {
		foreach my $i (0..2) {
			$mov_xlist->[$n][$i] -= $mov_com->[$i];
			$ref_xlist->[$n][$i] -= $ref_com->[$i];
		}
    }

	# initialize
	my $R = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $E0 = 0.0;

	foreach my $n (0..$nlist-1) {
		# E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n)
		foreach my $i (0..2) {
		  $E0 +=  $mov_xlist->[$n][$i] * $mov_xlist->[$n][$i]
				+ $ref_xlist->[$n][$i] * $ref_xlist->[$n][$i];
		}
	
		# R[i,j] = sum(over n): y(n,i) * x(n,j)
		foreach my $i (0..2) {
			foreach my $j (0..2) {
			   $R->[$i][$j] += $mov_xlist->[$n][$i] * $ref_xlist->[$n][$j];
			}
		}
	}
	$E0 *= 0.5;

	return ($nlist,$mov_com, $mov_to_ref, $R, $E0);
}

# helper funct
sub j_rotate {
	my ($a,$i,$j,$k,$l,$s,$tau) = @_;
	my $g = $a->[$i][$j];
	my $h = $a->[$k][$l];
	$a->[$i][$j] = $g-$s*($h+$g*$tau);
    $a->[$k][$l] = $h+$s*($g-$h*$tau); 
}

# ($d,$v,$nrot) = jacobi3($a)
#    computes eigenval and eigen_vec of a real 3x3
#    symmetric matrix. On output, elements of a that are above
#    the diagonal are destroyed. d[1..3] returns the
#    eigenval of a. v[1..3][1..3] is a matrix whose
#    columns contain, on output, the normalized eigen_vec of a. 
#    n_rot returns the number of Jacobi rotations that were required
sub jacobi3 {
	my $a = shift;

	my $v = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $b = [ $a->[0][0] , $a->[1][1] , $a->[2][2] ];
	my $d = [ $a->[0][0] , $a->[1][1] , $a->[2][2] ];
	my $z = [0,0,0];
	my $n_rot = 0;
	my $thresh = 0;

	# 50 tries!
	foreach my $count (0..49) {

		# sum off-diagonal elements
		my $sum = fabs($a->[0][1])+fabs($a->[0][2])+fabs($a->[1][2]);

		# if converged to machine underflow
		if ($sum == 0.0) {
	      return($d,$v,$n_rot);
		}

		# on 1st three sweeps..
		my $thresh = 0;
		if ($count < 3) {
			$thresh = $sum * 0.2 / 9.0;
		}

		foreach my $i (0,1) {
			foreach my $j ($i+1..2) {
				my $g = 100.0 * fabs($a->[$i][$j]);
		
				# after four sweeps, skip the rotation if
				# the off-diagonal element is small
				if ( $count > 3  
				      && fabs($d->[$i])+$g == fabs($d->[$i]) 
				      && fabs($d->[$j])+$g == fabs($d->[$j]) ) {
					$a->[$i][$j] = 0.0;
				} elsif (fabs($a->[$i][$j]) > $thresh) {
					my $h = $d->[$j] - $d->[$i];
					my ($t,$s,$tau,$theta);

					if (fabs($h)+$g == fabs($h)) {
						$t = $a->[$i][$j] / $h;
					} else {
						$theta = 0.5 * $h / ($a->[$i][$j]);
						$t = 1.0 / ( fabs($theta) + sqrt(1.0 + $theta*$theta) );
						if ($theta < 0.0) { $t = -$t; }
					}
			
					my $c = 1.0 / sqrt(1 + $t*$t);
					$s = $t * $c;
					$tau = $s / (1.0 + $c);
					$h = $t * $a->[$i][$j];
			
					$z->[$i] -= $h;
					$z->[$j] += $h;
					$d->[$i] -= $h;
					$d->[$j] += $h;
			
					$a->[$i][$j] = 0.0;

					foreach my $k (0..$i-1) {
						j_rotate($a, $k, $i, $k, $j, $s, $tau);
					}
					foreach my $k ($i+1..$j-1) {
						j_rotate($a, $i, $k, $k, $j, $s, $tau);
					}
					foreach my $k ($j+1..2) {
						j_rotate($a, $i, $k, $j, $k, $s, $tau);
					}
					foreach my $k (0..2) {
						j_rotate($v, $k, $i, $k, $j, $s, $tau);
					}
					$n_rot++;
				}
			}
		}
			
		foreach my $i (0..2) {
			$b->[$i] += $z->[$i];
			$d->[$i] = $b->[$i];
			$z->[$i] = 0.0;
		}
	}

	print STDERR "WARNING: Too many iterations in jacobi3!  You're bad and you should feel bad.\n";
	exit -1;
}


# ($eigen_vec, $eigenval) = diagonalize_symmetric( $matrix )
sub diagonalize_symmetric {
	my $matrix = shift;
	my ($eigenval,$vec,$n_rot) = jacobi3($matrix);

	# sort solutions by eigenval
	foreach my $i (0..2) {
		my $k = $i;
		my $val = $eigenval->[$i];

		foreach my $j ($i+1..2) {
			if ($eigenval->[$j] >= $val) {
				$k = $j;
				$val = $eigenval->[$k];
			}
		}

		if ($k != $i) {
			$eigenval->[$k] = $eigenval->[$i];
			$eigenval->[$i] = $val;
			foreach my $j (0..2) {
				$val = $vec->[$j][$i];
				$vec->[$j][$i] = $vec->[$j][$k];
				$vec->[$j][$k] = $val;
			}
		}
	}

	# transpose
	my $eigen_vec = [ [$vec->[0][0],$vec->[1][0],$vec->[2][0]] , 
	                  [$vec->[0][1],$vec->[1][1],$vec->[2][1]] , 
	                  [$vec->[0][2],$vec->[1][2],$vec->[2][2]] ];
	return ($eigen_vec, $eigenval);
}


# ($U,$residual) = calculate_rotation_matrix($R,$E0)
sub calculate_rotation_matrix {
	my ($R,$E0) = @_;

	my $RtR = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $left_eigenvec = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $right_eigenvec = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $eigenval = 0;

	 # Rt <- transpose of R
	my $Rt = [ [$R->[0][0],$R->[1][0],$R->[2][0]] , 
	           [$R->[0][1],$R->[1][1],$R->[2][1]] , 
	           [$R->[0][2],$R->[1][2],$R->[2][2]] ];

	# make symmetric RtR = Rt X R
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			$RtR->[$i][$j] = 0.0;
			foreach my $k (0..2) {
				$RtR->[$i][$j] += $Rt->[$k][$i] * $R->[$j][$k];
			}
		}
	}

	($right_eigenvec, $eigenval) = diagonalize_symmetric( $RtR );

	# right_eigenvec's should be an orthogonal system but could be left
	#   or right-handed. Let's force into right-handed system.
	$right_eigenvec->[2] = cross($right_eigenvec->[0], $right_eigenvec->[1]);

	# From the Kabsch algorithm, the eigenvec's of RtR
	#   are identical to the right_eigenvec's of R.
	#   This means that left_eigenvec = R x right_eigenvec
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			$left_eigenvec->[$i][$j] = dot($right_eigenvec->[$i], $Rt->[$j]);
		}
	}

	foreach my $i (0..2) {
		normalize($left_eigenvec->[$i]);
	}

	# Force left_eigenvec[2] to be orthogonal to the other vectors.
	# First check if the rotational matrices generated from the
	#   orthogonal eigenvectors are in a right-handed or left-handed
	#   coordinate system - given by sigma. Sigma is needed to
	#   resolve this ambiguity in calculating the RMSD.
	my $sigma = 1.0;
	my $v = cross($left_eigenvec->[0], $left_eigenvec->[1]);
	if (dot($v, $left_eigenvec->[2]) < 0.0) {
		$sigma = -1.0;
	}
	foreach my $i (0..2) {
	    $left_eigenvec->[2][$i] = $v->[$i];
	}

	# calc optimal rotation matrix U that minimises residual
	my $U = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			foreach my $k (0..2) {
				$U->[$i][$j] += $left_eigenvec->[$k][$i] * $right_eigenvec->[$k][$j];
    		}
		}
	}

	my $residual = $E0 - sqrt(fabs($eigenval->[0]))
	                   - sqrt(fabs($eigenval->[1]))
	                   - $sigma * sqrt(fabs($eigenval->[2]));

	return ($U,$residual);
}

#######
# end #
#######
