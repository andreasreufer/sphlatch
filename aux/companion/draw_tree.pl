#!/usr/bin/perl

open(fp_w,">tree.sm");
open(fp_r,"sys.hier");

sub draw_level {
	$l_idx = $_[0];
	$r_idx = $_[1];
	$x = $_[2];
	$y = $_[3];
	$w = $_[4];
	$h = $_[5];
	$l_siz = $_[6];
	$r_siz = $_[7];
#	printf "$l_idx $r_idx $x $y $w $h\n";
	printf "$l_siz $r_siz\n";
	printf fp_w ("relocate %e %e\n", $x-$w/2,$y);
	printf fp_w ("draw %e %e\n", $x+$w/2,$y);
	printf fp_w ("relocate %e %e\n", $x-$w/2,$y);
	printf fp_w ("draw %e %e\n", $x-$w/2,$y-$h);
	printf fp_w ("expand $l_siz\n");
	printf fp_w ("dot\n");
	printf fp_w ("expand 1.001\n");
#	printf fp_w ("label %i\n", $l_idx);
	printf fp_w ("relocate %e %e\n", $x+$w/2,$y);
	printf fp_w ("draw %e %e\n", $x+$w/2,$y-$h);
	printf fp_w ("expand $r_siz\n");
	printf fp_w ("dot\n");
	printf fp_w ("expand 1.001\n");
#	printf fp_w ("label %i\n", $r_idx);
	$n=$#part_idx+1;
	$part_idx[$n]=$l_idx;
	$part_x[$n]=$x-$w/2;
	$part_y[$n]=$y-$h;
	$part_w[$n]=0.45*$w;
	$n=$#part_idx+1;
	$part_idx[$n]=$r_idx;
	$part_x[$n]=$x+$w/2;
	$part_y[$n]=$y-$h;
	$part_w[$n]=0.45*$w;
}

$line = <fp_r>;
$line = <fp_r>;
$min = 2e33;
$min_per = 6;
$max_per = 0;
for ($i=0;$i<5;$i++) {
	$line = <fp_r>;
	@elem = split/\s+/,$line;
	$per = $elem[12];
#	printf "$per\n";
	$mass_p = 2e33*$elem[2]*$elem[5];
	if ($mass_p < $min) {
		$min = $mass_p;
	}
	if ($per < $min_per) {
		$min_per = $per;
		printf "$min_per\n";
	}
	if ($per > $max_per) {
		$max_per = $per;
		printf "$max_per\n";
	}
	$sum += $per;
}
printf "$min_per $max_per\n";
close(fp_r);
open(fp_r,"sys.hier");
$line = <fp_r>;
$line = <fp_r>;
$x_pos = 0.5;
$y_pos = 0.9;
$wid = 0.6;
#$h = 0.2;
printf fp_w ("dev x11\n");
printf fp_w ("expand 1.001\n");
printf fp_w ("limits 0 1 0 1\n");
#printf fp_w ("box\n");
printf fp_w ("ptype 9 3\n");
printf fp_w ("relocate 0.5 0.9\n");
#$sun_r = log(7e10)/5;
$sun_m = log(2e33);
printf ("$sun_m\n");
printf fp_w ("expand 5\n");
printf fp_w ("dot\n");
printf fp_w ("expand 1.001\n");
for ($i=0;$i<5;$i++) {
	$line = <fp_r>;
	@elem = split/\s+/,$line;
	$hei = .1+.3*(log($elem[12])-log($min_per))/(log($max_per)-log($min_per));
#	$hei = .2;
#	$exp_l = log($elem[4]*1.5e13)/5;
#	$exp_r = log($elem[7]*1.5e13)/5;
	$l_mass = $elem[2]*2e33;
	$exp_l = 1 + 4*(log($l_mass)-log($min))/($sun_m-log($min));
	$exp_r = 1 + 4*(log($elem[5]*$l_mass)-log($min))/($sun_m-log($min));
#	printf "rad_l = $rad_l rad_r = $rad_r\n";
#	printf "$elem[1] $elem[3] @elem[6]\n";
	if ($#part_idx >= 0) {
		for ($j=0;$j<=$#part_idx;$j++) {
#			printf "index = $part_idx[$j] width = $part_w[$j] elem[1] = $elem[1]\n";
			if ($part_idx[$j] == $elem[1]) {
				$x_pos = $part_x[$j];
				$y_pos = $part_y[$j];
				$wid = $part_w[$j];
#				printf "x_pos=  $x_pos\n"
			}
		}
	}
	&draw_level($elem[3],$elem[6],$x_pos,$y_pos,$wid,$hei,$exp_l,$exp_r);
}
