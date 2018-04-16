#!/usr/bin/perl -anl

sub hsv2rgb
{
	my ($h, $s, $v) = ($_[0] / 60, $_[1], int( $_[2] * 255 ) );

	my $i = int($h);
	my $f = ( $i % 2 ? $h - $i : 1 - ($h - $i) );
	my $m = int( $v * (1 - $s) );
	my $n = int( $v * ( 1 - $s * $f ) );
		
	my @tripel = ( [$v,$n,$m], [$n,$v,$m], [$m,$v,$n], [$m,$n,$v], [$n,$m,$v], [$v,$m,$n] );

	$i %= 6;
		
	return @{$tripel[$i]};
}

sub rgb2hsv
{
	# normalize colors
	@rgb = map { $_ / 255 } @_;

	# first index: color (r = 0, g = 1, b = 2); second index: value
	my @sorted = sort { $a->[1] <=> $b->[1] } map { [$i++, $_] } @rgb;

	# V = max. color
	my $v = $sorted[-1]->[1];

	# H: depends on which is the max. color
	my $h = 60 * (0 + ($rgb[1] - $rgb[2])/($sorted[-1]->[1]-$sorted[0]->[1]), 2 + ($rgb[2] - $rgb[0])/($sorted[-1]->[1]-$sorted[0]->[1]), 4 + ($rgb[0] - $rgb[1])/($sorted[-1]->[1]-$sorted[0]->[1]) )[$sorted[-1]->[0]];

	# S: normalized ratio between max and min colors
	my $s = ( $sorted[-1]->[1] == 0 ? 0 : ($sorted[-1]->[1] - $sorted[0]->[1]) / $sorted[-1]->[1] );

	return ($h, $s, $v);
}

/created by:/ and $name = (/created by: (.*)$/)[0];

/COLOR_MODEL/ and ( $colormodel = ($F[-1] =~ /\+?([A-Z]+)/)[0], print "<ColorMap name=\"$name\" space=\"$colormodel\">" );

/^[#BFN]|^$/ and next;

for my $i ( 0..1 )
{
	$val = $F[0 + 4 * $i];

	@color = ( $colormodel =~ /HSV/ ? map { $_ / 255 } hsv2rgb( @F[1 + 4 * $i .. 3 + 4 * $i] ) : map { $_ / 255 } @F[1 + 4 * $i .. 3 + 4 * $i] );

	print "\t<Point x=\"$val\" o=\"1\" r=\"$color[0]\" g=\"$color[1]\" b=\"$color[2]\"/>";
}

END { print "</ColorMap>"; }
