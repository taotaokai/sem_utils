#!/usr/bin/perl -anl

use constant PI => atan2( 0, -1 );
use POSIX "fmod";

$, = ' ';

sub deg2rad { $_[0] / 180 * PI }

sub rad2deg { $_[0] * 180 / PI }

sub acos
{
	my $square = $_[0] * $_[0];
	if( $square > 1 )
	{
		$square = 1;
	}

	return( atan2( sqrt(1 - $square), $_[0] ) );
}

# returns the great circle distance in radians
sub calc_distance
{
	my ($lat1, $lon1, $lat2, $lon2) = @_;

        my $a = deg2rad( 90 - $lat1 );
        my $b = deg2rad( 90 - $lat2 );
        my $gamma = deg2rad( abs( $lon2 - $lon1 ) );

        return acos( cos( $a ) * cos( $b ) + sin( $a ) * sin( $b ) * cos( $gamma ) );
}

sub interpolate_distance
{
	my ($lat1, $lon1, $lat2, $lon2, $ratio) = @_;

	my $epsilon2 = 1e-6;
	$lat1 -= $epsilon2 if( $lat1 == 90 );
	$lat2 -= $epsilon2 if( $lat2 == 90 );
	$lat1 += $epsilon2 if( $lat1 == -90 );
	$lat2 += $epsilon2 if( $lat2 == -90 );

        my $a = deg2rad( 90 - $lat1 );
        my $b = deg2rad( 90 - $lat2 );
        my $gamma = deg2rad( abs( $lon2 - $lon1 ) );
	$gamma = $epsilon2 if( 0 == fmod( $gamma, PI) );
	my $c = acos( cos( $a ) * cos( $b ) + sin( $a ) * sin( $b ) * cos( $gamma ) );
	my $beta = acos( ( cos( $b ) - cos( $a ) * cos( $c ) ) / ( sin( $a ) * sin( $c ) ) );

	$c *= $ratio;

	$b = acos( cos( $a ) * cos( $c ) + sin( $a ) * sin( $c ) * cos( $beta ) );
	$gamma = ($gamma <=> 0) * acos( ( cos( $c ) - cos( $a ) * cos( $b ) ) / ( sin( $a ) * sin( $b ) ) );

	$lat2 = 90 - rad2deg( $b );
	$lat2 = (abs($lat2) > 90 ? ($lat2 <=> 0) * 180 - $lat2 : $lat2 );
	#$lon2 = $lon1 + ( ($lon2 - $lon1) <=> 0 ) * rad2deg( $gamma );
	$lon2 = $lon1 + ( (abs($lon2 - $lon1) <= 180) ? ( ($lon2 - $lon1) <=> 0 ) : -( ($lon2 - $lon1) <=> 0 ) ) * rad2deg( $gamma );
	$lon2 = (fmod( $lon2 + ($lon2 <=> 0) * 180, 360) - ($lon2 <=> 0) * 180);

	return ( $lat2, $lon2 );
}

/^>/ and $segment++, next;

if( $plates{"$segment"} )
{
	my $prev_lat = $plates{"$segment"}[-1]{'lat'};
	my $prev_lon = $plates{"$segment"}[-1]{'lon'};

	my $distance = calc_distance( $F[0], $F[1], $prev_lat, $prev_lon );

	#my $points_to_add = int( 0.5 + $distance / deg2rad( 1 ) ) - 1;
	my $points_to_add = 0;
	
	for my $i (1..$points_to_add)
	{
		my ($new_lat, $new_lon) = interpolate_distance( $prev_lat, $prev_lon, $F[1], $F[0], ($i / ($points_to_add + 1) ) );
	
		push( @{ $plates{"$segment"} }, { 'lon' => $new_lon, 'lat' => $new_lat } );
	}
}

push( @{ $plates{"$segment"} }, { 'lon' => $F[0], 'lat' => $F[1] } );

END {

	my $points = 0;

	for $current_segment ( keys %plates )
	{
		$points += scalar( @{$plates{$current_segment}} );
	}

	print "# vtk DataFile Version 3.0";
	print "xy2vtk.pl";
	print "ASCII";
	print "DATASET POLYDATA";
	print "POINTS $points float";

	for $current_segment ( keys %plates )
	{
		for my $i (0..@{$plates{$current_segment}} - 1)
		{
			my $x = cos( deg2rad( $plates{$current_segment}[$i]{'lon'} ) ) * cos( deg2rad( $plates{$current_segment}[$i]{'lat'} ) );
			my $y = sin( deg2rad( $plates{$current_segment}[$i]{'lon'} ) ) * cos( deg2rad( $plates{$current_segment}[$i]{'lat'} ) );
			my $z = sin( deg2rad( $plates{$current_segment}[$i]{'lat'} ) );
		
			print $x, $y, $z;
		}
	}

	print "LINES " . scalar( keys %plates ) . " " . ($points + 1 * scalar( keys %plates ));

	my $offset = 0;

	for my $current_segment ( keys %plates )
	{
		my @arr = @{$plates{$current_segment}};
		print scalar(@{$plates{$current_segment}}), $offset..($#arr + $offset);	# do not close polygon, since it's an open polyline
			
		$offset += scalar(@{$plates{$current_segment}});
	}
}
