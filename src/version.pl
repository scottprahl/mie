#!/usr/bin/perl -w
#
# gets version from $Makefile and then creates a two new files
#    version.c
#    version.tex
# with the version information and a date stamp.
# if the makefile contains VERSION=1.0.0d (i.e., ends with a 'd' then the
#    version is assumed to be a development release and the time stamp in
#    the $header file is the current time otherwise it gets the last this is a development
#    version and the date becomes 
#
use POSIX qw(strftime);

my $makefile = '../Makefile';
my $header   = 'version.c';
my $tex      = 'version.tex';

my $version  = '0-0-0';

open IN, $makefile or die "could not open <$makefile>";
while (<IN>){
	last if m/VERSION\s*=/;
}
if (/=\s*([-0-9d]+)/) {$version=$1};
close IN;

$time='';
if ($version =~ /d/) {
	$time = strftime "%H:%M:%S %D", localtime;
} else {
	$time = strftime("%H:%M:%S %D", localtime((stat $makefile)[9]));
}

open OUT, ">$header" or die "could not open <$header> for output";
print OUT "/* AUTOMATICALLY generated version header by  <$0> */\n";
print OUT "/* version number is set in                   <$makefile>  */\n";
if ($version =~ /d/) {
	print OUT "/* date corresponds to last time program was compiled.       */\n";
} else {
	print OUT "/* date is the last modification date for     <$makefile>  */\n";
}
print OUT "\n";
print OUT "char *Version = \"$version ($time)\";\n";
close OUT;

open OUT, ">$tex" or die "could not open <$header> for output";
print OUT "%% AUTOMATICALLY generated version header by  <$0>\n";
print OUT "%% version number is set in                   <$makefile>\n";
if ($version =~ /d/) {
	print OUT "%% date corresponds to last time program was compiled.       \n";
} else {
	print OUT "%% date is the last modification date for     <$makefile>  \n";
}
print OUT "\n";
print OUT "\\def\\version{Version $version}\n";
close OUT;
