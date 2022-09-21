#!/usr/bin/perl

use warnings;
use strict;

#$! is a special variable that conveys the error message

while (my $line = <>) { #when my is used, the declared variable is local instead of global
chomp $line;
if ($line =~ /^\>/) { #if has fasta >
my $new_file = substr $line,1;
$new_file = (split /[ \t]/, $new_file)[0];
$new_file .= ".fa"; #append .fa to the variable
open (chr_file, ">$new_file") or die "error writing $new_file: $!";
}
print chr_file "$line\n";
}
close chr_file;