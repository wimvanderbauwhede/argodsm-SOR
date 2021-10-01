#!/usr/bin/env perl
use v5.20;
use strict;
use warnings;

if (scalar @ARGV < 2) {
    die "usage: $0 dsmNX dsmNY [niters im jm km]\n";
}
my ($dsmNX,$dsmNY,@rest) = @ARGV;
my $niters = 60;
if (@rest) {
    $niters = shift @rest;
}
my $im=64;
if (@rest) {
    $im = shift @rest;
}
my $jm=64;
if (@rest) {
    $jm = shift @rest;
}
my $km=64;
if (@rest) {
    $km = shift @rest;
}

my $sor_params_h = <<"ENDS";
#ifndef __SOR_PARAMS_H__
#define __SOR_PARAMS_H__
const int im = $im;
const int jm = $jm;
const int km = $km;
const int dsmNX = $dsmNX;
const int dsmNY = $dsmNY;
const int niters = $niters;
#endif

ENDS

open my $SPH, '>', 'sor_params.h' or die $!;
print $SPH $sor_params_h;
close $SPH;

if (not -d 'results') {
    mkdir 'results';
}

my $np = $dsmNX*$dsmNY;
my $oversubscribe = $np>8? '--oversubscribe' : '';
say('scons');
system('scons');
say("mpirun $oversubscribe -n $np ./sor-c");
# --mca mca_base_env_list \"ARGO_WRITE_BUFFER_WRITE_BACK_SIZE=128;ARGO_WRITE_BUFFER_SIZE=2048;ARGO_ALLOCATION_POLICY=1\" --mca btl vader,self --mca oob_tcp_if_include 127.0.0.1/
# ;ARGO_WRITE_BUFFER_SIZE=2048;ARGO_ALLOCATION_POLICY=1
# mpirun --mca btl vader,self --mca oob_tcp_if_include 127.0.0.1/ -n :
# --mca mca_base_env_list \"ARGO_WRITE_BUFFER_WRITE_BACK_SIZE=128\" 
system("mpirun  $oversubscribe -n $np ./sor-c") ;
