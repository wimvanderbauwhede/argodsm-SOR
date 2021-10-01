#!/usr/bin/perl
use warnings;
use strict;
our $W=1;
our $DBG=0;
=pod
This script generates specialised subroutines for a set of array dimensions, types and access modes because Fortran is not polymorphic.
The generated subroutines are

    oclMake${dim}D${type}Array${mode}Buffer
    oclWrite${dim}D${type}ArrayBuffer    
    oclRead${dim}D${type}ArrayBuffer
    
=cut

# Dimensions
my @dims = (1..7); # Fortran arrays are limited to 7 dimensions
# Types
my @c_types= qw(float double char short int long);
my @types = qw(Real Integer); # TODO: support Complex
my @kinds = qw(1 2 4 8 );# not supporting kind=16 as there is no standard C equivalent
# Modes
my @modes = qw(Read Write ReadWrite);

my %ftypes =(Float => 'real', Int => 'integer', Double => 'real(8)', Long => 'integer(8)');
my %wordsizes = (Float => 4, Double => 8, Int => 4, Long => 8);

sub gen_szstr {
	my $dim=shift;
	my @insts = map {"sz($_)" } (1..$dim);
	my $szstr=join(', ',@insts);
	return $szstr;
}

open my $IN, '<', 'dsm-C_TEMPL.h';
open my $OUT, '>', 'dsm-C.h';
        print $OUT "// !!! Don't edit this file!!! Edit dsmArray_TEMPL.f95 and run $0 !!!\n";

while (my $line = <$IN> ){
    print $OUT $line unless $line=~/^\s*#pragma GEN\s+WrapperSubs/i;
    if ($line=~/^\s*#pragma GEN\s+WrapperSubs/i) {
        print $OUT "\n#define GEN_WrapperSubs\n\n";
        print $OUT "\n// Generated API routines\n\n";
  
        for my $c_type (@c_types) {
                    my $code = "
// Scalar $c_type

void dsmalloc${c_type}f_(DSMScalarPtr,${c_type}*) ;
void dsmdelete${c_type}f_(DSMScalarPtr);
void dsmwrite${c_type}f_(DSMScalarPtr,${c_type}*);
void dsmread${c_type}f_(DSMScalarPtr,${c_type}*);

// Array of $c_type

void dsmalloc${c_type}arrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdelete${c_type}arrayf_(DSMArrayPtr) ;
void dsmwrite${c_type}arrayf_(DSMArrayPtr,int64_t*,${c_type}*);
void dsmread${c_type}arrayf_(DSMArrayPtr,int64_t*,${c_type}*);
// TODO! void dsmmemcpy${c_type}arrayf_(DSMArrayPtr, int64_t*, $c_type*, int64_t*);

// Local array of $c_type

void dsmalloc${c_type}localarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdelete${c_type}localarrayf_(DSMArrayPtr) ;
void dsmwrite${c_type}localarrayf_(DSMArrayPtr,int64_t*,${c_type}*);
void dsmread${c_type}localarrayf_(DSMArrayPtr,int64_t*,${c_type}*);
// TODO! void dsmmemcpy${c_type}localarrayf_(DSMArrayPtr, int64_t*, $c_type*, int64_t*);

        ";
                print $OUT $code;
        }

my $c_type = 'pointer';
print $OUT "
// Array of $c_type

void dsmalloc${c_type}arrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdelete${c_type}arrayf_(DSMArrayPtr) ;
void dsmwrite${c_type}arrayf_(DSMArrayPtr,int64_t*,int64_t*);
void dsmread${c_type}arrayf_(DSMArrayPtr,int64_t*,int64_t*);
// TODO! void dsmmemcpy${c_type}arrayf_(DSMArrayPtr, int64_t*, int64_t*, int64_t*);

// Local array of $c_type

void dsmalloc${c_type}localarrayf_(DSMArrayPtr, int64_t* ) ;
void dsmdelete${c_type}localarrayf_(DSMArrayPtr) ;
void dsmwrite${c_type}localarrayf_(DSMArrayPtr,int64_t*,int64_t*);
void dsmread${c_type}localarrayf_(DSMArrayPtr,int64_t*,int64_t*);
// TODO! void dsmmemcpy${c_type}localarrayf_(DSMArrayPtr, int64_t*, int64_t*, int64_t*);
";


    } # line
} # loop over template source
close $IN;
close $OUT;

sub create_idx_prod { (my $ndims) = @_;
    my @cs = map {'p%'.$_.'sz'} ('i' .. 'o');
    my @ncs = @cs[0..$ndims-1];
    return join('*', @ncs);
}

sub create_idx_set { (my $ndims) = @_;
    my @cs =  ('i' .. 'o');
    my @ncs = @cs[0..$ndims-1];
    return join(',', @ncs);
}
#(k+a%klh-1)*a%isz*a%jsz+(j+a%jlh)*a%isz+(i+a%ilh)
sub create_idx_expr { (my $ndims) = @_;
    my @idxs = ('i' .. 'o');
    my @cs =  map {'('.$_.'+a%'.$_.'lh)'} @idxs;
    my @ncs = @cs[0..$ndims-1];
    my @nncs=();
    for my $idx (0 .. $ndims-1) {
        my @szs=($ncs[$idx]);        
        if ($idx>0) {
            for my $jdx (0 .. $idx-1) {
                push @szs, 'a%'.$idxs[$jdx].'sz';
            }                             
        }
        push @nncs, join('*', @szs); 
    }


    return join(' + &'."\n".'            ', @nncs);
}

sub ctype { (my $ftype, my $kind, my $use_sized_types) = @_;
if (not defined $use_sized_types) {
    $use_sized_types=0;
}
$kind*=1;
if ($ftype eq 'Real') {
    if ($kind <=4 ) {
        return 'float'
    } elsif ($kind==8) {
        return 'double'        
    } elsif ($kind==16) {
        warn "Fortran 128-bit real to C long double may result in loss of precision!\n";
        return 'long double'; 
    } else {
        die "only 1,2,4,8,16 bytes supported.\n";
    }
} elsif  ($ftype eq 'Integer') {
    if ($kind ==1 ) {
        return $use_sized_types ? 'int8_t' : 'char'; 
    } elsif ($kind==2) {
        return $use_sized_types ? 'int16_t' : 'short';
    } elsif ($kind==4) {
        return $use_sized_types ? 'int32_t' : 'int';
    } elsif ($kind==8) {
        return $use_sized_types ? 'int64_t' : 'long';       
    } else {
        die "only 1,2,4,8 bytes supported.\n";
    }
}

}
