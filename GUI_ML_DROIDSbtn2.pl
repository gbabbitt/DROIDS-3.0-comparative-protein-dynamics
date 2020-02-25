#!/usr/bin/perl
use Tk;
use Tk::PNG;
use Tk::JPEG;
use Tk::Photo;
#use strict;
#use warnings;
use feature ":5.10";
use File::Copy;
use File::Path;
use List::Util qw( min );
use List::Util qw( max );
use List::Util qw(min max);
use Statistics::Descriptive();

# collect frame data
open (IN, "<"."MDframes.ctl") || die "could not open frames ctl file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "framenumber"){$framenumber = $value;}
      if ($header eq "framestep"){$framestep = $value;}
      if ($header eq "framegroups"){$framegroups = $value;}
      if ($header eq "framegrpfactor"){$framefactor = $value;}
}
close IN;

# specify the path to working directory for Chimera here
open(IN, "<"."paths.ctl") or die "could not find paths.txt file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $path = @INrow[1];
	 if ($header eq "chimera_path"){$chimera_path = $path;}
}
close IN;
print "path to Chimera .exe\t"."$chimera_path\n";

#### This uses a GUI to write the control files needed for the DROIDS scripts ####
print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
                          in Dynamic Simulations

- visual toolbox for functional evolutionary comparison
  of molecular dynamic simulation \n\n";

#print "Enter residue number at the start of both chains\n";
#print "(e.g. enter 389 if starts at THR 389.A) \n";
#print "(e.g. enter 1 if starts at MET 1.A) \n\n";
#my $startN = <STDIN>;
#chop($startN);

#### Declare variables ####
my $queryID = '';
my $refID = '';
my $lengthID = '';
my $testStr = '';
my $cutoffValue = '';
my $homology = '';
my $chainN = '';


# read control files
open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.ctl file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "query"){$queryID = $value;}
      if ($header eq "reference"){$refID = $value;}
      if ($header eq "length"){$lengthID = $value;}
      #if ($header eq "start"){$startN = $value;}
      if ($header eq "cutoff_value"){$cutoffValue = $value;}
      if ($header eq "test_type"){$testStr = $value;}
      if ($header eq "representations"){$repStr = $value;}
      if ($header eq "homology"){$homology = $value;}
      if ($header eq "num_chains"){$chainN = $value;}
}
close IN;

open(IN2, "<"."MDr.ctl") or die "could not find MDr.ctl file\n";
my @IN2 = <IN2>;
for (my $i = 0; $i < scalar @IN2; $i++){
	 my $INrow = $IN2[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "Number_Runs"){$number_runs = $value;}
}
close IN2;

#sleep(1);print "MASK LEARNING TO ONLY SIGNIFICANTLY DIFFERENT DYNAMICS? (type 'on' or 'off')\n\n";
#my $option = <STDIN>;
#chop($option);
#sleep(1);print "SELECT MOVIE VIEWER TYPE (type 'fix' or 'roll')\n\n";
#my $view = <STDIN>;
#chop($view);

#sub ctl {

##############################
#print("appending ctl file...\n");
#open(CTL, '>>', "DROIDS.ctl") or die "Could not open output file";
#print CTL "MLoption\t"."$option\t # mask method (apply to signif KS result...on or off)\n";
#print CTL "MLview\t"."$view\t # movie viewing type\n";
#close CTL;
#print("CTL file made\n");
#sleep(1);
##############################
mkdir("indAAtrain");
#mkdir("./trainingData_$queryID/indAA");
print("\nparsing training set time series for each amino acid in $refID...\n");
sleep(1);
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./indAAtrain/fluxtimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./trainingData_$refID/fluxtime_$refID"."_$tt".".txt")||die "could not open time series file "."fluxtime_$refID"."_$tt.txt\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $DATAseries = substr($INrow, 28, length($INrow) - 28);  # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	$aaID = int($a/4 - 0.1);
     $atomID = @INrow[1];
	$overallAVG = @INrow[2];
     if ($a == 0){  # create headers for file
          @DATAseries = split(/\s+/, $DATAseries);
          $framelistlength = scalar @DATAseries;
          $framelist = '';
          $frame = 0;
          for (my $aa = 0; $aa < $framelistlength; $aa++){$frame = $frame +1; $framelist = "$framelist"."$frame\t";}
          if ($tt == 0){print OUT "class\t"."$framelist\n";}
          next;}
     if ($aaID == $t){print OUT "0\t".$DATAseries;} #print "$tt\t"."$aaID\t"."$atomID\t"."$overallAVG\n";}
     }
   close IN;
  }
   close OUT;
   
 }

print("\nparsing training set time series for each amino acid in $queryID...\n");
sleep(1);
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">>"."./indAAtrain/fluxtimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./trainingData_$queryID/fluxtime_$queryID"."_$tt".".txt")||die "could not open time series file "."fluxtime_$queryID"."_$tt.txt\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $DATAseries = substr($INrow, 28, length($INrow) - 28); # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	$aaID = int($a/4 - 0.1);
     $atomID = @INrow[1];
	$overallAVG = @INrow[2];
     if($atomID eq "AtomicFlx"){next;} # skip header row if encountered
     if ($aaID == $t){print OUT "1\t".$DATAseries;} #print "$tt\t"."$aaID\t"."$atomID\t"."$overallAVG\n";}
     }
   close IN;
  }
   
 }
print("\nparsing is done\n");
sleep(1);
exit;
############################################
############################################

#}  # end control subroutine
###########################################################################################################
###########################################################################################################
