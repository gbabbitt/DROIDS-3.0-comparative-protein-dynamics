#!/usr/bin/perl
#use Tk;
#use Tk::PNG;
#use Tk::JPEG;
#use Tk::Photo;
#use strict;
#use warnings;
#use feature ":5.10";
#use File::Copy;
#use File::Path;
#use List::Util qw( min );
#use List::Util qw( max );
#use List::Util qw(min max);
use Statistics::Descriptive::Full();

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
#print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
#                         in Dynamic Simulations
#
#- visual toolbox for functional evolutionary comparison
#  of molecular dynamic simulation \n\n";

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
      if ($header eq "shape"){$vector_enter = $value;}
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
print("\nparsing training set time series for vibrational data for each amino acid in $refID...\n");
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

print("\nparsing training set time series for vibrational data for each amino acid in $queryID...\n");
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
  close OUT; 
 }

###########################################
###########################################

if ($vector_enter eq 'y'){

mkdir("indAAtrain_vector");

print("\nparsing training set time series for shape data for each amino acid in $refID...\n");
sleep(3);
 for (my $t = 0; $t < $lengthID; $t++){
  open(OUT, ">"."./indAAtrain_vector/vecttimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
  print OUT "class\t"."X\t"."Y\t"."Z\t"."XO\t"."YO\t"."ZO\n";
   for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./atomvect_$refID/vect_$refID"."_aa$t"."_$tt".".txt")||die "could not open vector file "."vect_$refID"."_aa$t"."_$tt".".txt\n";
    my @IN = <IN>;   
    $framecounter = 0;
    @frameXavg = (); @frameYavg = (); @frameZavg = (); @frameXOavg = (); @frameYOavg = (); @frameZOavg = ();
    for (my $a = 0; $a < scalar @IN; $a++) {
	$framecounter = $framecounter +1;
     $INrow = $IN[$a];
     #$DATAseries = substr($INrow, 28, length($INrow) - 28); # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	$frameID = @INrow[1];
     $Xvector = @INrow[2]; $Yvector = @INrow[3]; $Zvector = @INrow[4];
     $origXvector = @INrow[5]; $origYvector = @INrow[6]; $origZvector = @INrow[7];
     #print "AA $t\t"."run $tt\t"."frame $framenumber\t"."X $Xvector\t"."Y $Yvector\t"."Z $Zvector\t"."origX $origXvector\t"."origY $origYvector\t"."origZ $origZvector\n";
     if ($framecounter < $framestep){push (@frameXavg, $Xvector);push (@frameYavg, $Yvector);push (@frameZavg, $Zvector);push (@frameXOavg, $origXvector);push (@frameYOavg, $origYvector);push (@frameZOavg, $origZvector);}
     if ($framecounter == $framestep){
          $statSCORE = new Statistics::Descriptive::Full; # avg X vector
          $statSCORE->add_data (@frameXavg); $meanX = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Y vector
          $statSCORE->add_data (@frameYavg); $meanY = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Z vector
          $statSCORE->add_data (@frameZavg); $meanZ = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg X vector origin
          $statSCORE->add_data (@frameXOavg); $meanXO = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Y vector origin
          $statSCORE->add_data (@frameYOavg); $meanYO = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Z vector origin
          $statSCORE->add_data (@frameZOavg); $meanZO = $statSCORE->mean();
          print OUT "0\t"."$meanX\t"."$meanY\t"."$meanZ\t"."$meanXO\t"."$meanYO\t"."$meanZO\n";
          @frameXavg = (); @frameYavg = (); @frameZavg = (); @frameXOavg = (); @frameYOavg = (); @frameZOavg = ();
          }
     }

     close IN;
 
   }
   close OUT; 
 
 }

for (my $t = 0; $t < $lengthID; $t++){
  open(OUT, ">>"."./indAAtrain_vector/vecttimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
  #print OUT "class\t"."X\t"."Y\t"."X\t"."XO\t"."YO\t"."ZO\n";
   for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./atomvect_$queryID/vect_$queryID"."_aa$t"."_$tt".".txt")||die "could not open vector file "."vect_$queryID"."_aa$t"."_$tt".".txt\n";
    my @IN = <IN>;   
    $framecounter = 0;
    @frameXavg = (); @frameYavg = (); @frameZavg = (); @frameXOavg = (); @frameYOavg = (); @frameZOavg = ();
    for (my $a = 0; $a < scalar @IN; $a++) {
	$framecounter = $framecounter +1;
     $INrow = $IN[$a];
     #$DATAseries = substr($INrow, 28, length($INrow) - 28); # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	$frameID = @INrow[1];
     $Xvector = @INrow[2]; $Yvector = @INrow[3]; $Zvector = @INrow[4];
     $origXvector = @INrow[5]; $origYvector = @INrow[6]; $origZvector = @INrow[7];
     #print "AA $t\t"."run $tt\t"."frame $framenumber\t"."X $Xvector\t"."Y $Yvector\t"."Z $Zvector\t"."origX $origXvector\t"."origY $origYvector\t"."origZ $origZvector\n";
     if ($framecounter < $framestep){push (@frameXavg, $Xvector);push (@frameYavg, $Yvector);push (@frameZavg, $Zvector);push (@frameXOavg, $origXvector);push (@frameYOavg, $origYvector);push (@frameZOavg, $origZvector);}
     if ($framecounter == $framestep){
          $statSCORE = new Statistics::Descriptive::Full; # avg X vector
          $statSCORE->add_data (@frameXavg); $meanX = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Y vector
          $statSCORE->add_data (@frameYavg); $meanY = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Z vector
          $statSCORE->add_data (@frameZavg); $meanZ = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg X vector origin
          $statSCORE->add_data (@frameXOavg); $meanXO = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Y vector origin
          $statSCORE->add_data (@frameYOavg); $meanYO = $statSCORE->mean();
          $statSCORE = new Statistics::Descriptive::Full; # avg Z vector origin
          $statSCORE->add_data (@frameZOavg); $meanZO = $statSCORE->mean();
          print OUT "1\t"."$meanX\t"."$meanY\t"."$meanZ\t"."$meanXO\t"."$meanYO\t"."$meanZO\n";
          @frameXavg = (); @frameYavg = (); @frameZavg = (); @frameXOavg = (); @frameYOavg = (); @frameZOavg = ();
          }
     }

     close IN;
 
   }
   close OUT; 
 
 }

mkdir("indAAtrain_fluxvector");
print "\n\ncombining atom fluctuation and protein shape (i.e. vector ) information\n\n";
sleep(1);

for (my $t = 0; $t < $lengthID; $t++){
  open(OUT, ">"."./indAAtrain_fluxvector/fluxvectorAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
  #print OUT "class\t"."value\t"."datatype\t"."run\n";
  print OUT "class\t"."value\t"."datatype\t"."run\n";
  for (my $tt = 0; $tt < $framegroups; $tt++){
        open(IN1, "<"."./indAAtrain/fluxtimeAA_$refID"."_$t".".txt")||die "could not open flux file "."fluxtimeAA_$refID"."_$t".".txt\n";
        open(IN2, "<"."./indAAtrain_vector/vecttimeAA_$refID"."_$t".".txt")||die "could not open vector file "."vecttimeAA_$refID"."_$t".".txt\n";
        my @IN1 = <IN1>;
        my @IN2 = <IN2>;
        for (my $a = 0; $a < scalar @IN1; $a++) {
	       $IN1row = $IN1[$a];
            @IN1row = split(/\s+/, $IN1row);
	       $class1 = @IN1row[0];
            $flux = @IN1row[$tt+1];
            #print "AA "."$t"."  run "."$tt"."  class "."$class1"."  flux "."$flux\n";
            if($class1 eq '0' || $class1 eq '1'){print OUT "$class1\t"."$flux\t"."F\t"."$tt\n";}
            }
        for (my $aa = 0; $aa < scalar @IN2; $aa++){
               $IN2row = $IN2[$aa];
               @IN2row = split(/\s+/, $IN2row);
	          $class2 = @IN2row[0];
               $x = @IN2row[1];
               $y = @IN2row[2];
               $z = @IN2row[3];
               $xo = @IN2row[4];
               $yo = @IN2row[5];
               $zo = @IN2row[6];
               #$flux = @IN1row[$tt+1];
               #print "AA "."$t"."  run "."$tt"."  class "."$class2"."  x "."$x"."  y "."$y"."  z "."$z"."  xo "."$xo"."  yo "."$yo"."  zo "."$zo\n";
               if($class2 eq '0' || $class2 eq '1'){print OUT "$class2\t"."$x\t"."X\t"."$tt\n"."$class2\t"."$y\t"."Y\t"."$tt\n"."$class2\t"."$z\t"."Z\t"."$tt\n";}
             }
      close IN1;
      close IN2;
      }
   close OUT;
   }

}
############################################
print("\nparsing is done\n");
sleep(1);
exit;
############################################
############################################

#}  # end control subroutine
###########################################################################################################
###########################################################################################################
