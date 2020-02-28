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


# prompt user - choose best learning model to display
sleep(1);print "\nSELECT MOVIE VIEWER TYPE (type 'fix' or 'roll')\n\n";
my $view = <STDIN>;
chop($view);

sleep(1);print "\nSELECT INTERACTION TYPE (1=protein only | 2=protein-protein | 3=DNA-protein | 4=protein-ligand)\n\n";
my $stype_number = <STDIN>;
chop($stype_number);

if($stype_number == 1){$stype = "protein";}
if($stype_number == 2){$stype = "protprot";}
if($stype_number == 3){$stype = "dna";}
if($stype_number == 4){$stype = "ligand";}

if($stype eq "dna" || $stype eq "ligand"){$orig_queryID = $queryID; $queryID = $refID;}

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################


###########################################################################################################
###########################################################################################################
#sub play1 {

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

########################################

if($view eq "fix" && $stype ne "protprot") {
$attr = "class";
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

#$attr = "dRMSF";
#print("Preparing movie display...\n");
#print("close DROIDS movie windows to exit\n\n");
#my @movies;
#my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
##for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#my $movieStr = join(" ", @movies);
#system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

}

if($view eq "roll" && $stype ne "protprot") {
$attr = "class";     
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
for (my $i = 8; $i < 10; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

#$attr = "dRMSF";     
#print("Preparing movie display...\n");
#print("close DROIDS movie windows to exit\n\n");
#my @movies;
#my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
##for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 8; $i < 10; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#my $movieStr = join(" ", @movies);
#system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");
}

if($view eq "stereo" && $stype ne "protprot") {
$attr = "class";       
print("Preparing movie display...\n");
system("ffmpeg -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1L_6.mp4 -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1R_7.mp4 -filter_complex hstack Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewSTEREO.mp4");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

#$attr = "dRMSF";       
#print("Preparing movie display...\n");
#system("ffmpeg -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1L_6.mp4 -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1R_7.mp4 -filter_complex hstack Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewSTEREO.mp4");
#print("close DROIDS movie windows to exit\n\n");
#my @movies;
#my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
##for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#my $movieStr = join(" ", @movies);
#system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");
}
####################################
# reverse query and reference for prot-prot interaction
if($view eq "fix" && $stype eq "protprot") {
$attr = "class";
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

#$attr = "dRMSF";
#print("Preparing movie display...\n");
#print("close DROIDS movie windows to exit\n\n");
#my @movies;
#my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$orig_queryID"."_$refID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
##for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#my $movieStr = join(" ", @movies);
#system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

}

if($view eq "roll" && $stype eq "protprot") {
$attr = "class";     
print("Preparing movie display...\n");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
for (my $i = 8; $i < 10; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

#$attr = "dRMSF";     
#print("Preparing movie display...\n");
#print("close DROIDS movie windows to exit\n\n");
#my @movies;
#my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
##for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 8; $i < 10; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$orig_queryID"."_$refID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#my $movieStr = join(" ", @movies);
#system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");
}

if($view eq "stereo" && $stype eq "protprot") {
$attr = "class";       
print("Preparing movie display...\n");
system("ffmpeg -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1L_6.mp4 -i Videos/$orig_queryID"."_$refID"."_$repStr"."_$attr"."_$testStr"."_viewS1R_7.mp4 -filter_complex hstack Videos/$orig_queryID"."_$refID"."_$repStr"."_$attr"."_$testStr"."_viewSTEREO.mp4");
print("close DROIDS movie windows to exit\n\n");
my @movies;
my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
#for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
my $movieStr = join(" ", @movies);
system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");

#$attr = "dRMSF";       
#print("Preparing movie display...\n");
#system("ffmpeg -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1L_6.mp4 -i Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_viewS1R_7.mp4 -filter_complex hstack Videos/$orig_queryID"."_$refID"."_$repStr"."_$attr"."_$testStr"."_viewSTEREO.mp4");
#print("close DROIDS movie windows to exit\n\n");
#my @movies;
#my @axes = ('Z1', 'Z2', 'X2', 'X1', 'Y1', 'Y2', 'S1L', 'S1R', 'R1', 'R2');
##for (my $i = 0; $i < 6; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$refID"."_$queryID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#for (my $i = 6; $i < 8; $i++) { $axis = $axes[$i]; $movies[$i] = "Videos/$orig_queryID"."_$refID"."_$repStr"."_$attr"."_$testStr"."_view$axis"."_$i.mp4"; }
#my $movieStr = join(" ", @movies);
#system("x-terminal-emulator -e python DROIDS_gstreamer.py @movies");
}


#}

###########################################################################################################
###########################################################################################################
#####################################################



