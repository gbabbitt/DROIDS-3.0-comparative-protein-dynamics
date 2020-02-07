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
#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("maxDEMON - feature learning"); # Titles the main window
$mw->setPalette("gray75");
$image = $mw->Photo(-file => "demon.jpg");

# ML method Frame

my $methodFrame = $mw->Frame(	-label => "MACHINE LEARNING METHOD",
				-relief => "groove",
				-borderwidth => 2                    
				);
	my $bnpCheck = $methodFrame->Checkbutton( -text => "basic non-parametric (K nearest neighbors)",
						-foreground => 'navy',
                              #-value=>"bnp",
						-variable=>\$method_bnp
						);
     my $distCheck = $methodFrame->Checkbutton( -text => "probabilistic methods (naive Bayes, LDA, QDA)",
						-foreground => 'navy',
                              #-value=>"dist",
						-variable=>\$method_dist
						);
     my $kernCheck = $methodFrame->Checkbutton( -text => "kernel methods (support vector machine - SVM)",
						-foreground => 'navy',
                              #-value=>"kern",
						-variable=>\$method_kern
						);
     my $ensCheck = $methodFrame->Checkbutton( -text => "ensemble methods on trees (boosting, random forest)",
						-foreground => 'navy',
                              #-value=>"ens",
						-variable=>\$method_ens
						);
# ML option Frame
my $optionFrame = $mw->Frame(	-label => "MASK OPTIONS",
				-relief => "groove",
				-borderwidth => 2
				);
	my $offRadio = $optionFrame->Radiobutton( -text => "apply learners to entire protein",
						-foreground => 'navy',
                              -value=>"off",
						-variable=>\$option
						);
	my $onRadio = $optionFrame->Radiobutton( -text => "learn only where dynamics differ",
						-foreground => 'navy',
                              -value=>"on",
						-variable=>\$option
						);
	
# structure type option Frame
my $structureFrame = $mw->Frame(	-label => "STRUCTURE TYPE BEING DEPLOYED",
				-relief => "groove",
				-borderwidth => 2
				);
	my $proteinRadio = $structureFrame->Radiobutton( -text => "protein only (single or multi-chain)",
						-foreground => 'darkred',
                              -value=>"protein",
						-variable=>\$stype
						);
     
     my $dnaRadio = $structureFrame->Radiobutton( -text => "protein bound to DNA",
						-foreground => 'darkred',
                              -value=>"dna",
						-variable=>\$stype
						);
	my $ligandRadio = $structureFrame->Radiobutton( -text => "protein bound to small molecule ligand",
						-foreground => 'darkred',
                              -value=>"ligand",
						-variable=>\$stype
						);
     my $protprotRadio = $structureFrame->Radiobutton( -text => "protein bound to another protein",
						-foreground => 'darkred',
                              -value=>"protprot",
						-variable=>\$stype
						);

# movie viewing option Frame
my $viewFrame = $mw->Frame(	-label => "MOVIE VIEWING OPTIONS",
				-relief => "groove",
				-borderwidth => 2
				);
	my $fixRadio = $viewFrame->Radiobutton( -text => "monoscopic XYZ - normal 2D",
						-foreground => 'darkred',
                              -value=>"fix",
						-variable=>\$view
						);
     
     my $rollRadio = $viewFrame->Radiobutton( -text => "monoscopic ROLLING - normal 2D",
						-foreground => 'darkred',
                              -value=>"roll",
						-variable=>\$view
						);
	my $stereoRadio = $viewFrame->Radiobutton( -text => "3D stereoscopic - single perspective",
						-foreground => 'darkred',
                              -value=>"stereo",
						-variable=>\$view
						);
	
# Color Type Frame
my $seqFrame = $mw->Frame(	-label => "PROTEIN DIVERGENCE COLOR SCHEME",
				-relief => "groove",
				-borderwidth => 2
				);
     my $col1Radio = $seqFrame->Radiobutton(-text=>"'stoplight' color scheme",
						-foreground => 'navy',
                              -value=>"c1",
						-variable=>\$colorScheme
                              );
	my $col2Radio = $seqFrame->Radiobutton(-text=>"'temperature' color scheme",
						-foreground => 'navy',
                              -value=>"c2",
						-variable=>\$colorScheme
                              );
     
# frameCount Frame				

#my $frameFrame = $mw->Frame();
#		my $frameLabel = $frameFrame->Label(-text=>"number of movie frame sets for video: ");
#		my $frameEntry = $frameFrame->Entry(-borderwidth => 2,
#					-relief => "groove",
#					-textvariable=>\$frameCount
#                    );

$frameCount = $framenumber;


# variant movie Frame				

my $variantFrame = $mw->Frame();
		my $variantLabel = $variantFrame->Label(-text=>"PDB ID of variant to render: ");
		my $variantEntry = $variantFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$varselect
					);
# Buttons

my $stopButton = $mw -> Button(-text => "exit maxDEMON", 
				-command => \&stop,
                -background => 'gray45',
                -foreground => 'white'
				); # Creates a exit button
my $ctlButton = $mw -> Button(-text => "control file and training set processing", 
				-command => \&ctl,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $mdButton = $mw -> Button(-text => "run new MD simulation for comparative analysis", 
				-command => \&md,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $mlstatsButton = $mw -> Button(-text => "deploy machine learner(s) on MD simulations", 
				-command => \&mlstats,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $canonButton = $mw -> Button(-text => "identify conserved dynamics and variant impacts", 
				-command => \&canon,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $image1Button = $mw -> Button(-text => "show significantly conserved dynamics", 
				-command => \&image1,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $movie1Button = $mw -> Button(-text => "render movies showing machine learners in space and time", 
				-command => \&movie1,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $image2Button = $mw -> Button(-text => "show significant mutation and/or binding impact", 
				-command => \&image2,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
#my $movie2Button = $mw -> Button(-text => "render movies showing significant mutation and/or binding impact", 
#				-command => \&movie2,
 #               -background => 'gray45',
  #              -foreground => 'white'
   #             ); # Creates a ctl file button
my $play1Button = $mw -> Button(-text => "play machine learning classification movies", 
				-command => \&play1,
                -background => 'gray45',
                -foreground => 'white'
                ); # Creates a ctl file button
my $demonButton = $mw -> Button(
                -command => \&stop,
                -background => 'gray45',
                -foreground => 'white',
                -image => $image
				); # Creates an image
#### Organize GUI Layout ####
$demonButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$methodFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$structureFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$optionFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$mdButton->pack(-side=>"top",
			-anchor=>"s"
			);
$ctlButton->pack(-side=>"top",
			-anchor=>"s"
			);
$mlstatsButton->pack(-side=>"top",
			-anchor=>"s"
			);
$canonButton->pack(-side=>"top",
			-anchor=>"s"
			);
#$frameFrame->pack(-side=>"top",
#		-anchor=>"s"
#		);
#$frameLabel->pack(-side=>"left");
#$frameEntry->pack(-side=>"left");
$image1Button->pack(-side=>"top",
			-anchor=>"s"
			);
$movie1Button->pack(-side=>"top",
			-anchor=>"s"
			);
$seqFrame->pack(-side=>"top",
		-anchor=>"n"
		);

$col1Radio->pack(-anchor=>"w");
$col2Radio->pack(-anchor=>"w");

$variantFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$variantLabel->pack(-side=>"left");
$variantEntry->pack(-side=>"left");
$image2Button->pack(-side=>"top",
			-anchor=>"s"
			);
#$movie2Button->pack(-side=>"top",
#			-anchor=>"s"
#			);
$viewFrame->pack(-side=>"top",
		-anchor=>"s"
		);
$play1Button->pack(-side=>"top",
			-anchor=>"s"
			);
$bnpCheck->pack();
$distCheck->pack();
$kernCheck->pack();
$ensCheck->pack();
$offRadio->pack();
$onRadio->pack();
$fixRadio->pack();
$rollRadio->pack();
$stereoRadio->pack();
$proteinRadio->pack();
$protprotRadio->pack();
$dnaRadio->pack();
$ligandRadio->pack();
MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################
sub stop {exit;}
########################################################################################
########################################################################################
sub ctl {



##############################
print("appending ctl file...\n");
open(CTL, '>>', "DROIDS.ctl") or die "Could not open output file";
print CTL "MLmethod\t"."$method\t # machine learning method\n";
print CTL "MLoption\t"."$option\t # mask method (apply to signif KS result...on or off)\n";
print CTL "MLview\t"."$view\t # movie viewing type\n";
close CTL;
print("CTL file made\n");
sleep(1);
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

############################################
############################################

}  # end control subroutine
###########################################################################################################
###########################################################################################################
sub md {
if($stype eq 'protein'|| $stype eq 'protprot'){system("perl GUI_MLMD_DROIDS.pl");}    
if($stype eq 'dna'){system("perl GUI_MLMD_DROIDSdp.pl");}   
if($stype eq 'ligand'){system("perl GUI_MLMD_DROIDSlp.pl");}     
}
###########################################################################################################
###########################################################################################################
sub mlstats {

sleep(1);
print "BUG NOTE: if following error is encountered\n";
print "'cant open temp test file'\n";
print "then simply close the maxDemon GUI and reopen it\n";
print "at the terminal with 'perl GUI_ML_DROIDS.pl' and continue\n\n";
sleep(1);

if ($method_kern == 1 || $method_other == 1){
# prompt user - choose best learning model to display
sleep(1);print "CHOOSE KERNEL TYPE FOR SVM (type 'linear', 'polynomial', 'laplace' or 'RBF')\n\n";
my $kerntype_enter = <STDIN>;
chop($kerntype_enter);

if($kerntype_enter eq 'linear'){$kerntype = 'vanilladot';}
if($kerntype_enter eq 'polynomial'){$kerntype = 'polydot';}
if($kerntype_enter eq 'laplace'){$kerntype = 'laplacedot';}
if($kerntype_enter eq 'RBF'){$kerntype = 'rbfdot';}
elsif($kerntype_enter eq ''){$kerntype = 'vanilladot';}
sleep(1);
}


for (my $g = 0; $g < 2; $g++){
if ($g == 0) {open(MUT, "<"."copy_list.txt");}
if ($g == 1) {open(MUT, "<"."variant_list.txt");}
my @MUT = <MUT>;
#print @MUT;
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      if ($g == 1 && $p == 1){next;}
      if ($g == 1 && $p == 2){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      print "\nmaking control file for $fileIDq.pdb\n";
      sleep(2);
         
# initialize classpositionHISTO data files (zero values) - for R plots

# KNN method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOknn.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# NB method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOnb.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# LDA method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOlda.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# QDA method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOqda.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# SVM method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOsvm.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
# RFOR method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOrfor.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# ADABOOST method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOada.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;

############################################
############################################
if ($method_bnp == 1){
     print "running KNN classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassKNNtemp");     
 for (my $p = 1; $p<=$framefactor; $p++){
     $add = $p*$framegroups-$framegroups;
     $total = $framefactor*$framegroups;
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nKNN classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   #print Rinput "print(test)\n";
   # define subset of dataD to match size of dataT
   if ($p > 1){
   print Rinput "test <- test[,$add:$total]\n";
   #print Rinput "print(test)\n";
   }
   if ($p == 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassKNNtemp/classAA_$fileIDq"."_$r.txt')\n";
   }
   if ($p > 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassKNNtemp/classAA_$fileIDq"."_$r.txt', append = TRUE)\n";
   }
   print Rinput "for(i in 1:length(train)){
   train_slice <- train[,i, drop=FALSE]
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   lengthN <- length(ref_train_slice)
   Kval <- round(sqrt(lengthN))
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   knn.pred <- knn(train_slice, test_slice, class, k=Kval)
   myKNN1 <- as.numeric(knn.pred[1])
   myKNN2 <- as.numeric(knn.pred[2])
   myKNN3 <- as.numeric(knn.pred[3])
   myKNN4 <- as.numeric(knn.pred[4])
   my_avgKNN = (myKNN1+myKNN2+myKNN3+myKNN4)/4-1
   if(my_avgKNN > 0.5){my_KNN = 1}
   if(my_avgKNN < 0.5){my_KNN = 0}
   if(my_avgKNN == 0.5){my_KNN = 0.5}
   print (i)
   print(knn.pred)
   print (my_KNN)
   print (delta_rmsf)
   }\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
 } # end iterations 
 mkdir("./testingData_$fileIDq/indAAclassKNN");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassKNN/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassKNNtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns"){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns"){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassKNN/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassKNNtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:"){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
 
  
 
}

############################################
if ($method_dist == 1){ # METHOD IS PERMANENTLY DISABLED (it duplicates QDA)...change to 1 to enable it again
     print "running naive Bayes classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassNBtemp");     
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nNB classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(e1071)\n";
   print Rinput "library(doParallel)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   $c1 = "dataT\$c1";
   $c2 = "dataT\$c2";
   $c3 = "dataT\$c3";
   $c4 = "dataT\$c4";
   $c5 = "dataT\$c5";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n";
   print Rinput "test <- dataD\n";
   #print Rinput "print(train)\n";
   #print Rinput "print(test)\n";
   print Rinput "my_bind <- vector()
   classes <- vector()
   for(i in 1:length(test)){
     new <- class
     my_bind = cbind(new, classes)
     classes <- my_bind
     }
   classes <- as.data.frame(classes)
   stack_class <- stack(classes)
   names(stack_class)<-c(\"class\", \"slice\")
   stack_class_slice = stack_class[,1, drop=FALSE]
   stack_train <- stack(train)
   names(stack_train)<-c(\"c1\", \"slice\")
   stack_train_slice = stack_train[,1, drop=FALSE]\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassNBtemp/classAA_$fileIDq"."_$r.txt')\n";
   print Rinput "numCores <- detectCores()\n";
   print Rinput "print(numCores)\n";
   #print Rinput "cl <- makeCluster(numCores, type='FORK')\n";  
   print Rinput "registerDoParallel(numCores)\n";
   #print Rinput "for(i in 1:length(test)){
   print Rinput "foreach(i=1:length(test), .combine=rbind) %dopar% {
   train_slice = stack_train_slice
   class_slice = stack_class_slice
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   names(test_slice)<-c(\"c1\")
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   ref_train_slice <- train_slice[class == 0,]
   xy <- data.frame(class_slice, train_slice)
   nb.xy   <- naiveBayes(as.factor(class)~., data=xy)
   nb.pred <- predict(nb.xy, as.data.frame(test_slice), type='class')
   myNB1 <- as.numeric(nb.pred[1])
   myNB2 <- as.numeric(nb.pred[2])
   myNB3 <- as.numeric(nb.pred[3])
   myNB4 <- as.numeric(nb.pred[4])
   my_avgNB = (myNB1+myNB2+myNB3+myNB4)/4-1
   if(my_avgNB > 0.5){my_NB = 1}
   if(my_avgNB < 0.5){my_NB = 0}
   if(my_avgNB == 0.5){my_NB = 0.5}
   print (i)
   print(nb.pred)
   print(my_NB)
   print(delta_rmsf)
   }\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
  
 mkdir("./testingData_$fileIDq/indAAclassNB");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassNB/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassNBtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns"){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns"){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassNB/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassNBtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:"){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
 
}


############################################
if ($method_dist == 1){
     print "running linear discriminant analysis classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassLDAtemp");     
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nLDA classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(MASS)\n";
   print Rinput "library(doParallel)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   $c1 = "dataT\$c1";
   $c2 = "dataT\$c2";
   $c3 = "dataT\$c3";
   $c4 = "dataT\$c4";
   $c5 = "dataT\$c5";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n";
   print Rinput "test <- dataD\n";
   #print Rinput "print(train)\n";
   #print Rinput "print(test)\n";
   print Rinput "my_bind <- vector()
   classes <- vector()
   for(i in 1:length(test)){
     new <- class
     my_bind = cbind(new, classes)
     classes <- my_bind
     }
   classes <- as.data.frame(classes)
   stack_class <- stack(classes)
   names(stack_class)<-c(\"class\", \"slice\")
   stack_class_slice = stack_class[,1, drop=FALSE]
   stack_train <- stack(train)
   names(stack_train)<-c(\"c1\", \"slice\")
   stack_train_slice = stack_train[,1, drop=FALSE]\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassLDAtemp/classAA_$fileIDq"."_$r.txt')\n";
   print Rinput "numCores <- detectCores()\n";
   print Rinput "print(numCores)\n";
   #print Rinput "cl <- makeCluster(numCores, type='FORK')\n";  
   print Rinput "registerDoParallel(numCores)\n";
   #print Rinput "for(i in 1:length(test)){
   print Rinput "foreach(i=1:length(test), .combine=rbind) %dopar% {
   train_slice = stack_train_slice
   class_slice = stack_class_slice
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   names(test_slice)<-c(\"c1\")
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   ref_train_slice <- train_slice[class == 0,]
   xy <- data.frame(class_slice, train_slice)
   lda.xy   <- lda(as.factor(class)~., data=xy)
   lda.pred <- predict(lda.xy, as.data.frame(test_slice))\$class
   myLDA1 <- as.numeric(lda.pred[1])
   myLDA2 <- as.numeric(lda.pred[2])
   myLDA3 <- as.numeric(lda.pred[3])
   myLDA4 <- as.numeric(lda.pred[4])
   my_avgLDA = (myLDA1+myLDA2+myLDA3+myLDA4)/4-1
   if(my_avgLDA > 0.5){my_LDA = 1}
   if(my_avgLDA < 0.5){my_LDA = 0}
   if(my_avgLDA == 0.5){my_LDA = 0.5}
   print (i)
   print(lda.pred)
   print(my_LDA)
   print(delta_rmsf)
   }\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
  
 mkdir("./testingData_$fileIDq/indAAclassLDA");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassLDA/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassLDAtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns" && $classvalue <= 1){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassLDA/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassLDAtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
 
}

############################################
if ($method_dist == 1){
     print "running quadratic discriminant analysis classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassQDAtemp");     
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nQDA classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(MASS)\n";
   print Rinput "library(doParallel)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   $c1 = "dataT\$c1";
   $c2 = "dataT\$c2";
   $c3 = "dataT\$c3";
   $c4 = "dataT\$c4";
   $c5 = "dataT\$c5";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n";
   print Rinput "test <- dataD\n";
   #print Rinput "print(train)\n";
   #print Rinput "print(test)\n";
   print Rinput "my_bind <- vector()
   classes <- vector()
   for(i in 1:length(test)){
     new <- class
     my_bind = cbind(new, classes)
     classes <- my_bind
     }
   classes <- as.data.frame(classes)
   stack_class <- stack(classes)
   names(stack_class)<-c(\"class\", \"slice\")
   stack_class_slice = stack_class[,1, drop=FALSE]
   stack_train <- stack(train)
   names(stack_train)<-c(\"c1\", \"slice\")
   stack_train_slice = stack_train[,1, drop=FALSE]\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassQDAtemp/classAA_$fileIDq"."_$r.txt')\n";
   print Rinput "numCores <- detectCores()\n";
   print Rinput "print(numCores)\n";
   #print Rinput "cl <- makeCluster(numCores, type='FORK')\n";  
   print Rinput "registerDoParallel(numCores)\n";
   #print Rinput "for(i in 1:length(test)){
   print Rinput "foreach(i=1:length(test), .combine=rbind) %dopar% {
   train_slice = stack_train_slice
   class_slice = stack_class_slice
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   names(test_slice)<-c(\"c1\")
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   ref_train_slice <- train_slice[class == 0,]
   xy <- data.frame(class_slice, train_slice)
   qda.xy   <- qda(as.factor(class)~., data=xy)
   qda.pred <- predict(qda.xy,  as.data.frame(test_slice))\$class
   myQDA1 <- as.numeric(qda.pred[1])
   myQDA2 <- as.numeric(qda.pred[2])
   myQDA3 <- as.numeric(qda.pred[3])
   myQDA4 <- as.numeric(qda.pred[4])
   my_avgQDA = (myQDA1+myQDA2+myQDA3+myQDA4)/4-1
   if(my_avgQDA > 0.5){my_QDA = 1}
   if(my_avgQDA < 0.5){my_QDA = 0}
   if(my_avgQDA == 0.5){my_QDA = 0.5}
   print (i)
   print(qda.pred)
   print(my_QDA)
   print(delta_rmsf)
   }\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
  
 mkdir("./testingData_$fileIDq/indAAclassQDA");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassQDA/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassQDAtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns" && $classvalue <= 1){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassQDA/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassQDAtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
 
}


############################################
if ($method_kern == 1){
     print "running SVM classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassSVMtemp");
 for (my $p = 1; $p<=$framefactor; $p++){
     $add = $p*$framegroups-$framegroups;
     $total = $framefactor*$framegroups;
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nSVM classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(kernlab)\n";
   print Rinput "library(doParallel)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   $c1 = "dataT\$c1";
   $c2 = "dataT\$c2";
   $c3 = "dataT\$c3";
   $c4 = "dataT\$c4";
   $c5 = "dataT\$c5";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n";
   print Rinput "test <- dataD\n";
   #print Rinput "print(train)\n";
   #print Rinput "print(test)\n";
   #print Rinput "sink('./testingData_$fileIDq/indAAclassSVMtemp/classAA_$fileIDq"."_$r.txt')\n";
   #print Rinput "print(test)\n";
   # define subset of dataD to match size of dataT
   if ($p > 1){
   print Rinput "test <- test[,$add:$total]\n";
   #print Rinput "print(test)\n";
   print Rinput "for(i in 1:$framegroups){
   names(test)[i]<-paste('X',i, sep='')
   }\n";
   #print Rinput "print(test)\n";
   }
   #print Rinput "train <- train[-c(1)]\n"; # drops class column
   if ($p == 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassSVMtemp/classAA_$fileIDq"."_$r.txt')\n";
   }
   if ($p > 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassSVMtemp/classAA_$fileIDq"."_$r.txt', append = TRUE)\n";
   }
   print Rinput "numCores <- detectCores()\n";
   print Rinput "print(numCores)\n";
   #print Rinput "cl <- makeCluster(numCores, type='FORK')\n";  
   print Rinput "registerDoParallel(numCores)\n";
   #print Rinput "for(i in 1:length(train)){
   print Rinput "foreach(i=1:length(train), .combine=rbind) %dopar% {
      train_slice <- train[,i, drop=FALSE]
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   xy <- data.frame(class, train_slice)
   svm.xy   <- ksvm(as.factor(class)~., data=xy, kernel='$kerntype',C=2, cross=5)
   svm.pred <- predict(svm.xy, as.data.frame(test_slice), type='response')
   mySVM1 <- as.numeric(svm.pred[1])
   mySVM2 <- as.numeric(svm.pred[2])
   mySVM3 <- as.numeric(svm.pred[3])
   mySVM4 <- as.numeric(svm.pred[4])
   my_avgSVM = (mySVM1+mySVM2+mySVM3+mySVM4)/4-1
   if(my_avgSVM > 0.5){my_SVM = 1}
   if(my_avgSVM < 0.5){my_SVM = 0}
   if(my_avgSVM == 0.5){my_SVM = 0.5}
   print (i)
   print(svm.pred)
   print(my_SVM)
   print(delta_rmsf)
   }\n";
   #print Rinput "stopCluster(cl)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
 } # end iterations 
 mkdir("./testingData_$fileIDq/indAAclassSVM");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassSVM/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassSVMtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns" && $classvalue <= 1){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassSVM/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassSVMtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
 
 
}


###########################################
if ($method_ens == 1){
    print "running random forest classifier\n";
    mkdir("./testingData_$fileIDq/indAAclassRFORtemp");     
 for (my $p = 1; $p<=$framefactor; $p++){
     $add = $p*$framegroups-$framegroups;
     $total = $framefactor*$framegroups;
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nRFOR classifier learning residue $r on $fileIDq\n\n";
   #sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(randomForest)\n";
   print Rinput "library(doParallel)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class = dataT\$class\n"; # training class
   $c1 = "dataT\$c1";
   $c2 = "dataT\$c2";
   $c3 = "dataT\$c3";
   $c4 = "dataT\$c4";
   $c5 = "dataT\$c5";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n";
   print Rinput "test <- dataD\n";
   #print Rinput "sink('./testingData_$fileIDq/indAAclassRFORtemp/classAA_$fileIDq"."_$r.txt')\n";
   #print Rinput "print(test)\n";
   # define subset of dataD to match size of dataT
   if ($p > 1){
    print Rinput "test <- test[,$add:$total]\n";
    #print Rinput "print(test)\n";
    print Rinput "for(i in 1:$framegroups){
    names(test)[i]<-paste('X',i, sep='')
    }\n";
    #print Rinput "print(test)\n";
   }
   #print Rinput "train <- train[-c(1)]\n"; # drops class column
   if ($p == 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassRFORtemp/classAA_$fileIDq"."_$r.txt')\n";
   }
   if ($p > 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassRFORtemp/classAA_$fileIDq"."_$r.txt', append = TRUE)\n";
   }
   print Rinput "numCores <- detectCores()\n";
   print Rinput "print(numCores)\n";
   #print Rinput "cl <- makeCluster(numCores, type='FORK')\n";  
   print Rinput "registerDoParallel(numCores)\n";
   #print Rinput "for(i in 1:length(train)){
   print Rinput "foreach(i=1:length(train), .combine=rbind) %dopar% {
   train_slice <- train[,i, drop=FALSE]
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   ref_train_slice <- train_slice[class == 0,]
   xy <- data.frame(class, train_slice)
   rf.xy <- randomForest(as.factor(class)~., data=xy, ntree=500)
   rf.pred <- predict(rf.xy, as.data.frame(test_slice), type='response')
   myRF1 <- as.numeric(rf.pred[1])
   myRF2 <- as.numeric(rf.pred[2])
   myRF3 <- as.numeric(rf.pred[3])
   myRF4 <- as.numeric(rf.pred[4])
   my_avgRF = (myRF1+myRF2+myRF3+myRF4)/4-1
   if(my_avgRF > 0.5){my_RF = 1}
   if(my_avgRF < 0.5){my_RF = 0}
   if(my_avgRF == 0.5){my_RF = 0.5} 
   print (i)
   print(rf.pred)
   print(my_RF)
   print(delta_rmsf)
   }\n";
   #print Rinput "stopCluster(cl)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
 } # end iterations 
 mkdir("./testingData_$fileIDq/indAAclassRFOR");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassRFOR/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassRFORtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns" && $classvalue <= 1){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassRFOR/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassRFORtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
     
}


###########################################
if ($method_ens == 1){
    print "running adaptive boosting classifier\n";
    mkdir("./testingData_$fileIDq/indAAclassADAtemp");     
 for (my $p = 1; $p<=$framefactor; $p++){
     $add = $p*$framegroups-$framegroups;
     $total = $framefactor*$framegroups;
 for (my $r = 0; $r<$lengthID; $r++){
   print "\n\nADA classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load plotting libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(ada)\n";
   print Rinput "library(doParallel)\n";
   # read data into R
   print Rinput "dataT = read.table('indAAtrain/fluxtimeAA_$refID"."_$r.txt', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$r.txt', header = TRUE)\n";
   print Rinput "class = dataT\$class\n"; # training class
   $c1 = "dataT\$c1";
   $c2 = "dataT\$c2";
   $c3 = "dataT\$c3";
   $c4 = "dataT\$c4";
   $c5 = "dataT\$c5";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n";
   print Rinput "test <- dataD\n";
   #print Rinput "sink('./testingData_$fileIDq/indAAclassADAtemp/classAA_$fileIDq"."_$r.txt')\n";
   #print Rinput "print(test)\n";
   # define subset of dataD to match size of dataT
   if ($p > 1){
    print Rinput "test <- test[,$add:$total]\n";
    #print Rinput "print(test)\n";
    print Rinput "for(i in 1:$framegroups){
    names(test)[i]<-paste('X',i, sep='')
    }\n";
    #print Rinput "print(test)\n";
   }
   #print Rinput "train <- train[-c(1)]\n"; # drops class column
   if ($p == 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassADAtemp/classAA_$fileIDq"."_$r.txt')\n";
   }
   if ($p > 1){
   print Rinput "sink('./testingData_$fileIDq/indAAclassADAtemp/classAA_$fileIDq"."_$r.txt', append = TRUE)\n";
   }
   print Rinput "numCores <- detectCores()\n";
   print Rinput "print(numCores)\n";
   #print Rinput "cl <- makeCluster(numCores, type='FORK')\n";  
   print Rinput "registerDoParallel(numCores)\n";
   #print Rinput "for(i in 1:length(train)){
   print Rinput "foreach(i=1:length(train), .combine=rbind) %dopar% {
   train_slice <- train[,i, drop=FALSE]
   ref_train_slice <- train_slice[class == 0,]
   test_slice <- test[,i, drop=FALSE]
   mean_test_slice <- mean(test_slice[,1])
   mean_ref_train_slice <- mean(ref_train_slice)
   delta_rmsf = (mean_test_slice - mean_ref_train_slice)
   ref_train_slice <- train_slice[class == 0,]
   xy <- data.frame(class, train_slice)
   boost.xy <- ada(as.factor(class)~., data=xy, 18)
   boost.pred <- predict(boost.xy, as.data.frame(test_slice), type='vector')
   myBOOST1 <- as.numeric(boost.pred[1])
   myBOOST2 <- as.numeric(boost.pred[2])
   myBOOST3 <- as.numeric(boost.pred[3])
   myBOOST4 <- as.numeric(boost.pred[4])
   my_avgBOOST = (myBOOST1+myBOOST2+myBOOST3+myBOOST4)/4-1
   if(my_avgBOOST > 0.5){my_BOOST = 1}
   if(my_avgBOOST < 0.5){my_BOOST = 0}
   if(my_avgBOOST == 0.5){my_BOOST = 0.5}
   print (i)
   print(boost.pred)
   print(my_BOOST)
   print(delta_rmsf)
   }\n";
   
  #print Rinput "stopCluster(cl)\n";
  print Rinput "sink()\n";
  # write to output file and quit R
  print Rinput "q()\n";# quit R 
  print Rinput "n\n";# save workspace image?
  close Rinput;
  }
 } # end iterations 
 mkdir("./testingData_$fileIDq/indAAclassADA");
 mkdir("./testingData_$fileIDq/indAAdrmsf");
 if ($option eq "on"){  # apply mask...learn only where KS test is signif
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassADA/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n"; 
   open(IN, "<"."./testingData_$fileIDq/indAAclassADAtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
   open(IN2, "<"."./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/adjustKStests".".txt")||die "could not open adjustKStests.txt file\n";
    my @IN = <IN>;
    my @IN2 = <IN2>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     for (my $aa = 0; $aa < scalar @IN2; $aa++) {
        $IN2row = $IN2[$aa];
        @IN2row = split(/\s+/, $IN2row);
        $AApos = @IN2row[0];
        $ns = @IN2row[5];
        if ($AApos == $t && $search eq "Levels:" && $ns eq "ns" && $classvalue <= 1){print OUT "0\n";}
        if ($AApos == $t && $search eq "Levels:" && $ns ne "ns" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";}
        }
     }
   close IN;
   close IN2;
   close OUT;
   close OUT2;
 }
 }
 
 if ($option eq "off"){  # no mask option applied
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAclassADA/classAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(OUT2, ">"."./testingData_$fileIDq/indAAdrmsf/drmsfAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   open(IN, "<"."./testingData_$fileIDq/indAAclassADAtemp/classAA_$fileIDq"."_$t".".txt")||die "could not open temp time series file\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $nextINrow = $IN[$a+1];
     $nextnextINrow = $IN[$a+2];
     @INrow = split(/\s+/, $INrow);
	@nextINrow = split(/\s+/, $nextINrow);
     @nextnextINrow = split(/\s+/, $nextnextINrow);
     $search = @INrow[0];
     $classvalue = @nextINrow[1];
     $dRMSF = @nextnextINrow[1];
     if ($search eq "Levels:" && $classvalue <= 1){print OUT "$classvalue\n"; print OUT2 "$dRMSF\n";} 
     }
   close IN;
   close OUT;
   close OUT2;
 } 
     
 }
     
}


#############################################################################     
#############################################################################
sleep(1);
print "\n machine learning classification is completed\n\n";

print " copying flux profile training plots\n\n";
sleep(1);
my $oldfilename = "./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/DROIDSfluctuationAVGchain.txt";
my $newfilename = "./testingData_$fileIDq/DROIDSfluctuationAVGchain.txt";
copy($oldfilename, $newfilename);

##### make R Plots
if ($method_bnp == 1){
# KNN method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOknn.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassKNN/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}

if ($method_dist == 1){  # this is permanently disabled (replicates QDA)
# NB method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOnb.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassNB/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}

if ($method_dist == 1){
# LDA method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOlda.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassLDA/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}

if ($method_dist == 1){
# QDA method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOqda.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassQDA/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}

if ($method_kern == 1){
# SVM method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOsvm.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassSVM/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}


if ($method_ens == 1){
# RANDOM FOREST method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOrfor.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassRFOR/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}

if ($method_ens == 1){
# ADABOOST method
open (OUT, ">"."./testingData_$fileIDq/classpositionHISTOada.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 0; $r<$lengthID; $r++){
   $AApos = $r;
   @class_values = ();
   open (IN, "<"."./testingData_$fileIDq/indAAclassADA/classAA_$refID"."_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $class_value = @INrow[0];
      if ($class_value == 0 ||$class_value == 0.5 || $class_value == 1){push(@class_values, $class_value);}
      }
      close IN;
      $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
      $statSCORE->add_data (@class_values);
	 $sum_class = $statSCORE->mean();
      if ($AApos ne '' && $sum_class ne ''){print OUT "$AApos\t"."$sum_class\n";}
}
close OUT;
}

#####################################################
sleep(1);
print " plotting class position histogram and flux profiles\n\n";
sleep(1);

open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";

# read data into R
# KNN
print Rinput "datatableKNN = read.table('./testingData_$fileIDq/classpositionHISTOknn.txt', header = TRUE)\n"; 
$AAposition_knn = "datatableKNN\$position"; # AA position
$sum_classifiers_knn = "datatableKNN\$sum"; # sum of classifiers
print Rinput "dataframeKNN = data.frame(pos1=$AAposition_knn, Y1val=$sum_classifiers_knn)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./testingData_$fileIDq/classpositionHISTOnb.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
print Rinput "dataframeNB = data.frame(pos1=$AAposition_nb, Y1val=$sum_classifiers_nb)\n";
# LDA
print Rinput "datatableLDA = read.table('./testingData_$fileIDq/classpositionHISTOlda.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
print Rinput "dataframeLDA = data.frame(pos1=$AAposition_lda, Y1val=$sum_classifiers_lda)\n";
# QDA
print Rinput "datatableQDA = read.table('./testingData_$fileIDq/classpositionHISTOqda.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
print Rinput "dataframeQDA = data.frame(pos1=$AAposition_qda, Y1val=$sum_classifiers_qda)\n";
# SVM
print Rinput "datatableSVM = read.table('./testingData_$fileIDq/classpositionHISTOsvm.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
print Rinput "dataframeSVM = data.frame(pos1=$AAposition_svm, Y1val=$sum_classifiers_svm)\n";
# random forest
print Rinput "datatableRFOR = read.table('./testingData_$fileIDq/classpositionHISTOrfor.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
print Rinput "dataframeRFOR = data.frame(pos1=$AAposition_rfor, Y1val=$sum_classifiers_rfor)\n";
# adaboost
print Rinput "datatableADA = read.table('./testingData_$fileIDq/classpositionHISTOada.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
print Rinput "dataframeADA = data.frame(pos1=$AAposition_ada, Y1val=$sum_classifiers_ada)\n";
# atom flux profiles
print Rinput "datatable2 = read.table('./testingData_$fileIDq/DROIDSfluctuationAVGchain.txt', header = TRUE)\n"; 
$AAposition2 = "datatable2\$pos_ref"; # AA position
$AAlabel = "datatable2\$res_ref"; # AA  identity
$trainingflux_ref = "datatable2\$flux_ref_avg"; # flux profile ref training
$trainingflux_query = "datatable2\$flux_query_avg"; # flux profile query training
print Rinput "dataframe2 = data.frame(pos2=$AAposition2, Y2val=$trainingflux_ref)\n";
print Rinput "dataframe3 = data.frame(pos3=$AAposition2, Y3val=$trainingflux_query)\n";

# lineplots
print Rinput "myplot1 <- ggplot() + ggtitle('       learning performance along protein backbone for MD validation set') + labs(x = 'position (residue number)', y = 'avg classification over time intervals') + geom_line(data = dataframeKNN, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'K nearest neighbors')) + geom_line(data = dataframeNB, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'naive Bayes')) + geom_line(data = dataframeLDA, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'linear discriminant - LDA')) + geom_line(data = dataframeQDA, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'quadratic discriminant - QDA')) + geom_line(data = dataframeSVM, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'support vector machine-SVM')) + geom_line(data = dataframeRFOR, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'random forest')) + geom_line(data = dataframeADA, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'adaboost')) + scale_color_brewer(palette='Set1')\n"; 
print Rinput "myplot2 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = dataframe2, mapping = aes(x = pos2, y = Y2val, color = 'MD training set - reference protein=$refID')) + geom_line(data = dataframe3, mapping = aes(x = pos3, y = Y3val, color = 'MD training set - query protein=$queryID'))\n"; 

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(2, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$fileIDq/learningPLOT.pdf";
copy($oldfilename, $newfilename);	
print " machine learning is complete\n\n";

close MUT;
} #end for loop
} # end outer for loop
####################################################################
for (my $g = 0; $g < 2; $g++){
if ($g == 0) {open(MUT, "<"."copy_list.txt");}
if ($g == 1) {open(MUT, "<"."variant_list.txt");}
my @MUT = <MUT>;
#print @MUT;
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      if ($g == 1 && $p == 1){next;}
      if ($g == 1 && $p == 2){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      print " open PDF file for $fileIDq\n\n";
      sleep(2);
      print " close PDF viewer to continue\n\n";
      system "evince ./testingData_$fileIDq/learningPLOT.pdf\n";
      } #end for loop
close MUT;
print "\n machine learning is complete\n\n";
} # end outer for loop
} # end sub
###########################################################################################################
###########################################################################################################

sub image1 {

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

print "Enter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot'){
   print "Enter number position of N terminal on this chain (default = 0)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }
# create mask for signif Wilks lamda
open(OUT, ">"."./testingData_$queryID/adj_vertpvals_$queryID.txt")||die "could not create mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
open(IN, "<"."./testingData_$queryID/adj_pvals_$queryID.txt")||die "could not open mask file /testingData_$queryID/adj_pvals_$queryID.txt\n";  
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++) {
	$INrow = $IN[$i];
     @INrow = split(/] /, $INrow);
	$trunINrow = $INrow[1];
     #print OUT "$trunINrow";
     @trunINrow = split(/\s+/, $trunINrow);
     for (my $ii = 0; $ii < scalar(@trunINrow); $ii++){
     $value = $trunINrow[$ii];
     print OUT "$value\n";
     }
     }
close IN;
close OUT;

# make learned class attribute files for image rendering
  open(OUT, ">"."./attribute_files/classATTR_$refID"."_mask".".dat")||die "could not create ATTR time series file\n";
  #if ($stype ne "protprot"){open(OUT, ">"."./attribute_files/classATTR_$refID"."_mask".".dat")||die "could not create ATTR time series file\n";}  
  #if ($stype eq "protprot"){open(OUT, ">"."./attribute_files/classATTR_$orig_queryID"."_mask".".dat")||die "could not create ATTR time series file\n";}  
  print OUT "recipient: residues\n";
  print OUT "attribute: class\n";
  print OUT "\n";
  
   for (my $a = 0; $a < $lengthID+$offset+1; $a++){
   if ($a eq '' || $a <= $offset){next;}
   open(IN, "<"."./testingData_$queryID/adj_vertpvals_$queryID.txt")||die "could not open mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
    my @IN = <IN>;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $maskvalue = $INrow[0];
   
     #print "$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     if ($maskvalue == 1){print OUT "\t:"."$pos".".$chainMAP\t"."1\n";}
     #if ($f == $frame && $maskvalue == 1){print OUT "\t:"."$pos\t"."1\n";} # to test mask positioning
     if ($maskvalue == 0){print OUT "\t:"."$pos".".$chainMAP\t"."0\n";}
     
 
 close IN;
 }
 close OUT;

print "\n class attribute files for all frames are created\n";
sleep(1);


print("Preparing static display...\n");
print("close Chimera window to exit\n\n");
if ($homology eq "loose"){$mutType = "tan";}
if ($homology eq "strict"){$mutType = "tan";}
if ($colorScheme eq "c1" ){$colorType = "wb";}
if ($colorScheme eq "c2" ){$colorType = "wb";}
$attr = "class";
$min_val = 0;
$max_val = 1;

if ($stype eq "protein"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_tan.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_gray.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}

if ($stype eq "protprot"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_tan_pp.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_gray_pp.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}

if ($stype eq "dna"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_tan_dp.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_gray_dp.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}

if ($stype eq "ligand"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_tan_lp.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_gray_lp.py	--rep=$repStr --test=$testStr --qID=$orig_queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}    
     
}

###########################################################################################################
###########################################################################################################
sub image2 {

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

print "Enter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure (e.g. PDB 2oob = B)\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot'){
   
   print "Enter number position of N terminal on this chain (default = 0)\n";
   print "(e.g. PDB ID 200b = 45)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }

# create array of names of copies and variants
open(MUT, "<"."copy_list.txt");
my @MUT = <MUT>;
#print @MUT;
@copies = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@copies, $fileIDq);
      } #end for loop
close MUT;
print "\n copies are @copies\n\n";
sleep(1);

# create array of names of copies and variants
open(MUT, "<"."variant_list.txt");
my @MUT = <MUT>;
#print @MUT;
@variants = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      #if ($p == 1){next;}
      #if ($p == 2){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variants, $fileIDq);
      } #end for loop
close MUT;
print "\n variants are @variants\n\n";
sleep(1);



# select variant to render

$variantID = $varselect;

# create mask for non zero impacts on dynamics
open(OUT, ">"."./testingData_$queryID/impact_vert_$variantID.txt")||die "could not create mask file /testingData_$queryID/impact_vert_$variantID.txt\n";  
open(IN, "<"."./testingData_$queryID/impact_$variantID.txt")||die "could not open mask file /testingData_$variantID/impact_$variantID.txt\n";  
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++) {
	$INrow = $IN[$i];
     @INrow = split(/] /, $INrow);
	$trunINrow = $INrow[1];
     #print OUT "$trunINrow";
     @trunINrow = split(/\s+/, $trunINrow);
     for (my $ii = 0; $ii < scalar(@trunINrow); $ii++){
     $value = $trunINrow[$ii];
     print OUT "$value\n";
     }
     }
close IN;
close OUT;
 
# make relative temperature attribute files on impacted areas for image rendering

  open(OUT, ">"."./attribute_files/drmsfATTR_$refID"."_mask".".dat")||die "could not create ATTR time series file\n";  
  print OUT "recipient: residues\n";
  print OUT "attribute: dRMSF\n";
  print OUT "\n";
  @impacts = ();
  for (my $a = 0; $a < $lengthID+$offset+1; $a++){
   if ($a eq '' || $a <= $offset){next;}
   open(IN, "<"."./testingData_$queryID/impact_vert_$variantID.txt")||die "could not open mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
    my @IN = <IN>;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $maskvalue = $INrow[0];
     
     #print "$f\t"."$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     if ($maskvalue > 0){print OUT "\t:"."$pos".".$chainMAP\t"."$maskvalue\n";}
     #if ($f == $frame && $maskvalue > 0){print OUT "\t:"."$pos\t"."1\n";} # test mask positioning
     if ($maskvalue == 0){print OUT "\t:"."$pos".".$chainMAP\t"."0\n";}
     push (@impacts, $maskvalue); 
   close IN;
 }
close OUT;

print "\n class attribute files for all frames are created\n";
sleep(1);

  # render dynamic movies for local flux difference from MD training run averages
print("Rendering 10 movies on various axes...\n");
print("this may take several minutes...\n\n");
print("close Chimera window when 10 movie files appear in movies folder\n\n");
if ($homology eq "loose"){$mutType = "gray50";}
if ($homology eq "strict"){$mutType = "tan";}
if ($colorScheme eq "c1" ){$colorType = "wr";}
if ($colorScheme eq "c2" ){$colorType = "wr";}
$attr = "dRMSF";
$statSCORE = new Statistics::Descriptive::Full; # impact map
$statSCORE->add_data (@impacts);
$min_impact = $statSCORE->min();
$max_impact = $statSCORE->max();
$min_val = $min_impact;
$max_val = $max_impact;




sleep(1);
if ($stype eq "protein"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnrmsf_tan.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnrmsf_gray.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}

if ($stype eq "protprot"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_tan.py	--rep=$repStr --test=$testStr --qID=$refID --rID=$orig_queryID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnclass_gray.py	--rep=$repStr --test=$testStr --qID=$refID --rID=$orig_queryID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}

if ($stype eq "dna"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnrmsf_tan_dp.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnrmsf_gray_dp.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}

if ($stype eq "ligand"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"color_by_attr_learnrmsf_tan_lp.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"color_by_attr_learnrmsf_gray_lp.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
print("\n\n Display complete\n\n");
}    
   
    
     
}

###########################################################################################################
###########################################################################################################

sub movie1 {

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

print "number of frames to render = "."$frameCount\n\n";

# prompt user - choose best learning model to display
sleep(1);print "CHOOSE BEST LEARNER FOR VIDEO (type 'KNN', 'NB', 'LDA', 'QDA', 'SVM', 'ADA' or 'RFOR')\n\n";
my $learner = <STDIN>;
chop($learner);

sleep(1);print "SHOW LEARNERS ONLY WITHIN CONSERVED DYNAMIC REGIONS? (y or n)\n\n";
my $masktype = <STDIN>;
chop($masktype);

print "Enter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot'){
   
   print "Enter number position of N terminal on this chain (default = 0)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }

# create mask for signif Wilks lamda
open(OUT, ">"."./testingData_$queryID/adj_vertpvals_$queryID.txt")||die "could not create mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
open(IN, "<"."./testingData_$queryID/adj_pvals_$queryID.txt")||die "could not open mask file /testingData_$queryID/adj_pvals_$queryID.txt\n";  
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++) {
	$INrow = $IN[$i];
     @INrow = split(/] /, $INrow);
	$trunINrow = $INrow[1];
     #print OUT "$trunINrow";
     @trunINrow = split(/\s+/, $trunINrow);
     for (my $ii = 0; $ii < scalar(@trunINrow); $ii++){
     $value = $trunINrow[$ii];
     print OUT "$value\n";
     }
     }
close IN;
close OUT;

# make learned class attribute files for movie rendering
print "\n creating class attribute files for all frames\n";
#if($stype eq "protprot"){open(IN, "<"."./testingData_$refID/indAAclass$learner/classAA_$refID"."_1".".txt")||die "could not open time series file /testingData_$refID/indAAclass$learner/classAA_$refID"."_1".".txt\n";}
#if($stype ne "protprot"){open(IN, "<"."./testingData_$queryID/indAAclass$learner/classAA_$refID"."_1".".txt")||die "could not open time series file /testingData_$queryID/indAAclass$learner/classAA_$refID"."_1".".txt\n";}
open(IN, "<"."./testingData_$queryID/indAAclass$learner/classAA_$refID"."_1".".txt")||die "could not open time series file /testingData_$queryID/indAAclass$learner/classAA_$refID"."_1".".txt\n";
my @IN = <IN>;
$fgroups = scalar @IN;
close IN;
print "            number of groups of frames\t"."$fgroups\n";
mkdir("./testingData_$queryID/indAAclass_attr");
for (my $f = 1; $f <= $fgroups; $f++){
  open(OUT, ">"."./testingData_$queryID/indAAclass_attr/classATTR_$refID"."_$f".".dat")||die "could not create ATTR time series file\n";
  print OUT "recipient: residues\n";
  print OUT "attribute: class\n";
  print OUT "\n";
  
   for (my $a = 0; $a < $lengthID+$offset+1; $a++){
   if ($a eq '' || $a <= $offset){next;}
   open(IN, "<"."./testingData_$queryID/adj_vertpvals_$queryID.txt")||die "could not open mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
    my @IN = <IN>;
   $rowINDEX = $a-$offset-1;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $maskvalue = $INrow[0];
   #if($stype eq "protprot"){open(IN2, "<"."./testingData_$refID/indAAclass$learner/classAA_$refID"."_$a".".txt")||die "could not open time series file  /testingData_$refID/indAAclass$learner/classAA_$refID"."_$a".".txt\n";}
   #if($stype ne "protprot"){open(IN2, "<"."./testingData_$queryID/indAAclass$learner/classAA_$refID"."_$a".".txt")||die "could not open time series file  /testingData_$queryID/indAAclass$learner/classAA_$refID"."_$a".".txt\n";}
   open(IN2, "<"."./testingData_$queryID/indAAclass$learner/classAA_$refID"."_$rowINDEX".".txt")||die "could not open time series file  /testingData_$queryID/indAAclass$learner/classAA_$refID"."_$a".".txt\n";
    my @IN2 = <IN2>;
    $frame = 0;
    for (my $i = 0; $i < scalar @IN2; $i++) {
	$frame = $i+1;
     $IN2row = $IN2[$i];
     @IN2row = split(/\s+/, $IN2row);
	$value = $IN2row[0];
     #print "$f\t"."$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     if ($f == $frame && $maskvalue == 1){print OUT "\t:"."$pos".".$chainMAP\t"."$value\n";}
     #if ($f == $frame && $maskvalue == 1){print OUT "\t:"."$pos\t"."1\n";} # to test mask positioning
     if ($f == $frame && $maskvalue == 0 && $masktype eq "n"){print OUT "\t:"."$pos".".$chainMAP\t"."$value\n";}
     if ($f == $frame && $maskvalue == 0 && $masktype eq "y"){print OUT "\t:"."$pos".".$chainMAP\t"."0\n";}
     if ($f == $frame && $maskvalue == 0 && $masktype eq "N"){print OUT "\t:"."$pos".".$chainMAP\t"."$value\n";}
     if ($f == $frame && $maskvalue == 0 && $masktype eq "Y"){print OUT "\t:"."$pos".".$chainMAP\t"."0\n";}
     if ($f == $frame && $maskvalue == 0 && $masktype eq "no"){print OUT "\t:"."$pos".".$chainMAP\t"."$value\n";}
     if ($f == $frame && $maskvalue == 0 && $masktype eq "yes"){print OUT "\t:"."$pos".".$chainMAP\t"."0\n";}
     }
 close IN2;
 close IN;
 }
 close OUT;
}
print "\n class attribute files for all frames are created\n";
sleep(1);
 
 # render dynamic movies for local classification method given ML method
print("Rendering 10 movies on various axes...\n");
print("this may take several minutes...\n\n");
print("close Chimera window when 10 movie files appear in movies folder\n\n");
if ($homology eq "loose"){$mutType = "tan";}
if ($homology eq "strict"){$mutType = "tan";}
if ($colorScheme eq "c1" ){$colorType = "wb";}
if ($colorScheme eq "c2" ){$colorType = "wb";}
$attr = "class";
$min_val = 0;
$max_val = 1;
if ($stype ne "protprot"){ 
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"render_movies_learnclass_tan.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"render_movies_learnclass_gray.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
}
if ($stype eq "protprot"){ 
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"render_movies_learnclass_tan.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"render_movies_learnclass_gray.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$frameCount\"\n");}
}
print("\n\n Movies rendered\n");
 
     
}

###########################################################################################################
sub movie2 {

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

print "Enter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure (e.g. PDB 2oob = B)\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot'){
   
   print "Enter number position of N terminal on this chain (default = 0)\n";
   print "(e.g. PDB ID 200b = 45)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }


# create array of names of copies and variants
open(MUT, "<"."copy_list.txt");
my @MUT = <MUT>;
#print @MUT;
@copies = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@copies, $fileIDq);
      } #end for loop
close MUT;
print "\n copies are @copies\n\n";
sleep(1);

# create array of names of copies and variants
open(MUT, "<"."variant_list.txt");
my @MUT = <MUT>;
#print @MUT;
@variants = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      #if ($p == 1){next;}
      #if ($p == 2){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variants, $fileIDq);
      } #end for loop
close MUT;
print "\n variants are @variants\n\n";
sleep(1);



# select variant to render

$variantID = $varselect;

# create mask for non zero impacts on dynamics
open(OUT, ">"."./testingData_$queryID/impact_vert_$variantID.txt")||die "could not create mask file /testingData_$queryID/impact_vert_$variantID.txt\n";  
open(IN, "<"."./testingData_$queryID/impact_$variantID.txt")||die "could not open mask file /testingData_$variantID/impact_$variantID.txt\n";  
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++) {
	$INrow = $IN[$i];
     @INrow = split(/] /, $INrow);
	$trunINrow = $INrow[1];
     #print OUT "$trunINrow";
     @trunINrow = split(/\s+/, $trunINrow);
     for (my $ii = 0; $ii < scalar(@trunINrow); $ii++){
     $value = $trunINrow[$ii];
     print OUT "$value\n";
     }
     }
close IN;
close OUT;
 
# make relative temperature attribute files on impacted areas for movie rendering
print "\n creating class attribute files for all frames\n";
open(IN, "<"."./testingData_$queryID/indAAdrmsf/drmsfAA_$refID"."_1".".txt")||die "could not open time series file drmsfAA_$refID"."_1".".txt\n";
my @IN = <IN>;
$fgroups = scalar @IN;
close IN;
print "            number of groups of frames\t"."$fgroups\n";
mkdir("./testingData_$queryID/indAAdrmsf_attr");
for (my $f = 1; $f <= $fgroups; $f++){
  open(OUT, ">"."./testingData_$queryID/indAAdrmsf_attr/drmsfATTR_$refID"."_$f".".dat")||die "could not create ATTR time series file\n";  
  print OUT "recipient: residues\n";
  print OUT "attribute: dRMSF\n";
  print OUT "\n";
  for (my $a = 0; $a < $lengthID+$offset+1; $a++){
   if ($a eq '' || $a <= $offset){next;}
   open(IN, "<"."./testingData_$queryID/impact_vert_$variantID.txt")||die "could not open mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
    my @IN = <IN>;
   $rowINDEX = $a-$offset-1;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $maskvalue = $INrow[0];
   open(IN2, "<"."./testingData_$queryID/indAAdrmsf/drmsfAA_$refID"."_$rowINDEX".".txt")||die "could not open time series file  drmsfAA_$refID"."_$a".".txt\n";
    my @IN2 = <IN2>;
    $frame = 0;
    for (my $i = 0; $i < scalar @IN2; $i++) {
	$frame = $i+1;
     $IN2row = $IN2[$i];
     @IN2row = split(/\s+/, $IN2row);
	$value = $IN2row[0];
     #print "$f\t"."$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     if ($f == $frame && $maskvalue > 0){print OUT "\t:"."$pos\t"."$value\n";}
     #if ($f == $frame && $maskvalue > 0){print OUT "\t:"."$pos\t"."1\n";} # test mask positioning
     if ($f == $frame && $maskvalue == 0){print OUT "\t:"."$pos\t"."0\n";}
   
   }
   close IN2;
   close IN;
 }
close OUT;
}
print "\n class attribute files for all frames are created\n";
sleep(1);
 
  # render dynamic movies for local flux difference from MD training run averages
print("Rendering 10 movies on various axes...\n");
print("this may take several minutes...\n\n");
print("close Chimera window when 10 movie files appear in movies folder\n\n");
if ($homology eq "loose"){$mutType = "gray50";}
if ($homology eq "strict"){$mutType = "tan";}
if ($colorScheme eq "c1" ){$colorType = "rg";}
if ($colorScheme eq "c2" ){$colorType = "br";}
$attr = "dRMSF";
# find min and max values
print "\nfinding min and max dRMSF\n";
sleep(1);
@dRMSF=();
for (my $a = 0; $a < $lengthID; $a++){
   if ($a eq ''){next;}
   #$rowINDEX = $a-$offset-1;
   open(IN, "<"."./testingData_$queryID/indAAdrmsf/drmsfAA_$refID"."_$a".".txt")||die "could not open time series file  drmsfAA_$refID"."_$a".".txt\n";
    my @IN = <IN>;
    for (my $i = 0; $i < scalar @IN; $i++) {
	$INrow = $IN[$i];
     @INrow = split(/\s+/, $INrow);
	$value = $INrow[0];
     #print "value\t"."$value\n";
     if($value >= -2.5 && $value <= 2.5){push(@dRMSF, $value);}
    }
close IN;
}
#print "dRMSF\n";
#print @dRMSF;
#print "\n";
$statSCORE = new Statistics::Descriptive::Full; # avg abs delta flux
$statSCORE->add_data (@dRMSF);
$min_dRMSF = $statSCORE->min();
$max_dRMSF = $statSCORE->max();
$min_val = $min_dRMSF;
$max_val = $max_dRMSF;

print "color range = "."$min_val\t"."to\t"."$max_val\n";
sleep(1);
if ($stype ne "protprot"){ 
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"render_movies_learndrmsf_tan.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$framenumber\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"render_movies_learndrmsf_gray.py	--rep=$repStr --test=$testStr --qID=$queryID --rID=$refID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$framenumber\"\n");}
}
if ($stype eq "protprot"){
if ($mutType eq "tan"){system("$chimera_path"."chimera --script \"render_movies_learndrmsf_tan_pp.py	--rep=$repStr --test=$testStr --qID=$refID --rID=$orig_queryID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$framenumber\"\n");}
if ($mutType eq "gray50"){system("$chimera_path"."chimera --script \"render_movies_learndrmsf_gray_pp.py	--rep=$repStr --test=$testStr --qID=$refID --rID=$orig_queryID --lengthID=$lengthID --cutoff=$cutoffValue --colorType=$colorType --testType=$testStr --attr=$attr --minVal=$min_val --maxVal=$max_val --frameCount=$framenumber\"\n");}
}
print("\n\n Movies rendered\n");


}
###########################################################################################################
###########################################################################################################
sub play1 {

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


}

###########################################################################################################
###########################################################################################################
#####################################################
sub canon{

# prompt user - choose graphing option
sleep(1);
print "CHOOSE BAR PLOT TYPE - total mutational impact or just conserved regions\n";
print "NOTE: conserved regions option is intended only for very large proteins with proportionally small active sites and no allostery\n";
print "(type 'total', 'conserved' -note: default is 'total')\n\n";
my $bartype = <STDIN>;
chop($bartype);
if ($bartype eq ''){$bartype = 'total';}

print "Enter number position of N terminal on this chain\n";
print "(default = 0)\n\n";
my $startN = <STDIN>;
chop($startN);
if ($startN eq ''){$startN = 0;}

# create array of names of copies and variants
open(MUT, "<"."copy_list.txt");
my @MUT = <MUT>;
#print @MUT;
@copies = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@copies, $fileIDq);
      } #end for loop
close MUT;
print "\n copies are @copies\n\n";
sleep(1);

# create array of names of copies and variants
open(MUT, "<"."variant_list.txt");
my @MUT = <MUT>;
#print @MUT;
@variants = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variants, $fileIDq);
      } #end for loop
close MUT;
print "\n variants are @variants\n\n";
sleep(1);


# create array of names of copies and variants
open(MUT, "<"."variant_label_list.txt");
my @MUT = <MUT>;
#print @MUT;
@variant_labels = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variant_labels, $fileIDq);
      } #end for loop
close MUT;
print "\n variant labels are @variant_labels\n\n";
sleep(1);

print " plotting canonical correlation and flux profiles\n\n";

$window = '';
$cancor_threshold = '';
# prompt user - choose best learning model to display
sleep(1);print "CHOOSE SIZE OF SLIDING WINDOW (e.g. default = 20 = 20 AA window for r value)\n\n";
my $window = <STDIN>;
chop($window);
if ($window eq ''){$window = 20};

# set conserved threshold
$cancor_threshold = 0.3;
sleep(1);print "CHOOSE LEVEL OF SIGNIFICANCE FOR CONSERVED CAN CORR (default = 0.01)\n\n";
my $cancor_threshold = <STDIN>;
chop($cancor_threshold);
if ($cancor_threshold eq ''){$cancor_threshold = 0.01};

open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";
print Rinput "library(CCA)\n";
print Rinput "library(CCP)\n";
print Rinput "library(boot)\n";
#print Rinput "options(scipen=999)\n"; # eliminates printing of sci notation

#############################################
###  analyze copies for plots 1 and 2
#############################################
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
# read test MD data into R
# KNN
print Rinput "datatableKNN = read.table('./testingData_$copyID/classpositionHISTOknn.txt', header = TRUE)\n"; 
$AAposition_knn = "datatableKNN\$position"; # AA position
$sum_classifiers_knn = "datatableKNN\$sum"; # sum of classifiers
print Rinput "dataframeKNN_$v = data.frame(pos1=$AAposition_knn+$startN, Y1val=$sum_classifiers_knn)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./testingData_$copyID/classpositionHISTOnb.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
print Rinput "dataframeNB_$v = data.frame(pos1=$AAposition_nb+$startN, Y1val=$sum_classifiers_nb)\n";
# LDA
print Rinput "datatableLDA = read.table('./testingData_$copyID/classpositionHISTOlda.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
print Rinput "dataframeLDA_$v = data.frame(pos1=$AAposition_lda+$startN, Y1val=$sum_classifiers_lda)\n";
# QDA
print Rinput "datatableQDA = read.table('./testingData_$copyID/classpositionHISTOqda.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
print Rinput "dataframeQDA_$v = data.frame(pos1=$AAposition_qda+$startN, Y1val=$sum_classifiers_qda)\n";
# SVM
print Rinput "datatableSVM = read.table('./testingData_$copyID/classpositionHISTOsvm.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
print Rinput "dataframeSVM_$v = data.frame(pos1=$AAposition_svm+$startN, Y1val=$sum_classifiers_svm)\n";
# random forest
print Rinput "datatableRFOR = read.table('./testingData_$copyID/classpositionHISTOrfor.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
print Rinput "dataframeRFOR_$v = data.frame(pos1=$AAposition_rfor+$startN, Y1val=$sum_classifiers_rfor)\n";
# adaboost
print Rinput "datatableADA = read.table('./testingData_$copyID/classpositionHISTOada.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
print Rinput "dataframeADA_$v = data.frame(pos1=$AAposition_ada+$startN, Y1val=$sum_classifiers_ada)\n";

} #end for loop

# atom flux profiles
print Rinput "datatable2 = read.table('./testingData_$queryID/DROIDSfluctuationAVGchain.txt', header = TRUE)\n"; 
$AAposition2 = "datatable2\$pos_ref"; # AA position
$AAlabel = "datatable2\$res_ref"; # AA  identity
$trainingflux_ref = "datatable2\$flux_ref_avg"; # flux profile ref training
$trainingflux_query = "datatable2\$flux_query_avg"; # flux profile query training
print Rinput "dataframe2 = data.frame(pos2=$AAposition2+$startN, Y2val=$trainingflux_ref)\n";
print Rinput "dataframe3 = data.frame(pos3=$AAposition2+$startN, Y3val=$trainingflux_query)\n";

# canonical correlations
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$copyID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$copyID)\n";
} # end for loop

# output stats file
print Rinput "sink('./testingData_$queryID/cancor_stats.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and new variants\n\n')\n";
print Rinput "cat('there are 6 dimensions and 7 learners\n\n')\n";
print Rinput "cat('cor = overall correlation of all learners in each dimension\n')\n";
print Rinput "cat('coef = slope of each pairwise learner correlation for each dimension\n')\n";
print Rinput "cat('center = mean of each learner\n\n')\n";
print Rinput "cat('1 = KNN, 2 = naive Bayes, 3 = LDA, 4 = QDA, 5 = SVM, 6 = random forest, 7 = adaboost\n\n')\n";
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $copyID\n\n')\n";
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$copyID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$copyID)\n";

print Rinput "sink()\n";

# AA sliding window for canonical correlation
print Rinput "my_cancors_$copyID <- c()\n";
print Rinput "my_pvals_$copyID <- c()\n";
print Rinput "my_ccpos_$copyID <- c()\n";
print Rinput "my_bars_$copyID <- c()\n";
print Rinput "mylength <- length(newMD[,1])\n";
print Rinput "print(mylength)\n";
print Rinput "for(i in 1:(mylength-$window)){
   newMDslice <- newMD[i:(i+$window),1] 
   queryMDslice <- queryMD[i:(i+$window),1];
   queryMDrefslice <- queryMDref[i:(i+$window),1];
   mytest1 <- cancor(newMDslice, queryMDrefslice)
   mycor <- mytest1\$cor
   mytest2 <- cancor(queryMDslice, queryMDrefslice)
   myselfcor <- mytest2\$cor
   N = length(queryMDslice)
   p = 1
   q = 1
   mypvaltest <- p.asym(myselfcor, N, p, q, tstat = 'Wilks')
   mypval <- mypvaltest\$p.value[1]
   if(mypval <= $cancor_threshold){mybar = 1}
   if(mypval > $cancor_threshold){mybar = 0}
   if($v == 0){
   my_bars_$queryID <- c(my_bars_$queryID, mybar)
   }
   my_ccpos_$copyID <- c(my_ccpos_$copyID, i)
   my_cancors_$copyID <- c(my_cancors_$copyID, mycor)
   my_pvals_$copyID <- c(my_pvals_$copyID, mypval)
   }\n";
print Rinput "print(my_cancors_$copyID)\n";
print Rinput "print(my_pvals_$copyID)\n";
print Rinput "my_adjpvals_$copyID <- p.adjust(my_pvals_$copyID, method = 'fdr', n = length(my_pvals_$copyID))\n";
print Rinput "print(my_adjpvals_$copyID)\n";
print Rinput "print(my_ccpos_$copyID)\n";
print Rinput "print(my_bars_$queryID)\n";
print Rinput "my_adjbars_$queryID <- c()\n";
print Rinput "print(my_adjpvals_$copyID\[i])\n";
print Rinput "for(j in 1:(length(my_adjpvals_$copyID))){
   print(my_adjpvals_$copyID\[j])
   if (my_adjpvals_$copyID\[j] <= $cancor_threshold){myadjbar = 1}
   if (my_adjpvals_$copyID\[j] > $cancor_threshold){myadjbar = 0}
   my_adjbars_$queryID <- c(my_adjbars_$queryID, myadjbar) 
   }\n";
print Rinput "print(my_adjbars_$queryID)\n";  
print Rinput "sink('./testingData_$queryID/cancor_stats_$copyID.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $copyID ($window AA sliding window)\n\n')\n";
print Rinput "print(my_cancors_$copyID)\n";
print Rinput "sink()\n";
print Rinput "dataframe4_$copyID = data.frame(pos=my_ccpos_$copyID+$startN, cc=my_cancors_$copyID, bars=my_adjbars_$queryID)\n";
print Rinput "my_overall <- cancor(newMD, queryMD)\n"; # overall
print Rinput "my_ccor_$copyID <- my_overall\$cor[1]\n";
print Rinput "print(my_ccor_$copyID)\n";
# create output mask movie render file for Wilks lambda
print Rinput "sink('./testingData_$queryID/adj_pvals_$copyID.txt')\n";
print Rinput "print(my_adjbars_$queryID)\n";
print Rinput "sink()\n";
} # end for loop

# lineplots
# create plot 1
print Rinput "myplot1 <- ggplot() + ggtitle('learning performance - 2 MD validation sets') + labs(x = 'position (residue number)', y = 'avg class over time intervals') + theme(legend.position = 'none', plot.title = element_text(size = 10))\n"; 
print Rinput "mylist1 <- list()\n";
# create lines for plot 1
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     print Rinput "mylines_$v <- list(geom_line(data = dataframeKNN_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID', ymin = 0, ymax = 1)), geom_line(data = dataframeNB_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeLDA_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeQDA_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeSVM_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeRFOR_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeADA_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')))\n"; 
     print Rinput "mylist1 <- list(mylist1, mylines_$v)\n";
} # end for loop
print Rinput "myplot1final <- myplot1 + mylist1 + scale_color_brewer(palette='Set1')\n"; 
# create plot 2
print Rinput "myplot2 <- ggplot() + ggtitle('cancor for $window AA window + conserved dynamics zones -Wilks lambda') + geom_area(data = dataframe4_$copyID, mapping = aes(x = pos, y = bars, color = 'conserved dynamics'), alpha = 0.8) + labs(x = 'position (residue number)', y = 'local R value for cancor') + theme(legend.title = element_blank())\n"; 
print Rinput "mylist2 <- list()\n";
# create lines for plot 2
for (my $v = 0; $v < scalar(@copies); $v++){
     if ($v == 0) {next;} # skip first line
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     print Rinput "myline_$v <- geom_line(data = dataframe4_$copyID, mapping = aes(x = pos, y = cc, color = '$variantID_label', ymin = 0, ymax = 1))\n"; 
     print Rinput "mylist2 <- list(mylist2, myline_$v)\n";
} # end for loop
print Rinput "myplot2final <- myplot2 + mylist2 + scale_color_brewer(palette='Set2')\n"; 


# exit if no variant list to analyze

if (scalar(@variants) == 1 || scalar(@variants) == 0){
print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1final, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2final, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
#print Rinput "print(myplot2final, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$queryID/learningPLOTcompare.pdf";
copy($oldfilename, $newfilename);	
print " machine learning is complete\n\n";
print " close PDF and txt viewer to continue\n\n";
system "evince ./testingData_$queryID/learningPLOTcompare.pdf\n";

# open can corr results
system "gedit ./testingData_$queryID/cancor_stats.txt\n";
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     system "gedit ./testingData_$queryID/cancor_stats_$copyID.txt\n";
} # end for loop  
     
exit;
}


#############################################
###  analyze variants for plot 3
#############################################


for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     
     
# read test MD data into Rprint Rinput "my_pvals_$variantID <- c()\n";
# KNN
print Rinput "datatableKNN = read.table('./testingData_$variantID/classpositionHISTOknn.txt', header = TRUE)\n"; 
$AAposition_knn = "datatableKNN\$position"; # AA position
$sum_classifiers_knn = "datatableKNN\$sum"; # sum of classifiers
print Rinput "dataframeKNN_$v = data.frame(pos1=$AAposition_knn+$startN, Y1val=$sum_classifiers_knn)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./testingData_$variantID/classpositionHISTOnb.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
print Rinput "dataframeNB_$v = data.frame(pos1=$AAposition_nb+$startN, Y1val=$sum_classifiers_nb)\n";
# LDA
print Rinput "datatableLDA = read.table('./testingData_$variantID/classpositionHISTOlda.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
print Rinput "dataframeLDA_$v = data.frame(pos1=$AAposition_lda+$startN, Y1val=$sum_classifiers_lda)\n";
# QDA
print Rinput "datatableQDA = read.table('./testingData_$variantID/classpositionHISTOqda.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
print Rinput "dataframeQDA_$v = data.frame(pos1=$AAposition_qda+$startN, Y1val=$sum_classifiers_qda)\n";
# SVM
print Rinput "datatableSVM = read.table('./testingData_$variantID/classpositionHISTOsvm.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
print Rinput "dataframeSVM_$v = data.frame(pos1=$AAposition_svm+$startN, Y1val=$sum_classifiers_svm)\n";
# random forest
print Rinput "datatableRFOR = read.table('./testingData_$variantID/classpositionHISTOrfor.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
print Rinput "dataframeRFOR_$v = data.frame(pos1=$AAposition_rfor+$startN, Y1val=$sum_classifiers_rfor)\n";
# adaboost
print Rinput "datatableADA = read.table('./testingData_$variantID/classpositionHISTOada.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
print Rinput "dataframeADA_$v = data.frame(pos1=$AAposition_ada+$startN, Y1val=$sum_classifiers_ada)\n";

} #end for loop

$varcount = scalar(@variants);
# canonical correlations
for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$variantID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$variantID)\n";
} # end for loop

# output stats file
print Rinput "sink('./testingData_$queryID/cancor_stats.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and new variants\n\n')\n";
print Rinput "cat('there are 6 dimensions and 7 learners\n\n')\n";
print Rinput "cat('cor = overall correlation of all learners in each dimension\n')\n";
print Rinput "cat('coef = slope of each pairwise learner correlation for each dimension\n')\n";
print Rinput "cat('center = mean of each learner\n\n')\n";
print Rinput "cat('1 = KNN, 2 = naive Bayes, 3 = LDA, 4 = QDA, 5 = SVM, 6 = random forest, 7 = adaboost\n\n')\n";

print Rinput "my_impact_IDs <- c()\n";
print Rinput "my_impact_sums <- c()\n";
print Rinput "my_impact_sd <- c()\n";
print Rinput "my_impact_n <- c()\n";
print Rinput "my_impact_sums_upperSE <- c()\n";
print Rinput "my_impact_sums_lowerSE <- c()\n";
print Rinput "my_impact_cons_sums <- c()\n";
print Rinput "my_impact_cons_sums_upperSE <- c()\n";
print Rinput "my_impact_cons_sums_lowerSE <- c()\n";
print Rinput "my_impact_cons_sums <- c()\n";
print Rinput "my_impact_cons_sd <- c()\n";
print Rinput "my_impact_cons_n <- c()\n";

for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $copyID\n\n')\n";
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$variantID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$variantID)\n";

print Rinput "sink()\n";

# AA sliding window for canonical correlation
print Rinput "my_cancors_$variantID <- c()\n";
print Rinput "my_pvals_$variantID <- c()\n";
print Rinput "my_ccpos_$variantID <- c()\n";
print Rinput "my_impact_$variantID <- c()\n";
print Rinput "my_impact_sum_$variantID = 0\n";
print Rinput "my_impact_cons_sum_$variantID = 0\n";
print Rinput "my_impact_sd_$variantID = 0\n";
print Rinput "my_impact_cons_sd_$variantID = 0\n";
print Rinput "my_impact_n_$variantID = 0\n";
print Rinput "my_impact_length_$variantID = 0\n";
print Rinput "my_impact_cons_length_$variantID = 0\n";
print Rinput "my_bars_$variantID <- c()\n";
print Rinput "mylength <- length(newMD[,1])\n";
print Rinput "print(mylength)\n";
print Rinput "for(i in 1:(mylength-$window)){
   newMDslice <- newMD[i:(i+$window),1] 
   queryMDslice <- queryMD[i:(i+$window),1];
   queryMDrefslice <- queryMDref[i:(i+$window),1];
   mytest1 <- cancor(newMDslice, queryMDrefslice)
   mycor <- mytest1\$cor
   mytest2 <- cancor(queryMDslice, queryMDrefslice)
   myselfcor <- mytest2\$cor
   N = length(queryMDslice)
   p = 1
   q = 1
   mypvaltest <- p.asym(myselfcor, N, p, q, tstat = 'Wilks')
   mypval <- mypvaltest\$p.value[1]
   if(mypval <= $cancor_threshold){mybar = 1}
   if(mypval > $cancor_threshold){mybar = 0}
   if($v == 0){
   my_bars_$queryID <- c(my_bars_$queryID, mybar)
   }
   myimpact = abs(myselfcor*(log(mycor/myselfcor)))
   my_impact_sum_$variantID = my_impact_sum_$variantID + myimpact
   if(mybar == 1){my_impact_cons_sum_$variantID = my_impact_cons_sum_$variantID + myimpact}
   if(mybar == 0){my_impact_cons_sum_$variantID = my_impact_cons_sum_$variantID + 0}
   if ($v > 0){ 
     if (myimpact >= myupper){myimpact = abs(myselfcor*(log(mycor/myselfcor)))}
     if (myimpact < myupper){myimpact = 0}
   }
   my_impact_length_$variantID = my_impact_length_$variantID + 1
   if(mybar == 1){my_impact_cons_length_$variantID = my_impact_cons_length_$variantID + 1}
   
   my_impact_$variantID <- c(my_impact_$variantID, myimpact)
   my_ccpos_$variantID <- c(my_ccpos_$variantID, i)
   my_cancors_$variantID <- c(my_cancors_$variantID, mycor)
   my_pvals_$variantID <- c(my_pvals_$variantID, mypval)
   }\n";  # NOTE: can set impact to zero if less than upper bound
print Rinput "print(my_cancors_$variantID)\n";
print Rinput "print(my_pvals_$variantID)\n";
print Rinput "my_adjpvals_$variantID <- p.adjust(my_pvals_$variantID, method = 'fdr', n = length(my_pvals_$variantID))\n";
print Rinput "print(my_adjpvals_$variantID)\n";
print Rinput "print(my_ccpos_$variantID)\n";
print Rinput "print(my_impact_$variantID)\n";
print Rinput "custom.bootTTL <- function(times, data=my_impact_$variantID, length=my_impact_length_$variantID) { 
  boots <- rep(NA, times)
  for (i in 1:times) {
    boots[i] <- sd(sample(data, length, replace=TRUE))/sqrt(length)  
  }
  boots
}\n";
print Rinput "custom.bootCONS <- function(times, data=my_impact_$variantID, length=my_impact_cons_length_$variantID) { 
  boots <- rep(NA, times)
  for (i in 1:times) {
    boots[i] <- sd(sample(data, length, replace=TRUE))/sqrt(length)  
  }
  boots
}\n";
print Rinput "mySE <- sum(custom.bootTTL(times=1000))\n"; # bootstrap standard error of sum
print Rinput "print(mySE)\n";
print Rinput "myUpperSE <- (my_impact_sum_$variantID + mySE)\n";
print Rinput "myLowerSE <- (my_impact_sum_$variantID - mySE)\n";
print Rinput "myN = $varcount\n";
print Rinput "mySD = mySE\n";
print Rinput "print(mySD)\n";
print Rinput "print(my_impact_sum_$variantID)\n";
print Rinput "print(myUpperSE)\n";
print Rinput "print(myLowerSE)\n";
print Rinput "mySE_cons <- sum(custom.bootCONS(times=1000))\n"; # bootstrap standard error of sum
print Rinput "print(mySE_cons)\n";
print Rinput "myUpperSE_cons <- (my_impact_cons_sum_$variantID + mySE_cons)\n";
print Rinput "myLowerSE_cons <- (my_impact_cons_sum_$variantID - mySE_cons)\n";
print Rinput "mySD_cons = mySE_cons\n";
print Rinput "print(my_impact_cons_sum_$variantID)\n";
print Rinput "print(myUpperSE_cons)\n";
print Rinput "print(myLowerSE_cons)\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_IDs <- c('$variantID_label', my_impact_IDs)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sums <- c(my_impact_sum_$variantID, my_impact_sums)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sums_upperSE <- c(my_impact_sum_$variantID+mySE, my_impact_sums_upperSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sums_lowerSE <- c(my_impact_sum_$variantID-mySE, my_impact_sums_lowerSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sd <- c(mySD, my_impact_sd)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_n <- c(myN, my_impact_n)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sums <- c(my_impact_cons_sum_$variantID, my_impact_cons_sums)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sums_upperSE <- c(my_impact_cons_sum_$variantID+mySE_cons, my_impact_cons_sums_upperSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sums_lowerSE <- c(my_impact_cons_sum_$variantID-mySE_cons, my_impact_cons_sums_lowerSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sd <- c(mySD_cons, my_impact_cons_sd)}\n";
print Rinput "print(my_bars_$queryID)\n";
print Rinput "my_adjbars_$queryID <- c()\n";
print Rinput "print(my_adjpvals_$variantID\[i])\n";
print Rinput "for(j in 1:(length(my_adjpvals_$variantID))){
   print(my_adjpvals_$variantID\[j])
   if (my_adjpvals_$variantID\[j] <= $cancor_threshold){myadjbar = 1}
   if (my_adjpvals_$variantID\[j] > $cancor_threshold){myadjbar = 0}
   my_adjbars_$queryID <- c(my_adjbars_$queryID, myadjbar) 
   }\n";
print Rinput "print(my_adjbars_$queryID)\n";  
print Rinput "sink('./testingData_$queryID/cancor_stats_$variantID.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $variantID ($window AA sliding window)\n\n')\n";
print Rinput "print(my_cancors_$variantID)\n";
print Rinput "cat('\nMUTATIONAL IMPACT ANALYSIS OF $queryID and $variantID ($window AA sliding window)\n\n')\n";
print Rinput "print(my_impact_$variantID)\n";
print Rinput "sink()\n";
print Rinput "dataframe4_$variantID = data.frame(pos=my_ccpos_$variantID+$startN, cc=my_cancors_$variantID, bars=my_adjbars_$queryID)\n";
print Rinput "dataframe5_$variantID = data.frame(pos=my_ccpos_$variantID+$startN, cc=my_impact_$variantID)\n";
print Rinput "my_overall <- cancor(newMD, queryMD)\n"; # overall
print Rinput "my_ccor_$variantID <- my_overall\$cor[1]\n";
print Rinput "print(my_ccor_$variantID)\n";

# get bootstrap if first comparison
if($variantID eq $variants[0]){
    print Rinput "print(my_impact_$variantID)\n";
    print Rinput "mymean <- mean(my_impact_$variantID)\n";
    print Rinput "mysd <- sd(my_impact_$variantID)\n";
    print Rinput "myupper <- mymean + (2*mysd)\n";
    print Rinput "print(myupper)\n";
    }
# create output impact files for movie render
print Rinput "sink('./testingData_$queryID/impact_$variantID.txt')\n";
print Rinput "print(my_impact_$variantID)\n";
print Rinput "sink()\n";   
    
} # end for loop

# create plot 3
print Rinput "myplot3 <- ggplot() + ggtitle('signif local variant impact for $window AA window (+2sd of validation set)') + labs(x = 'position (residue number)', y = 'local index of impact') + theme(legend.title = element_blank())\n"; 
print Rinput "mylist3 <- list()\n";
# create lines for plot 3
for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     print Rinput "myline_$v <- geom_line(data = dataframe5_$variantID, mapping = aes(x = pos, y = cc, color = '$variantID_label'))\n"; 
     print Rinput "mylist3 <- list(mylist3, myline_$v)\n";
} # end for loop
print Rinput "myplot3final <- myplot3 + mylist3 + scale_color_brewer(palette='Set2')\n"; 


print Rinput "print(my_impact_IDs)\n";
print Rinput "print(my_impact_sums)\n";
print Rinput "print(my_impact_cons_sums)\n";
print Rinput "print(my_impact_sd)\n";
print Rinput "print(my_impact_cons_sd)\n";
print Rinput "print(my_impact_n)\n";

# create plot 4 
if ($bartype eq "total"){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_sums, upperSE = my_impact_sums_upperSE, lowerSE = my_impact_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, y = sum, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + scale_fill_brewer(palette = 'Set2') + ggtitle('total mutational impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "conserved"){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, upperSE = my_impact_cons_sums_upperSE, lowerSE = my_impact_cons_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + scale_fill_brewer(palette = 'Set2') + ggtitle('conserved region impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}

if ($bartype eq "total" && scalar(@variants) >= 4){ # ANOVA table from summary data
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_sums, sd = mySE, n = my_impact_n)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print(pval)\n";
print Rinput "sink('./maxDemon_results/ANOVAresult.txt')\n";
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_sums, sd = mySE, n = my_impact_n)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print('Fval')\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print('numerator DF')\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print('denominator DF')\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print('pval')\n";
print Rinput "print(pval)\n";
print Rinput "sink()\n"; 
}
if ($bartype eq "conserved" && scalar(@variants) >= 4){ # ANOVA table from summary data
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, sd = mySE_cons, n = my_impact_n)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print(pval)\n";
print Rinput "sink('./maxDemon_results/ANOVAresult.txt')\n";
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, sd = mySE_cons, n = my_impact_n)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print('Fval')\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print('numerator DF')\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print('denominator DF')\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print('pval')\n";
print Rinput "print(pval)\n";
print Rinput "sink()\n";
}
#####################################

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 2)))\n";
print Rinput "print(myplot2final, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))\n";
print Rinput "print(myplot3final, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))\n";
print Rinput "print(myplot1final, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
print Rinput "print(myplot4final, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))\n";
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
mkdir("maxDemon_results");
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$queryID/learningPLOTcompare.pdf";
copy($oldfilename, $newfilename);	
my $oldfilename2 = "Rplots.pdf";
my $newfilename2 = "./maxDemon_results/learningPLOTcompare.pdf";
copy($oldfilename2, $newfilename2);
print " machine learning is complete\n\n";
print " close PDF and txt viewer to continue\n\n";
system "evince ./testingData_$queryID/learningPLOTcompare.pdf\n";

# open can corr results
system "gedit ./testingData_$queryID/cancor_stats.txt\n";
     my $oldfilename3 = "./testingData_$queryID/cancor_stats.txt";
     my $newfilename3 = "./maxDemon_results/cancor_stats.txt";
     copy($oldfilename3, $newfilename3);
for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     system "gedit ./testingData_$queryID/cancor_stats_$variantID.txt\n";
     my $oldfilename4 = "./testingData_$queryID/cancor_stats_$variantID.txt";
     my $newfilename4 = "./maxDemon_results/cancor_stats_$variantID.txt";
     copy($oldfilename4, $newfilename4);
} # end for loop
if(scalar(@variants) >= 4){system "gedit ./maxDemon_results/ANOVAresult.txt\n";}
}
###########################################################################################################
###########################################################################################################





