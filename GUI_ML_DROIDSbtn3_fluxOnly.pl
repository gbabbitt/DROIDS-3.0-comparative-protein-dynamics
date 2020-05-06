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
# initialize ML methods
$method_bnp = 0;
$method_dist = 0;
$method_kern = 0;
$method_ens = 0;

open(IN3, "<"."MLmethods.txt") or die "could not find MLmethods.txt file\n";
my @IN3 = <IN3>;
for (my $i = 0; $i < scalar @IN3; $i++){
	 my $INrow = $IN3[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $method = @INrow[0];
	 if ($method eq "bnp"){$method_bnp = 1;}
      if ($method eq "dist"){$method_dist = 1;}
      if ($method eq "kern"){$method_kern = 1;}
      if ($method eq "ens"){$method_ens = 1;}
      if ($method eq "no_bnp"){$method_bnp = 0;}
      if ($method eq "no_dist"){$method_dist = 0;}
      if ($method eq "no_kern"){$method_kern = 0;}
      if ($method eq "no_ens"){$method_ens = 0;}
}
close IN3;

print ("\nMLmethods\n\n");
if ($method_bnp == 1){print"KNN activated\n";}
if ($method_bnp == 0){print"KNN deactivated\n";}
if ($method_dist == 1){print"NB/LDA/QDA activated\n";}
if ($method_dist == 0){print"NB/LDA/QDA deactivated\n";}
if ($method_kern == 1){print"SVM activated\n";}
if ($method_kern == 0){print"SVM deactivated\n";}
if ($method_ens == 1){print"random forest/adaboost activated\n";}
if ($method_ens == 0){print"random forest/adaboost deactivated\n";}
print ("\n\n\n");

###########################################################################################################
###########################################################################################################
system "x-terminal-emulator -e htop\n";
sleep(1);
print "\n\nINITIATING MACHINE LEARNING ON PROTEIN ATOM FLUCTUATION ONLY\n\n";
sleep(3);


if ($method_kern == 1 || $method_other == 1){
# prompt user - choose best learning model to display
sleep(1);print "CHOOSE KERNEL TYPE FOR SVM (type 'linear', 'polynomial', 'laplace' or 'RBF')\n\n";
my $kerntype_enter = <STDIN>;
chop($kerntype_enter);
if ($kerntype_enter eq ''){$kerntype_enter = 'linear';}

if($kerntype_enter eq 'linear'){$kerntype = 'vanilladot';}
if($kerntype_enter eq 'polynomial'){$kerntype = 'polydot';}
if($kerntype_enter eq 'laplace'){$kerntype = 'laplacedot';}
if($kerntype_enter eq 'RBF'){$kerntype = 'rbfdot';}
elsif($kerntype_enter eq ''){$kerntype = 'vanilladot';}
sleep(1);
}

sleep(1);print "MASK LEARNING ONLY TO SIGNIFICANTLY DIFFERENT DYNAMICS? (type 'on' or 'off' | default = 'off')\n\n";
my $option = <STDIN>;
chop($option);
if ($option eq ''){$option = 'off';}


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
if ($method_dist == 1){ 
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

} # end outer for loop
print "\n machine learning is complete\n\n";
exit;
###########################################################################################################
###########################################################################################################




