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


###########################################################################################################
###########################################################################################################
#####################################################
#sub canon{

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
#print Rinput "print(dataframeKNN_$v)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./testingData_$copyID/classpositionHISTOnb.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
print Rinput "dataframeNB_$v = data.frame(pos1=$AAposition_nb+$startN, Y1val=$sum_classifiers_nb)\n";
#print Rinput "print(dataframeNB_$v)\n";
# LDA
print Rinput "datatableLDA = read.table('./testingData_$copyID/classpositionHISTOlda.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
print Rinput "dataframeLDA_$v = data.frame(pos1=$AAposition_lda+$startN, Y1val=$sum_classifiers_lda)\n";
#print Rinput "print(dataframeLDA_$v)\n";
# QDA
print Rinput "datatableQDA = read.table('./testingData_$copyID/classpositionHISTOqda.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
print Rinput "dataframeQDA_$v = data.frame(pos1=$AAposition_qda+$startN, Y1val=$sum_classifiers_qda)\n";
#print Rinput "print(dataframeQDA_$v)\n";
# SVM
print Rinput "datatableSVM = read.table('./testingData_$copyID/classpositionHISTOsvm.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
print Rinput "dataframeSVM_$v = data.frame(pos1=$AAposition_svm+$startN, Y1val=$sum_classifiers_svm)\n";
#print Rinput "print(dataframeSVM_$v)\n";
# random forest
print Rinput "datatableRFOR = read.table('./testingData_$copyID/classpositionHISTOrfor.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
print Rinput "dataframeRFOR_$v = data.frame(pos1=$AAposition_rfor+$startN, Y1val=$sum_classifiers_rfor)\n";
#print Rinput "print(dataframeRFOR_$v)\n";
# adaboost
print Rinput "datatableADA = read.table('./testingData_$copyID/classpositionHISTOada.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
print Rinput "dataframeADA_$v = data.frame(pos1=$AAposition_ada+$startN, Y1val=$sum_classifiers_ada)\n";
#print Rinput "print(dataframeADA_$v)\n";

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

if (scalar(@variants) <= 8){print Rinput "myplot2final <- myplot2 + mylist2 + scale_color_brewer(palette='Set2')\n";} 
if (scalar(@variants) > 8){print Rinput "myplot2final <- myplot2 + mylist2\n";} 

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
if(scalar(@variants) <= 8){print Rinput "myplot3final <- myplot3 + mylist3 + scale_color_brewer(palette='Set2')\n";} 
if(scalar(@variants) > 8){print Rinput "myplot3final <- myplot3 + mylist3\n";} 

print Rinput "print(my_impact_IDs)\n";
print Rinput "print(my_impact_sums)\n";
print Rinput "print(my_impact_cons_sums)\n";
print Rinput "print(my_impact_sd)\n";
print Rinput "print(my_impact_cons_sd)\n";
print Rinput "print(my_impact_n)\n";

# create plot 4 
if ($bartype eq "total" && scalar(@variants) <= 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_sums, upperSE = my_impact_sums_upperSE, lowerSE = my_impact_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, y = sum, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + scale_fill_brewer(palette = 'Set2') + ggtitle('total mutational impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "conserved" && scalar(@variants) <= 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, upperSE = my_impact_cons_sums_upperSE, lowerSE = my_impact_cons_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + scale_fill_brewer(palette = 'Set2') + ggtitle('conserved region impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "total" && scalar(@variants) > 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_sums, upperSE = my_impact_sums_upperSE, lowerSE = my_impact_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, y = sum, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + ggtitle('total mutational impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "conserved" && scalar(@variants) > 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, upperSE = my_impact_cons_sums_upperSE, lowerSE = my_impact_cons_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + ggtitle('conserved region impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
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
print Rinput "print(dataframe6)\n";
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
print Rinput "print(dataframe6)\n";
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
#}



###########################################
###########################################

# prompt user - selest machine learner for MI matrix
sleep(1);print "CHOOSE BEST LEARNER FOR MUTUAL INFORMATION MATRIX (type 'KNN', 'NB', 'LDA', 'QDA', 'SVM', 'RFOR' or 'ADA')\n\n";
my $learner = <STDIN>;
chop($learner);

for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];

# initialize     
$residue1 = 0;
$residue2 = 0;
$MI = 0;
# open MI matrix files
open(MI, ">"."./testingData_$variantID/MIcolumn_$learner"."_$variantID.txt") or die "could not open MI matrix file for $variantID \n";    
print MI "pos1\t"."pos2\t"."MI\n";
open(MI2, ">"."./testingData_$variantID/MImatrix_$learner"."_$variantID.txt") or die "could not open MI matrix file for $variantID \n";    
#print MI2 "residue\t";
#for (my $n = 0; $n < $lengthID; $n++){print MI2 "$n\t";}
#print MI2 "\n";
 
   
print "\ncalculating mutual information matrix for $variantID \n";     
sleep(1);

for (my $i = 0; $i < $lengthID; $i++){
	 $residue1 = $i;
      print "calculating MI values for residue\t"."$residue1 for $variantID\n";
      open(POS1, "<"."./testingData_$variantID/indAAclass$learner/classAA_$refID"."_$i.txt") or die "could not open MI matrix file for ./testingData_$variantID/indAAclass$learner/classAA_$refID"."_$i.txt \n";    
      my @POS1 = <POS1>;
      if($i>0){print MI2 "\n";}
      #if($i==0){print MI2 "$residue1\t";}
      #if($i>0){print MI2 "\n"; print MI2 "$residue1\t";}
      for (my $j = 0; $j < $lengthID; $j++){
	     open(POS2, "<"."./testingData_$variantID/indAAclass$learner/classAA_$refID"."_$j.txt") or die "could not open MI matrix file for ./testingData_$variantID/indAAclass$learner/classAA_$refID"."_$j.txt \n";    
          my @POS2 = <POS2>;
          $residue2 = $j;
          
          #calculate freq A
          #print "calculating freq state A\t";
          $freqA = 0;
          $sumclassA = 0;
          $cntA = 0;
          for (my $ii = 0; $ii < scalar @POS1; $ii++){
              my $POS1row = $POS1[$ii];
              my @POS1row = split (/\s+/, $POS1row);
	         $class1 = $POS1row[0];
              $time1 = $ii;
              if($class1 == 0 || $class1 == 0.5 || $class1 == 1) {$sumclassA = $sumclassA+$class1; $cntA = $cntA+1;}
              #print "Residue1\t"."$residue1\t"."timeslice\t"."$time1\t"."class = "."$class1\n";
              }
              $freqA = $sumclassA/($cntA+0.0001);
              #print "freqA = "."$freqA\n";
          
          #calculate freq B
          #print "calculating freq state B\t";
          $freqB = 0;
          $sumclassB = 0;
          $cntB = 0;
          for (my $jj = 0; $jj < scalar @POS2; $jj++){
              my $POS2row = $POS2[$jj];
              my @POS2row = split (/\s+/, $POS2row);
	         $class2 = $POS2row[0];
              $time2 = $jj;
              if($class2 == 0 || $class2 == 0.5 || $class2 == 1) {$sumclassB = $sumclassB+$class2; $cntB = $cntB+1;}
              #print "Residue2\t"."$residue2\t"."timeslice\t"."$time2\t"."class = "."$class2\n";
              }
              $freqB = $sumclassB/($cntB+0.0001);
              #print "freqB = "."$freqB\n";
          
          #calculate freq A and B
          #print "calculating freq state A and B\t";
          $freqAB = 0;
          $sumclassAB = 0;
          $cntAB = 0;
          for (my $ii = 0; $ii < scalar @POS1; $ii++){
              my $POS1row = $POS1[$ii];
              my @POS1row = split (/\s+/, $POS1row);
	         $class1 = $POS1row[0];
              $time1 = $ii;
              #if($class1 == 0 || $class1 == 0.5 || $class1 == 1) {$sumclassA = $sumclassA+$class1; $cntA = $cntA+1;}
              #print "Residue1\t"."$residue1\t"."timeslice\t"."$time1\t"."class = "."$class1\n";
              for (my $jj = 0; $jj < scalar @POS2; $jj++){
              my $POS2row = $POS2[$jj];
              my @POS2row = split (/\s+/, $POS2row);
	         $class2 = $POS2row[0];
              $time2 = $jj;
              if($time1 ==$time2 && $class1 == $class2) {$sumclassAB = $sumclassAB+1; $cntAB = $cntAB+1;}
              if($time1 ==$time2 && $class1 != $class2) {$cntAB = $cntAB+1;}
              #print "Residue2\t"."$residue2\t"."timeslice\t"."$time2\t"."class = "."$class2\n";
              }
              }
              $freqAB = $sumclassAB/($cntAB+0.0001);
              #print "freqAB = "."$freqAB\n";
                    
          # calculate MI
          $MI = $freqAB*log(($freqAB+0.0001)/(($freqA*$freqB)+0.0001));
          if ($MI>1){$MI = 1;}
          #print "MImatrix\t"."$residue1\t"."$residue2\t"."$MI\n";
          print MI "$residue1\t"."$residue2\t"."$MI\n";
          print MI2 "$MI\t";
      }
}

close MI;
close MI2;

print "\nheatmapping mutual information matrix for $variantID (close .pdf to continue)\n";     
sleep(1);

open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "datamatrixMI = read.table('./testingData_$variantID/MImatrix_$learner"."_$variantID.txt', header = FALSE)\n";
#print Rinput "print(datamatrixMI)\n";
print Rinput "datamatrixMI<-as.matrix(datamatrixMI)\n";
print Rinput "myMImean <- mean(datamatrixMI)\n";
print Rinput "myMImean <- round(myMImean, digits=4)\n";
print Rinput "print(myMImean)\n";
print Rinput "myMIsd <- sd(datamatrixMI)\n";
print Rinput "myMIsd <- round(myMIsd, digits=4)\n";
print Rinput "print(myMIsd)\n";
#print Rinput "print(datamatrixMI)\n";
#print Rinput "datamatrixMI <- scale(datamatrixMI)\n";
#print Rinput "mymap1<-heatmap(datamatrixMI, Colv = 'Rowv', symm = TRUE, keep.dendro = FALSE)\n";
#print Rinput "print(mymap1)\n";
print Rinput "x <- (1:nrow(datamatrixMI))\n";
print Rinput "y <- (1:ncol(datamatrixMI))\n";
print Rinput "mymap2<-image(x+$startN, y+$startN, datamatrixMI, col = gray.colors(20), main = c('MUTUAL INFORMATION for $variantID (lighter = higher))', paste('mean MI = ', myMImean, 'sd MI = ', myMIsd)), xlab = 'residue position', ylab = 'residue position')\n";
print Rinput "print(mymap2)\n";      
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$variantID/MImatrix_$learner"."_$variantID.pdf";
my $newfilename2 = "./maxDemon_results/MImatrix_$learner"."_$variantID.pdf";
copy($oldfilename, $newfilename);
copy($oldfilename, $newfilename2);
my $oldfilename3 = "./testingData_$variantID/MImatrix_$learner"."_$variantID.txt";
my $newfilename3 = "./maxDemon_results/MImatrix_$learner"."_$variantID.txt";
copy($oldfilename3, $newfilename3);
my $oldfilename4 = "./testingData_$variantID/MIcolumn_$learner"."_$variantID.txt";
my $newfilename4 = "./maxDemon_results/MIcolumn_$learner"."_$variantID.txt";
copy($oldfilename4, $newfilename4);
print " mutual information matrix is complete\n\n";
print " close PDF and txt viewer to continue\n\n";


} # end for loop

for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     system "evince ./testingData_$variantID/MImatrix_$learner"."_$variantID.pdf\n";
}

###########################################################################################################
###########################################################################################################





