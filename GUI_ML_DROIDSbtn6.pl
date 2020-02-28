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

$frameCount = $framenumber;

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
      #if ($header eq "color_scheme"){$colorType = $value;}
}
close IN;

$colorScheme = "c1";

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




########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################


#sub movie1 {



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


sleep(1);print "\nSELECT INTERACTION TYPE (1=protein only | 2=protein-protein | 3=DNA-protein | 4=protein-ligand)\n\n";
my $stype_number = <STDIN>;
chop($stype_number);

if($stype_number == 1){$stype = "protein";}
if($stype_number == 2){$stype = "protprot";}
if($stype_number == 3){$stype = "dna";}
if($stype_number == 4){$stype = "ligand";}

if($stype eq "dna" || $stype eq "ligand"){$orig_queryID = $queryID; $queryID = $refID;}


# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

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
 
     
#}

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




