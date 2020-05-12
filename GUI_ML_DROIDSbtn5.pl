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

#sub image1 {

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
$orig_queryID = $queryID;  # create tag for original query in training set
$queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

print "\nEnter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot'){
   print "\nEnter number position of N terminal on this chain (default = 0)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }


print "\nDO YOU WANT TO MAP CONSERVED REGIONS OR MUTUAL INFORMATION VALUES?\n";
   print "(TYPE = 'cons' or 'mi')\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure\n";
   $typeMAP = <STDIN>;
   chop($typeMAP);
   if ($typeMAP eq ''){$typeMAP = 'cons';}
   if ($typeMAP eq 'CONS'){$typeMAP = 'cons';}
   if ($typeMAP eq 'MI'){$typeMAP = 'mi';}

if ($typeMAP eq 'mi'){   
# prompt user - selest machine learner for MI matrix
sleep(1);print "\nCHOOSE BEST LEARNER FOR MUTUAL INFORMATION MAPPING (type 'KNN', 'NB', 'LDA', 'QDA', 'SVM', 'RFOR' or 'ADA')\n\n";
$learner = <STDIN>;
chop($learner);
# prompt user - selest variant for MI matrix
sleep(1);print "\nSELECT VARIANT ID FOR MUTUAL INFORMATION MAPPING (e.g. type '1yet_unbound_2')\n\n";
$variantID = <STDIN>;
chop($variantID);
# prompt user - selest AA site for MI mapping
sleep(1);print "\nSELECT AMINO ACID REFERENCE SITE MUTUAL INFORMATION MAPPING (e.g. type '50' for interaction strengths with AA 50)\n\n";
$siteID = <STDIN>;
chop($siteID);
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

if ($typeMAP eq 'cons'){

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
}

if ($typeMAP eq 'mi'){

print "collecting MI values over all interaction sites for given amino acid\n";
   sleep(1);
   open(OUT, ">"."./maxDemon_results/MIsite_$learner"."_$variantID.txt")||die "could not open mask file /maxDemon_results/MIcolumn_$learner"."_$variantID.txt\n";  
   print OUT "position\t"."MIvalue\n";
   open(IN, "<"."./maxDemon_results/MIcolumn_$learner"."_$variantID.txt")||die "could not open mask file /maxDemon_results/MIcolumn_$learner"."_$variantID.txt\n";  
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++){
      $INrow = $IN[$a];
      @INrow = split(/\s+/, $INrow);
      $first = $INrow[0];
      $second = $INrow[1];
      $MIvalue = $INrow[2];
      if ($siteID == $first){print OUT "$second\t"."$MIvalue\n";}
    
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
   open(IN, "<"."./maxDemon_results/MIsite_$learner"."_$variantID.txt")||die "could not open mask file /maxDemon_results/MIcolumn_$learner"."_$variantID.txt\n";  
    my @IN = <IN>;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $MIvalue = $INrow[1];
   
     #print "$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     print OUT "\t:"."$pos".".$chainMAP\t"."$MIvalue\n";
     
 
 close IN;
 }
 close OUT;

print "\n class attribute files for all frames are created\n";
sleep(1);
}



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
     
#}

