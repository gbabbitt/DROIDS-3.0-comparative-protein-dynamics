#!/usr/bin/perl
use Tk;
#use strict;
#use warnings;
use feature ":5.10";
use Statistics::Descriptive();
use File::Copy;

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

# read control files
open(IN, "<"."DROIDS.ctl") or die "could not find DROIDS.txt file\n";
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
      if ($header eq "homology"){$homology = $value;}
      if ($header eq "num_chains"){$chainN = $value;}
      if ($header eq "shape"){$vector_enter = $value;}
}

close IN;

# read control files
open(IN, "<"."MDq.ctl") or die "could not find MDq.ctl file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "Force_Field"){$forceID = $value;}
      if ($header eq "DNA_Field" || $header eq "LIGAND_Field"){$dforceID = $value;}
      if ($header eq "PDB_ID"){$fileIDq = $value; $fileIDq = substr($fileIDq, 0, length ($fileIDq)-7);} # trims 'REDUCED' off end of ref label
      if ($header eq "LIGAND_ID"){$fileIDl = $value;}
      if ($header eq "Heating_Time"){$cutoffValueHeatFS = $value;}
      if ($header eq "Equilibration_Time"){$cutoffValueEqFS = $value;}
      if ($header eq "Production_Time"){$cutoffValueProdFS = $value;}
      if ($header eq "Temperature_Query"){$tempQ = $value;}
      if ($header eq "Solvation_Method"){$solvType = $value;}
      if ($header eq "Salt_Conc"){$cutoffValueSalt = $value;}
      
}
close IN;

if ($tempQ eq '#'){$tempQ = 300;}

open(IN, "<"."MDr.ctl") or die "could not find MDr.ctl file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "PDB_ID"){$fileIDr = $value; $fileIDr = substr($fileIDr, 0, length ($fileIDr)-7);}  # trims 'REDUCED' off end of ref label
      if ($header eq "Number_Runs"){$number_runs = $value;}
}
close IN;

#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### Declare variables ####
#my $fileIDq = '';
#my $fileIDr = '';
#my $forceID = '';
my $runsID = '';
my $implicit=0;
my $explicit=0;
my $solvType = '';
#my $cutoffValueHeat=100;
#my $cutoffValueEq=10;
#my $cutoffValueProd=10;
#my $cutoffValueSalt=0.0;
#my $cutoffValueHeatFS=0;
#my $cutoffValueEqFS=0;
#my $cutoffValueProdFS=0;
$cutoffValueHeat = $cutoffValueHeatFS/1000;
$cutoffValueEq = $cutoffValueEqFS/1000000;
$cutoffValueProd = $cutoffValueProdFS/1000;

$previousProdTime = $cutoffValueProdFS;

my @fullfile;
my @chainlen;
my @fullfile2;
my @chainlen2;


#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("MD control settings for ML deployment"); # Titles the main window
$mw->setPalette("gray");

my $MDheatScale = $mw->Scale(-label=>"Length of MD heating run (ps) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>1000,
			-variable=>\$cutoffValueHeat,
			-tickinterval=>200,
			-resolution=>10,
			-length=>205
			);

my $MDeqScale = $mw->Scale(-label=>"Length of MD equilibration run (ns) :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>20,
			-variable=>\$cutoffValueEq,
			-tickinterval=>5,
			-resolution=>1,
			-length=>205
			);

#my $MDprodScale = $mw->Scale(-label=>"Length of MD deployment run (ns) :",
#			-orient=>'h',
#			-digit=>3,
#			-from=>0,
#			-to=>100,
#			-variable=>\$cutoffValueProd,
#			-tickinterval=>20,
#			-resolution=>1,
#			-length=>205
#			);

my $MDsaltScale = $mw->Scale(-label=>"extra salt conc (M) (implicit only)  :",
			-orient=>'h',
			-digit=>3,
			-from=>0,
			-to=>0.6,
			-variable=>\$cutoffValueSalt,
			-tickinterval=>0.2,
			-resolution=>0.05,
			-length=>205
			);

# Solvation Frame
my $solnFrame = $mw->Frame(	-label => "METHOD OF SOLVATION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $implicitCheck = $solnFrame->Radiobutton( -text => "implicit - Generalized Born",
						-value=>"im",
						-variable=>\$solvType
						);
	my $explicitCheck = $solnFrame->Radiobutton( -text => "explicit - Particle Mesh Ewald",
						-value=>"ex",
						-variable=>\$solvType
						);

# Simulation Frame
my $simFrame = $mw->Frame(	-label => "MD SIMULATION ENGINE",
				-relief => "groove",
				-borderwidth => 2
				);
	my $amberCheck = $simFrame->Radiobutton( -text => "pmemd.cuda (amber) - licensed",
						-value=>"amber",
						-variable=>\$simType
						);
	my $openCheck = $simFrame->Radiobutton( -text => "OpenMM - open source",
						-value=>"open",
						-variable=>\$simType
						);

# MD production length Frame
my $mdprodFrame = $mw->Frame(	-label => "LENGTH OF PRODUCTION",
				-relief => "groove",
				-borderwidth => 2
				);
	my $oneCheck = $mdprodFrame->Radiobutton( -text => "1x training run",
						-value=>1,
						-variable=>\$prodLen
						);
	my $twoCheck = $mdprodFrame->Radiobutton( -text => "2x training run",
						-value=>2,
						-variable=>\$prodLen
						);
     my $threeCheck = $mdprodFrame->Radiobutton( -text => "3x training run",
						-value=>3,
						-variable=>\$prodLen
						);
	my $fiveCheck = $mdprodFrame->Radiobutton( -text => "5x training run",
						-value=>5,
						-variable=>\$prodLen
						);
     my $tenCheck = $mdprodFrame->Radiobutton( -text => "10x training run",
						-value=>10,
						-variable=>\$prodLen
						);

# PDB ID Frame				
my $pdbFrame = $mw->Frame();
#	my $QfileFrame = $pdbFrame->Frame();
#		my $QfileLabel = $QfileFrame->Label(-text=>"pdb ID query (for deploying learner) : ");
#		my $QfileEntry = $QfileFrame->Entry(-borderwidth => 2,
#					-relief => "groove",
#					-textvariable=>\$fileIDq
#					);
	my $forceFrame = $pdbFrame->Frame();
		my $forceLabel = $forceFrame->Label(-text=>"Force Field (from previous runs): ");
		my $forceEntry = $forceFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$forceID
					);
	my $dforceFrame = $pdbFrame->Frame();
		my $dforceLabel = $dforceFrame->Label(-text=>"additional DNA force field: ");
		my $dforceEntry = $dforceFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$dforceID
					);
    my $prodFrame = $pdbFrame->Frame();
		my $prodLabel = $prodFrame->Label(-text=>"length of training MD run (DO NOT CHANGE): ");
		my $prodEntry = $prodFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$cutoffValueProd
					);
     my $tempFrame = $pdbFrame->Frame();
		my $tempLabel = $tempFrame->Label(-text=>"temperature K on query PDB (DO NOT CHANGE): ");
		my $tempEntry = $tempFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$tempQ
					);
     my $runsFrame = $pdbFrame->Frame();
		my $runsLabel = $runsFrame->Label(-text=>"number of repeated MD sample runs: ");
		my $runsEntry = $runsFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$runsID
					);
          $runsID = 1;  # this is hard coded for deployment on single run
          
      my $chainFrame = $pdbFrame->Frame();
		my $chainLabel = $chainFrame->Label(-text=>"number of protein chains (e.g. 3 = A/B/C): ");
		my $chainEntry = $chainFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$chainN
					);
      #  create list of chain labels
      @alphabet = ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");
      @chainlist = ();
      if($chainN > 26) {print "warning - number of chains exceeds alphabet\n";}
      for(my $l = 0; $l < $chainN; $l++){
           $letter = $alphabet[$l];
           push(@chainlist, $letter);
           }
      print "chains in structure are...\n";
      print @chainlist;
      print "\n\n";
    
          
     my $startFrame = $pdbFrame->Frame();
		my $startLabel = $startFrame->Label(-text=>"start numbering AA's on chain at (e.g. 1): ");
		my $startEntry = $startFrame->Entry(-borderwidth => 2,
					-relief => "groove",
					-textvariable=>\$startN
					);
          $startN = 1; # this is now hard coded - Nov 2018 
          
# Buttons
my $varlistButton = $mw -> Button(-text => "create list of .pdb files for variants to analyze", 
				-command => \&varlist
				); # Creates a txt file button
my $controlButton = $mw -> Button(-text => "make control files for deployment run (.ctl)", 
				-command => \&control
				); # Creates a ctl file button
my $teLeapButton = $mw -> Button(-text => "generate topology and coordinate files (teLeap)", 
				-command => \&teLeap
				); # Creates a teLeap button
my $reduceButton = $mw -> Button(-text => "dry and reduce structure (run pdb4amber)", 
				-command => \&reduce
				); # Creates a pdb4amber button
my $alignButton = $mw -> Button(-text => "create sequence and structural alignment (UCSF Chimera)", 
				-command => \&align
				); # Creates a align button
my $launchButton = $mw -> Button(-text => "launch MD run", 
				-command => \&launch,
				-background => 'gray45',
				-foreground => 'white'
				); # Creates a launch button
my $infoButton = $mw -> Button(-text => "create atom info files", 
				-command => \&info
				); # Creates a file button

my $fluxButton = $mw -> Button(-text => "create atom fluctuation files", 
				-command => \&flux
				); # Creates a file button

my $doneButton = $mw -> Button(-text => "parse / prepare files atom fluctuation files", 
				-command => \&done
				); # Creates a file button
my $stopButton = $mw -> Button(-text => "exit back to machine learning", 
				-command => \&stop
				); # Creates a file button


#### Organize GUI Layout ####
$stopButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$doneButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$fluxButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$infoButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$launchButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$teLeapButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$alignButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$reduceButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$controlButton->pack(-side=>"bottom",
			-anchor=>"s"
			);
$varlistButton->pack(-side=>"bottom",
			-anchor=>"s"
			);

#$QfileLabel->pack(-side=>"left");
#$QfileEntry->pack(-side=>"left");
$forceLabel->pack(-side=>"left");
$forceEntry->pack(-side=>"left");
$dforceLabel->pack(-side=>"left");
$dforceEntry->pack(-side=>"left");
$prodLabel->pack(-side=>"left");
$prodEntry->pack(-side=>"left");
#$runsLabel->pack(-side=>"left");
#$runsEntry->pack(-side=>"left");
$tempLabel->pack(-side=>"left");
$tempEntry->pack(-side=>"left");
$chainLabel->pack(-side=>"left");
$chainEntry->pack(-side=>"left");
#$startLabel->pack(-side=>"left");
#$startEntry->pack(-side=>"left");

$forceFrame->pack(-side=>"top",
		-anchor=>"e");
$dforceFrame->pack(-side=>"top",
		-anchor=>"e");
#$QfileFrame->pack(-side=>"top",
#		-anchor=>"e");
$prodFrame->pack(-side=>"top",
		-anchor=>"e");
$tempFrame->pack(-side=>"top",
		-anchor=>"e");
#$runsFrame->pack(-side=>"top",
#		-anchor=>"e");
$chainFrame->pack(-side=>"top",
		-anchor=>"e");
#$startFrame->pack(-side=>"top",
#		-anchor=>"e");
$pdbFrame->pack(-side=>"top",
		-anchor=>"n");

$implicitCheck->pack();
$explicitCheck->pack();
$solnFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$amberCheck->pack();
$openCheck->pack();
$simFrame->pack(-side=>"top",
		-anchor=>"n"
		);
$MDheatScale->pack(-side=>"top");
$MDeqScale->pack(-side=>"top");
#$MDprodScale->pack(-side=>"top");
$MDsaltScale->pack(-side=>"top");

$mdprodFrame->pack(-side=>"top",
		-anchor=>"n");
$oneCheck->pack();
$twoCheck->pack();
$threeCheck->pack();
$fiveCheck->pack();
$tenCheck->pack();

MainLoop; # Allows Window to Pop Up


########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################
sub stop {exit;}

#####################################################################################################
sub varlist{
sleep(1);print "CHOOSE NUMBER OF REPLICATE RUNS ON IDENTICAL PDB (default = 2)\n\n";
my $num_copy = <STDIN>;
chop($num_copy);
if ($num_copy eq ''){$num_copy = 2;}

#NOTE: for DNA bound structures in DROIDS, the DNA bound form is the reference and so the reference homologs are simulated here

# make copies of query protein file for deploy 
for (my $c = 0; $c <= $num_copy ; $c++){
if ($c == 0){next;}
my $oldfilename = "$fileIDr".".pdb";
my $newfilename = "$fileIDr"."_$c.pdb";
copy($oldfilename, $newfilename);
}
open(MUT, ">"."copy_list.txt");
print MUT "PDB_IDs\n";
for (my $c = 0; $c <= $num_copy ; $c++){if ($c == 0){next;} print MUT "$fileIDr"."_$c\n";}
close MUT;
print "check PDB ID's for copies of $fileIDr then save and close\n\n";
system "gedit copy_list.txt\n";

# create mutate_protein.cmd script
open(MUT, ">"."variant_list.txt");
print MUT "PDB_IDs\n";
print MUT "$fileIDr"."_1\n";
print MUT "$fileIDr"."_2\n";
close MUT;
print "opening variant_list.txt using gedit\n\n";
print "type PDB ID's for additional target protein variants under 'PDB_IDs' then save and close\n\n";
print "for example\n\n";
print "PDB_IDs\n";
print "1cdw_bound_validation1\n";
print "1cdw_bound_validation2\n";
print "1cdw_bound_F34G\n";
print "1cdw_bound_H57D\n";
print "etc...\n\n\n";
print "LEAVE AS IS IF YOU ARE NOT ANALYZING ANY GENETIC OR DNA BINDING VARIANTS\n\n";
system "gedit variant_list.txt\n";

open(MUT, ">"."variant_label_list.txt");
print MUT "names\n";
print MUT "validation_run1\n";
print MUT "validation_run2\n";
close MUT;
print "opening variant_label_list.txt using gedit\n\n";
print "type names for additional variants as you want them to appear in plots then save and close\n\n";
print "for example\n\n";
print "names\n";
print "validation_run1\n";
print "validation_run2\n";
print "variant1\n";
print "variant2\n";
print "etc...\n\n\n";
print "LEAVE AS IS IF YOU ARE NOT ANALYZING ANY GENETIC OR DRUG BINDING VARIANTS\n\n";
system "gedit variant_label_list.txt\n";

print "\ncopy_list.txt, variant_label_list.txt, and variant_list.txt files were created\n";

}

########################################################################################
sub control { # Write a control file and then call appropriate scripts that reference control file

# append MDframes.ctl
print "appending MDframes.ctl\n";
sleep(1);
open (OUT, ">>"."MDframes.ctl") || die "could not open MDframes.ctl file\n";
print OUT "framegrpfactor\t"."$prodLen\n";
close OUT;

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
         
  ################  
     
     
     if ($solvType eq "im") {$repStr = "implicit";}
	if ($solvType eq "ex") {$repStr = "explicit";}
	
	# convert all times to femtosec
	$cutoffValueHeatFS = $cutoffValueHeat*1000;
	$cutoffValueEqFS = $cutoffValueEq*1000000;
	$cutoffValueProdFS = $cutoffValueProd*1000;

$currentProdTime = $cutoffValueProdFS*$prodLen;

   
### make query protein control file ###

### get chain information from PDB
## read in .pdb file
open(INPUTFILE, $fileIDq.".pdb");
    # load input into array
    #print"success";
    chomp(@fullfile = <INPUTFILE>);
    close(INPUTFILE);

my $count = 0;
for(my $line = 0; $line < scalar @fullfile; $line++){
    ## go down first column until you hit TER
    chomp($fullfile[$line]);
    my @entry = (split (/\s+/, $fullfile[$line]));
    if ($entry[0] eq "TER") {
        # get each chain length
        my $len = $entry[4];
        if ($len eq ''){$len = $entry[3]; $len =~ s/\D//g;}  #fixes concatenation of chain ID and residue number when > 1000
        $chainlen[$count] = $len;
        #print "$chainlen[$count]\n";
        $count++;
    }
}
### write control file
open(my $ctlFile1, '>', "MDq_deploy_$fileIDq.ctl") or die "Could not open output file";
print $ctlFile1 
"PDB_ID\t".$fileIDq."REDUCED\t# Protein Data Bank ID for MD run
Number_Chains\t$chainN\t# Number of chains on structure\n";
$allchainlen = 0;
for(my $ent = 0; $ent < scalar @chainlen; $ent++){
    my $chain = chr($ent + 65);
    print $ctlFile1 "length$chain\t$chainlen[$ent]\t #end of chain designated\n";
    print "MDq.ctl\n";
    print "length$chain\t$chainlen[$ent]\t #end of chain designated\n";
    $allchainlen = $allchainlen + $chainlen[$ent];
}

# define vector reference point (...as mid sequence in Chain A)
#$vectref = int(0.5*$allchainlen);
if ($vector_enter eq 'y'){
sleep(1);print "\nCHOOSE AN AMINO ACID RESIDUE AS REFERENCE POINT FOR VECTOR (i.e. shape) ANALYSIS (default = 1)\n\n";
my $vectref = <STDIN>;
chop($vectref);
if ($vectref eq ''){$vectref = 1;}
}

print $ctlFile1
"LIGAND_ID\t".$fileIDl."REDUCED\t# Protein Data Bank ID for MD run
Force_Field\t$forceID\t# AMBER force field to use in MD runs
DNA_Field\t$dforceID\t# AMBER force field to use in MD runs
Number_Runs\t$runsID\t# number of repeated samples of MD runs
Heating_Time\t$cutoffValueHeatFS\t# length of heating run (fs)
Equilibration_Time\t$cutoffValueEqFS\t# length of equilibration run (fs)
Production_Time\t$currentProdTime\t# length of production run (fs)
Solvation_Method\t$repStr\t# method of solvation (implicit or explicit)
Salt_Conc\t$cutoffValueSalt\t# salt concentration (implicit only, PME=O)
Temperature_Query\t$tempQ\t# temperature of query run (300K is same as ref run)";
close $ctlFile1;


print "MD control file is made (see MDq_deploy.ctl)\n";
################################
### make cpptraj ctl files
sleep(0.5);

if ($solvType eq "im") {$implicit = 1;}
if ($solvType eq "ex") {$explicit = 1;}

### make atom info control files ###	
open(ctlFile1, '>', "atominfo_$fileIDq"."_deploy.ctl") or die "Could not open output file";
my $parm_label1 = '';
if ($implicit == 1) {my $parm_label1 = "vac_"."$fileIDq"."REDUCED.prmtop"; print ctlFile1 "parm $parm_label1\n";}
if ($explicit == 1) {my $parm_label1 = "wat_"."$fileIDq"."REDUCED.prmtop"; print ctlFile1 "parm $parm_label1\n";}
my $traj_label1 = "prod_"."$fileIDq"."REDUCED_deploy".".nc";
print ctlFile1 "trajin $traj_label1\n";
print ctlFile1 "resinfo !(:WAT)\n"; # all residues but not water
print ctlFile1 "atominfo \@CA,C,O,N,H&!(:WAT)\n"; # mask for all protein backbone atoms eliminating water
close ctlFile1;

#for (my $i = 0; $i < $runsID; $i++){
### make atom flux control files ###	
open(ctlFile3, '>', "atomflux_$fileIDq"."_deploy.ctl") or die "Could not open output file";
my $parm_label3 = '';
if ($implicit == 1) {my $parm_label3 = "vac_"."$fileIDq"."REDUCED.prmtop"; print ctlFile3 "parm $parm_label3\n";}
if ($explicit == 1) {my $parm_label3 = "wat_"."$fileIDq"."REDUCED.prmtop"; print ctlFile3 "parm $parm_label3\n";}
my $traj_label3 = "prod_"."$fileIDq"."REDUCED_deploy".".nc";
print ctlFile3 "trajin $traj_label3\n";	
print ctlFile3 "rms first\n";
print ctlFile3 "average crdset MyAvg\n";
print ctlFile3 "run\n";
print ctlFile3 "rms ref MyAvg\n";
print ctlFile3 "atomicfluct out fluct_$fileIDq"."_deploy.txt \@CA,C,O,N,H&!(:WAT)\n";
#print ctlFile3 "byatom\n"; # hash out for avg atom flux, unhash for total atom flux
print ctlFile3 "run\n";
close ctlFile3;

if($vector_enter eq 'y'){
mkdir("atomvctl_$fileIDq"."_deploy");
mkdir("atomvect_$fileIDq"."_deploy");

for(my $j = 0; $j<$allchainlen; $j++){
open(ctlFile5, '>', "./atomvctl_$fileIDq"."_deploy/atomvector_$fileIDq"."_aa$j.ctl") or die "Could not open output file";
my $parm_label5 = '';
$jplus = $j+1;
if ($implicit == 1) {my $parm_label5 = "vac_"."$fileIDq"."REDUCED.prmtop"; print ctlFile5 "parm $parm_label5\n";}
if ($explicit == 1) {my $parm_label5 = "wat_"."$fileIDq"."REDUCED.prmtop"; print ctlFile5 "parm $parm_label5\n";}
my $traj_label5 = "prod_"."$fileIDq"."REDUCED_deploy".".nc";
print ctlFile5 "trajin $traj_label5\n";	
print ctlFile5 "vector out ./atomvect_$fileIDq"."_deploy/vect_$fileIDq"."_aa$j.txt :$jplus :$vectref\n";
print ctlFile5 "run\n";
close ctlFile5;
}
}
#  } # end per run loop 

my $prefix = "";
open(metafile, '>', "$fileIDr.meta") or die "Could not open output file";
if ($implicit == 1) {$prefix = "vac";}
if ($explicit == 1) {$prefix = "wat";}
print metafile "amber\n$prefix"."_$fileIDr.prmtop\nprod_$fileIDr"."_0.nc\n";
close metafile;

print "\n\ncpptraj control files is made\n\n";

##########################################
# make control file for DROIDS	
sleep(0.5);
print("Making DROIDS_deploy.ctl file...\n");
	if ($ribbon == 1 && $surface == 0) {$repStr = "ribbon";}  # opaque ribbon rep only
	if ($surface == 1 && $ribbon == 0) {$repStr =  "surface";} # opaque surface rep only
	if ($surface == 1 && $ribbon == 1) {$repStr =  "ribbonsurface";} # opaque ribbon with transparent surface

	$testStr = "flux"; $testStrLong = "fluctuation";  # file and folder labels
	
open(CTL, '>', "DROIDS_deploy_$fileIDq.ctl") or die "Could not open output file";
print CTL "query\t"."$fileIDq\t # Protein Data Bank ID for query structure\n";
#print CTL "reference\t"."$fileIDr\t # Protein Data Bank ID for reference structure (or neutral model)\n";
#print CTL "length\t"."$lengthID\t # number of amino acids on chain\n";
print CTL "num_chains\t"."$chainN\t # number of chains in structure\n";
$chainTTL = 0;
for(my $cnt = 0; $cnt < scalar @chainlen2; $cnt++){
    my $chain = chr($cnt + 65);
    #print "$cnt";
    #print "$chainlen2[$cnt]\n";
    #print "length$chain\t$chainlen2[$cnt]\n";
    print CTL "length$chain\t$chainlen2[$cnt]\t #end of chain designated\n";
    print "DROIDS_deploy.ctl\n";
    print "length$chain\t$chainlen2[$cnt]\t #end of chain designated\n";
    $chainTTL = $chainlen2[$cnt];
}
print CTL "length\t"."$chainTTL\t # total length of chain\n";
print CTL "start\t"."$startN\t # number of AA at start of chain\n";
#print CTL "cutoff_value\t"."$cutoffValue\t # p-value under which the KS comparison will be considered significant\n";
#print CTL "representations\t"."$repStr\t # methods of molecular representation in Chimera (ribbon and/or surface)\n";
#print CTL "test_type\t"."$testStr\t # test method (sequence = local Grantham dist, structure = RMSD, fluctuation = MD)\n";
#print CTL "color_scheme\t"."$colorType\t # output color scheme (red-green, yellow-blue, or orange-magenta)\n";
close CTL;
print("DROIDS ctl file is made\n");
##############################################
#  create list of chain labels
##############################################
@alphabet = ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z");
@chainlist = ();
if($chainN > 26) {print "warning - number of chains exceeds alphabet\n";}
for(my $l = 0; $l < $chainN; $l++){
     $letter = $alphabet[$l];
     push(@chainlist, $letter);
     }
print "chains in structure are...\n";
print @chainlist;
print "\n\n";

print "NOTE: if the chain designations look as if they have been calculated incorrectly\n";
print "you will need to edit...re-enter lengths manually in MDq_deploy.ctl, DROIDS_deploy.ctl\n\n";

close MUT;
} #end for loop
} # end outer for loop
print "\ncontrol files for all variants are made \n\n";

##############################################
}  # end sub


#####################################################################################################
sub align{

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
      print "\nmaking align file for $fileIDq.pdb\n";
      sleep(2);
      ############
      
print "STEP 1 - Here you will need to run MatchMaker in UCSF Chimera\n\n";
print "STEP 2 - Then run Match-Align in UCSF Chimera for each chain\n\n";
print "            if satisfied with alignment, save as a clustal file\n";
print "            (e.g. my_align.aln)\n\n";

print "continue? (y/n)\n";
my $go = <STDIN>;
chop($go);
if ($go eq "n") {exit;}
sleep(1);
print "            opening USCF Chimera and loading PDB ref structure\n\n";
print "            CREATE YOUR STRUCTURAL/SEQUENCE ALIGNMENT (.aln) NOW \n\n";
system("$chimera_path"."chimera $fileIDr"."REDUCED.pdb $fileIDq"."REDUCED.pdb\n");
# properly rename .aln file for DROIDS  
$chainlabel = '';
for (my $cl = 0; $cl < scalar @chainlist; $cl++){
     $chainlabel = $chainlist[$cl];
print "\nPlease enter name of your saved chain $chainlabel alignment file (e.g my$chainlabel"."_align.aln)\n";
my $align_name = "";
my $align_file = <STDIN>;
chop($align_file);
sleep(1);
# open and read first line header
open(IN, "<"."$align_file") or die "could not find alignment file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
      my $refINrow = $IN[$i+2];
	 my @INrow = split (/\s+/, $INrow);
	 my @refINrow = split (/\s+/, $refINrow);
      my $header = $INrow[0];
      my $ref_header = $refINrow[0];
      if ($header eq "CLUSTAL"){$align_name = $ref_header;}
      }
my @name_segment = split (/REDUCED/, $align_name);
if ($chainlabel eq "A"){$pdb_name = $name_segment[0];}
#$split_name = $name_segment[0]."_align".$chainlabel.".aln";
$split_name = $fileIDq."_align".$chainlabel.".aln";
$oldfilename = $align_file;
$newfilename = $split_name;
print "copying $align_file"." to $split_name\n";
# rename file with header
copy($oldfilename, $newfilename);
}
# create concatenated .aln master file
open (OUT, ">"."$pdb_name"."_align.aln");
for (my $ccl = 0; $ccl < scalar @chainlist; $ccl++){
     $chainlabel = $chainlist[$ccl];
     open (IN, "<"."$pdb_name"."_align".$chainlabel.".aln");
     my @IN = <IN>;
     print OUT @IN;
     print OUT "\n";
     #print @IN;
     #print "\n";
     close IN;
     }    
close OUT;

close MUT;
} #end for loop
} # end outer for loop
################
sleep(0.5);
print "\n\n alignment procedure is complete\n";
sleep(0.5);


	
}

#####################################################################################################

sub teLeap { # create topology and coordinate files 
system "perl teLeap_dnaproteinQuery_deploy.pl\n";
#system "perl teLeap_proteinReference.pl\n";
}

######################################################################################################
sub reduce { # create topology and coordinate files 

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
      print "\nmaking pdf4amber files for $fileIDq.pdb\n";
      sleep(2);
      ########
# copy + rename ref PDB before mutating 
system "pdb4amber -i $fileIDq.pdb -o ".$fileIDq."REDUCED.pdb --dry --reduce \n";
system "pdb4amber -i $fileIDr.pdb -o ".$fileIDr."REDUCED.pdb --dry --reduce \n";

close MUT;
} #end for loop
} # end outer for loop
}
######################################################################################################

sub launch { # launch MD run
if($simType eq "amber"){
    system "perl MD_proteinQuery_deploy.pl\n";
    sleep(2);
    print "\n\n";
    print "MD SIMULATIONS ARE COMPLETED\n\n";
    }
if($simType eq "open" && $solvType eq "ex"){
    system "conda config --set auto_activate_base true\n";
    system "x-terminal-emulator\n";
    print "\nRUN THE FOLLOWING SCRIPTS SEQUENTIALLY IN THE NEW TERMINAL\n";
    print "python MD_proteinQuery_deploy_openMM.py\n";
    sleep(2);
    print "\n\n";
    system "conda config --set auto_activate_base false\n";
    print "CLOSE TERMINAL WHEN BOTH MD SIMULATIONS ARE COMPLETED\n\n";
    }
if($simType eq "open" && $solvType eq "im"){
    print "Implicit solvent is not supported in OpenMM. Use explicit solvent\n\n";
    }
}

######################################################################################################

sub kill { # kill MD run
system "pkill pmemd\n";	
}

######################################################################################################

sub surv {
	### open job survalience terminals ######
system "x-terminal-emulator -e top\n";
system "x-terminal-emulator -e nvidia-smi -l 20\n";
}

######################################################################################################

###################################################################################################

sub info { # launch atom info

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
      print "\nmaking atom info file for $fileIDq.pdb\n";
      sleep(2);
      ##########
system("cpptraj "."-i ./atominfo_$fileIDq"."_deploy.ctl | tee cpptraj_atominfo_$fileIDq.txt");
#system("cpptraj "."-i ./atominfo_$fileIDr"."_0.ctl | tee cpptraj_atominfo_$fileIDr.txt");

close MUT;
} #end for loop
} # end outer for loop
}

###################################################################################################

sub flux { # launch atom fluctuation calc
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
      print "\nmaking atom flux file for $fileIDq.pdb\n";
      sleep(2);
      ########
      #for (my $i = 0; $i < $runsID; $i++){
      system("cpptraj "."-i ./atomflux_$fileIDq"."_deploy.ctl | tee cpptraj_atomflux_$fileIDq.txt");
      #system("cpptraj "."-i ./atomflux_$fileIDr"."_$i.ctl | tee cpptraj_atomflux_$fileIDr.txt");
      if($vector_enter eq 'y'){
      for(my $j = 0; $j<$allchainlen; $j++){
          system("cpptraj "."-i ./atomvctl_$fileIDq"."_deploy/atomvector_$fileIDq"."_aa$j.ctl | tee cpptraj_atomvector_$fileIDq.txt");
          }
      }
#  }

close MUT;
} #end for loop
} # end outer for loop
}
###################################################################################################

sub done {

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
      print "\nmaking final files for $fileIDq.pdb\n";
      sleep(2);
      ########

        # collect frame data
open (IN, "<"."MDframes.ctl") || die "could not open frames ctl file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $header = @INrow[0];
	 my $value = @INrow[1];
	 if ($header eq "framenumber"){$framenumber = $value*$prodLen;}
      if ($header eq "framestep"){$framestep = $value;}
      if ($header eq "framegroups"){$framegroups = $value;}
}
close IN;   
      
      
      
#print "Enter residue number at the start of both chains\n";
#print "(e.g. enter 389 if starts at THR 389.A) \n";
#print "(e.g. enter 1 if starts at MET 1.A) \n\n";
#my $startN = <STDIN>;
#chop($startN);

#my $fileIDr = $fileIDq;  # edit for retrieving alignment files for newly deployed simulations


sleep(2);
print "\n\n searching for atom info file = "."cpptraj_atominfo_$fileIDr.txt\n";
sleep(2);
print "\n\n creating atom_residue_list_$fileIDr.txt\n";
open(OUT, ">"."atom_residue_list_$fileIDr.txt") or die "could open output file\n";
print OUT "atomnumber\t"."atomlabel\t"."resnumber\t"."reslabel\n";
open(IN, "<"."cpptraj_atominfo_$fileIDr.txt") or die "could not find atom info file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 $atomnumber = @INrow[1];
	 $atomlabel = @INrow[2];
	 $resnumber = @INrow[3];
	 $resindex = $resnumber + ($startN - 1);
	 $reslabel = @INrow[4];
      if ($atomnumber eq "CA"|| $atomnumber eq "C" || $atomnumber eq "O" || $atomnumber eq "N"){ #finds correct whitespace frame when atomnumber > 10000
          $atomnumber = @INrow[0];
	     $atomlabel = @INrow[1];
	     $resnumber = @INrow[2];
	     $resindex = $resnumber + ($startN - 1);
	     $reslabel = @INrow[3];
          }
	 if ($atomlabel eq "CA"|| $atomlabel eq "C" || $atomlabel eq "O" || $atomlabel eq "N"){print OUT "$atomnumber\t $atomlabel\t $resindex\t $reslabel\n"}
   }
close IN;
close OUT;
sleep(2);
print "\n\n creating atom_residue_list_$fileIDq.txt\n";
open(OUT, ">"."atom_residue_list_unmodified_$fileIDq.txt") or die "could open output file\n";
print OUT "atomnumber\t"."atomlabel\t"."resnumber\t"."reslabel\n";
open(IN, "<"."cpptraj_atominfo_$fileIDq.txt") or die "could not find atom info file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
	 my @INrow = split (/\s+/, $INrow);
	 $atomnumber = @INrow[1];
	 $atomlabel = @INrow[2];
	 $resnumber = @INrow[3];
	 $resindex = $resnumber + ($startN - 1);
	 $reslabel = @INrow[4];
      if ($atomnumber eq "CA"|| $atomnumber eq "C" || $atomnumber eq "O" || $atomnumber eq "N"){ #finds correct whitespace frame when atomnumber > 10000
          $atomnumber = @INrow[0];
	     $atomlabel = @INrow[1];
	     $resnumber = @INrow[2];
	     $resindex = $resnumber + ($startN - 1);
	     $reslabel = @INrow[3];
          }
	 if ($atomlabel eq "CA"|| $atomlabel eq "C" || $atomlabel eq "O" || $atomlabel eq "N"){print OUT "$atomnumber\t $atomlabel\t $resindex\t $reslabel\n"}
   }
close IN;
close OUT;
sleep(2);
#################################################################
print "\n\n modifying alignment for color mapping to reference structure\n";
sleep(2);  # need to remove indels from ref sequence and any corresponding AA's in query

open(OUT1, ">"."$fileIDr"."_alignREFMAP.aln") or die "could not open output file\n";
open(IN1, "<"."$fileIDr"."_align.aln") or die "could not open alignment...did you save as (.aln)?\n";
print OUT1 "CLUSTAL W ALN saved from UCSF Chimera/MultAlignViewer\n\n";
my @IN1 = <IN1>;
my $position = 0;
for (my $i = 0; $i < scalar @IN1; $i++){
	my $IN1row = $IN1[$i];
	my $IN1nextrow = $IN1[$i+1];
     my $target = $fileIDr."REDUCED";
	if ($IN1row =~ m/$target/){my @IN1row = split(/\s+/, $IN1row); $header_ref = $IN1row[0]; $seq_ref =$IN1row[1]; print "$header_ref\t"."$seq_ref\n";
															my @IN1nextrow = split(/\s+/, $IN1nextrow); $header_query = $IN1nextrow[0]; $seq_query =$IN1nextrow[1]; print "$header_query\t"."$seq_query\n";
															my @seq_ref = split(//,$seq_ref);
															my @seq_query = split(//,$seq_query);
															my $new_seq_ref = "";
															my $new_seq_query = "";
															for (my $ii = 0; $ii < length $seq_ref; $ii++){
																      my $respos = $ii+1;
																			$position = $position+1;
																			my $AAref = @seq_ref[$ii]; 
																			my $AAquery = @seq_query[$ii];
																			if ($AAref ne "."){$new_seq_ref = $new_seq_ref.$AAref; $new_seq_query = $new_seq_query.$AAquery;}
															}
															print OUT1 "$header_ref\t"."$new_seq_ref\n";
															print OUT1 "$header_query\t"."$new_seq_query\n\n";
																													
															}
}
close OUT1;
close IN1;
sleep (2);

#################################################################
print "\n\n calculating AA sequence similarity and Grantham distances\n";
sleep(2);
print "\n\n creating vertical alignment files\n";
sleep(2);
open(OUT1, ">"."$fileIDr"."_vertalign_ref.aln") or die "could not open output file\n";
print OUT1 "respos\t"."seq_ref\n";
open(OUT2, ">"."$fileIDq"."_vertalign_query.aln") or die "could not open output file\n";
print OUT2 "respos\t"."seq_query\n";
open(OUT3, ">"."$fileIDr"."_vertalign_ref_indexed.aln") or die "could not open output file\n";
print OUT3 "respos\t"."seq_ref\n";
open(OUT4, ">"."$fileIDq"."_vertalign_query_indexed.aln") or die "could not open output file\n";
print OUT4 "respos\t"."seq_query\n";
open(OUT5, ">"."myGranthamDistances.txt") or die "could not open output file\n";
print OUT5 "respos\t"."distance\n";
open(IN1, "<"."$fileIDr"."_alignREFMAP.aln") or die "could not open alignment...did you save as (.aln)?\n";
my @IN1 = <IN1>;
my $position = 0;
my $positionINDEX = $startN-1;
my $AAsame_cnt = 0;
my $AA_cnt = 0;
my @gDISTS = ();
for (my $i = 0; $i < scalar @IN1; $i++){
	my $IN1row = $IN1[$i];
	my $IN1nextrow = $IN1[$i+1];
     my $target = $fileIDr."REDUCED";
	if ($IN1row =~ m/$target/){my @IN1row = split(/\s+/, $IN1row); $header_ref = $IN1row[0]; $seq_ref =$IN1row[1]; print "$header_ref\t"."$seq_ref\n";
															my @IN1nextrow = split(/\s+/, $IN1nextrow); $header_query = $IN1nextrow[0]; $seq_query =$IN1nextrow[1]; print "$header_query\t"."$seq_query\n";
															my @seq_ref = split(//,$seq_ref);
															my @seq_query = split(//,$seq_query);
															for (my $ii = 0; $ii < length $seq_ref; $ii++){
																      $gDIST = '';
																			my $respos = $ii+1;
																			my $AAref = @seq_ref[$ii]; 
																			my $AAquery = @seq_query[$ii];
																			$position = $position+1;
																			$positionINDEX = $positionINDEX+1;
																			if ($AAquery eq $AAref){$AAsame_cnt = $AAsame_cnt + 1;}
																			if ($AAquery ne $AAref || $AAquery eq $AAref){$AA_cnt = $AA_cnt + 1;}
																			open(IN2, "<"."amino1to3.txt") or die "could not open amino1to3.txt\n";
																			my @IN2 = <IN2>;
																			#print "AAref "."$AAref\n";
																			#print "AAquery "."$AAquery\n";
																			for (my $iii = 0; $iii < scalar @IN2; $iii++){
																					my $AArow = @IN2[$iii];
																					my @AArow = split(/\s+/, $AArow);
																					$AAone = @AArow[0]; $AAthree = @AArow[1];
																					if ($AAone eq $AAref){print OUT1 "$position\t"."$AAthree\n"}
																					if ($AAone eq $AAquery){print OUT2 "$position\t"."$AAthree\n"}
																					if ($AAone eq $AAref){print OUT3 "$positionINDEX\t"."$AAthree\n"}
																					if ($AAone eq $AAquery){print OUT4 "$positionINDEX\t"."$AAthree\n"}
																			  	}
																			# determine Grantham distance
																			if ($AAquery eq $AAref || $AAquery eq "." || $AAref eq "."){$gDIST = 0;}
																			if ($AAquery ne $AAref && $AAquery ne "." && $AAref ne "."){
																			open(IN3, "<"."GranthamScores.txt") or die "could not open amino1to3.txt\n";
																			my @IN3 = <IN3>;
																			for (my $iiii = 0; $iiii < scalar @IN3; $iiii++){
																					my $GDrow = @IN3[$iiii];
																					my @GDrow = split(/\s+/, $GDrow);
																					$AAqueryTEST = @GDrow[0]; $AArefTEST = @GDrow[1]; $gDISTtest = @GDrow[2];
																					if(uc $AAqueryTEST eq $AAquery && uc $AArefTEST eq $AAref){$gDIST = $gDISTtest;} # grantham matrix
																					elsif(uc $AAqueryTEST eq $AAref && uc $AArefTEST eq $AAquery){$gDIST = $gDISTtest;} # to cover other half of matrix
																					}
																			}
																			
																			print OUT5 "$position\t"."$gDIST\n";
																			if ($gDIST > 0) {push (@gDISTS, $gDIST);} # average only non-zero Grantham Distances
																																																									 
													        }
															}
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close IN1;
close IN2;

### whole sequence stats
$statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
           $statSCORE->add_data (@gDISTS);
					 $avg_gDIST = $statSCORE->mean();
					 $avg_gDIST = sprintf "%.2f", $avg_gDIST;

$AAseqsim = int(($AAsame_cnt/$AA_cnt+0.0001)*100);
$AAseq_matchfreq = ($AAsame_cnt/$AA_cnt+0.0001)*100;
$AAseq_matchfreq = sprintf "%.2f", $AAseq_matchfreq;
print "\n\nAA sequence similarity = "."$AAseqsim"."%\n";
print "avg Grantham Distance = "."$avg_gDIST"."\n";
open(OUT6, ">"."mySeqSTATS.txt") or die "could not open output file\n";
print OUT6 "label\t"."value\n";
print OUT6 "AAmatchFreq\t"."$AAseq_matchfreq\n";
print OUT6 "avgGranthamDist\t"."$avg_gDIST\n";
close OUT6;

sleep (2);


#################################################################################
print "\n\n adding gaps to query atom residue list (if needed)\n";
sleep(2);

open(IN1, "<"."atom_residue_list_unmodified_$fileIDq.txt") or die "could not open atom_residue_list.txt\n";
open(IN2, "<"."$fileIDq"."_vertalign_query_indexed.aln") or die "could not open output file\n";
open(OUT, ">"."atom_residue_list_$fileIDq.txt") or die "could not make atom_residue_list.txt\n";
#print OUT "atomnumber\t"."atomlabel\t"."resnumber\t"."reslabel\n";
my @IN1 = <IN1>;
my @IN2 = <IN2>;
$indelCount = 0;
for (my $i = 0; $i < scalar @IN1; $i++){ # scan residue list
	         my $IN1row = $IN1[$i];
			     my @IN1row = split(/\s+/, $IN1row);
			     my $atomnumber = $IN1row[0];
			     my $atomlabel = $IN1row[1];
			     my $resindex = $IN1row[2];
					 my $reslabel = $IN1row[3];
					 				 
					 for (my $j = 0; $j < scalar @IN2; $j++){ # scan alignment
			        my $IN2row = $IN2[$j];
			        my @IN2row = split(/\s+/, $IN2row);
			        my $pos_query = $IN2row[0] - $indelCount;
			        my $res_query = $IN2row[1];
							my $resindexgap = $resindex + $indelCount;
							#print "$pos_query\t"."$resindex\n";
				      if ($pos_query == $resindex && $res_query eq "xxx"){print OUT "na\t"."na\t"."na\t"."xxx\n"; print OUT "na\t"."na\t"."na\t"."xxx\n";print OUT "na\t"."na\t"."na\t"."xxx\n";print OUT "na\t"."na\t"."na\t"."xxx\n"; $indelCount = $indelCount+1}
							if ($pos_query == $resindex && $res_query ne "xxx"){print OUT "$atomnumber\t"."$atomlabel\t"."$resindexgap\t"."$reslabel\n";}		
			       
						 }
					 
					 #print "$reslabel\t".@skipped."\n";
					 
			}

close IN1;
close IN2;
close OUT;

#########################################################################################
##########  FLUX analysis     ###########################################################
#########################################################################################
print "\n\n collecting atomic fluctuation values (may take a minute)\n\n";
sleep(2);
open (OUT1, ">"."DROIDSfluctuation_deploy_$fileIDq.txt") or die "could not create output file\n";
print OUT1 "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
open(IN3, "<"."atom_residue_list_$fileIDr.txt") or die "could not open atom_residue_list.txt\n";
open(IN4, "<"."atom_residue_list_$fileIDq.txt") or die "could not open atom_residue_list.txt\n";
my @IN3 = <IN3>;
my @IN4 = <IN4>;
      for (my $i = 0; $i < scalar @IN3; $i++){ # scan atom type
			     my $IN3row = $IN3[$i];
	         my @IN3row = split(/\s+/, $IN3row); 
			     my $atomnumberR = $IN3row[0]; my $atomlabelR = $IN3row[1]; my $resnumberR = $IN3row[2]; my $reslabelR = $IN3row[3];
					 my $IN4row = $IN4[$i];
	         my @IN4row = split(/\s+/, $IN4row); 
			     my $atomnumberQ = $IN4row[0]; my $atomlabelQ = $IN4row[1]; my $resnumberQ = $IN4row[2]; my $reslabelQ = $IN4row[3];
					 #print "atom+res REF"."$atomnumberR\t"."$atomlabelR\t"."$resnumberR\t"."$reslabelR\n";	                  
					 #print "atom+res QUERY"."$atomnumberQ\t"."$atomlabelQ\t"."$resnumberQ\t"."$reslabelQ\n";
					 # assemble fluctuation data
			     for (my $ii = 0; $ii < $runsID; $ii++){  #scan flux data
	            $sample = $ii;
							open(IN5, "<"."fluct_$fileIDq"."_deploy.txt") or die "could not open fluct file for $fileIDq\n";
              open(IN6, "<"."fluct_$fileIDr"."_$ii.txt") or die "could not open fluct file for $fileIDr\n";
	            my @IN5 = <IN5>;
              my @IN6 = <IN6>;
			        $flux_query = '';
							$flux_ref = '';
							for (my $iii = 0; $iii < scalar @IN5; $iii++){
							    my $IN5row = $IN5[$iii];
									$IN5row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN5row = split(/\s+/, $IN5row);
									my $Qtest_atom_decimal = $IN5row[0];
									my $Qtest_atom = int($Qtest_atom_decimal);
							    #print "Q "."$Qtest_atom\t"."$atomnumberQ\n";
									if($atomnumberQ eq $Qtest_atom){$flux_query = $IN5row[1];}
							  }	
							for (my $iii = 0; $iii < scalar @IN6; $iii++){
									my $IN6row = $IN6[$iii];
									$IN6row =~ s/^\s+//;# need trim leading whitespace if present 
	                my @IN6row = split(/\s+/, $IN6row);
			            my $Rtest_atom_decimal = $IN6row[0];
									my $Rtest_atom = int($Rtest_atom_decimal);
									#print "R "."$Rtest_atom\t"."$atomnumberR\n";
									if($atomnumberR eq $Rtest_atom){$flux_ref = $IN6row[1];}
							  }
							
					    if($resnumberR =~/\d/ && $flux_query=~/\d/ && $reslabelQ ne "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."$flux_query\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."$flux_query\n";
							    }
							if($resnumberR =~/\d/ && $reslabelQ eq "xxx"){
							    #print "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."$flux_ref\t"."NA\n";
							    print OUT1 "$sample\t"."$resnumberR\t"."$reslabelR\t"."$reslabelQ\t"."$atomnumberR\t"."$atomlabelR\t"."NA\t"."NA\n";
							    }
							
							
							} 	
							
							close IN5;
              close IN6;
							
					
							}
					 
	
close IN3;
close IN4;
close OUT1;

########################################################################################
print "\n\n choose homology level for comparing backbone atom dynamics\n\n";
print " strict = collect only exact matching aligned residues\n";
print "          (e.g. position 5 -> LEU LEU)\n";
print "          (this will allow sites of mutations to be visualized later)\n\n";
print " loose  = collect any aligned residues\n";
print "          (e.g. position 5 -> LEU LEU or position 5 -> LEU ALA)\n"; 
print "          (this will NOT allow sites of mutations to be visualized later)\n\n";
## choose homology
#my $homology = <STDIN>;
#chop($homology);

$homology = "loose";
print "\nHOMOLOGY WILL BE LOOSE FOR THIS ANALYSIS\n\n";
sleep(2);

#$homology = "strict";
#print "\nHOMOLOGY WILL BE STRICT FOR THIS ANALYSIS\n\n";
#sleep(2);


open(CTL, '>>', "DROIDS.ctl") or die "Could not open output file";
print CTL "homology\t"."$homology\t # homology as 'strict' or 'loose'\n";
close CTL;

print "\n\n averaging DROIDSfluctuations by residue\n\n";
mkdir ("atomflux_deploy_$fileIDq") or die "please delete atomflux folder from previous run\n";
open (IN, "<"."DROIDSfluctuation_deploy_$fileIDq.txt") or die "could not create input file\n";
my @IN = <IN>;
open (OUT2, ">"."DROIDSfluctuationAVG_deploy_$fileIDq.txt") or die "could not create output file\n";
print OUT2 "pos_ref\t"."res_ref\t"."res_query\t"."flux_ref_avg\t"."flux_query_avg\t"."delta_flux\t"."abs_delta_flux\t"."KLdivergence\n";
@REFfluxAvg = ();
@QUERYfluxAvg = ();
$KL = 0;
for (my $j = 0; $j < scalar @IN; $j++){ # scan atom type
			     my $INrow = $IN[$j];
	         my @INrow = split(/\s+/, $INrow); 
			     my $sample = $INrow[0];
					 my $pos_ref = $INrow[1];
					 my $res_ref = $INrow[2];
					 my $res_query = $INrow[3];
					 my $atomnumber = $INrow[4];
					 my $atomlabel = $INrow[5];
					 my $flux_ref = $INrow[6];
					 my $flux_query = $INrow[7];
					 push(@REFfluxAvg, $flux_ref);
					 push(@QUERYfluxAvg, $flux_query);
					 my $INnextrow = $IN[$j+1];
	         my @INnextrow = split(/\s+/, $INnextrow); 
			     my $next_pos = $INnextrow[1];
					 print OUT "$sample\t"."$pos_ref\t"."$res_ref\t"."$res_query\t"."$atomnumber\t"."$atomlabel\t"."$flux_ref\t"."$flux_query\n";
					 
					 if ($homology eq "loose"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_query ne "xxx"){  # loose homology = collect all aligned residues  
           open (OUT, ">"."./atomflux_deploy_$fileIDq/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					                         
                          if ($pos_ref =~ m/\d/ && $j>1){
                              $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
                              $statSCORE->add_data (@REFfluxAvg);
					     $flux_ref_avg = $statSCORE->mean();
                              #$flux_ref_n = $statSCORE->count();
                              #print "flux_ref_n\t"."$flux_ref_n\n";
					     $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
                              $statSCORE->add_data (@QUERYfluxAvg);
					     $flux_query_avg = $statSCORE->mean();
                              #$flux_query_n = $statSCORE->count();
                              #print "flux_query_n\t"."$flux_query_n\n";
					     $delta_flux = ($flux_query_avg - $flux_ref_avg);
					     $abs_delta_flux = abs($flux_query_avg - $flux_ref_avg);
                              # calculate JS divergence
                              open (TMP1, ">"."flux_values_temp.txt") or die "could not create temp file\n";
                              print TMP1 "flux_ref\t"."flux_query\n";
                              for (my $t = 0; $t <= scalar @REFfluxAvg; $t++){print TMP1 "$REFfluxAvg[$t]\t"; print TMP1 "$QUERYfluxAvg[$t]\n";}
                              close TMP1;
                              open (TMP2, ">"."flux_values_KL.txt") or die "could not create temp file\n";
                              close TMP2;
                              #open (Rinput, "| R --vanilla")||die "could not start R command line\n";
                              #print Rinput "library('FNN')\n";
                              #print Rinput "data = read.table('flux_values_temp.txt', header = TRUE)\n"; 
                              #$flux_ref = "data\$flux_ref"; # flux on reference residue
                              #$flux_query = "data\$flux_query"; # flux on query residue
                              #print Rinput "d1 = data.frame(fluxR=$flux_ref, fluxQ=$flux_query)\n";
                              ##print Rinput "print(d1)\n";
                              #print Rinput "myKL<-KL.dist($flux_ref, $flux_query, k=10)\n";
                              #print Rinput "print(myKL[10])\n";
                              #print Rinput "sink('flux_values_KL.txt')\n";
                              #print Rinput "print(myKL[10])\n";
                              #print Rinput "sink()\n";
                              ## write to output file and quit R
                              #print Rinput "q()\n";# quit R 
                              #print Rinput "n\n";# save workspace image?
                              #close Rinput;
                              open (TMP3, "<"."flux_values_KL.txt") or die "could not create temp file\n";
                              my @TMP3 = <TMP3>;
                              for (my $tt = 0; $tt <= scalar @TMP3; $tt++){
                              $TMP3row = $TMP3[$tt];
                              @TMP3row = split (/\s+/, $TMP3row);
                              $header = $TMP3row[0];
                              $value = $TMP3row[1];
                              #print "$header\t"."$value\n";
                              if ($header eq "[1]"){$KL = $value;}
                              }
                              if ($delta_flux <= 0){$KL = -$KL;} # make KL value negative if dFLUX is negative
                              #print "my KL is "."$KL\n";
                              close TMP3;
                              print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\t"."$KL\n";
					     @REFfluxAvg = ();
                              @QUERYfluxAvg = ();
                              }
					 if ($next_pos eq ''){next;}
					 }}
					 					 
					 if ($homology eq "strict"){
					 if(($j == 1 || $pos_ref ne $next_pos) && $res_ref eq $res_query && $res_query ne "xxx"){ # strict homology = collect only exact matching residues  
           open (OUT, ">"."./atomflux_deploy_$fileIDq/DROIDSfluctuation_$next_pos.txt") or die "could not create output file\n";
           print OUT "sample\t"."pos_ref\t"."res_ref\t"."res_query\t"."atomnumber\t"."atomlabel\t"."flux_ref\t"."flux_query\n";
					         
                          if ($pos_ref =~ m/\d/ && $j>1){
                              $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - reference
                              $statSCORE->add_data (@REFfluxAvg);
					     $flux_ref_avg = $statSCORE->mean();
					     $statSCORE = new Statistics::Descriptive::Full; # residue avg flux - query
                              $statSCORE->add_data (@QUERYfluxAvg);
					     $flux_query_avg = $statSCORE->mean();
					     $delta_flux = ($flux_query_avg - $flux_ref_avg);
					     $abs_delta_flux = abs($flux_query_avg - $flux_ref_avg);
					     # calculate JS divergence
                              open (TMP1, ">"."flux_values_temp.txt") or die "could not create temp file\n";
                              print TMP1 "flux_ref\t"."flux_query\n";
                              for (my $t = 0; $t <= scalar @REFfluxAvg; $t++){print TMP1 "$REFfluxAvg[$t]\t"; print TMP1 "$QUERYfluxAvg[$t]\n";}
                              close TMP1;
                              open (TMP2, ">"."flux_values_KL.txt") or die "could not create temp file\n";
                              close TMP2;
                              #open (Rinput, "| R --vanilla")||die "could not start R command line\n";
                              #print Rinput "library('FNN')\n";
                              #print Rinput "data = read.table('flux_values_temp.txt', header = TRUE)\n"; 
                              #$flux_ref = "data\$flux_ref"; # flux on reference residue
                              #$flux_query = "data\$flux_query"; # flux on query residue
                              #print Rinput "d1 = data.frame(fluxR=$flux_ref, fluxQ=$flux_query)\n";
                              ##print Rinput "print(d1)\n";
                              #print Rinput "myKL<-KL.dist($flux_ref, $flux_query, k=10)\n";
                              #print Rinput "print(myKL[10])\n";
                              #print Rinput "sink('flux_values_KL.txt')\n";
                              #print Rinput "print(myKL[10])\n";
                              #print Rinput "sink()\n";
                              ## write to output file and quit R
                              #print Rinput "q()\n";# quit R 
                              #print Rinput "n\n";# save workspace image?
                              #close Rinput;
                              open (TMP3, "<"."flux_values_KL.txt") or die "could not create temp file\n";
                              my @TMP3 = <TMP3>;
                              for (my $tt = 0; $tt <= scalar @TMP3; $tt++){
                              $TMP3row = $TMP3[$tt];
                              @TMP3row = split (/\s+/, $TMP3row);
                              $header = $TMP3row[0];
                              $value = $TMP3row[1];
                              #print "$header\t"."$value\n";
                              if ($header eq "[1]"){$KL = $value;}
                              }
                              if ($delta_flux <= 0){$KL = -$KL;} # make KL value negative if dFLUX is negative
                              #print "my KL is "."$KL\n";
                              close TMP3;
                              print OUT2 "$pos_ref\t"."$res_ref\t"."$res_query\t"."$flux_ref_avg\t"."$flux_query_avg\t"."$delta_flux\t"."$abs_delta_flux\t"."$KL\n";
					     @REFfluxAvg = ();
                              @QUERYfluxAvg = ();
                              }
					 if ($next_pos eq ''){next;}
					 }}
					 
					 
																
}
close IN;
close OUT;
close OUT2;

sleep(2);

##################################################################################################
print "\n\n done parsing CPPTRAJ data files\n\n";
sleep(2);
#################################################################################################

# create chain ID column DROIDSfluctuationAVG.txt and make chain specific output data files

print " reading control file to get chain lengths\n\n";
@lengthlist = ();
$chainlabel = '';
for (my $cl = 0; $cl < scalar @chainlist; $cl++){
     $chainlabel = $chainlist[$cl];

my $AA_count = '';

open(IN, "<"."MDq_deploy_$fileIDq.ctl") or die "could not find CPPTRAJ input control file\n";
my @IN = <IN>;
for (my $c = 0; $c <= scalar @IN; $c++){
    my $INrow = $IN[$c];
    my @INrow = split (/\s+/, $INrow);
    my $header = $INrow[0];
    my $value = $INrow[1];
    #print "$header\t"."$value\n";
    if ($header eq "length$chainlabel") { $AA_count = $value; push (@lengthlist, $AA_count);}
}
close IN;
sleep(1);
}
print @chainlist;
print @lengthlist;
print "\n\n";
$pointer = 0;
$mychain = $chainlist[$pointer];
$mylength = $lengthlist[$pointer];
$prevlength = 0;
open(OUT, ">"."DROIDSfluctuationAVGchain_deploy_$fileIDq.txt") or die "could open DROIDS DATA file\n";
open(IN, "<"."DROIDSfluctuationAVG_deploy_$fileIDq.txt") or die "could not find DROIDS DATA file\n";
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++){
	 my $INrow = $IN[$i];
      chomp $INrow;
      my @INrow = split (/\s+/, $INrow);
	 if ($i == 0){print OUT "$INrow\t"."chain\n";}
      my $header = $INrow[0];
      #print "$header\t"."$mylength\t"."$mychain\n";
      if ($i > 0 && $header > $mylength){$pointer = $pointer+1; $mychain = $chainlist[$pointer]; $mylength = $lengthlist[$pointer];}
      if ($i > 0 && $header <= $mylength){print OUT "$INrow\t"."$mychain\n";}
      }
close IN;
sleep(1);
close OUT;
sleep(1);
print "chain lengths added to DROIDSfluctuationAVGchain_deploy_$fileIDq.txt file\n\n";
###############################################################
sleep(1);
my $frameCount = $framenumber;
my $stepsize = $framestep;
print "number of frames in time series = $frameCount\n";
sleep(1);
print "number of frames/step in time series = $framestep\n";
sleep(1);
if ($stepsize eq ''){$stepsize = 50;}

# choose topology file     
if($solvType eq "im"){$TOPfileQUERY = "vac_$fileIDq"."REDUCED.prmtop";}
if($solvType eq "ex"){$TOPfileQUERY = "wat_$fileIDq"."REDUCED.prmtop";}

# process .nc files on queryID
mkdir ("testingData_$fileIDq") or die "please delete testing data folder from previous run\n";
$TRAJfile = "prod_$fileIDq"."REDUCED_deploy.nc";
$OUTfile = "./testingData_$fileIDq/fluxtime_$fileIDq"."_deploy.txt";
$step = $stepsize;
$steplimit = $frameCount;
$start = 0;
$stop = $stepsize;
open (CPPTRAJ, "|"."cpptraj -p $TOPfileQUERY\n");
print CPPTRAJ "trajin $TRAJfile\n";
print CPPTRAJ "rms first\n";
print CPPTRAJ "average crdset MyAvg\n";
print CPPTRAJ "run\n";
print CPPTRAJ "rms ref MyAvg\n";
print CPPTRAJ "rms first average\n";
print CPPTRAJ "atomicfluct out $OUTfile \@CA,C,N,O,H&!(:WAT) start $start stop $steplimit\n";
print CPPTRAJ "run\n";
for(my $i = 0; $i<$steplimit/$step; $i++){
print CPPTRAJ "atomicfluct out $OUTfile \@CA,C,N,O,H&!(:WAT) start $start stop $stop\n";
print CPPTRAJ "run\n";
$start = $start + $step;
$stop = $stop + $step;
}
print CPPTRAJ "quit\n";
close CPPTRAJ;


sleep(1); print("time series file created\n\n"); sleep(1);
###############################################################
mkdir("./testingData_$fileIDq/indAAtest");
print("\nparsing testing set time series for each amino acid in $fileIDq...\n");
sleep(1);
open(OUT2, ">"."./testingData_$fileIDq/avgfluxtimeATOM_$fileIDq.txt")||die "could not create avg ATOM time series file\n";
print OUT2 "AA_position\t"."atomID\t"."ref_flux\n";
for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$t".".txt")||die "could not create AA time series file\n";  
   #for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./testingData_$fileIDq/fluxtime_$fileIDq"."_deploy".".txt")||die "could not open time series file "."fluxtime_$fileIDq"."_deploy.txt\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $DATAseries = substr($INrow, 28, length($INrow) - 28); # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	@DATAseries = split(/\s+/, $DATAseries);
     $aaID = int($a/4 - 0.1);
     $atomID = @INrow[1];
	$overallAVG = @INrow[2];
     if ($a == 0){  # create headers for file
          $framelistlength = scalar @DATAseries;
          $framelist = '';
          $framenumber = 0;
          for (my $aa = 0; $aa < $framelistlength; $aa++){$framenumber = $framenumber +1; $framelist = "$framelist"."$framenumber\t";}
          if ($tt == 0){print OUT "$framelist\n";}
          next;}
     if ($aaID == $t){print OUT $DATAseries;
                      $statSCORE = new Statistics::Descriptive::Full; # avg flux - reference
                      $statSCORE->add_data (@DATAseries);
	                 $ref_flux = $statSCORE->mean();
                      print OUT2 "$aaID\t"."$atomID\t"."$ref_flux\n";
                      }
     #print "$tt\t"."$aaID\t"."$atomID\t"."$overallAVG\n";}
     }
   close IN;
  #}
   close OUT;
 }
close OUT2;

#############################################################
if($vector_enter eq 'y'){

mkdir("./testingData_$fileIDq/indAAtest_vector");
print("\nparsing testing set time series for shape data for each amino acid in $fileIDq...\n");
sleep(3);
 for (my $t = 0; $t < $lengthID; $t++){
  open(OUT, ">"."./testingData_$fileIDq/indAAtest_vector/vecttimeAA_$fileIDq"."_$t".".txt")||die "could not create AA time series file\n";
  print OUT "X\t"."Y\t"."Z\t"."XO\t"."YO\t"."ZO\n";
   #for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./atomvect_$fileIDq"."_deploy/vect_$fileIDq"."_aa$t".".txt")||die "could not open vector file "."vect_$fileIDq"."_aa$t".".txt\n";
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
          print OUT "$meanX\t"."$meanY\t"."$meanZ\t"."$meanXO\t"."$meanYO\t"."$meanZO\n";
          @frameXavg = (); @frameYavg = (); @frameZavg = (); @frameXOavg = (); @frameYOavg = (); @frameZOavg = ();
          }
     }

     close IN;
 
   #}
   close OUT; 
 
 }
#################################

open(OUT3, ">"."./testingData_$fileIDq/avgfluxtimeAA_$fileIDq.txt")||die "could not open avg ATOM time series file\n";
print OUT3 "AA_position\t"."ref_flux\n";
@ref_flux = ();
open(IN3, "<"."./testingData_$fileIDq/avgfluxtimeATOM_$fileIDq.txt")||die "could not create avg AA time series file\n";
    my @IN3 = <IN3>;
    for (my $s = 0; $s < scalar @IN3; $s++) {
	if ($s == 0){next;}
     $IN3row = $IN3[$s];
     @IN3row = split(/\s+/, $IN3row);
	$atomPOS = @IN3row[0];
	$nextIN3row = $IN3[$s+1];
     @nextIN3row = split(/\s+/, $nextIN3row);
	$nextatomPOS = @nextIN3row[0];
     $ref_FLUX = @IN3row[2];
     if ($atomPOS == $nextatomPOS){push(@ref_flux, $ref_FLUX);}
     if ($atomPOS != $nextatomPOS){ # calc AA avg for ref flux
          push(@ref_flux, $ref_FLUX);
          $statSCORE = new Statistics::Descriptive::Full; # avg flux - reference
          $statSCORE->add_data (@ref_flux);
	     $avg_ref_FLUX = $statSCORE->mean();
          $len = scalar (@ref_flux);
          #print @ref_flux;
          #print "$len\n";
          print OUT3 "$atomPOS\t"."$avg_ref_FLUX\n";
          @ref_flux = ();}
    }

close IN2;
close OUT3;

#########################################################
mkdir("./testingData_$fileIDq/indAAtest_fluxvector");
print "\n\ncombining atom fluctuation and protein shape (i.e. vector ) information\n\n";
sleep(1);

for (my $t = 0; $t < $lengthID; $t++){
  open(OUT, ">"."./testingData_$fileIDq/indAAtest_fluxvector/fluxvectorAA_$fileIDq"."_$t".".txt")||die "could not create AA time series file\n";
  print OUT "value\t"."datatype\t"."run\n";
  for (my $tt = 0; $tt < $prodLen*$framegroups; $tt++){
        open(IN1, "<"."./testingData_$fileIDq/indAAtest/fluxtimeAA_$fileIDq"."_$t".".txt")||die "could not open flux file "."fluxtimeAA_$fileIDq"."_$t".".txt\n";
        open(IN2, "<"."./testingData_$fileIDq/indAAtest_vector/vecttimeAA_$fileIDq"."_$t".".txt")||die "could not open vector file "."vecttimeAA_$fileIDq"."_$t".".txt\n";
        my @IN1 = <IN1>;
        my @IN2 = <IN2>;
        for (my $a = 0; $a < scalar @IN1; $a++) {
	       $IN1row = $IN1[$a];
            @IN1row = split(/\s+/, $IN1row);
	       $flux = @IN1row[$tt];
            #print "AA "."$t"."  run "."$tt"."  class "."$class1"."  flux "."$flux\n";
            if($flux =~ m/\d/ && $flux != int($flux)){print OUT "$flux\t"."F\t"."$tt\n";}
            }
        for (my $aa = 0; $aa < scalar @IN2; $aa++){
               $IN2row = $IN2[$aa];
               @IN2row = split(/\s+/, $IN2row);
	          $x = @IN2row[0];
               $y = @IN2row[1];
               $z = @IN2row[2];
               $xo = @IN2row[3];
               $yo = @IN2row[4];
               $zo = @IN2row[5];
               #$flux = @IN1row[$tt+1];
               #print "AA "."$t"."  run "."$tt"."  class "."$class2"."  x "."$x"."  y "."$y"."  z "."$z"."  xo "."$xo"."  yo "."$yo"."  zo "."$zo\n";
               if($x =~ m/\d/){print OUT "$x\t"."X\t"."$tt\n"."$y\t"."Y\t"."$tt\n"."$z\t"."Z\t"."$tt\n";}
             }
      close IN1;
      close IN2;
      }
   close OUT;
   }

}

print("\nparsing is done\n");

close MUT;
} # end for loop
} # end outer for loop
###############################################################
} # end sub

##################################################################################################

