#!/usr/bin/perl
use File::Copy;

print "INSTALLATION SCRIPT FOR DROIDS v3.0 DEPENDENCIES\n";
sleep(2);

print "DROIDS v3.0 main dependencies include Amber16/18, R, and UCSF Chimera on a Debian\n";
print "system Linux build or virtual machine with one or two dedicated Nvidia GPUs + CUDA 9.0\n";
print "NOTE: avoid CUDA 7.0 and CUDA 10.0+ with Amber\n";
sleep(2);

print "You will need the following file\n";
print "DROIDS-3.0.tar.gz (free from our website or GitHub repo)\n";
print "You may need the following files as well for fresh system install\n";
print "Amber18.tar.gz (licensed from https://ambermd.org/)\n";
print "AmberTools18.tar.gz or AmberTools19.tar.gz (free from https://ambermd.org/)\n";
print "chimera-1.14-linux_x86_64 or preferred version (free from https://www.cgl.ucsf.edu/chimera/)\n";
print "Modeller .deb folder (e.g. modeller_9.20-1_amd64.deb from https://salilab.org/modeller/)";
print "\n\nDo you need to open web links to these programs?\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "y" || $yn eq "Y" || $yn eq "yes"){
sleep(1);
system('xdg-open http://ambermd.org/'); sleep(1);
system('xdg-open https://www.cgl.ucsf.edu/chimera/'); sleep(1);
system('xdg-open https://salilab.org/modeller/'); sleep(1);
print "you might also like ChimeraX (chimerax-daily.deb) for VR application\n\n";
system('xdg-open https://www.rbvi.ucsf.edu/chimerax/'); sleep(1);
};
sleep(1);
print "Please answer the following questions about your VM instance or Linux system\n";
sleep(1);

print "\nIs amber and UCSF Chimera already installed on your system? (e.g. 'yes','y' or 'n','no')\n\n";
  $skipping = <STDIN>; 
  chop($skipping);
  
print "\nEnter your admin user name for this computer or VM\n\n";
  $UserName = <STDIN>; 
  chop($UserName);

print "\nEnter chimera version you use or plan to install here (e.g. '1.11' or '1.13')\n\n";
  $ChimeraName = <STDIN>; 
  chop($ChimeraName);


sleep(1); print "\nchecking perl version - looking for Perl 5.0\n\n"; system('perl -v'); sleep(1);

#install Perl modules
sleep(1); print "\ninstalling perl modules\n\n"; sleep(1);system('sudo cpan App::cpanminus'); system ('sudo apt install cpanminus'); system('sudo cpanm Statistics::Descriptive'); sleep(1);

# install and update Debian packages
sleep(1); print "\nchecking Debian packages\n\n"; sleep(1);
print "\nDo you need to install the xfce4 desktop, xrdp, and VNC? (y/n) (this is needed for GCP VM instance)?\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "y" || $yn eq "Y" || $yn eq "yes"){sleep(1); print "\ninstalling xrdp and xfce4 desktop\n\n"; sleep(1); system('sudo apt-get install xrdp xfce4'); sleep(1);}
if($yn eq "y" || $yn eq "Y" || $yn eq "yes"){sleep(1); print "\ninstalling VNC\n\n"; sleep(1); system('sudo apt install tightvncserver'); sleep(1);}
sleep(1); print "\ninstalling htop\n\n"; sleep(1); system('sudo apt-get install htop'); sleep(1);
sleep(1); print "\ninstalling gedit\n\n"; sleep(1); system('sudo apt-get install gedit'); sleep(1);
sleep(1); print "\ninstalling gdebi\n\n"; sleep(1); system('sudo apt-get install gdebi'); sleep(1);
sleep(1); print "\ninstalling gparted\n\n"; sleep(1); system('sudo apt-get install gparted'); sleep(1);
sleep(1); print "\ninstalling vokoscreen\n\n"; sleep(1); system('sudo apt-get install vokoscreen'); sleep(1);
sleep(1); print "\ninstalling evince\n\n"; sleep(1); system('sudo apt-get install evince'); sleep(1);
sleep(1); print "\ninstalling grace\n\n"; sleep(1); system('sudo apt-get install grace'); sleep(1);
sleep(1); print "\ninstalling perl-tk\n\n"; sleep(1); system('sudo apt-get install perl-tk'); sleep(1);
sleep(1); print "\ninstalling python-tk\n\n"; sleep(1); system('sudo apt-get install python-tk'); sleep(1);
sleep(1); print "\ninstalling python-gi\n\n"; sleep(1); system('sudo apt-get install python-gi'); sleep(1);
sleep(1); print "\ninstalling gstreamer\n\n"; sleep(1); system('sudo apt-get install libgstreamer1.0-0 gstreamer1.0-plugins-base gstreamer1.0-plugins-good gstreamer1.0-plugins-bad gstreamer1.0-plugins-ugly gstreamer1.0-libav gstreamer1.0-doc gstreamer1.0-tools gstreamer1.0-x gstreamer1.0-alsa gstreamer1.0-gl gstreamer1.0-gtk3 gstreamer1.0-qt5 gstreamer1.0-pulseaudio'); sleep(1);
sleep(1); print "\ninstalling steam and VR dependencies\n\n"; sleep(1); system('sudo apt-get install steam steam-devices libvulkan1'); sleep(1);
sleep(1); print "\ninstalling Amber dependencies\n\n"; sleep(1); system('sudo apt-get install csh flex patch gfortran g++ make xorg-dev bison libbz2-dev'); sleep(1);
sleep(1); print "\nrunning updates\n\n"; sleep(1); system('sudo apt-get update'); sleep(1);

# skip Amber install option
if($skipping eq "y" || $skipping eq "Y" || $skipping eq "yes"){print "\n amber installation skipped\n\n"; goto Askip;}

sleep(1); print "\nyou should have tar.bz2 versions of Ambertools18/19 and Amber18 on desktop\n\n"; 
sleep(1); print "\ninstalling AmberTools18\n\n";
print "\nAre you using AmberTools 18 or 19? (type '18 or '19')\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "19"){sleep(1); print "\nunzipping ambertools\n\n"; sleep(1); system('tar jxvf AmberTools19.tar.bz2'); sleep(1);}
if($yn eq "18"){sleep(1); print "\nunzipping ambertools\n\n"; sleep(1); system('tar jxvf AmberTools18.tar.bz2'); sleep(1);}
system('export AMBERHOME=/home/'.$UserName.'/Desktop/amber18');
system('gnome-terminal');
sleep(1); print "\nin the new terminal run the following commands\n\n"; sleep(1);
print "cd amber18\n";
print "./configure -noX11 gnu\n";
print "make install\n";
#print "make test\n\n";
sleep(1); print "\nWhen AmberTools is installed, close secondary terminal\n\n"; sleep(1);
print "\nAre you ready to continue? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}

# update AMBERHOME on bashrc
sleep(1); print "\nupdating bashrc file\n"; sleep(1);
sleep(1); print "\nyou need to add the following lines to the bashrc file\n\n"; sleep(1);
print "source /home/".$UserName."/Desktop/amber18/amber.sh\n";
print "export AMBERHOME=/home/".$UserName."/Desktop/amber18\n";
print "export PATH="."\$PATH:"."\$AMBERHOME/bin\n";
#print "test -f /home/(host name)/Desktop/amber18/amber.sh\n";
system('gedit ~/.bashrc');
print "\nAre you ready to continue? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}

# run tests on AmberTools
#system('gnome-terminal');
#sleep(1); print "\nin the new terminal run the following commands\n\n"; sleep(1);
#print "cd amber18\n";
#print "make test\n\n";
#sleep(1); print "\nWhen AmberTools is done testing, close secondary terminal\n\n"; sleep(1);
#print "\nAre you ready to continue? (y/n)\n\n";
#  $yn = <STDIN>; 
#  chop($yn);
#if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}

# install CUDA
print "\nIs CUDA and cuda toolkit already installed? (y/n) (can probably skip if GCP cloud VM) \n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "y" || $yn eq "Y" || $yn eq "yes"){print "\n CUDA installation skipped\n\n"; goto CDskip;}
#install cuda
sleep(1); print "\ninstalling cuda tool kit\n\n"; sleep(1); system('sudo apt install nvidia-cuda-toolkit'); sleep(1);
sleep(1); print "\nyou will need to download Linux .deb version of CUDA\n\n"; sleep(1); 
sleep(1); print "\nIMPORTANT NOTE: double check Amber webpage for supported CUDA versions (avoid CUDA 7.0 and 10.0)\n\n";
print "/nDo you want to open CUDA download webpages (y/n)\n";
$yn = <STDIN>; 
  chop($yn);
if($yn eq "y" || $yn eq "Y"){system('xdg-open https://developer.nvidia.com/cuda-toolkit-archive'); sleep(1);}

# install CUDA
system('gnome-terminal');
sleep(1); print "\nin the new terminal run the following commands\n\n"; sleep(1);
print "sudo dpkg -i '.deb package name'\n";
print "sudo apt-get update\n";
print "sudo apt-get install cuda\n\n";
print "\ninstall patches manually if needed\n";
print "\nNOTE: double check that the downloaded .deb file permissions allow executable\n\n";
sleep(1); print "\nWhen CUDA is installed, close secondary terminal\n\n"; sleep(1);
print "\nAre you ready to continue? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}
CDskip:

# update cuda home
sleep(1); print "\nupdating bashrc file\n"; sleep(1);
print "\nyou need to add the following lines to the bashrc file\n\n";
print "export CUDA_HOME=/usr/local/cuda-(version e.g. 9.0)\n";
print "export PATH=\$CUDA_HOME/bin:/usr/local/cuda-(version e.g. 9.0)/bin:\$PATH\n";
print "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\${AMBERHOME}/lib:\$CUDA_HOME/lib64:\$CUDA_H\$\n";
print "test -f /home/".$UserName."/Desktop/amber18/amber.sh  && source /home/".$UserName."/Desktop/amber18/amber.sh \n";
system('gedit ~/.bashrc');
print "\nAre you ready to continue? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}


# install Amber18
sleep(1); print "\ninstalling Amber18 (pmemd.cuda)\n\n";
sleep(1); print "\nchecking nvcc and c compilers\n\n"; sleep(1);
system('nvcc -V');
system('gcc --version');
sleep(1);
sleep(1); print "\nunzipping amber18\n\n"; sleep(1); system('tar jxvf Amber18.tar.bz2'); sleep(1);
system('gnome-terminal');
sleep(1); print "\nin the new terminal run the following commands\n\n"; sleep(1);
print "cd amber18\n";
print "./configure -cuda gnu\n";
print "make install\n";
print "make test\n\n";

print "\nNOTE: if gcc compilers are too recent and 'make install' fails\n\n"; sleep(1);
print "sudo apt-get install gcc-5 g++-5 gfortran-5\n";
print "sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 60 --slave /usr/bin/g++ g++ /usr/bin/g++-5 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5\n";
print "sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 40 --slave /usr/bin/g++ g++ /usr/bin/g++-7 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-7 \n";
print "sudo update-alternatives --config gcc\n";

sleep(1); print "\nWhen Amber18 is installed, close second terminal\n\n"; sleep(1);
print "\nAre you ready to continue? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}

# copy and duplicate pmemd.cuda
sleep(1); print "copying and duplicating pmemd.cuda for dual GPUs\n"; sleep(1); 
copy("./amber18/bin/pmemd.cuda_SPFP", "./amber18/bin/pmemd0.cuda_SPFP")||die "could not find and copy pmemd.cuda\n";
copy("./amber18/bin/pmemd.cuda_SPFP", "./amber18/bin/pmemd1.cuda_SPFP")||die "could not find and copy pmemd.cuda\n";
copy("./amber18/bin/pmemd.cuda_DPFP", "./amber18/bin/pmemd0.cuda_DPFP")||die "could not find and copy pmemd.cuda\n";
copy("./amber18/bin/pmemd.cuda_DPFP", "./amber18/bin/pmemd1.cuda_DPFP")||die "could not find and copy pmemd.cuda\n";
system('chmod +x ./amber18/bin/pmemd0.cuda_SPFP');
system('chmod +x ./amber18/bin/pmemd1.cuda_SPFP');
system('chmod +x ./amber18/bin/pmemd0.cuda_DPFP');
system('chmod +x ./amber18/bin/pmemd1.cuda_DPFP');
Askip: # skip Amber install option


if($skipping eq "y" || $skipping eq "Y" || $skipping eq 'yes'){print "\n Chimera installations skipped\n\n"; goto Cskip;}

# install Chimera, ChimeraX and Modeller
sleep(1); print "\ninstalling UCSF Chimera binary file, and chimerax-daily tar.gz (option) and Modeller(option) .deb folders to your desktop\n\n"; sleep(1);
print "\n NOTE: make sure to have right clicked these, go to properties/permissions and allow them to run as executable\n\n"; sleep(1);
print "\nAre you ready? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}

print "\ntype name of Chimera bin file (e.g. chimera-1.13.1-linux_x86_64.bin)  NOTE:be sure file permission is set executable\n\n";
  $Cname = <STDIN>; 
  chop($Cname);
sleep(1); print "\ninstalling Chimera\n\n"; sleep(1); system('./'.$Cname); sleep(1);
#sleep(1); print "\nunzipping ChimeraX\n\n"; sleep(1); system('tar xvzf chimerax-daily.tar.gz'); sleep(1);
#sleep(1); print "\nIMPORTANT NOTE: ChimeraX runs from within extracted bin folder\n\n"; sleep(1);
print "\ntype name of Modeller .deb folder (e.g. modeller_9.20-1_amd64.deb)\n\n";
  $Mname = <STDIN>; 
  chop($Mname);
sleep(1); print "\ninstalling Modeller\n\n"; sleep(1); system('sudo apt install ./'.$Mname); sleep(1);
# install ChimeraX
print "\ninstall chimeraX \n";
system('gnome-terminal');
sleep(1); print "\nin the new terminal run the following commands\n\n"; sleep(1);
print "sudo dpkg -i chimerax-daily.deb\n";
print "sudo apt-get update\n";
print "sudo apt-get install chimerax-daily\n\n";

sleep(1); print "\nWhen chimeraX is installed, close secondary terminal\n\n"; sleep(1);
print "\nAre you ready to continue? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;}

Cskip:

# check R installation
sleep(1); print "\ninstalling R\n\n"; sleep(1); system('sudo apt-get install r-base-core r-base r-base-dev'); sleep(1);
sleep(1); print "\ninstalling R and R packages\n\n";

sleep(1); print "\ninstalling some R packages  (type 'y' if this hangs)\n\n";
print "\nIf script fails here, open terminal and install R (sudo apt-get install r-base-core r-base r-base-dev) then type 'R' at command line, then 'install.packages('package name') ...names = ggplot2 gridExtra dply caret FNN e1071 kernlab class MASS ada randomForest CCA CCP doParallel foreach...then 'q()'\n\n";

#install R and R packages
open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "chooseCRANmirror(graphics = getOption('menu.graphics'), ind = 81, local.only = TRUE)\n";
print Rinput "install.packages('ggplot2')\n";
print Rinput "install.packages('gridExtra')\n";
print Rinput "install.packages('dplyr')\n";
print Rinput "install.packages('caret')\n";
print Rinput "install.packages('FNN')\n";
print Rinput "install.packages('e1071')\n";
print Rinput "install.packages('kernlab')\n";
print Rinput "install.packages('class')\n";
print Rinput "install.packages('MASS')\n";
print Rinput "install.packages('ada')\n";
print Rinput "install.packages('randomForest')\n";
print Rinput "install.packages('CCA')\n";
print Rinput "install.packages('CCP')\n";
print Rinput "install.packages('doParallel')\n";
print Rinput "install.packages('foreach')\n";
print Rinput "install.packages('boot')\n";
# load some libraries to check installation
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";
print Rinput "library(lattice)\n";
print Rinput "library(FNN)\n";
print Rinput "library(MASS)\n";
print Rinput "library(CCA)\n";
print Rinput "library(CCP)\n";
print Rinput "library(e1071)\n";
print Rinput "library(kernlab)\n";
print Rinput "library(class)\n";
print Rinput "library(caret)\n";
print Rinput "library(dplyr)\n";
print Rinput "library('ada')\n";
print Rinput "library('randomForest')\n";
print Rinput "library('parallel')\n";
print Rinput "library('foreach')\n";
print Rinput "library('doParallel')\n";
print Rinput "library('boot')\n";
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";

#unzip and move DROIDS
sleep(1); print "\nlooking for DROIDS-3.0.tar.gz folder on desktop (rename download file if needed)\n\n"; sleep(1);
print "\nAre you ready? (y/n)\n\n";
  $yn = <STDIN>; 
  chop($yn);
if($yn eq "n" || $yn eq "N"){print "\ninstallation interrupted\n\n"; exit;} 
$Dname = 'DROIDS-3.0.tar.gz'; 
sleep(1); print "\nunzipping DROIDS\n\n"; sleep(1); system('tar xvzf '.$Dname); sleep(1);

# find chimera paths
sleep(1); print "\nlocating paths to steam, chimera and chimerax\n";
sleep(1); system ('which steam');
sleep(1); system ('locate /bin/ChimeraX | egrep ./ | grep bin');
sleep(1); system ('locate /bin/chimera | egrep ./ | grep bin');
# create paths.ctl
sleep(1); print "\ncreating paths.ctl file\n";
open (PTH, ">"."paths.ctl" || die "could not create paths.ctl file\n");
print PTH "amber_path	/home/$UserName/Desktop/amber18/	# path to amber home folder\n";
print PTH "chimera_path	~/.local/UCSF-Chimera64-$ChimeraName/bin/	# path to Chimera executable\n";
print PTH "chimerax_path	/home/$UserName/Desktop/chimerax-2019.01.19/bin/	# path to ChimeraX executable\n";
print PTH "teleap_path	/home/$UserName/Desktop/amber18/dat/leap/cmd/	# path to teLeap force field folder\n";
print PTH "steam_path	/usr/bin/steam 	# path to steam executable\n";
print "\ndouble check Chimera working paths and version numbers in paths.ctl file, edit if needed, then copy it manually into your DROIDS folder\n\n"; 
print "\nDO NOT INCLUDE FILENAME FOR EXECUTABLE IN THE PATH\n"; sleep(2);

system ('gedit paths.ctl');
print "\n\nDROIDS v3.0 INSTALLATION COMPLETE\n\n";

print "\n(optional) to install md4vr, remove md4vr.zip from DROIDS folder to desktop and unzip\n";
print "make sure steamOS and steamVR are working correctly,then follow README\n\n";

exit;  


