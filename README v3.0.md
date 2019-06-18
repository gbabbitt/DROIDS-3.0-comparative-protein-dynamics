DROIDS v3.0 starts with the command line call

perl DROIDS.pl

a GUI based script PATHS.pl will prompt paths to working directories for UCSF Chimera, Amber forcefields, and Amber Home directory and write a .ctl file. Then the main DROIDS menu will pop up. Requires Amber 16/18 license and some dependencies in R, python and perl (see user doc with download)

Dr. Gregory A. Babbitt1 and Dr. Ernest P. Fokoue2 
1Thomas H. Gosnell School of Life Sciences, Rochester Institute of Technology, Rochester NY, USA 14623
2 School of Mathematical Sciences, Rochester Institute of Technology, Rochester NY, USA 14623

The Babbitt and Fokoue Labs at RIT have developed DROIDS v3.0, a software package for comparative protein dynamics, which applies metrics of distributional divergence and statistical analysis to the root mean square fluctuations (rmsf) of protein backbone atoms and maps these results to both static and moving image of proteins. We have also developed maxDemon v1.0, a multi-method machine learning application that trains on the comparative protein dynamics, identifies functionally conserved dynamics, and deploys classifications of functional dynamic states to newly generated protein simulations. Nine different types of machine learners can be deployed on the dynamics of each amino acid, then the resulting classifications are rendered upon movie images of the novel MD runs. This results in movies of protein dynamics where the conserved functional states are identified in real time by color mapping, allowing users to see both when and where a novel MD simulation displays a specific functional state defined by the comparative training. DROIDS+maxDemon designed to compare impacts of genetic variants and drug binding variants on the functional aspects of protein dynamics. 

Our main goal is to try to visualize the impact of one of the longest time scale processes in the universe (i.e molecular evolution over 100s millions of years) on one of the shortest time scale processes (i.e. molecular motion over femtoseconds). To achieve this goal we use state-of-the-art biophysical simulations and graphics to design a gaming PC into a ‘computational microscope’ that is capable seeing how mutations and other molecular events like binding, bending and bonding affect the functioning of proteins and nucleic acids. DROIDS-1.0 (Detecting Relative Outlier Impacts in molecular Dynamic Simulation) is a GUI-based pipeline that works with AMBER16, Chimera 1.11 and CPPTRAJ to analyze and visualize comparative protein dynamics on GPU accelerated Linux graphics workstations.  DROIDS employs a statistical method (multiple test corrected KS tests on all backbone atoms of each amino acid) to detect significant changes in molecular dynamics simulated on two homologous PDB structures.  Quantitative differences in atom fluctuation are displayed graphically and mapped onto movie images of the protein dynamics at the level of individual residues.  P values indicating significant changes are also able to be similarly mapped.  DROIDS is useful for examining how mutations or binding interactions affect protein dynamics.DROIDS was produced by student effort at the Rochester Institute of Technology under the direction of Dr. Gregory A. Babbitt as a collaborative project between the Gosnell School of Life Sciences and the Biomedical Engineering Dept.  Visit our lab website (https://people.rit.edu/gabsbi/) and download DROIDS 1.0 from Github at https://github.com/gbabbitt/DROIDS-1.0. We will be posting video results periodically on our youtube channel at https://www.youtube.com/channel/UCJTBqGq01pBCMDQikn566Kw


please cite

Biophys J. 2018 Mar 13;114(5):1009-1017. doi: 10.1016/j.bpj.2018.01.020.
DROIDS 1.20: A GUI-Based Pipeline for GPU-Accelerated Comparative Protein Dynamics.
Babbitt GA1, Mortensen JS2, Coppola EE2, Adams LE3, Liao JK2.
