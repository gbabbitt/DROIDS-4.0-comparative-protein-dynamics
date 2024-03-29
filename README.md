# DROIDS-4.0-comparative-protein-dynamics
DROIDS is a machine learning assisted GUI based AMBER MD pipeline for statistical comparison of protein dynamics (python DROIDS.py). The DROIDS program runs multiple replicate MD runs in two comparative states (e.g. ligand bound vs unbound, DNA bound vs unbound, wild-type vs mutated, monomer vs dimer, or two temperatures) and statistically compares and color-maps local atom fluctuations across the protein. The subsequent maxDemon program trains a multi-agent machine learning algorithm upon the atom fluctuations in these two comparative states and then recognizes them in subsequent validation runs on molecular homologs in order to identify regions of functionally or evolutionary conserved dynamics. Information theoretics are also employed to compare genetic variant impacts on conserved dynamics and to detect coordinated protein motions. DROIDS 4.0 adds an option to train the machine learners on both atom fluctuation AND atom correlations with motions of adjacent amino acids. This potentially captures more complex intrinsic sequence dependent motions than atom fluctuation alone. This release also includes amberMDgui, a simple GUI for controlling replicate MD runs on a single PDB structure (python AMBERMD.py). This release also includes PDBmutator, a python extension of UCSF Chimera that creates all possible single substitution/replacement variants from a PDB structure of DNA/protein.  it also can create all variants represented in a CLUSTAL protein alignment file (.aln).  To run unzip the PDBmutator folder and run 'python PDBMUTATOR.py' from terminal.

Please see the following citations for more details
DROIDS 1.20: A GUI-Based Pipeline for GPU-Accelerated Comparative Protein Dynamics
Babbitt, Gregory A. et al. 2020
Biophysical Journal, Volume 114, Issue 5, 1009 - 1017

DROIDS 3.0—Detecting Genetic and Drug Class Variant Impact on Conserved Protein Binding Dynamics
Babbitt, Gregory A. et al. 2018
Biophysical Journal, Volume 118, Issue 3, 541 - 55
![image](/DROIDSgui.png)

To run DROIDS v4.0:
python DROIDS.py OR perl DROIDS.pl

a GUI based script PATHS.pl will prompt paths to working directories for UCSF Chimera, Amber forcefields, and Amber Home directory and write a .ctl file. Then the main DROIDS menu will pop up. Requires Amber 16/18 license and some dependencies in R, python and perl (see user doc with download)

Dr. Gregory A. Babbitt(1) and Dr. Ernest P. Fokoue(2) 
1Thomas H. Gosnell School of Life Sciences, Rochester Institute of Technology, Rochester NY, USA 14623
2 School of Mathematical Sciences, Rochester Institute of Technology, Rochester NY, USA 14623
![image](/DROIDSgui3.png)

The Babbitt and Fokoue Labs at RIT have developed DROIDS 4.0, a software package for comparative protein dynamics, which applies metrics of distributional divergence and statistical analysis to the root mean square fluctuations (rmsf) of protein backbone atoms and maps these results to both static and moving image of proteins. We have also developed maxDemon 2.0, a multi-method machine learning application that trains on the comparative protein dynamics, identifies functionally conserved dynamics, and deploys classifications of functional dynamic states to newly generated protein simulations. Nine different types of machine learners can be deployed on the dynamics of each amino acid, then the resulting classifications are rendered upon movie images of the novel MD runs. This results in movies of protein dynamics where the conserved functional states are identified in real time by color mapping, allowing users to see both when and where a novel MD simulation displays a specific functional state defined by the comparative training. DROIDS+maxDemon designed to compare impacts of genetic variants and drug binding variants on the functional aspects of protein dynamics. 
![image](/MAXDEMONgui.png)

Our main goal is to try to visualize the impact of one of the longest time scale processes in the universe (i.e molecular evolution over 100s millions of years) on one of the shortest time scale processes (i.e. molecular motion over femtoseconds). To achieve this goal we use state-of-the-art biophysical simulations and graphics to design a gaming PC into a ‘computational microscope’ that is capable seeing how mutations and other molecular events like binding, bending and bonding affect the functioning of proteins and nucleic acids. DROIDS (Detecting Relative Outlier Impacts in molecular Dynamic Simulation) is a GUI-based pipeline that works with AMBER18, Ambertools18, Chimera 1.11, and R to analyze and visualize comparative protein dynamics on GPU accelerated Linux graphics workstations.  DROIDS employs a statistical method (multiple test corrected KS tests on all backbone atoms of each amino acid) to detect significant changes in molecular dynamics simulated on two homologous PDB structures.  Quantitative differences in atom fluctuation are displayed graphically and mapped onto movie images of the protein dynamics at the level of individual residues.  P values indicating significant changes are also able to be similarly mapped.  DROIDS is useful for examining how mutations or binding interactions affect protein dynamics.DROIDS was produced by student effort at the Rochester Institute of Technology under the direction of Dr. Gregory A. Babbitt as a collaborative project between the Gosnell School of Life Sciences and the Biomedical Engineering Dept.  Visit our lab website (https://people.rit.edu/gabsbi/) or download DROIDS 3.0 from Github. We will be posting video results periodically on our youtube channel at https://www.youtube.com/channel/UCJTBqGq01pBCMDQikn566Kw

If you use any software included in the DROIDS project please cite

DROIDS 1.20: A GUI-Based Pipeline for GPU-Accelerated Comparative Protein Dynamics
Babbitt, Gregory A. et al. 2020
Biophysical Journal, Volume 114, Issue 5, 1009 - 1017

DROIDS 3.0—Detecting Genetic and Drug Class Variant Impact on Conserved Protein Binding Dynamics
Babbitt, Gregory A. et al. 2018
Biophysical Journal, Volume 118, Issue 3, 541 - 55
