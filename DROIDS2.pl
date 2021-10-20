#!/usr/bin/perl
use Tk;
use Tk::Text;
use Tk::Font;
use Tk::PNG;
use Tk::JPEG;
use Tk::Photo;
#use strict;
#use warnings;
use feature ":5.10";
use File::Copy;
use List::Util qw( min );
use List::Util qw( max );
use List::Util qw(min max);
use Statistics::Descriptive();

# specify the path to working directory for Chimera here

my $chimera_path = "/opt/UCSF/Chimera64-1.11/bin/";

#### Declare variables ####
my $testType = '';
my $gpuType = '';
my $mdType = '';

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("DROIDS 4.0 - free software for comparative protein dynamics"); # Titles the main window
$mw->setPalette("gray20");
$image1 = $mw->Photo(-file => "droids.png");
$image2 = $mw->Photo(-file => "droidprotein.jpeg");

# Analysis Pipeline Frame

my $pipeFrame = $mw->Frame(	-label => "CHOOSE MODE OF MOLECULAR DYNAMIC COMPARISON",
				-relief => "groove",
				-borderwidth => 2
                );
        
   	 my $sdmRadio = $pipeFrame->Radiobutton( -text => "(2) analyze impact of one or several mutations (site-directed mutagenesis) on a protein or complex			            (requires 1 PDB ID and list of (locations/types) of mutation)",
						-foreground => 'lightyellow',
                        -value=>"sdm",
						-variable=>\$testType
						);
	 my $ppRadio = $pipeFrame->Radiobutton( -text => "(1) analyze impact of protein-protein interaction upon binding   			                                                                    (requires 2 PDB ID's representing the complex and unbound protein)",
						-foreground => 'lightyellow',
                        -value=>"pp",
						-variable=>\$testType
						);
     my $dp2Radio = $pipeFrame->Radiobutton( -text => "(2) analyze impact of mutation(s) on DNA-protein interaction		 		                 (requires 2 PDB ID's as in #4 and (locations/types) of cis or trans regulatory mutation)",
						-foreground => 'lightgreen',
                        -value=>"dp2",
						-variable=>\$testType
						);
	 
	 my $lp2Radio = $pipeFrame->Radiobutton( -text => "(3)  analyze impact of mutation(s) on protein-ligand interaction		 			   	            (requires 3 PDB ID's as in #7 and list of (locations/types) of mutation)",
						-foreground => 'pink',
                        -value=>"lp2",
						-variable=>\$testType
						);
     
	 
# System Type Frame
my $gpuFrame = $mw->Frame(	-label => "HOW MANY EFFECTIVE GPU's IN SYSTEM?",
				-relief => "groove",
				-borderwidth => 2
                );
	my $gpu1Radio = $gpuFrame->Radiobutton( -text => "single GPU system - run MD sequentially (single GPU or multi GPU with SLI)",
						-foreground => 'lightblue',
                        -value=>"gpu1",
						-variable=>\$gpuType
						);
        my $gpu2Radio = $gpuFrame->Radiobutton( -text => "double GPU system - run query/ref MD simultaneously (i.e. no SLI)",
						-foreground => 'lightgreen',
                        -value=>"gpu2",
						-variable=>\$gpuType
						);

		
# Buttons

my $pipeButton = $mw -> Button(-text => "run DROIDS", 
				-command => \&go,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $exitButton = $mw -> Button(-text => "exit DROIDS", 
				-command => \&stop,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $droidsButton1 = $mw -> Button(
                -command => \&stop,
                -background => 'gray45',
                -foreground => 'white',
                -image => $image1
				); # Creates an image
my $droidsButton2 = $mw -> Button(
                -command => \&stop,
                -background => 'gray45',
                -foreground => 'white',
                -image => $image2
				); # Creates an image


#### Organize GUI Layout ####
$pipeFrame->pack(-side=>"top",
		-anchor=>"s");


$sdmRadio->pack(-anchor=>"w");
$dp2Radio->pack(-anchor=>"w");
$lp2Radio->pack(-anchor=>"w");
$droidsButton1->pack(-side=>"right",
			-anchor=>"n");
$droidsButton2->pack(-side=>"left",
			-anchor=>"n");

$gpuFrame->pack(-side=>"top",
		-anchor=>"s");
$gpu1Radio->pack(-anchor=>"w");
$gpu2Radio->pack(-anchor=>"w");

$pipeButton->pack(-side=>"top",
			-anchor=>"s");
$exitButton->pack(-side=>"top",
			-anchor=>"s");




MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################

sub stop {exit;}

########################################################################################
sub go {
# make control file for DROIDS	
print("\nlaunching DROIDS 4.0...\n");
  
     if ($testType eq "ss" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSss.pl\n";}  
	elsif ($testType eq "sdm" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSsdm.pl\n";}  
	elsif ($testType eq "sdm" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSsdm_dualGPU.pl\n";}
	elsif ($testType eq "dp2" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSdp2.pl\n";}
    elsif ($testType eq "dp2" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSdp2_dualGPU.pl\n";}
	elsif ($testType eq "lp2" && $gpuType eq "gpu1") {system "perl GUI_START_DROIDSlp2.pl\n";}
    elsif ($testType eq "lp2" && $gpuType eq "gpu2") {system "perl GUI_START_DROIDSlp2_dualGPU.pl\n";}
    else {print " PLEASE SELECT OPTIONS\n"}
        
  
  
  
  
}

########################################################################
