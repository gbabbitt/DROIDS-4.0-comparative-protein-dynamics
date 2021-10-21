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

#### Declare variables ####
my $testType = '';
my $gpuType = '';
my $mdType = '';

#### Create GUI ####
my $mw = MainWindow -> new; # Creates a new main window
$mw -> title("maxDemon 2.0 - genetic/drug impacts on functionally conserved dynamics"); # Titles the main window
$mw->setPalette("gray20");
$image1 = $mw->Photo(-file => "demon.jpg");
$image2 = $mw->Photo(-file => "demon.jpg");

# Analysis Pipeline Frame

my $pipeFrame = $mw->Frame(	-label => "MULTI-AGENT MACHINE LEARNING STACK",
				-relief => "groove",
				-borderwidth => 2
                );
        
   my $chkKNN = $pipeFrame->Checkbutton( -text => "KNN - K nearest neighbors",
						-foreground => 'lightblue',
                        -onvalue=>"Red",
						-offvalue=>"",
						-variable=>\$knnType
						);
	 my $chkPROB = $pipeFrame->Checkbutton( -text => "naiveBayes/LDA/QDA - probablistic learning",
						-foreground => 'lightyellow',
                        -onvalue=>"Red",
						-offvalue=>"",
						-variable=>\$probType
						);
	 my $chkSVM = $pipeFrame->Checkbutton( -text => "SVM - support vector machine",
						-foreground => 'lightgreen',
                        -onvalue=>"Red",
						-offvalue=>"",
						-variable=>\$svmType
						);
     my $chkENS = $pipeFrame->Checkbutton( -text => "randomForest/adaboost - ensemble learning",
						-foreground => 'pink',
                        -onvalue=>"Red",
						-offvalue=>"",
						-variable=>\$ensType
						);
   
	 
	 
# System Type Frame
my $vectFrame = $mw->Frame(	-label => "SELECT OPTION FOR MACHINE LEARNING FEATURE VECTOR?",
				-relief => "groove",
				-borderwidth => 2
                );
	my $vect1Radio = $vectFrame->Radiobutton( -text => "learn atom fluctuations and downstream correlations",
						-foreground => 'lightblue',
                        -value=>"vect1",
						-variable=>\$vectType
						);
    my $vect2Radio = $vectFrame->Radiobutton( -text => "learn atom fluctuations only ",
						-foreground => 'lightgreen',
                        -value=>"vect2",
						-variable=>\$vectType
						);

		
# Buttons

my $pipe1Button = $mw -> Button(-text => "(1) run MD (validation and variants)", 
				-command => \&go1,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe2Button = $mw -> Button(-text => "(2) process train/test data sets", 
				-command => \&go2,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe3Button = $mw -> Button(-text => "(3) deploy machine learning models", 
				-command => \&go3,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe4Button = $mw -> Button(-text => "(4) calc/map conserved dynamics", 
				-command => \&go4,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe5Button = $mw -> Button(-text => "(5) calc/map coordinated dynamics", 
				-command => \&go5,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe6Button = $mw -> Button(-text => "(6) calc/map genetic/drug impacts", 
				-command => \&go6,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe7Button = $mw -> Button(-text => "(7) render movies of learning process", 
				-command => \&go7,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button
my $pipe8Button = $mw -> Button(-text => "(8) play movies of learning process", 
				-command => \&go8,
                -background => 'gray35',
                -foreground => 'white'
				); # Creates a go button											
my $exitButton = $mw -> Button(-text => "exit MAXDEMON", 
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

$chkKNN->pack(-anchor=>"w");
$chkPROB->pack(-anchor=>"w");
$chkSVM->pack(-anchor=>"w");
$chkENS->pack(-anchor=>"w");
$droidsButton1->pack(-side=>"right",
			-anchor=>"n");
$droidsButton2->pack(-side=>"left",
			-anchor=>"n");

$vectFrame->pack(-side=>"top",
		-anchor=>"s");
$vect1Radio->pack(-anchor=>"w");
$vect2Radio->pack(-anchor=>"w");

$pipe1Button->pack(-side=>"top",
			-anchor=>"s");
$pipe2Button->pack(-side=>"top",
			-anchor=>"s");
$pipe3Button->pack(-side=>"top",
			-anchor=>"s");
$pipe4Button->pack(-side=>"top",
			-anchor=>"s");
$pipe5Button->pack(-side=>"top",
			-anchor=>"s");
$pipe6Button->pack(-side=>"top",
			-anchor=>"s");
$pipe7Button->pack(-side=>"top",
			-anchor=>"s");
$pipe8Button->pack(-side=>"top",
			-anchor=>"s");								
$exitButton->pack(-side=>"top",
			-anchor=>"s");




MainLoop; # Allows Window to Pop Up

########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################

sub stop {exit;}

########################################################################################
sub go1 {

# write MLmethods.txt file
open (ML, ">"."MLmethods.txt");
if ($vectType eq "vect1"){print ML "shape\n";}
if ($vectType eq "vect2"){print ML "no_shape\n";}
if ($knnType eq "Red"){print ML "bnp\n";}
if ($knnType ne "Red"){print ML "no_bnp\n";}
if ($probType eq "Red"){print ML "dist\n";}
if ($probType ne "Red"){print ML "no_dist\n";}
if ($svmType eq "Red"){print ML "kern\n";}
if ($svmType ne "Red"){print ML "no_kern\n";}
if ($ensType eq "Red"){print ML "ens\n";}
if ($ensType ne "Red"){print ML "no_ens\n";}
#system "gedit MLmethods.txt";
print "\n\nPlease enter type of structure (1=protein | 2=DNA+protein | 3=protein+ligand)\n\n";
  $structType = <STDIN>; 
  chop($structType);
print "\nrunning validation and variant MD simulations\n"; 
if($structType == 1){system "perl GUI_MLMD_DROIDS.pl";}
if($structType == 2){system "perl GUI_MLMD_DROIDSdp.pl";}
if($structType == 3){system "perl GUI_MLMD_DROIDSlp.pl";}
}
########################################################################
sub go2 {
print "\nmaking control file and processing training data sets\n";
system "perl GUI_ML_DROIDSbtn2.pl";
}
########################################################################
sub go3 {
print "\nrunning stacked machine learning model\n"; 
if ($vectType eq "vect1"){system "perl GUI_ML_DROIDSbtn3_fluxCorr.pl";}
if ($vectType eq "vect2"){system "perl GUI_ML_DROIDSbtn3_fluxOnly.pl";}
}
########################################################################
sub go4 {
print "\nrunning canonical correlation and variant impact analyses\n";
system "perl GUI_ML_DROIDSbtn4.pl";
}
########################################################################
sub go5 {
print "\nshowing conserved dynamics on protein structure\n";
system "perl GUI_ML_DROIDSbtn5.pl";
}
########################################################################
sub go6 {
print "\nrendering movies\n";
system "perl GUI_ML_DROIDSbtn6.pl";
}
########################################################################
sub go7 {
print "\nshowing variant impacts on protein structure\n";
system "perl GUI_ML_DROIDSbtn7.pl"; 
}
########################################################################
sub go8 {
print "\nplaying movies of learners on protein structure\n";	
system "perl GUI_ML_DROIDSbtn8.pl"; 
}
########################################################################