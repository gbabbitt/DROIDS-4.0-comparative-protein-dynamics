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

#### Introductory Message #############
system 'gedit READMEv4.0.md';


print " 

    DROIDS 4.0 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DROIDS 4.0 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DROIDS 4.0.  If not, see <http://www.gnu.org/licenses/>.

    Visit us on GitHub. \n\n";

print "continue to GUI ? (y/n)\n";
my $go = <STDIN>;
chop($go);

if ($go eq "n") {exit;}

#### This creates a GUI to write the control files needed for the GPU accelerated pmemd.cuda pipeline ####

#### This uses a GUI to write the control files needed for the DROIDS scripts ####
print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
                          in Dynamic Simulations

- analytical engine and visual toolbox for functional evolutionary
  comparison of molecular dynamic simulation \n\n";

#### open PATHS specification ####

system "perl PATHS.pl\n";

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

my $pipeFrame = $mw->Frame(	-label => "MACHINE LEARNING ASSISTED COMPARATIVE PROTEIN DYNAMICS",
				-relief => "groove",
				-borderwidth => 2
                );
        
   my $oneRadio = $pipeFrame->Radiobutton( -text => "DROIDS - direct pairwise comparative analyses",
						-foreground => 'lightblue',
                        -value=>"one",
						-variable=>\$testType
						);
	 my $twoRadio = $pipeFrame->Radiobutton( -text => "DROIDS - mutant model comparative analyses",
						-foreground => 'lightyellow',
                        -value=>"two",
						-variable=>\$testType
						);
	 my $threeRadio = $pipeFrame->Radiobutton( -text => "DROIDS + maxDemon - functional binding analyses",
						-foreground => 'lightgreen',
                        -value=>"three",
						-variable=>\$testType
						);
	 
	 

# Buttons

my $pipeButton = $mw -> Button(-text => "RUN DROIDS", 
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

$oneRadio->pack(-anchor=>"w");
$twoRadio->pack(-anchor=>"w");
$threeRadio->pack(-anchor=>"w");
#$droidsButton1->pack(-side=>"right",
#			-anchor=>"n");
#$droidsButton2->pack(-side=>"left",
#			-anchor=>"n");

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
  
     if ($testType eq "one" ) {system "perl DROIDS1.pl\n";}  
	elsif ($testType eq "two" ) {system "perl DROIDS2.pl\n";}
    elsif ($testType eq "three" ) {system "perl DROIDS3.pl\n";}  
	else {print " PLEASE SELECT OPTIONS\n"}
      
  
  
  
  
}

########################################################################
