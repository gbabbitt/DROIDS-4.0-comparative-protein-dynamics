#!/usr/bin/perl
#use Tk;
#use Tk::PNG;
#use Tk::JPEG;
#use Tk::Photo;
#use strict;
#use warnings;
#use feature ":5.10";
#use File::Copy;
#use File::Path;
#use List::Util qw( min );
#use List::Util qw( max );
#use List::Util qw(min max);
use Statistics::Descriptive::Full();

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
#print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
#                         in Dynamic Simulations
#
#- visual toolbox for functional evolutionary comparison
#  of molecular dynamic simulation \n\n";

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
      if ($header eq "shape"){$vector_enter = $value;}
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

#sleep(1);print "MASK LEARNING TO ONLY SIGNIFICANTLY DIFFERENT DYNAMICS? (type 'on' or 'off')\n\n";
#my $option = <STDIN>;
#chop($option);
#sleep(1);print "SELECT MOVIE VIEWER TYPE (type 'fix' or 'roll')\n\n";
#my $view = <STDIN>;
#chop($view);

#sub ctl {

##############################
#print("appending ctl file...\n");
#open(CTL, '>>', "DROIDS.ctl") or die "Could not open output file";
#print CTL "MLoption\t"."$option\t # mask method (apply to signif KS result...on or off)\n";
#print CTL "MLview\t"."$view\t # movie viewing type\n";
#close CTL;
#print("CTL file made\n");
#sleep(1);
##############################
mkdir("indAAtrain");
#mkdir("./trainingData_$queryID/indAA");
print("\nparsing training set time series for fluctuation data for each amino acid in $refID...\n");
sleep(1);
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">"."./indAAtrain/fluxtimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./trainingData_$refID/fluxtime_$refID"."_$tt".".txt")||die "could not open time series file "."fluxtime_$refID"."_$tt.txt\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $DATAseries = substr($INrow, 28, length($INrow) - 28);  # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	$aaID = int($a/4 - 0.1);
     $atomID = @INrow[1];
	$overallAVG = @INrow[2];
     if ($a == 0){  # create headers for file
          @DATAseries = split(/\s+/, $DATAseries);
          $framelistlength = scalar @DATAseries;
          $framelist = '';
          $frame = 0;
          for (my $aa = 0; $aa < $framelistlength; $aa++){$frame = $frame +1; $framelist = "$framelist"."$frame\t";}
          if ($tt == 0){print OUT "class\t"."$framelist\n";}
          next;}
     if ($aaID == $t){print OUT "0\t".$DATAseries;} #print "$tt\t"."$aaID\t"."$atomID\t"."$overallAVG\n";}
     }
   close IN;
  }
   close OUT;
   
 }

print("\nparsing training set time series for fluctuation data for each amino acid in $queryID...\n");
sleep(1);
 for (my $t = 0; $t < $lengthID; $t++){
   open(OUT, ">>"."./indAAtrain/fluxtimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
   for (my $tt = 0; $tt < $number_runs; $tt++){
    open(IN, "<"."./trainingData_$queryID/fluxtime_$queryID"."_$tt".".txt")||die "could not open time series file "."fluxtime_$queryID"."_$tt.txt\n";
    my @IN = <IN>;
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     $DATAseries = substr($INrow, 28, length($INrow) - 28); # to remove atomID and overall avg
     @INrow = split(/\s+/, $INrow);
	$aaID = int($a/4 - 0.1);
     $atomID = @INrow[1];
	$overallAVG = @INrow[2];
     if($atomID eq "AtomicFlx"){next;} # skip header row if encountered
     if ($aaID == $t){print OUT "1\t".$DATAseries;} #print "$tt\t"."$aaID\t"."$atomID\t"."$overallAVG\n";}
     }
   close IN;
  }
  close OUT; 
 }

###########################################
###########################################

if ($vector_enter eq 'y'){

mkdir("indAAtrain_corr");

print("\nparsing training set time series for correlation data for each amino acid in $refID...\n");
sleep(1);
 for (my $t = 0; $t < $lengthID; $t++){
  #if($t == 0){next;}
  $tt = $t+1;
  open(OUT, ">"."./indAAtrain_corr/corrtimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
  #print OUT "corr1\t"."corr3\t"."corr5\t"."corr9\n";
  $corrlabel1 = $framegroups + 1;
  $corrlabel3 = $framegroups + 2;
  $corrlabel5 = $framegroups + 3;
  $corrlabel9 = $framegroups + 4;
  print OUT "$corrlabel1\t"."$corrlabel3\t"."$corrlabel5\t"."$corrlabel9\n";
  open(IN, "<"."./atomcorr/DROIDScorrelation"."_$tt".".txt")||die "could not open corr file "."DROIDScorrelation"."_$t".".txt\n";
    my @IN = <IN>;   
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     @INrow = split(/\s+/, $INrow);
	$sample = @INrow[0];
     $ref_corr1 = @INrow[8]; $ref_corr3 = @INrow[10]; $ref_corr5 = @INrow[12]; $ref_corr9 = @INrow[14];
     if ($sample =~ /\d/) {print OUT "$ref_corr1\t"."$ref_corr3\t"."$ref_corr5\t"."$ref_corr9\n";}
     }
     close IN;
     close OUT; 
     }

print("\nparsing training set time series for correlation data for each amino acid in $queryID...\n");
sleep(1); 
 
for (my $t = 0; $t < $lengthID; $t++){
  #if($t == 0){next;}
  $tt = $t+1;
  open(OUT, ">>"."./indAAtrain_corr/corrtimeAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
  open(IN, "<"."./atomcorr/DROIDScorrelation"."_$tt".".txt")||die "could not open corr file "."DROIDScorrelation"."_$t".".txt\n";
    my @IN = <IN>;   
    for (my $a = 0; $a < scalar @IN; $a++) {
	$INrow = $IN[$a];
     @INrow = split(/\s+/, $INrow);
	$sample = @INrow[0];
     $query_corr1 = @INrow[9]; $query_corr3 = @INrow[11]; $query_corr5 = @INrow[13]; $query_corr9 = @INrow[15];
     if ($sample =~ /\d/) {print OUT "$query_corr1\t"."$query_corr3\t"."$query_corr5\t"."$query_corr9\n";}
     }
     close IN;
     close OUT; 
     }

     
mkdir("indAAtrain_flux+corr");
print "\n\ncombining atom fluctuation and atom correlation information\n\n";
sleep(1);

for (my $t = 0; $t < $lengthID; $t++){
  open(OUT, ">"."./indAAtrain_flux+corr/fluxcorrAA_$refID"."_$t".".txt")||die "could not create AA time series file\n";
  open(IN1, "<"."./indAAtrain/fluxtimeAA_$refID"."_$t".".txt")||die "could not open flux file "."fluxtimeAA_$refID"."_$t".".txt\n";
  open(IN2, "<"."./indAAtrain_corr/corrtimeAA_$refID"."_$t".".txt")||die "could not open corr file "."corrtimeAA_$refID"."_$t".".txt\n";
  my @IN1 = <IN1>;
  my @IN2 = <IN2>;
  for (my $a = 0; $a < scalar @IN1; $a++) {
	 $IN1row = $IN1[$a];
      $IN2row = $IN2[$a];
      chomp $IN1row;
      chomp $IN2row;
      if($IN2row == ''){$IN2row = "0.5\t"."0.5\t"."0.5\t"."0.5\t";} # avg if array is empty
      if($a == 0){print OUT "$IN1row"."$IN2row\n";}
      if($a != 0){print OUT "$IN1row\t"."$IN2row\n";}
      }
      close IN1;
      close IN2;
      close OUT;
      }
} # end if statement
############################################
print("\nparsing is done\n");
sleep(1);
exit;

###########################################################################################################
###########################################################################################################
