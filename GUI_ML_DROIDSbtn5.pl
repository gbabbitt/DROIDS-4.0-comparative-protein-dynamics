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
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variants, $fileIDq);
      } #end for loop
close MUT;
print "\n variants are @variants\n\n";
sleep(1);


# create array of names of copies and variants
open(MUT, "<"."variant_label_list.txt");
my @MUT = <MUT>;
#print @MUT;
@variant_labels = ();
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variant_labels, $fileIDq);
      } #end for loop
close MUT;
print "\n variant labels are @variant_labels\n\n";
sleep(1);

print "Enter number position of N terminal on this chain\n";
print "(default = 0)\n\n";
my $startN = <STDIN>;
chop($startN);
if ($startN eq ''){$startN = 0;}

print "Choose color palette for matrix (type 'heat' or 'wesanderson')\n";
print "NOTE: zissou1 color palette requires install.packages('wesanderson') in R \n";
print "(default = 'heat')\n\n";
my $pal = <STDIN>;
chop($pal);
if ($pal eq ''){$pal = 'heat';}

################################################################
################################################################
#print "\n THIS FEATURE IS UNDER RECONSTRUCTION....try again later\n\n"; exit;

# prompt user - selest machine learner for MI matrix
sleep(1);print "CHOOSE BEST LEARNER FOR MUTUAL INFORMATION MATRIX (type 'KNN', 'NB', 'LDA', 'QDA', 'SVM', 'RFOR' or 'ADA')\n\n";
my $learner = <STDIN>;
chop($learner);

for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];

# initialize     
$residue1 = 0;
$residue2 = 0;
$MI = 0;
# open MI matrix files
open(MI, ">"."./testingData_$variantID/MIcolumn_$learner"."_$variantID.txt") or die "could not open MI matrix file for $variantID \n";    
print MI "pos1\t"."pos2\t"."MI\n";
open(MI2, ">"."./testingData_$variantID/MImatrix_$learner"."_$variantID.txt") or die "could not open MI matrix file for $variantID \n";    
#print MI2 "residue\t";
#for (my $n = 0; $n < $lengthID; $n++){print MI2 "$n\t";}
#print MI2 "\n";

#  create vertical files of 0,1 for analysis
print "\n reading input files \n";     
sleep(1);

for (my $g = 1; $g < $lengthID; $g++){
open(POS1, "<"."./testingData_$variantID/indAAclass$learner/position_$g.txt") or die "could not open position_$g.txt \n";
	 open(VPOS1, ">"."./testingData_$variantID/indAAclass$learner/vertposition_$g.txt") or die "could not open vertposition_$g.txt \n";
	   my @POS1 = <POS1>;
	   # collect zeros and ones and print vertically
	   for (my $z = 0; $z < scalar @POS1; $z++){
              my $POS1row = $POS1[$z];
              my @POS1row = split (/\s+/, $POS1row);
	         $head = $POS1row[0];
		    $value = $POS1row;
		    #print "reading i\t"."$head\t"."$value\n";
	         if ($head ne "[1]" || $head ne "Levels:"){
			    $ZerosOnes = $value;
			    #print "ZerosOnes\t"."$ZerosOnes\n";
			    @ZerosOnes = split(/\s+/, $ZerosOnes);
			    for (my $zz = 0; $zz < scalar @ZerosOnes; $zz++){
				  $classvalue = $ZerosOnes[$zz];
				  if($classvalue eq "0" || $classvalue eq "1"){print VPOS1 "$classvalue\n";}
				  #print "classvalue\t"."$classvalue\n";
			       }
			    }
		    }
	 close POS1;
	 close VPOS1;
}
sleep(2);

print "\ncalculating mutual information matrix for $variantID \n";     
sleep(1);
for (my $i = 1; $i < $lengthID; $i++){
	 
      print "calculating MI values for residue\t"."$residue1 for $variantID\n";
      open(POS1, "<"."./testingData_$variantID/indAAclass$learner/vertposition_$i.txt") or die "could not open MI matrix file for ./testingData_$variantID/indAAclass$learner/position_$i.txt \n";    
      my @POS1 = <POS1>;
	 $residue1 = $i;
      if($i>0){print MI2 "\n";}
      #if($i==0){print MI2 "$residue1\t";}
      #if($i>0){print MI2 "\n"; print MI2 "$residue1\t";}
      for (my $j = 1; $j < $lengthID; $j++){
	     
	     open(POS2, "<"."./testingData_$variantID/indAAclass$learner/vertposition_$j.txt") or die "could not open MI matrix file for ./testingData_$variantID/indAAclass$learner/position_$j.txt \n";    
          my @POS2 = <POS2>;
          $residue2 = $j;
          
          #calculate freq A
          #print "calculating freq state A\t";
          $freqA = 0;
          $sumclassA = 0;
          $cntA = 0;
          for (my $ii = 0; $ii < scalar @POS1; $ii++){
              my $POS1row = $POS1[$ii];
              my @POS1row = split (/\s+/, $POS1row);
	         $class1 = $POS1row[0];
              $time1 = $ii;
              if($class1 == 0 || $class1 == 0.5 || $class1 == 1) {$sumclassA = $sumclassA+$class1; $cntA = $cntA+1;}
              #print "Residue1\t"."$residue1\t"."timeslice\t"."$time1\t"."class = "."$class1\n";
              }
              $freqA = $sumclassA/($cntA+0.0001);
              #print "freqA = "."$freqA\n";
          
          #calculate freq B
          #print "calculating freq state B\t";
          $freqB = 0;
          $sumclassB = 0;
          $cntB = 0;
          for (my $jj = 0; $jj < scalar @POS2; $jj++){
              my $POS2row = $POS2[$jj];
              my @POS2row = split (/\s+/, $POS2row);
	         $class2 = $POS2row[0];
              $time2 = $jj;
              if($class2 == 0 || $class2 == 0.5 || $class2 == 1) {$sumclassB = $sumclassB+$class2; $cntB = $cntB+1;}
              #print "Residue2\t"."$residue2\t"."timeslice\t"."$time2\t"."class = "."$class2\n";
              }
              $freqB = $sumclassB/($cntB+0.0001);
              #print "freqB = "."$freqB\n";
          
          #calculate freq A and B
          #print "calculating freq state A and B\t";
          $freqAB = 0;
          $sumclassAB = 0;
          $cntAB = 0;
          for (my $ii = 0; $ii < scalar @POS1; $ii++){
              my $POS1row = $POS1[$ii];
              my @POS1row = split (/\s+/, $POS1row);
	         $class1 = $POS1row[0];
              $time1 = $ii;
              #if($class1 == 0 || $class1 == 0.5 || $class1 == 1) {$sumclassA = $sumclassA+$class1; $cntA = $cntA+1;}
              #print "Residue1\t"."$residue1\t"."timeslice\t"."$time1\t"."class = "."$class1\n";
              for (my $jj = 0; $jj < scalar @POS2; $jj++){
              my $POS2row = $POS2[$jj];
              my @POS2row = split (/\s+/, $POS2row);
	         $class2 = $POS2row[0];
              $time2 = $jj;
              if($time1 ==$time2 && $class1 == $class2) {$sumclassAB = $sumclassAB+1; $cntAB = $cntAB+1;}
              if($time1 ==$time2 && $class1 != $class2) {$cntAB = $cntAB+1;}
              #print "Residue2\t"."$residue2\t"."timeslice\t"."$time2\t"."class = "."$class2\n";
              }
              }
              $freqAB = $sumclassAB/($cntAB+0.0001);
              #print "freqAB = "."$freqAB\n";
                    
          # calculate MI
          $MI = $freqAB*log(($freqAB+0.0001)/(($freqA*$freqB)+0.0001));
          if ($MI>1){$MI = 1;}
          #print "MImatrix\t"."$residue1\t"."$residue2\t"."$MI\n";
          print MI "$residue1\t"."$residue2\t"."$MI\n";
          print MI2 "$MI\t";
      }
}

close MI;
close MI2;

print "\nheatmapping mutual information matrix for $variantID (close .pdf to continue)\n";     
sleep(1);

open (Rinput, "| R --vanilla")||die "could not start R command line\n";
print Rinput "datamatrixMI = read.table('./testingData_$variantID/MImatrix_$learner"."_$variantID.txt', header = FALSE)\n";
#print Rinput "print(datamatrixMI)\n";
print Rinput "datamatrixMI<-as.matrix(datamatrixMI)\n";
print Rinput "myMImean <- mean(datamatrixMI)\n";
print Rinput "myMImean <- round(myMImean, digits=4)\n";
print Rinput "print(myMImean)\n";
print Rinput "myMIsd <- sd(datamatrixMI)\n";
print Rinput "myMIsd <- round(myMIsd, digits=4)\n";
print Rinput "print(myMIsd)\n";
#print Rinput "print(datamatrixMI)\n";
#print Rinput "datamatrixMI <- scale(datamatrixMI)\n";
#print Rinput "mymap1<-heatmap(datamatrixMI, Colv = 'Rowv', symm = TRUE, keep.dendro = FALSE)\n";
#print Rinput "print(mymap1)\n";
print Rinput "x <- (1:nrow(datamatrixMI))\n";
print Rinput "y <- (1:ncol(datamatrixMI))\n";
if($pal eq 'heat'){
print Rinput "mymap2<-image(x+$startN, y+$startN, datamatrixMI, col = hcl.colors(20, 'heat', rev = TRUE), main = c('MUTUAL INFORMATION for $variantID (red = significant))', paste('mean MI = ', myMImean, 'sd MI = ', myMIsd)), xlab = 'residue position', ylab = 'residue position')\n";
}
if($pal eq 'wesanderson'){
# for wes anderson color palette (install.packages('wesanderson'))
print Rinput "library(wesanderson)\n";
print Rinput "pal <- wes_palette('Zissou1', 100, type = 'continuous')\n";
print Rinput "mymap2<-image(x+$startN, y+$startN, datamatrixMI, col = pal, main = c('MUTUAL INFORMATION for $variantID (red = significant))', paste('mean MI = ', myMImean, 'sd MI = ', myMIsd)), xlab = 'residue position', ylab = 'residue position')\n";
}
print Rinput "print(mymap2)\n";      
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$variantID/MImatrix_$learner"."_$variantID.pdf";
my $newfilename2 = "./maxDemon_results/MImatrix_$learner"."_$variantID.pdf";
copy($oldfilename, $newfilename);
copy($oldfilename, $newfilename2);
my $oldfilename3 = "./testingData_$variantID/MImatrix_$learner"."_$variantID.txt";
my $newfilename3 = "./maxDemon_results/MImatrix_$learner"."_$variantID.txt";
copy($oldfilename3, $newfilename3);
my $oldfilename4 = "./testingData_$variantID/MIcolumn_$learner"."_$variantID.txt";
my $newfilename4 = "./maxDemon_results/MIcolumn_$learner"."_$variantID.txt";
copy($oldfilename4, $newfilename4);
print " mutual information matrix is complete\n\n";
print " close PDF and txt viewer to continue\n\n";


} # end for loop

for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     system "evince ./testingData_$variantID/MImatrix_$learner"."_$variantID.pdf\n";
}


################################################################
###  mapping significantly coordinated regions to PDB structure
################################################################

# select original structure or variant

# mapping

