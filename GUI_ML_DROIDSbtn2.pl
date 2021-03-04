#!/usr/bin/perl
#use Tk;
#use Tk::PNG;
#use Tk::JPEG;
#use Tk::Photo;
#use strict;
#use warnings;
#use feature ":5.10";
use File::Copy;
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

## training set
mkdir("indAAtrain_flux");
print("\nparsing training set time series for fluctuation data for each amino acid in $refID...\n");
sleep(1);
move("./trainingDataFLUX_$refID","./indAAtrain_flux/trainingDataFLUX_$refID"); 
print("\nparsing training set time series for fluctuation data for each amino acid in $queryID...\n");
sleep(1);
move("./trainingDataFLUX_$queryID","./indAAtrain_flux/trainingDataFLUX_$queryID"); 
## testing set
mkdir("indAAtest_flux");
for (my $g = 0; $g < 2; $g++){
if ($g == 0) {open(MUT, "<"."copy_list.txt");}
if ($g == 1) {open(MUT, "<"."variant_list.txt");}
my @MUT = <MUT>;
#print @MUT;
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      if ($g == 1 && $p == 1){next;}
      if ($g == 1 && $p == 2){next;}
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      print("\nparsing testing set time series for fluctuation data for each amino acid in $fileIDq...\n");
      sleep(1);
      move("./testingDataFLUX_$fileIDq","./indAAtest_flux/testingDataFLUX_$fileIDq"); 
close MUT;
} #end for loop
} # end outer for loop

####################################
if ($vector_enter eq 'y'){
## training set
mkdir("indAAtrain_corr");
mkdir("trash");
print("\nparsing training set time series for correlation data for each amino acid in $refID...\n");
sleep(1);
move("./trainingDataCORR_$refID"."_1res","./indAAtrain_corr/trainingDataCORR_$refID"."_1res"); 
move("./trainingDataCORR_$refID"."_3res","./indAAtrain_corr/trainingDataCORR_$refID"."_3res");
move("./trainingDataCORR_$refID"."_5res","./indAAtrain_corr/trainingDataCORR_$refID"."_5res");
move("./trainingDataCORR_$refID"."_9res","./indAAtrain_corr/trainingDataCORR_$refID"."_9res");
move("./trainingDataCORR_$refID"."_ALLPAIRS","./trash/trainingDataCORR_$refID"."_ALLPAIRS");
print("\nparsing training set time series for correlation data for each amino acid in $queryID...\n");
sleep(1); 
move("./trainingDataCORR_$queryID"."_1res","./indAAtrain_corr/trainingDataCORR_$queryID"."_1res"); 
move("./trainingDataCORR_$queryID"."_3res","./indAAtrain_corr/trainingDataCORR_$queryID"."_3res");
move("./trainingDataCORR_$queryID"."_5res","./indAAtrain_corr/trainingDataCORR_$queryID"."_5res");
move("./trainingDataCORR_$queryID"."_9res","./indAAtrain_corr/trainingDataCORR_$queryID"."_9res");
move("./trainingDataCORR_$queryID"."_ALLPAIRS","./trash/trainingDataCORR_$queryID"."_ALLPAIRS");
## testing set
mkdir("indAAtest_corr");
for (my $g = 0; $g < 2; $g++){
if ($g == 0) {open(MUT, "<"."copy_list.txt");}
if ($g == 1) {open(MUT, "<"."variant_list.txt");}
my @MUT = <MUT>;
#print @MUT;
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      if ($g == 1 && $p == 1){next;}
      if ($g == 1 && $p == 2){next;}
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      print("\nparsing testing set time series for correlation data for each amino acid in $fileIDq...\n");
      sleep(1);
      move("./testingDataCORR_$fileIDq"."_1res","./indAAtest_corr/testingDataCORR_$fileIDq"."_1res"); 
      move("./testingDataCORR_$fileIDq"."_3res","./indAAtest_corr/testingDataCORR_$fileIDq"."_3res");
      move("./testingDataCORR_$fileIDq"."_5res","./indAAtest_corr/testingDataCORR_$fileIDq"."_5res");
      move("./testingDataCORR_$fileIDq"."_9res","./indAAtest_corr/testingDataCORR_$fileIDq"."_9res");
      move("./testingDataCORR_$fileIDq"."_ALLPAIRS","./trash/testingDataCORR_$fileIDq"."_ALLPAIRS"); 
close MUT;
} #end for loop
} # end outer for loop
     
} # end if statement
############################################

# create final training data sets
 print("\nbuilding final training/testing data sets for machine learning\n");
      sleep(1);

#####################################################################
#####################################################################
if ($vector_enter ne 'y'){ # build data sets with local fluctuations only (0, 1, 3, 5, 9, sites downstream)
mkdir("trainingData");     

# collecting reference training data
for (my $i = 0; $i <= $lengthID; $i++){
print "collecting position $i for $refID\n";
open (OUT, ">"."./trainingData/"."position$i") || die " could not open output file\n";
print OUT "class\t"."flux0\t"."flux0\t"."flux3\t"."flux5\t"."flux9\n";
$class = 0;
$dir = "./indAAtrain_flux/trainingDataFLUX_$refID"."/";
opendir(DIR,$dir)or die "can't open directory $dir:$!";
while ($filename = readdir DIR){ # loop through files
    $filecount = $filecount + 1;
    $filelocation = "$dir"."$filename";
    open(IN, $filelocation) or die "Cannot open file";
    my @IN = <IN>;
    for (my $j = 0; $j < scalar @IN; $j++){
	 my $INrow = $IN[$j];
      my $INrow1 = $IN[$j+1];
      my $INrow3 = $IN[$j+3];
      my $INrow5 = $IN[$j+5];
      my $INrow9 = $IN[$j+9];
	 my @INrow = split (/\s+/, $INrow);
      my @INrow1 = split (/\s+/, $INrow1);
      my @INrow3 = split (/\s+/, $INrow3);
      my @INrow5 = split (/\s+/, $INrow5);
      my @INrow9 = split (/\s+/, $INrow9);
        $pos = $INrow[1];
        $flux = $INrow[2];
        $pos1 = $INrow1[1];
        $flux1 = $INrow1[2];
        $pos3 = $INrow3[1];
        $flux3 = $INrow3[2];
        $pos5 = $INrow5[1];
        $flux5 = $INrow5[2];
        $pos9 = $INrow9[1];
        $flux9 = $INrow9[2];
								if($flux1 eq ''){$flux1 = $flux;}
								if($flux3 eq ''){$flux3 = $flux;}
								if($flux5 eq ''){$flux5 = $flux;}
								if($flux9 eq ''){$flux9 = $flux;}
								#print "pos=$pos\t"."flux=$flux\t"."pos1=$pos1\t"."flux1=$flux1\t"."pos3=$pos3\t"."flux3=$flux3\t"."pos5=$pos5\t"."flux5=$flux5\t"."pos9=$pos9\t"."flux9=$flux9\n";
        if($pos == $i){print OUT "$class\t"."$flux\t"."$flux1\t"."$flux3\t"."$flux5\t"."$flux9\n";}
    } # end for
    close IN;
} # end while
close DIR;
close OUT;
} # end for
# collecting query training data
for (my $i = 0; $i <= $lengthID; $i++){
print "collecting position $i for $queryID\n";
open (OUT, ">>"."./trainingData/"."position$i") || die " could not open output file\n";
$class = 1;
$dir = "./indAAtrain_flux/trainingDataFLUX_$queryID"."/";
opendir(DIR,$dir)or die "can't open directory $dir:$!";
while ($filename = readdir DIR){ # loop through files
    $filecount = $filecount + 1;
    $filelocation = "$dir"."$filename";
    open(IN, $filelocation) or die "Cannot open file";
    my @IN = <IN>;
    for (my $j = 0; $j < scalar @IN; $j++){
	 my $INrow = $IN[$j];
      my $INrow1 = $IN[$j+1];
      my $INrow3 = $IN[$j+3];
      my $INrow5 = $IN[$j+5];
      my $INrow9 = $IN[$j+9];
	 my @INrow = split (/\s+/, $INrow);
      my @INrow1 = split (/\s+/, $INrow1);
      my @INrow3 = split (/\s+/, $INrow3);
      my @INrow5 = split (/\s+/, $INrow5);
      my @INrow9 = split (/\s+/, $INrow9);
        $pos = $INrow[1];
        $flux = $INrow[2];
        $pos1 = $INrow1[1];
        $flux1 = $INrow1[2];
        $pos3 = $INrow3[1];
        $flux3 = $INrow3[2];
        $pos5 = $INrow5[1];
        $flux5 = $INrow5[2];
        $pos9 = $INrow9[1];
        $flux9 = $INrow9[2];
								if($flux1 eq ''){$flux1 = $flux;}
								if($flux3 eq ''){$flux3 = $flux;}
								if($flux5 eq ''){$flux5 = $flux;}
								if($flux9 eq ''){$flux9 = $flux;}
        #print "pos=$pos\t"."flux=$flux\t"."pos1=$pos1\t"."flux1=$flux1\t"."pos3=$pos3\t"."flux3=$flux3\t"."pos5=$pos5\t"."flux5=$flux5\t"."pos9=$pos9\t"."flux9=$flux9\n";
        if($pos == $i){print OUT "$class\t"."$flux\t"."$flux1\t"."$flux3\t"."$flux5\t"."$flux9\n";}
    } # end for
    close IN;
} # end while
close DIR;
close OUT;
} # end for

# collecting testing data
for (my $g = 0; $g < 2; $g++){
if ($g == 0) {open(MUT, "<"."copy_list.txt");}
if ($g == 1) {open(MUT, "<"."variant_list.txt");}
my @MUT = <MUT>;
#print @MUT;
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      if ($g == 1 && $p == 1){next;}
      if ($g == 1 && $p == 2){next;}
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      mkdir("testingData_$fileIDq");
      for (my $i = 0; $i <= $lengthID; $i++){
      print "collecting position $i for $fileIDq\n";
      open (OUT, ">"."./testingData_$fileIDq/"."position$i") || die " could not open output file\n";
      print OUT "class\t"."flux0\t"."flux0\t"."flux3\t"."flux5\t"."flux9\n";
      $class = "unk";
      $dir = "./indAAtest_flux/testingDataFLUX_$fileIDq"."/";
      opendir(DIR,$dir)or die "can't open directory $dir:$!";
      while ($filename = readdir DIR){ # loop through files
          $filecount = $filecount + 1;
          $filelocation = "$dir"."$filename";
          open(IN, $filelocation) or die "Cannot open file";
          my @IN = <IN>;
          for (my $j = 0; $j < scalar @IN; $j++){
	       my $INrow = $IN[$j];
            my $INrow1 = $IN[$j+1];
            my $INrow3 = $IN[$j+3];
            my $INrow5 = $IN[$j+5];
            my $INrow9 = $IN[$j+9];
	       my @INrow = split (/\s+/, $INrow);
            my @INrow1 = split (/\s+/, $INrow1);
            my @INrow3 = split (/\s+/, $INrow3);
            my @INrow5 = split (/\s+/, $INrow5);
            my @INrow9 = split (/\s+/, $INrow9);
              $pos = $INrow[1];
              $flux = $INrow[2];
              $pos1 = $INrow1[1];
              $flux1 = $INrow1[2];
              $pos3 = $INrow3[1];
              $flux3 = $INrow3[2];
              $pos5 = $INrow5[1];
              $flux5 = $INrow5[2];
              $pos9 = $INrow9[1];
              $flux9 = $INrow9[2];
														if($flux1 eq ''){$flux1 = $flux;}
														if($flux3 eq ''){$flux3 = $flux;}
														if($flux5 eq ''){$flux5 = $flux;}
														if($flux9 eq ''){$flux9 = $flux;}
              #print "pos=$pos\t"."flux=$flux\t"."pos1=$pos1\t"."flux1=$flux1\t"."pos3=$pos3\t"."flux3=$flux3\t"."pos5=$pos5\t"."flux5=$flux5\t"."pos9=$pos9\t"."flux9=$flux9\n";
              if($pos == $i){print OUT "$class\t"."$flux\t"."$flux1\t"."$flux3\t"."$flux5\t"."$flux9\n";}
              } # end for
          close IN;
      } # end while
      close DIR;
      close OUT;
      } # end for
        
close MUT;
} #end for loop
} # end outer for loop

} # end
################################################################
################################################################
if ($vector_enter eq 'y'){  # build data sets with site-specific fluctuations and local correlations (at 1,3,5 and 9 AA sites downstream)
mkdir("trainingData");     

# collecting reference training data
for (my $i = 0; $i <= $lengthID; $i++){
print "collecting position $i for $refID\n";
open (OUT, ">"."./trainingData/"."position$i") || die " could not open output file\n";
print OUT "class\t"."flux0\t"."corr0\t"."corr3\t"."corr5\t"."corr9\n";
$class = 0;
$dir = "./indAAtrain_flux/trainingDataFLUX_$refID"."/";
$dir1 = "./indAAtrain_corr/trainingDataCORR_$refID"."_1res/";
$dir3 = "./indAAtrain_corr/trainingDataCORR_$refID"."_3res/";
$dir5 = "./indAAtrain_corr/trainingDataCORR_$refID"."_5res/";
$dir9 = "./indAAtrain_corr/trainingDataCORR_$refID"."_9res/";
opendir(DIR,$dir)or die "can't open directory $dir0:$!";
while ($filename = readdir DIR){ # loop through files
    $filecount = $filecount + 1;
    $filelocation = "$dir"."$filename";
    open(IN, $filelocation) or die "Cannot open file";
    my @IN = <IN>;
    # find matching corr file names to flux file name
    $subname = substr($filename,4);
    #print "$subname\n";
    $filename1 = "corr".$subname;
    $filename3 = "corr".$subname;
    $filename5 = "corr".$subname;
    $filename9 = "corr".$subname;
    #print "$filename1\t"."$filename3\t"."$filename5\t"."$filename9\n";
    $filelocation1 = "$dir1"."$filename1";
    open(IN1, $filelocation1);
    my @IN1 = <IN1>;
    $filelocation3 = "$dir3"."$filename3";
    open(IN3, $filelocation3);
    my @IN3 = <IN3>;
    $filelocation5 = "$dir5"."$filename5";
    open(IN5, $filelocation5);
    my @IN5 = <IN5>;
    $filelocation9 = "$dir9"."$filename9";
    open(IN9, $filelocation9);
    my @IN9 = <IN9>;
    
    for (my $j = 0; $j < scalar @IN; $j++){
	 my $INrow = $IN[$j];
      my $INrow1 = $IN1[$j];
      my $INrow3 = $IN3[$j];
      my $INrow5 = $IN5[$j];
      my $INrow9 = $IN9[$j];
	 my @INrow = split (/\s+/, $INrow);
      my @INrow1 = split (/\s+/, $INrow1);
      my @INrow3 = split (/\s+/, $INrow3);
      my @INrow5 = split (/\s+/, $INrow5);
      my @INrow9 = split (/\s+/, $INrow9);
        $pos = $INrow[1];
        $flux = $INrow[2];
        $pos1 = $INrow1[1];
        $corr1 = $INrow1[2];
        $pos3 = $INrow3[1];
        $corr3 = $INrow3[2];
        $pos5 = $INrow5[1];
        $corr5 = $INrow5[2];
        $pos9 = $INrow9[1];
        $corr9 = $INrow9[2];
        #print "pos=$pos\t"."flux=$flux\t"."pos1=$pos1\t"."corr1=$corr1\t"."pos3=$pos3\t"."corr3=$corr3\t"."pos5=$pos5\t"."corr5=$corr5\t"."pos9=$pos9\t"."corr9=$corr9\n";
        if($pos == $i){print OUT "$class\t"."$flux\t"."$corr1\t"."$corr3\t"."$corr5\t"."$corr9\n";}
        
    } # end for
    close IN;
    close IN1;
    close IN3;
    close IN5;
    close IN9;
} # end while
close DIR;
close OUT;
} # end for

# collecting query training data
for (my $i = 0; $i <= $lengthID; $i++){
print "collecting position $i for $queryID\n";
open (OUT, ">>"."./trainingData/"."position$i") || die " could not open output file\n";
$class = 1;
$dir = "./indAAtrain_flux/trainingDataFLUX_$queryID"."/";
$dir1 = "./indAAtrain_corr/trainingDataCORR_$queryID"."_1res/";
$dir3 = "./indAAtrain_corr/trainingDataCORR_$queryID"."_3res/";
$dir5 = "./indAAtrain_corr/trainingDataCORR_$queryID"."_5res/";
$dir9 = "./indAAtrain_corr/trainingDataCORR_$queryID"."_9res/";
opendir(DIR,$dir)or die "can't open directory $dir0:$!";
while ($filename = readdir DIR){ # loop through files
    $filecount = $filecount + 1;
    $filelocation = "$dir"."$filename";
    open(IN, $filelocation) or die "Cannot open file";
    my @IN = <IN>;
    # find matching corr file names to flux file name
    $subname = substr($filename,4);
    #print "$subname\n";
    $filename1 = "corr".$subname;
    $filename3 = "corr".$subname;
    $filename5 = "corr".$subname;
    $filename9 = "corr".$subname;
    #print "$filename1\t"."$filename3\t"."$filename5\t"."$filename9\n";
    $filelocation1 = "$dir1"."$filename1";
    open(IN1, $filelocation1);
    my @IN1 = <IN1>;
    $filelocation3 = "$dir3"."$filename3";
    open(IN3, $filelocation3);
    my @IN3 = <IN3>;
    $filelocation5 = "$dir5"."$filename5";
    open(IN5, $filelocation5);
    my @IN5 = <IN5>;
    $filelocation9 = "$dir9"."$filename9";
    open(IN9, $filelocation9);
    my @IN9 = <IN9>;
    
    for (my $j = 0; $j < scalar @IN; $j++){
	 my $INrow = $IN[$j];
      my $INrow1 = $IN1[$j];
      my $INrow3 = $IN3[$j];
      my $INrow5 = $IN5[$j];
      my $INrow9 = $IN9[$j];
	 my @INrow = split (/\s+/, $INrow);
      my @INrow1 = split (/\s+/, $INrow1);
      my @INrow3 = split (/\s+/, $INrow3);
      my @INrow5 = split (/\s+/, $INrow5);
      my @INrow9 = split (/\s+/, $INrow9);
        $pos = $INrow[1];
        $flux = $INrow[2];
        $pos1 = $INrow1[1];
        $corr1 = $INrow1[2];
        $pos3 = $INrow3[1];
        $corr3 = $INrow3[2];
        $pos5 = $INrow5[1];
        $corr5 = $INrow5[2];
        $pos9 = $INrow9[1];
        $corr9 = $INrow9[2];
        #print "pos=$pos\t"."flux=$flux\t"."pos1=$pos1\t"."corr1=$corr1\t"."pos3=$pos3\t"."corr3=$corr3\t"."pos5=$pos5\t"."corr5=$corr5\t"."pos9=$pos9\t"."corr9=$corr9\n";
        if($pos == $i){print OUT "$class\t"."$flux\t"."$corr1\t"."$corr3\t"."$corr5\t"."$corr9\n";}
        
    } # end for
    close IN;
    close IN1;
    close IN3;
    close IN5;
    close IN9;
} # end while
close DIR;
close OUT;
} # end for


# collecting testing data
for (my $g = 0; $g < 2; $g++){
if ($g == 0) {open(MUT, "<"."copy_list.txt");}
if ($g == 1) {open(MUT, "<"."variant_list.txt");}
my @MUT = <MUT>;
#print @MUT;
for (my $p = 0; $p < scalar @MUT; $p++){
	 if ($p == 0){next;}
      my $MUTrow = $MUT[$p];
      if ($g == 1 && $p == 1){next;}
      if ($g == 1 && $p == 2){next;}
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      mkdir("testingData_$fileIDq");
      for (my $i = 0; $i <= $lengthID; $i++){
      print "collecting position $i for $fileIDq\n";
      open (OUT, ">"."./testingData_$fileIDq/"."position$i") || die " could not open output file\n";
      print OUT "class\t"."flux0\t"."corr0\t"."corr3\t"."corr5\t"."corr9\n";
      $class = "unk";
      $dir = "./indAAtest_flux/testingDataFLUX_$fileIDq"."/";
      $dir1 = "./indAAtest_corr/testingDataCORR_$fileIDq"."_1res/";
      $dir3 = "./indAAtest_corr/testingDataCORR_$fileIDq"."_3res/";
      $dir5 = "./indAAtest_corr/testingDataCORR_$fileIDq"."_5res/";
      $dir9 = "./indAAtest_corr/testingDataCORR_$fileIDq"."_9res/";
      opendir(DIR,$dir)or die "can't open directory $dir0:$!";
      while ($filename = readdir DIR){ # loop through files
          $filecount = $filecount + 1;
          $filelocation = "$dir"."$filename";
          open(IN, $filelocation) or die "Cannot open file";
          my @IN = <IN>;
          # find matching corr file names to flux file name
          $subname = substr($filename,4);
          #print "$subname\n";
          $filename1 = "corr".$subname;
          $filename3 = "corr".$subname;
          $filename5 = "corr".$subname;
          $filename9 = "corr".$subname;
          #print "$filename1\t"."$filename3\t"."$filename5\t"."$filename9\n";
          $filelocation1 = "$dir1"."$filename1";
          open(IN1, $filelocation1);
          my @IN1 = <IN1>;
          $filelocation3 = "$dir3"."$filename3";
          open(IN3, $filelocation3);
          my @IN3 = <IN3>;
          $filelocation5 = "$dir5"."$filename5";
          open(IN5, $filelocation5);
          my @IN5 = <IN5>;
          $filelocation9 = "$dir9"."$filename9";
          open(IN9, $filelocation9);
          my @IN9 = <IN9>;
    
         for (my $j = 0; $j < scalar @IN; $j++){
        	 my $INrow = $IN[$j];
           my $INrow1 = $IN1[$j];
           my $INrow3 = $IN3[$j];
           my $INrow5 = $IN5[$j];
           my $INrow9 = $IN9[$j];
     	 my @INrow = split (/\s+/, $INrow);
           my @INrow1 = split (/\s+/, $INrow1);
           my @INrow3 = split (/\s+/, $INrow3);
           my @INrow5 = split (/\s+/, $INrow5);
           my @INrow9 = split (/\s+/, $INrow9);
             $pos = $INrow[1];
             $flux = $INrow[2];
             $pos1 = $INrow1[1];
             $corr1 = $INrow1[2];
             $pos3 = $INrow3[1];
             $corr3 = $INrow3[2];
             $pos5 = $INrow5[1];
             $corr5 = $INrow5[2];
             $pos9 = $INrow9[1];
             $corr9 = $INrow9[2];
             #print "pos=$pos\t"."flux=$flux\t"."pos1=$pos1\t"."corr1=$corr1\t"."pos3=$pos3\t"."corr3=$corr3\t"."pos5=$pos5\t"."corr5=$corr5\t"."pos9=$pos9\t"."corr9=$corr9\n";
             if($pos == $i){print OUT "$class\t"."$flux\t"."$corr1\t"."$corr3\t"."$corr5\t"."$corr9\n";}
        
          } # end for
          close IN;
          close IN1;
          close IN3;
          close IN5;
          close IN9;
     } # end while
     close DIR;
     close OUT;
     } # end for

close MUT;
} #end for loop
} # end outer for loop

} # end 



############################################


print("\nparsing is done\n");
sleep(1);
exit;

###########################################################################################################
###########################################################################################################
