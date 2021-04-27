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
# initialize ML methods
$method_bnp = 0;
$method_dist = 0;
$method_kern = 0;
$method_ens = 0;

open(IN3, "<"."MLmethods.txt") or die "could not find MLmethods.txt file\n";
my @IN3 = <IN3>;
for (my $i = 0; $i < scalar @IN3; $i++){
	 my $INrow = $IN3[$i];
	 my @INrow = split (/\s+/, $INrow);
	 my $method = @INrow[0];
	 if ($method eq "bnp"){$method_bnp = 1;}
      if ($method eq "dist"){$method_dist = 1;}
      if ($method eq "kern"){$method_kern = 1;}
      if ($method eq "ens"){$method_ens = 1;}
      if ($method eq "no_bnp"){$method_bnp = 0;}
      if ($method eq "no_dist"){$method_dist = 0;}
      if ($method eq "no_kern"){$method_kern = 0;}
      if ($method eq "no_ens"){$method_ens = 0;}
}
close IN3;

print ("\nMLmethods\n\n");
if ($method_bnp == 1){print"KNN activated\n";}
if ($method_bnp == 0){print"KNN deactivated\n";}
if ($method_dist == 1){print"NB/LDA/QDA activated\n";}
if ($method_dist == 0){print"NB/LDA/QDA deactivated\n";}
if ($method_kern == 1){print"SVM activated\n";}
if ($method_kern == 0){print"SVM deactivated\n";}
if ($method_ens == 1){print"random forest/adaboost activated\n";}
if ($method_ens == 0){print"random forest/adaboost deactivated\n";}
print ("\n\n\n");

###########################################################################################################
###########################################################################################################
system "x-terminal-emulator -e htop\n";
sleep(1);
print "\n\nINITIATING MACHINE LEARNING ON ATOM CORRELATION AND ATOM FLUCTUATION\n\n";
sleep(3);

mkdir("maxDemon_results");

sleep(1);print "CHOOSE TRAINING SET RANDOM SUBSAMPLING OPTION (type 'y' or 'n')\n";
print "if y, training sets will be fixed to n=20,000 time slices sampled with replacement\n";
print "note: use this option if machine learning is too slow (esp SVM and adaboost)\n\n";
my $sampletype_enter = <STDIN>;
chop($sampletype_enter);
if($sampletype_enter eq 'y' || $sampletype_enter eq 'Y'){$sampletype = 'y';}
if($sampletype_enter eq 'n' || $sampletype_enter eq 'Y'){$sampletype = 'n';}
if($sampletype_enter eq ''){$sampletype = 'n';}

if ($method_kern == 1 || $method_other == 1){
# prompt user - choose best learning model to display
sleep(1);print "CHOOSE KERNEL TYPE FOR SVM (type 'linear', 'polynomial', or 'RBF')\n\n";
my $kerntype_enter = <STDIN>;
chop($kerntype_enter);
if($kerntype_enter eq 'linear'){$kerntype = 'vanilladot';}
if($kerntype_enter eq 'polynomial'){$kerntype = 'polydot';}
if($kerntype_enter eq 'laplace'){$kerntype = 'laplacedot';}
if($kerntype_enter eq 'RBF'){$kerntype = 'rbfdot';}
elsif($kerntype_enter eq ''){$kerntype = 'vanilladot';}
sleep(1);
}


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
         
# initialize classpositionHISTO data files (zero values) - for R plots

# KNN method
open (OUT, ">"."./maxDemon_results/classpositionHISTOknn_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# NB method
open (OUT, ">"."./maxDemon_results/classpositionHISTOnb_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# LDA method
open (OUT, ">"."./maxDemon_results/classpositionHISTOlda_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# QDA method
open (OUT, ">"."./maxDemon_results/classpositionHISTOqda_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# SVM method
open (OUT, ">"."./maxDemon_results/classpositionHISTOsvm_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
# RFOR method
open (OUT, ">"."./maxDemon_results/classpositionHISTOrfor_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;
# ADABOOST method
open (OUT, ">"."./maxDemon_results/classpositionHISTOada_$fileIDq.txt\n");
print OUT "position\t"."sum\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $sum_class = 0;
   print OUT "$AApos\t"."$sum_class\n";
   }
close OUT;

############################################
############################################
if ($method_bnp == 1){
     print "running KNN classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassKNN");     
 for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nKNN classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   # read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "lengthN <- length(train)\n";
   print Rinput "Kval <- round(sqrt(lengthN))\n";
   print Rinput "knn.pred <- knn(train, test, class, k=Kval)\n";
   print Rinput "my_zero <- sum(knn.pred == 0)\n";
   print Rinput "my_one <- sum(knn.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassKNN/position"."_$r.txt')\n";
   print Rinput "print(knn.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
 
}


############################################
if ($method_dist == 1){ 
     print "running naive Bayes classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassNB");     
 for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nNB classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(e1071)\n";
			print Rinput "library(dplyr)\n";
   # read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "xy <- data.frame(class, train)\n";
			print Rinput "print('dimension of training set')\n";
			print Rinput "dim(xy)\n";
			if($sampletype eq 'y'){  # subsample training set option
			   print Rinput "xy <- sample_n(xy, 20000, replace = TRUE)\n";	  
			   print Rinput "print('new dimension of training set')\n";
			   print Rinput "dim(xy)\n";
			   }
   print Rinput "nb.xy   <- naiveBayes(as.factor(class)~., data=xy)\n";
   print Rinput "nb.pred <- predict(nb.xy, as.data.frame(test), type='class')\n";
   print Rinput "my_zero <- sum(nb.pred == 0)\n";
   print Rinput "my_one <- sum(nb.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassNB/position"."_$r.txt')\n";
   print Rinput "print(nb.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
 
}


############################################
if ($method_dist == 1){
     print "running linear discriminant analysis classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassLDA");     
 for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nLDA classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(MASS)\n";
			print Rinput "library(dplyr)\n";
   # read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "xy <- data.frame(class, train)\n";
			print Rinput "print('dimension of training set')\n";
			print Rinput "dim(xy)\n";
			if($sampletype eq 'y'){  # subsample training set option
			   print Rinput "xy <- sample_n(xy, 20000, replace = TRUE)\n";	  
			   print Rinput "print('new dimension of training set')\n";
			   print Rinput "dim(xy)\n";
			   }
   print Rinput "lda.xy   <- lda(as.factor(class)~., data=xy)\n";
   print Rinput "lda.pred <- predict(lda.xy, as.data.frame(test))\$class\n";
   print Rinput "my_zero <- sum(lda.pred == 0)\n";
   print Rinput "my_one <- sum(lda.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassLDA/position"."_$r.txt')\n";
   print Rinput "print(lda.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }
 
}

############################################
if ($method_dist == 1){
     print "running quadratic discriminant analysis classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassQDA");
  for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nQDA classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(MASS)\n";
			print Rinput "library(dplyr)\n";
   # read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "xy <- data.frame(class, train)\n";
			print Rinput "print('dimension of training set')\n";
			print Rinput "dim(xy)\n";
			if($sampletype eq 'y'){  # subsample training set option
			   print Rinput "xy <- sample_n(xy, 20000, replace = TRUE)\n";	  
			   print Rinput "print('new dimension of training set')\n";
			   print Rinput "dim(xy)\n";
			   }
   print Rinput "qda.xy   <- qda(as.factor(class)~., data=xy)\n";
   print Rinput "qda.pred <- predict(qda.xy, as.data.frame(test))\$class\n";
   print Rinput "my_zero <- sum(qda.pred == 0)\n";
   print Rinput "my_one <- sum(qda.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassQDA/position"."_$r.txt')\n";
   print Rinput "print(qda.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   }   

}


############################################
if ($method_kern == 1){
     print "running SVM classifier\n";
     mkdir("./testingData_$fileIDq/indAAclassSVM");
 for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nSVM classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(kernlab)\n";
			print Rinput "library(dplyr)\n";
			# read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "xy <- data.frame(class, train)\n";
			print Rinput "print('dimension of training set')\n";
			print Rinput "dim(xy)\n";
			if($sampletype eq 'y'){  # subsample training set option
			   print Rinput "xy <- sample_n(xy, 20000, replace = TRUE)\n";	  
			   print Rinput "print('new dimension of training set')\n";
			   print Rinput "dim(xy)\n";
			   }
			print Rinput "svm.xy   <- ksvm(as.factor(class)~., data=xy, kernel='$kerntype',C=2, cross=5)\n";
   print Rinput "svm.pred <- predict(svm.xy, as.data.frame(test), type='response')\n";
			print Rinput "my_zero <- sum(svm.pred == 0)\n";
   print Rinput "my_one <- sum(svm.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassSVM/position"."_$r.txt')\n";
   print Rinput "print(svm.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   } 

}


###########################################
if ($method_ens == 1){     
    print "running random forest classifier\n";
    mkdir("./testingData_$fileIDq/indAAclassRFOR");     
 for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nRFOR classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(randomForest)\n";
			print Rinput "library(dplyr)\n";
			# read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "xy <- data.frame(class, train)\n";
			print Rinput "print('dimension of training set')\n";
			print Rinput "dim(xy)\n";
			if($sampletype eq 'y'){  # subsample training set option
			   print Rinput "xy <- sample_n(xy, 20000, replace = TRUE)\n";	  
			   print Rinput "print('new dimension of training set')\n";
			   print Rinput "dim(xy)\n";
			   }
			print Rinput "rf.xy <- randomForest(as.factor(class)~., data=xy, ntree=500)\n";
   print Rinput "rf.pred <- predict(rf.xy, as.data.frame(test), type='response')\n";
			print Rinput "my_zero <- sum(rf.pred == 0)\n";
   print Rinput "my_one <- sum(rf.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassRFOR/position"."_$r.txt')\n";
   print Rinput "print(rf.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   } 
    
}


###########################################
if ($method_ens == 1){
    print "running adaptive boosting classifier\n";
    mkdir("./testingData_$fileIDq/indAAclassADA");     
 for (my $r = 1; $r<=$lengthID; $r++){
   print "\n\nADA classifier learning residue $r on $fileIDq\n\n";
   sleep(1);
   open (Rinput, "| R --vanilla")||die "could not start R command line\n";
   # load libraries
   print Rinput "library(ggplot2)\n";
   print Rinput "library(class)\n";
   print Rinput "library(ada)\n";
			print Rinput "library(dplyr)\n";
			# read data into R
   print Rinput "dataT = read.table('trainingData/position$r', header = TRUE)\n";
   print Rinput "dataD = read.table('testingData_$fileIDq/position$r', header = TRUE)\n";
   print Rinput "class <- dataT\$class\n"; # training class
   #print Rinput "print(dataD)\n";
   print Rinput "train <- dataT\n";
   print Rinput "train <- train[-c(1)]\n"; # drops class column
   print Rinput "test <- dataD\n";
   print Rinput "test <- test[-c(1)]\n"; # drops class column
   print Rinput "print(head(train))\n";
   print Rinput "print(head(test))\n";
   print Rinput "xy <- data.frame(class, train)\n";
			print Rinput "print('dimension of training set')\n";
			print Rinput "dim(xy)\n";
			if($sampletype eq 'y'){  # subsample training set option
			   print Rinput "xy <- sample_n(xy, 20000, replace = TRUE)\n";	  
			   print Rinput "print('new dimension of training set')\n";
			   print Rinput "dim(xy)\n";
			   }
			print Rinput "boost.xy <- ada(as.factor(class)~., data=xy, 18)\n";
   print Rinput "boost.pred <- predict(boost.xy, as.data.frame(test), type='vector')\n";
			print Rinput "my_zero <- sum(boost.pred == 0)\n";
   print Rinput "my_one <- sum(boost.pred == 1)\n";
   print Rinput "my_freq <- my_zero/(my_zero + my_one)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "mean_test <- mean(test[,1])\n";
   print Rinput "mean_train <- mean(train[,1])\n";
   print Rinput "delta_rmsf = (mean_test - mean_train)\n";
   print Rinput "sink('./testingData_$fileIDq/indAAclassADA/position"."_$r.txt')\n";
   print Rinput "print(boost.pred)\n";
   print Rinput "print(my_freq)\n";
   print Rinput "print(delta_rmsf)\n";
   print Rinput "sink()\n";
   # write to output file and quit R
   print Rinput "q()\n";# quit R 
   print Rinput "n\n";# save workspace image?
   close Rinput;
   } 
    
}

#############################################################################     
#############################################################################
sleep(1);
print "\n machine learning classification is completed\n\n";

print " copying flux profile training plots\n\n";
sleep(1);
my $oldfilename = "./DROIDS_results_$queryID"."_$refID"."_flux_$cutoffValue/DROIDSfluctuationAVGchain.txt";
my $newfilename = "./testingData_$fileIDq/DROIDSfluctuationAVGchain.txt";
copy($oldfilename, $newfilename);

##### make R Plots
if ($method_bnp == 1){
# KNN method
open (OUT, ">"."./maxDemon_results/classpositionHISTOknn_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassKNN/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}

if ($method_dist == 1){  # this is permanently disabled (replicates QDA)
# NB method
open (OUT, ">"."./maxDemon_results/classpositionHISTOnb_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassNB/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}

if ($method_dist == 1){
# LDA method
open (OUT, ">"."./maxDemon_results/classpositionHISTOlda_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassLDA/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}

if ($method_dist == 1){
# QDA method
open (OUT, ">"."./maxDemon_results/classpositionHISTOqda_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassQDA/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}

if ($method_kern == 1){
# SVM method
open (OUT, ">"."./maxDemon_results/classpositionHISTOsvm_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassSVM/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}

if ($method_ens == 1){ 
# RANDOM FOREST method
open (OUT, ">"."./maxDemon_results/classpositionHISTOrfor_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassRFOR/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}

if ($method_ens == 1){
# ADABOOST method
open (OUT, ">"."./maxDemon_results/classpositionHISTOada_$fileIDq.txt\n");
print OUT "position\t"."sum\t"."dRMSF\n";
for (my $r = 1; $r<=$lengthID; $r++){
   $AApos = $r;
   $freqML = 0.5;
   open (IN, "<"."./testingData_$fileIDq/indAAclassADA/position_$r.txt");
   my @IN = <IN>;
   for (my $i = 0; $i < scalar @IN; $i++) {
      $INrow = $IN[$i];
      @INrow = split(/\s+/, $INrow);
      $search = @INrow[0];
      $INrowNext = $IN[$i+1];
      @INrowNext = split(/\s+/, $INrowNext);
      $value1 = @INrowNext[1];
      $INrowNextNext = $IN[$i+2];
      @INrowNextNext = split(/\s+/, $INrowNextNext);
      $value2 = @INrowNextNext[1];
      if ($search eq "Levels:"){$freqML = $value1; $dRMSF = $value2; if($freqML eq 'NaN' || $freqML eq 'NA'){$freqML = 0.5;} }
      }
      close IN;
      if ($AApos ne ''){print OUT "$AApos\t"."$freqML\t"."$dRMSF\n";}
}
close OUT;
}
#####################################################
sleep(1);
print " plotting class position histogram and flux profiles\n\n";
sleep(1);

open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";

# read data into R
# KNNsleep(3);
print Rinput "datatableKNN = read.table('./maxDemon_results/classpositionHISTOknn_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_knn = "datatableKNN\$position"; # AA position
$sum_classifiers_knn = "datatableKNN\$sum"; # sum of classifiers
$dRMSF_knn = "datatableKNN\$dRMSF"; # avg dRMSF
print Rinput "dataframeKNN = data.frame(pos1=$AAposition_knn, Y1val=$sum_classifiers_knn)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./maxDemon_results/classpositionHISTOnb_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
$dRMSF_nb = "datatableNB\$dRMSF"; # avg dRMSF
print Rinput "dataframeNB = data.frame(pos1=$AAposition_nb, Y1val=$sum_classifiers_nb)\n";
# LDA
print Rinput "datatableLDA = read.table('./maxDemon_results/classpositionHISTOlda_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
$dRMSF_lda = "datatableLDA\$dRMSF"; # avg dRMSF
print Rinput "dataframeLDA = data.frame(pos1=$AAposition_lda, Y1val=$sum_classifiers_lda)\n";
# QDA
print Rinput "datatableQDA = read.table('./maxDemon_results/classpositionHISTOqda_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
$dRMSF_qda = "datatableQDA\$dRMSF"; # avg dRMSF
print Rinput "dataframeQDA = data.frame(pos1=$AAposition_qda, Y1val=$sum_classifiers_qda)\n";
# SVM
print Rinput "datatableSVM = read.table('./maxDemon_results/classpositionHISTOsvm_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
$dRMSF_svm = "datatableSVM\$dRMSF"; # avg dRMSF
print Rinput "dataframeSVM = data.frame(pos1=$AAposition_svm, Y1val=$sum_classifiers_svm)\n";
# random forest
print Rinput "datatableRFOR = read.table('./maxDemon_results/classpositionHISTOrfor_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
$dRMSF_rfor = "datatableRFOR\$dRMSF"; # avg dRMSF
print Rinput "dataframeRFOR = data.frame(pos1=$AAposition_rfor, Y1val=$sum_classifiers_rfor)\n";
# adaboost
print Rinput "datatableADA = read.table('./maxDemon_results/classpositionHISTOada_$fileIDq.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
$dRMSF_ada = "datatableADA\$dRMSF"; # avg dRMSF
print Rinput "dataframeADA = data.frame(pos1=$AAposition_ada, Y1val=$sum_classifiers_ada)\n";
# atom flux profiles
print Rinput "datatable2 = read.table('./testingData_$fileIDq/DROIDSfluctuationAVGchain.txt', header = TRUE)\n"; 
$AAposition2 = "datatable2\$pos_ref"; # AA position
$AAlabel = "datatable2\$res_ref"; # AA  identity
$trainingflux_ref = "datatable2\$flux_ref_avg"; # flux profile ref training
$trainingflux_query = "datatable2\$flux_query_avg"; # flux profile query training
print Rinput "dataframe2 = data.frame(pos2=$AAposition2, Y2val=$trainingflux_ref)\n";
print Rinput "dataframe3 = data.frame(pos3=$AAposition2, Y3val=$trainingflux_query)\n";

# lineplots
print Rinput "myplot1 <- ggplot() + ggtitle('       learning performance along protein backbone for MD validation set') + labs(x = 'position (residue number)', y = 'avg classification over time intervals') + geom_line(data = dataframeKNN, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'K nearest neighbors')) + geom_line(data = dataframeNB, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'naive Bayes')) + geom_line(data = dataframeLDA, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'linear discriminant - LDA')) + geom_line(data = dataframeQDA, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'quadratic discriminant - QDA')) + geom_line(data = dataframeSVM, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'support vector machine-SVM')) + geom_line(data = dataframeRFOR, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'random forest')) + geom_line(data = dataframeADA, mapping = aes(x = pos1, y = Y1val, ymin = 0, ymax = 1, color = 'adaboost')) + scale_color_brewer(palette='Set1')\n"; 
print Rinput "myplot2 <- ggplot() + labs(x = 'position (residue number)', y = 'avg FLUX') + geom_line(data = dataframe2, mapping = aes(x = pos2, y = Y2val, color = 'MD training set - reference protein=$refID')) + geom_line(data = dataframe3, mapping = aes(x = pos3, y = Y3val, color = 'MD training set - query protein=$queryID'))\n"; 

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(2, 1)))\n";
print Rinput "print(myplot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./maxDemon_results/fluxCorr_learningPLOT_$fileIDq.pdf";
copy($oldfilename, $newfilename);	
print " machine learning is complete\n\n";

close MUT;
} #end for loop
} # end outer for loop
####################################################################
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
      print " open PDF file for $fileIDq\n\n";
      sleep(2);
      print " close PDF viewer to continue\n\n";
      system "evince ./maxDemon_results/fluxCorr_learningPLOT_$fileIDq.pdf\n";
      } #end for loop
close MUT;

}# end outer for loop
########################################################################################################
sleep(1);
print "Do you want to graph the two-state molecular dynamics upon the classification space\n";
print "for each learning method at each amino acid site? (y or n)\n\n";
print "NOTE: requires python3, numpy, sklearn, and matplotlib installed\n\n";

my $YorN = <STDIN>;
chop($YorN);

if($YorN eq "y" || $YorN eq "Y" || $YorN eq "yes" || $YorN eq "YES"){
system "python3 classtraining_spaceplot_fluxCorr.py\n";     
sleep(3);
system "xdg-open ./trainingPlots/position1.png\n";
}
########################################################################################################
print "\n machine learning is complete\n\n";
exit;
###########################################################################################################
###########################################################################################################




