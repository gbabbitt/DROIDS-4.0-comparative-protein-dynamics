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
	 if ($header eq "chimerax_path"){$chimerax_path = $path;}
}
close IN;
print "path to Chimera .exe\t"."$chimera_path\n";
print "path to ChimeraX .exe\t"."$chimerax_path\n";

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


###########################################################################################################
###########################################################################################################

# prompt user - choose graphing option
sleep(1);
print "CHOOSE BAR PLOT TYPE - total mutational impact or just conserved regions\n";
print "NOTE: conserved regions option is intended only for very large proteins with proportionally small active sites and no allostery\n";
print "(type 'total', 'conserved' -note: default is 'total')\n\n";
my $bartype = <STDIN>;
chop($bartype);
if ($bartype eq ''){$bartype = 'total';}

print "Enter number position of N terminal on this chain\n";
print "(default = 0)\n\n";
my $startN = <STDIN>;
chop($startN);
if ($startN eq ''){$startN = 0;}

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

print " plotting canonical correlation and flux profiles\n\n";

$window = '';
$cancor_threshold = '';
# prompt user - choose best learning model to display
sleep(1);print "CHOOSE SIZE OF SLIDING WINDOW (e.g. default = 20 = 20 AA window for r value)\n\n";
my $window = <STDIN>;
chop($window);
if ($window eq ''){$window = 20};

# set conserved threshold
$cancor_threshold = 0.3;
sleep(1);print "CHOOSE LEVEL OF SIGNIFICANCE FOR CONSERVED CAN CORR (default = 0.01)\n\n";
my $cancor_threshold = <STDIN>;
chop($cancor_threshold);
if ($cancor_threshold eq ''){$cancor_threshold = 0.01};

open (Rinput, "| R --vanilla")||die "could not start R command line\n";

# load plotting libraries
print Rinput "library(ggplot2)\n";
print Rinput "library(gridExtra)\n";
print Rinput "library(CCA)\n";
print Rinput "library(CCP)\n";
print Rinput "library(boot)\n";
#print Rinput "options(scipen=999)\n"; # eliminates printing of sci notation

#############################################
###  analyze copies for plots 1 and 2
#############################################
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
# read test MD data into R
# KNN
print Rinput "datatableKNN = read.table('./maxDemon_results/classpositionHISTOknn_$copyID.txt', header = TRUE)\n"; 
$AAposition_knn = "datatableKNN\$position"; # AA position
$sum_classifiers_knn = "datatableKNN\$sum"; # sum of classifiers
print Rinput "dataframeKNN_$v = data.frame(pos1=$AAposition_knn+$startN, Y1val=$sum_classifiers_knn)\n";
#print Rinput "print(dataframeKNN_$v)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./maxDemon_results/classpositionHISTOnb_$copyID.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
print Rinput "dataframeNB_$v = data.frame(pos1=$AAposition_nb+$startN, Y1val=$sum_classifiers_nb)\n";
#print Rinput "print(dataframeNB_$v)\n";
# LDA
print Rinput "datatableLDA = read.table('./maxDemon_results/classpositionHISTOlda_$copyID.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
print Rinput "dataframeLDA_$v = data.frame(pos1=$AAposition_lda+$startN, Y1val=$sum_classifiers_lda)\n";
#print Rinput "print(dataframeLDA_$v)\n";
# QDA
print Rinput "datatableQDA = read.table('./maxDemon_results/classpositionHISTOqda_$copyID.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
print Rinput "dataframeQDA_$v = data.frame(pos1=$AAposition_qda+$startN, Y1val=$sum_classifiers_qda)\n";
#print Rinput "print(dataframeQDA_$v)\n";
# SVM
print Rinput "datatableSVM = read.table('./maxDemon_results/classpositionHISTOsvm_$copyID.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
print Rinput "dataframeSVM_$v = data.frame(pos1=$AAposition_svm+$startN, Y1val=$sum_classifiers_svm)\n";
#print Rinput "print(dataframeSVM_$v)\n";
# random forest
print Rinput "datatableRFOR = read.table('./maxDemon_results/classpositionHISTOrfor_$copyID.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
print Rinput "dataframeRFOR_$v = data.frame(pos1=$AAposition_rfor+$startN, Y1val=$sum_classifiers_rfor)\n";
#print Rinput "print(dataframeRFOR_$v)\n";
# adaboost
print Rinput "datatableADA = read.table('./maxDemon_results/classpositionHISTOada_$copyID.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
print Rinput "dataframeADA_$v = data.frame(pos1=$AAposition_ada+$startN, Y1val=$sum_classifiers_ada)\n";
#print Rinput "print(dataframeADA_$v)\n";

} #end for loop

# atom flux profiles
print Rinput "datatable2 = read.table('./testingData_$queryID/DROIDSfluctuationAVGchain.txt', header = TRUE)\n"; 
$AAposition2 = "datatable2\$pos_ref"; # AA position
$AAlabel = "datatable2\$res_ref"; # AA  identity
$trainingflux_ref = "datatable2\$flux_ref_avg"; # flux profile ref training
$trainingflux_query = "datatable2\$flux_query_avg"; # flux profile query training
print Rinput "dataframe2 = data.frame(pos2=$AAposition2+$startN, Y2val=$trainingflux_ref)\n";
print Rinput "dataframe3 = data.frame(pos3=$AAposition2+$startN, Y3val=$trainingflux_query)\n";


# canonical correlations
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$copyID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$copyID)\n";
} # end for loop

# output stats file
print Rinput "sink('./testingData_$queryID/cancor_stats.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and new variants\n\n')\n";
print Rinput "cat('there are 6 dimensions and 7 learners\n\n')\n";
print Rinput "cat('cor = overall correlation of all learners in each dimension\n')\n";
print Rinput "cat('coef = slope of each pairwise learner correlation for each dimension\n')\n";
print Rinput "cat('center = mean of each learner\n\n')\n";
print Rinput "cat('1 = KNN, 2 = naive Bayes, 3 = LDA, 4 = QDA, 5 = SVM, 6 = random forest, 7 = adaboost\n\n')\n";
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $copyID\n\n')\n";
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$copyID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$copyID)\n";

print Rinput "sink()\n";

# AA sliding window for canonical correlation
print Rinput "my_cancors_$copyID <- c()\n";
print Rinput "my_pvals_$copyID <- c()\n";
print Rinput "my_ccpos_$copyID <- c()\n";
print Rinput "my_bars_$copyID <- c()\n";
print Rinput "mylength <- length(newMD[,1])\n";
print Rinput "print(mylength)\n";
print Rinput "for(i in 1:(mylength-$window)){
   newMDslice <- newMD[i:(i+$window),1] 
   queryMDslice <- queryMD[i:(i+$window),1];
   queryMDrefslice <- queryMDref[i:(i+$window),1];
   mytest1 <- cancor(newMDslice, queryMDrefslice)
   mycor <- mytest1\$cor
   mytest2 <- cancor(queryMDslice, queryMDrefslice)
   myselfcor <- mytest2\$cor
   N = length(queryMDslice)
   p = 1
   q = 1
   mypvaltest <- p.asym(myselfcor, N, p, q, tstat = 'Wilks')
   mypval <- mypvaltest\$p.value[1]
   if(mypval <= $cancor_threshold){mybar = 1}
   if(mypval > $cancor_threshold){mybar = 0}
   if($v == 0){
   my_bars_$queryID <- c(my_bars_$queryID, mybar)
   }
   my_ccpos_$copyID <- c(my_ccpos_$copyID, i)
   my_cancors_$copyID <- c(my_cancors_$copyID, mycor)
   my_pvals_$copyID <- c(my_pvals_$copyID, mypval)
   }\n";
print Rinput "print(my_cancors_$copyID)\n";
print Rinput "print(my_pvals_$copyID)\n";
print Rinput "my_adjpvals_$copyID <- p.adjust(my_pvals_$copyID, method = 'fdr', n = length(my_pvals_$copyID))\n";
print Rinput "print(my_adjpvals_$copyID)\n";
print Rinput "print(my_ccpos_$copyID)\n";
print Rinput "print(my_bars_$queryID)\n";
print Rinput "my_adjbars_$queryID <- c()\n";
print Rinput "print(my_adjpvals_$copyID\[i])\n";
print Rinput "for(j in 1:(length(my_adjpvals_$copyID))){
   print(my_adjpvals_$copyID\[j])
   if (my_adjpvals_$copyID\[j] <= $cancor_threshold){myadjbar = 1}
   if (my_adjpvals_$copyID\[j] > $cancor_threshold){myadjbar = 0}
   my_adjbars_$queryID <- c(my_adjbars_$queryID, myadjbar) 
   }\n";
print Rinput "print(my_adjbars_$queryID)\n";  
print Rinput "sink('./testingData_$queryID/cancor_stats_$copyID.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $copyID ($window AA sliding window)\n\n')\n";
print Rinput "print(my_cancors_$copyID)\n";
print Rinput "sink()\n";
print Rinput "dataframe4_$copyID = data.frame(pos=my_ccpos_$copyID+$startN, cc=my_cancors_$copyID, bars=my_adjbars_$queryID)\n";
print Rinput "my_overall <- cancor(newMD, queryMD)\n"; # overall
print Rinput "my_ccor_$copyID <- my_overall\$cor[1]\n";
print Rinput "print(my_ccor_$copyID)\n";
# create output mask movie render file for Wilks lambda
print Rinput "sink('./testingData_$queryID/adj_pvals_$copyID.txt')\n";
print Rinput "print(my_adjbars_$queryID)\n";
print Rinput "sink()\n";
} # end for loop

# lineplots
# create plot 1
print Rinput "myplot1 <- ggplot() + ggtitle('learning performance - 2 MD validation sets') + labs(x = 'position (residue number)', y = 'avg class over time intervals') + theme(legend.position = 'none', plot.title = element_text(size = 10))\n"; 
print Rinput "mylist1 <- list()\n";
# create lines for plot 1
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     print Rinput "mylines_$v <- list(geom_line(data = dataframeKNN_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID', ymin = 0, ymax = 1)), geom_line(data = dataframeNB_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeLDA_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeQDA_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeSVM_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeRFOR_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')), geom_line(data = dataframeADA_$v, mapping = aes(x = pos1, y = Y1val, color = '$copyID')))\n"; 
     print Rinput "mylist1 <- list(mylist1, mylines_$v)\n";
} # end for loop
print Rinput "myplot1final <- myplot1 + mylist1 + scale_color_brewer(palette='Set1')\n"; 
# create plot 2
print Rinput "myplot2 <- ggplot() + ggtitle('cancor for $window AA window + conserved dynamics zones -Wilks lambda') + geom_area(data = dataframe4_$copyID, mapping = aes(x = pos, y = bars, color = 'conserved dynamics'), alpha = 0.8) + labs(x = 'position (residue number)', y = 'local R value for cancor') + theme(legend.title = element_blank())\n"; 
print Rinput "mylist2 <- list()\n";
# create lines for plot 2
for (my $v = 0; $v < scalar(@copies); $v++){
     if ($v == 0) {next;} # skip first line
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     print Rinput "myline_$v <- geom_line(data = dataframe4_$copyID, mapping = aes(x = pos, y = cc, color = '$variantID_label', ymin = 0, ymax = 1))\n"; 
     print Rinput "mylist2 <- list(mylist2, myline_$v)\n";
} # end for loop

if (scalar(@variants) <= 8){print Rinput "myplot2final <- myplot2 + mylist2 + scale_color_brewer(palette='Set2')\n";} 
if (scalar(@variants) > 8){print Rinput "myplot2final <- myplot2 + mylist2\n";} 

# exit if no variant list to analyze

if (scalar(@variants) == 1 || scalar(@variants) == 0){
print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 1)))\n";
print Rinput "print(myplot1final, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n";
print Rinput "print(myplot2final, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))\n";
#print Rinput "print(myplot2final, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";

# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$queryID/learningPLOTcompare.pdf";
copy($oldfilename, $newfilename);	
print " machine learning is complete\n\n";
print " close PDF and txt viewer to continue\n\n";
system "evince ./testingData_$queryID/learningPLOTcompare.pdf\n";

# open can corr results
system "gedit ./testingData_$queryID/cancor_stats.txt\n";
for (my $v = 0; $v < scalar(@copies); $v++){
     $copyID = $copies[$v];
     $queryID = $copies[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     system "gedit ./testingData_$queryID/cancor_stats_$copyID.txt\n";
} # end for loop  
     
exit;
}


#############################################
###  analyze variants for plot 3
#############################################


for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     
     
# read test MD data into Rprint Rinput "my_pvals_$variantID <- c()\n";
# KNN
print Rinput "datatableKNN = read.table('./maxDemon_results/classpositionHISTOknn_$variantID.txt', header = TRUE)\n"; 
$AAposition_knn = "datatableKNN\$position"; # AA position
$sum_classifiers_knn = "datatableKNN\$sum"; # sum of classifiers
print Rinput "dataframeKNN_$v = data.frame(pos1=$AAposition_knn+$startN, Y1val=$sum_classifiers_knn)\n";
# naive Bayes
print Rinput "datatableNB = read.table('./maxDemon_results/classpositionHISTOnb_$variantID.txt', header = TRUE)\n"; 
$AAposition_nb = "datatableNB\$position"; # AA position
$sum_classifiers_nb = "datatableNB\$sum"; # sum of classifiers
print Rinput "dataframeNB_$v = data.frame(pos1=$AAposition_nb+$startN, Y1val=$sum_classifiers_nb)\n";
# LDA
print Rinput "datatableLDA = read.table('./maxDemon_results/classpositionHISTOlda_$variantID.txt', header = TRUE)\n"; 
$AAposition_lda = "datatableLDA\$position"; # AA position
$sum_classifiers_lda = "datatableLDA\$sum"; # sum of classifiers
print Rinput "dataframeLDA_$v = data.frame(pos1=$AAposition_lda+$startN, Y1val=$sum_classifiers_lda)\n";
# QDA
print Rinput "datatableQDA = read.table('./maxDemon_results/classpositionHISTOqda_$variantID.txt', header = TRUE)\n"; 
$AAposition_qda = "datatableQDA\$position"; # AA position
$sum_classifiers_qda = "datatableQDA\$sum"; # sum of classifiers
print Rinput "dataframeQDA_$v = data.frame(pos1=$AAposition_qda+$startN, Y1val=$sum_classifiers_qda)\n";
# SVM
print Rinput "datatableSVM = read.table('./maxDemon_results/classpositionHISTOsvm_$variantID.txt', header = TRUE)\n"; 
$AAposition_svm = "datatableSVM\$position"; # AA position
$sum_classifiers_svm = "datatableSVM\$sum"; # sum of classifiers
print Rinput "dataframeSVM_$v = data.frame(pos1=$AAposition_svm+$startN, Y1val=$sum_classifiers_svm)\n";
# random forest
print Rinput "datatableRFOR = read.table('./maxDemon_results/classpositionHISTOrfor_$variantID.txt', header = TRUE)\n"; 
$AAposition_rfor = "datatableRFOR\$position"; # AA position
$sum_classifiers_rfor = "datatableRFOR\$sum"; # sum of classifiers
print Rinput "dataframeRFOR_$v = data.frame(pos1=$AAposition_rfor+$startN, Y1val=$sum_classifiers_rfor)\n";
# adaboost
print Rinput "datatableADA = read.table('./maxDemon_results/classpositionHISTOada_$variantID.txt', header = TRUE)\n"; 
$AAposition_ada = "datatableADA\$position"; # AA position
$sum_classifiers_ada = "datatableADA\$sum"; # sum of classifiers
print Rinput "dataframeADA_$v = data.frame(pos1=$AAposition_ada+$startN, Y1val=$sum_classifiers_ada)\n";

} #end for loop

$varcount = scalar(@variants)+1;
# canonical correlations
for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$variantID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$variantID)\n";
} # end for loop

# output stats file
print Rinput "sink('./testingData_$queryID/cancor_stats.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and new variants\n\n')\n";
print Rinput "cat('there are 6 dimensions and 7 learners\n\n')\n";
print Rinput "cat('cor = overall correlation of all learners in each dimension\n')\n";
print Rinput "cat('coef = slope of each pairwise learner correlation for each dimension\n')\n";
print Rinput "cat('center = mean of each learner\n\n')\n";
print Rinput "cat('1 = KNN, 2 = naive Bayes, 3 = LDA, 4 = QDA, 5 = SVM, 6 = random forest, 7 = adaboost\n\n')\n";

print Rinput "my_impact_IDs <- c()\n";
print Rinput "my_impact_sums <- c()\n";
print Rinput "my_impact_sd <- c()\n";
print Rinput "my_impact_n <- c()\n";
print Rinput "my_impact_sums_upperSE <- c()\n";
print Rinput "my_impact_sums_lowerSE <- c()\n";
print Rinput "my_impact_cons_sums <- c()\n";
print Rinput "my_impact_cons_sums_upperSE <- c()\n";
print Rinput "my_impact_cons_sums_lowerSE <- c()\n";
print Rinput "my_impact_cons_sums <- c()\n";
print Rinput "my_impact_cons_sd <- c()\n";
print Rinput "my_impact_cons_n <- c()\n";

for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $copyID\n\n')\n";
print Rinput "newMD <- cbind(dataframeKNN_".$v."[2], dataframeNB_".$v."[2], dataframeLDA_".$v."[2], dataframeQDA_".$v."[2], dataframeSVM_".$v."[2], dataframeRFOR_".$v."[2], dataframeADA_".$v."[2])\n";
print Rinput "queryMDref <- cbind(dataframeKNN_0[2], dataframeNB_0[2], dataframeLDA_0[2], dataframeQDA_0[2], dataframeSVM_0[2], dataframeRFOR_0[2], dataframeADA_0[2])\n";
print Rinput "queryMD <- cbind(dataframeKNN_1[2], dataframeNB_1[2], dataframeLDA_1[2], dataframeQDA_1[2], dataframeSVM_1[2], dataframeRFOR_1[2], dataframeADA_1[2])\n";
print Rinput "myCC_$variantID <- cancor(newMD, queryMD)\n"; # overall
print Rinput "print(myCC_$variantID)\n";

print Rinput "sink()\n";

# AA sliding window for canonical correlation
print Rinput "my_cancors_$variantID <- c()\n";
print Rinput "my_pvals_$variantID <- c()\n";
print Rinput "my_ccpos_$variantID <- c()\n";
print Rinput "my_impact_$variantID <- c()\n";
print Rinput "my_impact_sum_$variantID = 0\n";
print Rinput "my_impact_cons_sum_$variantID = 0\n";
print Rinput "my_impact_sd_$variantID = 0\n";
print Rinput "my_impact_cons_sd_$variantID = 0\n";
print Rinput "my_impact_n_$variantID = 0\n";
print Rinput "my_impact_length_$variantID = 0\n";
print Rinput "my_impact_cons_length_$variantID = 0\n";
print Rinput "my_bars_$variantID <- c()\n";
print Rinput "mylength <- length(newMD[,1])\n";
print Rinput "print(mylength)\n";
print Rinput "for(i in 1:(mylength-$window)){
   newMDslice <- newMD[i:(i+$window),1] 
   queryMDslice <- queryMD[i:(i+$window),1];
   queryMDrefslice <- queryMDref[i:(i+$window),1];
   mytest1 <- cancor(newMDslice, queryMDrefslice)
   mycor <- mytest1\$cor
   mytest2 <- cancor(queryMDslice, queryMDrefslice)
   myselfcor <- mytest2\$cor
   N = length(queryMDslice)
   p = 1
   q = 1
   mypvaltest <- p.asym(myselfcor, N, p, q, tstat = 'Wilks')
   mypval <- mypvaltest\$p.value[1]
   if(mypval <= $cancor_threshold){mybar = 1}
   if(mypval > $cancor_threshold){mybar = 0}
   if($v == 0){
   my_bars_$queryID <- c(my_bars_$queryID, mybar)
   }
   myimpact = abs(myselfcor*(log(mycor/myselfcor)))
   my_impact_sum_$variantID = my_impact_sum_$variantID + myimpact
   if(mybar == 1){my_impact_cons_sum_$variantID = my_impact_cons_sum_$variantID + myimpact}
   if(mybar == 0){my_impact_cons_sum_$variantID = my_impact_cons_sum_$variantID + 0}
   if ($v > 0){ 
     if (myimpact >= myupper){myimpact = abs(myselfcor*(log(mycor/myselfcor)))}
     if (myimpact < myupper){myimpact = 0}
   }
   my_impact_length_$variantID = my_impact_length_$variantID + 1
   if(mybar == 1){my_impact_cons_length_$variantID = my_impact_cons_length_$variantID + 1}
   
   my_impact_$variantID <- c(my_impact_$variantID, myimpact)
   my_ccpos_$variantID <- c(my_ccpos_$variantID, i)
   my_cancors_$variantID <- c(my_cancors_$variantID, mycor)
   my_pvals_$variantID <- c(my_pvals_$variantID, mypval)
   }\n";  # NOTE: can set impact to zero if less than upper bound
print Rinput "print(my_cancors_$variantID)\n";
print Rinput "print(my_pvals_$variantID)\n";
print Rinput "my_adjpvals_$variantID <- p.adjust(my_pvals_$variantID, method = 'fdr', n = length(my_pvals_$variantID))\n";
print Rinput "print(my_adjpvals_$variantID)\n";
print Rinput "print(my_ccpos_$variantID)\n";
print Rinput "print(my_impact_$variantID)\n";
print Rinput "custom.bootTTL <- function(times, data=my_impact_$variantID, length=my_impact_length_$variantID) { 
  boots <- rep(NA, times)
  for (i in 1:times) {
    boots[i] <- sd(sample(data, length, replace=TRUE))/sqrt(length)  
  }
  boots
}\n";
print Rinput "custom.bootCONS <- function(times, data=my_impact_$variantID, length=my_impact_cons_length_$variantID) { 
  boots <- rep(NA, times)
  for (i in 1:times) {
    boots[i] <- sd(sample(data, length, replace=TRUE))/sqrt(length)  
  }
  boots
}\n";
print Rinput "mySE <- sum(custom.bootTTL(times=1000))\n"; # bootstrap standard error of sum
print Rinput "print(mySE)\n";
print Rinput "myUpperSE <- (my_impact_sum_$variantID + mySE)\n";
print Rinput "myLowerSE <- (my_impact_sum_$variantID - mySE)\n";
print Rinput "myN = $varcount\n";
print Rinput "mySD = mySE\n";
print Rinput "print(mySD)\n";
print Rinput "print(my_impact_sum_$variantID)\n";
print Rinput "print(myUpperSE)\n";
print Rinput "print(myLowerSE)\n";
print Rinput "mySE_cons <- sum(custom.bootCONS(times=1000))\n"; # bootstrap standard error of sum
print Rinput "print(mySE_cons)\n";
print Rinput "myUpperSE_cons <- (my_impact_cons_sum_$variantID + mySE_cons)\n";
print Rinput "myLowerSE_cons <- (my_impact_cons_sum_$variantID - mySE_cons)\n";
print Rinput "mySD_cons = mySE_cons\n";
print Rinput "print(my_impact_cons_sum_$variantID)\n";
print Rinput "print(myUpperSE_cons)\n";
print Rinput "print(myLowerSE_cons)\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_IDs <- c('$variantID_label', my_impact_IDs)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sums <- c(my_impact_sum_$variantID, my_impact_sums)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sums_upperSE <- c(my_impact_sum_$variantID+mySE, my_impact_sums_upperSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sums_lowerSE <- c(my_impact_sum_$variantID-mySE, my_impact_sums_lowerSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_sd <- c(mySD, my_impact_sd)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_n <- c(myN, my_impact_n)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sums <- c(my_impact_cons_sum_$variantID, my_impact_cons_sums)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sums_upperSE <- c(my_impact_cons_sum_$variantID+mySE_cons, my_impact_cons_sums_upperSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sums_lowerSE <- c(my_impact_cons_sum_$variantID-mySE_cons, my_impact_cons_sums_lowerSE)}\n";
print Rinput "if(my_impact_sum_$variantID>0){my_impact_cons_sd <- c(mySD_cons, my_impact_cons_sd)}\n";
print Rinput "print(my_bars_$queryID)\n";
print Rinput "my_adjbars_$queryID <- c()\n";
print Rinput "print(my_adjpvals_$variantID\[i])\n";
print Rinput "for(j in 1:(length(my_adjpvals_$variantID))){
   print(my_adjpvals_$variantID\[j])
   if (my_adjpvals_$variantID\[j] <= $cancor_threshold){myadjbar = 1}
   if (my_adjpvals_$variantID\[j] > $cancor_threshold){myadjbar = 0}
   my_adjbars_$queryID <- c(my_adjbars_$queryID, myadjbar) 
   }\n";
print Rinput "print(my_adjbars_$queryID)\n";  
print Rinput "sink('./testingData_$queryID/cancor_stats_$variantID.txt')\n";
print Rinput "cat('CANONICAL CORRELATION ANALYSIS OF $queryID and $variantID ($window AA sliding window)\n\n')\n";
print Rinput "print(my_cancors_$variantID)\n";
print Rinput "cat('\nMUTATIONAL IMPACT ANALYSIS OF $queryID and $variantID ($window AA sliding window)\n\n')\n";
print Rinput "print(my_impact_$variantID)\n";
print Rinput "sink()\n";
print Rinput "dataframe4_$variantID = data.frame(pos=my_ccpos_$variantID+$startN, cc=my_cancors_$variantID, bars=my_adjbars_$queryID)\n";
print Rinput "dataframe5_$variantID = data.frame(pos=my_ccpos_$variantID+$startN, cc=my_impact_$variantID)\n";
print Rinput "my_overall <- cancor(newMD, queryMD)\n"; # overall
print Rinput "my_ccor_$variantID <- my_overall\$cor[1]\n";
print Rinput "print(my_ccor_$variantID)\n";

# get bootstrap if first comparison
if($variantID eq $variants[0]){
    print Rinput "print(my_impact_$variantID)\n";
    print Rinput "mymean <- mean(my_impact_$variantID)\n";
    print Rinput "mysd <- sd(my_impact_$variantID)\n";
    print Rinput "myupper <- mymean + (2*mysd)\n";
    print Rinput "print(myupper)\n";
    }
# create output impact files for movie render
print Rinput "sink('./testingData_$queryID/impact_$variantID.txt')\n";
print Rinput "print(my_impact_$variantID)\n";
print Rinput "sink()\n";   
    
} # end for loop

# create plot 3
print Rinput "myplot3 <- ggplot() + ggtitle('signif local variant impact for $window AA window (+2sd of validation set)') + labs(x = 'position (residue number)', y = 'local index of impact') + theme(legend.title = element_blank())\n"; 
print Rinput "mylist3 <- list()\n";
# create lines for plot 3
for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     $variantID_label = $variant_labels[$v];
     $queryID_label = $variant_labels[0];
     print Rinput "myline_$v <- geom_line(data = dataframe5_$variantID, mapping = aes(x = pos, y = cc, color = '$variantID_label'))\n"; 
     print Rinput "mylist3 <- list(mylist3, myline_$v)\n";
} # end for loop
if(scalar(@variants) <= 8){print Rinput "myplot3final <- myplot3 + mylist3 + scale_color_brewer(palette='Set2')\n";} 
if(scalar(@variants) > 8){print Rinput "myplot3final <- myplot3 + mylist3\n";} 

print Rinput "print(my_impact_IDs)\n";
print Rinput "print(my_impact_sums)\n";
print Rinput "print(my_impact_cons_sums)\n";
print Rinput "print(my_impact_sd)\n";
print Rinput "print(my_impact_cons_sd)\n";
print Rinput "print(my_impact_n)\n";

# create plot 4 
if ($bartype eq "total" && scalar(@variants) <= 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_sums, upperSE = my_impact_sums_upperSE, lowerSE = my_impact_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, y = sum, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + scale_fill_brewer(palette = 'Set2') + ggtitle('total mutational impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "conserved" && scalar(@variants) <= 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, upperSE = my_impact_cons_sums_upperSE, lowerSE = my_impact_cons_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + scale_fill_brewer(palette = 'Set2') + ggtitle('conserved region impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "total" && scalar(@variants) > 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_sums, upperSE = my_impact_sums_upperSE, lowerSE = my_impact_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, y = sum, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + ggtitle('total mutational impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}
if ($bartype eq "conserved" && scalar(@variants) > 8){
print Rinput "dataframe6 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, upperSE = my_impact_cons_sums_upperSE, lowerSE = my_impact_cons_sums_lowerSE)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myplot4 <- ggplot() + geom_col(data = dataframe6, aes(x = var, y = sum, fill = var)) + geom_errorbar(data = dataframe6, aes(x=var, ymin=lowerSE, ymax=upperSE), width=0.5, colour='gray', alpha=1, size=1) + ggtitle('conserved region impact (signif + ns)') + labs(x = 'variant', y = 'sum impact') + theme(legend.position = 'none', axis.text.x = element_text(angle = 30))\n"; 
print Rinput "myplot4final <- myplot4\n"; 
}

if ($bartype eq "total" && scalar(@variants) >= 4){ # ANOVA table from summary data
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_sums, sd = mySE, n = my_impact_n)\n";
print Rinput "print(dataframe7)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print(pval)\n";
print Rinput "sink('./maxDemon_results/ANOVAresult.txt')\n";
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_sums, sd = mySE, n = my_impact_n)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print('Fval')\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print('numerator DF')\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print('denominator DF')\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print('pval')\n";
print Rinput "print(pval)\n";
print Rinput "sink()\n"; 
}
if ($bartype eq "conserved" && scalar(@variants) >= 4){ # ANOVA table from summary data
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, sd = mySE_cons, n = my_impact_n)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print(pval)\n";
print Rinput "sink('./maxDemon_results/ANOVAresult.txt')\n";
print Rinput "library(rpsychi)\n";     
print Rinput "dataframe7 = data.frame(var = my_impact_IDs, sum = my_impact_cons_sums, sd = mySE_cons, n = my_impact_n)\n";
print Rinput "print(dataframe6)\n";
print Rinput "myANOVA <- with(dataframe7, ind.oneway.second(sum,sd,n))\n";
print Rinput "print(myANOVA)\n";
print Rinput "fval <- myANOVA\$anova.table[1,4]\n";
print Rinput "print('Fval')\n";
print Rinput "print(fval)\n";
print Rinput "numDF <- myANOVA\$anova.table[1,2]\n";
print Rinput "print('numerator DF')\n";
print Rinput "print(numDF)\n";
print Rinput "denomDF <- myANOVA\$anova.table[2,2]\n";
print Rinput "print('denominator DF')\n";
print Rinput "print(denomDF)\n";
print Rinput "pval <- pf(q=fval, df1=numDF, df2=denomDF, lower.tail=FALSE)\n";
print Rinput "print('pval')\n";
print Rinput "print(pval)\n";
print Rinput "sink()\n";
}
#####################################

print Rinput "library(grid)\n";
print Rinput "pushViewport(viewport(layout = grid.layout(3, 2)))\n";
print Rinput "print(myplot2final, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))\n";
print Rinput "print(myplot3final, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))\n";
print Rinput "print(myplot1final, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))\n";
print Rinput "print(myplot4final, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))\n";
# write to output file and quit R
print Rinput "q()\n";# quit R 
print Rinput "n\n";# save workspace image?
close Rinput;
print "\n\n";
print " copying plot\n\n";
sleep(1);
#mkdir("maxDemon_results");
my $oldfilename = "Rplots.pdf";
my $newfilename = "./testingData_$queryID/learningPLOTcompare.pdf";
copy($oldfilename, $newfilename);	
my $oldfilename2 = "Rplots.pdf";
my $newfilename2 = "./maxDemon_results/learningPLOTcompare.pdf";
copy($oldfilename2, $newfilename2);
print " machine learning is complete\n\n";
print " close PDF and txt viewer to continue\n\n";
system "evince ./testingData_$queryID/learningPLOTcompare.pdf\n";

# open can corr results
system "gedit ./testingData_$queryID/cancor_stats.txt\n";
     my $oldfilename3 = "./testingData_$queryID/cancor_stats.txt";
     my $newfilename3 = "./maxDemon_results/cancor_stats.txt";
     copy($oldfilename3, $newfilename3);
for (my $v = 0; $v < scalar(@variants); $v++){
     $variantID = $variants[$v];
     $queryID = $variants[0];
     system "gedit ./testingData_$queryID/cancor_stats_$variantID.txt\n";
     my $oldfilename4 = "./testingData_$queryID/cancor_stats_$variantID.txt";
     my $newfilename4 = "./maxDemon_results/cancor_stats_$variantID.txt";
     copy($oldfilename4, $newfilename4);
} # end for loop
if(scalar(@variants) >= 4){system "gedit ./maxDemon_results/ANOVAresult.txt\n";}
#}

###################################################
##  mapping conserved dynamics to PDB structure
###################################################

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
	 if ($header eq "chimerax_path"){$chimerax_path = $path;}
}
close IN;
print "path to Chimera .exe\t"."$chimera_path\n";
print "path to ChimeraX .exe\t"."$chimerax_path\n";

#### This uses a GUI to write the control files needed for the DROIDS scripts ####
print "\n\nWelcome to DROIDS- Detecting Relative Outlier Impacts
                          in Dynamic Simulations

- visual toolbox for functional evolutionary comparison
  of molecular dynamic simulation \n\n";

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




sleep(1);print "\nSELECT INTERACTION TYPE FOR MAPPING (1=protein only | 2=protein-protein | 3=DNA-protein | 4=protein-ligand)\n\n";
my $stype_number = <STDIN>;
chop($stype_number);

if($stype_number == 1){$stype = "protein";}
if($stype_number == 2){$stype = "protprot";}
if($stype_number == 3){$stype = "dna";}
if($stype_number == 4){$stype = "ligand";}

if($stype eq "dna" || $stype eq "ligand"){$orig_queryID = $queryID; $queryID = $refID;}

# enforce orig $queryID from training set
print "queryID label is "."$queryID\n";
if($queryID =~ m/_1_1/){$queryID = substr($queryID, 0, length($queryID)-4);}
if($queryID =~ m/_1/){$queryID = substr($queryID, 0, length($queryID)-2);}
print "queryID label is adjusted to "."$queryID\n";
$orig_queryID = $queryID;  # create tag for original query in training set
$queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

print "\nEnter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound
   structure\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot' || $stype eq 'protein'){
   print "\nEnter number position of N terminal on this chain (default = 0)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }



# create mask for signif Wilks lamda
open(OUT, ">"."./testingData_$queryID/adj_vertpvals_$queryID.txt")||die "could not create mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
open(IN, "<"."./testingData_$queryID/adj_pvals_$queryID.txt")||die "could not open mask file /testingData_$queryID/adj_pvals_$queryID.txt\n";  
my @IN = <IN>;
for (my $i = 0; $i < scalar @IN; $i++) {
	$INrow = $IN[$i];
     @INrow = split(/] /, $INrow);
	$trunINrow = $INrow[1];
     #print OUT "$trunINrow";
     @trunINrow = split(/\s+/, $trunINrow);
     for (my $ii = 0; $ii < scalar(@trunINrow); $ii++){
     $value = $trunINrow[$ii];
     print OUT "$value\n";
     }
     }
close IN;
close OUT;

# make learned class attribute files for image rendering
  open(OUT, ">"."./attribute_files/classATTR_$refID"."_mask".".dat")||die "could not create ATTR time series file\n";
  #if ($stype ne "protprot"){open(OUT, ">"."./attribute_files/classATTR_$refID"."_mask".".dat")||die "could not create ATTR time series file\n";}  
  #if ($stype eq "protprot"){open(OUT, ">"."./attribute_files/classATTR_$orig_queryID"."_mask".".dat")||die "could not create ATTR time series file\n";}  
  print OUT "recipient: residues\n";
  print OUT "attribute: class\n";
  print OUT "\n";
  
   for (my $a = 0; $a < $lengthID+$offset+1; $a++){
   if ($a eq '' || $a <= $offset){next;}
   open(IN, "<"."./testingData_$queryID/adj_vertpvals_$queryID.txt")||die "could not open mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
    my @IN = <IN>;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $maskvalue = $INrow[0];
   
     #print "$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     if ($maskvalue == 1){print OUT "\t:"."$pos\t"."1\n";}
     #if ($f == $frame && $maskvalue == 1){print OUT "\t:"."$pos\t"."1\n";} # to test mask positioning
     if ($maskvalue == 0){print OUT "\t:"."$pos\t"."0\n";}
     
 
 close IN;
 }
 close OUT;

print "\n class attribute files for all frames are created\n";
sleep(1);
###############################################################
print("Preparing display...\n");
print("close ChimeraX program to exit\n\n");
# copy visualization support files to folder
mkdir("ChimeraXvis");
copy("$queryID"."REDUCED.pdb", "ChimeraXvis/query.pdb");
copy("$refID"."REDUCED.pdb", "ChimeraXvis/reference.pdb");
copy("./attribute_files/classATTR_$refID"."_mask".".dat", "ChimeraXvis/attributeCLASS.dat");
# create control file for visualization
open(CTL, ">"."ChimeraXvis.ctl");
print CTL "model\t"."#1\n";
print CTL "chain\t"."/$chainMAP\n";
print CTL "structure\t"."ChimeraXvis/reference.pdb\n";
print CTL "structureADD\t"."ChimeraXvis/query.pdb\n";
print CTL "attr_file\t"."ChimeraXvis/attributeCLASS.dat\n";
print CTL "length\t"."$lengthID\n";
print CTL "attr\t"."class\n";
print CTL "minval\t"."0\n";
print CTL "maxval\t"."1\n";
print CTL "palette\t"."Greys-5\n";
print CTL "lighting\t"."simple\n";
print CTL "transparency\t"."50\n";
print CTL "background\t"."gray\n";
close CTL;
##############################################################
# system call ChimeraX
if ($stype eq "protein"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_protein.py\n";}
if ($stype eq "protprot"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_protprot.py\n";}
if ($stype eq "dna"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_dna.py\n";}	
if ($stype eq "ligand"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_ligand.py\n";}
	
###########################################################################################################
###########################################################################################################





