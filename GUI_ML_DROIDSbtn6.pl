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
		if ($header eq "chimerax_path"){$chimerax_path = $path;}
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




########################################################################################
######################     SUBROUTINES     #############################################
########################################################################################

#sub image2 {

#print "\n THIS FEATURE IS UNDER RECONSTRUCTION....try again later\n\n"; exit;

print "Enter chain ID where color map should start (default = A)\n";
   print "i.e. color map only to homologous parts of bound and unbound structure\n";
   $chainMAP = <STDIN>;
   chop($chainMAP);
   if ($chainMAP eq ''){$chainMAP = 'A';}

if($stype eq 'dna' || $stype eq 'ligand' || $stype eq 'protprot'){
   
   print "Enter number position of N terminal on this chain (default = 0)\n";
   $offset = <STDIN>;
   chop($offset);
   $offset = $offset-2;
   if ($offset eq ''){$offset = 0;}
   }

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
      #if ($p == 1){next;}
      #if ($p == 2){next;}
      my $MUTrow = $MUT[$p];
      my @MUTrow = split (/\s+/, $MUTrow);
	 $fileIDq = $MUTrow[0];
      push (@variants, $fileIDq);
      } #end for loop
close MUT;
print "\n variants are @variants\n\n";
sleep(1);

# prompt user - choose best learning model to display
sleep(1);print "SELECT VARIANT TO DISPLAY \n\n";
my $varselect = <STDIN>;
chop($varselect);

# select variant to render

$variantID = $varselect;

sleep(1);print "\nSELECT INTERACTION TYPE (1=protein only | 2=protein-protein | 3=DNA-protein | 4=protein-ligand)\n\n";
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
my $orig_queryID = $queryID;  # create tag for original query in training set
my $queryID = "$queryID"."_1"; # set query to first copy for this subroutine
print "queryID label is adjusted to "."$queryID\n\n\n";

# create mask for non zero impacts on dynamics
open(OUT, ">"."./testingData_$queryID/impact_vert_$variantID.txt")||die "could not create mask file /testingData_$queryID/impact_vert_$variantID.txt\n";  
open(IN, "<"."./testingData_$queryID/impact_$variantID.txt")||die "could not open mask file /testingData_$variantID/impact_$variantID.txt\n";  
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
 
# make relative temperature attribute files on impacted areas for image rendering

  open(OUT, ">"."./attribute_files/drmsfATTR_$refID"."_mask".".dat")||die "could not create ATTR time series file\n";  
  print OUT "recipient: residues\n";
  print OUT "attribute: dRMSF\n";
  print OUT "\n";
  @impacts = ();
  for (my $a = 0; $a < $lengthID+$offset+1; $a++){
   if ($a eq '' || $a <= $offset){next;}
   open(IN, "<"."./testingData_$queryID/impact_vert_$variantID.txt")||die "could not open mask file /testingData_$queryID/adj_vertpvals_$queryID.txt\n";  
    my @IN = <IN>;
   $INrow = $IN[$a-$offset];
    @INrow = split(/\s+/, $INrow);
    $maskvalue = $INrow[0];
     
     #print "$f\t"."$a\t"."$frame\t"."$value\n";
     $pos = $a+1;
     if ($maskvalue > 0){print OUT "\t:"."$pos\t"."$maskvalue\n";}
     #if ($f == $frame && $maskvalue > 0){print OUT "\t:"."$pos\t"."1\n";} # test mask positioning
     if ($maskvalue == 0){print OUT "\t:"."$pos\t"."0.00\n";}
     push (@impacts, $maskvalue); 
   close IN;
 }
close OUT;

print "\n class attribute files for all frames are created\n";
sleep(1);

# render dynamic movies for local flux difference from MD training run averages
print("close ChimeraX window when done\n\n");
if ($homology eq "loose"){$mutType = "gray50";}
if ($homology eq "strict"){$mutType = "tan";}
if ($colorScheme eq "c1" ){$colorType = "wr";}
if ($colorScheme eq "c2" ){$colorType = "wr";}
$attr = "dRMSF";

$statSCORE = new Statistics::Descriptive::Full; # impact map
$statSCORE->add_data (@impacts);
$min_impact = $statSCORE->min();
$max_impact = $statSCORE->max();
$min_val = $min_impact;
$max_val = $max_impact;
if($min_val eq ''){$min_val = 0;}
#print(@impacts);
print "\n";
print "min="."$min_val\t"."max="."$max_val\n";

###############################################################
print("Preparing display...\n");
print("close ChimeraX program to exit\n\n");
# copy visualization support files to folder
mkdir("ChimeraXvis");
copy("$queryID"."REDUCED.pdb", "ChimeraXvis/query.pdb");
copy("$refID"."REDUCED.pdb", "ChimeraXvis/reference.pdb");
copy("./attribute_files/drmsfATTR_$refID"."_mask".".dat", "ChimeraXvis/attributeDRMSF.dat");
# create control file for visualization
open(CTL, ">"."ChimeraXvis.ctl");
print CTL "model\t"."#1\n";
print CTL "chain\t"."/$chainMAP\n";
print CTL "structure\t"."ChimeraXvis/reference.pdb\n";
print CTL "structureADD\t"."ChimeraXvis/query.pdb\n";
print CTL "attr_file\t"."ChimeraXvis/attributeDRMSF.dat\n";
print CTL "length\t"."$lengthID\n";
print CTL "attr\t"."dRMSF\n";
print CTL "minval\t"."$min_val\n";
print CTL "maxval\t"."$max_val\n";
print CTL "palette\t"."YlOrRd\n";
print CTL "lighting\t"."simple\n";
print CTL "transparency\t"."50\n";
print CTL "background\t"."black\n";
close CTL;
##############################################################
##############################################################
# system call ChimeraX
if ($stype eq "protein"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_protein.py\n";}
if ($stype eq "protprot"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_protprot.py\n";}
if ($stype eq "dna"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_dna.py\n";}	
if ($stype eq "ligand"){system "$chimerax_path"."ChimeraX color_by_attr_chimerax_ligand.py\n";}

###########################################################################################################
###########################################################################################################





