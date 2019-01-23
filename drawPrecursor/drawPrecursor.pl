#!/usr/bin/perl

#package SrnaTools::Module::HpTool ;
#use base SrnaTools::Module ;
use strict ;
use warnings;
use File::Temp qw( tempfile tempdir );
use IO::String ;
use Cwd;
use File::Copy;


sub run{
# Commands 
my $RNAfold_cmd = 'RNAfold' ; 
my $RNAplot_cmd = 'RNAplot' ;
my $ps2pdf_cmd = 'ps2pdfwr';#'ps2pdf' ; 
my $ps2eps_cmd = 'ps2eps' ;
my $working_dir = getcwd;
my $cwd = getcwd;
($cwd) =($cwd=~/^(.+)$/);
# rgb colors for highlighting short seqs
my @rgb_colors = ('255 40 40', 
                  '0 191 255', 
                  '0 0 255', 
                  '0 255 255', 
                  '0 255 0', 
                  '255 255 0', 
                  '255 127 0',
                  '255 184 235',
                  '181 255 184',
                  '199 255 176',
                  '127 62 62',
                  '255 227 99',
                  '127 0 0', 
                  '127 0 127', 
                  '0 0 127',
                  '0 127 127',
                  '0 127 0',
                  '130 127 0',
                  '153 153 153',
                  '204 204 204') ;

  my $long_seq =  $_[0];#"TGTACTGCCCTTTCTCCATCCCCCAAATCTTTTGGAGTTTTTAACACTATAAATTGAGATACAGATGGAGATTTCTTGAGGCAGGGAGAGGAGGTCAGTGGCGGAGCTTGGAGATATCAGTAGCAGTGGAAGGGGCATGCAGAGGAGATTATATATGTTGATATGCTTCCTATGCTTCCTCTCTCCTCTGCCTGCCCCATCCACTCCTGCTGTTATCCCCTTCACGCGTCATACTGCGGATTAATCCCGTGTCCTCCTATATTTTTTTTCCAG";
   
  my $short_seq = $_[1];#"TGGAAGGGGCATGCAGAGGAGA";
  my $file = $short_seq;
  open (FH, "< $file") or die "Can't open $file for read: $!";
  my @lines;
  while (<FH>) {
    push (@lines, $_);
  }

  my $j=0;
  ### get match positions of short sequences and
  # construct postscript macro for RNAplot
  # in the format 
  # "start_pos end_pos 8 color omark" (8 is the thickness of the line)
  # start the macro by putting a circle around the
  # 5' end
  my $rnaplot_ps_macro = '1 cmark ' ;
  
  my $i = 0 ;
  foreach my $short_seq (@lines){
    my @positions = get_pos(\$long_seq,\$short_seq) ;
    foreach (@positions) {
     $rnaplot_ps_macro .= $_ .' 8 ' . convert_rgb_to_ps_format($rgb_colors[$i]) . ' omark ';
    }
    ++$i;
    last if $i > 19 ;
  }
 
### run RNAfold
  my $RNAfold_out_file = $_[2].'.RNAfold_out';#.'/data/RNAfold_out';
  my $cmd = "echo '$long_seq' | $RNAfold_cmd > $RNAfold_out_file" ; # run RNAfold
  my $out = `$cmd` ;
  
  
  ### run RNAplot
  # RNAplot always generates its output in a file named after the input sequence
  # or (if not fasta) "rna.ps". Therefore we have to generate another directory
  my $RNAplot_out_file = 'rna.ps' ; # default RNAplot file name
  my $RNAplot_ps_final = $_[2].'_rna.ps_final' ; # file after addition of label for pdf 

  # make temp dir
  my $RNAplot_dir = $working_dir.'/'.$_[2].'_RNAplot_out';#.'/data/RNAplot_dir';
  
  mkdir $RNAplot_dir;
  move($RNAfold_out_file,$RNAplot_dir);
  #copy($RNAfold_out_file,$RNAplot_dir);
  chmod 0777, $RNAplot_dir ; 
  chdir $RNAplot_dir ; 

  # RNAplot command
  $cmd = "$RNAplot_cmd --pre \"$rnaplot_ps_macro\" < $RNAfold_out_file " ; 
  $out = `$cmd` ;
  
  $rnaplot_ps_macro = '' ; # reset to empty for the next macro coomand

  ### modify and convert images
  # We add a label to the postscript file, then convert it to pdf
  # and also generate the jpg file to print to the result page
  
  # add label to postscript file
  # open the RNAplot outputfile in its temp dir and another file for writing the modified ps
  my $label ="Secondary structure for " .$_[2] ;

  open PS , '<', $RNAplot_dir .'/'. $RNAplot_out_file;
  open PS_FINAL, '>' , $RNAplot_dir . '/' . $RNAplot_ps_final;

  # add label to ps file
  # this is a bit of a hack:
  # the ps code is not inserted into the ideal position in the file
  # but by inserting just before "%%BeginProlog" we can print the label before the coordinate system is translated for plotting the RNA
  # post-script stack used:
  # Helvetica findfont -> change font
  # 12 scalefont -> set font size
  # 100 10 moveto -> set cursor position (bottom of page)
  # setfont
  # $label show -> print the label
  while (<PS>) {
    s/%%BeginProlog/\/Helvetica findfont\n12 scalefont\n100 10 moveto\nsetfont\n\($label\) show\n%%BeginProlog/ ;
    print PS_FINAL;
  }
  close PS ;
  close PS_FINAL ;

  unlink $working_dir.'/'.$RNAplot_out_file; # remove rna.ps file


  # convert to PDF and write to results file
  my $pdf_file = $_[2].'_Structure_plot.pdf';

  $cmd = $ps2pdf_cmd . ' '. $RNAplot_dir . '/' . $RNAplot_ps_final . ' ' . $pdf_file ;
  chmod 0666, $pdf_file ; 
  my $output = `$cmd` ;
  if ($output) {
    print "There was a problem running ps2pdf" ;
  }
  
  
 #chdir $working_dir;  
  
  

# get match positions of short sequence on long sequence
sub get_pos{
  my ($long_seq_ref, $short_seq_ref) = @_ ;
  my $seq;
  $seq = $$short_seq_ref;
  chomp $seq;
  $seq=~ s/\R//g;

  my @pos = () ;
  my $long_seq_len = length($$long_seq_ref) ;
  my $short_seq_len = length($seq) ;
  for (my $i = 1 ; $i <= $long_seq_len - $short_seq_len + 1; ++$i){#$long_seq_len - $short_seq_len + 1 ; ++$i) {
    my $sub=substr($$long_seq_ref, $i - 1, $short_seq_len);
    if (substr($$long_seq_ref, $i - 1, $short_seq_len) eq $seq) { # we have a match
    my $sub=substr($$long_seq_ref, $i - 1, $short_seq_len);
    my $start = $i ;
    my $stop = $i + $short_seq_len - 1 ;
    push @pos, $start . ' ' . $stop ;
    }
  }
  return @pos ;
}







# Convert space delimited RGB into notation for RNAplot macros
# # where 255 = 1
 sub convert_rgb_to_ps_format {
   my $input_rgb = shift ; # in format R G B
   my ($r,$g,$b)=($input_rgb=~/(\d+)/g) ;
   $r = $r / 255 ;
   $g = $g / 255 ;
   $b = $b / 255 ;
   return "$r $g $b" ;
 }





} # end of subroutine run

#Functional Call

my $long_SEQ=$ARGV[1];
my $short_SEQ=$ARGV[2];
my $miR=$ARGV[0];

run($long_SEQ,$short_SEQ,$miR);


