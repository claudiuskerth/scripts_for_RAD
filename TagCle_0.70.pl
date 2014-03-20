#!/usr/bin/env perl
use strict;
use warnings;
use String::Approx 'adist';
#################################################################################
my $version = "TagCle_0.70.pl (04/12/2013)";
my $usage = "
$version
-by Dr. Kang-Wook Kim (k.kim\@sheffield.ac.uk)

This program can delete illumina adapter and following sequences at both ends.
It can deal with the partial adapter sequences at the end of reads.
Importantly, this program uses fuzzy matcing algorithm to allow deletions, 
insertions and substitutions in the search of adapter sequence. To use fuzzy matching the Perl module String::Approx needs to be installed. 
Please install this package and add it to your Perl library path. Otherwise fuzzy matching will not be used.
Both fastq and fasta formats can be analyzed.
Both compressed and uncompressed data are recognized automatically and it allows flexible input format.

Usage: 
	(First you may want to change the permission to the program: chmod 700 remove_illumina_adapter_0.53.pl)
	./remove_illumina_adapter.pl -i input file -o output file
	(when -o is not specified it determines the output file name automatically) 
    -t adapter sequence (optional, default \'AGATCGGAAGAGCAC\')
    -f format (optional, fastq or fasta, default fastq)
    -fi, -fd, -fs specifies the level of fuzzy insertion, deletion and substitution (default fi=2, fd=2, fs=3 deactivated)
    	or just -fe for fuzziness (default fe = 3)
    -z (optional, output to compressed file)
 	-d (optional, print match stat file if y, default n)
    -me (optional, turn off fuzzy matching if added)
    -tp number of test print on STDIN to see how it is working (default 0) --> can be piped 
    
Example:
	./TagCle_0.70.pl -i input.fastq.gz (or just perl TagCle_0.70.pl -i input.fastq.gz)
	./TagCle_0.70.pl -i input.gz -o result.gz -z
	zcat test.fastq.gz | TagCle_0.70.pl
 	zcat test.fastq.gz | TagCle_0.70.pl -o output.fastq
	./TagCle_0.70.pl -i1 first-end-in.gz -i2 paired-end-in.gz -z -tp 1000 > result
	./TagCle_0.70.pl -i1 first-end-in.gz -i2 paired-end-in.gz -o1 first-end-out.gz -o2 paired-end-out.gz -z -tp 1000 > result
								.
								.
								
	cleaned_file, discarded_file, match_stat_file and log_file are produced as output
	Currently only paired-end mode will work properly.
\n";
#################################################################################
##
## To do list
## - switch to paired-end mode
## - control orphand read & stats
## - Reverse match ???
## - put shorter windows size (1~5bp) at the both ends of a sequence to solve the qualtiy drop.   


my $Infile1;
my $Infile2;
my $OUTfile1;
my $OUTfile2;
my $SE = {};
my $PE = {};
my $SE_PE_result = {};

my $compress_flag = 0;
my $PE_flag = 0;
my $reversematching_flag = 1;
my $fuzzymatching_flag = 1;
my $write_flag = "y";
my $log_write_flag = "y";
my $write_distance_flag = "n";
my $test_print_flag = 0;
my $fussiness_allowed = 3;
my $fussy_insertion = 2;
my $fuzzy_deletion = 2;
my $fuzzy_substitution = 3;
my $length_5end_match = 30;
my $Tag_match_distance_threshold = 4;
my $Tag_No_match_distance_threshold = 9;
my $min_output_length = 30;
#my $trim_length = 0;
my $score_limit = 20;
my $QCwin_len = 10;
my $Swindow_fraction = 0.04;
my $informat = "fastq";
my $outformat = "fastq";
my $SW_command = "smith-waterman-k";



# multiplexing adapter 1 AGATCGGAAGAGCACACGTCT
# multiplexing adapter 2 AGATCGGAAGAGCGTCGTGTA
#my $tag1 = "AGATCGGAAGAGCAC";
#my $tag2 = "AGATCGGAAGAGCGT";
my $tag1 = "AGATCGGAAGAGCG"; #Claudius
my $tag2 = "AGATCGGAAGAGCG"; #Claudius

my $taglength1 = length $tag1;	###
my $taglength2 = length $tag1;	###

### take command line arguments
while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-i1$/) { $Infile1 = shift @ARGV; }
	elsif ($_ =~ /^-i2$/) { $PE_flag = 1 if ($Infile2 = shift @ARGV)}
	elsif ($_ =~ /^-o1$/) { $OUTfile1 = shift @ARGV; }
	elsif ($_ =~ /^-o2$/) { $OUTfile2 = shift @ARGV; }
	elsif ($_ =~ /^-t1$/) { $tag1 = shift @ARGV; }
	elsif ($_ =~ /^-t2$/) { $tag2 = shift @ARGV; }
	elsif ($_ =~ /^-inf$/) { $informat = shift @ARGV; }
	elsif ($_ =~ /^-outf$/) { $outformat = shift @ARGV; }
	elsif ($_ =~ /^-fe$/) { $fussiness_allowed = shift @ARGV; }
	elsif ($_ =~ /^-fi$/) { $fussy_insertion = shift @ARGV; }
	elsif ($_ =~ /^-fd$/) { $fuzzy_deletion = shift @ARGV; }
	elsif ($_ =~ /^-fs$/) { $fuzzy_substitution = shift @ARGV; }
	#elsif ($_ =~ /^-tr$/) { $trim_length = shift @ARGV; }
	elsif ($_ =~ /^-tp$/) { $test_print_flag = shift @ARGV; }
	elsif ($_ =~ /^-me$/) { $fuzzymatching_flag = 0; }
	elsif ($_ =~ /^-z$/) { $compress_flag = 1; }
	elsif ($_ =~ /^-d$/) { $write_distance_flag = shift @ARGV; }
	elsif ($_ =~ /^-w$/) { $write_flag = shift @ARGV; }
	elsif ($_ =~ /^-l$/) { $log_write_flag = shift @ARGV; }
	elsif ($_ =~ /^-h$/) { die $usage; }
	elsif ($_ =~ /^-.*$/) { die $usage; }
}

eval {
	require String::Approx;
	use String::Approx qw (amatch aindex); };
die "remove_illumina_adapter.pl requires the Perl module String::Approx for fuzzy matching. Please install this package and add it to your Perl library path.
You can download it from http://search.cpan.org/~jhi/String-Approx/Approx.pm\n" if $@;

## check perl library
if ($@){
	print "remove_illumina_adapter.pl requires the Perl module String::Approx for fuzzy matching. Please install this package and add it to your Perl library path.
	Fuzzy matching will not be used this time.\n";
	$fuzzymatching_flag = 0;
}

my $revcom_tag1 = revcom($tag1);
my $revcom_tag2 = revcom($tag2);


### open input file(s) or get data from pipe
initialize_seq_hash ($Infile1, $OUTfile1, $SE);
initialize_seq_hash ($Infile2, $OUTfile2, $PE);
initialize_matchstat_hash ($SE_PE_result);

print "\n";
print "Processed SE file         :\t$Infile1\n";
print "Processed PE file         :\t$Infile2\n" if $PE_flag;
print "Adapter sequence 1        :\t$tag1\n";
print "Adapter sequence 2        :\t$tag2\n" if $PE_flag;
print "RevComp adapter sequence 1:\t$revcom_tag1\n";
print "RevComp adapter sequence 2:\t$revcom_tag2\n" if $PE_flag;
if (defined $fussiness_allowed){
	print "Fuzziness allowed         \t$fussiness_allowed\n";	
}
else{
	print "Insertion allowed         \t$fussy_insertion\n";
	print "Deletion allowed          \t$fuzzy_deletion\n";
	print "Substitution allowed      \t$fuzzy_substitution\n";
}
#print "length of 3' end trimming :\t$trim_length\n\n";



my $match_distance_fh;
if ($write_distance_flag eq "y"){
	my $mstat_file = $SE->{'match_distance_file'};
	open ($match_distance_fh, ">$mstat_file") or die $!;
}


### Main loop ###
while (defined read_next_seq($SE)){

	die $! unless (defined read_next_seq($PE));
	
	run_SWalignment($SE, $PE);

	### Foward adapter search 
	adapter_search ($tag1, $SE);
	adapter_search ($tag2, $PE);

	check_quality ($SE);
	check_quality ($PE);
	
	### stats for match mode
	$SE_PE_result->{$SE->{'match_mode'}}->{$PE->{'match_mode'}}->{'COUNT'}++;

	print $match_distance_fh "$SE->{'current_id'}\t$SE->{'match_mode'}\t$PE->{'match_mode'}\t$SE->{'SW_align_start'}\t$SE->{'SW_align_end'}\t$SE->{'SW_align_start'}\t$SE->{'SW_align_end'}\t$SE->{'SW_pct_identical'}\t$SE->{'SW_align_length'}\n" if $write_distance_flag eq "y";

	if($test_print_flag > 0){
		testprint($SE);  
		testprint($PE);  
		last if $SE->{'count_total_reads'} >= $test_print_flag;
	}

	select_output_seq($SE);
	select_output_seq($PE);

	### Select sequences to output
	if (check_minlength($SE) && check_minlength($PE)){
		write_seq($SE);
		write_seq($PE);
	}
#	else{
#		write_discarded($SE);
#		write_discarded($PE);
#	}

}

### print log files
unless ($log_write_flag eq "n") {
	print "\n--Single end reads\n";
	write_LOG($SE);
	print "\n--Paired end reads\n";
	write_LOG($PE);
}

### print match distance stats
if($write_distance_flag eq "y"){
	write_LOG_distance ($SE_PE_result, $SE);
}

close_files($SE);
close_files($PE);
exit;




###################################### Subroutines ######################################

sub testprint {
	return unless $test_print_flag;
#	my $seq_hash_ref = shift;
#	
#	print "Total reads              : $seq_hash_ref->{'count_total_reads'}\n";
#	print "seq id                   : $seq_hash_ref->{'current_id'}\n";
#	print "-------------------------\n";
#	print "SW align start           : $seq_hash_ref->{'SW_align_start'}\n";
#	print "SW align end             : $seq_hash_ref->{'SW_align_end'}\n";
#	print "SW align pcercent match  : $seq_hash_ref->{'SW_pct_identical'}\n";
#	print "SW align length          : $seq_hash_ref->{'SW_align_length'}\n";
#	print "-------------------------\n";
#	print "Adapter match mode       : $seq_hash_ref->{'match_mode'}\n";   # $seq_hash_ref->{'fuzzy_correction_log'}\n";
#	print "Adapter trimming start   : $seq_hash_ref->{'Atrim_start'}\n";
#	print "Quality trimming start   : $seq_hash_ref->{'Qtrim_start'}\n";
#	print "Final trimming start     : $seq_hash_ref->{'trim_start'}"; print "(Quality trimmed)" if $seq_hash_ref->{'Qtrimmed_flag'}==1; print "\n";
#	print "--------------------------\n";
#	print "Total matched            : $seq_hash_ref->{'count_total_matched_tags'}\n";
#	print "Exact match              : $seq_hash_ref->{'count_exact_match'}\n";
#	print "SWAlign match            : $seq_hash_ref->{'count_SW_match'}\n";
#	print "Fuzzy match              : $seq_hash_ref->{'count_fuzzy_match'}\n";
#	print "Non match                : $seq_hash_ref->{'count_non_match'}\n";
#	print "--------------------------\n";
#
#	print "<Before>\n$seq_hash_ref->{'current_seq'}\n$seq_hash_ref->{'current_qual'}\n";
#	if ($seq_hash_ref->{'Atrim_start'} > 0){
#		trim_seq($seq_hash_ref->{'Atrim_start'}-1, $seq_hash_ref);
#		print "<After adapter matching>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
#	}
#	else{
#		$seq_hash_ref->{'new_seq'} = $seq_hash_ref->{'current_seq'};
#		$seq_hash_ref->{'new_qual'} = $seq_hash_ref->{'current_qual'};
#		$seq_hash_ref->{'new_length'} = $seq_hash_ref->{'current_length'};
#		print "<After adapter matching>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
#	}
#
#
#	if (($seq_hash_ref->{'Qtrim_start'} > 0) && $seq_hash_ref->{'Qtrim_start'} < $seq_hash_ref->{'Atrim_start'}){
#		trim_seq($seq_hash_ref->{'Qtrim_start'}-1, $seq_hash_ref);
#		print "<After quality trimming>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
#	}
#	elsif($seq_hash_ref->{'Qtrim_start'} > 0 && $seq_hash_ref->{'Atrim_start'} < 0){
#		trim_seq($seq_hash_ref->{'Qtrim_start'}-1, $seq_hash_ref);
#		print "<After quality trimming>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
#	}
#	else{
#		print "<After quality trimming>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
#	}
#
#	print "------------------------------------------------------------------------------------------\n\n";

	my $seq_hash_ref = shift;
	
#	print "Total reads              : $seq_hash_ref->{'count_total_reads'}\n";
#	print "seq id                   : $seq_hash_ref->{'current_id'}\n";
#	print "-------------------------\n";
#	print "SW align start           : $seq_hash_ref->{'SW_align_start'}\n";
#	print "SW align end             : $seq_hash_ref->{'SW_align_end'}\n";
#	print "SW align pcercent match  : $seq_hash_ref->{'SW_pct_identical'}\n";
#	print "SW align length          : $seq_hash_ref->{'SW_align_length'}\n";
#	print "-------------------------\n";
#	print "Adapter match mode       : $seq_hash_ref->{'match_mode'}\n";   # $seq_hash_ref->{'fuzzy_correction_log'}\n";
#	print "Adapter trimming start   : $seq_hash_ref->{'Atrim_start'}\n";
#	print "Quality trimming start   : $seq_hash_ref->{'Qtrim_start'}\n";
#	print "Final trimming start     : $seq_hash_ref->{'trim_start'}"; print "(Quality trimmed)" if $seq_hash_ref->{'Qtrimmed_flag'}==1; print "\n";
#	print "--------------------------\n";
#	print "Total matched            : $seq_hash_ref->{'count_total_matched_tags'}\n";
#	print "Exact match              : $seq_hash_ref->{'count_exact_match'}\n";
#	print "SWAlign match            : $seq_hash_ref->{'count_SW_match'}\n";
#	print "Fuzzy match              : $seq_hash_ref->{'count_fuzzy_match'}\n";
#	print "Non match                : $seq_hash_ref->{'count_non_match'}\n";
#	print "--------------------------\n";

	if ($seq_hash_ref->{'Atrim_start'} > 0){
		print "seq id                   : $seq_hash_ref->{'current_id'}\n";
		print "SW align start           : $seq_hash_ref->{'SW_align_start'}\n";
		print "SW align end             : $seq_hash_ref->{'SW_align_end'}\n";
		print "SW align pcercent match  : $seq_hash_ref->{'SW_pct_identical'}\n";
		print "SW align length          : $seq_hash_ref->{'SW_align_length'}\n";
		print "<Before>\n$seq_hash_ref->{'current_seq'}\n$seq_hash_ref->{'current_qual'}\n";
		trim_seq($seq_hash_ref->{'Atrim_start'}-1, $seq_hash_ref);
		print "<After adapter matching>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
		print "------------------------------------------------------------------------------------------\n\n";
	}
	elsif (($seq_hash_ref->{'Qtrim_start'} > 0) && $seq_hash_ref->{'Qtrim_start'} < $seq_hash_ref->{'Atrim_start'}){
		print "seq id                   : $seq_hash_ref->{'current_id'}\n";
		print "Adapter match mode       : $seq_hash_ref->{'match_mode'}\n";   # $seq_hash_ref->{'fuzzy_correction_log'}\n";
		print "Adapter trimming start   : $seq_hash_ref->{'Atrim_start'}\n";
		print "Quality trimming start   : $seq_hash_ref->{'Qtrim_start'}\n";
		print "Final trimming start     : $seq_hash_ref->{'trim_start'}"; print "(Quality trimmed)" if $seq_hash_ref->{'Qtrimmed_flag'}==1; print "\n";
		print "<Before>\n$seq_hash_ref->{'current_seq'}\n$seq_hash_ref->{'current_qual'}\n";
		trim_seq($seq_hash_ref->{'Qtrim_start'}-1, $seq_hash_ref);
		print "<After quality trimming>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
		print "------------------------------------------------------------------------------------------\n\n";
	}
	elsif($seq_hash_ref->{'Qtrim_start'} > 0 && $seq_hash_ref->{'Atrim_start'} < 0){
		print "seq id                   : $seq_hash_ref->{'current_id'}\n";
		print "Adapter match mode       : $seq_hash_ref->{'match_mode'}\n";   # $seq_hash_ref->{'fuzzy_correction_log'}\n";
		print "Adapter trimming start   : $seq_hash_ref->{'Atrim_start'}\n";
		print "Quality trimming start   : $seq_hash_ref->{'Qtrim_start'}\n";
		print "Final trimming start     : $seq_hash_ref->{'trim_start'}"; print "(Quality trimmed)" if $seq_hash_ref->{'Qtrimmed_flag'}==1; print "\n";
		print "<Before>\n$seq_hash_ref->{'current_seq'}\n$seq_hash_ref->{'current_qual'}\n";
		trim_seq($seq_hash_ref->{'Qtrim_start'}-1, $seq_hash_ref);
		print "<After quality trimming>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
		print "------------------------------------------------------------------------------------------\n\n";
	}
	else{
		$seq_hash_ref->{'new_seq'} = $seq_hash_ref->{'current_seq'};
		$seq_hash_ref->{'new_qual'} = $seq_hash_ref->{'current_qual'};
		$seq_hash_ref->{'new_length'} = $seq_hash_ref->{'current_length'};
#		print "<After adapter matching>\n$seq_hash_ref->{'new_seq'}\n$seq_hash_ref->{'new_qual'}\n";
	}

}

sub initialize_seq_hash{
	my ($Infile, $OUTfile,$seq_hash_ref)=@_;
	my $Logfile;
	my $Discardedfile;
	my $Match_distance_file;
	my $in_fh;
	my $out_fh;
	my $log_fh;
	my $discarded_fh;
	my $match_distance_fh;

	### initialize sequence hash
	$seq_hash_ref->{'in_file'} = $Infile;
	$seq_hash_ref->{'out_file'} = "";
	$seq_hash_ref->{'log_file'} = "";
	$seq_hash_ref->{'in_FH'} = "";
	$seq_hash_ref->{'out_FH'} = "";
	$seq_hash_ref->{'log_FH'} = "";
	$seq_hash_ref->{'discarded_FH'} = "";
	$seq_hash_ref->{'match_distance_file'} = "";
	$seq_hash_ref->{'in_format'} = $informat;
	$seq_hash_ref->{'out_format'} = $outformat;
	$seq_hash_ref->{'count_total_matched_tags'} = 0;
	$seq_hash_ref->{'count_foward_tags'} = 0;
	$seq_hash_ref->{'count_reverse_tags'} = 0;
	$seq_hash_ref->{'count_start_tags'} = 0;
	$seq_hash_ref->{'count_end_tags'} = 0;
	$seq_hash_ref->{'count_3end_trimmed'} = 0;
	$seq_hash_ref->{'count_total_reads'} = 0;	
	$seq_hash_ref->{'count_exact_match'} = 0;	
	$seq_hash_ref->{'count_SW_match'} = 0;
	$seq_hash_ref->{'count_fuzzy_match'} = 0;	
	$seq_hash_ref->{'count_tail_match'} = 0;	
	$seq_hash_ref->{'count_non_match'} = 0;	
	$seq_hash_ref->{'count_fuzzy_correction'} = 0;	
	$seq_hash_ref->{'count_too_short'} = 0;	
	$seq_hash_ref->{'count_total_retained'} = 0;	
	$seq_hash_ref->{'count_total_discared'} = 0;	

	### open input file
	if (defined $Infile){
		if ($Infile =~ /.*\/(.*)\.gz$/){
			open($in_fh, "zcat $Infile |") or die $!;
			$seq_hash_ref->{'in_FH'} = $in_fh;
			$OUTfile = "cleaned_".$1 unless defined $OUTfile;
			$Logfile = "log_".$1;
			$Discardedfile = "discarded_".$1;
			$seq_hash_ref->{'match_distance_file'} = "match_stat_".$1;
		}
		else{
			open($in_fh, $Infile) or die $!;
			$seq_hash_ref->{'in_FH'} = $in_fh;
			my ($Infile_stub) = $Infile =~ /.*\/(.*)$/; #Claudius
			$OUTfile = "cleaned_".$Infile_stub unless defined $OUTfile;
			$Logfile = "log_".$Infile;
			$Discardedfile = "discarded_".$Infile;
			$seq_hash_ref->{'match_distance_file'} = "match_stat_".$Infile;
		}
	}
	else{
		if (-t STDIN and not defined $Infile){  #input attached to terminal
			die $usage;
		}
		else{  #input from pipe
			$seq_hash_ref->{'in_FH'} = *STDIN;
			$seq_hash_ref->{'in_file'} = "STDIN";
			my $Start_time = gmtime($^T);
			$Start_time =~ s/\s/\_/g;
			$Infile = "STDIN";
			$OUTfile = "cleaned_STDIN\_$Start_time" unless defined $OUTfile;
			$Logfile = "log_STDIN\_$Start_time";
			$Discardedfile = "discarded_STDIN\_$Start_time";
			$seq_hash_ref->{'match_distance_file'} = "match_distance_STDIN\_$Start_time";
		}
	}
	
	### open output file
	unless ($write_flag eq "n"){
		if ($compress_flag == 0){
			open ($out_fh, ">$OUTfile") or die $!, $OUTfile, "\n";
			$seq_hash_ref->{'out_FH'} = $out_fh;

			#Claudius
#			open ($discarded_fh, ">$Discardedfile") or die $!;
#			$seq_hash_ref->{'discarded_FH'} = $discarded_fh;
		}
		else{
			open $out_fh, "| gzip -c > $OUTfile\.gz" or die $!;
			$seq_hash_ref->{'out_FH'} = $out_fh;

			#Claudius
#			open $discarded_fh, "| gzip -c > $Discardedfile\.gz" or die $!;
#			$seq_hash_ref->{'discarded_FH'} = $discarded_fh;
		}
	}

	#Claudius
	### open log file
#	unless ($log_write_flag eq "n"){
#		open($log_fh,">$Logfile") or die $!;
#		$seq_hash_ref->{'log_FH'} = $log_fh;
#	}
}

sub initialize_matchstat_hash {
	my ($stathash_ref) = @_;
	my @matchmode_array = ('exact','SWAlign','fuzzy','non');
	my ($i, $j, $k);
	foreach $i (@matchmode_array) {
		foreach $j (@matchmode_array){
			$stathash_ref->{"$i"}->{"$j"}->{'COUNT'} = 0;
			for ($k=0; $k <= $length_5end_match; $k++){
				$stathash_ref->{"$i"}->{"$j"}->{'SE_distance'}->{$k}->{'COUNT'} = 0;
				$stathash_ref->{"$i"}->{"$j"}->{'PE_distance'}->{$k}->{'COUNT'} = 0;
			} 
		}
	}
}

sub read_next_seq {
	my $seq_hash_ref = shift;
	my $ifh = $seq_hash_ref->{'in_FH'};

	$seq_hash_ref->{'current_id'} = "";
	$seq_hash_ref->{'current_seq'} = "";
	$seq_hash_ref->{'current_qual_id'} = "";
	$seq_hash_ref->{'current_qual'} = "";
	$seq_hash_ref->{'current_length'} = 0;
	$seq_hash_ref->{'new_id'} = "";
	$seq_hash_ref->{'new_seq'} = "";
	$seq_hash_ref->{'new_qual_id'} = "";
	$seq_hash_ref->{'new_qual'} = "";
	$seq_hash_ref->{'new_length'} = 0;
	$seq_hash_ref->{'trim_flag'} = 0;
	$seq_hash_ref->{'match_mode'} = "";
	$seq_hash_ref->{'Atrim_start'} = -1;
	$seq_hash_ref->{'Qtrim_start'} = -1;
	$seq_hash_ref->{'trim_start'} = -1;
	$seq_hash_ref->{'Qtrimmed_flag'} = 0;

	unless ($seq_hash_ref->{'current_id'} = <$ifh>){
		return undef;	
	}
	chomp $seq_hash_ref->{'current_id'};
	$seq_hash_ref->{'current_id'} =~ s/\s/\_/g;

	$seq_hash_ref->{'current_seq'} = <$ifh>;
	chomp $seq_hash_ref->{'current_seq'};
	$seq_hash_ref->{'current_length'} = length $seq_hash_ref->{'current_seq'};

	if ($seq_hash_ref->{'in_format'} eq "fastq"){
		$seq_hash_ref->{'current_qual_id'} = <$ifh>;
		chomp $seq_hash_ref->{'current_qual_id'};
		$seq_hash_ref->{'current_qual'} = <$ifh>;
		chomp $seq_hash_ref->{'current_qual'}
	}
	$seq_hash_ref->{'count_total_reads'}++;
	return 1;
}

##################################################################


sub run_SWalignment {
	my ($seq1_hash_ref, $seq2_hash_ref)=@_;
	my $seq1 = $seq1_hash_ref->{'current_seq'};
	my $seq2 = $seq2_hash_ref->{'current_seq'};
	$seq2 = revcom($seq2);

 	my ($seq1_start,$seq1_end,$seq2_start,$seq2_end, $pctidentical,$align_length) = `$SW_command $seq1 $seq2`;
	chomp ($seq1_start,$seq1_end,$seq2_start,$seq2_end,$pctidentical,$align_length);

	$seq1_hash_ref->{'SW_align_start'} = $seq1_start;
	$seq1_hash_ref->{'SW_align_end'} = $seq1_end;
	$seq2_hash_ref->{'SW_align_start'} = $seq2_hash_ref->{'current_length'} - $seq2_end + 1;
	$seq2_hash_ref->{'SW_align_end'} = $seq2_hash_ref->{'current_length'} - $seq2_start + 1;
	$seq1_hash_ref->{'SW_pct_identical'}=$pctidentical;
	$seq1_hash_ref->{'SW_align_length'}=$align_length;
	$seq2_hash_ref->{'SW_pct_identical'}=$pctidentical;
	$seq2_hash_ref->{'SW_align_length'}=$align_length;

}


sub adapter_search {
	my ($tag_ref, $seq_hash_ref)=@_;
	
	if(defined exact_matching(\$tag_ref,$seq_hash_ref)){
	}
	elsif (defined SW_matching($seq_hash_ref) && $PE_flag == 1){
	}
#	elsif(defined fuzzy_matching(\$tag_ref,$seq_hash_ref) && $fuzzymatching_flag == 1){
#	}
	elsif( $fuzzymatching_flag == 1 && defined fuzzy_matching(\$tag_ref,$seq_hash_ref) ){ #Claudius: left precedence of '&&'
	}
	else{
		no_matching($seq_hash_ref);
	}
}



sub SW_matching {
	my ($seq_hash_ref) = @_;

	if ($seq_hash_ref->{'SW_pct_identical'} > 75 && $seq_hash_ref->{'SW_align_length'} > 30){
		if ($seq_hash_ref->{'SW_align_start'} < 4){ #Claudius: ???
			$seq_hash_ref->{'Atrim_start'} = $seq_hash_ref->{'SW_align_end'} + 1;
			$seq_hash_ref->{'trim_start'} = $seq_hash_ref->{'Atrim_start'};
#			$seq_hash_ref->{'new_length'} = $seq_hash_ref->{'SW_align_end'};
			$seq_hash_ref->{'count_total_matched_tags'}++;
			$seq_hash_ref->{'count_foward_tags'}++;
			$seq_hash_ref->{'count_SW_match'}++;
			$seq_hash_ref->{'trim_flag'} = 1;
			$seq_hash_ref->{'match_mode'} = "SWAlign";
			return 1;
		}
		else{
			return undef;
		}
	}
	else{
		return undef;
	}
}


sub exact_matching {
	my ($tag_ref, $seq_hash_ref)=@_;
	if($seq_hash_ref->{'current_seq'} =~ /$$tag_ref/i){
		$seq_hash_ref->{'Atrim_start'} = $-[0] + 1;
		$seq_hash_ref->{'trim_start'} = $seq_hash_ref->{'Atrim_start'};
		$seq_hash_ref->{'count_total_matched_tags'}++;
		$seq_hash_ref->{'count_foward_tags'}++;
		$seq_hash_ref->{'count_exact_match'}++;
		$seq_hash_ref->{'trim_flag'} = 1;
		$seq_hash_ref->{'match_mode'} = "exact";
		
		return 1;
	}
	else{
		return undef;
	}
}

sub fuzzy_matching {
	my ($tag_ref, $seq_hash_ref)=@_;
	my $index;
	
	if (defined $fussiness_allowed){ #Claudius: initialized on line 70, so always defined
		$index = aindex($$tag_ref,["i","$fussiness_allowed"],$seq_hash_ref->{'current_seq'});
	}
	else{
		$index = aindex($$tag_ref,["i","I$fussy_insertion","D$fuzzy_deletion","S$fuzzy_substitution"],$seq_hash_ref->{'current_seq'});
	}
	
	if($index > -1){

		$seq_hash_ref->{'Atrim_start'} = $index +1;
		$seq_hash_ref->{'trim_start'} = $index +1;
		$seq_hash_ref->{'count_total_matched_tags'}++;
		$seq_hash_ref->{'count_foward_tags'}++;
		$seq_hash_ref->{'count_fuzzy_match'}++;
		$seq_hash_ref->{'trim_flag'} = 1;
		$seq_hash_ref->{'match_mode'} = "fuzzy";
		return 1;
	}
	else {
		return undef;
	}
}



sub no_matching {
	my ($seq_hash_ref)=@_;
	$seq_hash_ref->{'count_non_match'}++;
	$seq_hash_ref->{'match_mode'} = "non";
	$seq_hash_ref->{'new_length'} = $seq_hash_ref->{'current_length'};
}


sub trim_seq {
	my ($index, $seq_hash_ref)=@_;
	$seq_hash_ref->{'new_seq'} = substr ($seq_hash_ref->{'current_seq'}, 0, $index);
	$seq_hash_ref->{'new_length'} = length $seq_hash_ref->{'new_seq'};
	#$seq_hash_ref->{'Atrim_start'} = $index +1;
	if ($seq_hash_ref->{'in_format'} eq "fastq"){
		$seq_hash_ref->{'new_qual'} = substr ($seq_hash_ref->{'current_qual'}, 0, $index); ## fq
	}
}

sub check_minlength {
	my ($seq_hash_ref)=@_;

	if ($seq_hash_ref->{'new_length'} >= $min_output_length){
		return 1;
	}
	else{
		$seq_hash_ref->{'count_too_short'}++;	
		write_LOG($seq_hash_ref,"Record $seq_hash_ref->{'current_id'} is shorter than $min_output_length (bp) after trimming!\n");
		return 0;
	}
}




sub select_output_seq {
	my $seq_hash_ref = shift;

	if ($seq_hash_ref->{'Atrim_start'} > 0 && $seq_hash_ref->{'Qtrim_start'} > 0){
		if ($seq_hash_ref->{'Atrim_start'} >= $seq_hash_ref->{'Qtrim_start'}){
			trim_seq($seq_hash_ref->{'Qtrim_start'}-1, $seq_hash_ref);
		} 
		else{
			trim_seq($seq_hash_ref->{'Atrim_start'}-1, $seq_hash_ref);
		}
	}
	elsif($seq_hash_ref->{'Atrim_start'} > 0 && $seq_hash_ref->{'Qtrim_start'} < 0){
			trim_seq($seq_hash_ref->{'Atrim_start'}-1, $seq_hash_ref);
	}
	elsif($seq_hash_ref->{'Atrim_start'} < 0 && $seq_hash_ref->{'Qtrim_start'} > 0){
			trim_seq($seq_hash_ref->{'Qtrim_start'}-1, $seq_hash_ref);
	}
	else {
		$seq_hash_ref->{'new_seq'} = $seq_hash_ref->{'current_seq'};
		$seq_hash_ref->{'new_qual'} = $seq_hash_ref->{'current_qual'};
		$seq_hash_ref->{'new_length'} = $seq_hash_ref->{'current_length'};
	}
}


sub write_seq {
	my $seq_hash_ref = shift;
	my $out_fh = $seq_hash_ref->{'out_FH'};
	$seq_hash_ref->{'count_total_retained'}++; 
	print $out_fh "$seq_hash_ref->{'current_id'}\n"; ###---> change?
	print $out_fh "$seq_hash_ref->{'new_seq'}\n";
	if ($seq_hash_ref->{'out_format'} eq "fastq"){
		print $out_fh "$seq_hash_ref->{'current_qual_id'}\n"; ###---> change?
		print $out_fh "$seq_hash_ref->{'new_qual'}\n"; ##fq
	}
}

sub write_discarded {
	my $seq_hash_ref = shift;
	my $discarded_fh = $seq_hash_ref->{'discarded_FH'};
	$seq_hash_ref->{'count_total_discared'}++;	

	print $discarded_fh "$seq_hash_ref->{'current_id'}\n"; ###---> change?
	print $discarded_fh "$seq_hash_ref->{'current_seq'}\n";
	if ($seq_hash_ref->{'out_format'} eq "fastq"){
		print $discarded_fh "$seq_hash_ref->{'current_qual_id'}\n"; ###---> change?
		print $discarded_fh "$seq_hash_ref->{'current_qual'}\n"; ##fq
	}
}



sub write_LOG {
	my ($seq_hash_ref, $msg) = @_;
	my $log_fh = $seq_hash_ref->{'log_FH'};
	if (defined $msg){
#		print $log_fh $msg; # Claudius
	}
	else{
		my $Start_time = gmtime($^T);
		my $End_time = gmtime (time());
		my $elapsed_time = time-$^T;
		my ($sec,$min,$hour,$mday) = gmtime($elapsed_time);
		$mday = $mday-1;

		
		
		
		my $percent_retained = get_percent ($seq_hash_ref->{'count_total_retained'}, $seq_hash_ref->{'count_total_reads'});
		my $percent_discarded = get_percent ($seq_hash_ref->{'count_total_discared'}, $seq_hash_ref->{'count_total_reads'});
		my $percent_total_matched = get_percent ($seq_hash_ref->{'count_total_matched_tags'}, $seq_hash_ref->{'count_total_reads'});

		#Claudius
#		print $log_fh "$version\n";
#		print $log_fh "\nProcessing\t$seq_hash_ref->{'in_file'}\n";
#		print $log_fh "Job started at: $Start_time\n";
#		print $log_fh "Job finished at: $End_time\n";
#		print $log_fh "Elapsed time: $mday day(s) $hour:$min:$sec\n\n";
#		print $log_fh "Adapter sequence 1        \t$tag1\n";
#		print $log_fh "Adapter sequence 2        \t$tag2\n";
#		print $log_fh "RevComp adapter sequence 1\t$revcom_tag1\n";
#		print $log_fh "RevComp adapter sequence 2\t$revcom_tag2\n";
#
#		if (defined $fussiness_allowed){
#			print $log_fh "Fuzziness allowed         \t$fussiness_allowed\n";	
#		}
#		else{
#			print $log_fh "Insertion allowed         \t$fussy_insertion\n";
#			print $log_fh "Deletion allowed          \t$fuzzy_deletion\n";
#			print $log_fh "Substitution allowed      \t$fuzzy_substitution\n";
#		}
#
#
#		#print $log_fh "length of 3' end trimming \t$trim_length\n";
#		print $log_fh "Total reads               \t$seq_hash_ref->{'count_total_reads'}\n";
#		print $log_fh "Total matched             \t$seq_hash_ref->{'count_total_matched_tags'}\t($percent_total_matched %)\n";
#		print $log_fh "Exact match               \t$seq_hash_ref->{'count_exact_match'}\n";
#		print $log_fh "SWAlign match             \t$seq_hash_ref->{'count_SW_match'}\n";
#		print $log_fh "Fuzzy match               \t$seq_hash_ref->{'count_fuzzy_match'}\n";
#		print $log_fh "Non-match                 \t$seq_hash_ref->{'count_non_match'}\n";
#		print $log_fh "Too_short after trimming  \t$seq_hash_ref->{'count_too_short'}\n";
#		print $log_fh "Total retained            \t$seq_hash_ref->{'count_total_retained'}\t($percent_retained %)\n";
#		print $log_fh "Total discared            \t$seq_hash_ref->{'count_total_discared'}\t($percent_discarded %)\n";
	
		print "Elapsed time              :\t$mday day(s) $hour:$min:$sec\n";
		print "Total reads               :\t$seq_hash_ref->{'count_total_reads'}\n";
		print "Total matched             :\t$seq_hash_ref->{'count_total_matched_tags'}\t($percent_total_matched %)\n";
		print "Exact match               :\t$seq_hash_ref->{'count_exact_match'}\n";
		print "SWAlign match             :\t$seq_hash_ref->{'count_SW_match'}\n";
		print "Fuzzy match               :\t$seq_hash_ref->{'count_fuzzy_match'}\n";
		print "Non-match                 :\t$seq_hash_ref->{'count_non_match'}\n";
		print "Too_short after trimming  :\t$seq_hash_ref->{'count_too_short'}\n";
		print "Total retained            :\t$seq_hash_ref->{'count_total_retained'} ($percent_retained %)\n";
		print "Total discared            :\t$seq_hash_ref->{'count_total_discared'} ($percent_discarded %)\n";
	}
}

sub write_LOG_distance {
	my ($mstat_hash_ref, $seq_hash_ref)=@_;
	my $mstat_file = $seq_hash_ref->{'match_distance_file'};
	open (my $match_distance_fh, ">$mstat_file") or die $!;
	if ($write_distance_flag eq "y"){
		### print match mode stats
		print $match_distance_fh "Read counts per matching mode\n";
		print $match_distance_fh "SE_mode\tPE_mode\tCounts\n";
		for my $SE_mode_temp (sort {$a cmp $b} keys %{$mstat_hash_ref}){
			for my $PEmode_temp (sort {$a cmp $b} keys %{$mstat_hash_ref->{$SE_mode_temp}}){
				my $Matchmode_Count = $mstat_hash_ref->{$SE_mode_temp}->{$PEmode_temp}->{'COUNT'};
				print $match_distance_fh "$SE_mode_temp\t$PEmode_temp\t$Matchmode_Count\n";
			}
		}
		
		### print match mode vs distance stats
		print $match_distance_fh "\nMatching distance between $length_5end_match bp of 5' end to paired sequence\n";
		print $match_distance_fh "SE_mode\tPE_mode\tEdit_distance\tSE2PE\t PE2SE\n";
		for my $SE_mode_temp (sort {$a cmp $b} keys %{$mstat_hash_ref}){
			for my $PEmode_temp (sort {$a cmp $b} keys %{$mstat_hash_ref->{$SE_mode_temp}}){
				for my $Edistance_temp (sort {$a <=> $b} keys %{$mstat_hash_ref->{$SE_mode_temp}->{$PEmode_temp}->{'SE_distance'}}){
					my $SEMatchmode_Count_perdistance = $mstat_hash_ref->{$SE_mode_temp}->{$PEmode_temp}->{'SE_distance'}->{$Edistance_temp}->{'COUNT'};
					my $PEMatchmode_Count_perdistance = $mstat_hash_ref->{$SE_mode_temp}->{$PEmode_temp}->{'PE_distance'}->{$Edistance_temp}->{'COUNT'};
					print $match_distance_fh "$SE_mode_temp\t$PEmode_temp\t$Edistance_temp\t$SEMatchmode_Count_perdistance\t$PEMatchmode_Count_perdistance\n";
				}
			}
		}
	}
	close $match_distance_fh;
}

sub close_files {
	my $seq_hash_ref = shift;
	unless ($seq_hash_ref->{'in_file'} eq "STDIN"){
		close $seq_hash_ref->{'in_FH'};	
	}

	close $seq_hash_ref->{'out_FH'} unless $write_flag eq "n";
	close $seq_hash_ref->{'log_FH'} unless $log_write_flag eq "n";		
}

sub revcom {
	my $inseq = shift;
	my $revcom_seq = reverse $inseq;
	$revcom_seq =~ tr/ATGCatgc/TACGtacg/;
	return $revcom_seq;
}

sub get_percent {
	my ($numerator, $denominator) = @_;
	my $answer = $numerator / $denominator * 100;
	$answer = sprintf "%.1f",$answer;
	return $answer;
}



sub check_quality {
	my ($seq_hash_ref)=@_;

    # Illumina 1.3+ encodes phred scores between ASCII values 64 (0 quality) and 104 (40 quality)
    #
    #   @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
    #   |         |         |         |         |
    #  64        74        84        94       104
    #   0        10(90%)   20(99%)   30(99.9%) 40(99.99%)
    #
    # Illumina 1.8+ encodes phred scores between ASCII values 33 (0 quality) and 74 (41 quality)
    # 
    #   !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
    #   |         |         |         |         |    
    #  33        43        53        63        73 
    #
    #   0        10(90%)   20(99%)   30(99.9%) 40(99.99%)
    
    my ($j, $k, @scores, @int_scores, $score, $len, $stop_pos, $mean, $int);

    $len      = $seq_hash_ref->{'current_length'};
    
    $stop_pos = $len - $QCwin_len;
    @scores   = split(//, $seq_hash_ref->{'current_qual'});

    foreach $j (0 .. $len - 1) {
        $int_scores[$j] = ord($scores[$j]) - 33;
    }
    foreach $j (0 .. $stop_pos) {
		$mean = 0;
	
		foreach $k ($j .. $j + $QCwin_len - 1) {
			$mean += $int_scores[$k];
		}
		$mean = $mean / $QCwin_len;
		if ($mean < $score_limit) {
			$seq_hash_ref->{'Qtrim_start'} = $j;
			$seq_hash_ref->{'trim_flag'} = 1;
			last;
		}
    }
    
    if($seq_hash_ref->{'trim_start'} > 0){
    	if ($seq_hash_ref->{'trim_start'} > $seq_hash_ref->{'Qtrim_start'} && $seq_hash_ref->{'Qtrim_start'} > 0){	
			$seq_hash_ref->{'trim_start'} = $seq_hash_ref->{'Qtrim_start'};
			$seq_hash_ref->{'Qtrimmed_flag'} = 1;
    	}
    }
	else{
		if ($seq_hash_ref->{'Qtrim_start'} > 0){
			$seq_hash_ref->{'trim_start'} = $seq_hash_ref->{'Qtrim_start'};
			$seq_hash_ref->{'Qtrimmed_flag'} = 1;
		}
	}


    return 0;
}


exit;
