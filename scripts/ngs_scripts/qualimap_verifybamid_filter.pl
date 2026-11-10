#!/usr/bin/perl
# Source: from the upstream variant-analysis pipeline scripts used for the
# China 100K Genomes Project cohort (internal deployment).
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );

my ($qcdir,$qualimapfile1,$qualimapfile2,$verifyfile,$outfile,$configure,$software);

GetOptions
(
		"fastqcdir=s"		=> \$qcdir,
		"qualimap1=s"		=> \$qualimapfile1,
		"qualimap2=s"		=> \$qualimapfile2,
		"verifybamid=s"				=> \$verifyfile,
		"tool=s"	=> \$software,
		"conf=s"	=> \$configure,
		"out=s"		=>	\$outfile
);


my $usage="
	perl $0  -fastqcdir <sample fastqc dir> -qualimap1 <genome_results.txt> -qualimap2 <Qualimap report_ BAM QC.html>  -verifybamid  <verifybamid.result>  -tool < qualimap || verifybamid || both>  -conf <configure.v1> -out <outfile>  

	-fastqcdir							the path of sample fastqc dir
	-qualimap1							the path of qualimap result file: genome_results.txt
	-qualimap2							the path of qualimap result file: Qualimap report_ BAM QC.html
	-verifybamid						the path of verifybamid result file
	-tool								fastqc || qualimap|| verifybamid || both
	-conf								configure.txt
	-out								outfile

	Usage : 
	[qualimap]				perl $0 -qualimap1 <genome_results.txt>  -qualimap2 <Qualimap report_ BAM QC.html> -fastqcdir  <sample fastqc dir> -tool < qualimap>  -conf <configure.txt> -out <outfile> 
	[verifybamid]			perl $0 -verifybamid  <LJD56.selfSM>   -tool <verifybamid>  -conf <configure.txt> -out <outfile> 

	[qualimap && verifybamid]  perl $0  -fastqcdir  <sample fastqc dir>  -qualimap1 <genome_results.txt>   -qualimap2 <Qualimap report_ BAM QC.html>  -verifybamid  <LJD56.selfSM>   -tool <both>  -conf <configure.txt> -out <outfile>
";

############################
if(!$software || !$outfile || !-f $configure)
{
	print "$usage\n";
	exit;
};

if( $software eq "both")
{
	if(  !-d $qcdir  ||  !-f $qualimapfile1 || !-f $qualimapfile2 || !-f $verifyfile )
	{
		print "$usage\n";
		exit;
	};
}elsif($software eq "qualimap")
{
	if( !-f $qualimapfile1 || !-f $qualimapfile2 ||!-d $qcdir )
	{
		print "$usage\n";
		exit;
	};
}elsif($software eq "verifybamid")
{
	if( !-f $verifyfile)
	{
		print "$usage\n";
		exit;
	};
}else
{
	print STDERR "please check your -tool !!!!\n";
	exit;
};

my %hash_config;  
&file_config($configure,\%hash_config,$software);

#print Dumper \%hash_config;

 my %hash_qualimap; my %hash_verifybamid; 

if($software eq "both" || $software eq "qualimap" )
{
	my @fastqcdirzip=glob "$qcdir/*fastqc.zip";
	
	my	$total_reads=&get_reads_num(\@fastqcdirzip);
	if($total_reads)
	{
		# Allow up to 1% read loss after alignment/deduplication (default)
		my $reads_ratio_cutoff = 0.99;
		my $cutoff_reads = int($total_reads * $reads_ratio_cutoff + 0.5);  # round to nearest
		$hash_config{'qualimap'}{'Total_reads_num_in_bam'}=$cutoff_reads;
		print STDERR "[INFO] Total reads from FASTQC: $total_reads, cutoff (>=99%): $cutoff_reads\n";
	}else
	{
		print STDERR "No Total Sequences in fastqc_data.txt !!\n ";
		die;
		exit(1);
	};
	&file_qualimap($qualimapfile1,$qualimapfile2,\%hash_qualimap)
};

if($software eq "both" || $software eq "verifybamid")
{

	&file_verifybamid($verifyfile,\%hash_verifybamid)
};
open O,">$outfile" or die "cannot write to file $outfile";

my $code=0;
if($software eq "both" || $software eq "qualimap")
{
	print O "[qualimap]\n";
	my %qualimap_result;
	&pre_qualimap(\%hash_config,\%hash_qualimap,\%qualimap_result);
	foreach my $final_result (sort keys %qualimap_result) 
	{
		
		if($qualimap_result{$final_result} =~/FAIL/)
		{
			$code=1;
			print STDERR "$final_result\t$qualimap_result{$final_result}\n";
		};
		print O "$final_result\t$qualimap_result{$final_result}\n";
	};
	print O "\n";	
};


if($software eq "both" || $software eq "verifybamid")
{
	print O "[verifybamid]\n";
	foreach my $index (sort keys %{$hash_config{'verifybamid'}}) 
	{
		my $result="PASS";
		if(exists $hash_verifybamid{$index})
		{
			if($hash_verifybamid{$index}>$hash_config{'verifybamid'}{$index})
			{
				$result="FAIL";
			};
			if($result eq "FAIL")
			{
				$code=1;
				print STDERR "$index\t$hash_config{'verifybamid'}{$index}\t$hash_verifybamid{$index}\t$result\n";
			};
			print O "$index\t$hash_config{'verifybamid'}{$index}\t$hash_verifybamid{$index}\t$result\n";
		}else
		{
			$result="NA";
			print O "$index\t$hash_config{'verifybamid'}{$index}\tNA\t$result\n";

		
		};
	};
};
close O;


if($code !=0)
{
	 print STDOUT "Unqualified \n";
	 exit(1);
}else
{
  print STDOUT "Successful!\n";
};




sub pre_qualimap
{
	my ($hash2_config,$hash2,$hash2_result)=@_;
	foreach my $key (sort keys %{$hash2_config->{'qualimap'}}) 
	{
		my $result="PASS";
		my $cutoff=$hash2_config->{'qualimap'}{$key};
		if(exists $hash2->{$key})
		{
			if($key eq "Total_reads_num_in_bam"  || $key eq "Coverage_depth_after_deduplication" || $key eq "Coverage_depth_on_1x" || $key eq "Coverage_depth_on_10x" || $key eq "Coverage_depth_on_20x" || $key eq "Mapping_rate" || $key eq "Average_read_length") 
			{
				$hash2->{$key}=~s/,//g;
				if($hash2->{$key}<$cutoff)
				{
					$result="FAIL";
				};
			}elsif($key eq "Gc_content")
			{
				my ($value1,$value2)=split/-/,$cutoff;
				if($hash2->{$key}<$value1 || $hash2->{$key}>$value2)
				{
					$result="FAIL";
				};
			}elsif( $key eq "Duplication_rate")
			{
				$hash2->{$key}=~s/,//g;
				if($hash2->{$key}>$cutoff)
				{
					$result="FAIL";
				};
			};
			$hash2_result->{$key}="$cutoff\t$hash2->{$key}\t$result";
		}else
		{
			$result="NA";
			$hash2_result->{$key}="$cutoff\tNA\t$result";
		};
	};
	
};
######($read1_fastqcfile,$read1_fastqcfile2,$read2_fastqcfile,$read2_fastqcfile2);


sub file_verifybamid
{
	my ($file,$hash)=@_;
	open F,"$file" or die "cannot open $file !!!";
	<F>;
	while(my $A=<F>)
	{
		chomp $A;
		next if($A=~/^$/);
		my @arr=split/\t/,$A;
		$hash->{'Contamination'}=$arr[6];
	};
	close F;
};



sub file_qualimap
{

	my ($file,$file2,$hash)=@_;
	open F,"$file" or die;
	while(my $A=<F>)
	{
		chomp $A;
		next if($A=~/^$/);
	
		if($A=~/number of reads = (.*)/)
		{
			$hash->{'Total_reads_num_in_bam'}=$1;
		};
		if($A=~/number of supplementary alignments = (.*)\(/)
		{
			my $sup=$1;
			$hash->{'Total_reads_num_in_bam'}=~s/,//g;
			$sup=~s/,//g;
			$hash->{'Total_reads_num_in_bam'}=$hash->{'Total_reads_num_in_bam'}-$sup;
		};
		if($A=~/mean coverageData = (.*)X/)
		{
			$hash->{'Coverage_depth_after_deduplication'}=$1;
		};
		if($A=~/There is a (.*)\% of reference with a coverageData >= 1X/)
		{
			$hash->{'Coverage_depth_on_1x'}=$1;
		};
		if($A=~/There is a (.*)\% of reference with a coverageData >= 10X/)
		{
			$hash->{'Coverage_depth_on_10x'}=$1;
		};
		if($A=~/There is a (.*)\% of reference with a coverageData >= 20X/)
		{
			$hash->{'Coverage_depth_on_20x'}=$1;
		};
		if($A=~/GC percentage = (.*)\%/)
		{
			$hash->{'Gc_content'}=$1;
		};
		if($A=~/number of mapped reads = .*?\((.*)\%\)/)   #########number of mapped reads = 126,158,440 (100%)
		{
			$hash->{'Mapping_rate'}=$1;
		};
		if($A=~/number of duplicated reads \(flagged\) = (.*)/)
		{
			my $dup_num=$1;
			$dup_num=~s/,//g;
			my $dup_rate=$dup_num/$hash->{'Total_reads_num_in_bam'};
			$dup_rate=sprintf("%.4f",$dup_rate);
			$hash->{'Duplication_rate'}=$dup_rate;
		};
	};
	close F;

	open F2,"$file2" or die;
	while(my $B=<F2>)
	{
		chomp $B;
		next if($B=~/^$/);
		if($B=~/Read.*\length/)
		{
			my $read_length=<F2>;
			if($read_length=~/.*\/.*\/(.*)\</)
			{
				my $read_mean=$1;
				$read_mean=~s/\s//g;
				$hash->{'Average_read_length'}=$read_mean;
			};
		};
	};
	close F2;
};


sub file_config                  #########&file_config($configure,\%hash_config,$tag,$software);
{
	my ($file,$hash,$qctool)=@_;
	my $key="NA";
	open M,"$file" or die;
	while (my $A=<M>)
	{
		chomp $A;
		next if ($A=~/^$/);
		my @arr=split/\t/,$A;
		if($A=~/\[(.*)\]/)
		{
			if($A!~/\[end\]/)
			{
				$key=$1;
			}else
			{
				next;
			};
			
		}else
		{
			if(($qctool eq $key )||($qctool eq "both"))
			{
				$hash->{$key}{$arr[0]}=$arr[1];
			};	
		};
	};
	close M;
};



sub get_reads_num 
{
	my ($filezip)=@_;
	my $read_num=0;
	for(my $i=0;$i<@$filezip;$i++)
	{
		my $zipfile = $filezip->[$i];
		my $zip=Archive::Zip->new($zipfile);
		my @members1=$zip->memberNames();
		foreach my $A (@members1)
		{
			if($A=~/fastqc\_data\.txt/)
			{
				my $neirongb= $zip->contents( $A ) . "\n";
				my @array_neirongb=split/\n/,$neirongb;
				foreach my $A_lineb (@array_neirongb)
				{
					if($A_lineb=~/^Total Sequences/)
					{
						 $read_num+=(split/\s/,$A_lineb)[-1];							
					};					
				};

			};
		};
	};
	return $read_num;
};
