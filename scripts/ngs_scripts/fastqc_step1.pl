# !/usr/bin/perl
# Source: from the upstream variant-analysis pipeline scripts used for the
# China 100K Genomes Project cohort (internal deployment).
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );


my ($fastqczip1,$fastqczip2,$outdir,$configure);

GetOptions
(
		"fastq1=s"		=> \$fastqczip1,
		"fastq2=s"		=> \$fastqczip2,
		"conf=s"	=> \$configure,
		"dirout=s"		=>	\$outdir
);


my $usage="
	perl $0 -fastq1 <read1_fastqc_result.zip>  -fastq2 <read2_fastqc_result.zip>  -conf <configure.v1> -dirout <outdir>  

	-fastq1								the path of read1  fastqc result file: *.zip
	-fastq2								the path of read2  fastqc result file: *.zip
	-conf								configure.txt
	-dirout								outfile

	Usage : 
	[fastqc]				perl $0 -fastq1 <read1_fastqc_result.zip> -fastq2 <read2_fastqc_result.zip>     -conf <configure.txt>  -dirout <outdir> 

";

if(!$fastqczip1|| !$fastqczip2 || !-f $configure || !$outdir)
{
	print STDERR  "$usage\n";
	exit;
};
if(!-d $outdir ){`mkdir $outdir `;};
$fastqczip1=`readlink -f $fastqczip1 `; chomp $fastqczip1;
$fastqczip2=`readlink -f $fastqczip2 `; chomp $fastqczip2;
$configure=`readlink -f $configure `; chomp $configure;
$outdir=`readlink -f $outdir `; chomp $outdir;


my $software="fastqc";
my %hash_config;  
&file_config($configure,\%hash_config,$software);
my %hash_fastqc;
&file_fastqc($fastqczip1,$fastqczip2, \%hash_fastqc );

my $readset=`readlink -f $fastqczip1 `; chomp $readset;
my $name=(split/\//,$readset)[-2];

open O,">$outdir/$name\.txt" or die "cannot write to file $outdir/$name";

foreach my $index (sort keys %{$hash_config{'fastqc'}}) 
{
	my $result="PASS";
	my ($result_tag1,$result_tag2)=("NA","NA");
	if(exists $hash_fastqc{$index})
	{
		my @index=@{$hash_fastqc{$index}};
		($result_tag1,$result_tag2)=($index[0],$index[1]);
		print O "$index\t$result_tag1|$result_tag2\n";

	}else
	{
		print O "$index\tNA|NA\n";

	};
};

print O "\n";	
close O;


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



sub file_fastqc
{
	my ($filezip1,$filezip2,$hash)=@_;

	my @filezip=($filezip1,$filezip2);
	
	for(my $i=0;$i<@filezip;$i++)
	{
		my $zipfile = $filezip[$i];
		my $zip=Archive::Zip->new($zipfile);
		my @members1=$zip->memberNames();
		foreach my $A (@members1)
		{
			if($A=~/summary\.txt/)
			{
				my $neirong= $zip->contents( $A ) . "\n";
				my @array_neirong=split/\n/,$neirong;
				foreach my $A_line (@array_neirong)
				{
					my @A_line_array=split/\t/,$A_line;
					push @{$hash->{$A_line_array[1]}}, $A_line_array[0];

				};
			};

			if($A=~/fastqc\_data\.txt/)
			{
				my %temp;
				my $neirongb= $zip->contents( $A ) . "\n";
				my @array_neirongb=split/\n/,$neirongb;
				my $key="NO";
				foreach my $A_lineb (@array_neirongb)
				{
					if($A_lineb=~/^\#Quality/)
					{
						$key="YES";
						next;
					};
					if($key eq "YES")
					{
						if($A_lineb!~/\>\>END_MODULE/)
						{
							my @A_lineb_array=split/\t/,$A_lineb;
							$temp{$A_lineb_array[0]}=$A_lineb_array[1];
						}else
						{
							last;
						};
					}else
					{
						next;
					};
				};
				my $all_1=0; my $Q20_1=0; my $Q30_1=0;
				foreach  (keys %temp) 
				{
					$all_1+=$temp{$_};
					if($_ <=20 )
					{
						$Q20_1+=$temp{$_};	
					};
					if($_ <=30)
					{
						$Q30_1+=$temp{$_};	
					};

				};
				$Q20_1=$all_1-$Q20_1;
				$Q30_1=$all_1-$Q30_1;
				#my $rounded_Q20_1 = sprintf("%.4f", $Q20_1);
				#my $rounded_Q30_1= sprintf("%.4f", $Q30_1);
				push @{$hash->{Qualification_rate_of_Q20}}, "$Q20_1;$all_1";
				push @{$hash->{Qualification_rate_of_Q30}}, "$Q30_1;$all_1";
			};
		};
	};
};
