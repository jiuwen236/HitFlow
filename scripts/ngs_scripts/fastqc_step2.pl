#!/usr/bin/perl
# Source: from the upstream variant-analysis pipeline scripts used for the
# China 100K Genomes Project cohort (internal deployment).
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename; 

my ($dir,$configure,$out,$names);

GetOptions
(
		"indir=s"		=> \$dir,
		"conf=s"	=> \$configure,
		"outfile=s"		=>	\$out,
		"filename=s" => \$names
);


my $usage="
	perl $0 -indir <the fastqc readset result dir>   -conf <configure.v1> -outfile <the  result>  
	-indir								the path of fastqc readset dir
	-conf								configure.txt
	-outfile							result out
	-filename							filename
	Usage : 
	[fastqc]				perl $0 -indir <the fastqc readset result dir>   -conf <configure.v1> -outfile <the  result>  -filename <names>

";
if(!$dir || !-f $configure || !$out)
{
	print STDERR "$usage\n";
	exit;
};
$configure=`readlink -f $configure `;chomp $configure;
my $dir_1=`readlink -f $dir `;chomp $dir_1;
my $software="fastqc";

my @resultfile=glob "$dir_1/readset*txt";

my %hash_config;  
&file_config($configure,\%hash_config,$software);
#print Dumper \%hash_config;
my %hash_fastqc;
for(my $i=0;$i<@resultfile;$i++)
{
	open F,"$resultfile[$i]" or die;
	while(my $A=<F>)
	{
		chomp $A;
		next if($A=~/^$/);
		my @arr=split/\t/,$A;
		if($arr[0] eq "Qualification_rate_of_Q20" || $arr[0] eq "Qualification_rate_of_Q30")
		{
			my ($Q1,$A1,$Q2,$A2)=split/[;|]/,$arr[-1];
			push @{$hash_fastqc{$arr[0]}},($Q1,$Q2);
			push @{$hash_fastqc{"All\_$arr[0]"}},($A1,$A2);
		}else
		{
			my ($tag1,$tag2)=split/\|/,$arr[-1];
			push @{$hash_fastqc{$arr[0]}},($tag1,$tag2);
		};
	};
	close F;
};

my $exit_code=0;
open O,">$out" or die;
open O1,">$out.json" or die;
print O1 "\{\"$names\":\{";
my $count=1;
foreach my $index_name (sort keys %{$hash_config{'fastqc'}}) 
{
	my $result="PASS";
	if(exists $hash_fastqc{$index_name})
	{
		my @index1=@{$hash_fastqc{$index_name}};
		if($index_name eq "Sequence Length Distribution")
		{	
			foreach my $tag (@index1) 
			{
				if($tag ne "PASS")
				{
					$result="FAIL";
				};
			};	
		}elsif($index_name eq "Qualification_rate_of_Q20" || $index_name  eq "Qualification_rate_of_Q30" )
		{
			my ($ratio,$part,$all)=(0,0,0);
			foreach my $tag (@index1) 
			{
				$part+=$tag;
			};	
			foreach my $tag (@{$hash_fastqc{"All\_$index_name"}}) 
			{
				$all+=$tag;
			};	
			
			$ratio=$part/$all;
			if($ratio < $hash_config{'fastqc'}{$index_name} )
			{
				$result="FAIL";
			};
		
		}else
		{
			foreach my $tag (@index1) 
			{
				if($tag eq $hash_config{'fastqc'}{$index_name})
				{
					$result="FAIL";
				};
			};	
		};	
		if($result eq "FAIL")
		{
			print STDERR "FAIL: $index_name\t$hash_config{'fastqc'}{$index_name}\t$result\n";
			 $exit_code=1;
		};
		print O "$index_name\t$hash_config{'fastqc'}{$index_name}\t$result\n";
		if($count==1){print O1 "\"$index_name\":\"$result\"";}
		else{print O1 ",\"$index_name\":\"$result\"";}
		$count++;
	}else
	{
		$result="NA";
		print O "$index_name\t$hash_config{'fastqc'}{$index_name}\t$result\n";
		if($count==1){print O1 "\"$index_name\":\"$result\"";}
		else{print O1 ",\"$index_name\":\"$result\"";}
		$count++;
	};
};
print O1 "\}\}";
close O;close O1;

if($exit_code!=0)
{
  #$exit_code=system("error error"); 
  print STDOUT "Unqualified \n";
  #exit($exit_code >> 8);
  exit(1);
}else
{
  print STDOUT "Successful!\n";
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


