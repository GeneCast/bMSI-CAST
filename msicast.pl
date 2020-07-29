#! /urs/bin/perl -w
use strict;
use List::Util;
use Getopt::Long;
use Statistics::Basic qw(:all);

my ($bin, $bam, $flank, $num, $minsize, $quality, $mode, $outdir, $sample, $help);
GetOptions(
                        "bin:s" => \$bin,
                        "bam:s" => \$bam,
			"flank:i"=> \$flank,
                        "num:i" => \$num,
			"minsize:i" => \$minsize,
			"quality:i"  => \$quality,
			"outdir:s" => \$outdir,
			"sample:s" => \$sample,
			"mode:s" => \$mode,
                        "h|?" => \$help,
);

if(!$bin || !$bam || $help){
        die &USAGE;
}

$flank ||= 2;
$num||= 2;
$minsize||=1;
$quality||=1;
$mode||="DelRatio";
my $e=0.0001;

my $baseline;
my $bed;
my $distribution;
if ($mode eq "DelRatio"){
	$baseline="$bin/delratio_baseline.txt";
	$bed="$bin/delratio_bed.txt";
}else{
	$baseline="$bin/kld_baseline.txt";
	$bed="$bin/kld_bed.txt";
	$distribution="$bin/distribution.txt";
}

my %baseline;
my %duprange;
open BASELINE,"$baseline"  or die "[ERROR] file $baseline open failed!\n";
while(<BASELINE>){
        chomp;
        next if (/^Name/);
        my @tmp=split /\t/;
	my $marker=shift @tmp;
	my $min=(split /\|/,$tmp[0])[1];
	my $max=(split /\|/,$tmp[-1])[1];
	my $info=join "\t",@tmp;
        $baseline{$marker}=$info;
	$duprange{$marker}="$min\t$max";
}
close BASELINE;

my %hash;
my %chr;
my %pos;
my %markerlen;
my %dp;
open BED,"$bed" or die "[ERROR] file $bed open failed!\n";
open BEDT,">$outdir/$sample\.bed.tmp" or die "[ERROR] file $outdir/$sample\.bed.tmp open failed!\n";
while(<BED>){
	chomp;
	next if (/Start/);
	my @tmp=split /\t/;
	$hash{$tmp[3]}=$_;
	$chr{$tmp[0]}=1;
	my $pos="$tmp[0]:$tmp[1]-$tmp[2]";
	$pos{$tmp[3]}=$pos;
	my $extends=$tmp[1];
	my $extende=$tmp[2];
	print BEDT "$tmp[0]\t$extends\t$extende\n";
	$markerlen{$tmp[3]}=$tmp[2]-$tmp[1];
	$dp{$tmp[3]}=$tmp[4];
}
close BED;
close BEDT;

my %dis4allele;
my %maxdisdup;
if ($mode ne "DelRatio" ){
	open DISTRIBUTION,"$distribution"  or die "[ERROR] file $distribution open failed!\n";
        while(<DISTRIBUTION>){
                chomp;
                next if (/^Name/ || /^$/ || /^#/);
                my @tmp=split /\t/;
                my $name=shift @tmp;
                my $maxdup=(split /\|/,$tmp[-1])[-1];

                $maxdisdup{$name}=$maxdup;
                foreach my $info(@tmp){
                        my @alleles=split /\|/,$info;
                        my $dup=pop @alleles;
                        my $max=0;
                        my $maxcut=0;
                        my $maxtype="";
                        foreach my $alleles(@alleles){
                                my($type,$ratio,$sd)=split /\:/,$alleles;
                                if ($ratio>$max){
                                        $maxtype=$type;
                                        $max=$ratio;
                                        $maxcut=$ratio+3*$sd;
                                }
                        }
                        $info=join "|",@alleles;
                        $dis4allele{$name}{$dup}="$maxtype\t$maxcut\t$info";
                }
        }
	close DISTRIBUTION;
}

my %ref2reads;
my %read;
my %insert;
my %dup;

open SAM,"samtools view -F 0xB00 -L $outdir/$sample\.bed.tmp -q $quality $bam|" or die "[ERROR] file $bam open failed!\n";### remove multi alignment
open OUT1,">$outdir/$sample\.result.txt" or die "[ERROR] file $outdir/$sample\.result.txt open failed!\n";
print OUT1 "ReadID\tMsequence\tReadSequence\n";
open OUT2,">$outdir/$sample\.stat.txt"or die "[ERROR] file $outdir/$sample\.stat.txt open failed!\n";
print OUT2 "Position\tName\tAverage_Depth\t$mode\tDup_Ratio\tIndelLength:AlleleFraction:SupportingCalls\n";
open OUT3,">$outdir/$sample\.msi_base.txt" or die "[ERROR] file $outdir/$sample\.msi_base.txt open failed!\n";
my $head;
foreach my $gene(sort keys %hash){
	$head.="\t$gene";
}
print OUT3 "sample_id\tunstable_loci\tpassing_loci\tmsi_score\tmsi_status$head\n";

my %all;
while(<SAM>){
	chomp;
	my @tmp=split /\t/;
	my $info="$tmp[0]\t$tmp[9]";
	my $StartAndDirection;
        my $flagall=$tmp[1];
	if(exists $chr{$tmp[2]} && !($flagall & 4)){
		%insert=();
		my $cigar=$tmp[5];$cigar=~s/S/S;/g;$cigar=~s/M/M;/g;$cigar=~s/D/D;/g;$cigar=~s/I/I;/g;$cigar=~s/H/H;/g;
                my @cigar=split /;/,$cigar;
		
		my $lenforend=0;
		my $endflag=0;
		foreach my $c(@cigar){
			if($c=~/(\d+)M/||$c=~/(\d+)D/){
				$lenforend+=$1;
			}
		}
		
		foreach my $gene(sort keys %hash){
                        my @str=split /\t/,$hash{$gene};
                        if ($tmp[2] eq $str[0] && ($str[1]+1-$tmp[3]>=$flank) && (($tmp[3]+$lenforend-1)-$str[2]>=$flank)) {
                        	$endflag=1        
			}
		}
		next if ($endflag==0);

                if ($tmp[6] eq "=" && $tmp[8]!=0){
                        if ($tmp[8]>0 ){
                                $StartAndDirection="$tmp[2]:$tmp[3]:$tmp[8]:+";
                        }elsif($tmp[8]<0){
                                $tmp[8]=0-$tmp[8];
                                $StartAndDirection="$tmp[2]:$tmp[7]:$tmp[8]:-";
                        }
                }else{  
                        my $cigar=$tmp[5];$cigar=~s/S/S;/g;$cigar=~s/M/M;/g;$cigar=~s/D/D;/g;$cigar=~s/I/I;/g;$cigar=~s/H/H;/g;
                        my @cigar=split /;/,$cigar;
                        my $lenforend=0;
                        foreach my $c(@cigar){
                                if($c=~/(\d+)M/||$c=~/(\d+)D/){
                                        $lenforend+=$1;
                                }
                        }
                	$StartAndDirection="$tmp[2]:$tmp[3]:$lenforend:$tmp[6]:0";
                }
		my $readsstart=0; my $lengthref=0;my $lengthreads=0;
                if ($cigar[0]=~/(\d+)S/){
                	$readsstart=$1;
		}elsif($cigar[0]=~/(\d+)M/){
			$lengthref=$1;
			$lengthreads=$1;
			for(my $j=1;$j<=$lengthref;$j++){
	 			$ref2reads{$j}=$j;
 			}
		}

                for (my $i=1;$i<@cigar;$i++){
               	        if ($cigar[$i]=~/(\d+)M/){
                       	        my $dis=$lengthreads-$lengthref;
				$lengthref+=$1;
				$lengthreads+=$1;
				my $lengthreftmp=$1;
				my $start=$lengthref-$lengthreftmp+1;
				for(my $j=$start;$j<=$lengthref;$j++){
					$ref2reads{$j}=$readsstart+$j+$dis;
				}
                        }
                        if ($cigar[$i]=~/(\d+)I/){
                        	my $dis=$lengthreads-$lengthref;
				$lengthreads+=$1;
				my $lengthreadstmp=$1;
				my $start=$lengthreads-$lengthreadstmp+1;
				for(my $j=$start;$j<=$lengthreads;$j++){
                                	$ref2reads{$lengthref}=$readsstart+$j+$dis;
					my $pos=$lengthref+$tmp[3]-1;
					$insert{$pos}=1
                                }
			}
			if ($cigar[$i]=~/(\d+)D/){
				my $dis=$lengthreads-$lengthref;
				$lengthref+=$1;
				my $lengthreftmp=$1;
				my $start=$lengthref-$lengthreftmp+1;
				for(my $j=$start;$j<=$lengthref;$j++){
                                	$ref2reads{$j}=$readsstart+$lengthreads+$dis;
					my $pos=$lengthref+$tmp[3]-1;
                                }
			}
                }
		foreach my $gene(sort keys %hash){
			next if (exists $read{$gene}{$tmp[0]}); 
			my @str=split /\t/,$hash{$gene};
			if ($tmp[2] eq $str[0]){
				if (($str[1]+1-$tmp[3]>=$flank) && (($tmp[3]+$lengthref-1)-$str[2]>=$flank)){
					if (!(%insert && !exists $insert{$str[1]})){
						my $sp=($str[1]-1)-$tmp[3]+1;
						my $ep=$str[2]-$tmp[3]+1; 
						my $rl=$ref2reads{$ep}-$ref2reads{$sp}-1;
						if ($str[1]+1-$tmp[3]==$flank){
							my $pre=substr($tmp[9],0,$flank+1);
							my @pre=split //,$pre;
							my $base=$pre[-1];
							my $flag=0;
							for (my $i=0;$i<@pre-1;$i++){
								if ($pre[$i] ne $base ){
									$flag=1;
								}
							}
							next if ($flag==0);
						}
						if (($tmp[3]+$lengthref-1)-$str[2]==$flank){
							my $post=substr($tmp[9],-($flank+1));
							my @post=split //,$post;
							my $base=$post[0];
							my $flag=0;
							for (my $i=1;$i<@post;$i++){
								if ($post[$i] ne $base ){
									$flag=1;
								}
							}
							next if ($flag==0);
						}
						my $substart=$ref2reads{$sp}+1; 
						my $subreads=substr($tmp[9],$substart,$rl);
						my $Len=length($subreads);
						$all{$gene}++;
						print OUT1 "$tmp[0]\t$subreads\t$tmp[9]\t$tmp[5]\t$gene\n";	
						$read{$gene}{$tmp[0]}=1;
						
						my @dupinfo=split /\:/,$StartAndDirection;
                                                my $position;
                                                if ($dupinfo[-1] eq "0"){
                                                        $position="$gene\t$dupinfo[0]\t$dupinfo[1]\t$dupinfo[2]:$dupinfo[3]:s";
                                                }else{
                                                        $position="$gene\t$dupinfo[0]\t$dupinfo[1]\t$dupinfo[2]";
                                                }
						$dup{$position}.="$Len ";		
					}
				}
			}	
		}
	}
}
close SAM;
close OUT1;

unlink "$outdir/$sample\.bed.tmp";
my %mode;
my %mult_mode;

open DUP, "> $outdir/$sample\.DUP.txt" or die "[ERROR] file $outdir/$sample\.DUP.txt open failed!\n";

foreach my $p( sort keys %dup){
	my @tmp=split /\t/,$p;
	my $gene=$tmp[0];
	my $pos="$tmp[1]\t$tmp[2]\t$tmp[3]";
	$dup{$p}=~s/\s+$//;
	my @len=split /\s+/,$dup{$p};
	next if (@len <$minsize );
	my $num4mode=mode(@len);
	$num4mode=~s/\[//;$num4mode=~s/\]$//;$num4mode=~s/\,//g;
	my @num4mode=split /\s+/,$num4mode;
	if (@num4mode ==1){
		$mode{$gene}{$num4mode}++;
		print DUP "$p\t$dup{$p}\n";
	}else{
		$mult_mode{$gene}{$pos}=$num4mode;
	}
}

close DUP;

my $positive=0;
my $passloci=0;
my %marker;

foreach my $gene(sort keys %mode){
	my $delratio=0;
	my $kld=0;
	my $value=0;
	my $peakallele;
	my $peak=0;
	my $hitnum=0;
	my $depth=0;
	my $line="";
	my $dup_base=0;
	my $dup_distr=0;	
	foreach my $length(sort keys %{$mode{$gene}}){
		if ($mode{$gene}{$length} >= $num ){
			$depth=$depth+$mode{$gene}{$length};
			if ($mode{$gene}{$length} >= $hitnum){
				$hitnum=$mode{$gene}{$length};
			}
		}
	}
	my $dup_ratio=sprintf("%.2f",($all{$gene}-$depth)/$all{$gene});
	$dup_base=$dup_ratio;
	my ($min,$max)=split /\t/,$duprange{$gene};
	
	if($dup_ratio>$max){
		$dup_base=sprintf("%.2f",$max);
	}

	if ($hitnum>0){
		foreach my $length(sort {$a<=>$b} keys %{$mode{$gene}}){
                	my $ratio=$mode{$gene}{$length}/$depth;
			my $len=$length-$markerlen{$gene};
			if ($ratio>$peak){
                                $peak=$ratio;
                                $peakallele=$len;
                        }
	
			if ($mode{$gene}{$length}>=$num){
		 		$line.="$len:$ratio:$mode{$gene}{$length} ";
		 		if ($mode eq "DelRatio" && $len <0){
					$delratio=$delratio+$ratio;
		 		}
        		}
		}
		$line=~s/\s+$//;
	
		if ($mode ne "DelRatio"){
			next if (!exists $dis4allele{$gene});
	                my %allele=();
        	        my %q=();
                	$dup_distr=$dup_ratio;
                	if ($dup_distr>$maxdisdup{$gene}){
                        	$dup_distr=$maxdisdup{$gene};
                	}
                	next if (!exists $dis4allele{$gene}{$dup_distr});
                	my @dis4allele=split /\t/,$dis4allele{$gene}{$dup_distr};
                	my @q=(split /\|/,$dis4allele[2]);
                	foreach my $q(@q){
                        	my $allele=(split /\:/,$q)[0];
                        	$q{$allele}=(split /\:/,$q)[1];
                        	$allele{$allele}=1;
                	}

               		 my %p=();
                	 my @p=split /\s+/,$line;
                	 my $depth=0;
                	 foreach my $p(@p){
                        	my $d=(split /\:/,$p)[2];
                        	$depth+= $d;
                	}
                	foreach my $p(@p){
                        	my $pro=sprintf("%0.4f",((split /\:/,$p)[2])/$depth);
                        	my $allele=(split /\:/,$p)[0];
                        	if ($pro>0){
                                	$p{$allele}=$pro;
                                	$allele{$allele}=1;
                        	}
                	}
                	my $all=keys %allele;
                	my $qnum=@q;
                	my $qnon=$all-$qnum;
                	my $pnum=keys %p;
                	my $pnon=$all-$pnum;

                	foreach my $allele(sort {$a<=>$b} keys %allele){
                        	if ($qnon>0){
                                	if (exists $q{$allele}){
						$q{$allele}=sprintf("%0.5f",$q{$allele}-$e/$qnum);
        	                        }elsif (!exists $q{$allele}){
                	                        $q{$allele}=sprintf("%0.5f",$e/$qnon);
                        	        }
                        	}

                        	if ($pnon>0){
                                	if (exists $p{$allele}){
                                        	 $p{$allele}=sprintf("%0.5f",$p{$allele}-$e/$pnum);
                                	}elsif (!exists $p{$allele}){
                                        	 $p{$allele}=sprintf("%0.5f",$e/$pnon);
                                	}
                        	}
                	}

                	my $maxq=0;
                	my $maxp=0;
                	my $maxqalle;
                	my $maxpalle;
                	my $qtmp=0;
                	my $ptmp=0;

                	foreach my $allele(sort {$a<=>$b} keys %allele){
                        	$qtmp+=$q{$allele};
                        	if ($q{$allele}>$maxq){
                                	$maxq=$q{$allele};
                                	$maxqalle=$allele;
                        	}
                        	$ptmp+=$p{$allele};
                        	if ($p{$allele}>$maxp){
                                	$maxp=$p{$allele};
                                	$maxpalle=$allele;
                        	}
                	}
                	my $splitq=($qtmp-1);
                	my $splitp=($ptmp-1);
                	$q{$maxqalle}=$q{$maxqalle}-$splitq;
                	$p{$maxpalle}=$p{$maxpalle}-$splitp;

                	foreach my $allele(sort {$a<=>$b} keys %allele){
                        	$kld+=$p{$allele}*((log($p{$allele})/log(2))-(log($q{$allele})/log(2)));
                	}
		}
		if ($mode eq "DelRatio"){
			$value=$delratio;
		}else{
			$value=$kld;
		}
		if ($depth >=$dp{$gene}){
                        $passloci++;
			my @baseline=split /\t/,$baseline{$gene};
			foreach my $info(@baseline){
				my ($cutoff,$dupbase)=split /\|/,$info;
				if ($dup_base == $dupbase){
					if ($value>=$cutoff){
						$marker{$gene}=1;
						$positive++;
						if ($mode ne "DelRatio"){
							my ($maxtype,$maxcut)=(split /\t/,$dis4allele{$gene}{$dup_distr})[0,1];
							if($peak >= $maxcut && $peakallele eq $maxtype){
								$marker{$gene}=0;
								$positive--;
							}
						}
					}else{
						$marker{$gene}=0;
					}
					last;
				}
			}                      
                }
                
		print OUT2 "$pos{$gene}\t$gene\t$depth\t$value\t$dup_ratio\t$line\n";
	}
}
close OUT2;

my $status;
my $pratio=sprintf("%.3f",$positive/$passloci);

if ($passloci>=15){
        if ($pratio >=0.1){
                $status="MSI-H";
        }else{$status="MSS";}
}else{
        $status="QNS";
}
print OUT3 "$sample\t$positive\t$passloci\t$pratio\t$status";

foreach my $gene(sort keys %hash){
        if (!exists $marker{$gene}){
                print OUT3 "\t";
        }else{
                print OUT3 "\t$marker{$gene}"
        }
}
print OUT3 "\n";
close OUT3;

sub USAGE{
        my $usage =<<"USAGE"

usage: perl $0 [options]
options:
      -bin <s> path of database
      -bam <s> alignment bam file
      -flank <i> flank size, default [2]
      -num <i> consistency numbern,default[2]
      -minsize <i> min read number of a family,default [1]
      -quality <i> mapping quality,default[1]
      -mode <s>	calling algorithm, default [DelRatio]
      -s <s> prefix of result files
      -o <s> output
example: perl $0 -bin ./bin  -bam result/sample20171207-B_I25.sorted.mkdup.realign.bam -flank 2 -num  2 -minsize 1 -quality 1 -mode DelRatio -o MSI -s sample20171207-B_I25
author : zhao.lili\@genecast.com.cn

USAGE
}
