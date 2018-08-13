#!/usr/bin/perl

# take uniprot-matches from stdin

use Getopt::Long;
use strict;

my $db_dir = "/net/bccs/cbu/berezovsky/babbit/db";

if ($#ARGV < 1) {
    my $b = `basename $0`;
    chomp $b;
    print "\nChecks which subgroups a profile matches based on a search on\n";
    print "either SCOP or uniref50.\n";
    print "\n\tusage: $b -n name-of-source -t type < matches\n\n";
    print "\ttype = 0, 1, 2 for scop, uniprot and pdb respectively\n";
    print "\tCode-mappings are taken from files in $db_dir/\n\n";
    exit;
}

my $type = -1;
my $code2sub_file = "";
my $input2code_file = "";
my $name = "";

GetOptions("t=i"=>\$type,
	   "n=s" => \$name);

if ($type == 0 ) { #scop
    $code2sub_file = "$db_dir/sfld.pdb.tab";
    $input2code_file = "$db_dir/dir.des.scop.txt_px_1.75.tab";
} elsif ($type == 1) { #uniprot
    $code2sub_file = "$db_dir/sfld.gi.tab";
    $input2code_file = "$db_dir/uniprot_gi_map.tab";
} elsif ($type == 2) { #pdb
    $code2sub_file = "$db_dir/sfld.pdb.tab";
} elsif ($type == -1) {
    die "Error: argument to -t must be specified\n";
} else {
    die "Error: type '$type' illegal\n";
}



my %input_matches;
my %input2subgroup; #gi/pdb to subgroup
my %code2subgroup; #uniprot/scop to subgroup
my %subgroupcount;


while (<>) {
    chomp;
    my @cols = split /\t/;
    my $input_id;
    if ($type == 0){
	$input_id = substr($cols[4], 0, 7);
    } elsif ($type == 1) {
	$cols[4] =~ m/^UniRef50_([\w-]+)/;
	$input_id = $1;
    } elsif ( $type == 2) {
	$input_id = lc substr($cols[4], 0, 4);
    }
    push @{$input_matches{$cols[0]}}, $input_id;
}

#creat hash from either gi or pdb codes to SFLD subgroups
open(INPUT, "< $code2sub_file") or 
    die "Error: failed to open file $code2sub_file\n";
while (<INPUT>) {
    chomp;
    my @cols = split /\t/;
    my $subgroup;
    my $code;
    if ($type == 0 || $type == 2) {
	$code = lc $cols[0];
	$subgroup = $cols[2];
	next if ($subgroup eq "?");
    } else {
	$subgroup = $cols[1];
	$code = $cols[5];
	next if ($subgroup eq "");
    } 
    $subgroup =~ s/\s/_/g;
    $code2subgroup{$code} = $subgroup;
}


# this part is actually superfluous for SCOP and PDB entries but
# in the future we might want to use information about chain and
# residues of the domains, so I'm keeping it this way for now...

close INPUT;
if ($type == 0 || $type == 1) {
    
    open(INPUT, "< $input2code_file") or 
	die "Error failed to open file $input2code_file\n";
    
    while (<INPUT>) {
	chomp;
	my @cols = split /\t/;
	my $code;
	my $input;
	if ($type == 0) {
	    $input = $cols[3];
	    $code = lc substr($cols[4], 0, 4);
	} elsif ($type == 1) {
	    $code = $cols[2];
	    $input = $cols[0];
	} 
	if (exists $code2subgroup{$code}) {
	    push @{$input2subgroup{$input}}, $code2subgroup{$code};
	}
   }
    
   close INPUT;
} else {
    # stupid, yes
    foreach my $code (keys %code2subgroup) {
	if (exists $code2subgroup{$code}) {
	    push @{$input2subgroup{$code}}, $code2subgroup{$code};
	}
    }
}

my %summary;
my %network;

my $sif_file;
my $prop_file;
if ($type == 0) {
    $sif_file = "scop_$name.sif";
    $prop_file = "scop_node_prop.tab"
} elsif ($type == 1) {
    $sif_file = "uniref50_$name.sif";
    $prop_file = "uniref50_node_prop.tab";
} elsif ($type == 2) {
    $sif_file = "pdb_$name.sif";
    $prop_file = "pdb_node_prop.tab";
}

open (SIF, ">$sif_file") or
    die "Error opening file sif_file\n";
open (NODE_PROP, ">$prop_file") or
    die "Error opening file $prop_file\n";
my @profiles = sort {$a <=> $b} keys %input_matches;
foreach my $profile (@profiles) {
    my @prof_matches = @{$input_matches{$profile}};
    my $N = $#prof_matches;
    for (my $i = 0; $i <= $N; ++$i) {
	my $matchi = $prof_matches[$i];
	my $group;
	if (exists $input2subgroup{$matchi}) {
	    my @subgroup = @{$input2subgroup{$matchi}};
	    $group = $subgroup[0];
	} else {
	    $group = "??";
	}

	print "# $profile\t$matchi\t$group\n";
	$summary{$profile}{$group}{$matchi} = 1;
	$subgroupcount{$group}{$matchi} = 1;
	print NODE_PROP "$matchi\t$group\n";

	for (my $j = $i+1; $j < $N; ++$j) {
	    print SIF "$matchi\t$profile\t$prof_matches[$j]\n";
	}	
    }
    
}
close SIF;
close NODE_PROP;

my %prof_prof_network;
foreach my $l1 (@profiles) {
    my @l1_matches = @{$input_matches{$l1}};
    foreach my $l2 (@profiles) {
	next if ($l1 ge $l2);
	my @l2_matches = @{$input_matches{$l2}};
	foreach my $p1 (@l1_matches) {
	    foreach my $p2 (@l2_matches) {
		if ($p1 eq $p2) {
		    my $G = $input2subgroup{$p1}->[0];
		    ++$prof_prof_network{$G}{$l1}{$l2};
		}
	    }
	}
    }
}

open (TAB, ">prof_prof.tab") or die "Error writing prof_prof.tab\n";

while (my ($group, $valuei) = each(%prof_prof_network)) {
    while (my ($l1, $valuej) = each(%{$valuei})) {
	while (my ($l2, $count) = each (%{$valuej})) {
	    $group =~ s/ /_/g;
	    $group =~ s/^$/\?\?/;
	    print TAB "$l1\t$group\t$l2\t$count\n";
	}
    }
}
close TAB;

my @subgroups = sort keys %subgroupcount;
my @prof_index = sort {$a <=> $b} keys %summary;

print "subgroup\tN";
foreach (@prof_index) {
    print "\t$_";
}
print "\n";
foreach my $group (@subgroups) {
    my $N = scalar keys %{$subgroupcount{$group}};
    print "$group\t$N";
    foreach my $prof (@prof_index) {
	printf("\t%1.3f",((scalar keys %{$summary{$prof}{$group}})/$N));
    }
    print "\n";
}

#print "\n# Database contents:\n";
#foreach (subgroups) {
#    print "$_ ";
#    print scalar keys %{$subgroupcount{$_}};
#    print "\n";
#}
#
#print "\n# Summary:\n";
#foreach my $prof_index (sort {$a <=> $b} keys %summary) {
#    my @subgroup = sort keys %{$summary{$prof_index}};
#    print "# $prof_index: "; 
#    foreach (@subgroup) {
#	#next if ($subgroupcount{$_} <= 0);
#	print "$_ ";
#	print ((scalar keys %{$summary{$prof_index}{$_}})/(scalar keys %{$subgroupcount{$_}}));
#	print ", ";
#    }
#    print "\n";
#}






