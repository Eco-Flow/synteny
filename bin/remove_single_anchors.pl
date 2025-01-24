use strict;
use warnings;

# Check for correct number of arguments
if (@ARGV != 2) {
    die "Usage: perl script.pl <input_file> <output_file>\n";
}

# Get input and output file names from command-line arguments
my ($input_file, $output_file) = @ARGV;

# Open input and output files
open(my $in, '<', $input_file) or die "Cannot open $input_file: $!";
open(my $out, '>', $output_file) or die "Cannot open $output_file: $!";

my @lines;
while (<$in>) {
    chomp;
    push @lines, $_;
}

close $in;

# Process lines
my @filtered_lines;
my $inside_block = 0;
my $temp_block   = []; # Temporary storage for a block

foreach my $line (@lines) {
    if ($line eq '###') {
        if ($inside_block) {
            # End of a block
            if (@$temp_block == 1) {
                # Singleton block, discard it
                @$temp_block = ();
            } else {
                # Add valid block to the result
                push @filtered_lines, @$temp_block;
            }
            @$temp_block = (); # Reset temporary block
        }
        # Add the '###' to the output, ensuring no duplicates
        push @filtered_lines, '###' if !@filtered_lines || $filtered_lines[-1] ne '###';
        $inside_block = 1;
    } else {
        if ($inside_block) {
            # Collect lines within a block
            push @$temp_block, $line;
        } else {
            # Outside any block, add directly
            push @filtered_lines, $line;
        }
    }
}

# Write the processed lines to the output file
foreach my $line (@filtered_lines) {
    print $out "$line\n";
}

close $out;

print "Processing complete. Results saved in $output_file\n";

