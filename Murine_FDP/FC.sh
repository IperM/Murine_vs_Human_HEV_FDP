#!/bin/bash
: <<COMMENT

# Directory containing your marked_dup_* files
input_directory="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Picard/Human"
output_directory="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/FeatureCounts/Human_BWA"

# Loop through each file matching the pattern marked_dup_*
for input_file in "$input_directory"/marked_dup_*; do
	# Check if there are matching files
	if [[ ! -e "$input_file" ]]; then
		echo "No files matching the pattern marked_dup_* found in $input_directory"
		break
	fi
					        
	# Get the base name of the input file (without directory)
	base_name=$(basename "$input_file")

	# Generate the output file name by adding 'count' at the start and changing the extension to .csv
	output_file="$output_directory/count${base_name%.*}.csv"

	# Run the featureCounts command with the generated output file name
	featureCounts -T 8 -F GTF -t exon -a /mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.111.gtf -g gene_id -o "$output_file" "$input_file"
done

COMMENT

input_directory="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Picard/Murine"
output_directory="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/FeatureCounts/Murine"_BWA

# Loop through each file matching the pattern marked_dup_*
for input_file in "$input_directory"/marked_dup_*; do
	# Check if there are matching files
	if [[ ! -e "$input_file" ]]; then
		echo "No files matching the pattern marked_dup_* found in $input_directory"
		break
	fi

	# Get the base name of the input file (without directory)
	base_name=$(basename "$input_file")

	# Generate the output file name by adding 'count' at the start and changing the extension to .csv
	output_file="$output_directory/count${base_name%.*}.csv"

	# Run the featureCounts command with the generated output file name
	featureCounts -T 8 -F GTF -t exon -a /mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.111.gtf -g gene_id -o "$output_file" "$input_file"
done
