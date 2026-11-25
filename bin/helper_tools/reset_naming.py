import os
import argparse
import glob
import re


def replace_underscores(value, replacements):
    """Replace selected underscores by their position in the original value."""
    position = 0
    result = []
    for character in value:
        if character == '_':
            position += 1
            character = replacements.get(position, character)
        result.append(character)
    return ''.join(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "This script prepares Bruker raw data before running "
            "1_PV2NIfTiConverter/pv_conv2Nifti.py. "
            "The raw data must follow the structure: projectfolder/subjects/ses/data. "
            "It automatically scans for all 'subject' files within the input folder and "
            "performs two modifications: "
            "1) Replaces selected underscores in SUBJECT_id and SUBJECT_study_name using "
            "--r1 TEXT through --r6 TEXT (use an empty string to delete); "
            "positions refer to each original value and unselected underscores stay unchanged. "
            "2) Replaces the word 'baseline' (case-insensitive) with 'PT0' in the study name. "
            "A corrected version of each 'subject' file is written back to disk. "
            "Example usage: "
            "python reset_naming.py -i /Volumes/Desktop/MRI/raw_data --r1 c"
        )
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the parent project folder containing the dataset, e.g. raw_data',
        type=str
    )

    for position in range(1, 7):
        parser.add_argument(
            f'--r{position}', metavar='TEXT',
            help=f'Replace original underscore {position} with TEXT; use "" to delete it.'
        )

    args = parser.parse_args()
    replacements = {}
    for position in range(1, 7):
        replacement = getattr(args, f'r{position}')
        if replacement is not None:
            replacements[position] = replacement

    # Get list of raw data folders or files in input directory
    list_of_raw = sorted([
        d for d in os.listdir(args.input)
        if os.path.isdir(os.path.join(args.input, d))
        or (os.path.isfile(os.path.join(args.input, d)) and (('zip' in d) or ('PvDataset' in d)))
    ])

    # Recursively find all files named "subject"
    subject_files = glob.glob(os.path.join(args.input, "**", "subject"), recursive=True)
    print(subject_files)
    print(list_of_raw)

    subject_id = "##$SUBJECT_id="
    session_id = "##$SUBJECT_study_name="

    for subject_file in subject_files:
        if not os.path.exists(subject_file):
            continue

        with open(subject_file, 'r') as infile:
            lines = infile.readlines()

        modified = False

        for idx, line in enumerate(lines):
            # Modify both subject_id and study_name entries
            if subject_id in line or session_id in line:
                if idx + 1 < len(lines):
                    original_next = lines[idx + 1]

                    # 1) Apply selections to original underscore positions.
                    new_next = replace_underscores(original_next, replacements)

                    # 2) For study_name only: replace "baseline" with "PT0"
                    if session_id in line:
                        new_next = re.sub(r'<\s*baseline[^\s>]*\s*>', 'PT0', new_next, flags=re.IGNORECASE)

                    # Apply only if something changed
                    if new_next != original_next:
                        lines[idx + 1] = new_next
                        modified = True

        # Write back the modified content if changes were made
        if modified:
            with open(subject_file, 'w') as outfile:
                outfile.writelines(lines)
            print(f"Modified: {subject_file}")
            print("Success")
