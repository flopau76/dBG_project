#!/bin/bash

# is fed a fasta file
# each time I see the > I dump the sequence into a file
# each time it is not a > at the beginning, I append the sequence (without newline) to the current file
# at the end I dump the last file

# filename is the first record after the '>' in the fasta
#!/bin/bash

# Get the first command-line argument passed to the shell script
PATH_TO_TEMP_DIRECTORY="$1"

# Check if a parameter was provided (optional but good practice)
if [ -z "$PATH_TO_TEMP_DIRECTORY" ]; then
    echo "Usage: $0 <PATH_TO_TEMP_DIRECTORY> < input_file.txt > output_file.txt"
    exit 1
fi

awk -v temp_directory="$PATH_TO_TEMP_DIRECTORY" '
BEGIN {
        file_count = 0;
        file_name = "";
}

/^>/ {  
        if (file_name != "") {
            close(file_name);
        }

        file_name = substr($1,2);
        gsub(/[^a-zA-Z0-9_.-]/, "_", file_name)
        file_count += 1;
        print $0 > temp_directory "/" file_name ".fa";
        print temp_directory "/" file_name ".fa";
}

!/^>/ { 
print $0 >> temp_directory "/" file_name ".fa";
}

END { close(file_name)}
'