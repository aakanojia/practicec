#include <stdlib.h>

#include "global.h"
#include "debug.h"

// Compare string function prototype
int compare(const char *str1, const char *str2);


/**
 * @brief Validates command line arguments passed to the program.
 * @details This function will validate all the arguments passed to the
 * program, returning 0 if validation succeeds and -1 if validation fails.
 * Upon successful return, the various options that were specified will be
 * encoded in the global variable 'global_options', where it will be
 * accessible elsewhere in the program.  For details of the required
 * encoding, see the assignment handout.
 *
 * @param argc The number of arguments passed to the program from the CLI.
 * @param argv The argument strings passed to the program from the CLI.
 * @return 0 if validation succeeds and -1 if validation fails.
 * @modifies global variable "global_options" to contain an encoded representation
 * of the selected program options.
 */
int validargs(int argc, char **argv)
{
    // Initialize global_options to 0
    global_options = 0;

    // If no arguments or only the program name is provided, return -1
    if (argc == 1) { return -1;}

    // Check for the -h flag
    for (int i = 1; i < argc; i++) {
        if (compare(argv[i], "-h") == 0) {
            global_options |= HELP_OPTION;
            return 0;
        }
    }

    // Check for the -m and -n options
    for (int i = 1; i < argc; i++) {
        if (compare(argv[i], "-m") == 0) {
            global_options |= MATRIX_OPTION;
        } else if (compare(argv[i], "-n") == 0) {
            global_options |= NEWICK_OPTION;
            // If -n is followed by -o, store the outlier name
            if (i + 2 < argc && compare(argv[i + 1], "-o") == 0) {
                outlier_name = argv[i + 2];
                i += 2; // Skip the next two arguments
            }
        }
    }

    // Check for invalid flag combinations
    if ((global_options & MATRIX_OPTION) && (global_options & NEWICK_OPTION)) { return -1; }

    // If neither MATRIX_OPTION nor NEWICK_OPTION is set, return -1
    if (!(global_options & MATRIX_OPTION) && !(global_options & NEWICK_OPTION)) { return -1; }

    // After parsing all flags, check if -o is provided without -n
    if (outlier_name && !(global_options & NEWICK_OPTION)) { return -1; }

    return 0;
}

// Custom string comparison and manipulation function
int compare(const char *str1, const char *str2) {
    while (*str1 && (*str1 == *str2)) {
        str1++;
        str2++;
    }
    return *(unsigned char *)str1 - *(unsigned char *)str2;
}
