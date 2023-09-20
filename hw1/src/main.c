#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "debug.h"

#define MAX_NODES 100 // Maximum number of nodes

double matrix[MAX_NODES][MAX_NODES]; //Pre-declared 2D array to store the matrix

int readAndValidateMatrix() {
    int numNodes = 0; // To keep track of the number of nodes read.
    char line[1024]; // Buffer to read each line
}

int main(int argc, char **argv)
{
    if(validargs(argc, argv))
        USAGE(*argv, EXIT_FAILURE);
    if(global_options == HELP_OPTION)
        USAGE(*argv, EXIT_SUCCESS);
    // Check for MATRIX_OPTION
    if(global_options & MATRIX_OPTION) {
        // Implement the logic for the matrix mode
        // Read from stdin, process, write to stdout
    }

    // Check for NEWICK_OPTION
    if(global_options & NEWICK_OPTION) {
        // Implement the logic for the Newick mode
        // If outlier_name is set, you also have the `-o` flag provided
        // Read from stdin, process, write to stdout
    }

    // If everything goes well
    return EXIT_SUCCESS;

    // If there's an error during processing, before the return statement above, you can:
    // fprintf(stderr, "Error message here");
    // return EXIT_FAILURE;
    return EXIT_FAILURE;
}

/*
 * Just a reminder: All non-main functions should
 * be in another file not named main.c
 */
