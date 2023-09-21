#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "debug.h"

#define MAX_TAXA 100
#define MAX_NEWICK_SIZE 4096

char taxa_names[MAX_TAXA][MAX_TAXA];
int distance_matrix[MAX_TAXA][MAX_TAXA];
int flags = 0; // Placeholder for flags
int num_all_nodes = MAX_TAXA;
char* outlier_name = NULL; // Placeholder for outlier name. This would be set based on user input

int compare(const char *str1, const char *str2);

int parse_integer(const char *str, int *value) {
    int i = 0;
    *value = 0;
    while (str[i] >= '0' && str[i] <= '9') {
        *value = (*value * 10) + (str[i] - '0');
        i++;
    }
    return i;
}

// void append_to_newick(char* dest, const char* src) {
//     while (*dest) {
//         dest++
//     }
//     while (*src) {
//         *dest = *src;
//         dest++;
//         src++;
//     }
//     *dest = '\0';
// }

// char* generate_newick_recursive(int node, char* buffer) {
//     // If the node is a leaf, append its name to the buffer
//     if (node < num_taxa) {
//         append_to_newick(buffer, taxa_names[node]);
//         return;
//     }

//     // Otherwise, it's an internal node. Recursively generate the Newick format for its children.
//     char left_newick[MAX_NEWICK_SIZE] = {0};
//     char right_newick[MAX_NEWICK_SIZE] = {0};
//     generate_newick_recursive(nodes[node].neighbors[0]->index, left_newick);
//     generate_newick_recursive(nodes[node].neighbors[1]->index, right_newick);

//     // Calculate the distances to the children
//     double left_distance = distances[node][nodes[node].neighbors[0]->index];
//     double right_distance = distances[node][nodes[node].neighbors[1]->index];

//     // Append the combined Newick formats of the children to the buffer
//     append_to_newick(buffer, "(");
//     append_to_newick(buffer, left_newick);
//     char left_dist_str[20]; // Assuming distance won't exceed this length
//     sprintf(left_dist_str, ":%.2f", left_distance);
//     append_to_newick(buffer, left_dist_str);
//     append_to_newick(buffer, ", ");
//     append_to_newick(buffer, right_newick);
//     char right_dist_str[20];
//     sprintf(right_dist_str, ":%.2f", right_distance);
//     append_to_newick(buffer, right_dist_str);
//     append_to_newick(buffer, ")");
// }

void generate_newick_recursive(NODE *node, char *buffer, int *pos) {
    if (node == NULL) return;

    // If it's a leaf node
    if (node->neighbors[1] == NULL && node->neighbors[2] == NULL) {
        // Append the name of the leaf node to the buffer
        int i = 0;
        while (node->name[i] != '\0') {
            buffer[*pos] = node->name[i];
            (*pos)++;
            i++;
        }
        return;
    }

    // If it's an internal node
    buffer[*pos] = '(';
    (*pos)++;

    // Recur for the first child
    generate_newick_recursive(node->neighbors[1], buffer, pos);

    buffer[*pos] = ',';
    (*pos)++;

    // Recur for the second child
    generate_newick_recursive(node->neighbors[2], buffer, pos);

    buffer[*pos] = ')';
    (*pos)++;
}


/**
 * @brief  Read genetic distance data and initialize data structures.
 * @details  This function reads genetic distance data from a specified
 * input stream, parses and validates it, and initializes internal data
 * structures.
 *
 * The input format is a simplified version of Comma Separated Values
 * (CSV).  Each line consists of text characters, terminated by a newline.
 * Lines that start with '#' are considered comments and are ignored.
 * Each non-comment line consists of a nonempty sequence of data fields;
 * each field iis terminated either by ',' or else newline for the last
 * field on a line.  The constant INPUT_MAX specifies the maximum number
 * of data characters that may be in an input field; fields with more than
 * that many characters are regarded as invalid input and cause an error
 * return.  The first field of the first data line is empty;
 * the subsequent fields on that line specify names of "taxa", which comprise
 * the leaf nodes of a phylogenetic tree.  The total number N of taxa is
 * equal to the number of fields on the first data line, minus one (for the
 * blank first field).  Following the first data line are N additional lines.
 * Each of these lines has N+1 fields.  The first field is a taxon name,
 * which must match the name in the corresponding column of the first line.
 * The subsequent fields are numeric fields that specify N "distances"
 * between this taxon and the others.  Any additional lines of input following
 * the last data line are ignored.  The distance data must form a symmetric
 * matrix (i.e. D[i][j] == D[j][i]) with zeroes on the main diagonal
 * (i.e. D[i][i] == 0).
 *
 * If 0 is returned, indicating data successfully read, then upon return
 * the following global variables and data structures have been set:
 *   num_taxa - set to the number N of taxa, determined from the first data line
 *   num_all_nodes - initialized to be equal to num_taxa
 *   num_active_nodes - initialized to be equal to num_taxa
 *   node_names - the first N entries contain the N taxa names, as C strings
 *   distances - initialized to an NxN matrix of distance values, where each
 *     row of the matrix contains the distance data from one of the data lines
 *   nodes - the "name" fields of the first N entries have been initialized
 *     with pointers to the corresponding taxa names stored in the node_names
 *     array.
 *   active_node_map - initialized to the identity mapping on [0..N);
 *     that is, active_node_map[i] == i for 0 <= i < N.
 *
 * @param in  The input stream from which to read the data.
 * @return 0 in case the data was successfully read, otherwise -1
 * if there was any error.  Premature termination of the input data,
 * failure of each line to have the same number of fields, and distance
 * fields that are not in numeric format should cause a one-line error
 * message to be printed to stderr and -1 to be returned.
 */

int read_distance_data(FILE *in) {
    char line[1024];
    int row = 0;

    // Read the first line for taxa names
    if (fgets(line, sizeof(line), in) != NULL) {
        int i = 0, j = 0, k = 0;
        while (line[i] != '\0') {
            if (line[i] == ',' || line[i] == '\n') {
                taxa_names[j][k] = '\0';
                j++;
                k = 0;
                i++;
            } else {
                taxa_names[j][k++] = line[i++];
            }
        }
    }

    // Read the distance matrix
    while (fgets(line, sizeof(line), in) != NULL) {
        int i = 0, col = 0, value;
        while (line[i] != '\0') {
            if (line[i] == ',' || line[i] == '\n') {
                i++;
            } else {
                i += parse_integer(&line[i], &value);
                distance_matrix[row][col++] = value;
            }
        }
        row++;
    }

    // TODO: Add error handling and validation

    return 0;
}

/**
 * @brief  Emit a representation of the phylogenetic tree in Newick
 * format to a specified output stream.
 * @details  This function emits a representation in Newick format
 * of a synthesized phylogenetic tree to a specified output stream.
 * See (https://en.wikipedia.org/wiki/Newick_format) for a description
 * of Newick format.  The tree that is output will include for each
 * node the name of that node and the edge distance from that node
 * its parent.  Note that Newick format basically is only applicable
 * to rooted trees, whereas the trees constructed by the neighbor
 * joining method are unrooted.  In order to turn an unrooted tree
 * into a rooted one, a root will be identified according by the
 * following method: one of the original leaf nodes will be designated
 * as the "outlier" and the unique node adjacent to the outlier
 * will serve as the root of the tree.  Then for any other two nodes
 * adjacent in the tree, the node closer to the root will be regarded
 * as the "parent" and the node farther from the root as a "child".
 * The outlier node itself will not be included as part of the rooted
 * tree that is output.  The node to be used as the outlier will be
 * determined as follows:  If the global variable "outlier_name" is
 * non-NULL, then the leaf node having that name will be used as
 * the outlier.  If the value of "outlier_name" is NULL, then the
 * leaf node having the greatest total distance to the other leaves
 * will be used as the outlier.
 *
 * @param out  Stream to which to output a rooted tree represented in
 * Newick format.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.  If the global variable "outlier_name" is
 * non-NULL, then it is an error if no leaf node with that name exists
 * in the tree.
 */
int emit_newick_format(FILE *out) {
    char newick_str[MAX_NEWICK_SIZE] = {0}; // Static buffer for Newick string
    int pos = 0;

    // Find the root node (it could be any node based on your logic)
    NODE *root = &nodes[0]; // Assuming nodes[0] is the root

    generate_newick_recursive(root, newick_str, &pos);

    newick_str[pos] = ';';
    newick_str[pos + 1] = '\0';

    fprintf(out, "%s\n", newick_str);

    return 0;
}

// int emit_newick_format(FILE *out) {
//     int outlier_node = -1;
//     char newick_str[MAX_NEWICK_SIZE] = {0}; // Static buffer for Newick string
//     if (outlier_name) {
//         // Find the node with the outlier_name
//         for (int i = 0; i < num_all_nodes; i++) {
//             if (compare(taxa_names[i], outlier_name) == 0) {
//                 outlier_node = i;
//                 break;
//             }
//         }
//         if (outlier_node == -1) return -1; // Error if no node with outlier_name exists
//     } else {
//         // Determine the outlier node based on the greatest total distance to other leaves
//         int max_distance = -1;
//         for (int i = 0; i < num_all_nodes; i++) {
//             int current_distance = 0;
//             for (int j = 0; j < num_all_nodes; j++) {
//                 current_distance += distance_matrix[i][j];
//             }
//             if (current_distance > max_distance) {
//                 max_distance = current_distance;
//                 outlier_node = i;
//             }
//         }
//     }

//     if (newick_str[0] != '\0') {
//         fprintf(out, "%s;\n", newick_str);
//     }
//     return 0;
// }

/**
 * @brief  Emit the synthesized distance matrix as CSV.
 * @details  This function emits to a specified output stream a representation
 * of the synthesized distance matrix resulting from the neighbor joining
 * algorithm.  The output is in the same CSV form as the program input.
 * The number of rows and columns of the matrix is equal to the value
 * of num_all_nodes at the end of execution of the algorithm.
 * The submatrix that consists of the first num_leaves rows and columns
 * is identical to the matrix given as input.  The remaining rows and columns
 * contain estimated distances to internal nodes that were synthesized during
 * the execution of the algorithm.
 *
 * @param out  Stream to which to output a CSV representation of the
 * synthesized distance matrix.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int emit_distance_matrix(FILE *out) {
    for (int i = 0; i < num_all_nodes; i++) {
        for (int j = 0; j < num_all_nodes; j++) {
            fprintf(out, "%d", distance_matrix[i][j]);
            if (j < num_all_nodes - 1) {
                fprintf(out, ",");  // Add comma for CSV format except for the last column
            }
        }
    }
    return 0;
}

/**
 * @brief  Build a phylogenetic tree using the distance data read by
 * a prior successful invocation of read_distance_data().
 * @details  This function assumes that global variables and data
 * structures have been initialized by a prior successful call to
 * read_distance_data(), in accordance with the specification for
 * that function.  The "neighbor joining" method is used to reconstruct
 * phylogenetic tree from the distance data.  The resulting tree is
 * an unrooted binary tree having the N taxa from the original input
 * as its leaf nodes, and if (N > 2) having in addition N-2 synthesized
 * internal nodes, each of which is adjacent to exactly three other
 * nodes (leaf or internal) in the tree.  As each internal node is
 * synthesized, information about the edges connecting it to other
 * nodes is output.  Each line of output describes one edge and
 * consists of three comma-separated fields.  The first two fields
 * give the names of the nodes that are connected by the edge.
 * The third field gives the distance that has been estimated for
 * this edge by the neighbor-joining method.  After N-2 internal
 * nodes have been synthesized and 2*(N-2) corresponding edges have
 * been output, one final edge is output that connects the two
 * internal nodes that still have only two neighbors at the end of
 * the algorithm.  In the degenerate case of N=1 leaf, the tree
 * consists of a single leaf node and no edges are output.  In the
 * case of N=2 leaves, then no internal nodes are synthesized and
 * just one edge is output that connects the two leaves.
 *
 * Besides emitting edge data (unless it has been suppressed),
 * as the tree is built a representation of it is constructed using
 * the NODE structures in the nodes array.  By the time this function
 * returns, the "neighbors" array for each node will have been
 * initialized with pointers to the NODE structure(s) for each of
 * its adjacent nodes.  Entries with indices less than N correspond
 * to leaf nodes and for these only the neighbors[0] entry will be
 * non-NULL.  Entries with indices greater than or equal to N
 * correspond to internal nodes and each of these will have non-NULL
 * pointers in all three entries of its neighbors array.
 * In addition, the "name" field each NODE structure will contain a
 * pointer to the name of that node (which is stored in the corresponding
 * entry of the node_names array).
 *
 * @param out  If non-NULL, an output stream to which to emit the edge data.
 * If NULL, then no edge data is output.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int build_taxonomy(FILE *out) {
    // Initialization
    int total_distances[MAX_TAXA] = {0};
    for (int i = 0; i < MAX_TAXA; i++) {
        for (int j = 0; j < MAX_TAXA; j++) {
            total_distances[i] += distance_matrix[i][j];
        }
    }

    // Main loop for the neighbor joining algorithm
    int N = MAX_TAXA; // Number of remaining taxa
    while (N > 2) {
        // Find the pair of taxa to join
        int min_i = -1, min_j = -1;
        double min_distance = 1.0e30; // A large initial value
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double reduced_distance = distance_matrix[i][j] - (double)(total_distances[i] + total_distances[j]) / (N - 2);
                if (reduced_distance < min_distance) {
                    min_distance = reduced_distance;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        // Calculate the branch lengths and update the distance matrix
        double branch_length_i = (distance_matrix[min_i][min_j] + (double)(total_distances[min_i] - total_distances[min_j]) / (N - 2)) / 2;
        double branch_length_j = distance_matrix[min_i][min_j] - branch_length_i;

        // Update the distance matrix
        for (int k = 0; k < N; k++) {
            if (k != min_i && k != min_j) {
                double new_distance = (distance_matrix[min_i][k] + distance_matrix[min_j][k] - distance_matrix[min_i][min_j]) / 2;
                distance_matrix[min_i][k] = new_distance;
                distance_matrix[k][min_i] = new_distance;
            }
        }
        // Remove min_j from the matrix by shifting rows and columns
        for (int k = min_j; k < N - 1; k++) {
            total_distances[k] = total_distances[k + 1];
            for (int l = 0; l < N; l++) {
                distance_matrix[k][l] = distance_matrix[k + 1][l];
                distance_matrix[l][k] = distance_matrix[l][k + 1];
            }
        }
        N--;
    }

    // Output the resulting tree or matrix based on flags
    if (flags & 0x1) { // 0x1 represents the -m flag
        //Output the matrix of estimated distances
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(out, "%d", distance_matrix[i][j]);
            }
            fprintf(out, "\n");
        }
    } else {
        // Output the edges (default behavior)
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                fprintf(out, "%s  %s %d\n", taxa_names[i], taxa_names[j], distance_matrix[i][j]);
            }
        }
    }

    // for (int i = 0; i < N; i++) {
    //     for (int j = i + 1; j < N; j++) {
    //         fprintf(out, "%s %s %d\n", taxa_names[i], taxa_names[j], distance_matrix[i][j]);
    //     }
    // }

    return 0;
}
