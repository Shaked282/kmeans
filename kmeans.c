#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 

double euclidean_distance(double vec1[], double vec2[], int num_of_coords){
    double d;
    int i;
    d = 0.0;
    for (i = 0; i < num_of_coords; i++){
        d += pow((vec1[i] - vec2[i]),2);
    }
    d = sqrt(d);
    return d;
}

void k_means(double** vectors, int num_of_vectors, int k, int num_of_coords, int iter){
    double** centroids;
    double*** clusters;
    int index_k, index_vector, index_coord;
    size_t index_size;
    int iter_cnt, min_cent_index;
    double min_dist, curr_dist;
    size_t* cluster_sizes;
    size_t new_size;
    double epsilon, curr_sum;
    double** new_centroids;
    int convergence;

    /* int l, m; */

    epsilon = 0.001;
    centroids = calloc(k, sizeof(double*));
    for (index_k = 0; index_k < k; index_k++) {
        centroids[index_k] = calloc(num_of_coords, sizeof(double));
        for (index_coord = 0; index_coord < num_of_coords; index_coord++) {
            centroids[index_k][index_coord] = vectors[index_k][index_coord];
        }
    }
    
    for (iter_cnt = 0; iter_cnt < iter; iter_cnt++) {
        new_centroids = calloc(k, sizeof(double*));
        for (index_k = 0; index_k < k; index_k++) {
            new_centroids[index_k] = calloc(num_of_coords, sizeof(double));
        }
        cluster_sizes = calloc(k, sizeof(size_t));
        clusters = calloc(k, sizeof(double**));
        for (index_k = 0; index_k < k; index_k++) {
            clusters[index_k] = calloc(1, sizeof(double*));
            clusters[index_k][0] = NULL;
        }

        printf("iter=%i\n",iter_cnt);
        for (index_vector = 0; index_vector < num_of_vectors; index_vector++) {
            min_dist = (double) INFINITY;
            min_cent_index = 0;
            for (index_k = 0; index_k < k; index_k++) {
                curr_dist = euclidean_distance(vectors[index_vector], centroids[index_k], num_of_coords);
                if (min_dist > curr_dist) {
                    min_cent_index = index_k;
                    min_dist = curr_dist;
                }
            }
            new_size = cluster_sizes[min_cent_index] + 1;
            if (clusters[min_cent_index][0] == NULL) {
                clusters[min_cent_index][0] = calloc(num_of_coords, sizeof(double));
                cluster_sizes[min_cent_index] = 1;
            } 
            else {
                clusters[min_cent_index] = realloc(clusters[min_cent_index], new_size*sizeof(double*));
                clusters[min_cent_index][new_size - 1] = calloc(num_of_coords, sizeof(double));
                cluster_sizes[min_cent_index] = new_size;
            }

            for (index_coord = 0; index_coord < num_of_coords; index_coord++) {
                clusters[min_cent_index][new_size-1][index_coord] = vectors[index_vector][index_coord];
            }
        }
        convergence = 0;
        for (index_k = 0; index_k < k; index_k++) {
            for (index_coord = 0; index_coord < num_of_coords; index_coord++) {
                curr_sum = 0.0;
                for (index_size = 0; index_size < cluster_sizes[index_k]; index_size++){
                    curr_sum += clusters[index_k][index_size][index_coord];
                }
                new_centroids[index_k][index_coord] = curr_sum / cluster_sizes[index_k];
            }
            /* for (l = 0; l<num_of_coords; l++){
                printf("old: %f new: %f\n",centroids[index_k][l], new_centroids[index_k][l]);
            } */
            if (euclidean_distance(centroids[index_k], new_centroids[index_k], num_of_coords) >= epsilon) {\
                convergence = 1;
            }
        }
/*         printf("new centroids!!!!\n");
        for (l = 0; l < k; l++){
            for (m = 0; m < num_of_coords; m++){
                printf("%f,",new_centroids[l][m]);
            }
            printf("\n");
        } */
        if (convergence == 0) {
            break;
        }
        if (iter_cnt != iter -1 && convergence != 0) {
            free(cluster_sizes);
        }
        for (index_k = 0; index_k < k; index_k++){
            for (index_coord = 0; index_coord < num_of_coords; index_coord++){
                centroids[index_k][index_coord] = new_centroids[index_k][index_coord];
            }
        }
        for (index_k = 0; index_k < k; index_k++) {
            free(new_centroids[index_k]);
        }
        free(new_centroids);

    }
    for (index_k = 0; index_k < k; index_k++){
        for (index_coord = 0; index_coord < num_of_coords; index_coord++) {
            printf("%.4f", centroids[index_k][index_coord]);
            if (index_coord != num_of_coords-1) {
                printf(",");
            }
        }
        printf("\n");
    }
    for (index_k = 0; index_k < k; index_k++) {
        for (index_size = 0; index_size < cluster_sizes[index_k]; index_size++) {
            if ( clusters[index_k][index_size] != NULL){
                free(clusters[index_k][index_size]);
            }
        }
        free(clusters[index_k]);
    }
    free(clusters);
    free(cluster_sizes);
    for (index_k = 0; index_k < k; index_k++) {
            free(centroids[index_k]);
        }
    free(centroids);
}

struct coord
{
    double value;
    struct coord *next;
};
struct vector
{
    struct vector *next;
    struct coord *coords;
};

void free_coords(struct coord *head)
{
    struct coord *tmp;
    while (head)
    {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

void free_vectors(struct vector *head)
{
    struct vector *tmp;
    while (head)
    {
        free_coords(head->coords);
        tmp = head;
        head = head->next;
        free(tmp);
    }
}

int main(int argc, char **argv)
{
    int K, iter;
    struct vector *head_vec, *curr_vec, *c_vec;
    struct coord *head_coord, *curr_coord, *c_coord;
    int num_of_vectors, num_of_coords;
    double n;
    char c;
    double kd, iterd;
    int vec_index, coor_index;
    double** vectors;

    if (argc < 2){
        printf("Invalid arguments!\n");
        exit(EXIT_FAILURE);
        return 1;
    }
    else if (argc == 2){
        K = atoi(argv[1]);
        kd = atof(argv[1]);
        iter = 200;
        iterd = 200;
    }
    else if (argc == 3){
        K = atoi(argv[1]);
        kd = atof(argv[1]);
        iter = atoi(argv[2]);
        iterd = atof(argv[2]);
    }
    printf("K=%i, iter=%i\n", K, iter);
    if (iter < 1 || iter > 1000 || iterd != iter) {
        printf("Invalid maximum iteration!\n");
        exit(EXIT_FAILURE);
        return 1;
    }
    head_coord = malloc(sizeof(struct coord));
    curr_coord = head_coord;
    curr_coord->next = NULL;
    head_vec = malloc(sizeof(struct vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;

    num_of_coords = 0;
    num_of_vectors = 0;
    c_vec = head_vec;
    while (scanf("%lf%c", &n, &c) == 2) {
        if (c == '\n') {
            if (num_of_vectors == 1) {
                ++num_of_coords;
            }
            num_of_vectors++;
            curr_coord->value = n;
            curr_vec->coords = head_coord;
            curr_vec->next = malloc(sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_coord = malloc(sizeof(struct coord));
            curr_coord = head_coord;
            curr_coord->next = NULL;
            continue;
        }
        if (num_of_vectors == 1){
            num_of_coords++;
        } 
        curr_coord->value = n;
        curr_coord->next = malloc(sizeof(struct coord));
        curr_coord = curr_coord->next;
        curr_coord->next = NULL;
    }

    if (K < 1 || K > num_of_vectors || K != kd) {
        printf("Invalid number of clusters!\n");
        exit(EXIT_FAILURE);
        return 1;
    }
    printf("num of vectors= %i\n", num_of_vectors);
    printf("num of coords= %i\n", num_of_coords);
    vectors = calloc(num_of_vectors, sizeof(double*));
    for (vec_index = 0; vec_index < num_of_vectors; vec_index++){
        c_coord = c_vec->coords;
        vectors[vec_index] = calloc(num_of_coords, sizeof(double));
        for (coor_index = 0; coor_index < num_of_coords; coor_index++){
            vectors[vec_index][coor_index] = c_coord->value;
            c_coord = c_coord->next;
        }
        c_vec = c_vec -> next;
    }

    k_means(vectors, num_of_vectors, K, num_of_coords, iter);
    for (vec_index = 0; vec_index < num_of_vectors; vec_index++){  /*!!!!!!!!!!!!*/
        free(vectors[vec_index]);
    }
    free(vectors);

    // free all linked list mallocs
    free_vectors(head_vec);
    return 0;
}
