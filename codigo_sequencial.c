#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#define N 2000  // Tamanho da grade
#define T 501 // Número de iterações no tempo
#define D 0.1  // Coeficiente de difusão
#define DELTA_T 0.01 // Intervalo de tempo entre iterações
#define DELTA_X 1.0 // Espaçamento entre os pontos da grade


// Função que resolve a equação de difusão
void diff_eq(double **C, double **C_new) { 
    for (int t = 0; t < T; t++) {
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
              // Atualiza a concentração de cada ponto usando a equação de difusão
                C_new[i][j] = C[i][j] + D * DELTA_T * (
                    (C[i+1][j] + C[i-1][j] + C[i][j+1] + C[i][j-1] - 4 * C[i][j]) / (DELTA_X * DELTA_X)
                );
            }
        }

        // // Calcula a diferença média entre as iterações e atualiza a matriz C
        double difmedio = 0.;
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                difmedio += fabs(C_new[i][j] - C[i][j]); // Soma da diferença absoluta
                C[i][j] = C_new[i][j]; // Atualiza a matriz C com os novos valores
            }
        }

        // Imprime a diferença média a cada 100 iterações
        if ((t%100) == 0)
          printf("interacao %d - diferenca=%g\n", t, difmedio/((N-2)*(N-2)));
    }
}

int main() {

    double tempo_inicial, tempo_final;

    tempo_inicial = omp_get_wtime();

    // Aloca memória para a matriz de concentração atual (C)
    double **C = (double **)malloc(N * sizeof(double *));
    if (C == NULL) { 
      fprintf(stderr, "Memory allocation failed\n");
      return 1;
    }
    for (int i = 0; i < N; i++) {
      C[i] = (double *)malloc(N * sizeof(double));
      if (C[i] == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
      }
    }
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        C[i][j] = 0.;
      }
    }

    // Concentração para a próxima iteração
    double **C_new = (double **)malloc(N * sizeof(double *));
    if (C_new == NULL) {
      fprintf(stderr, "Memory allocation failed\n");
      return 1;
    }
    for (int i = 0; i < N; i++) {
      C_new[i] = (double *)malloc(N * sizeof(double));
      if (C_new[i] == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
      }
    }
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        C_new[i][j] = 0.;
      }
    }

    // Configura a concentração inicial
    C[N/2][N/2] = 1.0;

    // Chama a função para resolver a equação de difusão
    diff_eq(C, C_new);

    tempo_final = omp_get_wtime();

    printf("Tempo total de execução: %f segundos\n", tempo_final - tempo_inicial);

    // Exibe o valor da concentração final no centro da matriz para verificação
    printf("Concentração final no centro: %f\n", C[N/2][N/2]);
    return 0;
}