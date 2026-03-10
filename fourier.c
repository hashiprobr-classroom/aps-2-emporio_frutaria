#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}

void fft(complex s[], complex t[], int n, int sign) {
    
    if (n==1){
        t[0] = s[0];
        return;
    }

    int metade = n / 2;
    complex sp[metade];
    complex si[metade];

    for (int i = 0; i < metade; i++){
        sp[i] = s[i*2];
        si[i] = s[i*2 + 1];

    }

    complex tp[metade];
    complex ti[metade];

    fft(sp, tp, metade, sign);
    fft(si, ti, metade, sign);

    for (int k = 0; k < metade; k++){
        double x = sign * 2 * PI * k / n;
        
        double term_a = ti[k].a * cos(x) - ti[k].b * sin(x);
        double term_b = ti[k].a * sin(x) + ti[k].b * cos(x);

        t[k].a = tp[k].a + term_a;
        t[k].b = tp[k].b + term_b;

        t[k + n / 2].a = tp[k].a - term_a;
        t[k + n / 2].b = tp[k].b - term_b;
    }

    
}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for (int i = 0; i < width; i++){
        complex s[height];
        complex t[height];

        for (int j = 0; j < height; j++){
            s[j] = matrix[j][i];
        }

        fft_forward(s, t, height);

        for (int k = 0; k < height; k++){
            matrix[k][i] = t[k];
        }
    }

    for (int l = 0; l < height; l++){
        complex s[width];
        complex t[width];

        for (int m = 0; m < width; m++){
            s[m] = matrix[l][m];
        }

        fft_forward(s, t, width);

        for (int n = 0; n < width; n++){
            matrix[l][n] = t[n];
        }

    }
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    for (int i = 0; i < height; i++){
        complex s[width];
        complex t[width];

        for (int j = 0; j < width; j++){
            s[j] = matrix[i][j];
        }

        fft_inverse(s, t, width);

        for (int k = 0; k < width; k++){
            matrix[i][k] = t[k];
        }
    }

    for (int l = 0; l < width; l++){
        complex s[height];
        complex t[height];

        for (int m = 0; m < height; m++){
            s[m] = matrix[m][l];
        }

        fft_inverse(s, t, height);

        for (int n = 0; n < height; n++){
            matrix[n][l] = t[n];
        }

    }
}

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
