#pragma once

// read points from input file
void read_points(const char* filename, int* nimg, int* npts, double** pts);

// read table from input file
void read_table(const char* filename, int* nrow, int* ncol, double** tab);
