#ifndef SBIF_DUMPING_H
#define SBIF_DUMPING_H

#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include "loading.h"
#include "check.h"
#include <libpq-fe.h>
#include <zlib.h>
#include <zip.h>

char* getOutfile(const char* out_folder, unsigned rep_id, const char* job_prefix);

// void insertDistanceData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix, const char *cell_line);

std::vector<float> computeDistanceVector(const my_chain &chain);

void insertDistanceDataFromVector(const char *conninfo, const char *cell_line, const char *chrid, unsigned rep_id, unsigned start, unsigned end, const std::vector<float> &distances_f);

std::vector<float> computeAvgVector(const std::vector<std::pair<unsigned, std::vector<float>>> &all_distances);

std::vector<float> computeFreqCondensed(const std::vector<std::pair<unsigned, std::vector<float>>> &all_distances, float threshold);

std::vector<float> squareformFullMatrix(const std::vector<float> &condensed);

void insertCalcDistance(const char *conninfo, const char *cell_line, const char *chrid, unsigned start, unsigned end, const std::vector<float> &avg_vec, const std::vector<float> &fq_full);

void dumpEnsemble(my_ensemble& chains, const char* out_folder, const char* job_prefix);

void insertSampleData(const char* conninfo, my_chain& chain, unsigned start, unsigned end, unsigned rep_id, const char* job_prefix, const char* cell_line);

double calculateDistance(const Node node1, const Node node2);

#endif //SBIF_DUMPING_H
