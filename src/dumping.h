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

void dumpSingleChainToZip(zip_t *zip_archive, my_chain& chain, unsigned rep_id, const char* job_prefix, const char* cell_line, unsigned start, unsigned end);

void dumpEnsemble(my_ensemble& chains, const char* out_folder, const char* job_prefix);

void insertSampleData(const char* conninfo, my_chain& chain, unsigned start, unsigned end, unsigned rep_id, const char* job_prefix, const char* cell_line);

double calculateDistance(const Node node1, const Node node2);

std::vector<std::vector<double>> allDistanceMatrix(const std::vector<Node>& chain);

#endif //SBIF_DUMPING_H
