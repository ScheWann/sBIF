#ifndef SBIF_DUMPING_H
#define SBIF_DUMPING_H

#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include "loading.h"
#include "check.h"
#include <libpq-fe.h>

char* getOutfile(const char* out_folder, unsigned rep_id, const char* job_prefix);

void dumpSingleChain(my_chain& chain, const char* out_folder, unsigned rep_id, const char* job_prefix);

void dumpEnsemble(my_ensemble& chains, const char* out_folder, const char* job_prefix);

void insertSampleData(const char* conninfo, my_chain& chain, unsigned start, unsigned end, unsigned rep_id, const char* job_prefix);

std::vector<std::vector<double>> allDistanceMatrix(const std::vector<Node>& chain);

// std::vector<char> createZipInMemory(const std::vector<my_chain>& chains, const std::string& job_prefix);
// void insertExampleData(my_chain& chain, unsigned start, unsigned end, unsigned rep_id, const char* job_prefix);
#endif //SBIF_DUMPING_H
