#include "dumping.h"
#include "check.h"
#include <libpq-fe.h>
#include <sstream>
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <zlib.h>
#include <zip.h>
#include <filesystem>
#include <numeric>
#include <cmath>

char *getOutfile(const char *out_folder, unsigned rep_id, const char *job_prefix)
{
    const char *seperator = "/";
    char rep[10];
    sprintf(rep, "%d", rep_id);
    char *out_file = (char *)malloc(sizeof(char) * MAX_CHAR);
    strcpy(out_file, out_folder);
    if (out_folder[strlen(out_folder) - 1] != *seperator)
    {
        strcat(out_file, seperator);
    }
    strcat(out_file, job_prefix);
    strcat(out_file, ".");
    strcat(out_file, rep);
    strcat(out_file, ".txt");
    return out_file;
}

// Euclidean distance calculation function
// double calculateDistance(const Node node1, const Node node2)
// {
//     return std::sqrt(
//         std::pow(node2.x - node1.x, 2) +
//         std::pow(node2.y - node1.y, 2) +
//         std::pow(node2.z - node1.z, 2));
// }

// insert into the database
std::string join(const std::vector<std::string> &elements, const std::string &delimiter)
{
    std::ostringstream oss;
    for (size_t i = 0; i < elements.size(); ++i)
    {
        oss << elements[i];
        if (i < elements.size() - 1)
        {
            oss << delimiter;
        }
    }
    return oss.str();
}

std::vector<float> computeDistanceVector(const my_chain &chain) {
    int n_beads = static_cast<int>(chain.size());
    size_t num_vals = static_cast<size_t>(n_beads) * (n_beads - 1) / 2;
    std::vector<float> distances;
    distances.reserve(num_vals);

    for (int i = 0; i < n_beads; ++i) {
        for (int j = i + 1; j < n_beads; ++j) {
            // double d = calculateDistance(chain[i], chain[j]);
            double d = calcuDist(chain[i], chain[j]);
            distances.push_back(static_cast<float>(d));
        }
    }
    return distances;
}

void insertDistanceDataFromVector(PGconn *conn, const char *cell_line, const char *chrid, unsigned rep_id, unsigned start, unsigned end, const std::vector<float> &distances_f)
{
    // Ensure the passed-in connection is valid
    if (PQstatus(conn) != CONNECTION_OK) {
        std::cerr << "connected databse failed " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    int n_beads = static_cast<int>((1 + std::sqrt(1 + 8.0 * distances_f.size())) / 2); 
    size_t bin_len = distances_f.size() * sizeof(float);
    std::string binBuf;
    binBuf.resize(bin_len);
    std::memcpy(&binBuf[0], distances_f.data(), bin_len);

    std::string str_rep     = std::to_string(rep_id);
    std::string str_start   = std::to_string(start);
    std::string str_end     = std::to_string(end);
    std::string str_n_beads = std::to_string(n_beads);

    // time stamp
    char time_buffer[20];
    {
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm *now_tm = std::localtime(&now_c);
        std::strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", now_tm);
    }
    std::string str_time = time_buffer;

    const int nParams = 8;
    const char *paramValues[nParams];
    int paramLengths[nParams];
    int paramFormats[nParams];

    // cell_line
    paramValues[0]  = cell_line;
    paramLengths[0] = static_cast<int>(std::strlen(cell_line));
    paramFormats[0] = 0;

    // chrid
    paramValues[1]  = chrid;
    paramLengths[1] = static_cast<int>(std::strlen(chrid));
    paramFormats[1] = 0;

    // sampleid
    paramValues[2]  = str_rep.c_str();
    paramLengths[2] = static_cast<int>(str_rep.size());
    paramFormats[2] = 0;

    // start_value
    paramValues[3]  = str_start.c_str();
    paramLengths[3] = static_cast<int>(str_start.size());
    paramFormats[3] = 0;

    // end_value
    paramValues[4]  = str_end.c_str();
    paramLengths[4] = static_cast<int>(str_end.size());
    paramFormats[4] = 0;

    // n_beads
    paramValues[5]  = str_n_beads.c_str();
    paramLengths[5] = static_cast<int>(str_n_beads.size());
    paramFormats[5] = 0;

    // distance_vector
    paramValues[6]  = binBuf.data();
    paramLengths[6] = static_cast<int>(bin_len);
    paramFormats[6] = 1;

    // insert_time
    paramValues[7]  = str_time.c_str();
    paramLengths[7] = static_cast<int>(str_time.size());
    paramFormats[7] = 0;

    PGresult *res = PQexecParams(conn,
        "INSERT INTO distance("
            "cell_line, chrid, sampleid, start_value, end_value, n_beads, distance_vector, insert_time"
        ") VALUES("
            "$1, $2, $3, $4, $5, $6, $7, $8"
        ");",
        nParams,
        nullptr,
        paramValues,
        paramLengths,
        paramFormats,
        0
    );

    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "insert distance failed: " << PQerrorMessage(conn) << std::endl;
    }
    PQclear(res);
}

// Compute the average distance vector from all distances collected
std::vector<float> computeAvgVector(const std::vector<std::pair<unsigned, std::vector<float>>> &all_distances)
{
    size_t n_samples = all_distances.size();
    if (n_samples == 0) return {};

    size_t L = all_distances[0].second.size();
    std::vector<double> sum_vec(L, 0.0);

    for (const auto &p : all_distances) {
        const std::vector<float> &dist = p.second;
        for (size_t i = 0; i < L; ++i) {
            sum_vec[i] += static_cast<double>(dist[i]);
        }
    }
    std::vector<float> avg_vec(L);
    for (size_t i = 0; i < L; ++i) {
        avg_vec[i] = static_cast<float>(sum_vec[i] / static_cast<double>(n_samples));
    }
    return avg_vec;
}

// Compute the fq condensed vector
std::vector<float> computeFreqCondensed(const std::vector<std::pair<unsigned, std::vector<float>>> &all_distances, float threshold)
{
    size_t n_samples = all_distances.size();
    if (n_samples == 0) return {};

    size_t L = all_distances[0].second.size();
    std::vector<int> count_leq(L, 0);

    for (const auto &p : all_distances) {
        const std::vector<float> &dist = p.second;
        for (size_t i = 0; i < L; ++i) {
            if (dist[i] <= threshold) {
                count_leq[i] += 1;
            }
        }
    }
    std::vector<float> freq_cond(L);
    for (size_t i = 0; i < L; ++i) {
        freq_cond[i] = static_cast<float>(count_leq[i]) / static_cast<float>(n_samples);
    }
    return freq_cond;
}

std::vector<float> squareformFullMatrix(const std::vector<float> &condensed)
{
    size_t L = condensed.size();
    size_t n_beads = static_cast<size_t>((1 + std::sqrt(1 + 8.0 * static_cast<double>(L))) / 2.0);

    std::vector<float> full_mat(n_beads * n_beads, 0.0f);
    for (size_t i = 0; i < n_beads; ++i) {
        for (size_t j = i + 1; j < n_beads; ++j) {
            size_t idx = i * n_beads - (i * (i + 1) / 2) + (j - i - 1);
            float v = condensed[idx];
            full_mat[i * n_beads + j] = v;
            full_mat[j * n_beads + i] = v;
        }
    }
    return full_mat;
}

std::pair<std::vector<float>, unsigned> computeBestVector(
    const std::vector<std::pair<unsigned,std::vector<float>>>& all_distances,
    const std::vector<float>& avg_vec)
{
    size_t n = all_distances.size(), L = avg_vec.size();
    // center the average vector
    double avg_mean = std::accumulate(avg_vec.begin(), avg_vec.end(), 0.0) / L;
    std::vector<double> avg_centered(L);
    for (size_t i = 0; i < L; ++i)
        avg_centered[i] = avg_vec[i] - avg_mean;
    double avg_norm_sq = std::accumulate(
        avg_centered.begin(), avg_centered.end(), 0.0,
        [](double s, double v){ return s + v*v; });
    double avg_norm = std::sqrt(avg_norm_sq);

    // calculate the best correlation
    size_t best_idx = 0;
    double best_corr = 0.0;
    for (size_t i = 0; i < n; ++i) {
        auto &vec = all_distances[i].second;
        double mean_i = std::accumulate(vec.begin(), vec.end(), 0.0) / L;

        double numer = 0, norm_i_sq = 0;
        for (size_t j = 0; j < L; ++j) {
            double x_cent = vec[j] - mean_i;
            numer    += x_cent * avg_centered[j];
            norm_i_sq += x_cent * x_cent;
        }
        double norm_i = std::sqrt(norm_i_sq);
        double corr = (norm_i * avg_norm == 0.0)
                      ? 0.0
                      : numer / (norm_i * avg_norm);
        if (std::abs(corr) > std::abs(best_corr)) {
            best_corr = corr;
            best_idx  = i;
        }
    }
    return std::make_pair(all_distances[best_idx].second, all_distances[best_idx].first);
}

void insertCalcDistance(const char *conninfo, const char *cell_line, const char *chrid, unsigned start, unsigned end, const std::vector<float> &avg_vec, const std::vector<float> &fq_full, const std::vector<float> &best_vec, unsigned best_sample_id)
{
    PGconn *conn = PQconnectdb(conninfo);
    if (PQstatus(conn) != CONNECTION_OK) {
        std::cerr << "connect calc_distance failed: " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    const int nParams = 8;
    const char *paramValues[nParams];
    int paramLengths[nParams];
    int paramFormats[nParams];

    char time_buffer[20];
    {
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm *now_tm = std::localtime(&now_c);
        std::strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", now_tm);
    }

    // (1) cell_line
    paramValues[0]  = cell_line;
    paramLengths[0] = static_cast<int>(std::strlen(cell_line));
    paramFormats[0] = 0;

    // (2) chrid
    paramValues[1]  = chrid;
    paramLengths[1] = static_cast<int>(std::strlen(chrid));
    paramFormats[1] = 0;

    // (3) start
    std::string str_start = std::to_string(start);
    paramValues[2]  = str_start.c_str();
    paramLengths[2] = static_cast<int>(str_start.size());
    paramFormats[2] = 0;

    // (4) end
    std::string str_end = std::to_string(end);
    paramValues[3]  = str_end.c_str();
    paramLengths[3] = static_cast<int>(str_end.size());
    paramFormats[3] = 0;

    // (5) avg_vec
    size_t avg_len = avg_vec.size() * sizeof(float);
    std::string avgBuf;
    avgBuf.resize(avg_len);
    std::memcpy(&avgBuf[0], avg_vec.data(), avg_len);
    paramValues[4]  = avgBuf.data();
    paramLengths[4] = static_cast<int>(avg_len);
    paramFormats[4] = 1;

    // (6) fq_full
    size_t fq_len = fq_full.size() * sizeof(float);
    std::string fqBuf;
    fqBuf.resize(fq_len);
    std::memcpy(&fqBuf[0], fq_full.data(), fq_len);
    paramValues[5]  = fqBuf.data();
    paramLengths[5] = static_cast<int>(fq_len);
    paramFormats[5] = 1;

    // (7) best_vec
    size_t best_len = best_vec.size() * sizeof(float);
    std::string bestBuf;
    bestBuf.resize(best_len);
    std::memcpy(&bestBuf[0], best_vec.data(), best_len);
    paramValues[6]  = bestBuf.data();
    paramLengths[6] = static_cast<int>(best_len);
    paramFormats[6] = 1;

    // (8) best_sample_id
    std::string str_best_sample_id = std::to_string(best_sample_id);
    paramValues[7]  = str_best_sample_id.c_str();
    paramLengths[7] = static_cast<int>(str_best_sample_id.size());
    paramFormats[7] = 0;

    PGresult *res = PQexecParams(conn,
        "INSERT INTO calc_distance ("
            "cell_line, chrid, start_value, end_value, avg_distance_vector, fq_distance_vector, best_vector, best_sample_id"
        ") VALUES ("
            "$1, $2, $3, $4, $5, $6, $7, $8"
        ");",
        nParams,
        nullptr,
        paramValues,
        paramLengths,
        paramFormats,
        0
    );

    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "insert calc_distance failed: " << PQerrorMessage(conn) << std::endl;
    }

    PQclear(res);
    PQfinish(conn);
}

void insertSampleData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix, const char *cell_line)
{
    // Establish a new connection for this thread
    PGconn *conn = PQconnectdb(conninfo);
    if (PQstatus(conn) != CONNECTION_OK)
    {
        std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    // Prepare the insert query
    std::string insertQuery = "INSERT INTO position (cell_line, chrID, sampleID, X, Y, Z, start_value, end_value, insert_time) VALUES ";
    std::vector<std::string> valueSets;

    // Get the current local time
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm *now_tm = std::localtime(&now_c);

    // Format the current time as YYYY-MM-DD HH:MM:SS using strftime
    char time_buffer[20]; // Buffer to hold formatted time
    std::strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", now_tm);
    std::string insertTime(time_buffer);

    for (Node node : chain)
    {
        char rep_id_str[12], x_str[32], y_str[32], z_str[32];
        snprintf(rep_id_str, sizeof(rep_id_str), "%u", rep_id);
        snprintf(x_str, sizeof(x_str), "%f", node.x);
        snprintf(y_str, sizeof(y_str), "%f", node.y);
        snprintf(z_str, sizeof(z_str), "%f", node.z);

        // Create value set for this node
        std::string valueSet = "('" + std::string(cell_line) + "', '" + std::string(job_prefix) + "', " +
                               std::string(rep_id_str) + ", " +
                               std::string(x_str) + ", " + std::string(y_str) + ", " +
                               std::string(z_str) + ", " + std::to_string(start) + ", " +
                               std::to_string(end) + ", '" + insertTime + "')";
        valueSets.push_back(valueSet);
    }

    // Join all value sets with commas
    insertQuery += join(valueSets, ", ");

    // Execute the batch insert
    PGresult *res = PQexec(conn, insertQuery.c_str());
    if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
        std::cerr << "Batch insert failed: " << PQerrorMessage(conn) << std::endl;
    }
    else
    {
        std::cout << "Inserted "
                  << cell_line << ": "
                  << job_prefix << "." << start << "-" << end
                  << " (" << valueSets.size() << " samples) successfully."
                  << std::endl;
    }
    PQclear(res);

    // Close the connection for this thread
    PQfinish(conn);
}

void dumpEnsemble(my_ensemble &chains, const char *out_folder, const char *job_prefix)
{
    auto n_samples = chains.size();
    for (unsigned i = 0; i != n_samples; ++i)
    {
        char *out_file = getOutfile(out_folder, i, job_prefix);
        FILE *output = fopen(out_file, "w");
        if (output == NULL)
        {
            fprintf(stderr, "The output file can not be opened!\n");
            exit(1);
        }
        for (Node node : chains[i])
        {
            fprintf(output, "%f\t%f\t%f\n", node.x, node.y, node.z);
        }
        fclose(output);
        // cout<<"Writing sample "<<i<<" ..."<<endl;
    }
}
