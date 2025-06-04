#include "dumping.h"
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
double calculateDistance(const Node node1, const Node node2)
{
    return std::sqrt(
        std::pow(node2.x - node1.x, 2) +
        std::pow(node2.y - node1.y, 2) +
        std::pow(node2.z - node1.z, 2));
}

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

// intert bead pair distance of each sample to the database
// void insertDistanceData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix, const char *cell_line)
// {
//     PGconn *conn = PQconnectdb(conninfo);
//     if (PQstatus(conn) != CONNECTION_OK)
//     {
//         std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
//         PQfinish(conn);
//         return;
//     }

//     std::string copyQuery = "COPY distance (cell_line, chrid, sampleid, start_value, end_value, n_beads, distance_vector, insert_time) FROM STDIN WITH (FORMAT csv)";
//     PGresult *res = PQexec(conn, copyQuery.c_str());
//     if (PQresultStatus(res) != PGRES_COPY_IN)
//     {
//         std::cerr << "COPY command failed: " << PQerrorMessage(conn) << std::endl;
//         PQclear(res);
//         PQfinish(conn);
//         return;
//     }
//     PQclear(res);

//     std::stringstream copyData;

//     // time stamp
//     char time_buffer[20];
//     auto now = std::chrono::system_clock::now();
//     std::time_t now_c = std::chrono::system_clock::to_time_t(now);
//     std::tm *now_tm = std::localtime(&now_c);
//     std::strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", now_tm);

//     int n_beads = chain.size();

//     // top triangle distance: vector<double>
//     std::vector<double> distances;
//     distances.reserve(n_beads * (n_beads - 1) / 2);
//     for (int i = 0; i < n_beads; ++i)
//     {
//         for (int j = i + 1; j < n_beads; ++j)
//         {
//             double d = calculateDistance(chain[i], chain[j]);
//             distances.push_back(d);
//         }
//     }

//     // like: {1.0,2.0,3.0,...}
//     std::stringstream arrStream;
//     arrStream << "{";
//     for (size_t i = 0; i < distances.size(); ++i)
//     {
//         arrStream << distances[i];
//         if (i != distances.size() - 1)
//             arrStream << ",";
//     }
//     arrStream << "}";
//     std::string arrStr = arrStream.str();

//     copyData << cell_line << ","              // cell_line
//              << job_prefix << ","             // chrID
//              << rep_id << ","                 // sampleID
//              << start << ","                  // start_value
//              << end << ","                    // end_value
//              << n_beads << ","                // n_beads
//              << "\"" << arrStr << "\","       // distance_vector
//              << time_buffer                   // insert_time
//              << "\n";

//     std::string copyDataStr = copyData.str();
//     if (PQputCopyData(conn, copyDataStr.c_str(), copyDataStr.size()) != 1)
//     {
//         std::cerr << "Failed to send data to PostgreSQL: " << PQerrorMessage(conn) << std::endl;
//     }

//     if (PQputCopyEnd(conn, NULL) != 1)
//     {
//         std::cerr << "Failed to complete COPY command: " << PQerrorMessage(conn) << std::endl;
//     }

//     PGresult *endRes = PQgetResult(conn);
//     if (PQresultStatus(endRes) != PGRES_COMMAND_OK)
//     {
//         std::cerr << "COPY failed: " << PQerrorMessage(conn) << std::endl;
//     }

//     PQclear(endRes);
//     PQfinish(conn);
// }

void insertDistanceData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix, const char *cell_line)
{
    PGconn *conn = PQconnectdb(conninfo);
    if (PQstatus(conn) != CONNECTION_OK) {
        std::cerr << "Connection failed: " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    int n_beads = static_cast<int>(chain.size());
    size_t num_vals = size_t(n_beads) * (n_beads - 1) / 2;
    std::vector<float> distances_f;
    distances_f.reserve(num_vals);
    for (int i = 0; i < n_beads; ++i) {
        for (int j = i + 1; j < n_beads; ++j) {
            double d = calculateDistance(chain[i], chain[j]);
            distances_f.push_back(static_cast<float>(d));
        }
    }
    size_t bin_len = distances_f.size() * sizeof(float);
    std::string binBuf;
    binBuf.resize(bin_len);
    std::memcpy(&binBuf[0], distances_f.data(), bin_len);

    std::string str_sampleid    = std::to_string(rep_id);
    std::string str_start_value = std::to_string(start);
    std::string str_end_value   = std::to_string(end);
    std::string str_n_beads     = std::to_string(n_beads);

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
    int         paramLengths[nParams];
    int         paramFormats[nParams];

    // (1) cell_line
    paramValues[0]  = cell_line;
    paramLengths[0] = std::strlen(cell_line);
    paramFormats[0] = 0;

    // (2) chrid
    paramValues[1]  = job_prefix;
    paramLengths[1] = std::strlen(job_prefix);
    paramFormats[1] = 0;

    // (3) sampleid
    paramValues[2]  = str_sampleid.c_str();
    paramLengths[2] = static_cast<int>(str_sampleid.size());
    paramFormats[2] = 0;

    // (4) start_value
    paramValues[3]  = str_start_value.c_str();
    paramLengths[3] = static_cast<int>(str_start_value.size());
    paramFormats[3] = 0;

    // (5) end_value
    paramValues[4]  = str_end_value.c_str();
    paramLengths[4] = static_cast<int>(str_end_value.size());
    paramFormats[4] = 0;

    // (6) n_beads
    paramValues[5]  = str_n_beads.c_str();
    paramLengths[5] = static_cast<int>(str_n_beads.size());
    paramFormats[5] = 0;

    // (7) distance_vector
    paramValues[6]  = binBuf.data();
    paramLengths[6] = static_cast<int>(bin_len);
    paramFormats[6] = 1;

    // (8) insert_time
    paramValues[7]  = str_time.c_str();
    paramLengths[7] = static_cast<int>(str_time.size());
    paramFormats[7] = 0;

    // 6. 执行 INSERT
    PGresult *res = PQexecParams(conn,
        "INSERT INTO distance("
          "cell_line, chrid, sampleid, start_value, end_value, n_beads, distance_vector, insert_time"
        ") VALUES("
          "$1,        $2,    $3,        $4,           $5,       $6,      $7,             $8"
        ");",
        nParams,
        nullptr,
        paramValues,
        paramLengths,
        paramFormats,
        0
    );

    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "INSERT 失败: " << PQerrorMessage(conn) << std::endl;
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
