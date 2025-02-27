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

void generateDistanceFile(my_chain &chain, int chain_index, zip_t *zip_archive, const char *job_prefix_char, const char *cell_line_char, unsigned start, unsigned end)
{
    char distance_filename[MAX_CHAR];
    snprintf(distance_filename, sizeof(distance_filename), "%s_%s_%u_%u_distance_%d.txt", job_prefix_char, cell_line_char, start, end, chain_index);

    // store the content of the distance file
    std::vector<char> distance_content;
    // reserve memory for the content to avoid reallocation
    distance_content.reserve(10000);

    int n_beads = chain.size();
    for (int i = 0; i < n_beads; ++i)
    {
        for (int j = i + 1; j < n_beads; ++j)
        {
            double dist = calculateDistance(chain[i], chain[j]);

            // Format the distance string and append it to the content vector
            char buffer[256];
            int written = snprintf(buffer, sizeof(buffer), "Distance between bead %d and bead %d: %.3f\n", i, j, dist);
            distance_content.insert(distance_content.end(), buffer, buffer + written);
        }
    }

    // store the content in a zip source
    zip_source_t *distance_source = zip_source_buffer(zip_archive, distance_content.data(), distance_content.size(), 0);
    if (distance_source == NULL)
    {
        fprintf(stderr, "Error creating zip source for distance file\n");
        exit(1);
    }

    if (zip_file_add(zip_archive, distance_filename, distance_source, ZIP_FL_OVERWRITE) < 0)
    {
        fprintf(stderr, "Error adding distance file to zip\n");
        exit(1);
    }
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
void insertDistanceData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix, const char *cell_line)
{
    // database connection
    PGconn *conn = PQconnectdb(conninfo);
    if (PQstatus(conn) != CONNECTION_OK)
    {
        std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    std::string insertQuery = "INSERT INTO distance (cell_line, chrID, sampleID, start_value, end_value, bead_i, bead_j, distance, insert_time) VALUES ";
    std::vector<std::string> valueSets;

    // get the current local time
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm *now_tm = std::localtime(&now_c);
    char time_buffer[20];
    std::strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", now_tm);
    std::string insertTime(time_buffer);

    // all bead pairs
    int n_beads = chain.size();
    for (int i = 0; i < n_beads; ++i)
    {
        for (int j = i + 1; j < n_beads; ++j)
        {
            // distance calculation
            double dist = calculateDistance(chain[i], chain[j]);

            char dist_str[32], i_str[12], j_str[12];
            snprintf(dist_str, sizeof(dist_str), "%.3f", dist);
            snprintf(i_str, sizeof(i_str), "%d", i);
            snprintf(j_str, sizeof(j_str), "%d", j);

            std::string valueSet = "('" + std::string(cell_line) + "', '" +
                                   std::string(job_prefix) + "', " +
                                   std::to_string(rep_id) + ", " +
                                   std::to_string(start) + ", " +
                                   std::to_string(end) + ", " +
                                   i_str + ", " + j_str + ", " +
                                   dist_str + ", '" + insertTime + "')";
            valueSets.push_back(valueSet);
        }
    }

    if (!valueSets.empty())
    {
        std::cout << "Inserting distances for "
                  << cell_line << ": "
                  << job_prefix << "." << start << "-" << end
                  << " (" << valueSets.size() << " pairs)..."
                  << std::endl;
        insertQuery += join(valueSets, ", ");

        PGresult *res = PQexec(conn, insertQuery.c_str());
        if (PQresultStatus(res) != PGRES_COMMAND_OK)
        {
            std::cerr << "Distance insert failed: " << PQerrorMessage(conn) << std::endl;
        }
        else
        {
            std::cout << "Inserted distances for "
                      << cell_line << ": "
                      << job_prefix << "." << start << "-" << end
                      << " (" << valueSets.size() << " pairs) successfully."
                      << std::endl;
        }
        PQclear(res);
    }

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

std::vector<std::vector<double>> allDistanceMatrix(const std::vector<Node> &chain)
{
    unsigned numNodes = chain.size();
    std::vector<std::vector<double>> distanceMatrix(numNodes, std::vector<double>(numNodes, 0.0));

    for (unsigned i = 0; i < numNodes; ++i)
    {
        for (unsigned j = i + 1; j < numNodes; ++j)
        {
            double distance = calculateDistance(chain[i], chain[j]);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
            cout << "Writing distance matrix " << distance << " ..." << endl;
        }
    }
    return distanceMatrix;
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
