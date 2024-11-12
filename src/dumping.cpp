#include "dumping.h"
#include <libpq-fe.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <zlib.h>
#include <zip.h>
#include <filesystem>

char* getOutfile(const char* out_folder, unsigned rep_id, const char* job_prefix)
{
    const char* seperator = "/";
    char rep[10];
    sprintf(rep, "%d", rep_id);
    char* out_file = (char*)malloc(sizeof(char) * MAX_CHAR);
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

void dumpSingleChainToZip(zip_t *zip_archive, my_chain& chain, unsigned rep_id, const char* job_prefix, const char* cell_line, unsigned start, unsigned end)
{
    // generate text for each chain
    char file_name[128];
    snprintf(file_name, sizeof(file_name), "%s.%s.%u.%u.%u.txt", cell_line, job_prefix, start, end, rep_id);

    // write the chain data to a string
    std::string chain_data;
    char header[128];
    snprintf(header, sizeof(header), "cell_line: %s, job_prefix: %s, rep_id: %u, start:%u, end: %u \n", cell_line, job_prefix, rep_id, start, end);
    chain_data.append(header);

    for (Node node : chain)
    {
        char line[128];
        snprintf(line, sizeof(line), "%f\t%f\t%f\n", node.x, node.y, node.z);
        chain_data.append(line);
    }

    // add the chain data to the zip archive
    zip_source_t *source = zip_source_buffer(zip_archive, chain_data.c_str(), chain_data.size(), 0);
    if (source == NULL)
    {
        fprintf(stderr, "Error creating zip source for chain %u\n", rep_id);
        exit(1);
    }

    if (zip_file_add(zip_archive, file_name, source, ZIP_FL_OVERWRITE) < 0)
    {
        fprintf(stderr, "Error adding file to zip archive: %s\n", zip_strerror(zip_archive));
        exit(1);
    }
}

// Euclidean distance calculation function
double calculateDistance(const Node& node1, const Node& node2) {
    return std::sqrt(
        std::pow(node2.x - node1.x, 2) +
        std::pow(node2.y - node1.y, 2) +
        std::pow(node2.z - node1.z, 2)
    );
}

// insert into the database
std::string join(const std::vector<std::string>& elements, const std::string& delimiter) {
    std::ostringstream oss;
    for (size_t i = 0; i < elements.size(); ++i) {
        oss << elements[i];
        if (i < elements.size() - 1) {
            oss << delimiter;
        }
    }
    return oss.str();
}

void insertSampleData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix, const char *cell_line)
{
    // Establish a new connection for this thread
    PGconn* conn = PQconnectdb(conninfo);
    if (PQstatus(conn) != CONNECTION_OK) {
        std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    // Prepare the insert query
    std::string insertQuery = "INSERT INTO position (cell_line, chrID, sampleID, X, Y, Z, start_value, end_value) VALUES ";
    std::vector<std::string> valueSets;

    for (Node node : chain) {
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
                                std::to_string(end) + ")";
        valueSets.push_back(valueSet);
    }

    // Join all value sets with commas
    insertQuery += join(valueSets, ", ");

    // Execute the batch insert
    PGresult* res = PQexec(conn, insertQuery.c_str());
    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "Batch insert failed: " << PQerrorMessage(conn) << std::endl;
    } else {
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

std::vector<std::vector<double>> allDistanceMatrix(const std::vector<Node>& chain) {
    unsigned numNodes = chain.size();
    std::vector<std::vector<double>> distanceMatrix(numNodes, std::vector<double>(numNodes, 0.0));

    for (unsigned i = 0; i < numNodes; ++i) {
        for (unsigned j = i + 1; j < numNodes; ++j) {
            double distance = calculateDistance(chain[i], chain[j]);
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
            cout<<"Writing distance matrix "<< distance <<" ..."<<endl;
        }
    }
    return distanceMatrix;
}

void dumpEnsemble(my_ensemble& chains, const char* out_folder, const char* job_prefix)
{
    auto n_samples = chains.size();
    for (unsigned i = 0; i != n_samples; ++i)
    {
        char* out_file = getOutfile(out_folder, i, job_prefix);
        FILE* output = fopen(out_file, "w");
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
        //cout<<"Writing sample "<<i<<" ..."<<endl;
    }
}
