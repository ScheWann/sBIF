#include "dumping.h"
#include <libpq-fe.h>
#include <sstream>
#include <zip.h> 

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

void dumpSingleChain(my_chain& chain, const char* out_folder, unsigned rep_id, const char* job_prefix)
{
    char* out_file = getOutfile(out_folder, rep_id, job_prefix);
    FILE* output = fopen(out_file, "w");
    if (output == NULL)
    {
        fprintf(stderr, "The output file can not be opened!\n");
        exit(1);
    }
    fprintf(output, "rep_id: %u, job_prefix: %s\n", rep_id, job_prefix);
    for (Node node : chain)
    {
        fprintf(output, "%f\t%f\t%f\n", node.x, node.y, node.z);
    }
    fclose(output);
    //cout<<"Writing sample "<<rep_id<<" ..."<<endl;
}

std::vector<char> createZipInMemory(const std::vector<my_chain>& chains, const std::string& job_prefix) {
    zip_source_t *source = nullptr;
    zip_t *archive = nullptr;
    zip_error_t error;

    // Initialize zip_error_t
    zip_error_init(&error);

    // Create the ZIP source buffer
    source = zip_source_buffer_create(nullptr, 0, 0, &error);
    if (!source) {
        std::string errorMsg = "Failed to create ZIP source buffer: " + std::string(zip_error_strerror(&error));
        zip_error_fini(&error);
        throw std::runtime_error(errorMsg);
    }

    // Open the archive from the source
    archive = zip_open_from_source(source, ZIP_CREATE, &error);
    if (!archive) {
        std::string errorMsg = "Failed to open ZIP archive: " + std::string(zip_error_strerror(&error));
        zip_source_free(source);
        zip_error_fini(&error);
        throw std::runtime_error(errorMsg);
    }

    // Add files to the archive
    for (size_t i = 0; i < chains.size(); ++i) {
        const my_chain& chain = chains[i];
        
        // Create the content for each chain
        std::ostringstream oss;
        oss << "rep_id: " << i << ", job_prefix: " << job_prefix << "\n";
        for (const Node& node : chain) {
            oss << node.x << "\t" << node.y << "\t" << node.z << "\n";
        }
        std::string data = oss.str();

        // Create a filename
        std::string filename = job_prefix + "_chain_" + std::to_string(i) + ".txt";

        // Add file to archive
        zip_source_t* data_source = zip_source_buffer_create(data.c_str(), data.size(), 0, &error);
        if (!data_source || zip_file_add(archive, filename.c_str(), data_source, ZIP_FL_OVERWRITE) == -1) {
            zip_source_free(data_source);  // Free data source if adding failed
            std::string errorMsg = "Failed to add file to ZIP archive: " + std::string(zip_error_strerror(&error));
            zip_discard(archive);  // Discard archive if we fail
            zip_error_fini(&error);
            throw std::runtime_error(errorMsg);
        }
    }

    // Close the archive
    if (zip_close(archive) == -1) {
        zip_source_free(source);
        zip_error_fini(&error);
        throw std::runtime_error("Failed to close ZIP archive");
    }

    // Get the data from the ZIP source buffer
    zip_stat_t zip_stat;
    zip_stat_init(&zip_stat);
    zip_stat_index(archive, 0, 0, &zip_stat);

    // Read the data from the ZIP source buffer
    std::vector<char> zip_data(zip_stat.size);
    zip_source_open(source);
    zip_source_read(source, zip_data.data(), zip_stat.size);
    zip_source_close(source);

    // Clean up the error object
    zip_error_fini(&error);

    return zip_data;
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

void insertSampleData(const char *conninfo, my_chain &chain, unsigned start, unsigned end, unsigned rep_id, const char *job_prefix)
{
    // Establish a new connection for this thread
    PGconn* conn = PQconnectdb(conninfo);
    if (PQstatus(conn) != CONNECTION_OK) {
        std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
        PQfinish(conn);
        return;
    }

    // Prepare the insert query
    std::string insertQuery = "INSERT INTO position (sampleID, chrID, X, Y, Z, start_value, end_value) VALUES ";
    std::vector<std::string> valueSets;

    for (Node node : chain) {
        char rep_id_str[12], x_str[32], y_str[32], z_str[32];
        snprintf(rep_id_str, sizeof(rep_id_str), "%u", rep_id);
        snprintf(x_str, sizeof(x_str), "%f", node.x);
        snprintf(y_str, sizeof(y_str), "%f", node.y);
        snprintf(z_str, sizeof(z_str), "%f", node.z);

        // Create value set for this node
        std::string valueSet = "(" + std::string(rep_id_str) + ", '" + std::string(job_prefix) + "', " +
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