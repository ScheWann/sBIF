#include "dumping.h"
#include <libpq-fe.h>
#include <sstream>

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