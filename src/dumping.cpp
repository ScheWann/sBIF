#include "dumping.h"
#include <libpq-fe.h>

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
void insertSampleData(PGconn* conn, my_chain& chain, unsigned start, unsigned end, unsigned rep_id, const char* job_prefix)
{
    if (conn == NULL)
    {
        std::cerr << "Connection to database failed!" << std::endl;
        return;
    }

    const char* insertQuery = "INSERT INTO position (sampleID, chrID, X, Y, Z, start_value, end_value) VALUES ($1, $2, $3, $4, $5, $6, $7);";

    for (Node node : chain)
    {
        const char* params[7];

        char rep_id_str[12];
        snprintf(rep_id_str, sizeof(rep_id_str), "%u", rep_id);

        char x_str[32], y_str[32], z_str[32];
        snprintf(x_str, sizeof(x_str), "%f", node.x);
        snprintf(y_str, sizeof(y_str), "%f", node.y);
        snprintf(z_str, sizeof(z_str), "%f", node.z);

        params[0] = rep_id_str;
        params[1] = job_prefix;
        params[2] = x_str;
        params[3] = y_str;
        params[4] = z_str;

        char start_str[12], end_str[12];
        snprintf(start_str, sizeof(start_str), "%u", start);
        snprintf(end_str, sizeof(end_str), "%u", end);
        
        params[5] = start_str;
        params[6] = end_str;

        PGresult* res = PQexecParams(conn, insertQuery, 7, NULL, params, NULL, NULL, 0);
        if (PQresultStatus(res) != PGRES_COMMAND_OK)
        {
            fprintf(stderr, "Insert failed: %s\n", PQerrorMessage(conn));
            std::cerr << "Insert failed for node (" << node.x << ", " << node.y << ", " << node.z 
                      << "): " << PQerrorMessage(conn) << std::endl;
        }
        else
        {
            std::cout << "Inserted sample with rep_id: " << rep_id 
                      << ", job_prefix: " << job_prefix 
                      << ", coordinates: (" << node.x << ", " << node.y << ", " << node.z << ")" << std::endl;
        }

        PQclear(res);
    }
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