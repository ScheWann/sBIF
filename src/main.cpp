#include <iostream>
#include <time.h>
#include "optimization.h"
#include <condition_variable>
#include <vector>
#include <ctime>
#include <mutex>
#include "dumping.h"
#include "parsingargs.h"
#include <string.h>
#include "help.h"
#include <sstream>
#include <thread>
#include <atomic>
#include <libpq-fe.h>

std::atomic<unsigned> global_sample_index(0);
int main(int argc, char *argv[])
{

    string tmpPara = "";
    for (int i = 1; i < argc; i++)
    {
        // cout << i << "=" << argv[i] <<"---"<< endl;
        if (strlen(argv[i]) == 0)
        {
            std::cout << "find NULL" << endl;
            tmpPara += char(31);
        }
        else
        {
            tmpPara += argv[i];
        }
        tmpPara += " ";
    }

    map<string, vector<string>> result;
    ParsingArgs pa;
    pa.AddArgType("i", "inter", ParsingArgs::MUST);
    pa.AddArgType("c", "chrom", ParsingArgs::MUST);
    pa.AddArgType("l", "chrlens", ParsingArgs::MUST);
    pa.AddArgType("s", "start", ParsingArgs::MUST);
    pa.AddArgType("e", "end", ParsingArgs::MUST);
    pa.AddArgType("cl", "cell_line", ParsingArgs::MUST);
    pa.AddArgType("o", "out", ParsingArgs::MAYBE);
    pa.AddArgType("r", "res", ParsingArgs::MAYBE);
    // pa.AddArgType("ex", "example", ParsingArgs::MAYBE);
    pa.AddArgType("d", "fibdens", ParsingArgs::MAYBE);
    pa.AddArgType("ns", "nsamp", ParsingArgs::MAYBE);
    pa.AddArgType("nr", "nruns", ParsingArgs::MAYBE);
    pa.AddArgType("n", "nsphere", ParsingArgs::MAYBE);
    pa.AddArgType("k", "kidist", ParsingArgs::MAYBE);
    pa.AddArgType("m", "maxtrial", ParsingArgs::MAYBE);
    pa.AddArgType("ni", "niter", ParsingArgs::MAYBE);
    pa.AddArgType("j", "jobpre", ParsingArgs::MAYBE);
    pa.AddArgType("p", "threads", ParsingArgs::MAYBE);
    pa.AddArgType("h", "help", ParsingArgs::NO);

    result.clear();
    // cout << "Input is:" << tmpPara << endl;
    std::string errPos;
    int iRet = pa.Parse(tmpPara, result, errPos);
    if ((0 > iRet) || (tmpPara.empty()))
    {
        std::cout << "Error: wrong options with flag " << iRet << endl;
        printHelp();
        return 0;
    }
    else
    {
        string inter_file;
        string chrom;
        string chrmfile;
        string cell_line;
        unsigned start;
        unsigned end;
        string out_folder;
        unsigned resolution = 2000;
        double fiber_density = 0.2368;
        unsigned n_samples = 50000;
        // unsigned n_examples = 5;
        unsigned n_samples_per_run = 100;
        unsigned n_sphere = 50;
        unsigned ki_dist = 80;
        unsigned max_trials = 100;
        unsigned n_iter = 100;
        unsigned threads = 1;
        string job_prefix = "test";

        if ((result.find("h") != result.end()) && (result.find("help") != result.end()))
        {
            printHelp();
            return 0;
        }
        if ((result.find("i") == result.end()) && (result.find("inter") == result.end()))
        {
            cout << "Error: missing required parameters..." << endl;
            printHelp();
            return 0;
        }
        if ((result.find("cl") == result.end()) && (result.find("cell_line") == result.end()))
        {
            cout << "Error: missing required parameters..." << endl;
            printHelp();
            return 0;
        }
        if ((result.find("c") == result.end()) && (result.find("chrom") == result.end()))
        {
            cout << "Error: missing required parameters..." << endl;
            printHelp();
            return 0;
        }
        if ((result.find("l") == result.end()) && (result.find("chrlens") == result.end()))
        {
            cout << "Error: missing required parameters..." << endl;
            printHelp();
            return 0;
        }
        if ((result.find("s") == result.end()) && (result.find("start") == result.end()))
        {
            cout << "Error: missing required parameters..." << endl;
            printHelp();
            return 0;
        }
        if ((result.find("e") == result.end()) && (result.find("end") == result.end()))
        {
            cout << "Error: missing required parameters..." << endl;
            printHelp();
            return 0;
        }
        map<std::string, std::vector<std::string>>::iterator it = result.begin();
        for (; it != result.end(); ++it)
        {
            if ((it->first.compare("i") == 0) || (it->first.compare("inter") == 0))
                inter_file = it->second[0];
            if ((it->first.compare("c") == 0) || (it->first.compare("chrom") == 0))
                chrom = it->second[0];
            if ((it->first.compare("l") == 0) || (it->first.compare("chrlens") == 0))
                chrmfile = it->second[0];
            if ((it->first.compare("s") == 0) || (it->first.compare("start") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> start;
            }
            if ((it->first.compare("e") == 0) || (it->first.compare("end") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> end;
            }
            if ((it->first.compare("cl") == 0) || (it->first.compare("cell_line") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> cell_line;
            }
            // if ((it->first.compare("ex") == 0) || (it->first.compare("example") == 0))
            // {
            //     std::stringstream item;
            //     item << it->second[0];
            //     item >> n_examples;
            // }
            if ((it->first.compare("o") == 0) || (it->first.compare("out") == 0))
                out_folder = it->second[0];
            if ((it->first.compare("r") == 0) || (it->first.compare("res") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> resolution;
            }
            if ((it->first.compare("d") == 0) || (it->first.compare("fibdens") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> fiber_density;
            }
            if ((it->first.compare("ns") == 0) || (it->first.compare("nsamp") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> n_samples;
            }
            if ((it->first.compare("nr") == 0) || (it->first.compare("nruns") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> n_samples_per_run;
            }
            if ((it->first.compare("n") == 0) || (it->first.compare("nsphere") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> n_sphere;
            }
            if ((it->first.compare("k") == 0) || (it->first.compare("kidist") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> ki_dist;
            }
            if ((it->first.compare("m") == 0) || (it->first.compare("maxtrial") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> max_trials;
            }
            if ((it->first.compare("ni") == 0) || (it->first.compare("niter") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> n_iter;
            }
            if ((it->first.compare("p") == 0) || (it->first.compare("threads") == 0))
            {
                std::stringstream item;
                item << it->second[0];
                item >> threads;
            }
            if ((it->first.compare("j") == 0) || (it->first.compare("jobpre") == 0))
                job_prefix = it->second[0];
            if ((it->first.compare("h")) == 0 || (it->first.compare("help") == 0))
            {
                printHelp();
            }
        }
        double diam = getDiam(resolution, fiber_density);
        unsigned n_runs = n_samples / n_samples_per_run;
        // string command;
        // command = "mkdir -p " + out_folder;
        // std::system(command.c_str());
        std::cout << "Parameters: " << endl;
        std::cout << "Interaction file :" << inter_file << endl;
        std::cout << "Chromosome :" << chrom << endl;
        std::cout << "Chrom lengths file :" << chrmfile << endl;
        std::cout << "Start position:" << start << endl;
        std::cout << "End position :" << end << endl;
        std::cout << "Cell line :" << cell_line << endl;
        std::cout << "Output folder :" << out_folder << endl;
        // std::cout << "Examples showing :" << n_examples << endl;
        std::cout << "Resolution :" << resolution << endl;
        std::cout << "Fiber density :" << fiber_density << endl;
        std::cout << "Number of samples :" << n_samples << endl;
        std::cout << "Number of samples per run :" << n_samples_per_run << endl;
        std::cout << "Number of potential sphere points :" << n_sphere << endl;
        std::cout << "Knock-in distance :" << ki_dist << endl;
        std::cout << "Maximum trials :" << max_trials << endl;
        std::cout << "Number of iteractions :" << n_iter << endl;
        std::cout << "Job prefix :" << job_prefix << endl;
        std::cout << "Number of threads :" << threads << endl;
        std::cout << "Bead diameter: " << diam << endl;
        std::cout << "Generating samples ..." << endl;
        unsigned region_size = end - start;
        unsigned n_nodes = (region_size % resolution == 0) ? (region_size / resolution) : (region_size / resolution + 1);
        vectord2d weights(n_nodes, vectord(n_nodes));
        const char *inter_file_char = inter_file.c_str();
        const char *chrmfile_char = chrmfile.c_str();
        const char *chrom_char = chrom.c_str();
        const char *cell_line_char = cell_line.c_str();
        const char *out_folder_char = out_folder.c_str();
        const char *job_prefix_char = job_prefix.c_str();

        // Instead of storing all distances, use online statistics
        unsigned expected_dist_size = 0;
        std::vector<double> sum_distances;
        std::vector<unsigned> count_condensed;
        std::vector<float> best_distances;
        unsigned best_sample_id = 0;
        double best_corr = 0.0;  // Use correlation instead of score
        std::mutex stats_mutex;

        vectord2d inter = readInterFiveCols(inter_file_char, weights, chrom_char, chrmfile_char, start, end, resolution);
        getInterNum(inter, n_samples_per_run, false, 1);

        // test
        // const char *conninfo = "host=localhost dbname=test user=siyuanzhao";

        const char *conninfo = "host=db port=5432 dbname=chromosome_db user=admin password=chromosome";

        clock_t begin, finish;
        double totaltime;
        begin = clock();


        auto t_position_start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel num_threads(threads)
        {
            // Create a thread-local database connection
            PGconn *thread_conn = PQconnectdb(conninfo);
            if (PQstatus(thread_conn) != CONNECTION_OK) {
                std::cerr << "Thread connection to database failed: " << PQerrorMessage(thread_conn) << std::endl;
                PQfinish(thread_conn);
                thread_conn = nullptr;
            }
            
            std::vector<float> local_dist;
            char local_job_prefix[64];
            char local_cell_line[64];
            strncpy(local_job_prefix, job_prefix_char, sizeof(local_job_prefix));
            strncpy(local_cell_line, cell_line_char, sizeof(local_cell_line));

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < n_runs; ++i)
            {
                if (thread_conn == nullptr) {
                    std::cerr << "Skipping run " << i << " due to database connection failure" << std::endl;
                    continue;
                }
                
                my_ensemble chains = SBIF(inter, weights, n_samples_per_run, n_sphere, diam, diam, ki_dist, max_trials, n_iter);
                    for (unsigned j = 0; j != n_samples_per_run; j++)
                    {
                        unsigned sample_id = i * n_samples_per_run + j;
                        insertSampleData(thread_conn, chains[j], start, end, sample_id, local_job_prefix, local_cell_line);
                        
                        // Compute distance vector
                        std::vector<float> local_dist = computeDistanceVector(chains[j]);
                        
                        // Insert distance data immediately to avoid memory buildup
                        insertDistanceDataFromVector(thread_conn, local_cell_line, local_job_prefix, sample_id, start, end, local_dist);
                        
                        // Update statistics online (thread-safe)
                        {
                            std::lock_guard<std::mutex> lock(stats_mutex);
                            if (sum_distances.empty()) {
                                expected_dist_size = local_dist.size();
                                sum_distances.resize(expected_dist_size, 0.0);
                                count_condensed.resize(expected_dist_size, 0);
                            }
                            
                            // Update running sum for average
                            for (size_t k = 0; k < local_dist.size(); ++k) {
                                sum_distances[k] += local_dist[k];
                                if (local_dist[k] <= 80.0f) {
                                    count_condensed[k]++;
                                }
                            }
                            
                            // Check if this is the best sample using Pearson correlation with current average
                            if (global_sample_index.load() > 10) { // Only start checking after some samples
                                unsigned current_count = global_sample_index.load();
                                
                                // Calculate current average vector
                                std::vector<double> current_avg(local_dist.size());
                                for (size_t k = 0; k < local_dist.size(); ++k) {
                                    current_avg[k] = sum_distances[k] / current_count;
                                }
                                
                                // Calculate means for centering
                                double avg_mean = std::accumulate(current_avg.begin(), current_avg.end(), 0.0) / current_avg.size();
                                double local_mean = std::accumulate(local_dist.begin(), local_dist.end(), 0.0) / local_dist.size();
                                
                                // Calculate Pearson correlation coefficient
                                double numer = 0.0, norm_avg_sq = 0.0, norm_local_sq = 0.0;
                                for (size_t k = 0; k < local_dist.size(); ++k) {
                                    double avg_centered = current_avg[k] - avg_mean;
                                    double local_centered = local_dist[k] - local_mean;
                                    numer += avg_centered * local_centered;
                                    norm_avg_sq += avg_centered * avg_centered;
                                    norm_local_sq += local_centered * local_centered;
                                }
                                
                                double norm_avg = std::sqrt(norm_avg_sq);
                                double norm_local = std::sqrt(norm_local_sq);
                                double corr = (norm_avg * norm_local == 0.0) ? 0.0 : numer / (norm_avg * norm_local);
                                
                                // Keep the sample with highest absolute correlation
                                if (std::abs(corr) > std::abs(best_corr)) {
                                    best_corr = corr;
                                    best_distances = local_dist;
                                    best_sample_id = sample_id;
                                }
                            }
                        }
                        
                        global_sample_index.fetch_add(1);
                    }
            }
            
            // Clean up thread-local database connection
            if (thread_conn != nullptr) {
                PQfinish(thread_conn);
            }
        }

        auto t_position_end = std::chrono::high_resolution_clock::now();
        double dur_position = std::chrono::duration<double>(t_position_end - t_position_start).count();
        std::cout << "[position and distance data inserted DONE] All samples' data inserted in " << dur_position << " seconds." << std::endl;

        // Compute final statistics from accumulated data
        auto t_calc_start = std::chrono::high_resolution_clock::now();
        
        // Calculate average vector
        std::vector<float> avg_vector(expected_dist_size);
        for (size_t i = 0; i < expected_dist_size; ++i) {
            avg_vector[i] = static_cast<float>(sum_distances[i] / n_samples);
        }
        
        // Calculate frequency condensed vector
        std::vector<float> freq_condensed(expected_dist_size);
        for (size_t i = 0; i < expected_dist_size; ++i) {
            freq_condensed[i] = static_cast<float>(count_condensed[i]) / n_samples;
        }
        
        // Convert to full matrix
        std::vector<float> freq_full = squareformFullMatrix(freq_condensed);
        
        auto t_calc_end = std::chrono::high_resolution_clock::now();
        double dur_calc = std::chrono::duration<double>(t_calc_end - t_calc_start).count();
        std::cout << "[statistics computed DONE] Computed all statistics in " << dur_calc << " seconds." << std::endl;

        // Insert average vector and frequency data into the database
        auto t_insert_calc_start = std::chrono::high_resolution_clock::now();
        insertCalcDistance(conninfo, cell_line_char, job_prefix_char, start, end, avg_vector, freq_full, best_distances, best_sample_id);
        auto t_insert_calc_end = std::chrono::high_resolution_clock::now();
        double dur_insert_calc = std::chrono::duration<double>(t_insert_calc_end - t_insert_calc_start).count();
        std::cout << "[calculated data inserted DONE] Inserted calculated data into the database in " << dur_insert_calc << " seconds." << std::endl;

        // Memory cleanup is no longer needed as we don't store all distances
        // Clear statistics vectors if needed
        sum_distances.clear();
        count_condensed.clear();
        avg_vector.clear();
        freq_condensed.clear();
        freq_full.clear();
        best_distances.clear();

        finish = clock();
        totaltime = (double)(finish - begin) / CLOCKS_PER_SEC;
        std::cout << "Total cost " << totaltime << " seconds!" << endl;
    }
    return 0;
}
