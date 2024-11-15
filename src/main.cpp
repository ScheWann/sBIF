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
    pa.AddArgType("do", "download", ParsingArgs::MAYBE);
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
        bool download = false;
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
            if ((it->first == "do") || (it->first == "download"))
            {
                std::string item = it->second[0];

                item = item.erase(0, item.find_first_not_of(" \t\n\r"));
                item = item.erase(item.find_last_not_of(" \t\n\r") + 1);

                std::cout << "Trimmed item value: '" << item << "'" << std::endl; // 输出修剪后的值
                download = (item == "true" || item == "1");
            }
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
        std::cout << "Download :" << download << endl;
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
        std::ostringstream zip_data;
        vectord2d inter = readInterFiveCols(inter_file_char, weights, chrom_char, chrmfile_char, start, end, resolution);
        getInterNum(inter, n_samples_per_run, false, 1);
        
        // test
        // const char* conninfo = "host=localhost dbname=test user=siyuanzhao";
        
        const char* conninfo = "host=db port=5432 dbname=chromosome_db user=admin password=chromosome";

        clock_t begin, finish;
        double totaltime;
        begin = clock();

        // create zip file path
        char out_file[MAX_CHAR];
        snprintf(out_file, sizeof(out_file), "%s/%s_%s_%u_%u_output.zip", out_folder_char, cell_line_char, job_prefix_char, start, end);

        int error = 0;
        zip_t *zip_archive = NULL;

#pragma omp parallel num_threads(threads)
        {
// only one thread opens the zip file
#pragma omp single
            {
                zip_archive = zip_open(out_file, ZIP_CREATE | ZIP_EXCL, &error);
                if (zip_archive == NULL)
                {
                    fprintf(stderr, "Error opening zip file for writing: %s\n", out_file);
                    exit(1);
                }
            }

// execute the main loop in parallel
#pragma omp for
            for (int i = 0; i < n_runs; ++i)
            {
                my_ensemble chains = SBIF(inter, weights, n_samples_per_run, n_sphere, diam, diam, ki_dist, max_trials, n_iter);

                if (download && zip_archive != NULL)
                {
                    // generate a txt file for each chain and add them to the zip file
                    for (unsigned j = 0; j != n_samples_per_run; ++j)
                    {
                        dumpSingleChainToZip(zip_archive, chains[j], i * n_samples_per_run + j, job_prefix_char, cell_line_char, start, end);
                    }
                }
                else
                {
                    // no download, insert the data into the database
                    for (unsigned j = 0; j != n_samples_per_run; j++)
                    {
                        insertSampleData(conninfo, chains[j], start, end, i * n_samples_per_run + j, job_prefix_char, cell_line_char);
                    }
                }
            }

// close the zip file
#pragma omp single
            {
                if (zip_close(zip_archive) < 0)
                {
                    fprintf(stderr, "Error closing zip file\n");
                    exit(1);
                }
                if (download) {
                    printf("Successfully created zip file: %s\n", out_file);
                }
            }
        }

        finish = clock();
        totaltime = (double)(finish - begin) / CLOCKS_PER_SEC;
        std::cout << "Total cost " << totaltime << " seconds!" << endl;
    }
    return 0;
}
