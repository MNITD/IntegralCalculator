#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <thread>
#include <math.h>
#include <mutex>
#include <atomic>
#include <map>
#include <cstring>


class Vertex {
public:
    double x{0};
    double y{0};

    Vertex(double x, double y) {
        this->x = x;
        this->y = y;
    }
};

class Figure {
private:
    std::vector<Vertex> vertexes;
    size_t max_size = 0;
public:
    Figure(size_t size) {
        this->max_size = size;
    }

    void add_vertex(Vertex &v) {
        if (vertexes.size() < max_size) {
            vertexes.push_back(v);
        }
    }

    Vertex get_vertex(unsigned int index) {
        return vertexes.at(index);
    }
};


double calculate_func(double x1, double x2) {
    const double FUN_CONST = 0.002;
    double temp_sum = 0.0;

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            temp_sum += 1 / (5 * (i + 2) + j + 3 + pow((x1 - 16 * j), 6.0) + pow((x2 - 16 * i), 6.0));
        }
    }
    temp_sum += FUN_CONST;
    temp_sum = pow(temp_sum, -1.0);
    return temp_sum;
}

void thread_handler(double x_start, double x_end, double y_start, double y_end, double step, double &sum) {

    double temp_sum = 0.0;
    double y_interval = y_end - y_start;
    double y_middle = (y_start + y_end) / 2;
    while (x_start < x_end) {
        double x_middle = x_start + step / 2;
        temp_sum += calculate_func(x_middle, y_middle) * y_interval * step;
        x_start += step;
    }
    sum += temp_sum;
}

void read(std::ifstream &file, int amount, std::vector<std::string> *text) {
    int i = 0;
    std::string line;

    while (i < amount and std::getline(file, line)) {
        text->push_back(line);
        i++;
    }
}

void reduce_symbols(std::string &line, const std::string &symbols) {

    for (std::string::iterator iter = line.begin(); iter != line.end(); iter++) {
        if (symbols.find(*iter) != std::string::npos) {
            *iter = ' ';
        }
    }
}

std::vector<std::string> split(std::string &line) {
    // std::operator>> separates by spaces
    std::vector<std::string> words;
    std::string word;
    std::istringstream split(line);

    while (split >> word) {
        words.push_back(word);
    }

    return words;
}

void start_threads(int threads_num, Figure &figure, double inc_x, double &sum) {
    std::vector<std::thread *> threads;
    double side = (figure.get_vertex(1).y - figure.get_vertex(0).y) / threads_num;

    for (int i = 0; i < threads_num; i++) {

        double temp_start_y = figure.get_vertex(0).y + side * i;
        double temp_end_y = temp_start_y + side;

//        threads.push_back(
//                new std::thread(thread_handler, figure.get_vertex(0).x, figure.get_vertex(2).x, temp_start_y, temp_end_y, inc_x, std::ref(sum),
//                                std::ref(sum_mutex)));
        thread_handler(figure.get_vertex(0).x, figure.get_vertex(2).x, temp_start_y, temp_end_y, inc_x, sum);
    }

//    std::for_each(threads.begin(), threads.end(), [](std::thread *thread) {
//        thread->join();
//    });
}


inline std::chrono::high_resolution_clock::time_point get_current_time_fenced() {
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto res_time = std::chrono::high_resolution_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return res_time;
}

template<class D>
inline long long to_us(const D &d) {
    return std::chrono::duration_cast<std::chrono::microseconds>(d).count();
}

int main() {
    const std::string symbols_c = "=\"";
    const std::string config_name("../addition/config.txt");
    const std::string abs_inaccuracy_name("abs_inaccuracy");
    const std::string rel_inaccuracy_name("rel_inaccuracy");
    const std::string x1_name("interval_start_x");
    const std::string y1_name("interval_start_y");
    const std::string x2_name("interval_end_x");
    const std::string y2_name("interval_end_y");
    const std::string result_name("res");
    const std::string threads_name("threads");
    char config[1024];
    double abs_inaccuracy;
    double rel_inaccuracy;

    double sendbuf[1], *resbuf;
    double sum = 0.0;
    double temp_sum = 0.0;
    double inc_x = 0.2;
    double actual_abs_inaccuracy = 0;
    int threads_num = 0;
    std::ifstream file_config(config_name);
    std::ofstream file_result;
    std::vector<std::string> config_text;

    std::map<std::string, std::string> config_dict;
    std::mutex sum_mutex;
    std::chrono::high_resolution_clock::time_point count_start_time;
    std::chrono::high_resolution_clock::time_point count_finish_time;

    Figure figure(4);


    int commsize, rank, len, number, max, interval;
    char procname[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &threads_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &len);

    /*
     * Read config file and create config dictionary
     */
    if (rank == 0) {
        std::string string_to_pass;

        read(file_config, 8, &config_text);
        std::for_each(config_text.begin(), config_text.end(), [&symbols_c, &](std::string &line) {
            reduce_symbols(line, symbols_c);
            string_to_pass+line+" ";
        });

        std::strncpy(config, string_to_pass.c_str(), sizeof(config));
        config[sizeof(config) - 1] = 0;

//        std::for_each(config_text.begin(), config_text.end(), [&config_dict](std::string &line) {
//            std::vector<std::string> words = split(line);
//            config_dict[words.at(0)] = words.at(1);
//        });

        resbuf = (double *) malloc(commsize * 1 * sizeof(double));


    }
    MPI_Bcast(config, sizeof(config),
              MPI_CHAR,
            0 , MPI_COMM_WORLD);

    if(rank != 0){
        std::string received_string(config);
        std::vector<std::string> words = split(received_string);
        for(int i = 0 ; i < words.size(); i+=2){
            config_dict[words.at(i)] = words.at(i+1);
        }
    }

    abs_inaccuracy = std::stod(config_dict[abs_inaccuracy_name]);
    rel_inaccuracy = std::stod(config_dict[rel_inaccuracy_name]);

    double x1 = std::stod(config_dict[x1_name]);
    double y1 = std::stod(config_dict[y1_name]);
    double x2 = std::stod(config_dict[x2_name]);
    double y2 = std::stod(config_dict[y2_name]);
    Vertex v1(x1, y1);
    Vertex v2(x1, y2);
    Vertex v3(x2, y2);
    Vertex v4(x2, y1);
    figure.add_vertex(v1);
    figure.add_vertex(v2);
    figure.add_vertex(v3);
    figure.add_vertex(v4);

    threads_num = std::stoi(config_dict[threads_name]);

    /*
     * Start of calculation 
     */

    count_start_time = get_current_time_fenced();

    start_threads(threads_num, figure, inc_x, sum);

    actual_abs_inaccuracy = abs(sum - temp_sum);
    while (actual_abs_inaccuracy > abs_inaccuracy) {
        sum = temp_sum;
        temp_sum = 0.0;
        inc_x /= 2;

        start_threads(threads_num, figure, inc_x, temp_sum);

        actual_abs_inaccuracy = abs(sum - temp_sum);
    }
//    sum = temp_sum;
    sendbuf[0] = temp_sum;


    MPI_Gather(sendbuf, 1,
               MPI_DOUBLE,
               resbuf, 1,
               MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        sum = 0.0;
        for (int i = 0; i < commsize; i++) {
            printf("Root: %f\n", resbuf[i]);
            sum += resbuf[i];
        }
    }



    count_finish_time = get_current_time_fenced();

    /*
     * Finish of calculation (above)
     */

    file_result.open(config_dict[result_name]);

    file_result << "Volume: " << sum << std::endl;
    file_result << "Abs Inaccuracy: " << actual_abs_inaccuracy << std::endl;
    file_result << "Time: " << to_us(count_finish_time - count_start_time) << std::endl;

    file_result.close();

    std::cout << "Time: " << to_us(count_finish_time - count_start_time) << std::endl;
    std::cout << "Volume: " << sum << std::endl;

    return 0;
}