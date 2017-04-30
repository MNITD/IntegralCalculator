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

void thread_handler(double x_start, double x_end, double y_start, double y_end,  double step, double & sum, std::mutex & sum_mutex);
void read(std::ifstream  & file, int amount, std::vector<std::string>*text);
void reduce_symbols(std::string & line, const std::string & symbols);
std::vector<std::string>  split(std::string & line);
double calculate_func(double x1, double x2);
void update_sum(double & sum, double & temp_sum, std::mutex & sum_mutex);

void thread_handler(double x_start, double x_end, double y_start, double y_end, double step, double & sum, std::mutex & sum_mutex){

    double temp_sum = 0.0;
    double y_interval = y_end - y_start;
    double y_middle = (y_start + y_end) / 2;
    while(x_start < x_end){
        double x_middle = x_start + step / 2;
        temp_sum += calculate_func(x_middle, y_middle) * y_interval * step;
        x_start += step;
    }
    update_sum(sum, temp_sum, sum_mutex);
}
void read(std::ifstream  & file, int amount, std::vector<std::string>*text){

    int i = 0;
    std::string line;

    while ( i < amount and std::getline(file, line) ){
        text->push_back(line);
        i++;
    }
}
void reduce_symbols(std::string & line,const std::string & symbols){

    for(std::string::iterator iter = line.begin(); iter != line.end(); iter++) {
        if ( symbols.find( *iter ) != std::string::npos ) {
            *iter = ' ';
        }
    }
}
std::vector<std::string> split(std::string & line){
    // std::operator>> separates by spaces
    std::vector<std::string> words;
    std::string word;
    std::istringstream split(line);

    while( split >> word ) {
        words.push_back(word);
        std::cout << word << std::endl;
    }

    return words;
}
void update_sum(double & sum, double & temp_sum, std::mutex & sum_mutex){
    std::lock_guard<std::mutex> lock(sum_mutex);
    sum += temp_sum;
}
double calculate_func(double x1, double x2){
    const double FUN_CONST = 0.002;
    double temp_sum = 0.0;

    for(int i = -2; i <= 2; i++){
        for(int j=-2; j <= 2; j++){
            temp_sum += 1 / (5 * ( i + 2 ) + j + 3 + pow((x1 - 16 * j), 6.0) + pow((x2 - 16 * i), 6.0) );
        }
    }
    temp_sum += FUN_CONST;
    temp_sum = pow(temp_sum, -1.0);
    return  temp_sum;
}

inline std::chrono::high_resolution_clock::time_point get_current_time_fenced(){
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto res_time = std::chrono::high_resolution_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return res_time;
}

template<class D>
inline long long to_us(const D& d) {
    return std::chrono::duration_cast<std::chrono::microseconds>(d).count();
}

int main() {

    const std::string symbols_c = "=\"";
    const std::string config_name("../addition/config.txt");
    const std::string start_x_name("interval_start_x");
    const std::string start_y_name("interval_start_y");
    const std::string end_x_name("interval_end_x");
    const std::string end_y_name("interval_end_y");
    const std::string result_name("res");
    const std::string threads_name("threads");
    double start_x;
    double start_y;
    double end_x;
    double end_y;
    double sum = 0.0;
    double step_x = 1.0;
    double step_y;
    int threads_num = 0;
    std::ifstream file_config(config_name);
    std::ofstream file_result;
    std::vector<std::string>config_text;
    std::vector<std::thread*> threads;
    std::map<std::string, std::string> config_dict;
    std::mutex sum_mutex;
    std::chrono::high_resolution_clock::time_point count_start_time;
    std::chrono::high_resolution_clock::time_point count_finish_time;


    /*
     * Read config file and create config dictionary
     */

    read(file_config, 8, &config_text);
    std::for_each(config_text.begin(), config_text.end(), [&symbols_c](std::string &line) {
        reduce_symbols(line, symbols_c);
    });
    std::for_each(config_text.begin(), config_text.end(), [&config_dict](std::string &line) {
        std::vector<std::string> words = split(line);
        config_dict[words.at(0)] = words.at(1);
    });

    start_x = std::stod(config_dict[start_x_name]);
    start_y = std::stod(config_dict[start_y_name]);
    end_x = std::stod(config_dict[end_x_name]);
    end_y = std::stod(config_dict[end_y_name]);
    threads_num = std::stoi(config_dict[threads_name]);

    /*
     * Start of calculation 
     */

    count_start_time = get_current_time_fenced();

    step_y = ( end_y - start_y ) / threads_num;

    for(int i = 0; i < threads_num; i++){

        double temp_start_y = start_y + step_y * i;
        double temp_end_y = temp_start_y + step_y;

        threads.push_back(new std::thread (thread_handler,  start_x,  end_x,  temp_start_y,  temp_end_y, step_x, std::ref(sum), std::ref(sum_mutex)));
    }

    std::for_each (threads.begin(), threads.end(), [](std::thread *thread){
        thread->join();
    });

    count_finish_time = get_current_time_fenced();

    /*
     * Finish of calculation (above)
     */

    std::cout << "Calculation: "<< to_us(count_finish_time - count_start_time) << std::endl;


    return 0;
}