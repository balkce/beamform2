/**
 * MCRA filter.
 * 
 * Based on:
 * Israel Cohen and Baruch Berdugo. Noise Estimation by Minima Controlled Recursive Averaging for Robust Speech Enhancement, IEEE SIGNAL PROCESSING LETTERS, VOL. 9, NO. 1, JANUARY 2002
 * 
 * And code from:
 * https://github.com/introlab/manyears/blob/master/manyears-C/dsplib/Preprocessing/mcra.c
 */

#include "rosjack/rosjack.hpp"
#include "util.h"

// Include FFTW header
#include <complex>
#include <fftw3.h>

// Eigen include
#include <Eigen/Eigen>

bool READY = false;

std::complex<double> *x_fft, *x_time, *y_fft, *y_time;
fftw_plan x_forward, y_inverse;

double *freqs;

// MCRA stuff
double smoothing_freq_window[3] = {0.25, 0.5, 0.25};
int smoothing_freq_window_pos[3] = {-1, 0, 1};
int smoothing_freq_window_size = 3;
int current_L = 0;
bool first_L = true;
double *S;
double *S_prev;
double *S_f;
double *S_tmp;
double *S_min;
double past_mag = 0;

// MCRA parameters
double alphaS = 0.95;
double alphaD = 0.95;
double alphaD2 = 0.97;
double delta = 0.001;
bool out_only_noise = false;
int L = 75; //this may need to be re-calculated: it's the number of windows between minima searches, and should represent between 0.5 and 1.5 s.
double out_amp = 2.0;

double *lambda; //noise estimation

//reused buffer
Eigen::MatrixXcd in_fft;
Eigen::MatrixXd in_fft_square;

double min (double a, double b){
    if(a > b){
        return b;
    }else{
        return a;
    }
}

void mcra (jack_ringbuffer_t **in, rosjack_data *out){
    
    int i;
    unsigned j,this_j;
    double this_mag, this_pha;
    
    // fft
    //only uses the first channel;
    overlap_and_add_prepare_input(in[0], x_time);
    fftw_execute(x_forward);
    for(j = 0; j < fft_win; j++){
        in_fft(0,j) = x_fft[j];
        in_fft_square(0,j) = norm(x_fft[j]);
    }
    
    //doing MCRA
    
    //smoothing frequency vector
    S_f[0] = abs(in_fft(0,0)); //passing on the DC component
    for(j = 1; j < fft_win; j++){
        S_f[j] = 0.0;
        for(i=0;i<smoothing_freq_window_size;i++){
            this_j = j+smoothing_freq_window_pos[i];
            if(this_j >= 1 && this_j < fft_win){
                S_f[j] += smoothing_freq_window[i]*in_fft_square(0,this_j);
            }
        }
    }
    
    //smoothing frequency vector through time
    for(j = 0; j < fft_win; j++){
        S[j] = (alphaS*S_prev[j]) + ((1-alphaS)*S_f[j]);
    }
    
    //carrying out minima search
    if(current_L > L){
        for(j = 0; j < fft_win; j++){
            S_min[j] = min(S_tmp[j],S[j]);
            S_tmp[j] = S[j];
        }
        current_L = 1;
        first_L = false;
    }else{
        for(j = 0; j < fft_win; j++){
            S_min[j] = min(S_min[j],S[j]);
            S_tmp[j] = min(S_tmp[j],S[j]);
        }
        current_L++;
    }
    
    //updating noise estimation
    for(j = 0; j < fft_win; j++){
        if (first_L || S[j] < S_min[j]*delta || lambda[j] > in_fft_square(0,j)){
            if (((first_L) && ((1.0f/(double)current_L) > alphaD))){
                lambda[j] = (1.0f/(double)current_L) * lambda[j] + (1.0f - (1.0f/(double)current_L)) * in_fft_square(0,j);
            }else{
                lambda[j] = alphaD2 * lambda[j] + (1.0f - alphaD) * in_fft_square(0,j);
            }
        }
    }
    
    //removing noise (lambda) from input
    y_fft[j] = in_fft(0,0);
    for(j = 1; j < fft_win; j++){
        this_pha = arg(in_fft(0,j));
        
        if (out_only_noise){
          this_mag = (sqrt(lambda[j]))*out_amp; //output only the noise
        }else{
          this_mag = (abs(in_fft(0,j))-sqrt(lambda[j]))*out_amp; //removing the noise magnitude
          if (this_mag < 0)
              this_mag = 0.0;
        }
            
        y_fft[j] = std::complex<double>(this_mag*cos(this_pha),this_mag*sin(this_pha));
        
        //if (j == inspect_j) 
        //    printf("%d \t noise mag at %f Hz: %f\n",current_L,freqs[inspect_j],(sqrt(lambda[j]))*out_amp); fflush(stdout);
    }
    
    //storing current S in S_prev for next iteration
    for(j = 0; j < fft_win; j++){
        S_prev[j] = S[j];
    }
    
    // ifft
    fftw_execute(y_inverse);

    // preparing output
    overlap_and_add_prepare_output(y_time,out);
}

int jack_callback (jack_nframes_t nframes, void *arg){
    //TimeVar t_bef = timeNow();
    jack_nframes_t i;
    
    //Inputing from ROS
    rosjack_data *out = new rosjack_data[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        do_overlap(in, out, nframes, mcra);
    }else{
        for (i = 0; i < nframes; i++){
            out[i] = 0.0;
        }
    }
    
    //Outputing to ROS
    output_to_rosjack (out, nframes, output_type);
    
    //std::cout << "Callback took: " << duration(timeNow()-t_bef)/1000000.0 << " ms.\n";
    return 0;
}

void mcra_handle_params(std::shared_ptr<rclcpp::Node> n){
    std::cout << "MCRA ROS parameters: " << std::endl;
    
    n->declare_parameter("alphaS",0.95);
    n->declare_parameter("alphaD",0.95);
    n->declare_parameter("alphaD2",0.97);
    n->declare_parameter("delta",0.001);
    n->declare_parameter("L",75);
    n->declare_parameter("out_amp",2.0);
    n->declare_parameter("out_only_noise",true);

    if (n->get_parameter("alphaS",alphaS)){
        RCLCPP_INFO(n->get_logger(),"Alpha for S: %f",alphaS);
    }else{
        alphaS = 0.95;
        RCLCPP_WARN(n->get_logger(),"Alpha for S argument not found in ROS param server, using default value (%f).",alphaS);
    }
    
    if (n->get_parameter("alphaD",alphaD)){
        RCLCPP_INFO(n->get_logger(),"Alpha for D: %f",alphaD);
    }else{
        alphaD = 0.95;
        RCLCPP_WARN(n->get_logger(),"Alpha for D argument not found in ROS param server, using default value (%f).",alphaD);
    }
    
    if (n->get_parameter("alphaD2",alphaD2)){
        RCLCPP_INFO(n->get_logger(),"Alpha2 for D: %f",alphaD2);
    }else{
        alphaD2 = 0.97;
        RCLCPP_WARN(n->get_logger(),"Alpha2 for D argument not found in ROS param server, using default value (%f).",alphaD2);
    }
    
    if (n->get_parameter("delta",delta)){
        RCLCPP_INFO(n->get_logger(),"Delta: %f",delta);
    }else{
        delta = 0.001;
        RCLCPP_WARN(n->get_logger(),"Delta argument not found in ROS param server, using default value (%f).",delta);
    }
    
    if (n->get_parameter("L",L)){
        RCLCPP_INFO(n->get_logger(),"Training Windows (L): %d",L);
    }else{
        L = 75;
        RCLCPP_WARN(n->get_logger(),"Training Windows (L) argument not found in ROS param server, using default value (%d).",L);
    }
    if (n->get_parameter("out_amp",out_amp)){
        RCLCPP_INFO(n->get_logger(),"Output Amplification: %f",out_amp);
    }else{
        out_amp = 2.0;
        RCLCPP_WARN(n->get_logger(),"Output Amplification argument not found in ROS param server, using default value (%f).",out_amp);
    }
    
    if (n->get_parameter("out_only_noise",out_only_noise)){
        RCLCPP_INFO(n->get_logger(),"Output only noise: %d",out_only_noise);
    }else{
        out_only_noise = true;
        RCLCPP_WARN(n->get_logger(),"Noise output argument not found in ROS param server, outputting filtered signal by default.");
    }
    
}

int main (int argc, char *argv[]) {
    /* ROS initialization*/
    rclcpp::init(argc, argv);
    std::shared_ptr<rclcpp::Node> rosjack_node = rclcpp::Node::make_shared("beamform2");
    handle_params(rosjack_node);
    mcra_handle_params(rosjack_node);
    
    
    /* create JACK agent */
    if(rosjack_create (ROSJACK_READ, rosjack_node, "jackaudio", rosjack_node->get_name(), number_of_microphones, jack_callback)){
        RCLCPP_ERROR(rosjack_node->get_logger(),"JACK agent could not be created.\n");
        rclcpp::shutdown();
        exit(1);
    }
    
    std::cout << "Pre-allocating space for internal buffers." << std::endl;
    prepare_overlap_and_add(); //fft_win is assigned here
    
    x_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    x_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    
    x_forward = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(x_time), reinterpret_cast<fftw_complex*>(x_fft), FFTW_FORWARD, FFTW_MEASURE);
    y_inverse = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(y_fft), reinterpret_cast<fftw_complex*>(y_time), FFTW_BACKWARD, FFTW_MEASURE);
    
    freqs = (double *)malloc(sizeof(double)*fft_win);
    calculate_frequency_vector(freqs,fft_win);
    
    S = (double *)malloc(sizeof(double)*fft_win);
    S_prev = (double *)malloc(sizeof(double)*fft_win);
    S_f = (double *)malloc(sizeof(double)*fft_win);
    S_tmp = (double *)malloc(sizeof(double)*fft_win);
    S_min = (double *)malloc(sizeof(double)*fft_win);
    lambda = (double *)malloc(sizeof(double)*fft_win);
    
    //initializing vectors
    for(unsigned int i = 0; i < fft_win; i++){
        S_prev[i] = 0.0;
        S_tmp[i] = 0.0;
        S_min[i] = 0.0;
        lambda[i] = 0.0;
    }
    
    in_fft.resize(number_of_microphones,fft_win);
    in_fft_square.resize(number_of_microphones,fft_win);
    
    READY = true;
    
    RCLCPP_INFO(rosjack_node->get_logger(),"Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    rclcpp::spin(rosjack_node);
    
    exit(0);
}
