/**
 * Mininimum-variance distortionless response, with frequency magnitude
 * thresholding so that it can be carried out online.
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
double *delays;
Eigen::MatrixXcd weights(1,1);
Eigen::MatrixXcd weights_h(1,1);

//MVDR config
unsigned int past_windows = 10;
double freq_mag_threshold = 1.5;
double freq_max = 4000;
double freq_min = 400;
double out_amp = 4.5;

//reused buffers
Eigen::MatrixXcd in_fft(1,1);
Eigen::MatrixXcd R(1,1);
Eigen::MatrixXcd invR(1,1);
Eigen::MatrixXcd whiteR(1,1);
Eigen::MatrixXcd MVDR_w(1,1);
std::vector<Eigen::MatrixXcd> past_ffts;

void update_weights(bool ini=false){
    calculate_delays(delays);
    
    long unsigned int i;
    unsigned int j;
    
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            if(ini){
                for(j = 0; j < fft_win; j++){
                    weights(i,j) = 1.0; 
                }
            }
        }else{
            for(j = 0; j < fft_win; j++){
                weights(i,j) = std::exp(-M_I*(double)2*PI*freqs[j]*delays[i]); 
            }
        }
    }
    weights_h = weights.adjoint();
}

void apply_weights (jack_ringbuffer_t **in, rosjack_data *out){
    int i;
    unsigned int j;
    double this_freq, this_mag;
    
    // fft
    for(i = 0; i < number_of_microphones; i++){
        overlap_and_add_prepare_input(in[i], x_time);
        fftw_execute(x_forward);
        for(j = 0; j < fft_win; j++){
            in_fft(i,j) = x_fft[j];
        }
    }
    
    //applying weights
    y_fft[0] = in_fft(0,0);
    for(j = 1; j < fft_win; j++){
        this_freq = abs(freqs[j]);
        this_mag = 0.0;
        for(i = 0; i < number_of_microphones; i++)
            this_mag += abs(in_fft(i,j));
        this_mag /= number_of_microphones*fft_win; //normalize the average against window length
        
        if(this_freq >= freq_min && this_freq <= freq_max){
            if(this_mag > freq_mag_threshold){
                //calculating sample covariance for this fft bin
                R = (past_ffts[j]*past_ffts[j].adjoint()).cwiseProduct(whiteR);
                invR = R.inverse();
                
                //calculating MVDR optimal weights for this fft bin
                MVDR_w = (invR*weights.col(j))/((weights_h.row(j)*invR*weights.col(j)));
                
                //applying weigths
                y_fft[j] = (MVDR_w.adjoint()*in_fft.col(j))(0,0);
            }else{
                y_fft[j] = in_fft(0,j)*0.01;
            }
            
            //shifting past_ffts and appending this window
            past_ffts[j].block(0,0,number_of_microphones,past_windows-1) = past_ffts[j].block(0,1,number_of_microphones,past_windows-1);
            past_ffts[j].col(past_windows-1) = in_fft.col(j);
        }else{
            y_fft[j] = 0.0;
        }
    }
    
    // ifft
    fftw_execute(y_inverse);
    
    // preparing output
    overlap_and_add_prepare_output(y_time,out);
    for (j = 0; j<fft_win; j++){
        out[j] *= out_amp;
    }
}

int jack_callback (jack_nframes_t nframes, void *arg){
    //TimeVar t_bef = timeNow();
    jack_nframes_t i;
    
    //Inputing from ROS
    rosjack_data *out = new rosjack_data[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        do_overlap(in, out, nframes, apply_weights);
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

void theta_roscallback(const std_msgs::msg::Float32 msg){
    std::cout << "Updating weights for angle: " << msg.data << std::endl;
    
    angle = msg.data;
    update_weights();
}

void mvdr_handle_params(std::shared_ptr<rclcpp::Node> n){
    std::cout << "MVDR ROS parameters: " << std::endl;
    
    n->declare_parameter("past_windows",10);
    n->declare_parameter("freq_mag_threshold",1.5);
    n->declare_parameter("freq_max",4000.0);
    n->declare_parameter("freq_min",400.0);
    n->declare_parameter("out_amp",4.5);
    
    if (n->get_parameter("past_windows",past_windows)){
        RCLCPP_INFO(n->get_logger(),"Number of Windows for Covariance Calculation: %d",past_windows);
    }else{
        past_windows = 10;
        RCLCPP_WARN(n->get_logger(),"Number of Windows for Covariance Calculation argument not found in ROS param server, using default value (%d).",past_windows);
    }
    
    if (n->get_parameter("freq_mag_threshold",freq_mag_threshold)){
        RCLCPP_INFO(n->get_logger(),"Frequency Magnitude Threshold: %f",freq_mag_threshold);
    }else{
        freq_mag_threshold = 1.5;
        RCLCPP_WARN(n->get_logger(),"Frequency Magnitude Threshold argument not found in ROS param server, using default value (%f).",freq_mag_threshold);
    }
    
    if (n->get_parameter("freq_max",freq_max)){
        RCLCPP_INFO(n->get_logger(),"Max Frequency: %f",freq_max);
    }else{
        freq_max = 4000;
        RCLCPP_WARN(n->get_logger(),"Max Frequency argument not found in ROS param server, using default value (%f).",freq_max);
    }
    
    if (n->get_parameter("freq_min",freq_min)){
        RCLCPP_INFO(n->get_logger(),"Min Frequency: %f",freq_min);
    }else{
        freq_min = 400;
        RCLCPP_WARN(n->get_logger(),"Min Frequency argument not found in ROS param server, using default value (%f).",freq_min);
    }
    
    if (n->get_parameter("out_amp",out_amp)){
        RCLCPP_INFO(n->get_logger(),"Output Amplification: %f",out_amp);
    }else{
        out_amp = 4.5;
        RCLCPP_WARN(n->get_logger(),"Output Amplification argument not found in ROS param server, using default value (%f).",out_amp);
    }
    
}


int main (int argc, char *argv[]) {
    /* ROS initialization*/
    rclcpp::init(argc, argv);
    std::shared_ptr<rclcpp::Node> rosjack_node = rclcpp::Node::make_shared("beamform2");
    handle_params(rosjack_node);
    mvdr_handle_params(rosjack_node);
    
    rclcpp::Subscription<std_msgs::msg::Float32>::SharedPtr theta_subscriber = rosjack_node->create_subscription<std_msgs::msg::Float32>("theta", 1000, theta_roscallback);
    
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
    
    delays = (double *) malloc (sizeof(double)*number_of_microphones);
    
    weights.resize(number_of_microphones,fft_win);
    weights_h.resize(fft_win,number_of_microphones);
    in_fft.resize(number_of_microphones,fft_win);
    
    for(unsigned int i = 0; i<fft_win;i++){
        Eigen::MatrixXcd this_past_fs(number_of_microphones,past_windows);
        this_past_fs.setZero();
        past_ffts.push_back(this_past_fs);
    }
    
    MVDR_w.resize(number_of_microphones,1);
    
    R.resize(number_of_microphones,number_of_microphones);
    invR.resize(number_of_microphones,number_of_microphones);
    
    whiteR.resize(number_of_microphones,number_of_microphones);
    whiteR.setOnes();
    for(int i = 0; i<number_of_microphones;i++){
        whiteR(i,i) = 1.001;
    }
    
    update_weights(true);
    
    READY = true;
    
    RCLCPP_INFO(rosjack_node->get_logger(),"Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    rclcpp::spin(rosjack_node);
    
    std::cout << "Closing theta subscriber." << std::endl;
    theta_subscriber.reset();
    
    exit(0);
}
