/**
 * Linearly-Constrained Mininimum Variance, with frequency magnitude
 * thresholding so that it can be carried out online.
 */

#include "rosjack/rosjack.hpp"
#include "beamform2/msg/interf_theta.hpp"
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
double **interf_delays;
std::vector<Eigen::MatrixXcd> weights;
std::vector<Eigen::MatrixXcd> weights_h;

//LCMV config
unsigned int past_windows = 10;
double freq_mag_threshold = 1.5;
double freq_max = 4000;
double freq_min = 400;
double out_amp = 4.5;
double interf_angle_threshold = 5.0;

//reused buffers
Eigen::MatrixXcd in_fft(1,1);
Eigen::MatrixXcd R(1,1);
Eigen::MatrixXcd invR(1,1);
Eigen::MatrixXcd whiteR(1,1);
Eigen::MatrixXcd LCMV_w(1,1);
std::vector<Eigen::MatrixXcd> past_ffts;

void update_weights(bool ini=false){
    calculate_delays(delays);
    
    long unsigned int i,k;
    unsigned int j;
    
    //weights for DOI
    for(i = 0; i < array_geometry.size(); i++){
        if (i == 0){
            if(ini){
                for(j = 0; j < fft_win; j++){
                    weights[j](i,0) = 1.0; 
                }
            }
        }else{
            for(j = 0; j < fft_win; j++){
                weights[j](i,0) = std::exp(-M_I*(double)2*PI*freqs[j]*delays[i]); 
            }
        }
    }
    
    //weights for interferences
    for(k = 0; k < interference_angles.size(); k++){
        calculate_interf_delays(interf_delays[k],k,interference_angles[k]);
        for(i = 0; i < array_geometry.size(); i++){
            if (i == 0){
                if(ini){
                    for(j = 0; j < fft_win; j++){
                        weights[j](i,k+1) = 1.0; 
                    }
                }
            }else{
                for(j = 0; j < fft_win; j++){
                    weights[j](i,k+1) = std::exp(-M_I*(double)2*PI*freqs[j]*interf_delays[k][i]); 
                }
            }
        }
    }
    
    //adjoint weights 
    for(j = 0; j < fft_win; j++){
        weights_h[j] = weights[j].adjoint();
    }
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
    for(j = 0; j < fft_win; j++){
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
                
                //calculating LCMV optimal weights for this fft bin
                LCMV_w = (invR*weights[j])*(weights_h[j]*invR*weights[j]).inverse();
                
                //applying weigths
                y_fft[j] = (LCMV_w.col(0).adjoint()*in_fft.col(j))(0,0);
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
    printf("Updating weights for angle: %f", msg.data);
    
    angle = msg.data;
    update_weights();
}

void lcmv_handle_params(std::shared_ptr<rclcpp::Node> n){
    std::cout << "LCMV ROS parameters: " << std::endl;
    
    n->declare_parameter("past_windows",10);
    n->declare_parameter("freq_mag_threshold",1.5);
    n->declare_parameter("freq_max",4000.0);
    n->declare_parameter("freq_min",400.0);
    n->declare_parameter("out_amp",4.5);
    n->declare_parameter("interf_angle_threshold",5.0);

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
    
    if (n->get_parameter("interf_angle_threshold",interf_angle_threshold)){
        RCLCPP_INFO(n->get_logger(),"Angular Threshold Between Interferences: %f",interf_angle_threshold);
    }else{
        interf_angle_threshold = 5.0;
        RCLCPP_WARN(n->get_logger(),"Angular Threshold Between Interferences argument not found in ROS param server, using default value (%f).",interf_angle_threshold);
    }
    
}

void free_interf_buffers(){
    long unsigned int i;
    
    for(i = 0; i<interference_angles.size();i++){
        free(interf_delays[i]);
    }
    free(interf_delays);
    
    std::vector<Eigen::MatrixXcd>().swap(weights);
    std::vector<Eigen::MatrixXcd>().swap(weights_h);
}

void allocate_interf_buffers(){
    //buffers to update: interf_delays, weights, weights_h (the last two could be recreated using update_weights), LCMV_w
    
    interf_delays = (double **) malloc (sizeof(double*)*interference_angles.size());
    for(long unsigned int i = 0; i<interference_angles.size();i++){
        interf_delays[i] = (double *) malloc (sizeof(double)*number_of_microphones);
    }
    
    
    for(unsigned int i = 0; i<fft_win;i++){
        Eigen::MatrixXcd this_weights(number_of_microphones,interference_angles.size()+1);
        Eigen::MatrixXcd this_weights_h(interference_angles.size()+1,number_of_microphones);
        
        this_weights.setZero();
        this_weights_h.setZero();
        
        weights.push_back(this_weights);
        weights_h.push_back(this_weights_h);
    }
    
    LCMV_w.resize(number_of_microphones,interference_angles.size()+1);
    
}

void interf_theta_roscallback(const beamform2::msg::InterfTheta msg){
    if (msg.id >= 1 && msg.id <= interference_angles.size()){
        printf("Updating weights for interference %d to angle %f", msg.id, msg.angle);
        
        interference_angles[msg.id-1] = msg.angle;
        
        //check if the new angle is too close to another interference
        long unsigned int i;
        for(i = 0; i<interference_angles.size();i++){
            if ((int)i != (msg.id-1) && abs(interference_angles[i]-msg.angle) < interf_angle_threshold){
                printf("New angle of interference %d is too close to the angle of interference %ld: %f", msg.id, i+1, interference_angles[i]);
                printf("Removing interference %d ", msg.id);
                
                READY = false;
                millisleep(30);
                free_interf_buffers();
                interference_angles.erase(interference_angles.begin()+msg.id-1);
                allocate_interf_buffers();
                READY = true;
                break;
            }
        }
        
        update_weights();
    }else if(msg.id > interference_angles.size()){
        printf("Interference id (%d) is greater than current number of interferences (%d)", msg.id, (int)interference_angles.size());
        printf("Assuming that it is a new interfernce...");
        
        long unsigned int i;
        for(i = 0; i<interference_angles.size();i++){
            if (abs(interference_angles[i]-msg.angle) < interf_angle_threshold){
                printf("Angle of new interference %d is too close to the angle of interference %ld: %f", ((int)interference_angles.size())+1, i+1, interference_angles[i]);
                printf("Not adding interference %d ", msg.id);
                break;
            }
        }
        
        if (i == interference_angles.size()){
            printf("Adding new interference %d in angle %f", ((int)interference_angles.size())+1, msg.angle);
            READY = false;
            millisleep(30);
            free_interf_buffers();
            interference_angles.push_back(msg.angle);
            allocate_interf_buffers();
            READY = true;
            
            update_weights();
        }
    }else{
        printf("Invalid interference id (%d), should be >= 1", msg.id);
    }
}

int main (int argc, char *argv[]) {
    /* ROS initialization*/
    rclcpp::init(argc, argv);
    std::shared_ptr<rclcpp::Node> rosjack_node = rclcpp::Node::make_shared("beamform2");
    handle_params(rosjack_node);
    lcmv_handle_params(rosjack_node);
    
    rclcpp::Subscription<std_msgs::msg::Float32>::SharedPtr theta_subscriber = rosjack_node->create_subscription<std_msgs::msg::Float32>("theta", 1000, theta_roscallback);
    rclcpp::Subscription<beamform2::msg::InterfTheta>::SharedPtr interf_theta_subscriber = rosjack_node->create_subscription<beamform2::msg::InterfTheta>("theta_interference", 1000, interf_theta_roscallback);
    
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
    
    in_fft.resize(number_of_microphones,fft_win);
    
    for(unsigned int i = 0; i<fft_win;i++){
        Eigen::MatrixXcd this_past_fs(number_of_microphones,past_windows);
        this_past_fs.setZero();
        past_ffts.push_back(this_past_fs);
    }
    
    R.resize(number_of_microphones,number_of_microphones);
    invR.resize(number_of_microphones,number_of_microphones);
    
    whiteR.resize(number_of_microphones,number_of_microphones);
    whiteR.setOnes();
    for(int i = 0; i<number_of_microphones;i++){
        whiteR(i,i) = 1.001;
    }
    
    allocate_interf_buffers();
    
    update_weights(true);
    
    READY = true;
    
    RCLCPP_INFO(rosjack_node->get_logger(),"Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    rclcpp::spin(rosjack_node);
    
    std::cout << "Closing theta subscriber." << std::endl;
    theta_subscriber.reset();
    interf_theta_subscriber.reset();
    
    exit(0);
}
