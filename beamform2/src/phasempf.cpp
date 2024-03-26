/**
 * Beamform that carries out phase-based frequency masking with a post-filter.
 * Based on:
 * Valin, J.M.; Yamamoto, S.; Rouat, Jean; Michaud, F.; Nakadai, K.; Okuno, H.G., "Robust Recognition of Simultaneous Speech by a Mobile Robot," IEEE Transactions on Robotics, vol.23, no.4, pp.742,752, Aug. 2007.
 * 
 * And code from:
 * https://github.com/introlab/manyears/blob/master/manyears-C/dsplib/Separation/postfilter.c
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
Eigen::MatrixXcd weights;

// Phase parameters
double min_mag = 0.0001;
double min_phase = 10;
double min_phase_diff_mean = 10*M_PI/180;
int smooth_size = 10;

// MCRA parameters
double MCRA_alphaS = 0.95;
double MCRA_alphaD = 0.95;
double MCRA_alphaD2 = 0.97;
double MCRA_delta = 0.001;
int MCRA_L = 75; //this may need to be re-calculated: it's the number of windows between minima searches, and should represent between 0.5 and 1.5 s.

// Other MCRA stuff
double smoothing_freq_window[3] = {0.25, 0.5, 0.25};
int smoothing_freq_window_pos[3] = {-1, 0, 1};
int smoothing_freq_window_size = 3;
int current_L = 0;
bool first_L = true;

//MPF parameters
double MPF_alphaS = 0.3;
double MPF_eta = 0.3;
double MPF_rev_gamma = 0.3;
double MPF_rev_delta = 1.0;

// Additional output arguments
double out_amp = 1.0;
double noise_floor = 0.001;
bool out_only_noise = false;
bool out_only_mcra = false;

//reused buffers
double *phases_aligned;
double *past_samples;
Eigen::MatrixXcd in_fft;
std::complex<double> *out_soi, *out_int;
double *out_soi_square, *out_int_square;
double *MCRA_S;
double *MCRA_S_prev;
double *MCRA_S_f;
double *MCRA_S_tmp;
double *MCRA_S_min;
double *MPF_Z;
double *MCRA_lambda_noise; //noise estimation
double *MPF_lambda_leak; //leak estimation
double **MPF_lambda_rev; //leak estimation
double *MPF_lambda; //noise estimation

void shift_data(double data, double *buf, int size){
  for(int i = 1; i < size; i++){
    buf[i-1] = buf[i];
  }
  buf[size-1] = data;
}

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
}

double get_overall_phase_diff(int min_i,int *num_i){
    if (min_i < number_of_microphones-1){
        double this_diff = 0;
        double this_diff_raw;
        for (int i = min_i+1; i < number_of_microphones; i++){
            this_diff_raw = abs(phases_aligned[min_i]-phases_aligned[i]);
            if (this_diff_raw > M_PI)
                this_diff_raw = 2*M_PI - this_diff_raw;
            this_diff += this_diff_raw;
            (*num_i)++;
        }
        return this_diff + get_overall_phase_diff(min_i+1,num_i);
    }else{
        return 0;
    }
}

double get_mean(double * data, int data_size){
    double data_sum = 0.0;
    
    for(int i = 0; i < data_size; i++){
        data_sum += data[i];
    }
    
    return data_sum/(double)data_size;
}

double min (double a, double b){
    if(a > b){
        return b;
    }else{
        return a;
    }
}

void mcra (){
    int i;
    unsigned int j,this_j;
    
    //smoothing frequency vector
    MCRA_S_f[0] = abs(out_soi[0]); //passing on the DC component
    for(j = 1; j < fft_win; j++){
        MCRA_S_f[j] = 0.0;
        for(i=0;i<smoothing_freq_window_size;i++){
            this_j = j+smoothing_freq_window_pos[i];
            if(this_j >= 1 && this_j < fft_win){
                MCRA_S_f[j] += smoothing_freq_window[i]*out_soi_square[j];
            }
        }
    }
    
    //smoothing frequency vector through time
    for(j = 0; j < fft_win; j++){
        MCRA_S[j] = (MCRA_alphaS*MCRA_S_prev[j]) + ((1-MCRA_alphaS)*MCRA_S_f[j]);
    }
    
    //carrying out minima search
    if(current_L > MCRA_L){
        for(j = 0; j < fft_win; j++){
            MCRA_S_min[j] = min(MCRA_S_tmp[j],MCRA_S[j]);
            MCRA_S_tmp[j] = MCRA_S[j];
        }
        current_L = 1;
        first_L = false;
    }else{
        for(j = 0; j < fft_win; j++){
            MCRA_S_min[j] = min(MCRA_S_min[j],MCRA_S[j]);
            MCRA_S_tmp[j] = min(MCRA_S_tmp[j],MCRA_S[j]);
        }
        current_L++;
    }
    
    //updating noise estimation
    for(j = 0; j < fft_win; j++){
        if (first_L || MCRA_S[j] < MCRA_S_min[j]*MCRA_delta || MCRA_lambda_noise[j] > out_soi_square[j]){
            if (((first_L) && ((1.0f/(double)current_L) > MCRA_alphaD))){
                MCRA_lambda_noise[j] = (1.0f/(double)current_L) * MCRA_lambda_noise[j] + (1.0f - (1.0f/(double)current_L)) * out_soi_square[j];
            }else{
                MCRA_lambda_noise[j] = MCRA_alphaD2 * MCRA_lambda_noise[j] + (1.0f - MCRA_alphaD) * out_soi_square[j];
            }
        }
    }
    
    //storing current S in S_prev for next iteration
    for(j = 0; j < fft_win; j++){
        MCRA_S_prev[j] = MCRA_S[j];
    }
}

void apply_weights (jack_ringbuffer_t **in, rosjack_data *out){
    int i;
    unsigned int j;
    double phase_diff_sum;
    int phase_diff_num;
    double phase_diff_mean;
    double mag_mean;
    double pha_mean;
    
    // fft
    for(i = 0; i < number_of_microphones; i++){
        overlap_and_add_prepare_input(in[i], x_time);
        fftw_execute(x_forward);
        for(j = 0; j < fft_win; j++){
            in_fft(i,j) = x_fft[j];
        }
    }
    
    out_soi[0] = in_fft(0,0);
    out_int[0] = in_fft(0,0);
    for(j = 1; j < fft_win; j++){
        //applying weights to align phases
        for(i = 0; i < number_of_microphones; i++){
          phases_aligned[i] = arg(conj(weights(i,j))*in_fft(i,j));
        }
        
        //getting the mean phase difference between all microphones
        phase_diff_num = 0;
        phase_diff_sum = get_overall_phase_diff(0,&phase_diff_num);
        phase_diff_mean = phase_diff_sum/(double)phase_diff_num;
        
        //creating new frequency data bin from mean magnitude
        mag_mean = 0;
        for(i = 0; i < number_of_microphones; i++){
            mag_mean += abs(in_fft(i,j));
        }
        mag_mean /= number_of_microphones;
        
        //and from the phase of the reference microphone
        pha_mean = arg(in_fft(0,j));
        
        
        if (phase_diff_mean < min_phase_diff_mean){
            out_soi[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
            
            mag_mean *= min_mag;
            out_int[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
        }else{
            out_int[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
            
            mag_mean *= min_mag;
            out_soi[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
        }
        
        out_soi_square[j] = norm(out_soi[j]);
        out_int_square[j] = norm(out_int[j]);
    }
    
    //estimate noise by MCRA
    mcra (); //carries out noise estimation by MCRA only with the SOI channel (out_soi)
             //stores the result in MCRA_lambda_noise
    
    //estimate total noise variance
    for(j = 0; j < fft_win; j++){
        //update Z with previous noise channel (out_int) and current noise channel power
        //Z_m(k,l) = alphaS * Z_m,(k,l-1) + (1 - alphaS) * |Y_m(k,l)|^2
        MPF_Z[j] = MPF_alphaS * MPF_Z[j] + (1 - MPF_alphaS) * out_int_square[j];
        
        //applying reduction factor to estimate leakage
        MPF_lambda_leak[j] = MPF_eta*MPF_Z[j];
        
        //updating reverberation estimation from both channels
        //lambda_rev(k,l) = gamma * lambda_rev(k,l-1) + ((1-gamma)/delta) * |Y_m(k,l)|^2
        MPF_lambda_rev[0][j] = MPF_rev_gamma * MPF_lambda_rev[0][j] + ((1-MPF_rev_gamma/MPF_rev_delta)) * out_soi_square[j];
        MPF_lambda_rev[1][j] = MPF_rev_gamma * MPF_lambda_rev[1][j] + ((1-MPF_rev_gamma/MPF_rev_delta)) * out_int_square[j];
        
        //updating noise variance
        MPF_lambda[j] = MCRA_lambda_noise[j] + MPF_lambda_leak[j] + MPF_lambda_rev[0][j] + MPF_lambda_rev[1][j];
        MPF_lambda[j] = sqrt(MPF_lambda[j]);
    }
    
    //removing noise/leakage/reverberation (MPF_lambda) from SOI (out_soi)
    y_fft[j] = out_soi[0];
    for(j = 1; j < fft_win; j++){
        pha_mean = arg(out_soi[j]);
        
        if (out_only_noise){
          mag_mean = MPF_lambda[j]*out_amp; //output only the noise
        }else{
          if (out_only_mcra){
            mag_mean = (abs(out_soi[j])-sqrt(MCRA_lambda_noise[j]))*out_amp; //removing the MCRA noise estimate
          }else{
            mag_mean = (abs(out_soi[j])-MPF_lambda[j])*out_amp; //removing the whole MPF noise estimate
          }
          
          if (mag_mean < 0)
              mag_mean = noise_floor;
        }
            
        y_fft[j] = std::complex<double>(mag_mean*cos(pha_mean),mag_mean*sin(pha_mean));
        
        //if (j == inspect_j) 
        //    printf("%d \t noise mag at %f Hz: %f\n",current_L,freqs[inspect_j],(sqrt(lambda[j]))*out_amp); fflush(stdout);
    }
    
    // ifft
    fftw_execute(y_inverse);
    
    // preparing output
    overlap_and_add_prepare_output(y_time,out);
}

void shift_data_rosjack(rosjack_data data, rosjack_data *buf, int size){
  for(int i = 1; i < size; i++){
    buf[i-1] = buf[i];
  }
  buf[size-1] = data;
}

rosjack_data calculate_power(rosjack_data *buf, int size){
  rosjack_data p = 0;
  for(int i = 0; i < size; i++){
    p += buf[i]*buf[i];
  }
  return sqrt(p/(rosjack_data)size);
}

int jack_callback (jack_nframes_t nframes, void *arg){
    //TimeVar t_bef = timeNow();
    jack_nframes_t i;
    
    //Inputing from ROS
    rosjack_data *out = new rosjack_data[nframes];
    if(READY){
        rosjack_data **in = input_from_rosjack (nframes);
        do_overlap(in, out, nframes, apply_weights);
        
        // this type of masking requires to smooth the output since it
        // inserts many discontinuities in the time domain
        for (i = 0; i<nframes; i++){
            shift_data(out[i],past_samples,smooth_size);
            out[i] = get_mean(past_samples,smooth_size);
        }
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

void phasempf_handle_params(std::shared_ptr<rclcpp::Node> n){
    std::cout << "PhaseMPF ROS parameters: " << std::endl;
    
    n->declare_parameter("min_phase",10.0);
    n->declare_parameter("min_mag",10.0);
    n->declare_parameter("smooth_size",20);
    n->declare_parameter("MCRA_alphaS",0.95);
    n->declare_parameter("MCRA_alphaD",0.95);
    n->declare_parameter("MCRA_alphaD2",0.97);
    n->declare_parameter("MCRA_delta",0.001);
    n->declare_parameter("MCRA_L",75);
    n->declare_parameter("MPF_alphaS",0.3);
    n->declare_parameter("MPF_eta",0.3);
    n->declare_parameter("MPF_rev_gamma",0.3);
    n->declare_parameter("MPF_rev_delta",1.0);
    n->declare_parameter("out_amp",2.0);
    n->declare_parameter("noise_floor",0.001);
    n->declare_parameter("out_only_noise",false);
    n->declare_parameter("out_only_mcra",false);

    if (n->get_parameter("min_phase",min_phase)){
        RCLCPP_INFO(n->get_logger(),"Min Phase Threshold: %f",min_phase);
    }else{
        min_phase = 10.0;
        RCLCPP_WARN(n->get_logger(),"Min Phase Threshold argument not found in ROS param server, using default value (%f).",min_phase);
    }
    min_phase_diff_mean = min_phase*M_PI/180;
    
    if (n->get_parameter("min_mag",min_mag)){
        RCLCPP_INFO(n->get_logger(),"Min Mag Threshold: %f",min_mag);
    }else{
        min_mag = 10.0;
        RCLCPP_WARN(n->get_logger(),"Min Mag Threshold argument not found in ROS param server, using default value (%f).",min_mag);
    }
    
    if (n->get_parameter("smooth_size",smooth_size)){
        RCLCPP_INFO(n->get_logger(),"Smooth Filter Size: %d",smooth_size);
        if(smooth_size < 1){
          smooth_size = 20;
          RCLCPP_WARN(n->get_logger(),"Invalid Smooth Filter Size argument, using default value (%d).",smooth_size);
        }
    }else{
        smooth_size = 20;
        RCLCPP_WARN(n->get_logger(),"Smooth Filter Size argument not found in ROS param server, using default value (%d).",smooth_size);
    }
    
    if (n->get_parameter("MCRA_alphaS",MCRA_alphaS)){
        RCLCPP_INFO(n->get_logger(),"MCRA Alpha for S: %f",MCRA_alphaS);
    }else{
        MCRA_alphaS = 0.95;
        RCLCPP_WARN(n->get_logger(),"MCRA Alpha for S argument not found in ROS param server, using default value (%f).",MCRA_alphaS);
    }
    
    if (n->get_parameter("MCRA_alphaD",MCRA_alphaD)){
        RCLCPP_INFO(n->get_logger(),"MCRA Alpha for D: %f",MCRA_alphaD);
    }else{
        MCRA_alphaD = 0.95;
        RCLCPP_WARN(n->get_logger(),"MCRA Alpha for D argument not found in ROS param server, using default value (%f).",MCRA_alphaD);
    }
    
    if (n->get_parameter("MCRA_alphaD2",MCRA_alphaD2)){
        RCLCPP_INFO(n->get_logger(),"MCRA Alpha2 for D: %f",MCRA_alphaD2);
    }else{
        MCRA_alphaD2 = 0.97;
        RCLCPP_WARN(n->get_logger(),"MCRA Alpha2 for D argument not found in ROS param server, using default value (%f).",MCRA_alphaD2);
    }
    
    if (n->get_parameter("MCRA_delta",MCRA_delta)){
        RCLCPP_INFO(n->get_logger(),"MCRA Delta: %f",MCRA_delta);
    }else{
        MCRA_delta = 0.001;
        RCLCPP_WARN(n->get_logger(),"MCRA Delta argument not found in ROS param server, using default value (%f).",MCRA_delta);
    }
    
    if (n->get_parameter("MCRA_L",MCRA_L)){
        RCLCPP_INFO(n->get_logger(),"MCRA Training Windows (L): %d",MCRA_L);
    }else{
        MCRA_L = 10;
        RCLCPP_WARN(n->get_logger(),"MCRA Training Windows (L) argument not found in ROS param server, using default value (%d).",MCRA_L);
    }
    
    if (n->get_parameter("MPF_alphaS",MPF_alphaS)){
        RCLCPP_INFO(n->get_logger(),"MPF Alpha: %f",MPF_alphaS);
    }else{
        MPF_alphaS = 0.3;
        RCLCPP_WARN(n->get_logger(),"MPF Alpha argument not found in ROS param server, using default value (%f).",MPF_alphaS);
    }
    
    if (n->get_parameter("MPF_eta",MPF_eta)){
        RCLCPP_INFO(n->get_logger(),"MPF Reduction Factor: %f",MPF_eta);
    }else{
        MPF_eta = 0.3;
        RCLCPP_WARN(n->get_logger(),"MPF Reduction Factor argument not found in ROS param server, using default value (%f).",MPF_eta);
    }
    
    if (n->get_parameter("MPF_rev_gamma",MPF_rev_gamma)){
        RCLCPP_INFO(n->get_logger(),"MPF Reverberation Gamma: %f",MPF_rev_gamma);
    }else{
        MPF_rev_gamma = 0.3;
        RCLCPP_WARN(n->get_logger(),"MPF Reverberation Gamma argument not found in ROS param server, using default value (%f).",MPF_rev_gamma);
    }
    
    if (n->get_parameter("MPF_rev_delta",MPF_rev_delta)){
        RCLCPP_INFO(n->get_logger(),"MPF Reverberation Delta: %f",MPF_rev_delta);
    }else{
        MPF_rev_delta = 1.0;
        RCLCPP_WARN(n->get_logger(),"MPF Reverberation Delta argument not found in ROS param server, using default value (%f).",MPF_rev_delta);
    }
    
    if (n->get_parameter("out_amp",out_amp)){
        RCLCPP_INFO(n->get_logger(),"Output Amplification: %f",out_amp);
    }else{
        out_amp = 2.0;
        RCLCPP_WARN(n->get_logger(),"Output Amplification argument not found in ROS param server, using default value (%f).",out_amp);
    }
    
    if (n->get_parameter("noise_floor",noise_floor)){
        RCLCPP_INFO(n->get_logger(),"Noise Floor: %f",noise_floor);
    }else{
        noise_floor = 0.001;
        RCLCPP_WARN(n->get_logger(),"Noise Floor argument not found in ROS param server, using default value (%f).",noise_floor);
    }
    
    if (n->get_parameter("out_only_noise",out_only_noise)){
        RCLCPP_INFO(n->get_logger(),"Output only noise: %d",out_only_noise);
    }else{
        out_only_noise = false;
        RCLCPP_WARN(n->get_logger(),"Noise output argument not found in ROS param server, outputting filtered signal by default.");
    }
    
    if (n->get_parameter("out_only_mcra",out_only_mcra)){
        RCLCPP_INFO(n->get_logger(),"Only filter with MCRA: %d",out_only_mcra);
    }else{
        out_only_mcra = false;
        RCLCPP_WARN(n->get_logger(),"Only filter with MCRA argument not found in ROS param server, outputting with multi-channel post-filter.");
    }
}

int main (int argc, char *argv[]) {
    /* ROS initialization*/
    rclcpp::init(argc, argv);
    std::shared_ptr<rclcpp::Node> rosjack_node = rclcpp::Node::make_shared("beamform2");
    handle_params(rosjack_node);
    phasempf_handle_params(rosjack_node);
    
    rclcpp::Subscription<std_msgs::msg::Float32>::SharedPtr theta_subscriber = rosjack_node->create_subscription<std_msgs::msg::Float32>("theta", 1000, theta_roscallback);
    
    /* create JACK agent */
    if(rosjack_create (ROSJACK_READ, rosjack_node, "jackaudio", rosjack_node->get_name(), number_of_microphones, jack_callback)){
        RCLCPP_ERROR(rosjack_node->get_logger(),"JACK agent could not be created.\n");
        rclcpp::shutdown();
        exit(1);
    }
    
    std::cout << "Pre-allocating space for internal buffers." << std::endl;
    prepare_overlap_and_add(); //fft_win is assigned here, with 2 outputs
    
    x_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    x_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_time = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    y_fft = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    
    x_forward = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(x_time), reinterpret_cast<fftw_complex*>(x_fft), FFTW_FORWARD, FFTW_MEASURE);
    y_inverse = fftw_plan_dft_1d(fft_win, reinterpret_cast<fftw_complex*>(y_fft), reinterpret_cast<fftw_complex*>(y_time), FFTW_BACKWARD, FFTW_MEASURE);
    
    freqs = (double *)malloc(sizeof(double)*fft_win);
    calculate_frequency_vector(freqs,fft_win);
    
    past_samples = (double *) calloc (smooth_size,sizeof(double));
    
    delays = (double *) malloc (sizeof(double)*number_of_microphones);
    phases_aligned = (double *) malloc (sizeof(double)*number_of_microphones);
    
    out_soi = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    out_int = (std::complex<double>*) fftw_malloc(sizeof(std::complex<double>) * fft_win);
    out_soi_square = (double *)malloc(sizeof(double)*fft_win);
    out_int_square = (double *)malloc(sizeof(double)*fft_win);
    
    MCRA_S = (double *)malloc(sizeof(double)*fft_win);
    MCRA_S_prev = (double *)malloc(sizeof(double)*fft_win);
    MCRA_S_f = (double *)malloc(sizeof(double)*fft_win);
    MCRA_S_tmp = (double *)malloc(sizeof(double)*fft_win);
    MCRA_S_min = (double *)malloc(sizeof(double)*fft_win);
    MCRA_lambda_noise = (double *)malloc(sizeof(double)*fft_win);
    
    MPF_Z = (double *)malloc(sizeof(double)*fft_win);
    MPF_lambda_leak = (double *)malloc(sizeof(double)*fft_win);
    MPF_lambda_rev = (double **)malloc(sizeof(double*)*2);
    MPF_lambda_rev[0] = (double *)malloc(sizeof(double)*fft_win);
    MPF_lambda_rev[1] = (double *)malloc(sizeof(double)*fft_win);
    MPF_lambda = (double *)malloc(sizeof(double)*fft_win);
    
    //initializing vectors
    for(unsigned int i = 0; i < fft_win; i++){
        MCRA_S_prev[i] = 0.0;
        MCRA_S_tmp[i] = 0.0;
        MCRA_S_min[i] = 0.0;
        MCRA_lambda_noise[i] = 0.0;
        MPF_Z[i] = 0.0;
        MPF_lambda_leak[i] = 0.0;
        MPF_lambda_rev[0][i] = 0.0;
        MPF_lambda_rev[1][i] = 0.0;
        MPF_lambda[i] = 0.0;
    }
    
    weights.resize(number_of_microphones,fft_win);
    in_fft.resize(number_of_microphones,fft_win);
    
    update_weights(true);
    
    READY = true;
    
    RCLCPP_INFO(rosjack_node->get_logger(),"Beamform ROS node started.");
    
    /* keep running until stopped by the user */
    rclcpp::spin(rosjack_node);
    
    std::cout << "Closing theta subscriber." << std::endl;
    theta_subscriber.reset();
    
    exit(0);
}
