#include <math.h>
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "m_pd.h"
#include <mysofa.h>
#include <fftw3.h>

static t_class *mysofa_tilde_class;

#define MAX_BLOCKSIZE 16384
#define MAX_N_POINTS 3000


typedef struct _mysofa_tilde {
    t_object x_obj;  //
    t_float rightIR[MAX_BLOCKSIZE];
    t_float leftIR[MAX_BLOCKSIZE];
    t_float f;
    t_float azimuth;
    t_float elevation;
    t_float distance;
    t_float values[3];
    t_float x,y,z;
    t_float leftDelay;
    t_float rightDelay;
    t_float fftsizeArg;
    t_float fftsize;
    t_float convsize;
    t_float nbins;
    t_float flag;
    t_symbol *filenameArg;
    t_float filter_length;
    t_float err;
    t_float delaysize;
    t_float sr;
    
    t_sample x_conv_before; //original data before fft
    t_sample x_ir_before; //ir data before fft
    t_sample x_conv_after;//original data after fft
    t_sample x_ir_after; //ir data after fft
    t_sample buffer_right[MAX_BLOCKSIZE];
    t_sample buffer_left[MAX_BLOCKSIZE];
    t_sample previousImpulse[MAX_N_POINTS];
    t_float crossCoef[MAX_BLOCKSIZE];
    t_int bufferPin;
    t_int nPts;//It is used to convert a float type message into a signal when it enters the signal inlet.
    
    
    t_inlet *x_in2;
    t_inlet *x_in3;
    t_inlet *x_in4;
    t_outlet *x_out1;
    t_outlet *x_out2;
    
    t_int M;    /* number of filters stored in sofa file */
    t_int N;    /* numer of sample points per stored filter */
    t_int R;    /* number of receivers */
    
    struct MYSOFA_EASY *sofa;
    char filename[1000];
    float *in_signal_fftwf;
    float *in_left_ir_fftwf, *in_right_ir_fftwf;
    float *out_left_final, *out_right_final;
    float *out_left_final_p, *out_right_final_p;

    
    fftwf_complex *out_signal_fftwf;
    fftwf_complex *out_left_ir_fftwf,*out_right_ir_fftwf;
    fftwf_complex *in_left_final, *in_right_final;

    
    fftwf_plan plan1, plan2, plan3, plan4,plan5;
    
    
} t_mysofa_tilde;

// MAX_BLOCKSIZE 8192
//blockScale = MAX_BLOCKSIZE / blocksize;
//scaledBlocksize = blocksize * blockScale;
// blocksizeDelta = MAX_BLOCKSIZE -1 - scaledBlocksize






t_int *mysofa_tilde_perform(t_int *w) {
    t_mysofa_tilde *x = (t_mysofa_tilde *)(w[1]);
    t_sample  *in =    (t_sample *)(w[2]);
    t_sample  *out1 =    (t_sample *)(w[3]);
    t_sample  *out2 =    (t_sample *)(w[4]);
    int       n =           (int)(w[5]);
    
    
    if(x->err!=0){
        error("The sofa file is not loaded");
    }

    
    
    else{
        int i = 0;
        int j = 0;
        float realD,imagD,realS,imagS;
        float mux = 1.0/x->fftsize;
        x->nbins = x->fftsize/2 + 1;
        
        
        x->values[0] = x->azimuth;//x->azimuth;
        x->values[1] = x->elevation;//x->elevation;
        x->values[2] = x->distance;//x->distance;
        
        mysofa_s2c(x->values);
        
        
        //get leftIR and rightIR
        if(x->x != x->values[0] || x->y != x->values[1] || x->z != x->values[2]){
            
            x->x = x->values[0];
            x->y = x->values[1];
            x->z = x->values[2];
            
            mysofa_getfilter_float(x->sofa,x->x,x->y,x->z,x->leftIR,x->rightIR,&x->leftDelay,&x->rightDelay);
            x->delaysize = x->rightDelay + x->leftDelay;
        }
        
        
        
        //signal
        for(i = 0; i<x->fftsize ; i++){
            if(i<n){
                x->in_signal_fftwf[i] = in[i];
            }
            else{
                x->in_signal_fftwf[i] = 0.0;
            }
        }
        
        
        
        //leftIR
        for(i =0;i<x->fftsize;i++){
            if(i<x->filter_length){
                x->in_left_ir_fftwf[i] = x->leftIR[i];
            }
            else{
            x->in_left_ir_fftwf[i] = 0.0;
            }
            
        }
        
        
        //right IR
        for(i =0;i<x->fftsize;i++){
            if(i<x->filter_length){
                x->in_right_ir_fftwf[i] = x->rightIR[i];
            }
            else{
                x->in_right_ir_fftwf[i] = 0.0;
            }
            
        }

        
        //fft
        fftwf_execute(x->plan1);
        fftwf_execute(x->plan2);
        fftwf_execute(x->plan3);
        
        
        
        for( i = 0; i < x->nbins; i++ ){
            //left convolution
            realD = x->out_signal_fftwf[i][0];
            imagD = x->out_signal_fftwf[i][0];
            realS = x->out_left_ir_fftwf[i][0];
            imagS = x->out_left_ir_fftwf[i][1];
            x->in_left_final[i][0] = (realD * realS - imagD * imagS)*mux;
            x->in_left_final[i][1] = (realD * imagS + imagD * realS)*mux;
            
            //right convolution
            realD = x->out_signal_fftwf[i][0];
            imagD = x->out_signal_fftwf[i][0];
            realS = x->out_right_ir_fftwf[i][0];
            imagS = x->out_right_ir_fftwf[i][1];
            x->in_right_final[i][0] = (realD * realS - imagD * imagS)*mux;
            x->in_right_final[i][1] = (realD * imagS + imagD * realS)*mux;
        }
        
        
        
        //inverse fft
        fftwf_execute(x->plan4);
        fftwf_execute(x->plan5);
        
        
        
        //delay panning
        x->delaysize = x->fftsize + x->leftDelay + x->rightDelay;
        
        j=0;
        for(i = 0; i<x->delaysize; i++){
            if(i<x->leftDelay || i > (x->leftDelay + x->fftsize)){
                x->out_left_final_p[i] = 0.0;
            }
            
            else{
                x->out_left_final_p[i] = x->out_left_final[j];
                j++;
            }
        }
        
        j=0;
        for(i = 0; i<x->delaysize; i++){
            if(i<x->rightDelay || i > (x->rightDelay + x->fftsize)){
                x->out_right_final_p[i] = 0.0;
            }
            else{
                x->out_right_final_p[i] = x->out_right_final[j];
                j++;
            }
        }
        
        
        
        //output and storage buffer
        j = 0;
        for(i=0;i<x->delaysize;i++){
            if(i<n){
                out1[i] = x->buffer_left[i] + x->out_left_final_p[i];
                out2[i] = x->buffer_right[i] + x->out_right_final_p[i];
            }
            
            else{
                x->buffer_left[j] = x->buffer_left[i] + x->out_left_final_p[i];
                x->buffer_right[j] = x->buffer_right[i] + x->out_right_final_p[i];
                j++;
            }
        }
        
    }
    //out1 is left
    //out2 is right
    
    return (w+6);
    
}





void mysofa_tilde_dsp(t_mysofa_tilde *x, t_signal **sp) {
    dsp_add(mysofa_tilde_perform,
            5,
            x,
            sp[0]->s_vec, //in_signal
            sp[1]->s_vec, //out_signal_r
            sp[2]->s_vec, //out_signal_l
            sp[0]->s_n);
    
    int filter_length, err;
    int i=0;
    int size[8] = {128, 256, 512, 1024, 2048, 4096, 8192,16384};
    
    x->err = 100.0;
    x->sr = sp[0]->s_sr;
    x->sofa = mysofa_open(x->filenameArg->s_name, x->sr, &filter_length, &err);
    //mysofa_tilde_open(x, x->filenameArg);
    x->filter_length = filter_length;
    x->convsize = x->filter_length + sp[0]->s_n - 1;
    x->err = err;
    x->azimuth = 0;
    x->elevation = 0;
    x->distance = 60;


    if(x->err != 0){
        error( "The file could not be read.");
    }
    
    
    else{
    
        while(1){
            if(i==8){
                post("blocksize is too large");
                break;
            }
            if(x->convsize <= size[i]){
                x->fftsize = size[i];
                break;
            }
            i++;
        }
        
        post("filter_length : %f",x->filter_length);
        post("convsize : %f",x->convsize);
        post("fftsize : %f",x->fftsize);
        
        
        if(!x->sofa) {
            post("Error opening the SOFA file");
        }

        x->in_signal_fftwf = fftwf_alloc_real(x->fftsize);
        x->out_signal_fftwf = fftwf_alloc_complex(x->fftsize);
        x->plan1 = fftwf_plan_dft_r2c_1d(x->fftsize, x->in_signal_fftwf, x->out_signal_fftwf,FFTW_ESTIMATE);
        
        x->in_left_ir_fftwf = fftwf_alloc_real(x->fftsize);
        x->out_left_ir_fftwf = fftwf_alloc_complex(x->fftsize);
        x->plan2 = fftwf_plan_dft_r2c_1d(x->fftsize, x->in_left_ir_fftwf, x->out_left_ir_fftwf, FFTW_ESTIMATE);
        
        x->in_right_ir_fftwf = fftwf_alloc_real(x->fftsize);
        x->out_right_ir_fftwf = fftwf_alloc_complex(x->fftsize);
        x->plan3 = fftwf_plan_dft_r2c_1d(x->fftsize, x->in_right_ir_fftwf, x->out_right_ir_fftwf, FFTW_ESTIMATE);
        
        
        x->in_left_final = fftwf_alloc_complex(x->fftsize);
        x->out_left_final= fftwf_alloc_real(x->fftsize);
        x->plan4 = fftwf_plan_dft_c2r_1d(x->fftsize,x->in_left_final,x->out_left_final, FFTW_ESTIMATE);
        
        x->in_right_final = fftwf_alloc_complex(x->fftsize);
        x->out_right_final= fftwf_alloc_real(x->fftsize);
        x->plan5 = fftwf_plan_dft_c2r_1d(x->fftsize,x->in_right_final,x->out_right_final, FFTW_ESTIMATE);
        
        x->out_right_final_p = (float*)malloc(sizeof(float)*x->fftsize + 10000);
        x->out_left_final_p = (float*)malloc(sizeof(float)*x->fftsize + 10000);

        
        for(i = 0; i<x->fftsize + 10000; i++){
            x->buffer_right[i] = 0.0;
            x->buffer_left[i] = 0.0;
        }
        
    }
    
}



void mysofa_tilde_symbol(t_mysofa_tilde *x, t_symbol *s){
    x->filenameArg = s;
}





void mysofa_tilde_free(t_mysofa_tilde *x) {
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    inlet_free(x->x_in4);
    outlet_free(x->x_out1);
    outlet_free(x->x_out2);
    
    fftwf_destroy_plan(x->plan1);
    fftwf_free(x->in_signal_fftwf);
    fftwf_free(x->out_signal_fftwf);
    fftwf_destroy_plan(x->plan2);
    fftwf_destroy_plan(x->plan3);
    fftwf_free(x->in_left_ir_fftwf);
    fftwf_free(x->in_right_ir_fftwf);
    fftwf_free(x->out_left_ir_fftwf);
    fftwf_free(x->out_right_ir_fftwf);
    fftwf_destroy_plan(x->plan4);
    fftwf_destroy_plan(x->plan5);
    fftwf_free(x->in_left_final);
    fftwf_free(x->in_right_final);
    fftwf_free(x->out_left_final);
    fftwf_free(x->out_right_final);
    
    free(x->out_left_final_p);
    free(x->out_right_final_p);
    
    mysofa_close(x->sofa);
}

float *in_signal_fftwf;
float *in_left_ir_fftwf, *in_right_ir_fftwf;
float *out_left_final, *out_right_final;
float *out_left_final_p, *out_right_final_p;


fftwf_complex *out_signal_fftwf;
fftwf_complex *out_left_ir_fftwf,*out_right_ir_fftwf;
fftwf_complex *in_left_final, *in_right_final;




void *mysofa_tilde_new(void) {
    t_mysofa_tilde *x = (t_mysofa_tilde *)pd_new(mysofa_tilde_class);
    x->x_in2 = floatinlet_new(&x->x_obj, &x->azimuth);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->elevation);
    x->x_in4 = floatinlet_new(&x->x_obj, &x->distance);
    
    
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    x->err = 100.0;
    post("err : %d",x->err);
    

    return (void *)x;
}





void mysofa_tilde_setup(void) {
    
    mysofa_tilde_class=class_new(gensym("mysofa~"),
                                 (t_newmethod)mysofa_tilde_new,
                                 (t_method)mysofa_tilde_free,
                                 sizeof(t_mysofa_tilde),
                                 CLASS_DEFAULT,
                                 0);
    
    class_addmethod(mysofa_tilde_class,
                    (t_method)mysofa_tilde_dsp,
                    gensym("dsp"),A_CANT,0);
    
    class_addsymbol(mysofa_tilde_class, mysofa_tilde_symbol);
    CLASS_MAINSIGNALIN(mysofa_tilde_class, t_mysofa_tilde, f);
    
}
