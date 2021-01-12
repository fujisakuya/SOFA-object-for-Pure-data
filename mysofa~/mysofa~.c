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

#define MAX_BLOCKSIZE 8192
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
    t_float flag;
    t_symbol *filenameArg;
    t_float filter_length;
    t_int err;
    
    
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
    t_inlet *x_in5;
    t_outlet *x_out1;
    t_outlet *x_out2;
    
    t_int M;    /* number of filters stored in sofa file */
    t_int N;    /* numer of sample points per stored filter */
    t_int R;    /* number of receivers */
    
    struct MYSOFA_EASY *sofa;
    char filename[1000];
    fftwf_complex *in_signal_fftwf,*out_signal_fftwf;
    fftwf_complex *in_left_ir_fftwf, *out_left_ir_fftwf,*in_right_ir_fftwf,*out_right_ir_fftwf;
    fftwf_complex *in_left_final, *out_left_final, *in_right_final, *out_right_final;
    
    fftwf_plan plan1, plan2, plan3, plan4,plan5;
    
    
} t_mysofa_tilde;

// MAX_BLOCKSIZE 8192
 //blockScale = MAX_BLOCKSIZE / blocksize;
 //scaledBlocksize = blocksize * blockScale;
// blocksizeDelta = MAX_BLOCKSIZE -1 - scaledBlocksize

// t_symbol *filenameArg


/*
 // This function is copied from https://github.com/sofacoustics/SOFAlizer-for-pd/blob/master/SOFAlizer~/SOFAlizer~.c
 static void mysofa_tilde_open (t_mysofa_tilde *x, t_symbol *filenameArg) {
 strcpy(x->filename, filenameArg->s_name);//->s_name);
 if (x->filename[0] == '\0') {
 error("No SOFA file specified");
 return;
 } else {
 post("SOFA file %s will be opened..", x->filename);
 }
 int filter_length, err;
 x->sofa = mysofa_open(x->filename, 48000, &filter_length, &err); // change this
 if(x->sofa == NULL) {
 error("SOFA file %s couldn't be opened.", x->filename);
 return;
 }
 x->M = x->sofa->hrtf->M;    // number of filters stored in sofa file
 x->N = x->sofa->hrtf->N;    // numer of sample points per stored filter
 x->R = x->sofa->hrtf->R;    //number of receivers
 
 post("Metadata: %s %u HRTFs, %u samples @ %u Hz. %s, %s, %s.\n",
 x->filename,(unsigned int)x->M, (unsigned int)x->N,
 (unsigned int)(x->sofa->hrtf->DataSamplingRate.values[0]),
 mysofa_getAttribute(x->sofa->hrtf->attributes, "DatabaseName"),
 mysofa_getAttribute(x->sofa->hrtf->attributes, "Title"),
 mysofa_getAttribute(x->sofa->hrtf->attributes, "ListenerShortName"));
 
 // Calculate distance of first measurement
 int i;
 for (i = 0; i < 3; i++) {    // read in first three coordinates
 x->values[i] = x->sofa->hrtf->SourcePosition.values[i];
 }
 x->distance = sqrtf(powf(x->values[0], 2.f) + powf(x->values[1], 2.f) + powf(x->values[2], 2.f)); // calculate distance
 post("Distance is %f m. \n", x->distance);
 
 // Check for IR delays
 float delay = 0;
 int warn = 0;
 for (i = 0; i < x->R; i++) {
 delay = x->sofa->hrtf->DataDelay.values[i];
 if (delay != 0.0) {
 warn = 1;
 }
 }
 if (warn == 1) {
 error(" Warning: This SOFA file will be incorrectly processed besause of non zero IR delays!"); // To change so we can use th delays too
 }
 }
 
 */



t_int *mysofa_tilde_perform(t_int *w) {
    t_mysofa_tilde *x = (t_mysofa_tilde *)(w[1]);
    t_sample  *in =    (t_sample *)(w[2]);
    t_sample  *out1 =    (t_sample *)(w[3]);
    t_sample  *out2 =    (t_sample *)(w[4]);
    int       n =           (int)(w[5]);
    
    float blocksize = n;
    int i = 0;
    float realD,imagD,realS,imagS;
    float mux = 1.0/blocksize;
    
    
    x->leftDelay = 0;
    x->rightDelay = 0;


    
    
    
    x->values[0] = x->azimuth;//x->azimuth;
    x->values[1] = x->elevation;//x->elevation;
    x->values[2] = x->distance;//x->distance;
    
    mysofa_s2c(x->values);
    
    
    /*
    post("%f, ",x->values[0]);
    post("%f, ",x->values[1]);
    post("%f\n",x->values[2]);
     */
    
    x->x = x->values[0];
    x->y = x->values[1];
    x->z = x->values[2];
    
     mysofa_getfilter_float(x->sofa,x->x,x->y,x->z,x->leftIR,x->rightIR,&x->leftDelay,&x->rightDelay);
    
    

    
     
     //signal
     for( i = 0; i<n ; i++){
     x->in_signal_fftwf[i][0] = in[i];
     x->in_signal_fftwf[i][1] = 0.0;
     }
     
     
     for( i = n; i<x->fftsize ; i++){
     x->in_signal_fftwf[i][0] = 0.0;
     x->in_signal_fftwf[i][1] = 0.0;
     }
    
    
    
    //leftIR
    for(i =0;i<x->filter_length;i++){
        x->in_left_ir_fftwf[i][0] = x->leftIR[i];
        x->in_left_ir_fftwf[i][1] = 0.0;
    }
    
    for(i =x->filter_length; i<x->fftsize; i++){
        x->in_left_ir_fftwf[i][0] = 0.0;
        x->in_left_ir_fftwf[i][1] = 0.0;
    }
    
    
    //right IR
    for(i =0;i<x->filter_length;i++){
        x->in_right_ir_fftwf[i][0] = x->rightIR[i];
        x->in_right_ir_fftwf[i][1] = 0.0;
    }
    
    for(i =x->filter_length; i<x->fftsize; i++){
        x->in_right_ir_fftwf[i][0] = 0.0;
        x->in_right_ir_fftwf[i][1] = 0.0;
    }
    
    
    
    
    
    
    
    /*
     //left IR
     for(i =0;i<x->filter_length;i++){
     scaledBlocksize = (x->filter_length-i)*blockScale;
     x->in_left_ir_fftwf[i][0] = x->leftIR[i] * x->crossCoef[scaledBlocksize] + x->previousImpulse[i] * x->crossCoef[MAX_BLOCKSIZE - scaledBlocksize];
     x->in_left_ir_fftwf[i][1] = 0.0;
     x->previousImpulse[i] = x->in_left_ir_fftwf[i][0];
     }
     
     for(i =x->filter_length; i<x->fftsize; i++){
     x->in_left_ir_fftwf[i][0] = 0.0;
     x->in_left_ir_fftwf[i][1] = 0.0;
     }
     
     
     //right IR
     for(i =0;i<x->filter_length;i++){
     scaledBlocksize = (x->filter_length-i)*blockScale;
     x->in_right_ir_fftwf[i][0] = x->leftIR[i] * x->crossCoef[scaledBlocksize] + x->previousImpulse[i] * x->crossCoef[MAX_BLOCKSIZE - scaledBlocksize];
     x->in_right_ir_fftwf[i][1] = 0.0;
     x->previousImpulse[i] = x->in_right_ir_fftwf[i][0];
     }
     
     for(i=x->filter_length; i<x->fftsize; i++){
     x->in_left_ir_fftwf[i][0] = 0.0;
     x->in_left_ir_fftwf[i][1] = 0.0;
     }
     */

    

     
     fftwf_execute(x->plan1);
     fftwf_execute(x->plan2);
     fftwf_execute(x->plan3);
    

    
    
    
     for( i = 0; i < x->convsize; i++ ){
     realD = x->out_signal_fftwf[i][0];
     imagD = x->out_signal_fftwf[i][1];
     realS = x->out_left_ir_fftwf[i][0];
     imagS = x->out_left_ir_fftwf[i][1];
     x->out_left_final[i][0] = (realD * realS - imagD * imagS)*mux;
     x->out_left_final[i][1] = (realD * imagS + imagD * realS)*mux;
     }
     
     
     
     for( i = 0; i < x->convsize; i++ ){
     realD = x->out_signal_fftwf[i][0];
     imagD = x->out_signal_fftwf[i][1];
     realS = x->out_right_ir_fftwf[i][0];
     imagS = x->out_right_ir_fftwf[i][1];
     x->out_right_final[i][0] = (realD * realS - imagD * imagS)*mux;
     x->out_right_final[i][1] = (realD * imagS + imagD * realS)*mux;
     }
    

    
  /*
    for(i=0;i<n;i++){
        post("right %f\n",x->out_right_final[i][0]);
    }
    
    for(i=0;i<n;i++){
        post("left %f\n", x->out_left_final[i][0] );
    }

     */
    
     fftwf_execute(x->plan4);
     fftwf_execute(x->plan5);
    
    
    /*
    for(i=0; i<n; i++){
        post("x->in_right_final[i][0] : %lf\n",x->in_right_final[i][0]);
    }
    
    for(i=0; i<n; i++){
        post("x->in_right_final[i][1] : %lf\n",x->in_right_final[i][1]);
    }
     */
    
    
    
    
     for(i=x->convsize;i<x->fftsize;i++){
         x->in_left_final[i][0] = 0.0;
         x->in_right_final[i][0] = 0.0;
     }

    
    
    

    /*
     if(x->flag == 0){
     
         for(i=0;i<n;i++){
             out1[i] = x->in_right_final[i][0];
             out2[i] = x->in_left_final[i][0];
         }
         
         int j = 0;
         for(i=n; i<x->convsize; i++){
             x->buffer_right[j] = x->in_right_final[i][0];
             x->buffer_left[j] = x->in_left_final[i][0];
             j++;
         }
     
         for(i=0; i<n; i++){
             post("out1[i] : %lf\n",out1[i]);
         }
         for(i=0; i<n; i++){
             post("out2[i] : %lf\n",out2[i]);
         }

         
            x->flag = 1;
         
        }
    */

    
    for(i = x->filter_length; i<x->convsize; i++){
            x->buffer_right[i] = 0.0;
            x->buffer_left[i] = 0.0;
    }
    
    
    //in_left_finalが1を超えていた
    
    for(i = 0;i<x->convsize;i++){
        x->buffer_right[i] = x->buffer_right[i] + x->in_right_final[i][0];
        x->buffer_left[i] = x->buffer_left[i] + (x->in_left_final[i][0] * mux);
    }

    
     
    for(i=0;i<n;i++){
        out1[i] = x->buffer_right[i];
        out2[i] = x->buffer_left[i];
    }
    

         
    int j = 0;
    for(i=n; i<x->convsize; i++){
        x->buffer_right[j] = x->in_right_final[i][0];
        x->buffer_left[j] = x->in_left_final[i][0];
        j++;
    }
    
    
         
    /*
    for(i=0; i<n; i++){
        post("out1[i] : %lf\n",out1[i]);
    }
    for(i=0; i<n; i++){
        post("out2[i] : %lf\n",out2[i]);
    }
    */
         
    
   //  post("%f\n",x->flag);
    
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
    int size[8] = {128, 256, 512, 1024, 2048, 4096, 8192};
    
    // x->sofa = mysofa_open("/Users/fujisawasakuwataru/Downloads/Pd/mysofa/RIEC_hrtf_all/RIEC_hrir_subject_001.sofa", 48000, &filter_length, &err);
    
    x->sofa = mysofa_open(x->filenameArg->s_name, 48000, &filter_length, &err);
    //mysofa_tilde_open(x, x->filenameArg);
    x->filter_length = filter_length;
    x->convsize = x->filter_length + sp[0]->s_n;
    x->flag = 0;
    
    while(1){
        if(i==8){
            post("blocksize is too large");
            break;
        }
        if(x->fftsize<=size[i]){
            x->fftsize = size[i];
            break;
        }
        i++;
    }
    
    
    
    
    if(!x->sofa) {
        post("Error opening the SOFA file");
    }
    post("%d",err);
    
    
    x->in_signal_fftwf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->out_signal_fftwf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->plan1 = fftwf_plan_dft_1d(x->fftsize, x->in_signal_fftwf, x->out_signal_fftwf, FFTW_FORWARD, FFTW_ESTIMATE);
    
    x->in_left_ir_fftwf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->out_left_ir_fftwf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->plan2 = fftwf_plan_dft_1d(x->fftsize, x->in_left_ir_fftwf, x->out_left_ir_fftwf, FFTW_FORWARD, FFTW_ESTIMATE);
    
    x->in_right_ir_fftwf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->out_right_ir_fftwf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->plan3 = fftwf_plan_dft_1d(x->fftsize, x->in_left_ir_fftwf, x->out_left_ir_fftwf, FFTW_FORWARD, FFTW_ESTIMATE);
    
    
    x->in_left_final = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->out_left_final= (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->plan4 = fftwf_plan_dft_1d(x->fftsize,x->out_left_final,x->in_left_final, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    x->in_right_final = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->out_right_final= (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * x->fftsize);
    x->plan5 = fftwf_plan_dft_1d(x->fftsize,x->out_right_final,x->in_right_final, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    
    for(i = 0; i<x->fftsize; i++){
        x->previousImpulse[i] = 0;
    }
    
    for (i = 0; i < MAX_BLOCKSIZE; i++) {
        x->crossCoef[i] = 1.0 * i / (MAX_BLOCKSIZE-1.0);
    }
    
    for(i = 0; i<x->fftsize; i++){
        x->buffer_right[i] = 0.0;
        x->buffer_left[i] = 0.0;
    }
    
}






void mysofa_tilde_free(t_mysofa_tilde *x) {
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    inlet_free(x->x_in4);
    inlet_free(x->x_in5);
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

    
    mysofa_close(x->sofa);
}




void *mysofa_tilde_new(void) {
    t_mysofa_tilde *x = (t_mysofa_tilde *)pd_new(mysofa_tilde_class);
    //  This is called every time an object is created.
    x->x_in2 = symbolinlet_new(&x->x_obj, &x->filenameArg);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->azimuth);
    x->x_in4 = floatinlet_new(&x->x_obj, &x->elevation);
    x->x_in5 = floatinlet_new(&x->x_obj, &x->distance);
    
    
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    
    
    //Additional signal outlet. This is made in same way.
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
    CLASS_MAINSIGNALIN(mysofa_tilde_class, t_mysofa_tilde, f);
    
    // class_addmethod(mysofa_tilde_class, (t_method)mysofa_tilde_open, gensym("open"), A_DEFSYMBOL, 0); /* method to open new sofa file */
}

//Name object "mysofa~"
//Instructor。Means that it is initialized every time
//Destructor。 we need destructor because we created an additional inlet.
//Allocate static memory.
//Graphic settings for objects. Make it the default.
//Use A_GIMME if you need more than 6 arguments for the symbolic object. 0 means End.
//By attaching a dsp selector, it becomes a signal class, and by having A_CANT, invalid dsp cannot be sent.
//In CLASS_MAINSIGNALIN, the inlet that exists from the beginning becomes the signal inlet.
