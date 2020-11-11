#include <math.h>
#include <string.h>
#include "fftw3.h"
#include "m_pd.h"

static t_class *conv_tilde_class;
//two inlets
//First inlet is original signal
//Second inlet is IR signal
#define MAX_BLOCKSIZE 8192
#define MAX_N_POINTS 3000
#define fftsize 1024

typedef struct _conv_tilde{
    t_object x_obj;  //

    t_sample x_conv_before; //original data before fft
    t_sample x_ir_before; //ir data before fft
    t_sample x_conv_after;//original data after fft
    t_sample x_ir_after; //ir data after fft
    t_sample f;
    t_sample previousImpulse[MAX_N_POINTS];
    t_float crossCoef[MAX_BLOCKSIZE];
    t_int bufferPin;
    t_int nPts;//It is used to convert a float type message into a signal when it enters the signal inlet.

    fftw_complex *in_signal_fftw,*out_signal_fftw, *in_ir_fftw, *out_ir_fftw, *in_final, *out_final;
    fftw_plan plan1, plan2, plan3;

    t_inlet *x_in2;
    t_outlet *x_out1;
    t_outlet *x_out2;


}t_conv_tilde;


       // MAX_BLOCKSIZE 8192
       //blockScale = MAX_BLOCKSIZE / blocksize;
       //scaledBlocksize = blocksize * blockScale;
       //blocksizeDelta = MAX_BLOCKSIZE -1 - scaledBlocksize


t_int *conv_tilde_perform(t_int *w)
{

    t_conv_tilde *x = (t_conv_tilde *)(w[1]);
    t_sample  *in1 =    (t_sample *)(w[2]);
    t_sample  *in2 =    (t_sample *)(w[3]);
    t_sample  *out1 =    (t_sample *)(w[4]);
    t_sample  *out2 =    (t_sample *)(w[5]);
    int       n =           (int)(w[6]);


    int blocksize = n;
    int blockScale  = MAX_BLOCKSIZE / blocksize;
    unsigned scaledBlocksize;
    int i = 0;
   // int fftsize = 44100;
    double realD,imagD,realS,imagS;
    double mux = 1.0/blocksize;
    
    
    for (int i = 0; i < MAX_BLOCKSIZE; i++) {
      x->crossCoef[i] = 1.0 * i / (MAX_BLOCKSIZE-1.0);
    }
    

    //for (int i = 0; i < MAX_BLOCKSIZE; i++) {
      //  x->crossCoef[i] = 1.0 * i / (MAX_BLOCKSIZE-1.0);
        //x->crossCoef[i] = cos((MyPI/2) * ((float)i/(MAX_BLOCKSIZE-1.0))); ///cosine-sine cross-fading
    //}



     for( i = 0; i<n ; i++){
        x->in_signal_fftw[i][0] = in1[i];
        x->in_signal_fftw[i][1] = 0.0;
    }



  //  for( i = n; i<(n*2) ; i++){
   //     x->in_signal_fftw[i][0] = 0.0;
   //     x->in_signal_fftw[i][1] = 0.0;
  //  }
    
    

     for(i =0;i<n;i++){
        scaledBlocksize = (n-i)*blockScale;
        x->in_ir_fftw[i][0] = in2[i] * x->crossCoef[scaledBlocksize] + x->previousImpulse[i] * x->crossCoef[MAX_BLOCKSIZE - scaledBlocksize];
        x->in_ir_fftw[i][1] = 0.0;
        x->previousImpulse[i] = x->in_ir_fftw[i][0];
         }
    


 //   for( i = n; i<(n*2) ; i++){
  //      x->in_ir_fftw[i][0] = 0.0;
  //      x->in_ir_fftw[i][1] = 0.0;
   // }



    fftw_execute(x->plan1);
    fftw_execute(x->plan2);
    
    


    /*

    for( i = 0; i < n; i++ ){
        realD = x->out_signal_fftw[i][0];
        imagD = x->out_signal_fftw[i][1];
        realS = x->out_ir_fftw[i][0];
        imagS = x->out_ir_fftw[i][1];
        x->out_final[i][0] = (realD * realS - imagD * imagS)*mux;
        x->out_final[i][1] = (realD * imagS + imagD * realS)*mux;
    }
*/


     for( i = 0; i < (n*2); i++ ){
     realD = x->out_signal_fftw[i][0];
     imagD = x->out_signal_fftw[i][1];
     realS = x->out_ir_fftw[i][0];
     imagS = x->out_ir_fftw[i][1];
     x->out_final[i][0] = (realD * realS - imagD * imagS)*mux;
     x->out_final[i][1] = (realD * imagS + imagD * realS)*mux;
     }


    fftw_execute(x->plan3);


    for(i=0; i<n; i++){
        out1[i] = x->in_final[i][0] * (0.5 - 0.5 * cos(2 * M_PI * ((float)i/n)));
    }

//*(0.5 - 0.5 * cos(2 * M_PI * ((float)i/n)))
    

    for(i=0; i<n; i++){
        out2[i] = x->in_final[i][0] * (0.5 - 0.5 * sin(2 * M_PI * ((float)i/n)));
    }

    return (w+7);

}



void conv_tilde_dsp(t_conv_tilde *x, t_signal **sp)
{
    dsp_add(conv_tilde_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
    /*
    x->in_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n));
    x->out_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *(sp[0]->s_n));
    x->plan1 = fftw_plan_dft_1d((sp[0]->s_n), x->in_signal_fftw, x->out_signal_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

    x->in_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n));
    x->out_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n));
    x->plan2 = fftw_plan_dft_1d((sp[0]->s_n), x->in_ir_fftw, x->out_ir_fftw, FFTW_FORWARD, FFTW_ESTIMATE);


    x->in_final = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n));
    x->out_final= (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n));
    x->plan3 = fftw_plan_dft_1d((sp[0]->s_n),x->out_final,x->in_final, FFTW_BACKWARD, FFTW_ESTIMATE);
    */





    x->in_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n*2));
    x->out_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *(sp[0]->s_n*2));
    x->plan1 = fftw_plan_dft_1d((sp[0]->s_n*2), x->in_signal_fftw, x->out_signal_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

    x->in_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n*2));
    x->out_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n*2));
    x->plan2 = fftw_plan_dft_1d((sp[0]->s_n*2), x->in_ir_fftw, x->out_ir_fftw, FFTW_FORWARD, FFTW_ESTIMATE);


    x->in_final = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n*2));
    x->out_final= (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (sp[0]->s_n*2));
    x->plan3 = fftw_plan_dft_1d((sp[0]->s_n*2),x->out_final,x->in_final, FFTW_BACKWARD, FFTW_ESTIMATE);


    for(int i = 0; i<sp[0]->s_n; i++){
        x->previousImpulse[i] = 0;
    }

}




void conv_tilde_free(t_conv_tilde *x)
{
    inlet_free(x->x_in2);
    outlet_free(x->x_out1);
    outlet_free(x->x_out2);

    fftw_destroy_plan(x->plan1);
    fftw_free(x->in_signal_fftw);
    fftw_free(x->out_signal_fftw);
    fftw_destroy_plan(x->plan2);
    fftw_free(x->in_ir_fftw);
    fftw_free(x->out_ir_fftw);
    fftw_destroy_plan(x->plan3);
    fftw_free(x->in_final);
    fftw_free(x->out_final);
}



void *conv_tilde_new(t_floatarg f)
{
    t_conv_tilde *x = (t_conv_tilde *)pd_new(conv_tilde_class);
    // This is called every time an object is created.
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    //Additional signal inlet. The argument is a reference to the lookup tail symbolic selector "signal".
    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);
    //Additional signal outlet. This is made in same way.
    x->nPts = f;
    x->bufferPin = 0;



    return (void *)x;
}





void conv_tilde_setup(void){


    conv_tilde_class=class_new(gensym("conv~"),
                               (t_newmethod)conv_tilde_new,
                               (t_method)conv_tilde_free,
                               sizeof(t_conv_tilde),
                               CLASS_DEFAULT,A_DEFFLOAT,0);

    class_addmethod(conv_tilde_class,
                    (t_method)conv_tilde_dsp,gensym("dsp"),A_CANT,0);
    CLASS_MAINSIGNALIN(conv_tilde_class, t_conv_tilde, f);

}
//Name object "conv~"
//Instructor。Means that it is initialized every time
//Destructor。 we need destructor because we created an additional inlet.
//Allocate static memory.
//Graphic settings for objects. Make it the default.
//Use A_GIMME if you need more than 6 arguments for the symbolic object. 0 means End.
//By attaching a dsp selector, it becomes a signal class, and by having A_CANT, invalid dsp cannot be sent.
//In CLASS_MAINSIGNALIN, the inlet that exists from the beginning becomes the signal inlet.
