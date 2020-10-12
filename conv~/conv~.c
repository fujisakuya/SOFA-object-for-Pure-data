#include <math.h>
#include <string.h>
#include "fftw3.h"
#include "m_pd.h"

static t_class *conv_tilde_class;
//two inlets
//First inlet is original signal
//Second inlet is IR signal


typedef struct _conv_tilde{
    t_object x_obj;  //
    
    t_sample x_conv_before; //original data before fft
    t_sample x_ir_before; //ir data before fft
    t_sample x_conv_after;//original data after fft
    t_sample x_ir_after; //ir data after fft
    t_float f;//It is used to convert a float type message into a signal when it enters the signal inlet.
    
    fftw_complex *in_signal_fftw,*out_signal_fftw, *in_ir_fftw, *out_ir_fftw, *in_final, *out_final;
    fftw_plan plan1, plan2, plan3;
    
    t_inlet *x_in2;
    t_outlet *x_out;
    
    int fftsize;
    
}t_conv_tilde;






t_int *conv_tilde_perform(t_int *w)
{
    
    t_conv_tilde *x = (t_conv_tilde *)(w[1]);
    t_sample  *in1 =    (t_sample *)(w[2]);
    t_sample  *in2 =    (t_sample *)(w[3]);
    t_sample  *out =    (t_sample *)(w[4]);
    int          n =           (int)(w[5]);
    
    int i;
   // int fftsize = 44100;
    double realD,imagD,realS,imagS;
    double mux = 1.0/n;


    
    

    
    for(i=0; i<n; i++){
        x->in_signal_fftw[i][0] = in1[i];
        x->in_signal_fftw[i][1] = 0.0;
        
        x->in_ir_fftw[i][0] = in2[i];
        x->in_ir_fftw[i][1] = 0.0;
    }
    
    
    fftw_execute(x->plan1);
    fftw_execute(x->plan2);
    
    
    for(i=0; i<n; i++){
        x->out_final[i][0] =  x->out_signal_fftw[i][0] * x->out_ir_fftw[i][0];
        x->out_final[i][1] =  x->out_signal_fftw[i][1] * x->out_ir_fftw[i][1];
    }
    
    
    
    
     for( i = 0; i < n; i++ )
     {
     realD = x->out_signal_fftw[i][0];
     imagD = x->out_signal_fftw[i][1];
     realS = x->out_ir_fftw[i][0];
     imagS = x->out_ir_fftw[i][1];
     x->out_final[i][0] = (realD * realS - imagD * imagS)*mux;
     x->out_final[i][1] = (realD * imagS + imagD * realS)*mux;
     }
    
    
    
    
    fftw_execute(x->plan3);
    
    for(i=0; i<n; i++){
        out[i] = x->in_final[i][0];
    }
    
    
    return (w+6);
}



void conv_tilde_dsp(t_conv_tilde *x, t_signal **sp)
{
    dsp_add(conv_tilde_perform, 5, x,
            sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
    
    
    x->in_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp[0]->s_n);
    x->out_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *sp[0]->s_n);
    x->plan1 = fftw_plan_dft_1d(sp[0]->s_n, x->in_signal_fftw, x->out_signal_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    
    x->in_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp[0]->s_n);
    x->out_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp[0]->s_n);
    x->plan2 = fftw_plan_dft_1d(sp[0]->s_n, x->in_ir_fftw, x->out_ir_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    
    
    x->in_final = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp[0]->s_n);
    x->out_final= (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sp[0]->s_n);
    x->plan3 = fftw_plan_dft_1d(sp[0]->s_n,x->out_final,x->in_final, FFTW_BACKWARD, FFTW_ESTIMATE);
    
}






void conv_tilde_free(t_conv_tilde *x)
{
    inlet_free(x->x_in2);
    outlet_free(x->x_out);
    
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



void *conv_tilde_new(void)
{
    t_conv_tilde *x = (t_conv_tilde *)pd_new(conv_tilde_class);
    // This is called every time an object is created.
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    //Additional signal inlet. The argument is a reference to the lookup tail symbolic selector "signal".
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    //Additional signal outlet. This is made in same way.
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
    CLASS_MAINSIGNALIN(conv_tilde_class,t_conv_tilde,f);
    
}
//Name object "conv~"
//Instructor。Means that it is initialized every time
//Destructor。 we need destructor because we created an additional inlet.
//Allocate static memory.
//Graphic settings for objects. Make it the default.
//Use A_GIMME if you need more than 6 arguments for the symbolic object. 0 means End.
//By attaching a dsp selector, it becomes a signal class, and by having A_CANT, invalid dsp cannot be sent.
//In CLASS_MAINSIGNALIN, the inlet that exists from the beginning becomes the signal inlet.



