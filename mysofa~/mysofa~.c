#include <math.h>
#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "m_pd.h"
#include "mysofa.h"

static t_class *mysofa_tilde_class;

#define MAX_BLOCKSIZE 8192

typedef struct _mysofa_tilde{
    t_object x_obj;  //
    t_float rightIR[MAX_BLOCKSIZE];
    t_float leftIR[MAX_BLOCKSIZE];
    t_int filter_length;
    t_int err;
    t_float f;
    t_float azimuth;
    t_float elevation;
    t_float distance;
    t_float values[3];
    t_float x,y,z;
    t_float leftDelay;
    t_float rightDelay;




    t_inlet *x_in2;
    t_inlet *x_in3;
    t_inlet *x_in4;
    t_outlet *x_out1;
    t_outlet *x_out2;
    struct MYSOFA_EASY *hrtf;
}t_mysofa_tilde;


// MAX_BLOCKSIZE 8192
//blockScale = MAX_BLOCKSIZE / blocksize;
//scaledBlocksize = blocksize * blockScale;
//blocksizeDelta = MAX_BLOCKSIZE -1 - scaledBlocksize


t_int *mysofa_tilde_perform(t_int *w)
{

    t_mysofa_tilde *x = (t_mysofa_tilde *)(w[1]);
    t_sample  *in =    (t_sample *)(w[2]);
    t_sample  *out1 =    (t_sample *)(w[3]);
    t_sample  *out2 =    (t_sample *)(w[4]);
    int       n =           (int)(w[5]);

    
    int i;
    x->leftDelay = 0;
    x->rightDelay = 0;




    x->values[0] = x->azimuth;
    x->values[1] = x->elevation;
    x->values[2] = x->distance;

    mysofa_s2c(x->values);

    x->x = x->values[0];
    x->y = x->values[1];
    x->z = x->values[2];

s

    mysofa_getfilter_float(x->hrtf, x->x, x->y, x->z,x->leftIR,x->rightIR,&x->leftDelay,&x->rightDelay);

   for(i=0;i<n;i++){
        out1[i] = x->leftIR[i];
        out2[i] = x->rightIR[i];
    }


    return (w+6);

}



void mysofa_tilde_dsp(t_mysofa_tilde *x, t_signal **sp)
{
    
    dsp_add(mysofa_tilde_perform,
            5,
            x,
            sp[0]->s_vec, //in_signal
            sp[1]->s_vec, //out_signal_r
            sp[2]->s_vec, //out_signal_l
            sp[0]->s_n);
    
    
    x->hrtf = (struct MYSOFA_EASY*)malloc(sizeof(struct MYSOFA_EASY));
    x->hrtf = mysofa_open("RIEC_hrir_subject_001.sofa", 48000, &x->filter_length, &x->err);
    

}




void mysofa_tilde_free(t_mysofa_tilde *x)
{
    inlet_free(x->x_in2);
    inlet_free(x->x_in3);
    inlet_free(x->x_in4);
    outlet_free(x->x_out1);
    outlet_free(x->x_out2);
    mysofa_close(x->hrtf);

}




void *mysofa_tilde_new(t_floatarg f)
{
    t_mysofa_tilde *x = (t_mysofa_tilde *)pd_new(mysofa_tilde_class);
    // This is called every time an object is created.
    x->x_in2 = floatinlet_new(&x->x_obj, &x->azimuth);
    x->x_in3 = floatinlet_new(&x->x_obj, &x->elevation);
    x->x_in4 = floatinlet_new(&x->x_obj, &x->distance);


    x->x_out1 = outlet_new(&x->x_obj, &s_signal);
    x->x_out2 = outlet_new(&x->x_obj, &s_signal);

    //Additional signal outlet. This is made in same way.
    return (void *)x;
}





void mysofa_tilde_setup(void){

    mysofa_tilde_class=class_new(gensym("mysofa~"),
                                 (t_newmethod)mysofa_tilde_new,
                                 (t_method)mysofa_tilde_free,
                                 sizeof(t_mysofa_tilde),
                                 CLASS_DEFAULT,
                                 A_DEFFLOAT,0);

    class_addmethod(mysofa_tilde_class,
                    (t_method)mysofa_tilde_dsp,
                    gensym("dsp"),A_CANT,0);
    CLASS_MAINSIGNALIN(mysofa_tilde_class, t_mysofa_tilde, f);

    post("singen~ (c) 2020 hoge hoge");

}
//Name object "mysofa~"
//Instructor。Means that it is initialized every time
//Destructor。 we need destructor because we created an additional inlet.
//Allocate static memory.
//Graphic settings for objects. Make it the default.
//Use A_GIMME if you need more than 6 arguments for the symbolic object. 0 means End.
//By attaching a dsp selector, it becomes a signal class, and by having A_CANT, invalid dsp cannot be sent.
//In CLASS_MAINSIGNALIN, the inlet that exists from the beginning becomes the signal inlet.
