#include <math.h>
#include <string.h>
#include "fftw3.h"
#include "m_pd.h"

static t_class *conv_tilde_class;
//インレットが2つ
//1つめが元となる信号
//2つめがIRのデータ


typedef struct _conv_tilde{
    t_object x_obj;  //
    
    t_sample x_conv_before; //fft前の元データ
    t_sample x_ir_before; //fft前のIRデータ
    t_sample x_conv_after;//fft後の元データ
    t_sample x_ir_after; //fft後のIRデータ
    t_float f;//シグナルインレットにfloat型のメッセージが入った時それをシグナルに変えるときに使います。
    
    
    fftw_complex *in_signal_fftw,*out_signal_fftw, *in_ir_fftw, *out_ir_fftw, *in_final, *out_final;
    fftw_plan plan1, plan2, plan3;
    
    
    t_inlet *x_in2;
    t_outlet *x_out;
    
    
}t_conv_tilde;






t_int *conv_tilde_perform(t_int *w)
{
    
    t_conv_tilde *x = (t_conv_tilde *)(w[1]);
    t_sample  *in1 =    (t_sample *)(w[2]);
    t_sample  *in2 =    (t_sample *)(w[3]);
    t_sample  *out =    (t_sample *)(w[4]);
    int          n =           (int)(w[5]);
    
    int i;
    int fftsize = 44100;
    //int l = 44100;
    
    
    
    
    x->in_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    x->out_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *fftsize);
    x->plan1 = fftw_plan_dft_1d(fftsize, x->in_signal_fftw, x->out_signal_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    
    x->in_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    x->out_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *fftsize);
    x->plan2 = fftw_plan_dft_1d(fftsize, x->in_ir_fftw, x->out_ir_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
    
    
    x->in_final = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    x->out_final= (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *fftsize);
    x->plan3 = fftw_plan_dft_1d(fftsize,x->out_final,x->in_final, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    
    
    for(i=0; i<fftsize; i++){
        x->in_signal_fftw[i][0] = in1[i];
        x->in_signal_fftw[i][1] = 0.0;
        
        x->in_ir_fftw[i][0] = in2[i];
        x->in_ir_fftw[i][1] = 0.0;
    }
    
    
    /*
     for(; i<l; i++){
     x->in_signal_fftw[i][0] = 0.0;
     x->in_signal_fftw[i][1] = 0.0;
     x->in_ir_fftw[i][0] = 0.0;
     x->in_ir_fftw[i][1] = 0.0;
     }
     */
    
    
    fftw_execute(x->plan1);
    fftw_execute(x->plan2);
    
    
    
    
    
    for(i=0; i<fftsize; i++){
        x->out_final[i][0] =  x->out_signal_fftw[i][0] * x->out_ir_fftw[i][0];
        x->out_final[i][1] =  x->out_signal_fftw[i][1] * x->out_ir_fftw[i][1];
    }
    
    
    
    
    
    /*
     double mux = 1.0/fftsize;
     
     for( i = 0; i < fftsize; i++ )
     {
     double realD = x->in_signal_fftw[i][0]; //freq_dataがin_dataのoutで、実数部分を入れます。
     double imagD = x->in_signal_fftw[i][1]; //freq_dataがin_dataのoutで、虚数部分を入れます。
     double realS = x->in_ir_fftw[i][0]; //freq_sequenceがin_sequenceのoutで、実数部分を入れます。
     double imagS = x->in_ir_fftw[i][1]; //freq_sequenceがin_sequenceのoutで、実数部分を入れます。
     x->out_final[i][0] = (realD * realS - imagD * imagS)*mux;
     x->out_final[i][1] = (realD * imagS + imagD * realS)*mux;
     
     }
     */
    
    
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
    //オブジェクトが生成される度に呼び出される。
    
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    //追加のシグナルインレット。引数はルックアップテーウルのシンボリックセレクター「シグナル」への参照。
    x->x_out = outlet_new(&x->x_obj, &s_signal);
    //アウトレットも同じように作られる。
    
    
    
    
    
    
    
    
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
//名前をconvにしたよ。
//インストラクタ。毎回初期化されることを意味する
//デスストラクタ。追加のインレットを使用したため必要
//静的メモリの確保を行います。
//オブジェクトのグラフィック設定。デフォにします。
//シンボリックオブジェクトの引数が6つ以上必要ならばA_GIMMEを使用。終わりの０
//dspセレクターをつけておくことでシグナルクラスになりA_CANTを持つことで不正なdspは送れなくなる
//CLASS_MAINSIGNALINにより最初からあるインレットがシグナルインレットとなります。



