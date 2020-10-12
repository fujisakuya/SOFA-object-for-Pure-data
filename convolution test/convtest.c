#include <stdio.h>
#include <fftw3.h>
#include <math.h>// max amplitude for a given sample
FILE  *filename1,*filename2, *filename3, *plotfile1,*plotfile2;

int main(int argc, const char * argv[]) {

    const int fftsize = 1024;
    int l = 44100;
    double re,im;
    float amplitude[fftsize];


    fftw_complex *in_signal_fftw , *out_signal_fftw;
    fftw_complex *in_ir_fftw , *out_ir_fftw;
    fftw_complex *in_final, *out_final;
    fftw_plan plan1, plan2, plan3;
    double realD,imagD,realS,imagS;
    double mux = 1.0/fftsize;


    int i,s1,s2,s3;




    if((filename1 = fopen("conv1.txt","w"))==NULL){
      printf("FILE not open\n");
      return -1;
    }

    if((filename2 = fopen("conv2.txt","w"))==NULL){
      printf("FILE not open\n");
      return -1;
    }

    if((filename3 = fopen("conv3.txt","w"))==NULL){
      printf("FILE not open\n");
      return -1;
    }

    if((plotfile1 = fopen("convplot1.txt","w"))==NULL){
      printf("FILE not open\n");
      return -1;
    }

    if((plotfile2 = fopen("convplot2.txt","w"))==NULL){
      printf("FILE not open\n");
      return -1;
    }









    in_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    out_signal_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *fftsize);
    plan1 = fftw_plan_dft_1d(fftsize, in_signal_fftw, out_signal_fftw, FFTW_FORWARD, FFTW_ESTIMATE);


    in_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    out_ir_fftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *fftsize);
    plan2 = fftw_plan_dft_1d(fftsize, in_ir_fftw, out_ir_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

    in_final = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    out_final= (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize);
    plan3 = fftw_plan_dft_1d(fftsize, out_final, in_final, FFTW_BACKWARD, FFTW_ESTIMATE);



    for(i=0; i<fftsize; i++){
      in_signal_fftw[i][0] = sin((float)400.0 * 2 * M_PI * i / 44100.0 );
      in_signal_fftw[i][1] = 0.0;

      in_ir_fftw[i][0] = sin((float)440.0 * 2 * M_PI * i / 44100.0 );
      in_ir_fftw[i][1] = 0.0;
    }


/*
    for(i=fftsize; i<l; i++){
      in_signal_fftw[i][0] = 0.0;
      in_signal_fftw[i][1] = 0.0;
      in_ir_fftw[i][0] = 0.0;
      in_ir_fftw[i][1] = 0.0;
    }
*/

    fftw_execute(plan1);
    fftw_execute(plan2);








/*
    plan1 = fftw_plan_dft_1d(fftsize, out_signal_fftw, in_signal_fftw, FFTW_BACKWARD, FFTW_ESTIMATE);
     fftw_execute(plan1);
       for(i = 0; i<fftsize; i++){
           fprintf(filename5," bin:%d , re.out:%f,   im.out:%f,  re.in:%f  im.in:%f  \n",i,out_signal_fftw[i][0],out_signal_fftw[i][1],in_signal_fftw[i][0],in_signal_fftw[i][1]);
     }
*/



      for( i = 0; i < fftsize; i++ )  {
        realD = out_signal_fftw[i][0]; //freq_dataがin_dataのoutで、実数部分を入れます。
        imagD = out_signal_fftw[i][1]; //freq_dataがin_dataのoutで、虚数部分を入れます。
        realS = out_ir_fftw[i][0]; //freq_sequenceがin_sequenceのoutで、実数部分を入れます。
        imagS = out_ir_fftw[i][1]; //freq_sequenceがin_sequenceのoutで、実数部分を入れます。
        out_final[i][0] = (realD * realS - imagD * imagS)*mux;
        out_final[i][1] = (realD * imagS + imagD * realS)*mux;
      }



      for(i = 0; i<fftsize; i++){
          fprintf(filename2," bin:%d , re.out:%f,   im.out:%f,  re.in:%f  im.in:%f  \n",i,out_final[i][0],out_final[i][1],in_final[i][0],in_final[i][1]);
    }


/*

      for(i=0;i<fftsize;i++){
        re = in_final[i][0];         // real complex number
        im = in_final[i][1];         // Imaginary complex number
        amplitude[i] = sqrt(re*re + im*im);
      }
  */


      /*writing file to plot */
      for(i = 0; i<fftsize; i++){
        fprintf(plotfile1,"%d, %lf\n",i,amplitude[i]);
      }

      /*Inverse fft*/
    fftw_execute(plan3);



   /*writing file to inverse fft*/
     for(i = 0; i<fftsize; i++){
         fprintf(filename3," bin:%d , re.out:%f,   im.out:%f,  re.in:%f  im.in:%f  \n",i,out_final[i][0],out_final[i][1],in_final[i][0],in_final[i][1]);
   }



   for(i = 0; i<fftsize; i++){
       fprintf(plotfile2,"%d, %f\n",i,in_final[i][0]);
   }




		fftw_destroy_plan(plan1);
    fftw_destroy_plan(plan2);
    fftw_destroy_plan(plan3);
    fftw_free(out_signal_fftw);
    fftw_free(out_ir_fftw);
    fftw_free(out_final);
    fftw_free(in_signal_fftw);
    fftw_free(in_ir_fftw);
    fftw_free(in_final);
    return 0;
}