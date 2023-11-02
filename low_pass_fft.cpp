#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>

#include <cstdlib>
#include <cmath>
#include "complexo.h"
#include "mat.h"
#include "array.h"


bool FFT(int dir,int m,Array <complexo> &c);

bool P2(int n, int &count,int &aux){

    aux   =  1;
    count =  0;
    while (aux < n){
        aux<<=1;
        count++;
    }
    //aux >>= 1;
    return (bool)(n?!(n&(n - 1)):0);
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
bool FFT2D(matriz <complexo> &c,int nx,int ny,int dir){
   int i,j;
   int m,twopm;
   
   /* Transform the rows */
   Array <complexo> cc(nx);
   
   if (!P2(nx,m,twopm) || twopm != nx)
      return(false);
   for (j=0;j<ny;j++) {    
      for (i=0;i<nx;i++) {
         cc[i] = c[i][j];
      }
      FFT(dir,m,cc);
      for (i=0;i<nx;i++) {
         c[i][j] = cc[i];
      }
   }
   
   cc.setDim(ny);
   /* Transform the columns */
   if (!P2(ny,m,twopm) || twopm != ny)
      return(false);
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         cc[j] = c[i][j];
      }
      FFT(dir,m,cc);
      for (j=0;j<ny;j++) {
         c[i][j] = cc[j];
      }
   }
   
   return(true);
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/

void swap(complexo &a, complexo &b){
	complexo aux = a;
	
	a = b;
	b = aux;
}

bool FFT(int dir,int m,Array <complexo> &c){
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i = 0;i < m; i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn - 1;i++) {
      if (i < j)
         swap(c[i], c[j]);
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j = 0;j < l1; j++) {
         for (i = j;i < nn; i+=l2) {
            i1 = i + l1;
            complexo t = complexo(u1*c[i1].real() - u2*c[i1].imag(),u1*c[i1].imag() + u2*c[i1].real());
            c[i1] = c[i] - t;
            c[i] += t;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++)
         c[i] /= (double)nn;
         
   }

   return(true);
}

void create_highpass_filter(cv::Mat &dft_Filter, int limite){
	cv::Mat tmp = cv::Mat(dft_Filter.rows, dft_Filter.cols, CV_64F);
	cv::Point center = cv::Point(dft_Filter.rows / 2, dft_Filter.cols / 2);
	double radius;


	for(int i = 0; i < dft_Filter.rows; i++){
		for(int j = 0; j < dft_Filter.cols; j++){
			radius = (double) sqrt(pow((i - center.x), 2.0) + pow((double) (j - center.y), 2.0));
			if(radius < limite)			
				tmp.at<double>(i,j) = 0;
			else
				tmp.at<double>(i,j) = 1;
		}
	}
	dft_Filter = tmp;
}

void create_lowpass_filter(cv::Mat &dft_Filter, int limite){
	cv::Mat tmp = cv::Mat(dft_Filter.rows, dft_Filter.cols, CV_64F);
	cv::Point center = cv::Point(dft_Filter.rows / 2, dft_Filter.cols / 2);
	double radius;

	for(int i = 0; i < dft_Filter.rows; i++){
		for(int j = 0; j < dft_Filter.cols; j++){
			radius = (double) sqrt(pow((i - center.x), 2.0) + pow((double) (j - center.y), 2.0));
			if(radius < limite)			
				tmp.at<double>(i,j) = 1;
			else
				tmp.at<double>(i,j) = 0;
		}
	}
	dft_Filter = tmp;
}

using namespace cv;
using namespace std;

int main( int argc, char** argv ){
	char* imageName = argv[1];

 	Mat I;
 	
	
	if(argc==2){
		I = imread(imageName, 1);
		if(!I.data ){
   			cout << "não consegui abrir a imagem" <<endl;
   			return -1;
 		}
	}else return -1;
		
 	namedWindow(imageName, cv::WINDOW_AUTOSIZE);
	imshow(imageName,I);

	Mat BW;
	cvtColor(I,BW, cv::COLOR_BGR2GRAY);
	
	namedWindow("Tons de cinza", cv::WINDOW_AUTOSIZE);
	imshow("Tons de cinza",BW);
 
	int Ny=BW.rows,Nx=BW.cols;
	int NFFTx = (int)pow(2.0, ceil(log((double)Nx)/log(2.0)));  
	int NFFTy = (int)pow(2.0, ceil(log((double)Ny)/log(2.0)));
	cout << Nx <<" "<<Ny<<endl;	
	
	Mat   fourier[] = {Mat::zeros(NFFTy,NFFTx, CV_64F),Mat::zeros(NFFTy,NFFTx, CV_64F)};
	for(int i = 0;i < Ny;i++)
		for(int j = 0;j < Nx;j++){
			fourier[0].at<double>(i,j) = (double)BW.at<unsigned char>(i,j);		
		}
	
       matriz <complexo> A(NFFTy,NFFTx);
       
	for(int i = 0;i < NFFTy;i++)
		for(int j = 0;j < NFFTx;j++){
			A[i][j].real() = 0;
			A[i][j].imag() = 0;		
		}
	
	for(int i = 0;i < NFFTy;i++)
		for(int j = 0;j < NFFTx;j++){
			A[i][j].real() = fourier[0].at<double>(i,j);
			A[i][j].imag() = 0;		
		}
		
	
	FFT2D(A,NFFTy,NFFTx,1);
	
	//-----------------------------------------
	complexo tmp13,tmp24;
	const int m2 = NFFTy/2, n2 = NFFTx/2;  	
	for(int i = 0; i < m2; i++){
		for(int k = 0; k < n2; k++){
          		tmp13         = A[i][k];
          		A[i][k]       = A[i+m2][k+n2];
          		A[i+m2][k+n2] = tmp13;

          		tmp24         = A[i+m2][k];
          		A[i+m2][k]    = A[i][k+n2];
          		A[i][k+n2]    = tmp24;
     		}
	}
	
	for(int i = 0;i < NFFTy;i++)
		for(int j = 0;j < NFFTx;j++){
			fourier[0].at<double>(i,j) = A[i][j].real();
			fourier[1].at<double>(i,j) = A[i][j].imag();		
		}
	
	///-------------------------------------------------------
	Mat filtro = fourier[0].clone();
	create_lowpass_filter(filtro,15);	
	fourier[0] = fourier[0].mul(filtro);
	fourier[1] = fourier[1].mul(filtro);	
	namedWindow( "Espectro", cv::WINDOW_AUTOSIZE );
	Mat espectro;	
	normalize(fourier[0],espectro, 0, 1, NORM_MINMAX);	
        imshow( "Espectro",espectro);
	///-------------------------------------------------------
	for(int i = 0;i < NFFTy;i++)
		for(int j = 0;j < NFFTx;j++){
			A[i][j].real() = fourier[0].at<double>(i,j);
			A[i][j].imag() = fourier[1].at<double>(i,j);		
		}

	for(int i = 0; i < m2; i++){
		for(int k = 0; k < n2; k++){
          		tmp13         = A[i][k];
          		A[i][k]       = A[i+m2][k+n2];
          		A[i+m2][k+n2] = tmp13;

          		tmp24         = A[i+m2][k];
          		A[i+m2][k]    = A[i][k+n2];
          		A[i][k+n2]    = tmp24;
     		}
	}
	
	//-----------------------------------------
	FFT2D(A,NFFTy,NFFTx,-1);
	
	Mat out(Ny,Nx, CV_64F);
	for(int i = 0;i < Ny;i++)
		for(int j = 0;j < Nx;j++){
			out.at<double>(i,j) = A[i][j].real();	
		}
	
	namedWindow( "passa baixa FFT", cv::WINDOW_AUTOSIZE );
	normalize(out,out, 0, 1, NORM_MINMAX);
	
	imshow( "Test image",out);
	waitKey(0);

 	return 0;
}
