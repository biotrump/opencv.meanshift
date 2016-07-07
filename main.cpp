//http://qiankanglai.me/2012/03/19/meanshift/
#include <iostream>
#include <stdio.h>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgproc/imgproc_c.h"

#include "MeanShift.h"
using namespace std;
using namespace cv;

int gui_spatialRad=SPATIAL_RADIUS;
int gui_color_radius=COLOR_RADIUS*10;
int gui_range_radius=RANGE_RADIUS/10;
int gui_shift_threshold=30;

int spatial_radius;
float color_radius;
int range_radius;
float shift = 5;

IplImage *img=NULL;

void on_Meanshift(int arg )  //the callback function
{
	//spatialRad=15,color_radius=20,range_radius=20;
	spatial_radius=gui_spatialRad;
	color_radius=gui_color_radius/10.0f;
	range_radius=gui_range_radius*10;
	shift =gui_shift_threshold/10.0f;
	printf("spatial_radius=%d, color_radius=%f,range_radius=%d,shift=%f\n",
		   spatial_radius, color_radius,range_radius, shift);

	//ilabels : label to every pixel after meanshift
	int **ilabels = new int *[img->height];
	for(int i=0;i<img->height;i++) ilabels[i] = new int [img->width];

	int regionCount = MeanShiftRGB(img, ilabels);
	vector<int> color(regionCount);
	CvRNG rng= cvRNG(cvGetTickCount());
	for(int i=0;i<regionCount;i++)
		color[i] = cvRandInt(&rng);

	// Draw random color
	for(int i=0 ; i<img->height ; i++)
		for(int j=0 ; j<img->width ; j++)
		{
			int cl = ilabels[i][j];
			((uchar *)(img->imageData + i*img->widthStep))[j*img->nChannels + 0] = (color[cl])&255;
			((uchar *)(img->imageData + i*img->widthStep))[j*img->nChannels + 1] = (color[cl]>>8)&255;
			((uchar *)(img->imageData + i*img->widthStep))[j*img->nChannels + 2] = (color[cl]>>16)&255;
		}
	cvShowImage("MeanShift",img);
}

int main(int argc, char* argv[])
{
//	IplImage *img = cvLoadImage("input.png");

	if(argc > 1)
		img = cvLoadImage(argv[1]);
	else
		img = cvLoadImage("input.png");
	if(img==NULL){
		cout << "fail to load image file" << endl;
		exit(1);
	}

	// Mean shift
	cvNamedWindow("MeanShift",CV_WINDOW_AUTOSIZE);
	cvShowImage("MeanShift",img);

	//cvNamedWindow("ms_dst",CV_WINDOW_AUTOSIZE);
	cvCreateTrackbar("spatialRad","MeanShift",&gui_spatialRad, 100, on_Meanshift);
    cvCreateTrackbar("colorRad","MeanShift",&gui_color_radius, 100, on_Meanshift);
	cvCreateTrackbar("rangeRad","MeanShift",&gui_range_radius, 100, on_Meanshift);
	cvCreateTrackbar("shiftThre","MeanShift",&gui_shift_threshold, 100, on_Meanshift);

	cvWaitKey();

	cvDestroyWindow("MeanShift");
	cvReleaseImage(&img);

	return 0;
}
