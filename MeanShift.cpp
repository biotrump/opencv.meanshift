//http://qiankanglai.me/2012/03/19/meanshift/
#include <stdio.h>
#include "MeanShift.h"

//https://vcansimplify.wordpress.com/2014/08/15/604/
//DRAWING REGION ADJACENCY GRAPHS
//RAList : Region Adjacency List
//RAM : Region Adajacency Matrix
//RAG : Region Adjacency Graph
RAList::RAList( void )
{
	label			= -1;
	next			= 0;	//NULL
}

RAList::~RAList( void )
{}

int RAList::Insert(RAList *entry)
{
	if(!next)
	{
		next		= entry;
		entry->next = 0;
		return 0;
	}
	if(next->label > entry->label)
	{
		entry->next	= next;
		next		= entry;
		return 0;
	}
	exists	= 0;
	cur		= next;
	while(cur)
	{
		if(entry->label == cur->label)
		{
			exists = 1;
			break;
		}
		else if((!(cur->next))||(cur->next->label > entry->label))
		{
			entry->next	= cur->next;
			cur->next	= entry;
			break;
		}
		cur = cur->next;
	}
	return (int)(exists);
}

int MeanShiftRGB(const IplImage* img, int **labels)
{
	DECLARE_TIMING(timer);
	START_TIMING(timer);

	int level = 1;
	double COLOR_RADIUS2=color_radius*color_radius;
	int minRegion = 50;

	// use Lab rather than L*u*v!
	// since Luv may produce noise points
	IplImage *result = cvCreateImage(cvGetSize(img),img->depth,img->nChannels);
	cvCvtColor(img, result, CV_RGB2Lab);

	// Step One. Filtering stage of meanshift segmentation
	// http://rsbweb.nih.gov/ij/plugins/download/Mean_Shift.java
	for(int i=0;i<img->height;i++)
		for(int j=0;j<img->width;j++)
		{
			int ic = i;
			int jc = j;
			int icOld, jcOld;
			float LOld, UOld, VOld;
			//Lab->LUV
			float L = (float)((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 0];
			float U = (float)((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 1];
			float V = (float)((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 2];
			// in the case of 8-bit and 16-bit images R, G and B are converted to floating-point format and scaled to fit 0 to 1 range
			// http://opencv.willowgarage.com/documentation/c/miscellaneous_image_transformations.html
			L = L*100/255;
			U = U-128;
			V = V-128;
			double shift = 5;
			for (int iters=0;shift > 3 && iters < 100;iters++)
			{
				icOld = ic;
				jcOld = jc;
				LOld = L;
				UOld = U;
				VOld = V;

				float mi = 0;
				float mj = 0;
				float mL = 0;
				float mU = 0;
				float mV = 0;
				int num=0;

				int i2from = max(0,i-spatial_radius), i2to = min(img->height, i+spatial_radius+1);
				int j2from = max(0,j-spatial_radius), j2to = min(img->width, j+spatial_radius+1);
				for (int i2=i2from; i2 < i2to;i2++) {
					for (int j2=j2from; j2 < j2to; j2++) {
						float L2 = (float)((uchar *)(result->imageData + i2*img->widthStep))[j2*result->nChannels + 0],
							U2 = (float)((uchar *)(result->imageData + i2*img->widthStep))[j2*result->nChannels + 1],
							V2 = (float)((uchar *)(result->imageData + i2*img->widthStep))[j2*result->nChannels + 2];
						L2 = L2*100/255;
						U2 = U2-128;
						V2 = V2-128;

						double dL = L2 - L;
						double dU = U2 - U;
						double dV = V2 - V;
						if (dL*dL+dU*dU+dV*dV <= COLOR_RADIUS2) {
							mi += i2;
							mj += j2;
							mL += L2;
							mU += U2;
							mV += V2;
							num++;
						}
					}
				}
				float num_ = 1.f/num;
				L = mL*num_;
				U = mU*num_;
				V = mV*num_;
				ic = (int) (mi*num_+0.5);
				jc = (int) (mj*num_+0.5);
				int di = ic-icOld;
				int dj = jc-jcOld;
				double dL = L-LOld;
				double dU = U-UOld;
				double dV = V-VOld;

				shift = di*di+dj*dj+dL*dL+dU*dU+dV*dV;
			}

			L = L*255/100;
			U = U+128;
			V = V+128;
			((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 0] = (uchar)L;
			((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 1] = (uchar)U;
			((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 2] = (uchar)V;
		}

		IplImage *tobeshow = cvCreateImage(cvGetSize(img),img->depth,img->nChannels);
		cvCvtColor(result, tobeshow, CV_Lab2RGB);
		cvSaveImage("filtered.png", tobeshow);
		cvReleaseImage(&tobeshow);

		// Step Two. Cluster/segmentation
		// Connect
		int regionCount = 0;
		int *modePointCounts = new int[img->height*img->width];
		memset(modePointCounts, 0, img->width*img->height*sizeof(int));
		float *mode = new float[img->height*img->width*3];
		{
			int label = -1;
			for(int i=0;i<img->height;i++)
				for(int j=0;j<img->width;j++)
					labels[i][j] = -1;
			for(int i=0;i<img->height;i++)
				for(int j=0;j<img->width;j++)
					if(labels[i][j]<0)
					{
						int cnts=0;
						labels[i][j] = ++label;
						modePointCounts[label]=1;
						float L = (float)((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 0],
							U = (float)((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 1],
							V = (float)((uchar *)(result->imageData + i*img->widthStep))[j*result->nChannels + 2];
						mode[label*3+0] = L*100/255;
						mode[label*3+1] = 354*U/255-134;
						mode[label*3+2] = 256*V/255-140;
						// Fill
						std::stack<CvPoint> neighStack;
						neighStack.push(cvPoint(i,j));
						const int dxdy[][2] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};
						while(!neighStack.empty())
						{
							CvPoint p = neighStack.top();
							neighStack.pop();
							for(int k=0;k<8;k++)
							{
								int i2 = p.x+dxdy[k][0], j2 = p.y+dxdy[k][1];
								if(i2>=0 && j2>=0 && i2<img->height && j2<img->width && labels[i2][j2]<0
									&& color_distance(result, i,j,i2,j2)<COLOR_RADIUS2)
								{
									labels[i2][j2] = label;
									neighStack.push(cvPoint(i2,j2));
									modePointCounts[label]++;
									L = (float)((uchar *)(result->imageData + i2*img->widthStep))[j2*result->nChannels + 0];
									U = (float)((uchar *)(result->imageData + i2*img->widthStep))[j2*result->nChannels + 1];
									V = (float)((uchar *)(result->imageData + i2*img->widthStep))[j2*result->nChannels + 2];
									mode[label*3+0] += L*100/255;
									mode[label*3+1] += 354*U/255-134;
									mode[label*3+2] += 256*V/255-140;
									cnts++;
								}
							}
						}
						if(cnts){
							mode[label*3+0] /= modePointCounts[label];
							mode[label*3+1] /= modePointCounts[label];
							mode[label*3+2] /= modePointCounts[label];
						}else{
							std::cout << "orphan mode:" << label << ",@(" << i <<"," << j <<")" <<std::endl;
							std::cout << "modePointCounts:" << modePointCounts[label] <<std::endl;
							
						}
					}
					//current Region count
					regionCount = label+1;
		}
		std::cout<<"Mean Shift(Connect):"<<regionCount<<std::endl;
		int oldRegionCount = regionCount;

		// TransitiveClosure : merge nearby cluster
		for(int counter = 0, deltaRegionCount = 1; counter<5 && deltaRegionCount>0; counter++)
		{
			// 1.Build RAM using classifiction structure
			RAList *raList = new RAList [regionCount], *raPool = new RAList [10*regionCount];	//10 is hard coded!
			for(int i = 0; i < regionCount; i++)
			{
				raList[i].label = i;
				raList[i].next = NULL;
			}
			for(int i = 0; i < regionCount*10-1; i++)
			{
				raPool[i].next = &raPool[i+1];
			}
			raPool[10*regionCount-1].next = NULL;
			RAList	*raNode1, *raNode2, *oldRAFreeList, *freeRAList = raPool;
			for(int i=0;i<img->height;i++)
				for(int j=0;j<img->width;j++)
				{
					if(i>0 && labels[i][j]!=labels[i-1][j])
					{
						// Get 2 free node
						raNode1			= freeRAList;
						raNode2			= freeRAList->next;
						oldRAFreeList	= freeRAList;
						freeRAList		= freeRAList->next->next;
						// connect the two region
						raNode1->label	= labels[i][j];
						raNode2->label	= labels[i-1][j];
						if(raList[labels[i][j]].Insert(raNode2))	//already exists!
							freeRAList = oldRAFreeList;
						else
							raList[labels[i-1][j]].Insert(raNode1);
					}
					if(j>0 && labels[i][j]!=labels[i][j-1])
					{
						// Get 2 free node
						raNode1			= freeRAList;
						raNode2			= freeRAList->next;
						oldRAFreeList	= freeRAList;
						freeRAList		= freeRAList->next->next;
						// connect the two region
						raNode1->label	= labels[i][j];
						raNode2->label	= labels[i][j-1];
						if(raList[labels[i][j]].Insert(raNode2))
							freeRAList = oldRAFreeList;
						else
							raList[labels[i][j-1]].Insert(raNode1);
					}
				}

				// 2.Treat each region Ri as a disjoint set
				for(int i = 0; i < regionCount; i++)
				{
					RAList	*neighbor = raList[i].next;
					while(neighbor)
					{
						if(color_distance(&mode[3*i], &mode[3*neighbor->label])<COLOR_RADIUS2)
						{
							int iCanEl = i, neighCanEl	= neighbor->label;
							while(raList[iCanEl].label != iCanEl) iCanEl = raList[iCanEl].label;
							while(raList[neighCanEl].label != neighCanEl) neighCanEl = raList[neighCanEl].label;
							if(iCanEl<neighCanEl)
								raList[neighCanEl].label = iCanEl;
							else
							{
								//raList[raList[iCanEl].label].label = iCanEl;
								raList[iCanEl].label = neighCanEl;
							}
						}
						neighbor = neighbor->next;
					}
				}
				// 3. Union Find
				for(int i = 0; i < regionCount; i++)
				{
					int iCanEl	= i;
					while(raList[iCanEl].label != iCanEl) iCanEl	= raList[iCanEl].label;
					raList[i].label	= iCanEl;
				}
				// 4. Traverse joint sets, relabeling image.
				int *modePointCounts_buffer = new int[regionCount];
				memset(modePointCounts_buffer, 0, regionCount*sizeof(int));
				float *mode_buffer = new float[regionCount*3];
				int	*label_buffer = new int[regionCount];

				for(int i=0;i<regionCount; i++)
				{
					label_buffer[i]	= -1;
					mode_buffer[i*3+0] = 0;
					mode_buffer[i*3+1] = 0;
					mode_buffer[i*3+2] = 0;
				}
				for(int i=0;i<regionCount; i++)
				{
					int iCanEl	= raList[i].label;
					modePointCounts_buffer[iCanEl] += modePointCounts[i];
					for(int k=0;k<3;k++)
						mode_buffer[iCanEl*3+k] += mode[i*3+k]*modePointCounts[i];
				}
				int	label = -1;
				for(int i = 0; i < regionCount; i++)
				{
					int iCanEl	= raList[i].label;
					if(label_buffer[iCanEl] < 0)
					{
						label_buffer[iCanEl]	= ++label;

						for(int k = 0; k < 3; k++)
							mode[label*3+k]	= (mode_buffer[iCanEl*3+k])/(modePointCounts_buffer[iCanEl]);

						modePointCounts[label]	= modePointCounts_buffer[iCanEl];
					}
				}
				regionCount = label+1;
				for(int i = 0; i < img->height; i++)
					for(int j = 0; j < img->width; j++)
						labels[i][j]	= label_buffer[raList[labels[i][j]].label];

				delete [] mode_buffer;
				delete [] modePointCounts_buffer;
				delete [] label_buffer;

				//Destroy RAM
				delete[] raList;
				delete[] raPool;

				deltaRegionCount = oldRegionCount - regionCount;
				oldRegionCount = regionCount;
				std::cout<<"Mean Shift(TransitiveClosure):"<<regionCount<<std::endl;
		}

		// Prune : small size clusters are removed
		{
			int *modePointCounts_buffer = new int[regionCount];
			float *mode_buffer = new float[regionCount*3];
			int	*label_buffer = new int [regionCount];
			int minRegionCount;

			do{
				minRegionCount = 0;
				// Build RAM again
				RAList *raList = new RAList [regionCount], *raPool = new RAList [10*regionCount];	//10 is hard coded!
				for(int i = 0; i < regionCount; i++)
				{
					raList[i].label = i;
					raList[i].next = NULL;
				}
				for(int i = 0; i < regionCount*10-1; i++)
				{
					raPool[i].next = &raPool[i+1];
				}
				raPool[10*regionCount-1].next = NULL;
				RAList	*raNode1, *raNode2, *oldRAFreeList, *freeRAList = raPool;
				for(int i=0;i<img->height;i++)
					for(int j=0;j<img->width;j++)
					{
						if(i>0 && labels[i][j]!=labels[i-1][j])
						{
							// Get 2 free node
							raNode1			= freeRAList;
							raNode2			= freeRAList->next;
							oldRAFreeList	= freeRAList;
							freeRAList		= freeRAList->next->next;
							// connect the two region
							raNode1->label	= labels[i][j];
							raNode2->label	= labels[i-1][j];
							if(raList[labels[i][j]].Insert(raNode2))	//already exists!
								freeRAList = oldRAFreeList;
							else
								raList[labels[i-1][j]].Insert(raNode1);
						}
						if(j>0 && labels[i][j]!=labels[i][j-1])
						{
							// Get 2 free node
							raNode1			= freeRAList;
							raNode2			= freeRAList->next;
							oldRAFreeList	= freeRAList;
							freeRAList		= freeRAList->next->next;
							// connect the two region
							raNode1->label	= labels[i][j];
							raNode2->label	= labels[i][j-1];
							if(raList[labels[i][j]].Insert(raNode2))
								freeRAList = oldRAFreeList;
							else
								raList[labels[i][j-1]].Insert(raNode1);
						}
					}
					// Find small regions
					for(int i = 0; i < regionCount; i++)
						if(modePointCounts[i] < minRegion)
						{
							minRegionCount++;
							RAList *neighbor = raList[i].next;
							int candidate = neighbor->label;
							float minDistance = color_distance(&mode[3*i], &mode[3*candidate]);
							neighbor = neighbor->next;
							while(neighbor)
							{
								float minDistance2 = color_distance(&mode[3*i], &mode[3*neighbor->label]);
								if(minDistance2<minDistance)
								{
									minDistance = minDistance2;
									candidate = neighbor->label;
								}
								neighbor = neighbor->next;
							}
							int iCanEl = i, neighCanEl	= candidate;
							while(raList[iCanEl].label != iCanEl) iCanEl = raList[iCanEl].label;
							while(raList[neighCanEl].label != neighCanEl) neighCanEl = raList[neighCanEl].label;
							if(iCanEl < neighCanEl)
								raList[neighCanEl].label	= iCanEl;
							else
							{
								//raList[raList[iCanEl].label].label	= neighCanEl;
								raList[iCanEl].label = neighCanEl;
							}
						}
						for(int i = 0; i < regionCount; i++)
						{
							int iCanEl	= i;
							while(raList[iCanEl].label != iCanEl)
								iCanEl	= raList[iCanEl].label;
							raList[i].label	= iCanEl;
						}
						memset(modePointCounts_buffer, 0, regionCount*sizeof(int));
						for(int i = 0; i < regionCount; i++)
						{
							label_buffer[i]	= -1;
							mode_buffer[3*i+0]	= 0;
							mode_buffer[3*i+1]	= 0;
							mode_buffer[3*i+2]	= 0;
						}
						for(int i=0;i<regionCount; i++)
						{
							int iCanEl	= raList[i].label;
							modePointCounts_buffer[iCanEl] += modePointCounts[i];
							for(int k=0;k<3;k++)
								mode_buffer[iCanEl*3+k] += mode[i*3+k]*modePointCounts[i];
						}
						int	label = -1;
						for(int i = 0; i < regionCount; i++)
						{
							int iCanEl	= raList[i].label;
							if(label_buffer[iCanEl] < 0)
							{
								label_buffer[iCanEl]	= ++label;

								for(int k = 0; k < 3; k++)
									mode[label*3+k]	= (mode_buffer[iCanEl*3+k])/(modePointCounts_buffer[iCanEl]);

								modePointCounts[label]	= modePointCounts_buffer[iCanEl];
							}
						}
						regionCount = label+1;
						for(int i = 0; i < img->height; i++)
							for(int j = 0; j < img->width; j++)
								labels[i][j]	= label_buffer[raList[labels[i][j]].label];

						//Destroy RAM
						delete[] raList;
						delete[] raPool;
						std::cout<<"Mean Shift(Prune):"<<regionCount<<std::endl;
			}while(minRegionCount > 0);

			delete [] mode_buffer;
			delete [] modePointCounts_buffer;
			delete [] label_buffer;
		}

		// Output
		STOP_TIMING(timer);
		std::cout<<"Mean Shift(ms):"<<GET_TIMING(timer)<<std::endl;

		cvReleaseImage(&result);
		delete []mode;
		delete []modePointCounts;
		return regionCount;
}

static IplImage *tobeshow16SC3=NULL;
/** @brief A spatial-range joint domain feature space. The range domain is a generic feature space
 * with 3 channels. The range domain is NDI normalized between -32767 to 32767.
 * spatial (x,y) is the image dimension.
 * so there are still 5d , joint spatial-range domain feature space.
 * @param[IN,OUT] labels[][] : preallocated table to store the label of every pixel
 * after meanshift
 * @param[IN] img : cvCreateImageHeader(IPL_DEPTH_16S, 3), 3 channel feature space with
 * 16bit per channel. Each channel coresponds to its NDI.
 * @return : regionCount, cluster numbers
 */
int MeanShiftGeneric3DS16(int16_t *ndi_buf, const int width, const int height,
						  const int nChannels, const uint32_t shift_thrshold,
						  const int16_t max_iters, int **labels)
{
	DECLARE_TIMING(timer);
	START_TIMING(timer);

	int level = 1;
	float RANGE_RADIUS2=range_radius*range_radius;
	int minRegion = 50;

	// Step One. Filtering stage of meanshift segmentation
	// http://rsbweb.nih.gov/ij/plugins/download/Mean_Shift.java
	for(int i=0;i<height;i++)
		for(int j=0;j<width;j++)
		{
			//start pixel as a center
			int ic = i;	int jc = j;

			//NDI feature space of pixel(ic,jc)
			int16_t ND0 = ndi_buf [ (i*width+j)*nChannels + 0];
			int16_t ND1 = ndi_buf [ (i*width+j)*nChannels + 1];
			int16_t ND2 = ndi_buf [ (i*width+j)*nChannels + 2];

			//mean shift smoothing
			int icOld, jcOld;
			int16_t ND0_Old, ND1_Old, ND2_Old;

			uint32_t shift = SHIFT_THRESHOLD_MAX;//mean shift distance
			for (int iters=0;shift > shift_thrshold && iters < max_iters;iters++)
			{
				icOld = ic;
				jcOld = jc;
				ND0_Old = ND0;
				ND1_Old = ND1;
				ND2_Old = ND2;

				float mi = 0;
				float mj = 0;
				float mND0 = 0;
				float mND1 = 0;
				float mND2 = 0;
				int num=0;

				//(hyper) cube as a window instead of a (hyper) ball to
				//save the distance square boundary checks
				//here is just a 2D rectangle bounding box around the pixel (ic,jc)
				int i2from = max(0,i-spatial_radius), i2to = min(height, i+spatial_radius+1);
				int j2from = max(0,j-spatial_radius), j2to = min(width, j+spatial_radius+1);
				//finding the pixels inside the box and get the weight center of this bounding box
				for (int i2=i2from; i2 < i2to;i2++) {
					for (int j2=j2from; j2 < j2to; j2++) {
						//the NDI value of 3 channels @(i2,j2)
						int16_t ND0_2 =	ndi_buf [ (i2*width+j2)*nChannels + 0];
						int16_t ND1_2 =	ndi_buf [ (i2*width+j2)*nChannels + 1];
						int16_t ND2_2 =	ndi_buf [ (i2*width+j2)*nChannels + 2];
						//nearby pixel distance to center pixel in NDI feature space
						float d0_2 = (ND0_2 - ND0) * (ND0_2 - ND0);
						float d1_2 = (ND1_2 - ND1) * (ND1_2 - ND1);
						float d2_2 = (ND2_2 - ND2) * (ND2_2 - ND2);
						float r_d=d0_2 + d1_2 + d2_2;
						//check if the nearby pixel in NDI space is within the bandwidth/radius?
						if (r_d <= RANGE_RADIUS2) {
							//yes, the pixel is inside the radius in featur space,
							//so sum up spatial and range to calculate the weight center
							//of the bounding box later.
							//spatial-range joint domain feature space
							mi += i2;
							mj += j2;
							mND0 += ND0_2;
							mND1 += ND1_2;
							mND2 += ND2_2;
							num++;
						}
					}
				}
				//get the weight center of the bounding box by dividing the pixels
				//inside the boudning box
				float num_ = 1.f/num;
				ND0 = mND0*num_;
				ND1 = mND1*num_;
				ND2 = mND2*num_;
				//move to next center by the mean shift
				ic = (int) (mi*num_+0.5f);
				jc = (int) (mj*num_+0.5f);
				//mean shift vector
				int di = ic-icOld;
				int dj = jc-jcOld;
				float dnd0 = ND0-ND0_Old;
				float dnd1 = ND1-ND1_Old;
				float dnd2 = ND2-ND2_Old;
				//mean shift vector length
				shift = di*di+dj*dj+dnd0*dnd0+dnd1*dnd1+dnd2*dnd2;
			}
			//smoothing current pixel(i,j) by replace the range value
			//with its end point of the mode

			ndi_buf[ (i*width+j)*nChannels + 0] = ND0;
			ndi_buf[ (i*width+j)*nChannels + 1] = ND1;
			ndi_buf[ (i*width+j)*nChannels + 2] = ND2;
		}

		if(!tobeshow16SC3){
			tobeshow16SC3 = cvCreateImageHeader(cvSize(width,height), IPL_DEPTH_16S, nChannels);
		}
		cvSetData(tobeshow16SC3, ndi_buf, width*nChannels*sizeof(int16_t));//width*nChannels*sizeof(IPL_DEPTH_16S) bytes per row
		cvShowImage("filter16SC3", tobeshow16SC3);

		// Step Two. Cluster
		// Connect
		int regionCount = 0;
		int *modePointCounts = new int[height*width];
		memset(modePointCounts, 0, width*height*sizeof(int));
		int16_t *mode = new int16_t[height*width*nChannels];
		{
			int label = -1;
			for(int i=0;i<height;i++)
				for(int j=0;j<width;j++)
					labels[i][j] = -1;
			for(int i=0;i<height;i++)
				for(int j=0;j<width;j++)
					if(labels[i][j]<0)
					{
						labels[i][j] = ++label;
						int16_t ND0 = ndi_buf[ (i*width+j)*nChannels + 0];
						int16_t	ND1 = ndi_buf[ (i*width+j)*nChannels + 1];
						int16_t	ND2 = ndi_buf[ (i*width+j)*nChannels + 2];
						mode[label*nChannels+0] = ND0;
						mode[label*nChannels+1] = ND1;
						mode[label*nChannels+2] = ND2;
						// Fill
						std::stack<CvPoint> neighStack;
						neighStack.push(cvPoint(i,j));
						//8-connected neighborhood
						const int dxdy[][2] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,1},{1,-1},{1,0},{1,1}};
						while(!neighStack.empty())
						{
							CvPoint p = neighStack.top();
							neighStack.pop();
							for(int k=0;k<8;k++)
							{
								int i2 = p.x+dxdy[k][0], j2 = p.y+dxdy[k][1];
								if(i2>=0 && j2>=0 && i2<height && j2<width && labels[i2][j2]<0 &&
									range_distance(ndi_buf, width, nChannels, i,j,i2,j2)<RANGE_RADIUS2)
								{
									labels[i2][j2] = label;
									neighStack.push(cvPoint(i2,j2));
									modePointCounts[label]++;
									ND0 = ndi_buf[ (i2*width+j2)*nChannels + 0];
									ND1 = ndi_buf[ (i2*width+j2)*nChannels + 1];
									ND2 = ndi_buf[ (i2*width+j2)*nChannels + 2];
									mode[label*nChannels+0] += ND0;
									mode[label*nChannels+1] += ND1;
									mode[label*nChannels+2] += ND2;
								}
							}
						}
						mode[label*3+0] /= modePointCounts[label];
						mode[label*3+1] /= modePointCounts[label];
						mode[label*3+2] /= modePointCounts[label];
					}
					//current Region count
					regionCount = label+1;
		}
		std::cout<<"Mean Shift(Connect):"<<regionCount<<std::endl;
		int oldRegionCount = regionCount;

		// TransitiveClosure By Region Adjacency Graph/Matrix
		for(int counter = 0, deltaRegionCount = 1; counter<5 && deltaRegionCount>0; counter++)
		{
			// 1.Build RAM using classifiction structure
			RAList *raList = new RAList [regionCount], *raPool = new RAList [10*regionCount];	//10 is hard coded!
			for(int i = 0; i < regionCount; i++)
			{
				raList[i].label = i;
				raList[i].next = NULL;
			}
			for(int i = 0; i < regionCount*10-1; i++)
			{
				raPool[i].next = &raPool[i+1];
			}
			raPool[10*regionCount-1].next = NULL;
			RAList	*raNode1, *raNode2, *oldRAFreeList, *freeRAList = raPool;
			for(int i=0;i<height;i++)
				for(int j=0;j<width;j++)
				{
					if(i>0 && labels[i][j]!=labels[i-1][j])
					{
						// Get 2 free node
						raNode1			= freeRAList;
						raNode2			= freeRAList->next;
						oldRAFreeList	= freeRAList;
						freeRAList		= freeRAList->next->next;
						// connect the two region
						raNode1->label	= labels[i][j];
						raNode2->label	= labels[i-1][j];
						if(raList[labels[i][j]].Insert(raNode2))	//already exists!
							freeRAList = oldRAFreeList;
						else
							raList[labels[i-1][j]].Insert(raNode1);
					}
					if(j>0 && labels[i][j]!=labels[i][j-1])
					{
						// Get 2 free node
						raNode1			= freeRAList;
						raNode2			= freeRAList->next;
						oldRAFreeList	= freeRAList;
						freeRAList		= freeRAList->next->next;
						// connect the two region
						raNode1->label	= labels[i][j];
						raNode2->label	= labels[i][j-1];
						if(raList[labels[i][j]].Insert(raNode2))
							freeRAList = oldRAFreeList;
						else
							raList[labels[i][j-1]].Insert(raNode1);
					}
				}

				// 2.Treat each region Ri as a disjoint set
				for(int i = 0; i < regionCount; i++)
				{
					RAList	*neighbor = raList[i].next;
					while(neighbor)
					{
						if(range_distance(&mode[3*i], &mode[3*neighbor->label])<RANGE_RADIUS2)
						{
							int iCanEl = i, neighCanEl	= neighbor->label;
							while(raList[iCanEl].label != iCanEl) iCanEl = raList[iCanEl].label;
							while(raList[neighCanEl].label != neighCanEl) neighCanEl = raList[neighCanEl].label;
							if(iCanEl<neighCanEl)
								raList[neighCanEl].label = iCanEl;
							else
							{
								//raList[raList[iCanEl].label].label = iCanEl;
								raList[iCanEl].label = neighCanEl;
							}
						}
						neighbor = neighbor->next;
					}
				}
				// 3. Union Find
				for(int i = 0; i < regionCount; i++)
				{
					int iCanEl	= i;
					while(raList[iCanEl].label != iCanEl) iCanEl	= raList[iCanEl].label;
					raList[i].label	= iCanEl;
				}
				// 4. Traverse joint sets, relabeling image.
				int *modePointCounts_buffer = new int[regionCount];
				memset(modePointCounts_buffer, 0, regionCount*sizeof(int));
				float *mode_buffer = new float[regionCount*3];
				int	*label_buffer = new int[regionCount];

				for(int i=0;i<regionCount; i++)
				{
					label_buffer[i]	= -1;
					mode_buffer[i*nChannels+0] = 0;
					mode_buffer[i*nChannels+1] = 0;
					mode_buffer[i*nChannels+2] = 0;
				}
				for(int i=0;i<regionCount; i++)
				{
					int iCanEl	= raList[i].label;
					modePointCounts_buffer[iCanEl] += modePointCounts[i];
					for(int k=0;k<nChannels;k++)
						mode_buffer[iCanEl*nChannels+k] += mode[i*nChannels+k]*modePointCounts[i];
				}
				int	label = -1;
				for(int i = 0; i < regionCount; i++)
				{
					int iCanEl	= raList[i].label;
					if(label_buffer[iCanEl] < 0)
					{
						label_buffer[iCanEl]	= ++label;

						for(int k = 0; k < nChannels; k++)
							mode[label*nChannels+k]	= (mode_buffer[iCanEl*nChannels+k])/(modePointCounts_buffer[iCanEl]);

						modePointCounts[label]	= modePointCounts_buffer[iCanEl];
					}
				}
				regionCount = label+1;
				for(int i = 0; i < height; i++)
					for(int j = 0; j < width; j++)
						labels[i][j]	= label_buffer[raList[labels[i][j]].label];

				delete [] mode_buffer;
				delete [] modePointCounts_buffer;
				delete [] label_buffer;

				//Destroy RAM
				delete[] raList;
				delete[] raPool;

				deltaRegionCount = oldRegionCount - regionCount;
				oldRegionCount = regionCount;
				std::cout<<"Mean Shift(TransitiveClosure):"<<regionCount<<std::endl;
		}

		// Prune small region
		{
			int *modePointCounts_buffer = new int[regionCount];
			int *mode_buffer = new int[regionCount*3];
			int	*label_buffer = new int [regionCount];
			int minRegionCount;

			do{
				minRegionCount = 0;
				// Build RAM again
				RAList *raList = new RAList [regionCount], *raPool = new RAList [10*regionCount];	//10 is hard coded!
				for(int i = 0; i < regionCount; i++)
				{
					raList[i].label = i;
					raList[i].next = NULL;
				}
				for(int i = 0; i < regionCount*10-1; i++)
				{
					raPool[i].next = &raPool[i+1];
				}
				raPool[10*regionCount-1].next = NULL;
				RAList	*raNode1, *raNode2, *oldRAFreeList, *freeRAList = raPool;
				for(int i=0;i<height;i++)
					for(int j=0;j<width;j++)
					{
						if(i>0 && labels[i][j]!=labels[i-1][j])
						{
							// Get 2 free node
							raNode1			= freeRAList;
							raNode2			= freeRAList->next;
							oldRAFreeList	= freeRAList;
							freeRAList		= freeRAList->next->next;
							// connect the two region
							raNode1->label	= labels[i][j];
							raNode2->label	= labels[i-1][j];
							if(raList[labels[i][j]].Insert(raNode2))	//already exists!
								freeRAList = oldRAFreeList;
							else
								raList[labels[i-1][j]].Insert(raNode1);
						}
						if(j>0 && labels[i][j]!=labels[i][j-1])
						{
							// Get 2 free node
							raNode1			= freeRAList;
							raNode2			= freeRAList->next;
							oldRAFreeList	= freeRAList;
							freeRAList		= freeRAList->next->next;
							// connect the two region
							raNode1->label	= labels[i][j];
							raNode2->label	= labels[i][j-1];
							if(raList[labels[i][j]].Insert(raNode2))
								freeRAList = oldRAFreeList;
							else
								raList[labels[i][j-1]].Insert(raNode1);
						}
					}
					// Find small regions
					for(int i = 0; i < regionCount; i++)
						if(modePointCounts[i] < minRegion)
						{
							minRegionCount++;
							RAList *neighbor = raList[i].next;
							int candidate = neighbor->label;
							float minDistance = range_distance(&mode[3*i], &mode[3*candidate]);
							neighbor = neighbor->next;
							while(neighbor)
							{
								float minDistance2 = range_distance(&mode[3*i], &mode[3*neighbor->label]);
								if(minDistance2<minDistance)
								{
									minDistance = minDistance2;
									candidate = neighbor->label;
								}
								neighbor = neighbor->next;
							}
							int iCanEl = i, neighCanEl	= candidate;
							while(raList[iCanEl].label != iCanEl) iCanEl = raList[iCanEl].label;
							while(raList[neighCanEl].label != neighCanEl) neighCanEl = raList[neighCanEl].label;
							if(iCanEl < neighCanEl)
								raList[neighCanEl].label	= iCanEl;
							else
							{
								//raList[raList[iCanEl].label].label	= neighCanEl;
								raList[iCanEl].label = neighCanEl;
							}
						}
						for(int i = 0; i < regionCount; i++)
						{
							int iCanEl	= i;
							while(raList[iCanEl].label != iCanEl)
								iCanEl	= raList[iCanEl].label;
							raList[i].label	= iCanEl;
						}
						memset(modePointCounts_buffer, 0, regionCount*sizeof(int));
						for(int i = 0; i < regionCount; i++)
						{
							label_buffer[i]	= -1;
							mode_buffer[nChannels*i+0]	= 0;
							mode_buffer[nChannels*i+1]	= 0;
							mode_buffer[nChannels*i+2]	= 0;
						}
						for(int i=0;i<regionCount; i++)
						{
							int iCanEl	= raList[i].label;
							modePointCounts_buffer[iCanEl] += modePointCounts[i];
							for(int k=0;k<nChannels;k++)
								mode_buffer[iCanEl*nChannels+k] += mode[i*nChannels+k]*modePointCounts[i];
						}
						int	label = -1;
						for(int i = 0; i < regionCount; i++)
						{
							int iCanEl	= raList[i].label;
							if(label_buffer[iCanEl] < 0)
							{
								label_buffer[iCanEl]	= ++label;

								for(int k = 0; k < nChannels; k++)
									mode[label*nChannels+k]	= (mode_buffer[iCanEl*nChannels+k])/(modePointCounts_buffer[iCanEl]);

								modePointCounts[label]	= modePointCounts_buffer[iCanEl];
							}
						}
						regionCount = label+1;
						for(int i = 0; i < height; i++)
							for(int j = 0; j < width; j++)
								labels[i][j]	= label_buffer[raList[labels[i][j]].label];

						//Destroy RAM
						delete[] raList;
						delete[] raPool;
						std::cout<<"Mean Shift(Prune):"<<regionCount<<std::endl;
			}while(minRegionCount > 0);

			delete [] mode_buffer;
			delete [] modePointCounts_buffer;
			delete [] label_buffer;
		}

		// Output
		STOP_TIMING(timer);
		std::cout<<"Mean Shift(ms):"<<GET_TIMING(timer)<<std::endl;

		for(int l=0; l<regionCount;l++){
			printf("mode %d:(%d,%d,%d)\n",l, mode[l*nChannels+0],mode[l*nChannels+1],mode[l*nChannels+2]);
			printf("\tcounts=%d\n",modePointCounts[l]);
		}
		delete []mode;
		delete []modePointCounts;

		cvReleaseImage(&tobeshow16SC3);
		tobeshow16SC3=NULL;
		return regionCount;
}
