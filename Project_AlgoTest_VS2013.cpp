// Project_AlgoTest_VS2013.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "time.h"

//darklib
//#include"darknet_invoke.h"



//
#include "commondata.h"
#include "commonfile.h"
#include "CommonFuncs.h"

#include "panorama.hpp"
//#include "geotiff.h"
//#include "relativepose.hpp"
#include "register.hpp"
#include "cvInvoke.hpp"
#include "ba.hpp"
#include "ceresba.hpp"
#include "panorama.hpp"
#include "bundlerio.hpp"
#include "panoba.hpp"
#include "sift.hpp"
#include "triangulate.hpp"
#include "pos.hpp"
#include "funcs.hpp"
#include "pano2plane.h"
#include "pos.hpp"
#include "detect.hpp"
#include "absOri.hpp"
#include "test.h"
#include "LatLong-UTMconversion.h"
#include "readJpegHeader.hpp"
#include "sfm.hpp"
#include "lane.hpp"






//rslib
#include "Aerosol.h"
#include "scatter.hpp"
#include "Polder.h"
#include "bpdf.h"


#include <opencv2/stitching/stitcher.hpp>

//#include"darknet.h"


int ReadPosData(char* posfile, char* imagepath,
	int startindex, int nImageNumber,
	CameraType camType, vector<CameraPara>& cameras,
	double& agx, double& agy, double& agz){

	CReadPosBase* posdata = new CReadCarPos();
	posdata->ReadPOSData(posfile);

	char** filenames = NULL;
	int n = 0;
	int nfile = 0;

	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 0);
	printf("%d \n", nfile);

	filenames = f2c(nfile, 256);
	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 1);

	if (nfile == 0)
	{
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 0);
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 1);
	}

	printf("image number: %d \n", nfile);


	int nrange = nImageNumber;
	//int startindex = 0;
	//int endindex = min(startindex + nrange, nfile - 1);
	int endindex = startindex + nrange;

	char** selectedFiles = f2c(endindex - startindex, 256);
	int id = 0;
	for (int i = startindex; i < endindex; i++)
	{
		strcpy(selectedFiles[id], filenames[i]);
		id++;
	}
	int nimage = endindex - startindex;

	//set the camera parameters based on POS
	int ht = 4096;
	int wd = 8192;

	cameras.resize(nimage);

	srand(NULL);

	double minx = 1000000000;
	double maxx = -100000000;
	double miny = 1000000000;
	double maxy = -1000000000;
	double minz = 1000000000;
	double maxz = -1000000000;
	for (int i = 0; i < nimage; i++)
	{

		double rv = (rand() / double(RAND_MAX) - 0.5) * 2;

		POSInfo pi;
		for (int k = 0; k < posdata->GetPOSNum(); k++)
		{
			posdata->GetPOS(k, pi);
			if (strstr(selectedFiles[i], pi.filetitle) != NULL)
				break;
		}

		strcpy(cameras[i].title, pi.filetitle);

		cameras[i].rows = ht;
		cameras[i].cols = wd;
		cameras[i].camtype = camType;

		agx += pi.gx;
		agy += pi.gy;
		agz += pi.height;

		if (minx > pi.gx) minx = pi.gx;
		if (maxx < pi.gx) maxx = pi.gx;
		if (miny > pi.gy) miny = pi.gy;
		if (maxy < pi.gy) maxy = pi.gy;
		if (minz > pi.height) minz = pi.height;
		if (maxz < pi.height) maxz = pi.height;

		cameras[i].T[0] = pi.gx;
		cameras[i].T[1] = pi.gy;
		cameras[i].T[2] = pi.height;

		cameras[i].gx = pi.gx;
		cameras[i].gy = pi.gy;
		cameras[i].gz = pi.height;

		cameras[i].lon = pi.lon;
		cameras[i].lat = pi.lat;

		cameras[i].yaw = pi.yaw;
		cameras[i].pitch = -pi.pitch; //for orbit
		cameras[i].roll =  -pi.roll;  //for orbit

		cameras[i].bIsExplicit = true;

		eular2rot(cameras[i].R, cameras[i].pitch, cameras[i].roll, cameras[i].yaw);
	}

	agx /= double(nimage);
	agy /= double(nimage);
	agz /= double(nimage);

	double wx = maxx - minx;
	double wy = maxy - miny;
	double wz = maxz - minz;

	wx = 1;
	wy = 1;
	wz = 1;

	printf("%lf %lf %lf \n", agx, agy, agz);

	for (int i = 0; i < cameras.size(); i++)
	{
		/*cameras[i].gx -= agx;
		cameras[i].gy -= agy;
		cameras[i].gz -= agz;*/

		cameras[i].T[0] -= agx;
		cameras[i].T[1] -= agy;
		cameras[i].T[2] -= agz;

		cameras[i].T[0] /= wx;
		cameras[i].T[1] /= wy;
		cameras[i].T[2] /= wz;
	}

	return 0;
}



Mat thinImage(const cv::Mat & src, const int maxIterations = -1)  
{
   assert(src.type() == CV_8UC1);  
   cv::Mat dst;  
   int width = src.cols;  
   int height = src.rows;  
   src.copyTo(dst);  
   int count = 0;  //记录迭代次数    
   while (true)  
   {  
       count++;  
       if (maxIterations != -1 && count > maxIterations) //限制次数并且迭代次数到达    
           break;  
       std::vector<uchar *> mFlag; //用于标记需要删除的点    
       //对点标记    
       for (int i = 0; i < height; ++i)  
       {  
           uchar * p = dst.ptr<uchar>(i);  
           for (int j = 0; j < width; ++j)  
           {  
               //如果满足四个条件，进行标记    
               //  p9 p2 p3    
               //  p8 p1 p4    
               //  p7 p6 p5    
               uchar p1 = p[j];  
               if (p1 != 1) continue;  
               uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);  
               uchar p8 = (j == 0) ? 0 : *(p + j - 1);  
               uchar p2 = (i == 0) ? 0 : *(p - dst.step + j);  
               uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - dst.step + j + 1);  
               uchar p9 = (i == 0 || j == 0) ? 0 : *(p - dst.step + j - 1);  
               uchar p6 = (i == height - 1) ? 0 : *(p + dst.step + j);  
               uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + dst.step + j + 1);  
               uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + dst.step + j - 1);  
               if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)  
               {  
                   int ap = 0;  
                   if (p2 == 0 && p3 == 1) ++ap;  
                   if (p3 == 0 && p4 == 1) ++ap;  
                   if (p4 == 0 && p5 == 1) ++ap;  
                   if (p5 == 0 && p6 == 1) ++ap;  
                   if (p6 == 0 && p7 == 1) ++ap;  
                   if (p7 == 0 && p8 == 1) ++ap;  
                   if (p8 == 0 && p9 == 1) ++ap;  
                   if (p9 == 0 && p2 == 1) ++ap;  
 
                   if (ap == 1 && p2 * p4 * p6 == 0 && p4 * p6 * p8 == 0)  
                   {  
                       //标记    
                       mFlag.push_back(p + j);  
                   }  
               }  
           }  
       }  
 
       //将标记的点删除    
       for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)  
       {  
           **i = 0;  
       }  
 
       //直到没有点满足，算法结束    
       if (mFlag.empty())  
       {  
           break;  
       }  
       else  
       {  
           mFlag.clear();//将mFlag清空    
       }  
 
       //对点标记    
       for (int i = 0; i < height; ++i)  
       {  
           uchar * p = dst.ptr<uchar>(i);  
           for (int j = 0; j < width; ++j)  
           {  
               //如果满足四个条件，进行标记    
               //  p9 p2 p3    
               //  p8 p1 p4    
               //  p7 p6 p5    
			   uchar p1 = p[j];
			   if (p1 != 1) continue;
			   uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);
			   uchar p8 = (j == 0) ? 0 : *(p + j - 1);
			   uchar p2 = (i == 0) ? 0 : *(p - dst.step + j);
			   uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - dst.step + j + 1);
			   uchar p9 = (i == 0 || j == 0) ? 0 : *(p - dst.step + j - 1);
			   uchar p6 = (i == height - 1) ? 0 : *(p + dst.step + j);
			   uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + dst.step + j + 1);
			   uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + dst.step + j - 1);
  
                if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 
					&& (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)  
                {  
                    int ap = 0;  
                    if (p2 == 0 && p3 == 1) ++ap;  
                    if (p3 == 0 && p4 == 1) ++ap;  
                    if (p4 == 0 && p5 == 1) ++ap;  
                    if (p5 == 0 && p6 == 1) ++ap;  
                    if (p6 == 0 && p7 == 1) ++ap;  
                    if (p7 == 0 && p8 == 1) ++ap;  
                    if (p8 == 0 && p9 == 1) ++ap;  
                    if (p9 == 0 && p2 == 1) ++ap;  
  
                    if (ap == 1 && p2 * p4 * p8 == 0 && p2 * p6 * p8 == 0)  
                    {  
                        //标记    
                        mFlag.push_back(p + j);  
                    }  
                }  
            }  
        }  
  
        //将标记的点删除    
        for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)  
        {  
            **i = 0;  
        }  
  
        //直到没有点满足，算法结束    
        if (mFlag.empty())  
        {  
            break;  
        }  
        else  
        {  
            mFlag.clear();//将mFlag清空    
        }  
    }  
    return dst;  
}  

/** 
* @brief 对骨骼化图数据进行过滤，实现两个点之间至少隔一个空白像素 
* @param thinSrc为输入的骨骼化图像,8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白 
*/  
void filterOver(cv::Mat thinSrc)  
{  
    assert(thinSrc.type() == CV_8UC1);  
    int width = thinSrc.cols;  
    int height = thinSrc.rows;  
    for (int i = 0; i < height; ++i)  
    {  
        uchar * p = thinSrc.ptr<uchar>(i);  
        for (int j = 0; j < width; ++j)  
        {  
            // 实现两个点之间至少隔一个像素  
            //  p9 p2 p3    
            //  p8 p1 p4    
            //  p7 p6 p5    
            uchar p1 = p[j];  
            if (p1 != 1) continue;  
            uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);  
            uchar p8 = (j == 0) ? 0 : *(p + j - 1);  
            uchar p2 = (i == 0) ? 0 : *(p - thinSrc.step + j);  
            uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - thinSrc.step + j + 1);  
            uchar p9 = (i == 0 || j == 0) ? 0 : *(p - thinSrc.step + j - 1);  
            uchar p6 = (i == height - 1) ? 0 : *(p + thinSrc.step + j);  
            uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + thinSrc.step + j + 1);  
            uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + thinSrc.step + j - 1);  
            if (p2 + p3 + p8 + p9 >= 1)  
            {  
                p[j] = 0;  
            }  
        }  
    }  
}  

/** 
* @brief 从过滤后的骨骼化图像中寻找端点和交叉点 
* @param thinSrc为输入的过滤后骨骼化图像,8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白 
* @param raudis卷积半径，以当前像素点位圆心，在圆范围内判断点是否为端点或交叉点 
* @param thresholdMax交叉点阈值，大于这个值为交叉点 
* @param thresholdMin端点阈值，小于这个值为端点 
* @return 为对src细化后的输出图像,格式与src格式相同，元素中只有0与1,1代表有元素，0代表为空白 
*/  
std::vector<cv::Point> getPoints(const cv::Mat &thinSrc, unsigned int raudis = 4, unsigned int thresholdMax = 6, unsigned int thresholdMin = 4)  
{  
    assert(thinSrc.type() == CV_8UC1);  
    int width = thinSrc.cols;  
    int height = thinSrc.rows;  
    cv::Mat tmp;  
    thinSrc.copyTo(tmp);  
    std::vector<cv::Point> points;  
    for (int i = 0; i < height; ++i)  
    {  
        for (int j = 0; j < width; ++j)  
        {  
            if (*(tmp.data + tmp.step * i + j) == 0)  
            {  
                continue;  
            }  
            int count=0;  
            for (int k = i - raudis; k < i + raudis+1; k++)  
            {  
                for (int l = j - raudis; l < j + raudis+1; l++)  
                {  
                    if (k < 0 || l < 0||k>height-1||l>width-1)  
                    {  
                        continue;  
                          
                    }  
                    else if (*(tmp.data + tmp.step * k + l) == 1)  
                    {  
                        count++;  
                    }  
                }  
            }  
  
            if (count > thresholdMax||count<thresholdMin)  
            {  
                Point point(j, i);  
                points.push_back(point);  
            }  
        }  
    }  
    return points;  
}  


void drawLine(cv::Mat &image, double theta, double rho, cv::Scalar color)
{
	if (theta < PI / 4. || theta > 3.*PI / 4.)// ~vertical line
	{
		cv::Point pt1(rho / cos(theta), 0);
		cv::Point pt2((rho - image.rows * sin(theta)) / cos(theta), image.rows);
		cv::line(image, pt1, pt2, cv::Scalar(255), 1);
	}
	else
	{
		cv::Point pt1(0, rho / sin(theta));
		cv::Point pt2(image.cols, (rho - image.cols * cos(theta)) / sin(theta));
		cv::line(image, pt1, pt2, color, 1);
	}
}

typedef struct stSegment
{
	int left, right;
	int row;
};

int TestGrabcut()
{

	//the object rectangle
	cv::Rect rect;
	rect.x = 0;
	rect.y = 140;
	rect.width  = 720;
	rect.height = 180;

	//string filename = "C:\\Work\\Data\\deeplearning\\jyz\\quexian\\qx2.jpg";
	//string filename = "C:\\Work\\Data\\segmentation\\qx97-crop.jpg";
	//string filename = "C:\\Work\\Data\\segmentation\\qx3-crop.jpg";
	//string filename = "C:\\Work\\Data\\segmentation\\qx59-crop.jpg";
	//string filename = "C:\\Work\\Data\\segmentation\\nw000522-crop.jpg";
	string filename = "C:\\Work\\Data\\segmentation\\nw001449-crop.jpg";



	//detect the jyz object


	
	//
	Mat Img = imread(filename);
	blur(Img, Img, Size(3, 3));

	
	
	//resize 



	//grab
	//rect.x = 0;
	//rect.y = 0;
	//rect.width  = Img.cols-1;
	//rect.height = Img.rows-1;
	//cv::Mat bgModel, fgModel, fgmask;
	//cv::grabCut(Img, fgmask, rect, bgModel, fgModel, 3, cv::GC_INIT_WITH_RECT);
	//cv::compare(fgmask, cv::GC_PR_FGD, fgmask, cv::CMP_EQ);
	//imwrite("c:\\temp\\grabcut-mask.jpg", fgmask);


	//canny edge
	//Mat edge, grayImage;
	//cvtColor(Img, grayImage, COLOR_BGR2GRAY);
	//blur(grayImage, edge, Size(3, 3));
	//Canny(edge, edge, 16, 128, 3);
	//imwrite("c:\\temp\\canny.jpg", edge);



	/*
	//otsu
	Mat imageOtsu;
	//threshold(grayImage, binImg, 128, 1, cv::THRESH_BINARY);
	//threshold(grayImage, imageOtsu, 0, 255, CV_THRESH_OTSU); //Opencv Otsu算法 
	adaptiveThreshold(grayImage, imageOtsu, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 41, 0);
	imwrite("c:\\temp\\otsu.jpg", imageOtsu);

	Mat binImg = imageOtsu/255;
	Mat erodeStruct = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
	erode(binImg, binImg, erodeStruct);
	dilate(binImg, binImg, erodeStruct);
	Mat dilateStruct = getStructuringElement(MORPH_RECT, Size(5, 1));
	dilate(binImg, binImg, dilateStruct);
	imwrite("c:\\temp\\binimg.jpg", binImg*255);


	//skeleton
	Mat skeImg = thinImage(binImg);
	imwrite("c:\\temp\\skeleton.jpg", skeImg*255);
	*/

	//segmentation using Hue
	//hue
	int hmin = 150;
	int hmin_Max = 360;
	int hmax = 200;
	int hmax_Max = 360;
	//light  
	int lmin = 50;
	int lmin_Max = 255;
	int lmax = 255;
	int lmax_Max = 255;
	//sature  
	int smin = 0;
	int smin_Max = 255;
	int smax = 255;
	int smax_Max = 255;

	Mat bgr,hls;
	Img.convertTo(bgr, CV_32FC3, 1.0 / 255, 0);
	cvtColor(bgr, hls, COLOR_BGR2HLS);

	//Mat dst = Mat::zeros(bgr.size(), CV_32FC3);
	Mat mask;
	inRange(hls, Scalar(hmin, lmin / float(lmin_Max), smin / float(smin_Max)), 
		         Scalar(hmax, lmax / float(lmax_Max), smax / float(smax_Max)), 
				 mask);

	/*for (int r = 0; r < bgr.rows; r++)
	{
		for (int c = 0; c < bgr.cols; c++)
		{
			if (mask.at<uchar>(r, c) == 255)
			{
				dst.at<Vec3f>(r, c) = bgr.at<Vec3f>(r, c);
			}
		}
	}
	dst.convertTo(dst, CV_8UC3, 255, 0);*/
	imwrite("c:\\temp\\seg-color.jpg", mask);  

	
	//linear regression
	vector<Point> forePoints;
	for (int r = 0; r < bgr.rows; r++)
	{
		for (int c = 0; c < bgr.cols; c++)
		{
			if (mask.at<uchar>(r, c) == 255)
			{
				Point p;
				p.x = c;
				p.y = r;
				forePoints.push_back(p);
			}
		}
	}
	cv::Vec4f line;
	cv::fitLine(forePoints,
		line,
		CV_DIST_HUBER,
		0,
		0.01,
		0.01);

	double cos_theta = line[0];
	double sin_theta = line[1];
	double x0 = line[2], y0 = line[3];
	double phi = atan2(sin_theta, cos_theta) + PI / 2.0;
	double rho = y0 * cos_theta - x0 * sin_theta;
	std::cout << "phi = " << phi / PI * 180 << std::endl;
	std::cout << "rho = " << rho << std::endl;
	//drawLine(Img, phi, rho, cv::Scalar(0));
	//double k = sin_theta / cos_theta;
	//imwrite("c:\\temp\\line.jpg", Img);


	//rotation
	double angle = phi/PI*180 - 90;
	cv::Point2f center(mask.cols / 2, mask.rows / 2);
	cv::Mat rot = cv::getRotationMatrix2D(center, angle, 1);
	cv::Rect bbox = cv::RotatedRect(center, mask.size(), angle).boundingRect();
	rot.at<double>(0, 2) += bbox.width / 2.0 - center.x;
	rot.at<double>(1, 2) += bbox.height / 2.0 - center.y;
	cv::Mat dst;
	cv::warpAffine(mask, dst, rot, bbox.size(), INTER_NEAREST);
	Mat rotateImg;
	cv::warpAffine(Img, rotateImg, rot, bbox.size());
	

	//morphology
	Mat erodeStruct = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
	erode(dst, dst, erodeStruct);
	dilate(dst, dst, erodeStruct);
	Mat dilateStruct = getStructuringElement(MORPH_RECT, Size(1, 5));
	dilate(dst, dst, dilateStruct);
	imwrite("c:\\temp\\rotate.jpg", dst);

	//scan the right edge point 
	vector<vector<int>> allCross;
	vector<vector<int>> crossMat;
	vector<int> segNumber;
	for (int r = 0; r < dst.rows; r++)
	{
		uchar *p = dst.ptr<uchar>(r);
		vector<int> lineCross;
		for (int c = 1; c < dst.cols-1; c++)
		{
			if (p[c] == 0)
				continue;
			if ( p[c - 1] == 0 && p[c]>0 )
			{
				lineCross.push_back(c);				
			}
		}
		allCross.push_back(lineCross);
		if (lineCross.size() < 3)
			continue;
		crossMat.push_back(lineCross);
		segNumber.push_back( lineCross.size() );
	}
	sort(segNumber.begin(), segNumber.end());
	int medianSegNumber = segNumber[ segNumber.size()/2 ];
	printf("median seg: %d \n", medianSegNumber);

	//calculate the distance of segment
	vector<vector<int>> crossDis;
	vector<int> segDis;
	for (int j = 0; j < crossMat.size(); j++)
	{
		if (crossMat[j].size()<(0.5*medianSegNumber) || crossMat[j].size() >(1.5*medianSegNumber))
			continue;

		vector<int> lineCrossDis;
		for (int i = 0; i < crossMat[j].size() - 1; i++)
		{
			int dis = crossMat[j][i + 1] - crossMat[j][i];
			lineCrossDis.push_back(dis);
			//printf("%d ", dis);
			segDis.push_back(dis);
		}
		//printf("\n");
		crossDis.push_back(lineCrossDis);
	}
	sort(segDis.begin(), segDis.end());
	int medianDis = segDis[segDis.size() / 2];
	double thresh = segDis[segDis.size() / 2] * 1.5;
	printf("medianDis  thresh: %d  %lf \n", medianDis, thresh);

	//generate histogram
	int nHist = rotateImg.cols / medianDis + 1;
	vector<double> lossHist(nHist, 0);
	vector<stSegment> validSeg;
	for (int j = 0; j < allCross.size(); j++)
	{
		int ncp = allCross[j].size();

		//1.segment constraint
		if (ncp < (medianSegNumber*0.5) || ncp>(medianSegNumber*1.5))
			continue;

		for (int i = 0; i < ncp - 1; i++)
		{
			int dis = allCross[j][i + 1] - allCross[j][i];

			//2.distance constrait
			if (dis>thresh && dis<thresh*3)
			{
				//in the loss area, background is more than foreground
				int fgSum = 0;
				int bgSum = 0;
				uchar *dp = dst.ptr<uchar>(j);
				for (int k = allCross[j][i]; k < allCross[j][i + 1]; k++)
				{
					if (dp[k]>0)
						fgSum++;
					else
						bgSum++;

				}

				//3.background & foreground constaint
				if (bgSum > fgSum && fgSum>4 )
				{
					//generate histogram
					double hv = double(allCross[j][i] + allCross[j][i + 1])*0.5 / double(medianDis);
					int    ih = hv;
					if (ih >= (nHist - 1)) ih = nHist - 2;
					double dh = hv - ih;
					lossHist[ih] += (1 - dh);
					lossHist[ih + 1] += dh;

					stSegment seg;
					seg.left = allCross[j][i];
					seg.right = allCross[j][i+1];
					seg.row = j;
					validSeg.push_back(seg);
				}
			}
		}
	}

	printf("histogram... \n");
	for (int i = 0; i < lossHist.size(); i++)
	{
		printf("%lf ", lossHist[i]);
	}
	printf("\n");

	//determin the loss area by histogram
	vector<Rect> lossRects;
	for (int i = 0; i < lossHist.size()-1; i++)
	{
		if ( (lossHist[i]+lossHist[i+1])>medianDis)
		{
			Rect rec;
			int id = 0;
			for (int k = 0; k < validSeg.size(); k++)
			{
				double cx = double(validSeg[k].left + validSeg[k].right)*0.5;
				double dIndex = cx / double(medianDis);

				if (fabs(dIndex - i - 0.5) < 1)
				{
					if (id == 0){
						rec.x = validSeg[k].left;
						rec.y = validSeg[k].row;
						rec.width = validSeg[k].right - validSeg[k].left;
						rec.height = 1;
					}
					else{
						Rect curRec;
						curRec.x = validSeg[k].left;
						curRec.y = validSeg[k].row;
						curRec.width = validSeg[k].right - validSeg[k].left;
						curRec.height = 1;

						//cluster
						rec = rec | curRec;
					}
					id++;
				}
			}
			lossRects.push_back(rec);
		}
	}

	printf("detect loss number: %d \n", lossRects.size());
	for (int i = 0; i < lossRects.size(); i++)
	{
		rectangle(rotateImg, lossRects[i], CV_RGB(255, 0, 0), 2);
	}

	imwrite("c:\\temp\\detecion.jpg", rotateImg);


	//Mat lossObj(dst.rows, dst.cols, CV_8UC1, Scalar::all(0));
	//for (int j = 0; j < allCross.size(); j++)
	//{
	//	int ncp = allCross[j].size();

	//	//the line is invalid when its segments are few or much
	//	if (ncp < (medianSegNumber*0.5) || ncp>(medianSegNumber*1.5))
	//		continue;

	//	for (int i = 0; i < ncp - 1; i++)
	//	{
	//		int dis = allCross[j][i + 1] - allCross[j][i];
	//		if (dis>thresh)
	//		{
	//			//in the loss area, background is more than foreground
	//			int fgSum = 0;
	//			int bgSum = 0;
	//			uchar *dp = dst.ptr<uchar>(j);
	//			for (int k = allCross[j][i]; k < allCross[j][i + 1]; k++)
	//			{
	//				if (dp[k]>0)
	//					fgSum++;
	//				else
	//					bgSum++;
	//			}

	//			uchar *p = lossObj.ptr<uchar>(j);
	//			if ( bgSum>fgSum && fgSum>4)
	//			{
	//				//
	//				double hv = double(allCross[j][i] + allCross[j][i + 1])*0.5 / double(medianDis);
	//				int    ih = hv;
	//				if (ih >= (nHist - 1)) ih = nHist - 2;
	//				double dh = hv - ih;

	//				int lossHit =   lossHist[ih]*(1-dh) + lossHist[ih+1]*dh;

	//				//the probability
	//				if (lossHit < (medianDis*0.5) )
	//					continue;

	//				//fill the area
	//				for (int k = allCross[j][i]; k < allCross[j][i + 1]; k++)
	//				{
	//					p[k] = 255;
	//				}

	//				//draw the cross point
	//				Point p1, p2;
	//				p1.x = allCross[j][i];
	//				p1.y = j;
	//				p2.x = allCross[j][i + 1];
	//				p2.y = j;
	//				circle(rotateImg, p1, 1, CV_RGB(255, 0, 0));
	//				circle(rotateImg, p2, 1, CV_RGB(255, 0, 0));
	//			}
	//		}
	//	}
	//}


	//detect based on hist


	/*erodeStruct = getStructuringElement(MORPH_ELLIPSE, Size(2, 2));
	erode(lossObj, lossObj, erodeStruct);
	dilate(lossObj, lossObj, erodeStruct);
	dilateStruct = getStructuringElement(MORPH_RECT, Size(1, 5));
	dilate(lossObj, lossObj, dilateStruct);*/


	////foreground cluster
	//MyRect* pRect = (MyRect*)malloc(512 * sizeof(MyRect));
	//int  nRect = 0;
	//SegLabel(lossObj.data, lossObj.rows, lossObj.cols, pRect, &nRect);

	//for (int i = 0; i < nRect; i++)
	//{
	//	Rect foreRect;
	//	foreRect.x = pRect[i].left;
	//	foreRect.y = pRect[i].top;
	//	foreRect.width = pRect[i].right - pRect[i].left;
	//	foreRect.height = pRect[i].bottom - pRect[i].top;

	//	printf("rect height: %d  %d \n", foreRect.height, int(medianDis*0.5));

	//	if (foreRect.height > medianDis*0.5)
	//	{
	//		rectangle(rotateImg, foreRect, CV_RGB(255,0,0));
	//	}
	//}

	//free(pRect);

	//imwrite("c:\\temp\\jyz-loss.jpg", lossObj);

	/*printf("histogram... \n");
	for (int i = 0; i < lossHist.size(); i++)
	{
		printf("%lf ", lossHist[i]);
	}
	printf("\n");*/

	return 0;
}


int TestPano()
{
	char* filepath = "C:\\Work\\Data\\Panorama\\2016.10.13-yizhuang\\L10_1013\\indoor\\jpg\\ladybug_panoramic_000000.jpg";

	SphereToCilinder(filepath, "c:\\temp\\cilinder.jpg");

	return 0;
}

int TestGeotiff()
{
	//char* filepath = "c:\\work\\data\\images\\timg.jpg";
	//PrintImageInfo(filepath);

	test_gdal();

	return 0;
}


//radiative transfer validation for pure molecular, 2018.1.30, xiedonghai
int TestRT()
{
	double w = 0.9;
	double g = 0.7;
	double aod = 0.0;
	double surfRef = 0.0;

	//rayleight optical depth
    double wav = 0.47; // 0.543;
	double rayleighOd = RayleighOD(wav, 0);
	printf("rayleigh aod: %lf %lf \n", wav, rayleighOd);

	//geometry condition
    double sza = 64.69; // 36.8699;   // 23.0739;
    double az = 51; // 180;       // 180;
	double vza[16] = { 0.0000, 11.4783, 16.2602, 23.0739, 32.8599, 43.9455, 50.2082, 58.6677,
		66.4218, 71.3371, 73.7398, 78.4630, 80.7931, 84.2608, 86.5602, 88.8540 };

	printf("wav: %lf  surf: %lf  sza: %lf  az: %lf  raylei_od: %lf  aod: %lf \n",
		wav, surfRef, sza, az, rayleighOd, aod);

	//single scattering for rayleight
	//interpolation based on LUT
    char*  lutpath0  = "C:\\Work\\Data\\luts-validation\\lut_threesurface\\2\\blue\\0.dat";// "C:\\Work\\Data\\luts-validation\\Double_543_208_619_130_2241_531_1039_1537_2_0.dat";
    char*  lutpath10 = "C:\\Work\\Data\\luts-validation\\lut_threesurface\\2\\blue\\10.dat";;// "C:\\Work\\Data\\luts-validation\\Double_543_208_619_130_2241_531_1039_1537_2_10.dat";
    char*  lutpath25 = "C:\\Work\\Data\\luts-validation\\lut_threesurface\\2\\blue\\25.dat";;// "C:\\Work\\Data\\luts-validation\\Double_543_208_619_130_2241_531_1039_1537_2_25.dat";
	CAtmosphereCorrectionRt3 cac;
	cac.LoadLut(lutpath0, lutpath10, lutpath25, wav);

	CLutBinCont lut;
	lut.Load(lutpath0);
	//printf("Rayleigh: lut interpolation \n");

	FILE* fp = fopen("c:\\temp\\rt-compare-surf-0.25.txt", "w");
	for (int i = 0; i < 16; i++)
	{
		//1. single scattering simulation
		//mocular contribution
		double pf = RayleighPhaseFunction(sza, vza[i], 0, az);
		double rayleighContribution = rayleighOd*pf / (4 * cos(sza*DPI)*cos(vza[i] * DPI));

		//surface contribution
		double tdown = Transmission(sza, rayleighOd, 0.0);
		double tup = Transmission(vza[i], rayleighOd, 0.0);
		double bs = BackScatteringRatio(rayleighOd, 0.0, 0.7);
		double surfContribution = tdown*tup*surfRef / (1 - surfRef*bs);

		//aerosol contribution
		double aerosolPf = AerosolFunction(sza, vza[i], 0, az, g);
		double aerosolContribution = w*aod*aerosolPf / (4 * cos(sza*DPI)*cos(vza[i] * DPI));

		double singleSimulatedValue = rayleighContribution + aerosolContribution + surfContribution;

		//calculate the T and S based on lut
		double T, S;
		cac.CalculateTAndS(sza, vza[i], az, 0.0, T, S);
		//printf("T and S -- single: %lf %lf  lut: %lf %lf \n", tdown*tup, bs, T, S);


		//2.lut interpolation based on rt3 lut
		double r, rp;
		lut.GetValue(sza, vza[i], az, 0, r, rp);
		double lutSim = r + T*surfRef / (1 - surfRef*S);


		//
		printf("%lf %lf %lf  %lf \n", vza[i], singleSimulatedValue, lutSim);
		fprintf(fp, "%lf %lf \n", singleSimulatedValue, lutSim);

	}
	fclose(fp);

	////rayleigh function
	//for (int i = 0; i <= 180; i++)
	//{
	//	double pf = RayleighPhaseFunction(double(i));
	//	printf("%d %lf \n", i, pf);
	//}
	return 0;
}


//int SimulateParasol(double inputaod, stMultiView mdata)
//{
//
//
//	return 0;
//}



int TestRTUsingParasol(double inputaod)
{
	double wav = 0.865;
	double aod = inputaod;

	printf("input aod: %lf \n", aod);
		
	/*char*  lutpath0  = "C:\\Work\\Data\\Luts\\EastAsia-Industry1\\Double_865_219_531_153_2724_583_131_1485_8_0.dat";
	char*  lutpath10 = "C:\\Work\\Data\\Luts\\EastAsia-Industry1\\Double_865_219_531_153_2724_583_131_1485_8_10.dat";
	char*  lutpath25 = "C:\\Work\\Data\\Luts\\EastAsia-Industry1\\Double_865_219_531_153_2724_583_131_1485_8_25.dat";*/
	
	/*char*  lutpath0  = "C:\\Work\\Data\\Luts\\EastAsia-Industry2\\Double_865_257_535_269_2580_568_192_1483_7_0.dat";
	char*  lutpath10 = "C:\\Work\\Data\\Luts\\EastAsia-Industry2\\Double_865_257_535_269_2580_568_192_1483_7_10.dat";
	char*  lutpath25 = "C:\\Work\\Data\\Luts\\EastAsia-Industry2\\Double_865_257_535_269_2580_568_192_1483_7_25.dat";*/

	/*char*  lutpath0  = "C:\\Work\\Data\\Luts\\EastAsia-Rural1\\Double_865_192_504_79_2915_618_75_1468_10_0.dat";
	char*  lutpath10 = "C:\\Work\\Data\\Luts\\EastAsia-Rural1\\Double_865_192_504_79_2915_618_75_1468_10_10.dat";
	char*  lutpath25 = "C:\\Work\\Data\\Luts\\EastAsia-Rural1\\Double_865_192_504_79_2915_618_75_1468_10_25.dat";
	*/
	//char*  lutpath0  = "C:\\Work\\Data\\Luts\\Haze-Beijing-Fine-NA\\Double_869_219_510_127_2731_602_80_1434_5_0.dat";
	//char*  lutpath10 = "C:\\Work\\Data\\Luts\\Haze-Beijing-Fine-NA\\Double_869_219_510_127_2731_602_80_1434_5_10.dat";
	//char*  lutpath25 = "C:\\Work\\Data\\Luts\\Haze-Beijing-Fine-NA\\Double_869_219_510_127_2731_602_80_1434_5_25.dat";

	//char*  lutpath0  = "C:\\Work\\Data\\Luts\\Haze-Beijing-Fine-MA\\Double_869_183_498_103_2673_625_112_1491_10_0.dat";
	//char*  lutpath10 = "C:\\Work\\Data\\Luts\\Haze-Beijing-Fine-MA\\Double_869_183_498_103_2673_625_112_1491_10_10.dat";
	//char*  lutpath25 = "C:\\Work\\Data\\Luts\\Haze-Beijing-Fine-MA\\Double_869_183_498_103_2673_625_112_1491_10_25.dat";

	//2015.12.09
	//char*  lutpath0  = "C:\\Work\\Data\\luts-validation\\20051209\\Double_870_212_565_111_3078_653_96_1576_20_0.dat";
	//char*  lutpath10 = "C:\\Work\\Data\\luts-validation\\20051209\\Double_870_212_565_111_3078_653_96_1576_20_10.dat";
	//char*  lutpath25 = "C:\\Work\\Data\\luts-validation\\20051209\\Double_870_212_565_111_3078_653_96_1576_20_25.dat";

	//
	char*  lutpath0  = "C:\\Work\\Data\\luts-validation\\20091028\\Double_870_279_540_117_2410_551_148_1538_0_0.dat";
	char*  lutpath10 = "C:\\Work\\Data\\luts-validation\\20091028\\Double_870_279_540_117_2410_551_148_1538_0_10.dat";
	char*  lutpath25 = "C:\\Work\\Data\\luts-validation\\20091028\\Double_870_279_540_117_2410_551_148_1538_0_25.dat";


	CAtmosphereCorrectionRt3 cac;
	cac.LoadLut(lutpath0, lutpath10, lutpath25, wav);

	//bpdf model
	double alfa = 0.025;
	double beita = 45;

	//char* parasolFile = "C:\\Work\\Data\\parasol\\Haze\\2005-12-29\\P3L1TBG1025076JD";
	char* parasolFile = "C:\\Work\\Data\\parasol\\Haze\\2009-10-28\\P3L1TBG1112178KD";
	//char* parasolFile = "C:\\Work\\Data\\parasol\\Haze\\2010-11-17\\P3L1TBG1136197KD";


	CReasParasolLevel1File  readParasol;
	readParasol.Load(parasolFile);

	//39.977N, 116.381E
	double lon = 116.381;
	double lat = 39.977;
	double surfRef = 0.1;

	stMultiView mdata = readParasol.GetPtMultiViewData(lon, lat);

	double errorp = 0;
	double error = 0;
	int nsum = 0;

	double sumrp = 0;

	//FILE* fp = fopen("c:\\temp\\parasol-2005-12-29.txt", "w");
	FILE* fp = fopen("c:\\temp\\parasol-2009-10-28.txt", "w");
	//FILE* fp = fopen("c:\\temp\\parasol-2010-11-17.txt", "w");
	for (int i = 0; i < mdata.nView; i++)
	{
		if (mdata.mViewData[i].rp8 < 0)
			continue;

		double sza = mdata.mViewData[i].sza;
		double vza = mdata.mViewData[i].vza;
		double azi = mdata.mViewData[i].azi;
		double sca = mdata.mViewData[i].sca;

		//if (sca>140)
		//	continue;

		//simulate the data
		double T, S;
		cac.CalculateTAndS(sza, vza, azi, aod, T, S);
		double r, rp;
		cac.GetPathRef(sza, vza, azi, aod, r, rp);

		//polarized reflectance
		double surfPR = BPDF_Nadal(sza, vza, azi, alfa, beita);
		double simulateTOAP = rp + T*surfPR / (1 - S*surfPR);

		//reflectance
		double simulateTOA = r + T*surfRef / (1 - S*surfRef);

		printf("%lf %lf %lf %lf   %lf  %lf  %lf  %lf  \n",
			sza, vza, azi, sca, simulateTOA, mdata.mViewData[i].r8, simulateTOAP, mdata.mViewData[i].rp8);

		fprintf(fp, "%lf %lf %lf %lf    %lf %lf  %lf %lf  \n",
			sza, vza, azi, sca, simulateTOA, mdata.mViewData[i].r8, simulateTOAP, mdata.mViewData[i].rp8);

		//error += fabs(simulateTOA - mdata.mViewData[i].rp8) / ((simulateTOA + mdata.mViewData[i].rp8)*0.5);
		/*double ep1 = fabs(simulateTOAP - mdata.mViewData[i].rp8) / ((simulateTOAP + mdata.mViewData[i].rp8)*0.5);
		printf("%lf \n", ep1);
		errorp += fabs(simulateTOAP - mdata.mViewData[i].rp8) / ((simulateTOAP + mdata.mViewData[i].rp8)*0.5);
		error  += fabs(simulateTOA - mdata.mViewData[i].r8) / ((simulateTOA + mdata.mViewData[i].r8)*0.5);*/
		
		errorp += fabs(simulateTOAP - mdata.mViewData[i].rp8);
		sumrp += simulateTOAP;

		nsum++;
	}
	fclose(fp);

	printf("%lf \n", errorp / sumrp);
	//printf("error: %lf  %lf \n", error / nsum, errorp/nsum);

	return 0;
}


void TestRand()
{
	int m = 50000;
	int n = 16;
	for (int i = 0; i < n; ++i)
	{
		int index = rand() % m;
		printf("%d ", index);
	}
	printf("\n");
}

int invokeCeres()
{

	test_ceres();
	//test_robust_fitting();

	return 0;
}


//ceres ba, new apis, 2017.12.31
int main_realimages(int argc, char* argv[])
{
	time_t  start = time(NULL);
	//double t = (double)getTickCount();

	//set the type of camera
	//CameraType camType = PerspectiveCam; //;PanoramCam;
	CameraType camType = PanoramCam;
	printf("SFM integration .... \n");

	//perspective data
	//char imagepath[256] = "C:\\Work\\Data\\ba1";

	//panorama data
	char* imagepath = "C:\\Work\\Data\\Panorama\\beijing\\image_all";
	//char imagepath[256] = "C:\\Work\\Data\\Panorama\\20180226_test_xian";
	//char imagepath[256] = "C:\\Work\\Programs\\PgStation\\CMVS\\bin\\pmvs\\visualize";
	//char imagepath[256] = "C:\\Work\\Data\\panorama\\2016.10.13-yizhuang\\L10_1013\\indoor\\jpg";
	//char imagepath[256] = "C:\\Work\\Data\\panorama\\cheku20161011\\garageoutput";
	//char imagepath[256] = "C:\\Work\\Data\\panorama\\cheku2016.10.26\\Output";
	//char imagepath[256] = "C:\\Work\\Data\\UAV\\binchuan";
	//char imagepath[256] = "C:\\Work\\Programs\\PgStation\\Bundler_V0.4\\examples\\kermit";
	char* outpath = "C:\\Work\\Programs\\PgStation\\CMVS\\bin\\pmvs";

	char* posfile = "C:\\Work\\Data\\Panorama\\beijing\\noloss.cam";

	//if (argc == 2)
	//{
	//	strcpy(imagepath, argv[1]);
	//	printf("path: %s \n", imagepath);
	//}

	char** filenames = NULL;
	int n = 0;
	int nfile = 0;

	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 0);
	printf("%d \n", nfile);

	filenames = f2c(nfile, 256);
	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 1);

	if (nfile == 0)
	{
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 0);
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 1);
	}


	printf("image number: %d \n", nfile);

	if (nfile<2)
	{
		printf("images are less than 2 ! \n");
		return -1;
	}

	int nrange = 50;
	int startindex = 130;
	int endindex = min(startindex + nrange, nfile - 1);
	int nimage = endindex - startindex + 1;
	char** selectedFiles = f2c(nimage, 256);
	int id = 0;
	for (int i = startindex; i <= endindex; i++)
	{
		strcpy(selectedFiles[id], filenames[i]);
		id++;
	}

	vector<CameraPara> cameras;
	cameras.resize(nimage);
	double agx = 0, agy = 0, agz = 0;
	ReadPosData(posfile, imagepath, startindex, nimage, PanoramCam, cameras, agx, agy, agz);


	//1. generating the image feature points
	vector<ImgFeature> imgFeatures;
	double dProgress = 0;
	DetectFileFeaturePts(selectedFiles, nimage, imgFeatures, 900, dProgress);

	//2. matching images 
	vector<PairMatchRes> matchRes;
	MatchImageFiles(imgFeatures, matchRes, camType, 3);// imgFeatures.size());

	//3.generating the tracks
	vector<TrackInfo> tracks;
	CGenerateTracksBase* pGenerateTrack = new CFastGenerateTrack();
	pGenerateTrack->GenerateTracks(imgFeatures, matchRes, tracks);

	//output the tracks
	FILE* fp = fopen("c:\\temp\\mytrack.txt", "w");
	fprintf(fp, "%d \n", tracks.size());
	for (int i = 0; i<tracks.size(); i++)
	{
		fprintf(fp, "%d  ", tracks[i].views.size());
		for (int j = 0; j<tracks[i].views.size(); j++)
		{
			fprintf(fp, " %d %d ", tracks[i].views[j].first, tracks[i].views[j].second);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);


	//4. bundle adjustment
	int numImage = imgFeatures.size();
	
	double focalLen = (imgFeatures[0].ht + imgFeatures[0].wd) * 0.5;
	for (int i = 0; i<numImage; i++)
	{
		cameras[i].focalLen = focalLen; //initialize the focal length 
		memset(cameras[i].R, 0, sizeof(double)* 9);
		cameras[i].R[0] = 1;
		cameras[i].R[4] = 1;
		cameras[i].R[8] = 1;
		cameras[i].rows = imgFeatures[0].ht;
		cameras[i].cols = imgFeatures[0].wd;
		cameras[i].camtype = camType;
	}

	CBABase* pBA = NULL;

	if (camType == PerspectiveCam)
	{
		pBA = new CCeresBA();
	}
	else if (camType == PanoramCam)
	{
		pBA = new CPanoBA();
	}

	pBA->BundleAdjust(cameras.size(), cameras, imgFeatures, matchRes, tracks, false);

	delete pBA;

	/*
	//generate a 3D point for debug
	Point2DDouble pl, pr;
	vector<Point2DDouble> lpts, rpts;
	int ht = 4000;
	int wd = 8000;
	double radius = double(wd) / (2 * PI);
	pl.p[0] = 5130 - wd*0.5;
	pl.p[1] = ht*0.5 - 1431;
	lpts.push_back(pl);
	pr.p[0] = 5184 - wd*0.5;
	pr.p[1] = ht*0.5 - 1413;
	rpts.push_back(pr);
	CTriangulateBase* pTri = new CTriangulateCV();
	vector<Point3DDouble> gpts;
	vector<double> errorarray;
	pTri->Triangulate(lpts, rpts, cameras[0], cameras[1], gpts, errorarray);
	printf("Generating 3D point... \n");
	for (int i = 0; i < gpts.size(); i++)
	{
		double gx = gpts[i].p[0];
		double gy = gpts[i].p[1];
		double gz = gpts[i].p[2];
		printf("ground point: %lf %lf %lf \n", gpts[i].p[0], gpts[i].p[1], gpts[i].p[2]);

		double p[3];
		p[0] = gx - cameras[0].T[0];
		p[1] = gy - cameras[0].T[1];
		p[2] = gz - cameras[0].T[2];
		double gp[3];
		mult(cameras[0].R, p, gp, 3, 3, 1);

		//from 3D to panoram point
		double ix, iy;
		GrdToSphere_center(gp[0], gp[1], gp[2], radius, ix, iy);
		printf("image point: %lf %lf \n", ix + wd*0.5, ht*0.5 - iy);
	}
	*/

	//getchar();

	//output the BA  to the pmvs files, added by xiedonghai, 2018.1.21
	if (0)
	{

		//6. generate projection images, save them and the corresponding projective matrix into the files
		int nLevel = 3;
		if (camType == PerspectiveCam)
		{

			for (int k = 0; k < numImage; k++)
			{
				char projFile[256];
				sprintf(projFile, "%s\\txt\\%.8d.txt", outpath, k);
				WritePMVSCamFile(projFile, cameras[k]);

				//Mat = imread(filenames[i]);
				char imageFile[256];
				sprintf(imageFile, "%s\\visualize\\%.8d.jpg", outpath, k);
				IplImage* pimage = cvLoadImage(filenames[k]);
				cvSaveImage(imageFile, pimage);
				cvReleaseImage(&pimage);
			}
		}

		if (camType == PanoramCam)
		{
			nLevel = 1;

			vector<stCameraPara> camParasPerspective;
			bool bIsBundlerOut = false;
			double vangle = 60, hangle = 60;
			double anglestep = 60;
			double focalratio = 1;

			for (int k = 0; k < numImage; k++)
			{
				double it[3];
				mult(cameras[k].R, cameras[k].T, it, 3, 3, 1);
				for (int m = 0; m < 3; m++)
				{
					it[m] = -it[m];
				}

				PanoToPlanes(k, filenames[k], anglestep, vangle, hangle, focalratio,
					cameras[k].R, it, camParasPerspective, outpath);
			}

			//
		}

		//generate the bundler.out file 
		//WriteBundlerOutFile("bundler.out", camParas);		
		printf("generate pmvs_options.txt ... \n");
		//int nImageNum = 2 * (360/anglestep);
		//int nImageNum = 2;
		/* Write the options file */
		char buf[512];
		sprintf(buf, "%s\\pmvs_options.txt", outpath);
		FILE *f_opt = fopen(buf, "w");
		fprintf(f_opt, "level %d \n", nLevel);
		fprintf(f_opt, "csize 2\n");
		fprintf(f_opt, "threshold 0.7\n");
		fprintf(f_opt, "wsize 7\n");
		fprintf(f_opt, "minImageNum 2\n");
		fprintf(f_opt, "CPU 4\n");
		fprintf(f_opt, "setEdge 0\n");
		fprintf(f_opt, "useBound 0\n");
		fprintf(f_opt, "useVisData 0\n");
		fprintf(f_opt, "sequence -1\n");
		fprintf(f_opt, "timages -1 0 %d \n", numImage);
		fprintf(f_opt, "oimages -3\n");
		fclose(f_opt);
	}

	//t = ((double)getTickCount() - t) / getTickFrequency();
	time_t  end = time(NULL);

	printf("\n running time: %d \n", end-start);

	//distance between the camera centers
	printf("Camera Distance .... \n");
	for (int i = 0; i < cameras.size() - 1; i++){
		double gx1 = cameras[i].T[0];
		double gy1 = cameras[i].T[1];
		double gz1 = cameras[i].T[2];

		double gx2 = cameras[i + 1].T[0];
		double gy2 = cameras[i + 1].T[1];
		double gz2 = cameras[i + 1].T[2];

		double distOpt = sqrt((gx1 - gx2)*(gx1 - gx2) +
			(gy1 - gy2)*(gy1 - gy2) + (gz1 - gz2)*(gz1 - gz2));

		//gx1 = cameras[i].gx;
		//gy1 = cameras[i].gy;
		//gz1 = cameras[i].gz;

		//gx2 = cameras[i + 1].gx;
		//gy2 = cameras[i + 1].gy;
		//gz2 = cameras[i + 1].gz;

		//double distOri = sqrt((gx1 - gx2)*(gx1 - gx2) +
		//	(gy1 - gy2)*(gy1 - gy2) + (gz1 - gz2)*(gz1 - gz2));

		//printf("%.2lf %.2lf  \n", distOri, distOpt);
		printf("%.6lf  \n", distOpt);
	}

	//absolute orientation
	stAbsPOS absparas;
	AbsOriOrthogonal(absparas, cameras);

	printf("absolute orientation ... \n");
	for (int i = 0; i < cameras.size(); i++)
	{
		printf("%.4lf %.4lf %.4lf \n",
			cameras[i].gx - cameras[i].xs,
			cameras[i].gy - cameras[i].ys,
			cameras[i].gz - cameras[i].zs
			);
	}


	return 0;
}

int BAWithPos()
{
	CameraType camType = PanoramCam;

	//reading pos data of format orbit
	
	//char* posfile = "C:\\Work\\Data\\Panorama\\20180226_test_xian\\444_new_orbit.txt";
	//CReadPosBase* posdata = new CReadCarPosOrbit();
	//posdata->ReadPOSData(posfile);
	////reading image data
	//char* imagepath = "C:\\Work\\Data\\Panorama\\20180226_test_xian";
	

	char* posfile = "C:\\Work\\Data\\Panorama\\beijing\\withloss.cam";
	CReadPosBase* posdata = new CReadCarPos();
	posdata->ReadPOSData(posfile);
	char* imagepath = "C:\\Work\\Data\\Panorama\\beijing\\image_all";
	
	int nrange = 22;
	
	char** filenames = NULL;
	int n = 0;
	int nfile = 0;

	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 0);
	printf("%d \n", nfile);

	filenames = f2c(nfile, 256);
	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 1);

	if (nfile == 0)
	{
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 0);
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 1);
	}
	
	printf("image number: %d \n", nfile);

	int startindex = 0;
	int endindex = min(startindex + nrange, nfile - 1);

	char** selectedFiles = f2c(endindex - startindex + 1, 256);
	int id = 0;
	for (int i = startindex; i <= endindex; i++)
	{
		strcpy(selectedFiles[id], filenames[i]);
		id++;
	}
	int nimage = endindex - startindex + 1;
	
	//set the camera parameters based on POS
	int ht = 4096;
	int wd = 8192;
	vector<CameraPara> cameras;
	cameras.resize(nimage);
	double agx = 0;
	double agy = 0;
	double agz = 0;

	srand(NULL);

	double minx = 1000000000;
	double maxx = -100000000;
	double miny = 1000000000;
	double maxy = -1000000000;
	double minz = 1000000000;
	double maxz = -1000000000;
	for (int i = 0; i < nimage; i++)
	{

		double rv = (rand() / double(RAND_MAX) - 0.5) * 2;

		POSInfo pi;
		for (int k = 0; k < posdata->GetPOSNum(); k++)
		{
			posdata->GetPOS(k, pi);
			if (strstr(selectedFiles[i], pi.filetitle) != NULL)
				break;
		}

		strcpy(cameras[i].title, pi.filetitle);

		cameras[i].rows = ht;
		cameras[i].cols = wd;
		cameras[i].camtype = camType;

		agx += pi.gx;
		agy += pi.gy;
		agz += pi.height;

		if (minx > pi.gx) minx = pi.gx;
		if (maxx < pi.gx) maxx = pi.gx;
		if (miny > pi.gy) miny = pi.gy;
		if (maxy < pi.gy) maxy = pi.gy;
		if (minz > pi.height) minz = pi.height;
		if (maxz < pi.height) maxz = pi.height;

		cameras[i].T[0] = pi.gx;
		cameras[i].T[1] = pi.gy;
		cameras[i].T[2] = pi.height;

		cameras[i].gx = pi.gx;
		cameras[i].gy = pi.gy;
		cameras[i].gz = pi.height;

		cameras[i].yaw   = pi.yaw;
		cameras[i].pitch = -pi.pitch; //for orbit
		cameras[i].roll  = -pi.roll;  //for orbit

		cameras[i].bIsExplicit = true;

		eular2rot(cameras[i].R, cameras[i].pitch, cameras[i].roll, cameras[i].yaw);
	}

	agx /= double(nimage);
	agy /= double(nimage);
	agz /= double(nimage);

	double wx = maxx - minx;
	double wy = maxy - miny;
	double wz = maxz - minz;

	wx = 1;
	wy = 1;
	wz = 1;

	printf("%lf %lf %lf \n", agx, agy, agz);

	for (int i = 0; i < cameras.size(); i++)
	{
		cameras[i].T[0] -= agx;
		cameras[i].T[1] -= agy;
		cameras[i].T[2] -= agz;

		cameras[i].T[0] /= wx;
		cameras[i].T[1] /= wy;
		cameras[i].T[2] /= wz;
	}

	//1.generating the image feature points
	vector<ImgFeature> imgFeatures;
	double dProgress = 0;
	DetectFileFeaturePts(selectedFiles, nimage, imgFeatures, 900, dProgress);

	//2. matching images 
	vector<PairMatchRes> matchRes;
	MatchImageFiles(imgFeatures, matchRes, camType, 3);

	//3.generating the tracks
	vector<TrackInfo> tracks;
	CGenerateTracksBase* pGenerateTrack = new CFastGenerateTrack();
	pGenerateTrack->GenerateTracks(imgFeatures, matchRes, tracks);
	
	/////////////////////  bundle adjustment based on POS ////////////////////////
	CPanoBAWithPos* pBA = new CPanoBAWithPos();
	pBA->SetCenter(agx, agy, agz);
	pBA->BundleAdjust(cameras.size(), cameras, imgFeatures, tracks);
	delete pBA;

	//output the trajectory
	FILE* fp = fopen("c:\\temp\\trajectory.txt", "w");
	for (int i = 0; i < cameras.size(); i++)
	{
		fprintf(fp, "%s,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
			cameras[i].title,
			cameras[i].gx, cameras[i].gy, cameras[i].gz,
			cameras[i].yaw, cameras[i].pitch, cameras[i].roll,
			(cameras[i].T[0]*wx + agx), (cameras[i].T[1]*wy + agy), (cameras[i].T[2]*wz + agz),
			cameras[i].az, cameras[i].ax, cameras[i].ay);
	}
	fclose(fp);

	//generate simulated data, added by xdh, 2018.3.12
	int ncam = cameras.size();
	
	//vector<vector<POINT2>> imgFeaturePts;
	//imgFeaturePts.resize(ncam);

	FILE* fp1 = fopen("c:\\temp\\tracks.txt", "w");
	for (int i = 0; i < tracks.size(); i++)
	{
		int view = tracks[i].GetViews();
		if (view == 0)
			continue;

		Point3DDouble gp = tracks[i].GetGround();
		fprintf(fp1, "%lf %lf %lf   ", gp.p[0]+agx, gp.p[1]+agy, gp.p[2]+agz);
		for (int k = 0; k < view; k++){
			ImageKey key = tracks[i].GetImageKey(k);
			int imageid = key.first;
			int ptid = key.second;

			//projection
			Point3DDouble gp = tracks[i].GetGround();
			Point2DDouble ip = imgFeatures[imageid].GetCenteredPt(ptid);
			Point2DDouble projPt;
			GrdToImg(gp, projPt, cameras[imageid]);

			fprintf(fp1, "%d %d %6.2lf %6.2lf %6.2lf %6.2lf   ", imageid, ptid, 
				projPt.p[0]+wd*0.5, ht*0.5-projPt.p[1],
				ip.p[0]+wd*0.5, ht*0.5-ip.p[1]);
		}
		fprintf(fp1, "\n");
	}
	fclose(fp1);


	return 0;
}

int TestPos()
{

	printf("\n\n tesing pos ... \n");

	int ht = 4096;
	int wd = 8192;
	CameraType camType = PanoramCam;

	//reading pos data
	//char* posfile = "C:\\Work\\Data\\Panorama\\20180226_test_xian\\444_new_orbit.txt";
	//CReadPosBase* posdata = new CReadCarPosOrbit();
	//posdata->ReadPOSData(posfile);

	char* posfile = "C:\\Work\\Data\\Panorama\\beijing\\noloss.cam";
	CReadPosBase* posdata = new CReadCarPos();
	posdata->ReadPOSData(posfile);

	vector<CameraPara> cameras;
	cameras.resize(2);
	int nstart = 132;
	for (int i = 0; i < 2; i++)
	{
		POSInfo pi;
		posdata->GetPOS(nstart+i, pi);

		cameras[i].rows = ht;
		cameras[i].cols = wd;
		cameras[i].camtype = camType;

		cameras[i].T[0] = pi.gx;
		cameras[i].T[1] = pi.gy;
		cameras[i].T[2] = pi.height;
		cameras[i].yaw   = pi.yaw;
		cameras[i].pitch = -pi.pitch;
		cameras[i].roll  = -pi.roll;

		printf("%lf %lf %lf %lf %lf %lf \n", pi.gx, pi.gy, pi.height, pi.pitch, pi.roll, pi.yaw);

		cameras[i].bIsExplicit = true;

		eular2rot(cameras[i].R, cameras[i].pitch, cameras[i].roll, cameras[i].yaw);
	}

	//241: 1920,2066
    //242: 1480,2059
	vector<Point2DDouble> pts;
	Point3DDouble gps;
	double ferror;
	CTriangulateBase* pTriangulate = new CTriangulateCV();

	double x1 = 4809.84;// 5164; //4428;
	double y1 = 2282.86;// 2095; //2195;
	double x2 = 5232.16;// 5223; //4461;
	double y2 = 2328.75;// 2092; //2207;

	printf("%lf %lf %lf %lf \n", x1, y1, x2, y2);

	pts.resize(2);
	pts[0].p[0] = x1 - wd*0.5;
	pts[0].p[1] = ht*0.5 - y1;
	pts[1].p[0] = x2 - wd*0.5;
	pts[1].p[1] = ht*0.5 - y2;

	pTriangulate->Triangulate(pts, cameras, gps, true, ferror);

	printf("grd: %lf %lf %lf  error: %lf \n", gps.p[0], gps.p[1], gps.p[2], ferror);

	//calculate projections
	
	//GrdToPanoImageCenter(gps.p[0], gps.p[1], gps.p[2], radius, ix, iy);
	Point2DDouble ip1, ip2;
	GrdToImg(gps, ip1, cameras[0]);
	GrdToImg(gps, ip2, cameras[1]);

	printf("image point: %lf %lf  %lf %lf \n", ip1.p[0]+wd*0.5, ht*0.5-ip1.p[1],
		ip2.p[0]+wd*0.5, ht*0.5-ip2.p[1]);

	return 0;
}

int TestJYZDetect(int argc, char* argv[])
{
	char* cfgfile = "";
	char* weightfile = "";

	//network* net = load_network(cfgfile, weightfile, 0);

	return 0;
}


//
//int TestJYZDetect(int argc, char* argv[])
//{
//	//string filename = "C:\\Work\\Data\\segmentation\\qx97-crop.jpg";
//	//string filename = "C:\\Work\\Data\\segmentation\\qx3-crop.jpg";
//	//string filename = "C:\\Work\\Data\\segmentation\\qx59-crop.jpg";
//	//string filename = "C:\\Work\\Data\\segmentation\\nw000522-crop.jpg";
//	//string filename = "C:\\Work\\Data\\segmentation\\nw001449-crop.jpg";
//
//	/*
//	//for uav cars
//	char* imagepath = "C:\\Work\\Data\\deeplearning\\uav-car\\DL_XML";
//	//char* imagepath = "C:\\Work\\Data\\deeplearning\\ir\\images";
//	char* cfgfile    = "C:\\Work\\Programs\\DeepLearning\\yolo\\yolov2\\darknet-master\\cfg\\uav\\yolo-uav.cfg";
//	////char* labelfile  = "C:\\Work\\Programs\\DeepLearning\\darknet-master\\data\\uavobj.names";
//	char* weightfile = "C:\\Work\\Data\\deeplearning\\uav-car\\backup\\yolo-uav_1600.weights";
//	//char* outpath = "c:\\temp\\uav-car";
//	*/
//
//
//	/*
//	//for jyz
//	char* cfgfile    = "C:\\Work\\Programs\\DeepLearning\\darknet-master\\cfg\\jyz\\yolo-jyz.cfg";
//	char* weightfile = "C:\\Work\\Data\\deeplearning\\jyz\\backup\\yolo-jyz_final.weights";
//
//	//for jyz-loss
//	//char* cfgfile    = "C:\\Work\\Programs\\DeepLearning\\darknet-master\\cfg\\jyz\\yolo-jyz.cfg";
//	//char* weightfile = "C:\\Work\\Data\\deeplearning\\jyz-loss\\yolo-jyz-loss_1000.weights";
//	//image path
//	//char* imagepath = "C:\\Work\\Data\\deeplearning\\jyz-loss\\quexian";
//	//char* imagepath = "C:\\Work\\Data\\deeplearning\\jyz\\2018.03.14_test";
//	char* imagepath = "C:\\Work\\Data\\test";
//	*/
//
//	
//	//for infra-light car
//	char* cfgfile    = "C:\\Work\\Programs\\DeepLearning\\yolo\\yolov2\\darknet-master\\cfg\\jyz\\yolo-jyz.cfg";
//	char* weightfile = "C:\\Work\\Data\\deeplearning\\ir\\yolo-jyz_3500.weights";
//	char* imagepath  = "C:\\Work\\Data\\deeplearning\\ir\\images";
//	
//
//
//	//voc
//	/*char* cfgfile    = "C:\\Work\\Programs\\DeepLearning\\yolo\\yolov2\\darknet-master\\cfg\\yolo-voc.cfg";
//	char* weightfile = "C:\\Work\\Data\\deeplearning\\yolo\\yolo.weights";
//	char* imagepath  = "C:\\Work\\Data\\deeplearning\\ir\\images";*/
//
//
//	char* outpath    = "C:\\temp";
//
//
//	/*if (argc == 2){
//		filename = argv[1];
//	}
//	*/
//
//	char** filenames = NULL;
//	int n=0, nfile=0;
//	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 0);
//	printf("%d \n", nfile);
//	filenames = f2c(nfile, 256);
//	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 1);
//	if (nfile == 0)
//	{
//		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 0);
//		filenames = f2c(nfile, 256);
//		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 1);
//	}
//	printf("%d \n", nfile);
//
//	//namedWindow("image");
//	//startWindowThread();
//
//	//batch detection using YOLO2
//	//DetectJYZ(cfgfile, weightfile, labelfile, filenames, nfile, outpath, 0.07);
//
//	network* net = load_network(cfgfile, weightfile, 0);
//
//	box* pboxes = new box[MAX_BOX_NUM];
//	int  nbox;
//	for (int i = 0; i < nfile; i++)
//	{
//		printf("file: %s %d \n", filenames[i], i);
//
//		//char* errorfile = "C:\\Work\\Data\\deeplearning\\jyz-loss\\quexian\\qx138.jpg";
//		//strcpy(filenames[i], errorfile);
//
//		//waitKey(0);
//
//		//draw box
//		Mat img = imread(filenames[i]);
//
//		nbox = 0;
//		DetectObject(net, filenames[i], 0.05, pboxes, &nbox);
//		printf("box number: %d \n", nbox);
//
//		
//		if (nbox < 1)
//			continue;
//
//		Mat outimage = img;
//		for (int i = 0; i < nbox; i++){
//			box b = pboxes[i];
//			rectangle(outimage, Rect(b.x, b.y, b.w, b.h),CV_RGB(255, 0, 0), 11);
//		}
//		char outfile[256];
//		sprintf(outfile, "%s\\%d.jpg", outpath, i);
//		printf("%s \n", outfile);
//		imwrite(outfile, outimage);
//		
//
//		//imshow("image", outimage);
//		
//
//
//		//loss detection
//		if (0){
//			bool bIsLoss = false;
//			vector<vector<Rect>> lossRects;
//			lossRects.resize(nbox);
//			for (int k = 0; k < nbox; k++){
//				box b = pboxes[k];
//				//crop image
//				Mat cropimage = img(Rect(b.x, b.y, b.w, b.h));
//				//vector<Rect> lossRects;
//				DetectJYZLoss(cropimage, lossRects[k]);
//				if (lossRects[k].size()>0)
//					bIsLoss = true;
//				/*char file[256];
//				sprintf(file, "c:\\temp\\crop-%d-%d.jpg", i, k);
//				imwrite(file, cropimage);	*/
//			}
//
//			//draw rects
//			for (int i = 0; i < lossRects.size(); i++){
//				box b = pboxes[i];
//				printf("%d \n", lossRects[i].size());
//				for (int ki = 0; ki < lossRects[i].size(); ki++)
//				{
//					//printf("%d %d %d %d \n", int(b.x) + lossRects[i][ki].x, b.y + lossRects[i][ki].y,
//					//	lossRects[i][ki].width, lossRects[i][ki].height);
//					rectangle(img, Rect(b.x + lossRects[i][ki].x, b.y + lossRects[i][ki].y,
//						lossRects[i][ki].width, lossRects[i][ki].height),
//						CV_RGB(255, 0, 0), 21);
//				}
//			}
//
//			if (bIsLoss){
//				char* pdes = strrchr(filenames[i], '\\');
//				int   index = pdes - filenames[i];
//				char title[512];
//				strcpy(title, filenames[i] + index + 1);
//				char outfile[256];
//				sprintf(outfile, "%s\\%s", outpath, title);
//				printf("%s \n", outfile);
//				imwrite(outfile, img);
//
//				/*
//				//output rectangle
//				char boxfile[256];
//				sprintf(boxfile, "%s\\%s.txt", outpath, title);
//				printf("%s \n", boxfile);
//
//				FILE* fp = fopen(boxfile, "w");
//				fprintf(fp, "%s \n", "insulatord");
//				for (int i = 0; i < lossRects.size(); i++){
//				box b = pboxes[i];
//				//printf("%d \n", lossRects[i].size());
//				for (int ki = 0; ki < lossRects[i].size(); ki++)
//				{
//				fprintf(fp, "%d %d %d %d \n",
//				int(b.x) + lossRects[i][ki].x,
//				int(b.y) + lossRects[i][ki].y,
//				int(b.x) + lossRects[i][ki].x + lossRects[i][ki].width,
//				int(b.y) + lossRects[i][ki].y + lossRects[i][ki].height);
//				}
//				}
//				fclose(fp);
//				*/
//			}
//		}			
//	}
//	delete pboxes;
//
//	////detect the loss area using CV method
//	//Mat Img = imread(filename);
//	//vector<Rect> lossRects;
//	//DetectJYZLoss(Img, lossRects);
//	//printf("detect loss number: %d \n", lossRects.size());
//	//for (int i = 0; i < lossRects.size(); i++)
//	//{
//	//	rectangle(Img, lossRects[i], CV_RGB(255, 0, 0), 2);
//	//}
//	//imwrite("c:\\temp\\detecion.jpg", Img);
//
//	return 0;
//}
//
//
int TestJYZLoss()
{

	//char* file = "C:\\temp\\crop-19-1.jpg";
	char* file = "C:\\Work\\Data\\segmentation\\qx3-crop.jpg";

	Mat img = imread(file);

	vector<Rect> lossRects;
	DetectJYZLoss(img, lossRects);

	for (int i = 0; i < lossRects.size(); i++){
		rectangle(img, Rect(lossRects[i].x, lossRects[i].y, lossRects[i].width, lossRects[i].height),
			CV_RGB(255, 0, 0), 5);
	}

	imwrite("c:\\temp\\detection-ori.jpg", img);


	return 0;
}

int SimulateDataBA()
{
	//read the feature point coordinates
	int ht = 4096;
	int wd = 8192;
	int nCamera = 4;
	char *datapath = "C:\\Work\\Programs\\python\\ba";
	vector<ImgFeature> imgFeatures;
	imgFeatures.resize(nCamera);
	for (int i = 0; i < nCamera; i++)
	{
		char filename[256];
		sprintf(filename, "%s\\cam%d.txt", datapath, i);
		printf("%s \n", filename);
		FILE* fp = fopen(filename, "r");

		int rows = GetFileRows(filename);
		stPtFeature pt;
		for (int k = 0; k < rows; k++){
			double x, y;
			fscanf(fp, "%lf,%lf", &x, &y);
			printf("%lf %lf \n", x, y);
			pt.cx = x-wd*0.5;
			pt.cy = ht*0.5-y;
			pt.extra = k;    //track index
			imgFeatures[i].AddFeatPt(pt);
		}
		fclose(fp);
	}

	//read the POS data
	char* posfile = "C:\\Work\\Data\\Panorama\\beijing\\image_all\\0226.cam";
	char* posfile_right = "C:\\Work\\Data\\Panorama\\beijing\\noloss.cam";
	char* imagepath = "C:\\Work\\Data\\Panorama\\beijing\\image_all";
	
	CameraType camType = PanoramCam;
	vector<CameraPara> cameras;
	double agx = 0, agy = 0, agz = 0;
	int nstart = 0;
	int nrange = nCamera;
	ReadPosData(posfile, imagepath, nstart, nrange, camType, cameras, agx, agy, agz);

	double agx1 = 0, agy1 = 0, agz1 = 0;
	vector<CameraPara> camerasRight;
	ReadPosData(posfile_right, imagepath, nstart, nrange, camType, camerasRight, agx1, agy1, agz1);



	//generate the tracks
	int nTrack = 24;
	vector<TrackInfo> tracks;
	tracks.resize(nTrack);
	for (int i = 0; i < nTrack; i++){
		ImageKey key;
		for (int k = 0; k < nCamera; k++){
			key.first  = k; //camera index
			key.second = i; //point index
			tracks[i].views.push_back(key);
		}
	}

	//calculate the control points
	vector<TrackInfo> tracksRight = tracks; 
	CaculateTrackSeqGrd(imgFeatures, tracksRight, camerasRight, true);
	printf("control points: \n");
	for (int i = 0; i<tracksRight.size(); i++)
	{
		printf("%d  %lf %lf %lf   %lf \n",
			i,
			tracksRight[i].grd.p[0] + agx1, 
			tracksRight[i].grd.p[1] + agy1,
			tracksRight[i].grd.p[2] + agz1,
			tracksRight[i].GetProjectionError());
	}

	//////calculate the 3D coordinates for optimization
	//CaculateTrackSeqGrd(imgFeatures, tracks, cameras, true);
	//printf("Track grd for optimization: \n");
	//for (int i = 0; i<tracks.size(); i++)
	//{
	//	printf("%lf %lf %lf  %lf \n",
	//		tracks[i].grd.p[0] + agx, 
	//		tracks[i].grd.p[1] + agy,
	//		tracks[i].grd.p[2] + agz,
	//		tracks[i].GetProjectionError());
	//}

	//set control points
	vector<int> ctrlTrackIndex;
	ctrlTrackIndex.push_back(5);
	ctrlTrackIndex.push_back(6);
	ctrlTrackIndex.push_back(14);
	ctrlTrackIndex.push_back(15);
	//ctrlTrackIndex.push_back(13);
	//ctrlTrackIndex.push_back(16);
	//ctrlTrackIndex.push_back(3);


	for (int i = 0; i < ctrlTrackIndex.size(); i++)
	{
		int id = ctrlTrackIndex[i];
		tracks[id].SetAsCtrlPt();
		tracks[id].grd = tracksRight[id].grd;
	}


	///////////////////////  bundle adjustment based on POS ////////////////////////
	CPanoBAWithPos* pBA = new CPanoBAWithPos();
	pBA->SetCenter(agx, agy, agz);
	pBA->BundleAdjust(cameras.size(), cameras, imgFeatures, tracks);
	delete pBA;

	//output the trajectory
	FILE* fp = fopen("c:\\temp\\trajectory.txt", "w");
	for (int i = 0; i < cameras.size(); i++)
	{
		fprintf(fp, "%s,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf\n",
			cameras[i].title,
			cameras[i].gx, cameras[i].gy, cameras[i].gz,
			cameras[i].yaw, cameras[i].pitch, cameras[i].roll,
			(cameras[i].T[0] + agx), (cameras[i].T[1] + agy), (cameras[i].T[2] + agz),
			cameras[i].az, cameras[i].ax, cameras[i].ay);
	}
	fclose(fp);

	//output the projection error
	printf("Projection Error .... \n");
	for (int i = 0; i < tracks.size(); i++)
	{
		Point3DDouble grd = tracks[i].grd;

		int nview = tracks[i].GetViews();
		for (int k = 0; k < nview; k++)
		{
			ImageKey ik = tracks[i].GetImageKey(k);
			int camid = ik.first;
			int ptid  = ik.second;
			
			Point2DDouble ip;
			GrdToImg(grd, ip, cameras[camid]);
			double dx = ip.p[0] - imgFeatures[camid].featPts[ptid].cx;
			double dy = ip.p[1] - imgFeatures[camid].featPts[ptid].cy;
			printf("%.2lf %.2lf   ", fabs(dx), fabs(dy));
		}
		printf("\n");
	}

	//distance between the camera centers
	/*printf("Camera Distance .... \n");
	for (int i = 0; i < cameras.size() - 1; i++){
		double gx1 = cameras[i].T[0];
		double gy1 = cameras[i].T[1];
		double gz1 = cameras[i].T[2];

		double gx2 = cameras[i+1].T[0];
		double gy2 = cameras[i+1].T[1];
		double gz2 = cameras[i+1].T[2];

		double distOpt = sqrt((gx1 - gx2)*(gx1 - gx2) + 
			(gy1 - gy2)*(gy1 - gy2) + (gz1 - gz2)*(gz1 - gz2));

		gx1 = cameras[i].gx;
		gy1 = cameras[i].gy;
		gz1 = cameras[i].gz;

		gx2 = cameras[i + 1].gx;
		gy2 = cameras[i + 1].gy;
		gz2 = cameras[i + 1].gz;

		double distOri = sqrt((gx1 - gx2)*(gx1 - gx2) +
			(gy1 - gy2)*(gy1 - gy2) + (gz1 - gz2)*(gz1 - gz2));

		printf("%.2lf %.2lf  \n", distOri, distOpt);
	}*/

	printf("trajectory error ... \n");
	for (int i = 0; i < cameras.size(); i++){
		double ox = cameras[i].T[0] + agx;
		double oy = cameras[i].T[1] + agy;
		double oz = cameras[i].T[2] + agz;

		double px = cameras[i].gx;
		double py = cameras[i].gy;
		double pz = cameras[i].gz;
	
		printf("%.4lf %.4lf %.4lf  \n", ox-px, oy-py, oz-pz);
	}

	return 0;
}

int TestCeres(){

	test_robust_fitting();
	
	return 0;
}


int GenerateCtrlPts()
{
	//reading image points



	return 0;
}

int BAWithControlPt()
{
	int ht = 4096;
	int wd = 8192;


	//read the POS data
	char* posfile = "C:\\Work\\Data\\Panorama\\beijing\\image_all\\0226.cam";
	char* posfile_right = "C:\\Work\\Data\\Panorama\\beijing\\noloss.cam";
	char* imagepath = "C:\\Work\\Data\\Panorama\\beijing\\image_all";


	//reading image filenames
	char** filenames = NULL;
	int n = 0;
	int nfile = 0;
	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 0);
	printf("%d \n", nfile);
	filenames = f2c(nfile, 256);
	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 1);
	if (nfile == 0)
	{
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 0);
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 1);
	}
	printf("image number: %d \n", nfile);

	int nrange = 4;
	int start = 0;
	int endindex = min(start + nrange, nfile - 1);

	char** selectedFiles = f2c(endindex - start, 256);
	int id = 0;
	for (int i = start; i <endindex; i++)
	{
		strcpy(selectedFiles[id], filenames[i]);
		id++;
	}
	
	int nimage = endindex - start;


	//reading cameras
	CameraType camType = PanoramCam;
	vector<CameraPara> cameras;
	double agx = 0, agy = 0, agz = 0;
	ReadPosData(posfile, imagepath, start, nrange, camType, cameras, agx, agy, agz);

	double agx1 = 0, agy1 = 0, agz1 = 0;
	vector<CameraPara> camerasRight;
	ReadPosData(posfile_right, imagepath, start, nrange, camType, camerasRight, agx1, agy1, agz1);

	//reading the projections of control point
	char* ctrlPtfile = "C:\\Work\\Programs\\python\\ba\\ctrlpt.txt";
	int nrows = GetFileRows(ctrlPtfile);
	vector<stCtrlPt> ctrlPts;
	FILE* fp = fopen(ctrlPtfile, "r");
	for (int i = 0; i < nrows; i++)
	{
		int id;
		double number;
		int camid;
		double cx, cy;

		stCtrlPt cpt;
		fscanf(fp, "%d %lf", &id, &number);
		for (int k = 0; k < number; k++)
		{
			fscanf(fp, "%d %lf %lf ", &camid, &cx, &cy);
			Point2DDouble ip;
			ip.p[0] = cx - wd*0.5;
			ip.p[1] = ht*0.5 - cy;
			cpt.camIds.push_back(camid);
			cpt.projPts.push_back(ip);
		}
		ctrlPts.push_back(cpt);
	}

	//calculate the 3D coordinates of control points
	CalculateCtrlPt(ctrlPts, camerasRight, true);
	for (int i = 0; i<ctrlPts.size(); i++)
	{
		printf("%d  %lf %lf %lf   %lf \n",
			i,
			ctrlPts[i].grd.p[0] + agx1,
			ctrlPts[i].grd.p[1] + agy1,
			ctrlPts[i].grd.p[2] + agz1,
			ctrlPts[i].derror);
	}

	//1.generating the image feature points
	vector<ImgFeature> imgFeatures;
	double dProgress = 0;
	DetectFileFeaturePts(selectedFiles, nimage, imgFeatures, 900, dProgress);

	//2. matching images 
	vector<PairMatchRes> matchRes;
	MatchImageFiles(imgFeatures, matchRes, camType, 3);

	//3.generating the tracks
	vector<TrackInfo> tracks;
	CGenerateTracksBase* pGenerateTrack = new CFastGenerateTrack();
	pGenerateTrack->GenerateTracks(imgFeatures, matchRes, tracks);

	//4.ba
	CPanoGlobalBA gba;
	gba.Init(cameras, imgFeatures, tracks);
	gba.AddFreePtBlock();
	gba.Run();
	gba.Output(cameras, tracks);

	printf("trajectory error ... \n");
	for (int i = 0; i < cameras.size(); i++){
		double ox = cameras[i].T[0] + agx;
		double oy = cameras[i].T[1] + agy;
		double oz = cameras[i].T[2] + agz;

		double px = camerasRight[i].gx;
		double py = camerasRight[i].gy;
		double pz = camerasRight[i].gz;

		printf("%.4lf %.4lf %.4lf  \n", ox - px, oy - py, oz - pz);
	}

	return 0;
}

int TestPosConvertion()
{

	//char* posfile   = "C:\\Work\\Data\\Panorama\\zibo\\0403\\2\\0331_am_new.cam";
	//char* imagepath = "C:\\Work\\Data\\Panorama\\zibo\\0403\\2";
	//double x1 = 5827.65; 
	//double y1 = 2083.2;  
	//double x2 = 6309.81; 
	//double y2 = 2095.93; 
	//double x3 = 6760.05;
	//double y3 = 2104.14;

	char* inputfile = "C:\\Work\\Data\\Panorama\\zibo\\0408-rtk\\4\\input.txt";
	//char* inputfile = "C:\\Work\\Data\\Panorama\\zibo\\0408-rtk\\10\\input.txt";

	char   posfile[256];
	char   imagepath[256];
	double x[3];
	double y[3];
	int    numpt=0;

	FILE* fp = fopen(inputfile, "r");
	fscanf(fp, "%s", imagepath);
	fscanf(fp, "%s", posfile);
	fscanf(fp, "%d", &numpt);
	for (int i = 0; i < numpt; i++)
	{
		fscanf(fp, "%lf %lf ", &(x[i]), &(y[i]) );
	}
	fclose(fp);

	int startindex = 0;
	int nimage = numpt;

	int ht = 4096;
	int wd = 8192;
	CameraType camType = PanoramCam;

	vector<CameraPara> cameras;
	cameras.resize(nimage);
	double agx = 0, agy = 0, agz = 0;
	ReadPosData(posfile, imagepath, startindex, nimage, PanoramCam, cameras, agx, agy, agz);


	vector<Point2DDouble> pts;
	Point3DDouble gps;
	double ferror;
	CTriangulateBase* pTriangulate = new CTriangulateCV();

	//printf("input image point: %lf %lf %lf %lf \n", x1, y1, x2, y2);

	pts.resize(numpt);
	for (int i = 0; i < numpt; i++)
	{
		pts[i].p[0] = x[i] - wd*0.5;
		pts[i].p[1] = ht*0.5 - y[i];
	}

	//pts[0].p[0] = x1 - wd*0.5;
	//pts[0].p[1] = ht*0.5 - y1;
	//pts[1].p[0] = x2 - wd*0.5;
	//pts[1].p[1] = ht*0.5 - y2;
	//pts[2].p[0] = x3 - wd*0.5;
	//pts[2].p[1] = ht*0.5 - y3;

	pTriangulate->Triangulate(pts, cameras, gps, true, ferror);
	printf("grd: %lf %lf %lf  error: %lf \n", gps.p[0]+agx, gps.p[1]+agy, gps.p[2]+agz, ferror);

	fp = fopen("c:\\temp\\testgrd.txt", "w");
	fprintf(fp, "%lf %lf %lf \n", gps.p[0] + agx, gps.p[1] + agy, gps.p[2] + agz);
	fclose(fp);

	//GrdToPanoImageCenter(gps.p[0], gps.p[1], gps.p[2], radius, ix, iy);
	Point2DDouble ip1, ip2;
	GrdToImg(gps, ip1, cameras[0]);
	GrdToImg(gps, ip2, cameras[1]);
	printf("reprojected image point: %lf %lf  %lf %lf \n", ip1.p[0] + wd*0.5, ht*0.5 - ip1.p[1],
		ip2.p[0] + wd*0.5, ht*0.5 - ip2.p[1]);

	/*printf("input lon, lat :   \n");
	double lat, lon;
	UTMtoLL(23, gps.p[1] + agy, gps.p[0] + agx, 50, lat, lon);
	printf("%lf %lf \n", lat, lon);

	int nd = int(lat);
	int ns = int((lat - nd) * 60);
	double nm = ((lat - nd) * 60 - ns) * 60;
	
	printf("%d %d %.8lf  %lf \n", nd, ns, nm);

	nd = int(lon);
	ns = int((lon - nd) * 60);
	nm = ((lon - nd) * 60 - ns) * 60;

	printf("%d %d %.8lf  %lf \n", nd, ns, nm);*/

	//double gx, gy;
	//LLtoUTM(23, 36.83668533, 118.02269237, gy, gx, 50);
	////LLtoUTM(23, lat, lon, gy, gx, 50);
	//printf("%lf %lf \n", gx, gy);


	return 0;
}



int BAFreeNet()
{
	time_t  start = time(NULL);
	//double t = (double)getTickCount();

	//set the type of camera
	CameraType camType = PerspectiveCam; //;PanoramCam;
	//CameraType camType = PanoramCam;
	printf("SFM integration .... \n");

	//panorama data
	//char* imagepath = "C:\\Work\\Data\\Panorama\\beijing\\image_all";
	//char imagepath[256] = "C:\\Work\\Data\\Panorama\\20180226_test_xian";
	//char imagepath[256] = "C:\\Work\\Data\\panorama\\2016.10.13-yizhuang\\L10_1013\\indoor\\jpg";
	//char imagepath[256] = "C:\\Work\\Data\\panorama\\cheku20161011\\garageoutput";
	//char imagepath[256] = "C:\\Work\\Data\\panorama\\cheku2016.10.26\\Output";
	char imagepath[256] = "C:\\Work\\Data\\UAV\\binchuan";
	//char imagepath[256] = "C:\\Work\\Programs\\PgStation\\Bundler_V0.4\\examples\\kermit";
	char* outpath = "C:\\Work\\Programs\\PgStation\\CMVS\\bin\\pmvs";

	char** filenames = NULL;
	int n = 0;
	int nfile = 0;

	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 0);
	printf("%d \n", nfile);

	filenames = f2c(nfile, 256);
	GetDirFileName(filenames, imagepath, &n, &nfile, "JPG", 1);

	if (nfile == 0)
	{
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 0);
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, imagepath, &n, &nfile, "jpg", 1);
	}
	printf("image number: %d \n", nfile);

	if (nfile<2)
	{
		printf("images are less than 2 ! \n");
		return -1;
	}

	//reading pos from the JPEG file header
	vector<CameraPara> cameras;
	cameras.resize(nfile);
	ReadingPosFromJpegImages(filenames, nfile, cameras);


	//1. generating the image feature points
	vector<ImgFeature> imgFeatures;
	double dProgress = 0;
	DetectFileFeaturePts(filenames, nfile, imgFeatures, 900, dProgress);

	//2. matching images 
	vector<PairMatchRes> matchRes;
	MatchImageFiles(imgFeatures, matchRes, camType, 3);// imgFeatures.size());

	//3.generating the tracks
	vector<TrackInfo> tracks;
	CGenerateTracksBase* pGenerateTrack = new CFastGenerateTrack();
	pGenerateTrack->GenerateTracks(imgFeatures, matchRes, tracks);

	//output the tracks
	FILE* fp = fopen("c:\\temp\\mytrack.txt", "w");
	fprintf(fp, "%d \n", tracks.size());
	for (int i = 0; i<tracks.size(); i++)
	{
		fprintf(fp, "%d  ", tracks[i].views.size());
		for (int j = 0; j<tracks[i].views.size(); j++)
		{
			fprintf(fp, " %d %d ", tracks[i].views[j].first, tracks[i].views[j].second);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);


	//4. bundle adjustment
	int numImage = imgFeatures.size();

	double focalLen = (imgFeatures[0].ht + imgFeatures[0].wd) * 0.5;
	for (int i = 0; i<numImage; i++)
	{
		cameras[i].focalLen = focalLen; //initialize the focal length 
		memset(cameras[i].R, 0, sizeof(double)* 9);
		cameras[i].R[0] = 1;
		cameras[i].R[4] = 1;
		cameras[i].R[8] = 1;
		cameras[i].rows = imgFeatures[0].ht;
		cameras[i].cols = imgFeatures[0].wd;
		cameras[i].camtype = camType;
	}

	CBABase* pBA = NULL;

	if (camType == PerspectiveCam)
	{
		pBA = new CCeresBA();
	}
	else if (camType == PanoramCam)
	{
		pBA = new CPanoBA();
	}

	pBA->BundleAdjust(cameras.size(), cameras, imgFeatures, matchRes, tracks, false);

	delete pBA;


	//t = ((double)getTickCount() - t) / getTickFrequency();
	time_t  end = time(NULL);

	printf("\n running time: %d \n", end - start);


	//absolute orientation
	stAbsPOS absparas;
	AbsOriOrthogonal(absparas, cameras);

	//generate dem


	//orthoimage




	return 0;
}


int TestG2o()
{
	return 0;
}


//new api for free bundle adjustment
int MosaicTest()
{
	CSFMSystem sfmSystem;

	//string imagepath = "C:\\Work\\Data\\UAV\\bingchuan";
	//string imagepath   = "C:\\Work\\Data\\uav\\heilongjiang1";
	//string projectpath = "C:\\Work\\Data\\uav\\project1";

	//string imagepath   = "F:\\Data\\uav\\shanxi\\temp1";
	//string projectpath = "F:\\Data\\uav\\shanxi\\temp1\\work";
    
	//string imagepath   = "C:\\Work\\Data\\uav\\heilongjiang\\images";
	//string projectpath = "C:\\Work\\Data\\uav\\heilongjiang\\work";

    //string imagepath   = "C:\\Work\\Data\\UAV\\zhengshe";
    //string projectpath = "C:\\Work\\Data\\UAV\\zhengshe\\work";

    string imagepath = "C:\\Work\\Data\\UAV\\zhong";
    string projectpath = "C:\\Work\\Data\\UAV\\zhong\\work";
    
    //string imagepath   = "C:\\Work\\Data\\UAV\\bingchuan1";
    //string projectpath = "C:\\Work\\Data\\UAV\\bingchuan1\\work";


	string mosaicfile = projectpath + "\\" + "mosaic.tif";

	CameraType camType = PerspectiveCam;

	
	time_t start, end;
	double cost;
	time(&start);

	//sparse points reconstruction
	sfmSystem.SetProjectPath(projectpath);
	sfmSystem.LoadImages(imagepath);
	sfmSystem.DetectFeatures(640);
	sfmSystem.ImageMatch(camType);
	sfmSystem.GenerateTracks();

	if (sfmSystem.BA(camType) < 0)
	{
		printf("BA failed! \n");
		return -1;;
	}
	//generate orthoimage and mosaic
	if (sfmSystem.AbsOrientation() >= 0){
		sfmSystem.DEMInterpolation();
		sfmSystem.GenerateOrthoImages();
		sfmSystem.Fusion(0.1, mosaicfile);
	}
	else{
		printf("absolute orientation failed! \n");
	}

	time(&end);
	cost = difftime(end, start);
	
	printf("time: %lf (minutes) \n", cost/60.0);

	return 0;
}

int ABSOriTest()
{
	string imagepath = "C:\\Work\\Data\\uav\\bingchuan\\images";
	string projectpath = "C:\\Work\\Data\\uav\\bingchuan\\work";

	string baout = projectpath + "\\" + "ba.out";
	string absout = projectpath + "\\" + "aor.out";

	CSFMSystem sfmSystem;
	sfmSystem.SetProjectPath(projectpath);
	sfmSystem.LoadImages(imagepath);


	sfmSystem.ReadCamsAndGrdpts(baout);
	sfmSystem.AbsOrientation();
	sfmSystem.DEMInterpolation();
	//sfmSystem.ReadCamsAndGrdpts(absout);
	//sfmSystem.GenerateOrthoImages();

	return 0;
}

void testRotation()
{
    double R[9] = { 0.421251525498582, -0.598691389531382 ,-0.681260429179694,
        0.862850635906912 ,  0.495923109239094 ,  0.097719239854121,
        0.279349122748140, -0.628990373395024  , 0.725489612466114 };

    //double iR[9];

    invers_matrix(R, 3);

    //rotation to rodrigues angles
    double aa[3];
    rot2aa(R, aa);

    printf("%lf %lf %lf \n", aa[0], aa[1], aa[2]);

}

Vec3b RandomColor(int value)
{
    value = value % 255;  //生成0~255的随机数
    RNG rng;
    int aa = rng.uniform(0, value);
    int bb = rng.uniform(0, value);
    int cc = rng.uniform(0, value);
    return Vec3b(aa, bb, cc);
}

void LaneDetect()
{

    char* filename = "C:\\Work\\Data\\pointcloud\\sanding_process\\sanding_process\\161103_065409_0_original.png";

    Mat ori = imread(filename);
    //medianBlur(ori, ori, 3);
        

    //interpolate the black holes
    Mat dilateImage;    
    int erosion_size = 2;
    Mat element = getStructuringElement(MORPH_ELLIPSE,
        Size(2 * erosion_size + 1, 2 * erosion_size + 1),
        Point(erosion_size, erosion_size));

    //erode(ori, ori, element);
    dilate(ori, dilateImage, element);
    imwrite("c:\\temp\\lane_dilate.jpg", dilateImage);
    
    //smooth filter
    Mat smoothImage;
    medianBlur(dilateImage, smoothImage, 3);
    //GaussianBlur(smoothImage, smoothImage, Size(5, 5), 1);
    imwrite("c:\\temp\\lane_filter.jpg", smoothImage);

    //rgb to gray
    Mat grayImage;
    cvtColor(smoothImage, grayImage, CV_RGB2GRAY);

    DetectLane(grayImage);


    /*
    //segmentation
    Mat segImage; 
    adaptiveThreshold(grayImage, segImage, 255,ADAPTIVE_THRESH_GAUSSIAN_C, 
        THRESH_BINARY, 5, 0);
    //threshold(grayImage, segImage, 0, 255, CV_THRESH_BINARY|CV_THRESH_OTSU);
    imwrite("c:\\temp\\lane_seg.jpg", segImage);


    //edge
    Canny(grayImage, grayImage, 80, 150);
    //imshow("Canny Image", grayImage);
    imwrite("c:\\temp\\lane_edge.jpg", grayImage);
    */

    /*
    //find contours
    vector<vector<Point>> contours;
    vector<Vec4i> hierarchy;
    findContours(grayImage, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point());
    Mat imageContours = Mat::zeros(grayImage.size(), CV_8UC1);
    Mat marks(grayImage.size(), CV_32S);   //Opencv分水岭第二个矩阵参数
    marks = Scalar::all(0);
    int index = 0;
    int compCount = 0;
    for (; index >= 0; index = hierarchy[index][0], compCount++)
    {
        drawContours(marks, contours, index, Scalar::all(compCount + 1), 1, 8, hierarchy);
        drawContours(imageContours, contours, index, Scalar(255), 1, 8, hierarchy);
    }
    
    Mat marksShows;
    convertScaleAbs(marks, marksShows);
    //imshow("marksShow", marksShows);
    //imshow("contours", imageContours);
    watershed(smoothImage, marks);

    Mat afterWatershed;
    convertScaleAbs(marks, afterWatershed);
    //imshow("After Watershed", afterWatershed);

    //对每一个区域进行颜色填充	
    Mat PerspectiveImage=Mat::zeros(smoothImage.size(),CV_8UC3);	
    for(int i=0;i<marks.rows;i++)	
    {		
        for(int j=0;j<marks.cols;j++)
        {
            int index=marks.at<int>(i,j);			
            if(marks.at<int>(i,j)==-1)			
            {				
                PerspectiveImage.at<Vec3b>(i,j)=Vec3b(255,255,255);			
            }			 			
            else			
            {				
                PerspectiveImage.at<Vec3b>(i,j) = RandomColor(index);			
            }		
        }	
    }	
    //imshow("After ColorFill",PerspectiveImage);

    //分割并填充颜色的结果跟原始图像融合
    Mat wshed;
    addWeighted(smoothImage, 0.4, PerspectiveImage, 0.6, 0, wshed);
    //imshow("AddWeighted Image", wshed);

    waitKey();
    */

}

//image: gray 
void fillholes(Mat& image)
{
    Mat fillimage = image.clone();
    int rows = image.rows;
    int cols = image.cols;
    for (int j = 1; j < rows-1; j++)
    {
        uchar* p  = image.ptr<uchar>(j);
        uchar* dp = fillimage.ptr<uchar>(j);
        for (int i = 1; i < cols-1; i++)
        {
            int value = p[i];
            if (value > 0)
                continue;

            double sum = 0;
            int    nsum = 0;
            for (int m = -1; m <= 1; m++)
            {
                uchar* pl = image.ptr<uchar>(j+m);
                for (int n = -1; n <= 1; n++)
                {
                    int g = pl[i + n];
                    if (g > 0)
                    {
                        sum += g;
                        nsum++;
                    }
                }
            }

            if (nsum > 0)
            {
                double ave = sum / double(nsum);
                dp[i] = int(ave+0.5);
            }
        }
    }

    image.release();
    image = fillimage.clone();
}


void LaneDetect1()
{
    char* filename = "C:\\Work\\Data\\pointcloud\\sanding_process\\sanding_process\\161103_065409_0_original.png";
    Mat image = imread(filename);
    Mat gray;
    cvtColor(image, gray, CV_RGB2GRAY);
    imwrite("c:\\temp\\lane_ori.jpg", gray);
    
    //1. fill holes
    //int erosion_size = 1;
    //Mat element = getStructuringElement(MORPH_RECT,
    //    Size(2 * erosion_size + 1, 2 * erosion_size + 1),
    //    Point(erosion_size, erosion_size));
    ////erode(ori, ori, element);
    //dilate(gray, gray, element);
    fillholes(gray);
    fillholes(gray);
    //equalizeHist(gray, gray);
    imwrite("c:\\temp\\lane_fill.jpg", gray);
    int rows = gray.rows;
    int cols = gray.cols;

    /*
    float* pbuffer = new float[rows*cols];
    for (int j = 0; j < rows; j++)
    {
        uchar* p = gray.ptr<uchar>(j);
        for (int i = 0; i < cols; i++)
        {
            pbuffer[j*cols + i] = p[i];
        }
    }
    AodInterpolation(pbuffer, rows, cols);
    for (int j = 0; j < rows; j++)
    {
        uchar* p = gray.ptr<uchar>(j);
        for (int i = 0; i < cols; i++)
        {
            p[i] = (uchar)(pbuffer[j*cols + i]);
        }
    }
    delete[] pbuffer;
    imwrite("c:\\temp\\lane_interpolate.jpg", gray);
    */
    
    //canny edge
   /* Mat edgeimage = gray.clone();
    medianBlur(edgeimage, edgeimage, 3);
    Canny(edgeimage, edgeimage, 8, 16);
    imwrite("c:\\temp\\lane_edge.jpg", edgeimage);*/
    
    //2. binary image
    int lanewid = 12; //the width of lane 
    printf("%d %d \n", rows, cols);    
    //Mat binaryImage(rows, cols, CV_8UC1, 0);
    Mat binaryImage = gray.clone();
    //imwrite("c:\\temp\\lane_fore.jpg", binaryImage);
    for (int j = lanewid; j < rows-lanewid; j++)
    {
        printf(".");
        uchar *p = gray.ptr<uchar>(j);
        uchar *dp = binaryImage.ptr<uchar>(j);
        for (int i = lanewid; i < cols-lanewid; i++)
        {
            dp[i] = 0;

            if (p[i] == 0)
            {                
                continue;
            }

            double dx_left  = p[i] - p[i - lanewid];
            double dx_right = p[i] - p[i + lanewid];

            if (dx_left > 4 && dx_right > 4)
            {
                dp[i] = 255;
            }
        }
    }
    imwrite("c:\\temp\\lane_binary.jpg", binaryImage);


    //remove the noise
   /* int erosion_size = 1;
    Mat erose_element = getStructuringElement(MORPH_RECT,
        Size(2 * erosion_size + 1, 2 * erosion_size + 1),
        Point(erosion_size, erosion_size));
    erode(binaryImage, binaryImage, erose_element);
    
  */
    medianBlur(binaryImage, binaryImage, 3);
    int dilate_size = 2;
    Mat dilate_element = getStructuringElement(MORPH_ELLIPSE,
        Size(2 * dilate_size + 1, 2 * dilate_size + 1),
        Point(dilate_size, dilate_size));
    dilate(binaryImage, binaryImage, dilate_element);
    imwrite("c:\\temp\\lane_fore.jpg", binaryImage);

    //3. contours
    vector<vector<Point>> contours;
    std::vector<Vec4i> hierarchy;
    findContours(binaryImage, contours, hierarchy,
        CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
    printf("contour number: %d \n", contours.size());
    
    //contour approxy
    vector<vector<Point>> poly(contours.size());
    for (int i = 0; i < contours.size(); i++)
    {
        approxPolyDP(Mat(contours[i]), poly[i], 2, true);
    }


    Mat drawing = Mat::zeros(gray.size(), CV_8UC3);
    ////imwrite("c:\\temp\\lane_contours.jpg", drawing);

    RNG rng(12345);
    for (int i = 0; i < contours.size(); i++)
    {
        //if (contourArea(contours[i]) > 8 && arcLength(contours[i], false) > 16)
        //if (contourArea(contours[i]) > 8 && arcLength(contours[i], false) > 16)
        {
            printf("%d ", i);

            Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            //drawContours(drawing, contours, i, color, 1, 8, hierarchy);
            drawContours(drawing, poly, i, color, 2, 8, vector<Vec4i>(), 0, Point());
        }        
    }

    //boundingRect

    /*int index = 0;
    int compCount = 0;
    for (; index >= 0; index = hierarchy[index][0], compCount++)
    {
        Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
        drawContours(drawing, contours, index, Scalar(255), 1, 8, hierarchy);
    }*/
    imwrite("c:\\temp\\lane_contours.jpg", drawing);
}

void DllLaneDetecct()
{
    //string filename = "C:\\Work\\Data\\pointcloud\\sanding_process\\sanding_process\\161103_065409_0_original.png";
    string filename = "C:\\Work\\Data\\pointcloud\\sanding_process\\sanding_process\\161103_065409_1_original.png";
    //char* filename = "C:\\Work\\Data\\lane\\test.jpg";
    //Mat image = imread(filename);
    
    CVDetectLane(filename);
}

void testfusion()
{
    //char* imagepath = "C:\\Work\\Data\\UAV\\test1\\work";
    char* imagepath = "C:\\Work\\Data\\UAV\\zhong\\work";

    char** filenames = NULL;
    int n = 0, nfile = 0;

    GetDirFileName(filenames, imagepath, &n, &nfile, "jpeg", 0);
    filenames = f2c(nfile, 256);
    GetDirFileName(filenames, imagepath, &n, &nfile, "jpeg", 1);
    
    vector<stGeoInfo> mValidDOMGeoData;
    for (int i = 0; i < nfile; i++)
    {
        //printf("%s \n", filenames[i]);
        char geofile[256];
        strcpy(geofile, filenames[i]);
        char* p = strrchr(geofile, '.');
        strcpy(p + 1, "geo");
        printf("%s \n", geofile);

        //read the geoinfo
        stGeoInfo geoinfo;
        FILE* fp = fopen(geofile, "r");
        fscanf(fp, "%d %d %d  %lf %lf %lf %lf ", &(geoinfo.zoneNumber),
            &(geoinfo.ht), &(geoinfo.wd),
            &(geoinfo.dx), &(geoinfo.dy),
            &(geoinfo.left), &(geoinfo.top)
        );
        fclose(fp);

        mValidDOMGeoData.push_back(geoinfo);
    }

    //mosaic and fusion	
    double minx, maxx, miny, maxy;
    double maxrx, maxry;
    double rx, ry;
    minx = 5000000000;	maxx = -5000000000;
    miny = 5000000000;	maxy = -5000000000;
    maxrx = 0;
    maxry = 0;
    int zoneNumber = mValidDOMGeoData[0].zoneNumber;
    for (int i = 0; i<mValidDOMGeoData.size(); i++)
    {
        //if (zoneNumber != geoArray[i].zoneNumber)
        //	continue;		
        //resolution
        rx = mValidDOMGeoData[i].dx;
        ry = fabs(mValidDOMGeoData[i].dy);
        maxrx = max(rx, maxrx);
        maxry = max(ry, maxry);
        //position
        minx = min(minx, mValidDOMGeoData[i].left);
        maxx = max(maxx, mValidDOMGeoData[i].left + mValidDOMGeoData[i].wd*rx);
        miny = min(miny, mValidDOMGeoData[i].top - mValidDOMGeoData[i].ht*ry);
        maxy = max(maxy, mValidDOMGeoData[i].top);
    }
    printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

    double outResolution = 0.1;
    double resolution = outResolution;
    int oht = (maxy - miny) / resolution;
    int owd = (maxx - minx) / resolution;
    printf("%d %d \n", oht, owd);
    Mat oimage(oht, owd, CV_16SC3, Scalar(0,0,0));
    
    std::vector<cv::Mat> out_pyr_laplace;
    int nlevel = 4;
    cv::detail::createLaplacePyr(oimage, nlevel, out_pyr_laplace);
    printf("size of pyramid: %d \n", out_pyr_laplace.size());

    std::vector<cv::Mat> out_pyr_weights(nlevel + 1);
    out_pyr_weights[0].create(oht, owd, CV_32FC1);
    out_pyr_weights[0] = Scalar(0);
    for (int i = 0; i < nlevel; ++i)
        cv::pyrDown(out_pyr_weights[i], out_pyr_weights[i + 1]);

    for (int fi = 0; fi < nfile; fi++)
    {
        //convert the image type
        Mat image = imread(filenames[fi]);
        Mat srcimage;
        image.convertTo(srcimage, CV_16SC3);

        
        //create laplacian images
        std::vector<cv::Mat> pyr_laplace;
        cv::detail::createLaplacePyr(srcimage, nlevel, pyr_laplace);
        printf("size of pyramid: %d \n", pyr_laplace.size());

        
        //create weight image
        Mat weightImage;
        int w = image.cols;
        int h = image.rows;
        weightImage.create(h, w, CV_32FC1);
        float *p = (float*)weightImage.data;
        float x_center = w / 2;
        float y_center = h / 2;
        float dis_max = sqrt(x_center*x_center + y_center * y_center);
        int weightType = 0; // svar.GetInt("Map2D.WeightType", 0);
        //uchar* pbuffer = (uchar*)(image.data);
        for (int i = 0; i < h; i++)
        {
            uchar* pbuffer = image.ptr<uchar>(i);
            for (int j = 0; j < w; j++)
            {
                int r = pbuffer[j*3];
                int g = pbuffer[j*3 + 1];
                int b = pbuffer[j*3 + 2];
                if ((r + g + b) == 0)
                {
                    *p = 0;
                }
                else
                {
                    float dis = (i - y_center)*(i - y_center) + (j - x_center)*(j - x_center);
                    dis = 1 - sqrt(dis) / dis_max;
                    if (0 == weightType)
                        *p = dis;
                    else *p = dis * dis;
                    if (*p <= 1e-5) *p = 1e-5;
                }        
                p++;
            }
        }
        
        //create weight pyramid
        std::vector<cv::Mat> pyr_weights(nlevel + 1);
        pyr_weights[0] = weightImage;
        for (int i = 0; i < nlevel; ++i)
            cv::pyrDown(pyr_weights[i], pyr_weights[i + 1]);

        
        //fill the pyramid
        int l = (mValidDOMGeoData[fi].left - minx) / outResolution;
        int t = (maxy - mValidDOMGeoData[fi].top ) / outResolution;
        for (int i = 0; i <= nlevel; i++)
        {
            int pht = pyr_laplace[i].rows;
            int pwd = pyr_laplace[i].cols;
            printf("%d %d \n", pht, pwd);
            for (int m = 0; m < pht; m++)
            {
                //weight data
                int mt = m + t;
                if (mt >= oht) mt = oht - 1;
                short* outp = out_pyr_laplace[i].ptr<short>(mt);
                short* inp  = pyr_laplace[i].ptr<short>(m);

                //image data
                float* outw = out_pyr_weights[i].ptr<float>(mt);
                float* inw  = pyr_weights[i].ptr<float>(m);
                
                for (int n = 0; n < pwd; n++)
                {
                    //ignore the background
                    //int sg = inp[n * 3] + inp[n * 3 + 1] + inp[n * 3 + 2];
                    //if (sg == 0)  continue;
                    if (inw[n] == 0)
                        continue;

                    //update the image and weight
                    int nl = n + l;
                    if (nl >= owd) nl = owd - 1;
                    if (inw[n] > outw[nl])
                    {
                        outw[nl] = inw[n];
                        outp[(nl)*3]   = inp[n*3];
                        outp[(nl)*3+1] = inp[n*3+1];
                        outp[(nl)*3+2] = inp[n*3+2];
                    }
                }
            }
            //zoom
            l=0.5*l;
            t=0.5*t;
        }
    }  

    cv::detail::restoreImageFromLaplacePyr(out_pyr_laplace);
    Mat result = out_pyr_laplace[0].clone();
    if (result.type() == CV_16SC3)
        result.convertTo(result, CV_8UC3);
    result.setTo(cv::Scalar::all(0), out_pyr_weights[0] == 0);

    imwrite("c:\\temp\\fusion.jpg", result);
}


void  MainCall(void* c, double progress, string name)
{
    printf("%lf \n", progress);
    printf("%s \n", name);

    setTrackbarPos("progress", "mosaic", progress);

    cv::Mat image = cv::Mat::zeros(cv::Size(480, 240), CV_8UC3);
    image.setTo(cv::Scalar(100, 0, 0));
    int font_face = cv::FONT_HERSHEY_COMPLEX;
    double font_scale = 0.5;
    int thickness = 2;
    int baseline;
    cv::Size text_size = cv::getTextSize(name, font_face, font_scale, thickness, &baseline);
    cv::Point origin;
    origin.x = image.cols / 2 - text_size.width / 2;
    origin.y = image.rows / 2 + text_size.height / 2;
    cv::putText(image, name, origin, font_face, font_scale, cv::Scalar(0, 255, 255), thickness, 8, 0);
    cv::imshow("mosaic", image);
}

CSFMSystem sfm;


void mouseHandler(int event, int x, int y, int flags, void* param)
{
    switch (event)
    {
    case CV_EVENT_LBUTTONDOWN:

        sfm.Stop();

        printf("Left button down \n");

        //cvWaitKey();

        break;
    case CV_EVENT_LBUTTONUP:


        printf("Left button up \n");

        break;
    }
}

void TestMosaicThread()
{
    //string imagepath = "F:\\Data\\UAV\\byler\\1-sub";
    //string projectpath = "F:\\Data\\UAV\\byler\\1-sub\\work";
    //string imagepath = "F:\\Data\\UAV\\byler\\2-sub";
    //string projectpath = "F:\\Data\\UAV\\byler\\2-sub\\work";
    //string imagepath   = "C:\\Work\\Data\\uav\\heilongjiang\\images";
    //string projectpath = "C:\\Work\\Data\\uav\\heilongjiang\\work";
    string imagepath = "C:\\Work\\Data\\UAV\\zhong";
    string projectpath = "C:\\Work\\Data\\UAV\\zhong\\work";
    string mosaicfile = projectpath + "\\" + "mosaic.tif";
    
    //thread 1
    SFMParameters paras;
    paras.imagePath = imagepath;
    paras.projectPath = projectpath;
    paras.outfilePath = mosaicfile;
    paras.outResolution = 0.1;
    paras.camtype = PerspectiveCam;
    paras.maxHt = 640;

    //CSFMSystem sfm(paras);
    sfm.LoadParams(paras);
    sfm.RegCallBack(MainCall);
    sfm.Run();

    namedWindow("mosaic");
    int position = 0;
    createTrackbar("progress", "mosaic", &position, 100, 0);
    setMouseCallback("mosaic", mouseHandler);
    
    waitKey();

    printf("finished! \n");
}

int _tmain(int argc, char* argv[])
{
	printf("Test... \n");

    //TestMosaicThread();

    //MosaicTest();
    
    testfusion();

    //TestRT();

    //DllLaneDetecct();

    //LaneDetect1();    
	////test memory
	//vector<TrackInfo> tracks;
	//tracks.resize(5000000);
	//printf("%d \n", sizeof(tracks));
	//vector<TrackInfo> newtracks = tracks;
	//vector<TrackInfo>& newtracks1 = tracks;

    //testRotation();
    

	//ABSOriTest();

	//BAFreeNet();

	//TestPosConvertion();

	//BAWithControlPt();


	//TestCeres();

	//SimulateDataBA();

	//TestJYZLoss();

	//detection using yolo
	//TestJYZDetect(argc, argv);

	//TestGrabcut();

	//TestPos();
	
	//SimulateDataBA();

	//BAWithPos();

	//main_realimages(argc, argv);

	//invokeCeres();

	//test_gdal();

	//TestRand();
	//TestRand();
	//TestRand();
	
	//TestPano();
	//TestGeotiff();

	//if (argc > 1)
	//{
	//	//printf("%s \n", argv[0]);

	//	double aod = atof(argv[1]);
	//	TestRTUsingParasol(aod);
	//}

	printf("Finished ! \n");

	return 0;
}

