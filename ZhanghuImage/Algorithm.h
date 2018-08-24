#pragma once
#include"math.h"  
#include <complex>
#define PI 3.1415926535897932384626  
using namespace std;
class Algorithm
{
public:
	Algorithm(Bitmap* pCurBitmap);
	~Algorithm(void);

	Bitmap* m_pCurBitmap;
	BYTE*  m_pGrayData;
	int    m_Width;
	int    m_Height;


	bool ImageSmoothingSharping(double* pTemp, int sz, Bitmap*& pBitmap);
	bool HistEqual(Bitmap*& pBitmap);

	bool FFT2D(Bitmap*& pBitmap);

	bool RunlengthCode(CString filename);

	bool RunlengthRecode(CString filename, Bitmap*& pBitmap);

private:
	bool GrayData2Bitmap(BYTE* pGrayData, Bitmap*& pBitmap);

	bool Hist(double*& hist);

	bool RGBData2Bitmap(BYTE* pRGBData, Bitmap*& pBitmap);

	complex<double> * DataFitFormat(unsigned char *data, int lWidth, int lHeight);//将数组转换为适合FFT处理的数据（数组长度为2的整数次幂）,填充的数据补零操作.当lHeight=1时表示为对一维数组处理.data为对二维数据的一维表示，是按照从左到右，从上到下。 
	void InitTDAndFD(complex<double> *&TD, complex<double> *&FD, unsigned char *data, int lWidth, int lHeight, int &w, int &h);//初始化TD和FD(FFT操作)
	void FFT_1D(complex<double> *&TD, complex<double>*&FD, int Len);//一维FFT运算，len为一维数组的真实长度。而TD和FD数组的长度都是经过 InitTDAndFD得到的适合FFT处理的数组长度为2的整数次幂的数组 。
	void FFT_2D(complex<double>*&TD, complex<double> *&FD, int lWidth, int lHeight);//由一维FFT推算出二维FFT。lWidth，lHeight分别为要处理数据的宽和高。TD数组的长度为2的整数

};

