#include "stdafx.h"
#include "Algorithm.h"




Algorithm::Algorithm(Bitmap * pCurBitmap)
{
	m_pCurBitmap = pCurBitmap;
	m_Width = m_pCurBitmap->GetWidth();
	m_Height = m_pCurBitmap->GetHeight();

	BYTE* scan0;
	int stride;
	BitmapData* pBitmapData = new BitmapData();
	Rect rc(0, 0, m_Width, m_Height);
	m_pCurBitmap->LockBits(&rc, ImageLockModeRead | ImageLockModeWrite,
		PixelFormat24bppRGB, pBitmapData);
	scan0 = (BYTE*)pBitmapData->Scan0;
	stride = pBitmapData->Stride;

	m_pGrayData = new BYTE[m_Width * m_Height];

	for (int i = 0; i < m_Height; i++) {
		for (int j = 0; j < m_Width; j++) {
			m_pGrayData[i * m_Width + j] = (scan0[i * stride + 3 * j]
				+ scan0[i * stride + 3 * j + 1] + scan0[i * stride + 3 * j + 2]) / 3;
		}
	}


	m_pCurBitmap->UnlockBits(pBitmapData);

}

Algorithm::~Algorithm(void)
{
	if (NULL != m_pGrayData) {
		delete[] m_pGrayData;
	}

}

bool Algorithm::ImageSmoothingSharping(double * pTemp, int sz, Bitmap *& pBitmap)
{
	BYTE* pOlData = new BYTE[m_Height * m_Width];
	memcpy(pOlData, m_pGrayData, m_Height * m_Width * sizeof(BYTE));

	int i2, j2;
	int sz2 = 2 * sz + 1;
	double tmp;
	for (int i = 0; i < m_Height; i++) {
		for (int j = 0; j < m_Width; j++) {
			m_pGrayData[i * m_Width + j] = 0;
			tmp = 0;
			for (int i1 = -sz; i1 <= sz; i1++) {
				for (int j1 = -sz; j1 <= sz; j1++) {
					i2 = i + i1;
					j2 = j + j1;
					if (i2 < 0 || i2 >= m_Height || j2 < 0 || j2 >= m_Width)
						continue;

					tmp += (pOlData[i2 * m_Width + j2] * pTemp[(i1 + sz) * sz2 + j1 + sz]);
				}
			}
			if (tmp < 0)
				tmp = -tmp;
			m_pGrayData[i * m_Width + j] = (BYTE)tmp;

		}
	}
	GrayData2Bitmap(m_pGrayData, pBitmap);

	return true;
}

bool Algorithm::HistEqual(Bitmap *& pBitmap)
{
	double* hist = NULL;
	Hist(hist);
	BYTE* pGrayData = new BYTE[m_Height*m_Width];
	for (int i = 0; i<m_Height*m_Width; i++) {
		pGrayData[i] = BYTE(hist[m_pGrayData[i]] * 255);
	}
	GrayData2Bitmap(pGrayData, pBitmap);


	delete[] pGrayData;
	delete[] hist;
	return true;
}

bool Algorithm::FFT2D(Bitmap*& pBitmap) {
	complex<double> *TD = NULL;
	complex<double> *FD = NULL;
	int w;
	int h;
	complex<double> ft;
	double tmp;

	InitTDAndFD(TD, FD, m_pGrayData, m_Width, m_Height, w, h);

	FFT_2D(TD, FD, m_Width, m_Height);

	BYTE* pGrayData = new BYTE[m_Height * m_Width];
	double* pGrayData2 = new double[m_Height * m_Width];
	double maxgray = 0;
	int a = (w - m_Width) / 2;
	int b = (h - m_Height) / 2;
	for (int i = 0; i < m_Height; i++) {
		for (int j = 0; j < m_Width; j++)
		{
			ft = FD[(i + b)*w + (j + a)];
			tmp = sqrt(ft.real()*ft.real() + ft.imag()*ft.imag());

			pGrayData2[i*m_Width + j] = tmp;
			if (tmp > maxgray)
				maxgray = tmp;

		}
	}

	maxgray = maxgray / 10;
	for (int i = 0; i < m_Height * m_Width; i++) {
		pGrayData[i] = 255;
		if (pGrayData2[i] <= maxgray)
			pGrayData[i] = (BYTE)(pGrayData2[i] * 255 / maxgray);

	}
	GrayData2Bitmap(pGrayData, pBitmap);
	delete[] pGrayData;
	delete[] FD;
	delete[] TD;


	return true;

}

bool Algorithm::RunlengthCode(CString filename)
{
	CStdioFile listCon;

	listCon.Open(filename, CFile::modeWrite | CFile::typeBinary | CFile::modeCreate | CFile::modeNoTruncate);


	listCon.Write(&m_Width, sizeof(int));
	listCon.Write(&m_Height, sizeof(int));


	Color pixel;
	m_pCurBitmap->GetPixel(0, 0, &pixel);

	CString strR, strG, strB;
	strR = _T("");
	strG = _T("");
	strB = _T("");

	int k = 1;

	BYTE r, g, b;
	r = pixel.GetR();
	g = pixel.GetG();
	b = pixel.GetB();

	int rcnt = 1;
	int gcnt = 1;
	int bcnt = 1;

	int i, j;

	while (k < m_Width * m_Height) {
		i = k / m_Width;
		j = k % m_Width;
		m_pCurBitmap->GetPixel(j, i, &pixel);
		if (pixel.GetR() == r)
			rcnt++;
		else
		{
			listCon.Write(&r, 1);
			listCon.Write(&rcnt, sizeof(int));
			r = pixel.GetR();
			rcnt = 1;
		}

		if (k == m_Width * m_Height - 2) {
			k = k;
		}

		k++;
	}

	listCon.Write(&r, 1);
	listCon.Write(&rcnt, sizeof(int));

	k = 1;

	while (k < m_Width * m_Height) {
		i = k / m_Width;
		j = k % m_Width;
		m_pCurBitmap->GetPixel(j, i, &pixel);

		if (pixel.GetG() == g)
			gcnt++;
		else
		{
			listCon.Write(&g, 1);
			listCon.Write(&gcnt, sizeof(int));
			g = pixel.GetG();
			gcnt = 1;
		}


		k++;
	}
	listCon.Write(&g, 1);
	listCon.Write(&gcnt, sizeof(int));

	k = 1;

	while (k < m_Width * m_Height) {
		i = k / m_Width;
		j = k % m_Width;
		m_pCurBitmap->GetPixel(j, i, &pixel);

		if (pixel.GetB() == b)
			bcnt++;
		else
		{
			listCon.Write(&b, 1);
			listCon.Write(&bcnt, sizeof(int));
			b = pixel.GetB();
			bcnt = 1;
		}

		k++;
	}
	listCon.Write(&b, 1);
	listCon.Write(&bcnt, sizeof(int));

	listCon.Close();
	return true;

	/*
	CString str = _T("");

	CString str2;

	str2.Format(_T("%d,%d"),m_pCurBitmap->GetWidth(), m_pCurBitmap->GetHeight());

	str += str2;

	Color pixel;
	m_pCurBitmap->GetPixel(0,0,&pixel);

	CString strR, strG, strB;
	strR = _T("");
	strG = _T("");
	strB = _T("");

	int k = 1;

	BYTE r,g,b;
	r = pixel.GetR();
	g = pixel.GetG();
	b = pixel.GetB();

	int rcnt = 1;
	int gcnt = 1;
	int bcnt = 1;

	int i,j;

	while ( k < m_Width * m_Height){
	i  = k  / m_Width;
	j = k % m_Width;
	m_pCurBitmap->GetPixel(i,j,&pixel);
	if ( pixel.GetR() == r )
	rcnt++;
	else
	{
	str2 = "";
	str2.Format(_T("%d,%d;"),r,rcnt);
	strR += str2;
	r = pixel.GetR();
	rcnt = 1;
	}

	if ( pixel.GetG() == g )
	gcnt++;
	else
	{
	str2 = "";
	str2.Format(_T("%d,%d;"),g,gcnt);
	strG += str2;
	g = pixel.GetG();
	gcnt = 1;
	}

	if ( pixel.GetB() == b )
	bcnt++;
	else
	{
	str2 = "";
	str2.Format(_T("%d,%d;"),b,bcnt);
	strB += str2;
	b = pixel.GetB();
	bcnt = 1;
	}

	k++;
	}

	str += _T("R:")  + strR + _T("G:") + strG + _T("B:") + strB;

	CStdioFile listCon;

	listCon.Open(filename,CFile::modeWrite | CFile::typeText | CFile::modeCreate | CFile::modeNoTruncate);


	listCon.WriteString(str);


	return true;

	*/


}

bool Algorithm::RunlengthRecode(CString filename, Bitmap *& pBitmap)
{
	CStdioFile listCon;

	listCon.Open(filename, CFile::modeRead | CFile::typeBinary);


	listCon.Read(&m_Width, sizeof(int));
	listCon.Read(&m_Height, sizeof(int));

	BYTE* pData = new BYTE[m_Width*m_Height * 3];



	BYTE r, g, b;

	int rcnt;
	int gcnt;
	int bcnt;

	int i;
	int k = 0;
	while (true) {
		listCon.Read(&r, 1);
		listCon.Read(&rcnt, sizeof(int));
		for (i = 0; i < rcnt; i++)
		{
			pData[3 * (k + i)] = r;
		}
		k += rcnt;
		if (k == m_Width * m_Height - 2) {
			k = k;
		}

		if (k >= m_Width * m_Height)
			break;

	}

	k = 0;
	while (true) {
		listCon.Read(&g, 1);
		listCon.Read(&gcnt, sizeof(int));
		for (i = 0; i < gcnt; i++)
		{
			pData[3 * (k + i) + 1] = g;
		}
		k += gcnt;
		if (k >= m_Width * m_Height)
			break;

	}

	k = 0;
	while (true) {
		listCon.Read(&b, 1);
		listCon.Read(&bcnt, sizeof(int));
		for (i = 0; i < bcnt; i++)
		{
			pData[3 * (k + i) + 2] = b;
		}
		k += bcnt;
		if (k >= m_Width * m_Height)
			break;

	}


	this->RGBData2Bitmap(pData, pBitmap);

	listCon.Close();
	delete[] pData;





	return true;
}

bool Algorithm::GrayData2Bitmap(BYTE * pGrayData, Bitmap *& pBitmap)
{
	if (NULL != pBitmap) {
		delete pBitmap;
		pBitmap = NULL;
	}
	pBitmap = new Bitmap(m_Width, m_Height, PixelFormat24bppRGB);

	BYTE* scan0;
	int stride;
	BitmapData* pBitmapData = new BitmapData();
	Rect rc(0, 0, m_Width, m_Height);
	pBitmap->LockBits(&rc, ImageLockModeRead | ImageLockModeWrite,
		PixelFormat24bppRGB, pBitmapData);
	scan0 = (BYTE*)pBitmapData->Scan0;
	stride = pBitmapData->Stride;
	for (int i = 0; i < m_Height; i++) {
		for (int j = 0; j < m_Width; j++) {
			scan0[i * stride + 3 * j] = pGrayData[i * m_Width + j];
			scan0[i * stride + 3 * j + 1] = pGrayData[i * m_Width + j];
			scan0[i * stride + 3 * j + 2] = pGrayData[i * m_Width + j];
		}
	}
	pBitmap->UnlockBits(pBitmapData);
	return true;
}

bool Algorithm::Hist(double *& hist)
{
	if (NULL != hist) {
		delete[] hist;
	}
	const int dim = 256;
	hist = new double[dim];
	memset(hist, 0, dim * sizeof(double));
	for (int i = 0; i < m_Height; i++) {
		for (int j = 0; j < m_Width; j++) {
			hist[m_pGrayData[i * m_Width + j]]++;
		}
	}
	for (int i = 0; i<dim; i++) {
		hist[i] = hist[i] / (m_Height*m_Width);
	}
	for (int i = 1; i<dim; i++) {
		hist[i] = hist[i] + hist[i - 1];
	}
	return true;
}

bool Algorithm::RGBData2Bitmap(BYTE * pRGBData, Bitmap *& pBitmap)
{
	if (NULL == pRGBData)
		return false;

	if (NULL != pBitmap) {
		delete pBitmap;
		pBitmap = NULL;
	}
	pBitmap = new Bitmap(m_Width, m_Height, PixelFormat24bppRGB);

	BYTE* scan0;
	int stride;
	BitmapData* pBitmapData = new BitmapData();
	Rect rc(0, 0, m_Width, m_Height);
	pBitmap->LockBits(&rc, ImageLockModeRead | ImageLockModeWrite,
		PixelFormat24bppRGB, pBitmapData);
	scan0 = (BYTE*)pBitmapData->Scan0;
	stride = pBitmapData->Stride;
	for (int i = 0; i < m_Height; i++) {
		for (int j = 0; j < m_Width; j++) {
			scan0[i * stride + 3 * j + 2] = pRGBData[3 * i * m_Width + 3 * j];
			scan0[i * stride + 3 * j + 1] = pRGBData[3 * i * m_Width + 3 * j + 1];
			scan0[i * stride + 3 * j] = pRGBData[3 * i * m_Width + 3 * j + 2];
		}
	}
	pBitmap->UnlockBits(pBitmapData);



	return true;
}

complex<double>* Algorithm::DataFitFormat(unsigned char * data, int lWidth, int lHeight)
{
	complex<double> *TD;
	int w = 1;
	int h = 1;
	int wp = 0;//存储w的2的幂数
	int hp = 0;//存储h的2的幂数
			   //////计算刚好大于或等于lWidth，lHeight的2的整数次幂，和相应的幂数///////////////
	while (w<lWidth)
	{
		w = w * 2;
		wp++;
	}
	while (h<lHeight)
	{
		h = h * 2;
		hp++;
	}
	TD = new complex<double>[w*h];
	////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i<h; i++)
	{
		if (i<lHeight)
		{
			for (int j = 0; j<w; j++)
			{
				if (j< lWidth)
				{
					TD[i*w + j] = complex<double>(data[i*lWidth + j], 0);
					if ((i + j) % 2 != 0)
						TD[i*w + j] = complex<double>(data[i*lWidth + j] * (-1), 0);//将char数据，准换为实数为data数据，虚数为0的复数
				}
				else
				{
					TD[i*w + j] = complex<double>(0, 0);//对于超出原数据的数据进行补零操作
				}

			}
		}
		else
		{
			for (int j = 0; j<w; j++)
			{
				TD[i*w + j] = complex<double>(0, 0);//对于超出原数据的数据进行补零操作

			}
		}
	}

	return TD;
}

void Algorithm::InitTDAndFD(complex<double>*& TD, complex<double>*& FD, unsigned char * data, int lWidth, int lHeight, int & w, int & h)
{
	w = 1;
	h = 1;
	int wp = 0;//存储w的2的幂数
	int hp = 0;//存储h的2的幂数
	complex<double> *TmpFD;
	//////计算刚好大于或等于lWidth，lHeight的2的整数次幂，和相应的幂数///////////////
	while (w<lWidth)
	{
		w = w * 2;
		wp++;
	}
	while (h<lHeight)
	{
		h = h * 2;
		hp++;
	}

	TmpFD = new complex<double>[w*h];

	for (int i = 0; i<w*h; i++)
	{
		TmpFD[i] = complex<double>(0.0, 0.0);//FD初值设为0
	}


	TD = DataFitFormat(data, lWidth, lHeight);//调用已有函数DataFitFormat初始化TD
	FD = TmpFD;

}

void Algorithm::FFT_1D(complex<double>*& TD, complex<double>*& FD, int Len)
{
	//long i,j,k;
	int l = 1;
	int lp = 0;
	int p = 0;
	double angle = 0;//中间变量及角度
	complex<double> *W, *X1, *X2, *X;

	while (l<Len)
	{
		l = l * 2;
		lp++;
	}
	int r = lp;

	long N = 1 << r;//快速傅里叶变换点数 2的r次幂;
	W = new complex<double>[N / 2];//存放旋转系数
	X1 = new complex<double>[N];//
	X2 = new complex<double>[N];//分配运算的存储器
	for (long i = 0; i<N / 2; i++)
	{
		angle = -i * PI * 2 / N;
		W[i] = complex<double>(cos(angle), sin(angle));
	}

	memcpy(X1, TD, sizeof(complex<double>)*N);//将TD所在的内存数据拷贝到X1中

											  ///////////////////////////核心算法：蝶形运算/////////////
	for (long k = 0; k<r; k++)
	{
		for (long j = 0; j<(1 << k); j++)
		{
			for (long i = 0; i<(1 << (r - k - 1)); i++)
			{
				p = j * (1 << (r - k));
				X2[i + p] = X1[i + p] + X1[i + p + (int)(1 << (r - k - 1))];
				X2[i + p + (int)(1 << (r - k - 1))] = (X1[i + p] - X1[i + p + (int)(1 << (r - k - 1))]) *W[i*(1 << k)];
			}
		}
		X = X1;
		X1 = X2;
		X2 = X;
	}

	/////////////////重新排序，将反序->正序//////////////
	for (int j = 0; j<N; j++)
	{
		p = 0;
		for (int i = 0; i<r; i++)
		{
			if (j&(1 << i))
			{
				p += 1 << (r - i - 1);
			}
		}
		FD[j] = X1[p];
	}
	delete W;
	delete X1;
	delete X2;

}

void Algorithm::FFT_2D(complex<double>*& TD, complex<double>*& FD, int lWidth, int lHeight)
{
	int w = 1;
	int h = 1;
	int wp = 0;//存储w的2的幂数
	int hp = 0;//存储h的2的幂数
	complex<double> *TmpTD, *TmpFD;//存放临时的一列的数据

								   //////计算刚好大于或等于lWidth，lHeight的2的整数次幂，和相应的幂数///////////////
	while (w<lWidth)
	{
		w = w * 2;
		wp++;
	}
	while (h<lHeight)
	{
		h = h * 2;
		hp++;
	}
	////////////////////////先按y方向进行快速的一维FFT运算 //////////////////////////////////
	TmpTD = new complex<double>[h];
	TmpFD = new complex<double>[h];

	for (int i = 0; i<w; i++)
	{
		//先按y方向进行快速的一维FFT运算 
		for (int j = 0; j<h; j++)
		{
			TmpTD[j] = TD[j*w + i];
		}
		FFT_1D(TmpTD, TmpFD, lHeight);
		//保存结果
		for (int j = 0; j<h; j++)
		{
			TD[j*w + i] = TmpFD[j];
		}
	}
	delete[] TmpTD;
	delete[] TmpFD;
	///////////////////////再按x方向进行快速的一维FFT运算///////////////////////////
	TmpTD = new complex<double>[w];
	TmpFD = new complex<double>[w];

	for (int i = 0; i<h; i++)
	{
		////再按x方向进行快速的一维FFT运算
		for (int j = 0; j<w; j++)
		{
			TmpTD[j] = TD[i*w + j];
		}
		FFT_1D(TmpTD, TmpFD, lWidth);
		//保存结果
		for (int j = 0; j<w; j++)
		{
			FD[i*w + j] = TmpFD[j];
		}
	}
	delete[] TmpTD;
	delete[] TmpFD;

}

