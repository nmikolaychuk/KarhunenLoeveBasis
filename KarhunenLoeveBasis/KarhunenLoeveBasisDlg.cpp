
// KarhunenLoeveBasisDlg.cpp: файл реализации
//

#include "pch.h"
#include "framework.h"
#include "KarhunenLoeveBasis.h"
#include "KarhunenLoeveBasisDlg.h"
#include "afxdialogex.h"

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <random>
#include <cassert>
#include <string.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define DOTS(x,y) (xp*((x)-xmin)),(yp*((y)-ymax))
#define DOTSSING(x,y) (xp_sing*((x)-xmin_sing)),(yp_sing*((y)-ymax_sing))
#define DOTSV(x,y) (xp_v*((x)-xmin_v)),(yp_v*((y)-ymax_v))

using namespace std;

// Диалоговое окно CKarhunenLoeveBasisDlg

struct PRNG
{
	std::mt19937 engine;
};

void initGenerator(PRNG& generator)
{
	// Создаём псевдо-устройство для получения случайного зерна.
	std::random_device device;
	// Получаем случайное зерно последовательности
	generator.engine.seed(device());
}

CKarhunenLoeveBasisDlg::CKarhunenLoeveBasisDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_KARHUNENLOEVEBASIS_DIALOG, pParent)
//part 1
	, result_part1(_T(""))
	, nevyaz_obr(0)
	, nevyaz_psevdo(0)
	, determ(0)
	, dim_m(0)
	, dim_n(0)
//part 2	
	, ampl1(0.05)
	, ampl2(0.01)
	, ampl3(0.02)
	, omega1(0.05)
	, omega2(0.04)
	, omega3(0.00006)
	, fi1(0)
	, fi2(0)
	, fi3(0.4)
	, sigma1(25.0)
	, sigma2(15.0)
	, sigma3(10.5)
	, t1(80.0)
	, t2(250.0)
	, t3(410.0)
	, Length(512)
	, dim_akm(50)
	, number_of_vector(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CKarhunenLoeveBasisDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_RESULT_TEXT, result_part1);
	DDX_Text(pDX, IDC_EDIT_NEV_OBR_MATR, nevyaz_obr);
	DDX_Text(pDX, IDC_EDIT_NEV_PSEVDO_MATR, nevyaz_psevdo);
	DDX_Text(pDX, IDC_EDIT_DETERMINANT, determ);
	DDX_Text(pDX, IDC_EDIT_RAZMER_M, dim_m);
	DDX_Text(pDX, IDC_EDIT_RAZMER_N, dim_n);
	DDX_Text(pDX, IDC_EDIT_A1, ampl1);
	DDX_Text(pDX, IDC_EDIT_A2, ampl2);
	DDX_Text(pDX, IDC_EDIT_A3, ampl3);
	DDX_Text(pDX, IDC_EDIT_OMEGA_1, omega1);
	DDX_Text(pDX, IDC_EDIT_OMEGA_2, omega2);
	DDX_Text(pDX, IDC_EDIT_OMEGA_3, omega3);
	DDX_Text(pDX, IDC_EDIT_FI_1, fi1);
	DDX_Text(pDX, IDC_EDIT_FI_2, fi2);
	DDX_Text(pDX, IDC_EDIT_FI_3, fi3);
	DDX_Text(pDX, IDC_EDIT_SIGMA_1, sigma1);
	DDX_Text(pDX, IDC_EDIT_SIGMA_2, sigma2);
	DDX_Text(pDX, IDC_EDIT_SIGMA_3, sigma3);
	DDX_Text(pDX, IDC_EDIT_T1, t1);
	DDX_Text(pDX, IDC_EDIT_T2, t2);
	DDX_Text(pDX, IDC_EDIT_T3, t3);
	DDX_Text(pDX, IDC_EDIT_LENGTH, Length);
	DDX_Text(pDX, IDC_EDIT_RAZM_AKM, dim_akm);
	DDX_Control(pDX, IDC_RADIO_POLIGARM, m_radio_poligarm_signal);
	DDX_Control(pDX, IDC_RADIO_GAUSS, m_radio_gauss_signal);
	DDX_Control(pDX, IDC_RADIO_EXP, m_radio_exp_signal);
	//  DDX_Text(pDX, IDC_EDIT_COUNT_VECTOR, number_of_vector);
	DDX_Text(pDX, IDC_EDIT_COUNT_VECTOR, number_of_vector);
}

BEGIN_MESSAGE_MAP(CKarhunenLoeveBasisDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_SVD_PART1, &CKarhunenLoeveBasisDlg::OnBnClickedButtonSvdPart1)
	ON_BN_CLICKED(IDC_BUTTON_SIGNAL, &CKarhunenLoeveBasisDlg::OnBnClickedButtonSignal)
END_MESSAGE_MAP()


// Обработчики сообщений CKarhunenLoeveBasisDlg

BOOL CKarhunenLoeveBasisDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Задает значок для этого диалогового окна.  Среда делает это автоматически,
	//  если главное окно приложения не является диалоговым
	SetIcon(m_hIcon, TRUE);			// Крупный значок
	SetIcon(m_hIcon, FALSE);		// Мелкий значок

	// TODO: добавьте дополнительную инициализацию
	PicWnd = GetDlgItem(IDC_DRAW_SIGNAL);			//связываем с ID окон
	PicDc = PicWnd->GetDC();
	PicWnd->GetClientRect(&Pic);

	PicWndSingular = GetDlgItem(IDC_DRAW_SINGULAR);			//связываем с ID окон
	PicDcSingular = PicWndSingular->GetDC();
	PicWndSingular->GetClientRect(&PicSingular);

	PicWndV = GetDlgItem(IDC_DRAW_SING_FUNC);			//связываем с ID окон
	PicDcV = PicWndV->GetDC();
	PicWndV->GetClientRect(&PicV);

	// перья
	setka_pen.CreatePen(		//для сетки
		PS_DOT,					//пунктирная
		-1,						//толщина 1 пиксель
		RGB(0, 0, 0));			//цвет  черный

	osi_pen.CreatePen(			//координатные оси
		PS_SOLID,				//сплошная линия
		2,						//толщина 2 пикселя
		RGB(0, 0, 0));			//цвет черный

	signal_pen.CreatePen(			//график
		PS_SOLID,				//сплошная линия
		-1,						//толщина -1 пикселя
		RGB(255, 0, 0));			//цвет синий

	singular_value_pen.CreatePen(			//график
		PS_SOLID,				//сплошная линия
		-1,						//толщина -1 пикселя
		RGB(0, 0, 255));			//цвет синий

	own_vector_pen.CreatePen(			//график
		PS_SOLID,				//сплошная линия
		-1,						//толщина -1 пикселя
		RGB(200, 0, 255));			//цвет синий

	RedrawSignal(ymin, ymax);
	RedrawSingular(ymin_sing, ymax_sing);
	RedrawOwnVector(ymin_v, ymax_v);

	UpdateData(FALSE);
	return TRUE;  // возврат значения TRUE, если фокус не передан элементу управления
}

// При добавлении кнопки свертывания в диалоговое окно нужно воспользоваться приведенным ниже кодом,
//  чтобы нарисовать значок.  Для приложений MFC, использующих модель документов или представлений,
//  это автоматически выполняется рабочей областью.

void CKarhunenLoeveBasisDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // контекст устройства для рисования

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Выравнивание значка по центру клиентского прямоугольника
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Нарисуйте значок
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// Система вызывает эту функцию для получения отображения курсора при перемещении
//  свернутого окна.
HCURSOR CKarhunenLoeveBasisDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CKarhunenLoeveBasisDlg::Mashtab(double arr[], int dim, double* mmin, double* mmax)		//определяем функцию масштабирования
{
	*mmin = *mmax = arr[0];

	for (int i = 0; i < dim; i++)
	{
		if (*mmin > arr[i]) *mmin = arr[i];
		if (*mmax < arr[i]) *mmax = arr[i];
	}
}

void CKarhunenLoeveBasisDlg::RedrawSignal(double ymin, double ymax)
{
	UpdateData(TRUE);
	PicDc->FillSolidRect(&Pic, RGB(250, 250, 250));			//закрашиваю фон 
	PicDc->SelectObject(&osi_pen);		//выбираем перо

	double window_signal_width = Pic.Width();
	double window_signal_height = Pic.Height();
	xp = (window_signal_width / (xmax - xmin));			//Коэффициенты пересчёта координат по Х
	yp = -(window_signal_height / (ymax - ymin));			//Коэффициенты пересчёта координат по У

	//создаём Ось Y
	PicDc->MoveTo(DOTS(0, ymax));
	PicDc->LineTo(DOTS(0, ymin));
	//создаём Ось Х
	PicDc->MoveTo(DOTS(xmin, 0));
	PicDc->LineTo(DOTS(xmax, 0));

	PicDc->SelectObject(&setka_pen);

	//отрисовка сетки по y
	for (double x = 0; x <= xmax; x += Length / 10)
	{
		if (x != 0) {
			PicDc->MoveTo(DOTS(x, ymax));
			PicDc->LineTo(DOTS(x, ymin));
		}
	}
	//отрисовка сетки по x
	for (double y = -ymax; y < ymax; y += ymax / 5)
	{
		if (y != 0) {
			PicDc->MoveTo(DOTS(xmin, y));
			PicDc->LineTo(DOTS(xmax, y));
		}
	}


	//подпись точек на оси
	CFont font;
	font.CreateFontW(14.5, 0, 0, 0, FW_REGULAR, 0, 0, 0, ANSI_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS || CLIP_LH_ANGLES, DEFAULT_QUALITY, DEFAULT_PITCH, _T("Century Gothic"));
	PicDc->SelectObject(font);

	//подпись осей
	PicDc->TextOutW(DOTS(0.02 * Length, 0.98 * ymax), _T("x")); //Y
	PicDc->TextOutW(DOTS(xmax - Length / 15, 0.2 * ymax), _T("t")); //X

	//по Y с шагом 5
	for (double i = -ymax; i <= ymax; i += ymax / 5)
	{
		CString str;
		if (i != 0)
		{
			str.Format(_T("%.3f"), i);
			PicDc->TextOutW(DOTS(xmin + Length / 50, i + 0.03 * ymax), str);
		}
	}
	//по X с шагом 0.5
	for (double j = 0; j <= xmax; j += Length / 10)
	{
		CString str;
		if (j != 0) {
			str.Format(_T("%.f"), j);
			PicDc->TextOutW(DOTS(j - Length / 100, 0), str);
		}
	}
}

void CKarhunenLoeveBasisDlg::RedrawSingular(double ymin_sing, double ymax_sing)
{
	UpdateData(TRUE);
	PicDcSingular->FillSolidRect(&PicSingular, RGB(250, 250, 250));			//закрашиваю фон 
	PicDcSingular->SelectObject(&osi_pen);		//выбираем перо

	double window_signal_width = PicSingular.Width();
	double window_signal_height = PicSingular.Height();
	xp_sing = (window_signal_width / (xmax_sing - xmin_sing));			//Коэффициенты пересчёта координат по Х
	yp_sing = -(window_signal_height / (ymax_sing - ymin_sing));			//Коэффициенты пересчёта координат по У

	//создаём Ось Y
	PicDcSingular->MoveTo(DOTSSING(0, ymax_sing));
	PicDcSingular->LineTo(DOTSSING(0, ymin_sing));
	//создаём Ось Х
	PicDcSingular->MoveTo(DOTSSING(xmin_sing, 0));
	PicDcSingular->LineTo(DOTSSING(xmax_sing, 0));

	PicDcSingular->SelectObject(&setka_pen);

	//отрисовка сетки по y
	for (double x = 0; x <= xmax_sing; x += dim_akm / 10)
	{
		if (x != 0) {
			PicDcSingular->MoveTo(DOTSSING(x, ymax_sing));
			PicDcSingular->LineTo(DOTSSING(x, ymin_sing));
		}
	}
	//отрисовка сетки по x
	for (double y = 0; y < ymax_sing; y += ymax_sing / 9)
	{
		if (y != 0) {
			PicDcSingular->MoveTo(DOTSSING(xmin_sing, y));
			PicDcSingular->LineTo(DOTSSING(xmax_sing, y));
		}
	}


	//подпись точек на оси
	CFont font;
	font.CreateFontW(14.5, 0, 0, 0, FW_REGULAR, 0, 0, 0, ANSI_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS || CLIP_LH_ANGLES, DEFAULT_QUALITY, DEFAULT_PITCH, _T("Century Gothic"));
	PicDcSingular->SelectObject(font);

	//подпись осей
	PicDcSingular->TextOutW(DOTSSING(0.02 * dim_akm, 0.98 * ymax_sing), _T("Sigma")); //Y
	PicDcSingular->TextOutW(DOTSSING(xmax_sing - dim_akm / 15, -ymin_sing - 0.05 * ymax_sing), _T("N")); //X

	//по Y с шагом 5
	for (double i = 0; i <= ymax_sing; i += ymax_sing / 9)
	{
		CString str;
		if (i != 0)
		{
			str.Format(_T("%.3f"), i);
			PicDcSingular->TextOutW(DOTSSING(xmin_sing + dim_akm / 50, i + 0.03 * ymax_sing), str);
		}
	}
	//по X с шагом 0.5
	for (double j = 0; j <= xmax_sing; j += dim_akm / 10)
	{
		CString str;
		if (j != 0) {
			str.Format(_T("%.f"), j);
			PicDcSingular->TextOutW(DOTSSING(j - dim_akm * 0.01, -ymin_sing - 0.16 * ymax_sing), str);
		}
	}
}

void CKarhunenLoeveBasisDlg::RedrawOwnVector(double ymin_v, double ymax_v)
{
	UpdateData(TRUE);
	PicDcV->FillSolidRect(&PicV, RGB(250, 250, 250));			//закрашиваю фон 
	PicDcV->SelectObject(&osi_pen);		//выбираем перо

	double window_signal_width = PicV.Width();
	double window_signal_height = PicV.Height();
	xp_v = (window_signal_width / (xmax_v - xmin_v));			//Коэффициенты пересчёта координат по Х
	yp_v = -(window_signal_height / (ymax_v - ymin_v));			//Коэффициенты пересчёта координат по У

	//создаём Ось Y
	PicDcV->MoveTo(DOTSV(0, ymax_v));
	PicDcV->LineTo(DOTSV(0, ymin_v));
	//создаём Ось Х
	PicDcV->MoveTo(DOTSV(xmin_v, 0));
	PicDcV->LineTo(DOTSV(xmax_v, 0));

	PicDcV->SelectObject(&setka_pen);

	//отрисовка сетки по y
	for (double x = 0; x <= xmax_v; x += dim_akm / 10)
	{
		if (x != 0) {
			PicDcV->MoveTo(DOTSV(x, ymax_v));
			PicDcV->LineTo(DOTSV(x, ymin_v));
		}
	}
	//отрисовка сетки по x
	for (double y = -ymax_v; y < ymax_v; y += ymax_v / 5)
	{
		if (y != 0) {
			PicDcV->MoveTo(DOTSV(xmin_v, y));
			PicDcV->LineTo(DOTSV(xmax_v, y));
		}
	}


	//подпись точек на оси
	CFont font;
	font.CreateFontW(14.5, 0, 0, 0, FW_REGULAR, 0, 0, 0, ANSI_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS || CLIP_LH_ANGLES, DEFAULT_QUALITY, DEFAULT_PITCH, _T("Century Gothic"));
	PicDcV->SelectObject(font);

	//подпись осей
	PicDcV->TextOutW(DOTSV(0.02 * dim_akm, 0.98 * ymax_v), _T("Sigma")); //Y
	PicDcV->TextOutW(DOTSV(xmax_v - dim_akm / 15, -ymin_v - 0.05 * ymax_v), _T("N")); //X

	//по Y с шагом 5
	for (double i = -ymax_v; i <= ymax_v; i += ymax_v / 5)
	{
		CString str;
		if (i != 0)
		{
			str.Format(_T("%.3f"), i);
			PicDcV->TextOutW(DOTSV(xmin_v + dim_akm / 50, i + 0.03 * ymax_v), str);
		}
	}
	//по X с шагом 0.5
	for (double j = 0; j <= xmax_v; j += dim_akm / 10)
	{
		CString str;
		if (j != 0) {
			str.Format(_T("%.f"), j);
			PicDcV->TextOutW(DOTSV(j - dim_akm * 0.01, 0 - 0.05 * ymax_v), str);
		}
	}
}

int CKarhunenLoeveBasisDlg::svd_hestenes(int m_m, int n_n, double* a, double* u, double* v, double* sigma)
{
	double thr = 1.E-4f, nul = 1.E-16f;
	int n, m, i, j, l, k, lort, iter, in, ll, kk;
	double alfa, betta, hamma, eta, t, cos0, sin0, buf, s;
	n = n_n;
	m = m_m;
	for (i = 0; i < n; i++)
	{
		in = i * n;
		for (j = 0; j < n; j++)
			if (i == j) v[in + j] = 1.;
			else v[in + j] = 0.;
	}
	for (i = 0; i < m; i++)
	{
		in = i * n;
		for (j = 0; j < n; j++)
		{
			u[in + j] = a[in + j];
		}
	}

	iter = 0;
	while (1)
	{
		lort = 0;
		iter++;
		for (l = 0; l < n - 1; l++)
			for (k = l + 1; k < n; k++)
			{
				alfa = 0.; betta = 0.; hamma = 0.;
				for (i = 0; i < m; i++)
				{
					in = i * n;
					ll = in + l;
					kk = in + k;
					alfa += u[ll] * u[ll];
					betta += u[kk] * u[kk];
					hamma += u[ll] * u[kk];
				}

				if (sqrt(alfa * betta) < nul)	continue;
				if (fabs(hamma) / sqrt(alfa * betta) < thr) continue;

				lort = 1;
				eta = (betta - alfa) / (2.f * hamma);
				t = (double)((eta / fabs(eta)) / (fabs(eta) + sqrt(1. + eta * eta)));
				cos0 = (double)(1. / sqrt(1. + t * t));
				sin0 = t * cos0;

				for (i = 0; i < m; i++)
				{
					in = i * n;
					buf = u[in + l] * cos0 - u[in + k] * sin0;
					u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
					u[in + l] = buf;

					if (i >= n) continue;
					buf = v[in + l] * cos0 - v[in + k] * sin0;
					v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
					v[in + l] = buf;
				}
			}

		if (!lort) break;
	}

	for (i = 0; i < n; i++)
	{
		s = 0.;
		for (j = 0; j < m; j++)	s += u[j * n + i] * u[j * n + i];
		s = (double)sqrt(s);
		sigma[i] = s;
		if (s < nul)	continue;
		for (j = 0; j < m; j++)	u[j * n + i] /= s;
	}
	//======= Sortirovka ==============
	for (i = 0; i < n - 1; i++)
		for (j = i; j < n; j++)
			if (sigma[i] < sigma[j])
			{
				s = sigma[i]; sigma[i] = sigma[j]; sigma[j] = s;
				for (k = 0; k < m; k++)
				{
					s = u[i + k * n]; u[i + k * n] = u[j + k * n]; u[j + k * n] = s;
				}
				for (k = 0; k < n; k++)
				{
					s = v[i + k * n]; v[i + k * n] = v[j + k * n]; v[j + k * n] = s;
				}
			}

	return iter;
}

double getRandomdouble(PRNG& generator, double minValue, double maxValue)
{
	// Проверяем корректность аргументов
	assert(minValue < maxValue);

	// Создаём распределение
	std::uniform_real_distribution<double> distribution(minValue, maxValue);

	// Вычисляем псевдослучайное число: вызовем распределение как функцию,
	//  передав генератор произвольных целых чисел как аргумент.
	return distribution(generator.engine);
}

// ----------------------------- Определитель --------------------------
// Вычисляет определитель матрицы T размерностью N
double CKarhunenLoeveBasisDlg::det(double** T, UINT32 N)
{
	double det__;
	int sub_j, s;
	double** subT;    // Субматрица как набор ссылок на исходную матрицу
	switch (N)
	{
	case 1:
		return T[0][0];
	case 2:
		return T[0][0] * T[1][1] - T[0][1] * T[1][0];
	default:
		if (N < 1)  // Некорректная матрица
		{
			MessageBox(L"Некорректная матрица!", L"Ошибка", MB_ICONERROR | MB_OK);
		}
		subT = new double* [N - 1];  // Массив ссылок на столбцы субматрицы
		det__ = 0;
		s = 1;        // Знак минора
		for (int i = 0; i < N; i++)  // Разложение по первому столбцу
		{
			sub_j = 0;
			for (int j = 0; j < N; j++)// Заполнение субматрицы ссылками на исходные столбцы
				if (i != j)      // исключить i строку
					subT[sub_j++] = T[j] + 1;  // здесь + 1 исключает первый столбец

			det__ = det__ + s * T[i][0] * det(subT, N - 1);
			s = -s;
		};
		delete[] subT;
		return det__;
	};
}

// Функция, производящая обращение матрицы.
// Принимает:
//     matrix - матрица для обращения
//     result - матрица достаточного размера для вмещения результата
//     size   - размерность матрицы
// Возвращает:
//     true в случае успешного обращения, false в противном случае
bool inverse(double** matrix, double** result, int size)
{
	// Изначально результирующая матрица является единичной
	// Заполняем единичную матрицу
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
			result[i][j] = 0.0;

		result[i][i] = 1.0;
	}

	// Копия исходной матрицы
	double** copy = new double* [size]();

	// Заполняем копию исходной матрицы
	for (int i = 0; i < size; ++i)
	{
		copy[i] = new double[size];

		for (int j = 0; j < size; ++j)
			copy[i][j] = matrix[i][j];
	}

	// Проходим по строкам матрицы (назовём их исходными)
	// сверху вниз. На данном этапе происходит прямой ход
	// и исходная матрица превращается в верхнюю треугольную
	for (int k = 0; k < size; ++k)
	{
		// Если элемент на главной диагонали в исходной
		// строке - нуль, то ищем строку, где элемент
		// того же столбца не нулевой, и меняем строки
		// местами
		if (fabs(copy[k][k]) < 1e-8)
		{
			// Ключ, говорязий о том, что был произведён обмен строк
			bool changed = false;

			// Идём по строкам, расположенным ниже исходной
			for (int i = k + 1; i < size; ++i)
			{
				// Если нашли строку, где в том же столбце
				// имеется ненулевой элемент
				if (fabs(copy[i][k]) > 1e-8)
				{
					// Меняем найденную и исходную строки местами
					// как в исходной матрице, так и в единичной
					std::swap(copy[k], copy[i]);
					std::swap(result[k], result[i]);

					// Взводим ключ - сообщаем о произведённом обмене строк
					changed = true;

					break;
				}
			}

			// Если обмен строк произведён не был - матрица не может быть
			// обращена
			if (!changed)
			{
				// Чистим память
				for (int i = 0; i < size; ++i)
					delete[] copy[i];

				delete[] copy;

				// Сообщаем о неудаче обращения
				return false;
			}
		}

		// Запоминаем делитель - диагональный элемент
		double div = copy[k][k];

		// Все элементы исходной строки делим на диагональный
		// элемент как в исходной матрице, так и в единичной
		for (int j = 0; j < size; ++j)
		{
			copy[k][j] /= div;
			result[k][j] /= div;
		}

		// Идём по строкам, которые расположены ниже исходной
		for (int i = k + 1; i < size; ++i)
		{
			// Запоминаем множитель - элемент очередной строки,
			// расположенный под диагональным элементом исходной
			// строки
			double multi = copy[i][k];

			// Отнимаем от очередной строки исходную, умноженную
			// на сохранённый ранее множитель как в исходной,
			// так и в единичной матрице
			for (int j = 0; j < size; ++j)
			{
				copy[i][j] -= multi * copy[k][j];
				result[i][j] -= multi * result[k][j];
			}
		}
	}

	// Проходим по вернхней треугольной матрице, полученной
	// на прямом ходе, снизу вверх
	// На данном этапе происходит обратный ход, и из исходной
	// матрицы окончательно формируется единичная, а из единичной -
	// обратная
	for (int k = size - 1; k > 0; --k)
	{
		// Идём по строкам, которые расположены выше исходной
		for (int i = k - 1; i + 1 > 0; --i)
		{
			// Запоминаем множитель - элемент очередной строки,
			// расположенный над диагональным элементом исходной
			// строки
			double multi = copy[i][k];

			// Отнимаем от очередной строки исходную, умноженную
			// на сохранённый ранее множитель как в исходной,
			// так и в единичной матрице
			for (int j = 0; j < size; ++j)
			{
				copy[i][j] -= multi * copy[k][j];
				result[i][j] -= multi * result[k][j];
			}
		}
	}

	// Чистим память
	for (int i = 0; i < size; ++i)
		delete[] copy[i];

	delete[] copy;

	// Сообщаем об успехе обращения
	return true;
}

void CKarhunenLoeveBasisDlg::OnBnClickedButtonSvdPart1()
{
	// TODO: добавьте свой код обработчика уведомлений
	UpdateData(TRUE);

	PRNG generator;
	initGenerator(generator);

	double** A = new double*[size];
	for (int i = 0; i < size; i++)
	{
		A[i] = new double[size];
	}

	double** Aobr = new double* [size];
	for (int i = 0; i < size; i++)
	{
		Aobr[i] = new double[size];
	}
	
	double* b = new double[size];
	double* x = new double[size];
	double* x_psevdo = new double[size];

	for (int i = 0; i < size; i++)
	{
		b[i] = 0;
		x[i] = 0;
		x_psevdo[i] = 0;
		for (int j = 0; j < size; j++)
		{
			A[i][j] = 0;
			Aobr[i][j] = 0;
		}
	}

	A[0][0] = 1.01;
	A[0][1] = 2.01;
	A[0][2] = 3.01;
	A[1][0] = 4.01;
	A[1][1] = 5.01;
	A[1][2] = 6.01;
	A[2][0] = 7.01;
	A[2][1] = 8.01;
	A[2][2] = 9.01;

	for (int i = 0; i < size; i++)
	{
		b[i] = getRandomdouble(generator, -1, 1);
	}

	ofstream out("result_part_1.txt");

	out << "\t      ИНФОРМАЦИОННЫЕ ТЕХНОЛОГИИ\n";
	out << "Задание 3. Построение и визуализация базиса Каруена-Лоэва. Часть 1.";

	out << "\n\nМатрица коэффициентов А:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << A[i][j] << "\t";
		}
		out << "\n";
	}

	double determinant = det(A, size);

	out << "\nОпределитель матрицы: " << setprecision(3) << determinant;

	out << "\n\nВектор свободных членов: ";
	for (int i = 0; i < size; i++)
	{
		out << setprecision(3) << b[i] << "      ";
	}

	out << "\n\n********* ОБРАТНАЯ МАТРИЦА (КЛАССИЧЕСКИЙ МЕТОД) ********";

	if (determinant != 0)
	{
		//транспонирование А
		inverse(A, Aobr, size);
	
		out << "\n\nОбратная матрица к матрице А:\n\n";
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				out << setprecision(3) << Aobr[i][j] << "\t";
			}
			out << "\n";
		}

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				x[i] += Aobr[i][j] * b[j];
			}
		}
	}

	double** edinichnaya = new double* [size];
	for (int i = 0; i < size; i++)
	{
		edinichnaya[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			edinichnaya[i][j] = 0;
			for (int k = 0; k < size; k++)
			{
				edinichnaya[i][j] += Aobr[i][k] * A[k][j];
			}
		}
	}

	out << "\n\nПроверка обратной матрицы:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << edinichnaya[i][j] << "\t";
		}
		out << "\n";
	}

	//решение методом псевдообратной матрицы Мура-Пенроуза

	double* MatrA = new double[size * size];
	double* MatrU = new double[size * size];
	double* MatrV = new double[size * size];
	double* VectorSigm = new double[size];
	double* MatrSigm = new double[size * size];
	double* TranspMatrSigm = new double[size * size];

	double* Buffer = new double[size * size];
	double* Buffer1 = new double[size];
	double* Buffer2 = new double[size * size];
	double* Buffer3 = new double[size * size];

	for (int i = 0; i < size; i++)
	{
		VectorSigm[i] = 0;
		Buffer1[i] = 0;
		for (int j = 0; j < size; j++)
		{
			MatrA[i * size + j] = A[i][j];
			MatrU[i * size + j] = 0;
			MatrV[i * size + j] = 0;
			MatrSigm[i * size + j] = 0;
			TranspMatrSigm[i * size + j] = 0;
			Buffer[i * size + j] = 0;
			Buffer2[i * size + j] = 0;
			Buffer3[i * size + j] = 0;
		}
	}

	out << "\n************** ПСЕВДООБРАТНАЯ МАТРИЦА (СВД) **************";

	svd_hestenes(size, size, MatrA, MatrU, MatrV, VectorSigm);

	double porog = VectorSigm[0] / 1000000;

	//транспонируем матрицу U
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			Buffer[i * size + j] = MatrU[i + j * size];
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
				MatrSigm[i * size + j] = VectorSigm[i];
			else
				MatrSigm[i * size + j] = 0;
		}
	}

	//транспонируем матрицу Sigma
	for (int i = 0; i < size; i++)
	{
		if (VectorSigm[i] <= porog)
			Buffer1[i] = 0;
		else
			Buffer1[i] = 1 / VectorSigm[i];
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
				TranspMatrSigm[i * size + j] = Buffer1[i];
			else
				TranspMatrSigm[i * size + j] = 0;
		}
	}

	//умножаем V на транспонированную Sigma
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				Buffer2[i * size + j] += MatrV[i * size + k] * TranspMatrSigm[k * size + j];
			}
		}
	}

	//умножаем результат предыдущего действия на транспонированную матрицу U (определение псевдообратной матрицы A-крест)
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				Buffer3[i * size + j] += Buffer2[i * size + k] * Buffer[k * size + j];
			}
		}
	}

	out << "\n\nВектор Сигма:\n\n";
	for (int i = 0; i < size; i++)
	{
		out << setprecision(3) << VectorSigm[i] << "\t";
	}

	out << "\n\nМатрица Сигма:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << setprecision(3) << MatrSigm[i * size + j] << "\t";
		}
		out << "\n";
	}

	out << "\nТранспонированная матрица Сигма:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << setprecision(3) << TranspMatrSigm[i * size + j] << "\t";
		}
		out << "\n";
	}

	out << "\nМатрица U:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << setprecision(3) << MatrU[i * size + j] << "\t";
		}
		out << "\n";
	}

	out << "\nМатрица V:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << setprecision(3) << MatrV[i * size + j] << "\t";
		}
		out << "\n";
	}

	out << "\nПроизведение VSigma^-1:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << setprecision(3) << Buffer2[i * size + j] << "\t";
		}
		out << "\n";
	}

	out << "\nПсевдообратная матрица Мура-Пенроуза (A+):\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << setprecision(3) << Buffer3[i * size + j] << "\t";
		}
		out << "\n";
	}

	//число обусловленности

	double cond_rand = 0;
	cond_rand = VectorSigm[0] / VectorSigm[size - 1];

	out << "\nЧисло обусловленности матрицы А: " << cond_rand;

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			x_psevdo[i] += Buffer3[i * size + j] * b[j];
		}
	}

	out << "\n\n************************** РЕШЕНИЕ **************************";

	out << "\n\nРешение СЛАУ методом обратной матрицы:\n\n";
	for (int i = 0; i < size; i++)
	{
		out << setprecision(3) << x[i] << '\n';
	}

	out << "\n\nРешение X через псевдообратную матрицу:\n\n";
	for (int i = 0; i < size; i++)
	{
		out << setprecision(3) << x_psevdo[i] << "\n";
	}
	
	double nevyazki_obr = 0;
	double nevyazki_psevdo = 0;
	double* sum_ax_obr = new double[size];
	double* sum_ax_psevdo = new double[size];
	for (int i = 0; i < size; i++)
	{
		sum_ax_obr[i] = 0;
		sum_ax_psevdo[i] = 0;
		for (int j = 0; j < size; j++)
		{
			sum_ax_obr[i] += A[i][j] * x[j];
			sum_ax_psevdo[i] += A[i][j] * x_psevdo[j];
		}
	}

	for (int i = 0; i < size; i++)
	{
		nevyazki_obr += (b[i] - sum_ax_obr[i]) * (b[i] - sum_ax_obr[i]);
		nevyazki_psevdo += (b[i] - sum_ax_psevdo[i]) * (b[i] - sum_ax_psevdo[i]);
	}
	nevyazki_obr /= size;
	nevyazki_psevdo /= size;

	out << "\nНевязки системы для обратной матрицы: " << nevyazki_obr;
	out << "\n\nНевязки системы для псевдообратной матрицы: " << nevyazki_psevdo;

	nevyaz_obr = nevyazki_obr;
	nevyaz_psevdo = nevyazki_psevdo;

	determ = determinant;

	dim_m = size;
	dim_n = size;

	out.close();

	CFile file;
	file.Open(L"result_part_1.txt", CFile::modeRead);
	CStringA str;
	LPSTR pBuf = str.GetBuffer(file.GetLength() + 1);
	file.Read(pBuf, file.GetLength() + 1);
	pBuf[file.GetLength()] = NULL;
	CStringA decodedText = str;
	file.Close();
	str.ReleaseBuffer();

	result_part1 = "";
	result_part1 = str;

	UpdateData(FALSE);

	for (int i = 0; i < size; ++i)
	{
		delete[] A[i];
	}
	delete[] A;

	for (int i = 0; i < size; ++i)
	{
		delete[] Aobr[i];
	}
	delete[] Aobr;

	for (int i = 0; i < size; ++i)
	{
		delete[] edinichnaya[i];
	}
	delete[] edinichnaya;

	delete[] MatrA;
	delete[] MatrU;
	delete[] MatrSigm;
	delete[] MatrV;
	delete[] Buffer;
	delete[] Buffer1;
	delete[] Buffer2;
	delete[] Buffer3;
	delete[] VectorSigm;
	delete[] TranspMatrSigm;
	delete[] b;
	delete[] x;
	delete[] x_psevdo;

	delete[] sum_ax_obr;
	delete[] sum_ax_psevdo;
}

double CKarhunenLoeveBasisDlg::GaussSignal(int t)
{
	UpdateData(TRUE);
	double result = 0.0;
	double ampl[] = { ampl1, ampl2, ampl3};
	double disp[] = { sigma1, sigma2, sigma3, };
	double t0[] = { t1, t2, t3 };
	for (int i = 0; i < 3; i++)
	{
		result += ampl[i] * exp(-((t - t0[i]) / disp[i]) * ((t - t0[i]) / disp[i]));
	}
	return result;
}

double CKarhunenLoeveBasisDlg::PoligarmSignal(int t)
{
	UpdateData(TRUE);
	double result = 0.0;
	double ampl[] = { ampl1, ampl2, ampl3 };
	double frec[] = { omega1, omega2, omega3, };
	double phase[] = { fi1, fi2, fi3 };
	for (int i = 0; i < 3; i++)
	{
		result += ampl[i] * sin(2 * Pi * frec[i] * t + phase[i]);
	}
	return result;
}

double CKarhunenLoeveBasisDlg::ExpSignal(int t)
{
	UpdateData(TRUE);
	double result = 0.0;
	double ampl[] = { ampl1, ampl2, ampl3 };
	double frec[] = { omega1, omega2, omega3, };
	double phase[] = { fi1, fi2, fi3 };
	for (int i = 0; i < 3; i++)
	{
		result += exp(-ampl[i] * t) * sin(2 * Pi * frec[i] * t + phase[i]);
	}
	return result;
}

void CKarhunenLoeveBasisDlg::OnBnClickedButtonSignal()
{
	// TODO: добавьте свой код обработчика уведомлений
	UpdateData(TRUE);

	memset(DefSignal, 0, Length * sizeof(double));

	if (m_radio_poligarm_signal.GetCheck() == BST_CHECKED)
	{
		double* signal = new double[Length];
		memset(signal, 0, Length * sizeof(double));

		for (int i = 0; i < Length; ++i)
		{
			signal[i] = PoligarmSignal(i);
			DefSignal[i] = signal[i];
		}

		Mashtab(signal, Length, &ymin, &ymax);
		xmax = Length;
		RedrawSignal(ymin, ymax);

		PicDc->SelectObject(&signal_pen);
		PicDc->MoveTo(DOTS(0, signal[0]));

		for (int i = 0; i < Length; ++i)
		{
			PicDc->LineTo(DOTS(i, signal[i]));
		}

		delete[] signal;
	}

	if (m_radio_gauss_signal.GetCheck() == BST_CHECKED)
	{
		double* signal = new double[Length];
		memset(signal, 0, Length * sizeof(double));

		for (int i = 0; i < Length; ++i)
		{
			signal[i] = GaussSignal(i);
			DefSignal[i] = signal[i];
		}

		Mashtab(signal, Length, &ymin, &ymax);
		RedrawSignal(-0.15 * ymax, ymax);

		PicDc->SelectObject(&signal_pen);
		PicDc->MoveTo(DOTS(0, signal[0]));

		for (int i = 0; i < Length; ++i)
		{
			PicDc->LineTo(DOTS(i, signal[i]));
		}

		delete[] signal;
	}

	if (m_radio_exp_signal.GetCheck() == BST_CHECKED)
	{
		double* signal = new double[Length];
		memset(signal, 0, Length * sizeof(double));

		for (int i = 0; i < Length; ++i)
		{
			signal[i] = ExpSignal(i);
			DefSignal[i] = signal[i];
		}

		Mashtab(signal, Length, &ymin, &ymax);
		RedrawSignal(-ymax, ymax);

		PicDc->SelectObject(&signal_pen);
		PicDc->MoveTo(DOTS(0, signal[0]));

		for (int i = 0; i < Length; ++i)
		{
			PicDc->LineTo(DOTS(i, signal[i]));
		}

		delete[] signal;
	}

	if (m_radio_poligarm_signal.GetCheck() == BST_UNCHECKED && m_radio_gauss_signal.GetCheck() == BST_UNCHECKED && m_radio_exp_signal.GetCheck() == BST_UNCHECKED)
	{
		MessageBox(L"Пожалуйста, выберите сигнал для построения!", L"Information", MB_ICONINFORMATION | MB_OK);
	}

	if (m_radio_poligarm_signal.GetCheck() == BST_CHECKED || m_radio_gauss_signal.GetCheck() == BST_CHECKED || m_radio_exp_signal.GetCheck() == BST_CHECKED)
	{

		//AKП

		double* Rxx = new double[dim_akm];
		double* R = new double[dim_akm * dim_akm];
		memset(Rxx, 0, dim_akm * sizeof(double));
		double summ;

		for (int m = 0; m < dim_akm; m++)
		{
			summ = 0;
			Rxx[m] = 0;
			for (int k = 0; k < dim_akm - m; k++)
			{
				R[m * dim_akm + k] = 0;
				summ += DefSignal[k] * DefSignal[k + m];
			}
			Rxx[m] = summ / (dim_akm - m);
		}

		//AKM

		int index = 0;
		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				if (i - j >= 0)
				{
					R[i * dim_akm + j] = Rxx[i - j];
				}
				else
				{
					index = dim_akm + (i - j);
					R[i * dim_akm + j] = Rxx[index];
				}
			}
		}

		ofstream out("matrix_part_2.txt");

		out << "АКМ(размерность " << dim_akm << ")\n\n";
		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				out << setprecision(2) << R[i * dim_akm + j] << "\t";
			}
			out << "\n\n";
		}

		double* MatrA = new double[dim_akm * dim_akm];
		double* MatrU = new double[dim_akm * dim_akm];
		double* MatrV = new double[dim_akm * dim_akm];
		double* VectorSigm = new double[dim_akm];
		double* MatrSigm = new double[dim_akm * dim_akm];
		double* TranspMatrSigm = new double[dim_akm * dim_akm];

		double* Buffer = new double[dim_akm * dim_akm];
		double* Buffer1 = new double[dim_akm];
		double* Buffer2 = new double[dim_akm * dim_akm];
		double* Buffer3 = new double[dim_akm * dim_akm];

		for (int i = 0; i < dim_akm; i++)
		{
			VectorSigm[i] = 0;
			Buffer1[i] = 0;
			for (int j = 0; j < dim_akm; j++)
			{
				MatrA[i * dim_akm + j] = R[i * dim_akm + j];
				MatrU[i * dim_akm + j] = 0;
				MatrV[i * dim_akm + j] = 0;
				MatrSigm[i * dim_akm + j] = 0;
				TranspMatrSigm[i * dim_akm + j] = 0;
				Buffer[i * dim_akm + j] = 0;
				Buffer2[i * dim_akm + j] = 0;
				Buffer3[i * dim_akm + j] = 0;
			}
		}

		svd_hestenes(dim_akm, dim_akm, MatrA, MatrU, MatrV, VectorSigm);

		//транспонируем матрицу U
		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				Buffer[i * dim_akm + j] = MatrU[i + j * dim_akm];
			}
		}

		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				if (i == j)
					MatrSigm[i * dim_akm + j] = VectorSigm[i];
				else
					MatrSigm[i * dim_akm + j] = 0;
			}
		}

		//транспонируем матрицу Sigma
		for (int i = 0; i < dim_akm; i++)
		{
			Buffer1[i] = 1 / VectorSigm[i];
		}

		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				if (i == j)
					TranspMatrSigm[i * dim_akm + j] = Buffer1[i];
				else
					TranspMatrSigm[i * dim_akm + j] = 0;
			}
		}

		//умножаем V на транспонированную Sigma
		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				for (int k = 0; k < dim_akm; k++)
				{
					Buffer2[i * dim_akm + j] += MatrV[i * dim_akm + k] * TranspMatrSigm[k * dim_akm + j];
				}
			}
		}

		//умножаем результат предыдущего действия на транспонированную матрицу U (определение псевдообратной матрицы A-крест)
		for (int i = 0; i < dim_akm; i++)
		{
			for (int j = 0; j < dim_akm; j++)
			{
				for (int k = 0; k < dim_akm; k++)
				{
					Buffer3[i * dim_akm + j] += Buffer2[i * dim_akm + k] * Buffer[k * dim_akm + j];
				}
			}
		}

		out << "Вектор Сигма:\n\n";
		for (int i = 0; i < dim_akm; i++)
		{
			out << setprecision(3) << VectorSigm[i] << "\t";
		}

		Mashtab(VectorSigm, dim_akm, &ymin_sing, &ymax_sing);
		RedrawSingular(-0.15 * ymax_sing, ymax_sing);

		PicDcSingular->SelectObject(&singular_value_pen);
		PicDcSingular->MoveTo(DOTSSING(0, VectorSigm[0]));

		for (int i = 0; i < dim_akm; ++i)
		{
			PicDcSingular->LineTo(DOTSSING(i, VectorSigm[i]));
		}

		double* Vector = new double[dim_akm];
		for (int i = 0; i < dim_akm; i++)
		{
			Vector[i] = 0;
			for (int j = 0; j < dim_akm; j++)
			{
				Vector[i] = Buffer[i + number_of_vector * j];
			}
		}

		Mashtab(Vector, dim_akm, &ymin_v, &ymax_v);
		RedrawOwnVector(-ymax_v, ymax_v);

		PicDcV->SelectObject(&own_vector_pen);
		PicDcV->MoveTo(DOTSV(0, Vector[0]));

		for (int i = 0; i < dim_akm; ++i)
		{
			PicDcV->LineTo(DOTSV(i, Vector[i]));
		}

		out << "Собственный вектор:\n\n";
		for (int i = 0; i < dim_akm; i++)
		{
			out << Vector[i] << "\t";
		}

		delete[] Rxx;
		delete[] R;

		delete[] Vector;
		delete[] MatrA;
		delete[] MatrU;
		delete[] MatrSigm;
		delete[] MatrV;
		delete[] Buffer;
		delete[] Buffer1;
		delete[] Buffer2;
		delete[] Buffer3;
		delete[] VectorSigm;
		delete[] TranspMatrSigm;
	}
}
