
// KarhunenLoeveBasisDlg.cpp: файл реализации
//

#include "pch.h"
#include "framework.h"
#include "KarhunenLoeveBasis.h"
#include "KarhunenLoeveBasisDlg.h"
#include "afxdialogex.h"

#include <math.h>
#include <iostream>
#include <fstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace std;

// Диалоговое окно CKarhunenLoeveBasisDlg



CKarhunenLoeveBasisDlg::CKarhunenLoeveBasisDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_KARHUNENLOEVEBASIS_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CKarhunenLoeveBasisDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CKarhunenLoeveBasisDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON_SVD_PART1, &CKarhunenLoeveBasisDlg::OnBnClickedButtonSvdPart1)
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

float CKarhunenLoeveBasisDlg::Psi()
{
	float r = 0;
	for (int i = 1; i <= 12; i++)
	{
		r += ((rand() % 100) / (100 * 1.0) * 2) - 1;		// [-1;1]
	}
	return r / 12;
}

// ----------------------------- Определитель --------------------------
// Вычисляет определитель матрицы T размерностью N
float CKarhunenLoeveBasisDlg::det(float** T, UINT32 N)
{
	float det__;
	int sub_j, s;
	float** subT;    // Субматрица как набор ссылок на исходную матрицу
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
		subT = new float* [N - 1];  // Массив ссылок на столбцы субматрицы
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
};

void CKarhunenLoeveBasisDlg::OnBnClickedButtonSvdPart1()
{
	// TODO: добавьте свой код обработчика уведомлений
	const int M = 9;
	float* koeff = new float[M];
	srand(time(0));
	int size = 3;
	float** matrA;
	matrA = new float* [size];
	for (int i = 0; i < size; ++i)
	{
		matrA[i] = new float[size];
	}
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			matrA[i][j] = Psi();
		}
	}

	ofstream out("Матрица коэффициентов А.txt");
	out << "Матрица коэффициентов А:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << matrA[i][j] << " ";
		}
		out << "\n";
	}

	float determinant = det(matrA, 3);

	out << "\n\nОпределитель матрицы: " << determinant;

	out.close();

	for (int i = 0; i < size; ++i)
	{
		delete[] matrA[i];
	}
	delete[] matrA;
}
