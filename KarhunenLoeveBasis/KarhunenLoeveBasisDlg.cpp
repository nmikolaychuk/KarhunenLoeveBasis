
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

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

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
	, result_part1(_T(""))
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CKarhunenLoeveBasisDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_RESULT_TEXT, result_part1);
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

int CKarhunenLoeveBasisDlg::svd_hestenes(int m_m, int n_n, float* a, float* u, float* v, float* sigma)
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

float getRandomFloat(PRNG& generator, float minValue, float maxValue)
{
	// Проверяем корректность аргументов
	assert(minValue < maxValue);

	// Создаём распределение
	std::uniform_real_distribution<float> distribution(minValue, maxValue);

	// Вычисляем псевдослучайное число: вызовем распределение как функцию,
	//  передав генератор произвольных целых чисел как аргумент.
	return distribution(generator.engine);
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
}

// Функция, производящая обращение матрицы.
// Принимает:
//     matrix - матрица для обращения
//     result - матрица достаточного размера для вмещения результата
//     size   - размерность матрицы
// Возвращает:
//     true в случае успешного обращения, false в противном случае
bool inverse(float** matrix, float** result, int size)
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
	float** copy = new float* [size]();

	// Заполняем копию исходной матрицы
	for (int i = 0; i < size; ++i)
	{
		copy[i] = new float[size];

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
			float multi = copy[i][k];

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

	float** A = new float*[size];
	for (int i = 0; i < size; i++)
	{
		A[i] = new float[size];
	}

	float** Aobr = new float* [size];
	for (int i = 0; i < size; i++)
	{
		Aobr[i] = new float[size];
	}
	
	float* b = new float[size];
	float* x = new float[size];
	float* x_psevdo = new float[size];

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
		b[i] = getRandomFloat(generator, -1, 1);
	}

	ofstream out("result_part_1.txt");

	out << "\t\t\tИНФОРМАЦИОННЫЕ ТЕХНОЛОГИИ\n";
	out << "\tЗадание 3. Построение и визуализация базиса Каруена-Лоэва. Часть 1.";

	out << "\n\nМатрица коэффициентов А:\n\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			out << A[i][j] << "\t";
		}
		out << "\n";
	}

	float determinant = det(A, size);

	out << "\nОпределитель матрицы: " << determinant;

	out << "\n\nВектор свободных членов: ";
	for (int i = 0; i < size; i++)
	{
		out << setprecision(3) << b[i] << "\t";
	}

	out << "\n\n******************* ОБРАТНАЯ МАТРИЦА (КЛАССИЧЕСКИЙ МЕТОД) *******************";

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

	//решение методом псевдообратной матрицы Мура-Пенроуза

	float* MatrA = new float[size * size];
	float* MatrU = new float[size * size];
	float* MatrV = new float[size * size];
	float* VectorSigm = new float[size];
	float* MatrSigm = new float[size * size];
	float* TranspMatrSigm = new float[size * size];

	float* Buffer = new float[size * size];
	float* Buffer1 = new float[size];
	float* Buffer2 = new float[size * size];
	float* Buffer3 = new float[size * size];

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

	out << "\n************************ ПСЕВДООБРАТНАЯ МАТРИЦА (СВД) ************************";

	svd_hestenes(size, size, MatrA, MatrU, MatrV, VectorSigm);

	float porog = VectorSigm[0] / 1000;

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

	out << "\n\nМатрица U:\n\n";
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

	float cond_rand = 0;
	cond_rand = VectorSigm[0] / VectorSigm[size - 1];

	out << "\nЧисло обусловленности матрицы А: " << cond_rand;

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			x_psevdo[i] += Buffer3[i * size + j] * b[j];
		}
	}

	out << "\n\n*********************************** РЕШЕНИЕ ***********************************";

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
	
	float nevyazki_obr = 0;
	float nevyazki_psevdo = 0;
	float* sum_ax_obr = new float[size];
	float* sum_ax_psevdo = new float[size];
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
}
