
// KarhunenLoeveBasisDlg.h: файл заголовка
//

#pragma once


// Диалоговое окно CKarhunenLoeveBasisDlg
class CKarhunenLoeveBasisDlg : public CDialogEx
{
	// Создание
public:
	CKarhunenLoeveBasisDlg(CWnd* pParent = nullptr);	// стандартный конструктор

	int svd_hestenes(int, int, double*, double*, double*, double*);
	double det(double** T, UINT32 N);
	void Mashtab(double arr[], int dim, double* mmin, double* mmax);
	void RedrawSignal(double ymin, double ymax);
	void RedrawSingular(double ymin_sing, double ymax_sing);
	void RedrawOwnVector(double ymin_v, double ymax_v);
	double PoligarmSignal(int t);
	double GaussSignal(int t);
	double ExpSignal(int t);

	// Данные диалогового окна
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_KARHUNENLOEVEBASIS_DIALOG };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// поддержка DDX/DDV


// Реализация
protected:
	HICON m_hIcon;

	// Созданные функции схемы сообщений
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonSvdPart1();

	CWnd* PicWnd;				//области рисования
	CDC* PicDc;
	CRect Pic;

	CWnd* PicWndSingular;				//области рисования
	CDC* PicDcSingular;
	CRect PicSingular;

	CWnd* PicWndV;				//области рисования
	CDC* PicDcV;
	CRect PicV;

	CPen osi_pen;				//ручки
	CPen setka_pen;
	CPen signal_pen;
	CPen singular_value_pen;
	CPen own_vector_pen;

	int size = 3;
	int Length;
	double nevyaz_obr;
	double nevyaz_psevdo;
	double determ;
	int dim_m;
	int dim_n;
	double ampl1;
	double ampl2;
	double ampl3;
	double omega1;
	double omega2;
	double omega3;
	double fi1;
	double fi2;
	double fi3;
	double sigma1;
	double sigma2;
	double sigma3;
	double t1;
	double t2;
	double t3;
	int dim_akm;

	double* DefSignal = new double[Length];

	double xp = 0, yp = 0,				//коэфициенты пересчета
		xmin = -Length * 0.1, xmax = Length,				//максисимальное и минимальное значение х 
		ymin = -0.001, ymax = 0.006;				//максисимальное и минимальное значение y

	double xp_sing = 0, yp_sing = 0,				//коэфициенты пересчета
		xmin_sing = -dim_akm * 0.1, xmax_sing = dim_akm,				//максисимальное и минимальное значение х 
		ymin_sing = -0.001, ymax_sing = 0.006;				//максисимальное и минимальное значение y

	double xp_v = 0, yp_v = 0,				//коэфициенты пересчета
		xmin_v = -dim_akm * 0.1, xmax_v = dim_akm,				//максисимальное и минимальное значение х 
		ymin_v = -0.001, ymax_v = 0.006;

	double Pi = 3.141592653589;

	CString result_part1;

	CButton m_radio_poligarm_signal;
	CButton m_radio_gauss_signal;
	CButton m_radio_exp_signal;
	afx_msg void OnBnClickedButtonSignal();
//	double number_of_vector;
	int number_of_vector;
};
