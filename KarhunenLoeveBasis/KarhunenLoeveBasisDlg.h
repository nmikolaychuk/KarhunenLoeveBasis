
// KarhunenLoeveBasisDlg.h: файл заголовка
//

#pragma once


// Диалоговое окно CKarhunenLoeveBasisDlg
class CKarhunenLoeveBasisDlg : public CDialogEx
{
// Создание
public:
	CKarhunenLoeveBasisDlg(CWnd* pParent = nullptr);	// стандартный конструктор

	int svd_hestenes(int, int, float*, float*, float*, float*);
	float det(float** T, UINT32 N);

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
	int size = 3;
	CString result_part1;
};
