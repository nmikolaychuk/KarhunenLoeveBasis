
// KarhunenLoeveBasis.h: главный файл заголовка для приложения PROJECT_NAME
//

#pragma once

#ifndef __AFXWIN_H__
	#error "включить pch.h до включения этого файла в PCH"
#endif

#include "resource.h"		// основные символы


// CKarhunenLoeveBasisApp:
// Сведения о реализации этого класса: KarhunenLoeveBasis.cpp
//

class CKarhunenLoeveBasisApp : public CWinApp
{
public:
	CKarhunenLoeveBasisApp();

// Переопределение
public:
	virtual BOOL InitInstance();

// Реализация

	DECLARE_MESSAGE_MAP()
};

extern CKarhunenLoeveBasisApp theApp;
