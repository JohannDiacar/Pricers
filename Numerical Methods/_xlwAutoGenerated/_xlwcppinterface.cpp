//// 
//// Autogenerated by xlw 
//// Do not edit this file, it will be overwritten 
//// by InterfaceGenerator 
////

#include "xlw/MyContainers.h"
#include <xlw/CellMatrix.h>
#include "..\cppinterface.h"
#include <xlw/xlw.h>
#include <xlw/XlFunctionRegistration.h>
#include <stdexcept>
#include <xlw/XlOpenClose.h>
#include <xlw/HiResTimer.h>
using namespace xlw;

namespace {
const char* LibraryName = "MyTestLibrary";
};


// registrations start here


namespace
{
XLRegistration::Arg
EchoShortArgs[]=
{
{ "x"," number to be echoed ","XLF_OPER"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerEchoShort("xlEchoShort",
"EchoShort",
" echoes a short ",
LibraryName,
EchoShortArgs,
1
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlEchoShort(
LPXLFOPER xa)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper xb(
	(xa));
short x(
	xb.AsShort("x"));

short result(
	EchoShort(
		x)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

//////////////////////////
// Methods that will get registered to execute in AutoOpen
//////////////////////////

//////////////////////////
// Methods that will get registered to execute in AutoClose
//////////////////////////

